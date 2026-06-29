# ==============================================================================
# C2_param_p1_R0check_reverse.R
# Purpose:
#   Given target R0 values per transmission route,
#   back-calculate the per-contact transmission probabilities (p_inf_*)
#
# Network structure (new):
#   Layer 1: household (full clique, mass-action)
#            → p_inf_household_close    × close_only_home (age-weighted)
#            → p_inf_household_physical × phys_only_home  (age-weighted)
#   Layer 2: community — 6 sub-layers (close/physical × daily/weekly/monthly)
#            → per-edge, frequency dependent
#   Layer 3: hospital
#            → patient→HCW: mass-action (÷ n_hcw)
#            → HCW→HCW:     mass-action (÷ n_hcw)
#   Funeral: NegBin(funeral_avg, funeral_k) attendees
#            → household attendees (all): p_inf_funeral_household
#            → community attendees (funeral_avg - mean_hh_deg): p_inf_funeral_community
# ==============================================================================

library(dplyr)

# ==============================================================================
# [Configuration]
# ==============================================================================

case_tag    <- "case1_1M"
network_dir <- "output/network"
mats_path   <- "output/MPMmat/DRC_network_input_matrices.rds"

# ── Target R0 per route ──────────────────────────────────────────────────────
R0_hh_close_target    <- 0.0
R0_hh_physical_target <- 0.5
R0_comm_target        <- 0.2
R0_hosp_target        <- 881/(28600 - 881)*0.3
R0_fun_target         <- 0.5

# ── Ratio of physical to close p_inf within community ────────────────────────
# Used to split R0_comm into close and physical components.
# physical contacts are assumed to be more infectious than close.
phys_to_close_ratio   <- 10^6  # p_inf_physical = phys_to_close_ratio × p_inf_close

# ── Fixed severity parameters ─────────────────────────────────────────────────
prob_hospitalised_genPop <- 0.30
prob_death_comm          <- 0.50
prob_death_hosp          <- 0.20
p_unsafe_funeral         <- 0.50
funeral_unsafe_multiplier <- 1.0
funeral_safe_multiplier   <- 0.1

# ── Funeral ───────────────────────────────────────────────────────────────────
funeral_avg <- 20   # NegBin mean attendee count

# ── Natural history means (days) ──────────────────────────────────────────────
T_incub       <- 9.0
T_onset_hosp  <- 5.0
T_onset_death <- 12.0
T_onset_recov <- 18.0
T_hosp_death  <- 8.0
T_hosp_recov  <- 13.0

# Generation time: Gamma(shape, rate)
Tg_shape <- 2.5
Tg_rate  <- Tg_shape / 12.0  # mean = 12 days

# ==============================================================================
# [Section 1] Load network data
# ==============================================================================

message("=== Loading network data ===")

nodes   <- readRDS(file.path(network_dir, sprintf("%s_nodes.rds",            case_tag)))
layer1  <- readRDS(file.path(network_dir, sprintf("%s_layer1_household.rds", case_tag)))
layer2d <- readRDS(file.path(network_dir, sprintf("%s_layer2_daily.rds",     case_tag)))
layer2w <- readRDS(file.path(network_dir, sprintf("%s_layer2_weekly.rds",    case_tag)))
layer2m <- readRDS(file.path(network_dir, sprintf("%s_layer2_monthly.rds",   case_tag)))
layer3a <- readRDS(file.path(network_dir, sprintf("%s_layer3_admission.rds", case_tag)))
mats    <- readRDS(mats_path)

close_only_home <- mats$close_only_home  # 16x16
phys_only_home  <- mats$phys_only_home   # 16x16

age_lookup <- nodes %>% select(person_id, age_group)

# ==============================================================================
# [Section 2] Compute network degree statistics
# ==============================================================================

message("=== Computing degree statistics ===")

# ── Household degree ──────────────────────────────────────────────────────────
hh_deg_df <- bind_rows(
  layer1 %>% select(pid = from),
  layer1 %>% select(pid = to)
) %>% group_by(pid) %>% summarise(hh_deg = n(), .groups = "drop")

deg_df <- data.frame(pid = nodes$person_id) %>%
  left_join(hh_deg_df, by = "pid")
deg_df$hh_deg[is.na(deg_df$hh_deg)] <- 0L

mean_hh_deg <- mean(deg_df$hh_deg)
hh_sizes    <- deg_df$hh_deg + 1L
hh_sizes    <- hh_sizes[hh_sizes > 1L]

# E[(n-1)/n]: fraction of household members each person can infect
# (mass-action: contacts / household size)
# mean_freq_hh <- mean((hh_sizes - 1L) / hh_sizes)
mean_freq_hh <- mean((hh_sizes - 1L))

# ── Community degree by sub-layer ─────────────────────────────────────────────
compute_mean_deg <- function(layer) {
  d <- bind_rows(
    layer %>% select(pid = from),
    layer %>% select(pid = to)
  ) %>% group_by(pid) %>% summarise(deg = n(), .groups = "drop")
  all_pids <- data.frame(pid = nodes$person_id) %>%
    left_join(d, by = "pid")
  all_pids$deg[is.na(all_pids$deg)] <- 0L
  mean(all_pids$deg)
}

mean_deg_close_daily    <- compute_mean_deg(layer2d %>% filter(is_physical == 0L))
mean_deg_close_weekly   <- compute_mean_deg(layer2w %>% filter(is_physical == 0L))
mean_deg_close_monthly  <- compute_mean_deg(layer2m %>% filter(is_physical == 0L))
mean_deg_phys_daily     <- compute_mean_deg(layer2d %>% filter(is_physical == 1L))
mean_deg_phys_weekly    <- compute_mean_deg(layer2w %>% filter(is_physical == 1L))
mean_deg_phys_monthly   <- compute_mean_deg(layer2m %>% filter(is_physical == 1L))

mean_deg_comm_total <- mean_deg_close_daily  + mean_deg_close_weekly  +
  mean_deg_close_monthly + mean_deg_phys_daily   +
  mean_deg_phys_weekly   + mean_deg_phys_monthly

# ── Household contact matrix scaling factors (population-average) ─────────────
# For mass-action household transmission:
#   E[p_eff_hh] = p_inf_hh_close × mean(close_only_home) + p_inf_hh_phys × mean(phys_only_home)
# Use mean over all nonzero age-group pairs as scaling factor
mean_close_home <- mean(close_only_home[close_only_home > 0])
mean_phys_home  <- mean(phys_only_home[phys_only_home > 0])

# ── Hospital ──────────────────────────────────────────────────────────────────
mean_adm_deg <- mean(sapply(layer3a$hcw_list, length))

# ── Funeral ───────────────────────────────────────────────────────────────────
# E[hh attendees]   ≈ mean_hh_deg  (all attend)
# E[comm attendees] ≈ funeral_avg - mean_hh_deg
E_hh_att   <- mean_hh_deg
E_comm_att <- max(0, funeral_avg - mean_hh_deg)

message(sprintf("  mean HH deg         : %.3f", mean_hh_deg))
message(sprintf("  E[(n-1)/n] HH       : %.4f", mean_freq_hh))
message(sprintf("  mean_close_home     : %.4f", mean_close_home))
message(sprintf("  mean_phys_home      : %.4f", mean_phys_home))
message(sprintf("  mean deg close daily  : %.3f", mean_deg_close_daily))
message(sprintf("  mean deg close weekly : %.3f", mean_deg_close_weekly))
message(sprintf("  mean deg close monthly: %.3f", mean_deg_close_monthly))
message(sprintf("  mean deg phys daily   : %.3f", mean_deg_phys_daily))
message(sprintf("  mean deg phys weekly  : %.3f", mean_deg_phys_weekly))
message(sprintf("  mean deg phys monthly : %.3f", mean_deg_phys_monthly))
message(sprintf("  mean comm deg total : %.3f", mean_deg_comm_total))
message(sprintf("  mean adm HCWs       : %.2f", mean_adm_deg))
message(sprintf("  E[hh att funeral]   : %.2f", E_hh_att))
message(sprintf("  E[comm att funeral] : %.2f", E_comm_att))

# ==============================================================================
# [Section 3] Transmission window fractions via generation time CDF
# ==============================================================================

# Phase 1 end: onset → hospitalisation (if hosp) or onset → outcome (if not)
t_phase1_hosp    <- T_incub + T_onset_hosp
t_phase1_no_hosp <- T_incub + (prob_death_comm * T_onset_death +
                                 (1 - prob_death_comm) * T_onset_recov)

# P(generation time falls in Phase 1) — weighted by hospitalisation prob
p_tg_phase1 <- prob_hospitalised_genPop *
  pgamma(t_phase1_hosp,    Tg_shape, Tg_rate) +
  (1 - prob_hospitalised_genPop) *
  pgamma(t_phase1_no_hosp, Tg_shape, Tg_rate)

# Phase 2 window (hospitalisation → outcome)
t_hosp_outcome <- T_hosp_death * prob_death_hosp +
  T_hosp_recov  * (1 - prob_death_hosp)
p_tg_phase2 <- pgamma(t_phase1_hosp + t_hosp_outcome, Tg_shape, Tg_rate) -
  pgamma(t_phase1_hosp,                   Tg_shape, Tg_rate)

# Funeral: one-shot event at time of death
p_tg_funeral <- 1.0

message(sprintf("\n  P(Tg in Phase 1)    : %.4f", p_tg_phase1))
message(sprintf("  P(Tg in Phase 2)    : %.4f", p_tg_phase2))

# ==============================================================================
# [Section 4] Reverse calculation
# ==============================================================================

# ── Household ─────────────────────────────────────────────────────────────────
# R0_hh_close = p_inf_hh_close × mean_close_home × mean_freq_hh × p_tg_phase1
# R0_hh_phys  = p_inf_hh_phys  × mean_phys_home  × mean_freq_hh × p_tg_phase1
p_inf_household_close    <- R0_hh_close_target /
  (mean_close_home * mean_freq_hh * p_tg_phase1)
p_inf_household_physical <- R0_hh_physical_target /
  (mean_phys_home  * mean_freq_hh * p_tg_phase1)

# ── Community ─────────────────────────────────────────────────────────────────
# R0_comm = (p_close × Σ deg_close + p_phys × Σ deg_phys) × p_tg_phase1
# where p_phys = phys_to_close_ratio × p_close (for all strata)
#
# R0_comm = p_close × (Σ deg_close + phys_to_close_ratio × Σ deg_phys) × p_tg_phase1
deg_close_total <- mean_deg_close_daily + mean_deg_close_weekly + mean_deg_close_monthly
deg_phys_total  <- mean_deg_phys_daily  + mean_deg_phys_weekly  + mean_deg_phys_monthly

p_inf_community_close <- R0_comm_target /
  ((deg_close_total + phys_to_close_ratio * deg_phys_total) * p_tg_phase1)
p_inf_community_physical <- phys_to_close_ratio * p_inf_community_close

# Daily/weekly/monthly use same per-contact probability (frequency already in degree)
p_inf_community_close_daily    <- p_inf_community_close
p_inf_community_close_weekly   <- p_inf_community_close
p_inf_community_close_monthly  <- p_inf_community_close
p_inf_community_physical_daily   <- p_inf_community_physical
p_inf_community_physical_weekly  <- p_inf_community_physical
p_inf_community_physical_monthly <- p_inf_community_physical

# ── Hospital ──────────────────────────────────────────────────────────────────
# R0_hosp = prob_hosp × p_inf_patient_to_hcw × p_tg_phase2
# (mass-action: dividing by n_hcw is done inside simulation,
#  so p_inf_patient_to_hcw here is the "total" probability before division)
p_inf_patient_to_hcw <- R0_hosp_target /
  (prob_hospitalised_genPop * p_tg_phase2)

# ── Funeral ───────────────────────────────────────────────────────────────────
# R0_fun = prob_death × funeral_mult_eff ×
#          (p_inf_funeral_household × E_hh_att +
#           p_inf_funeral_community × E_comm_att)
# Assume p_inf_funeral_household = p_inf_funeral_community × hh_to_comm_ratio
# For simplicity, solve with equal p (set hh_to_comm_ratio = 1 or adjust)
hh_to_comm_funeral_ratio <- 2.0  # household attendees assumed 2x more infectious

funeral_mult_eff <- p_unsafe_funeral    * funeral_unsafe_multiplier +
  (1-p_unsafe_funeral) * funeral_safe_multiplier

p_inf_funeral_community <- R0_fun_target /
  (prob_death_comm * funeral_mult_eff *
     (hh_to_comm_funeral_ratio * E_hh_att + E_comm_att))
p_inf_funeral_household <- hh_to_comm_funeral_ratio * p_inf_funeral_community

# ==============================================================================
# [Section 5] Output
# ==============================================================================

R0_total <- R0_hh_close_target + R0_hh_physical_target +
  R0_comm_target + R0_hosp_target + R0_fun_target

cat(sprintf("\n=== R0 inverse calculator: %s ===\n", case_tag))

cat(sprintf("\n  Target R0 breakdown:\n"))
cat(sprintf("    R0_hh_close      : %.3f\n", R0_hh_close_target))
cat(sprintf("    R0_hh_physical   : %.3f\n", R0_hh_physical_target))
cat(sprintf("    R0_community     : %.3f\n", R0_comm_target))
cat(sprintf("    R0_hospital      : %.3f\n", R0_hosp_target))
cat(sprintf("    R0_funeral       : %.3f\n", R0_fun_target))
cat(sprintf("    R0_total         : %.3f\n", R0_total))

cat(sprintf("\n  ── Implied transmission probabilities ──\n"))
cat(sprintf("    p_inf_household_close    : %.5f\n", p_inf_household_close))
cat(sprintf("    p_inf_household_physical : %.5f\n", p_inf_household_physical))
cat(sprintf("    p_inf_community_close    : %.5f  (daily/weekly/monthly)\n",
            p_inf_community_close))
cat(sprintf("    p_inf_community_physical : %.5f  (daily/weekly/monthly)\n",
            p_inf_community_physical))
cat(sprintf("    p_inf_patient_to_hcw     : %.5f\n", p_inf_patient_to_hcw))
cat(sprintf("    p_inf_funeral_household  : %.5f\n", p_inf_funeral_household))
cat(sprintf("    p_inf_funeral_community  : %.5f\n", p_inf_funeral_community))

# ── Verification: forward R0 calculation ─────────────────────────────────────
R0_hh_close_check <- p_inf_household_close    * mean_close_home *
  mean_freq_hh * p_tg_phase1
R0_hh_phys_check  <- p_inf_household_physical * mean_phys_home  *
  mean_freq_hh * p_tg_phase1

R0_comm_check <- (p_inf_community_close    * deg_close_total +
                    p_inf_community_physical * deg_phys_total) * p_tg_phase1

R0_hosp_check <- prob_hospitalised_genPop * p_inf_patient_to_hcw * p_tg_phase2

R0_fun_check  <- prob_death_comm * funeral_mult_eff *
  (p_inf_funeral_household * E_hh_att +
     p_inf_funeral_community * E_comm_att)

R0_total_check <- R0_hh_close_check + R0_hh_phys_check +
  R0_comm_check + R0_hosp_check + R0_fun_check

cat(sprintf("\n  ── Verification (forward R0) ──\n"))
cat(sprintf("    R0_hh_close  (check): %.4f  target: %.3f\n",
            R0_hh_close_check, R0_hh_close_target))
cat(sprintf("    R0_hh_phys   (check): %.4f  target: %.3f\n",
            R0_hh_phys_check,  R0_hh_physical_target))
cat(sprintf("    R0_community (check): %.4f  target: %.3f\n",
            R0_comm_check,     R0_comm_target))
cat(sprintf("    R0_hospital  (check): %.4f  target: %.3f\n",
            R0_hosp_check,     R0_hosp_target))
cat(sprintf("    R0_funeral   (check): %.4f  target: %.3f\n",
            R0_fun_check,      R0_fun_target))
cat(sprintf("    R0_total     (check): %.4f  target: %.3f\n",
            R0_total_check,    R0_total))
