# ==============================================================================
# C1_network_R0_inverse.R
# Purpose:
#   Given target R0 values per transmission route,
#   back-calculate the per-contact transmission probabilities (p)
# ==============================================================================

library(dplyr)

case_tag <- "case1_1M"

# ==============================================================================
# [Section 1] Load network degree stats
# ==============================================================================

network_dir <- "output/network"
nodes   <- readRDS(file.path(network_dir, sprintf("%s_nodes.rds",            case_tag)))
layer1  <- readRDS(file.path(network_dir, sprintf("%s_layer1_household.rds", case_tag)))
layer2  <- readRDS(file.path(network_dir, sprintf("%s_layer2_community.rds", case_tag)))
layer3a <- readRDS(file.path(network_dir, sprintf("%s_layer3_admission.rds", case_tag)))

# HH degree
hh_deg_df <- bind_rows(
  layer1 %>% select(pid=from),
  layer1 %>% select(pid=to)
) %>% group_by(pid) %>% summarise(hh_deg=n(), .groups="drop")

# Community degree
comm_deg_df <- bind_rows(
  layer2 %>% select(pid=from),
  layer2 %>% select(pid=to)
) %>% group_by(pid) %>% summarise(comm_deg=n(), .groups="drop")

deg_df <- data.frame(pid=nodes$person_id) %>%
  left_join(hh_deg_df,   by="pid") %>%
  left_join(comm_deg_df, by="pid")
deg_df$hh_deg[is.na(deg_df$hh_deg)]     <- 0L
deg_df$comm_deg[is.na(deg_df$comm_deg)] <- 0L

mean_hh_deg   <- mean(deg_df$hh_deg)
mean_comm_deg <- mean(deg_df$comm_deg)

hh_sizes     <- deg_df$hh_deg + 1L
hh_sizes     <- hh_sizes[hh_sizes > 1]
mean_freq_hh <- mean((hh_sizes - 1L) / hh_sizes)  # E[(n-1)/n]

mean_adm_deg <- mean(sapply(layer3a$hcw_list, length))

# ==============================================================================
# [Configuration] — Change ONLY these
# ==============================================================================

# Target R0 per route
R0_hh_target   <- 0.4*1.3
R0_comm_target  <- 0.3*1.3
R0_hosp_target  <- 0.1
R0_fun_target   <- 0.3*1.3

# Fixed severity parameters
prob_hospitalised_genPop <- 0.50
prob_death_comm          <- 0.70
prob_death_hosp          <- 0.50
p_unsafe_funeral         <- 0.70
safe_funeral_efficacy    <- 0.80

# Natural history means
T_incub       <- 8.5
T_onset_hosp  <- 1.0
T_onset_death <- 9.3
T_onset_recov <- 13.0
T_hosp_death  <- 4.5
T_hosp_recov  <- 8.0

# Generation time: Gamma(shape=2.5, mean=15.4)
Tg_shape <- 2.5
Tg_rate  <- Tg_shape / 15.4

# ==============================================================================
# [Section 2] Compute transmission window fractions
# ==============================================================================

# Phase 1 end times
t_phase1_hosp    <- T_incub + T_onset_hosp
t_phase1_no_hosp <- T_incub + (prob_death_comm * T_onset_death +
                                 (1 - prob_death_comm) * T_onset_recov)

p_tg_phase1 <- prob_hospitalised_genPop * pgamma(t_phase1_hosp,    Tg_shape, Tg_rate) +
  (1 - prob_hospitalised_genPop) * pgamma(t_phase1_no_hosp, Tg_shape, Tg_rate)

# Phase 2 window
t_hosp_outcome <- T_hosp_death * prob_death_hosp + T_hosp_recov * (1 - prob_death_hosp)
p_tg_phase2 <- pgamma(t_phase1_hosp + t_hosp_outcome, Tg_shape, Tg_rate) -
  pgamma(t_phase1_hosp,                   Tg_shape, Tg_rate)

p_tg_funeral <- 1.0  # Funeral = one-shot event

# ==============================================================================
# [Section 3] Inverse calculation
# ==============================================================================

# Household (frequency dependent):
#   R_hh = p_hh × E[(n-1)/n] × p_tg_phase1
p_inf_household <- R0_hh_target / (mean_freq_hh * p_tg_phase1)

# Community (density dependent):
#   R_comm = p_comm × mean_comm_deg × p_tg_phase1
p_inf_community <- R0_comm_target / (mean_comm_deg * p_tg_phase1)

# Hospital (frequency dependent, cancels):
#   R_hosp = prob_hosp × p_hosp × p_tg_phase2
p_inf_patient_to_hcw <- R0_hosp_target / (prob_hospitalised_genPop * p_tg_phase2)

# Funeral:
#   R_fun = prob_death × p_eff_fun × (mean_hh_deg + mean_comm_deg) × p_tg_funeral
#   p_eff_fun = p_unsafe × p_fun_unsafe + (1-p_unsafe) × p_fun_unsafe × (1-eff)
#             = p_fun_unsafe × [p_unsafe + (1-p_unsafe)×(1-eff)]
funeral_multiplier <- p_unsafe_funeral + (1 - p_unsafe_funeral) * (1 - safe_funeral_efficacy)
p_inf_funeral_unsafe <- R0_fun_target /
  (prob_death_comm * funeral_multiplier * (mean_hh_deg + mean_comm_deg) * p_tg_funeral)
p_inf_funeral_safe   <- p_inf_funeral_unsafe * (1 - safe_funeral_efficacy)

# ==============================================================================
# [Section 4] Output
# ==============================================================================

R0_total <- R0_hh_target + R0_comm_target + R0_hosp_target + R0_fun_target

cat(sprintf("=== R0 inverse calculator: %s ===\n", case_tag))
cat(sprintf("\n  Target R0 breakdown:\n"))
cat(sprintf("    R0_household : %.3f\n", R0_hh_target))
cat(sprintf("    R0_community : %.3f\n", R0_comm_target))
cat(sprintf("    R0_hospital  : %.3f\n", R0_hosp_target))
cat(sprintf("    R0_funeral   : %.3f\n", R0_fun_target))
cat(sprintf("    R0_total     : %.3f\n", R0_total))

cat(sprintf("\n  Network stats:\n"))
cat(sprintf("    mean HH size     : %.2f\n", mean(hh_sizes)))
cat(sprintf("    E[(n-1)/n] HH    : %.4f\n", mean_freq_hh))
cat(sprintf("    mean comm degree : %.2f\n", mean_comm_deg))
cat(sprintf("    P(Tg in phase1)  : %.4f\n", p_tg_phase1))
cat(sprintf("    P(Tg in phase2)  : %.4f\n", p_tg_phase2))

cat(sprintf("\n  ── Implied transmission probabilities ──\n"))
cat(sprintf("    p_inf_household      : %.4f\n", p_inf_household))
cat(sprintf("    p_inf_community      : %.4f\n", p_inf_community))
cat(sprintf("    p_inf_patient_to_hcw : %.4f\n", p_inf_patient_to_hcw))
cat(sprintf("    p_inf_funeral_unsafe : %.4f\n", p_inf_funeral_unsafe))
cat(sprintf("    p_inf_funeral_safe   : %.4f\n", p_inf_funeral_safe))

# Sanity check: verify by forward calculation
p_eff_funeral <- p_unsafe_funeral * p_inf_funeral_unsafe +
  (1 - p_unsafe_funeral) * p_inf_funeral_safe
R0_check <- p_inf_household      * mean_freq_hh * p_tg_phase1 +
  p_inf_community       * mean_comm_deg * p_tg_phase1 +
  prob_hospitalised_genPop * p_inf_patient_to_hcw * p_tg_phase2 +
  prob_death_comm * p_eff_funeral * (mean_hh_deg + mean_comm_deg) * p_tg_funeral

cat(sprintf("\n  ── Verification (forward R0) : %.4f ──\n", R0_check))
