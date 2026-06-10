# ==============================================================================
# R0 approximation from network structure + transmission parameters
# ==============================================================================

library(dplyr)

case_tag <- "case1_1M"

nodes   <- readRDS(sprintf("output/network/%s_nodes.rds",            case_tag))
layer1  <- readRDS(sprintf("output/network/%s_layer1_household.rds", case_tag))
layer2  <- readRDS(sprintf("output/network/%s_layer2_community.rds", case_tag))
layer3a <- readRDS(sprintf("output/network/%s_layer3_admission.rds", case_tag))

N <- nrow(nodes)

# ==============================================================================
# Compute per-person degrees
# ==============================================================================

# HH degree: number of household contacts per person (= hh_size - 1)
hh_deg_df <- bind_rows(
  layer1 %>% select(pid = from),
  layer1 %>% select(pid = to)
) %>%
  group_by(pid) %>%
  summarise(hh_deg = n(), .groups = "drop")

# Community degree: number of community contacts per person
comm_deg_df <- bind_rows(
  layer2 %>% select(pid = from),
  layer2 %>% select(pid = to)
) %>%
  group_by(pid) %>%
  summarise(comm_deg = n(), .groups = "drop")

# Join to nodes
deg_df <- data.frame(pid = nodes$person_id) %>%
  left_join(hh_deg_df,   by = "pid") %>%
  left_join(comm_deg_df, by = "pid")
deg_df$hh_deg[is.na(deg_df$hh_deg)]     <- 0L
deg_df$comm_deg[is.na(deg_df$comm_deg)] <- 0L

mean_hh_deg   <- mean(deg_df$hh_deg)
mean_comm_deg <- mean(deg_df$comm_deg)

# Admission degree: mean number of HCWs a patient is connected to
mean_adm_deg <- mean(sapply(layer3a$hcw_list, length))


# Transmission parameters (from fiber defaults)
p_inf_household        <- 0.20
p_inf_community        <- 0.02
p_inf_patient_to_hcw   <- 0.03
p_unsafe_funeral       <- 0.70
p_inf_funeral_unsafe   <- 0.01
p_inf_funeral_safe     <- p_inf_funeral_unsafe * (1 - 0.80)  # safe_funeral_efficacy=0.80
prob_hospitalised_genPop <- 0.50
prob_death_comm          <- 0.70
prob_death_hosp          <- 0.50


# ==============================================================================
# R0 components (one infected individual → expected secondary cases)
# Each edge = one Bernoulli trial
# ==============================================================================

# ==============================================================================
# R0 components (one infected individual → expected secondary cases)
# Each edge = one Bernoulli trial, weighted by P(Tg falls in transmission window)
# ==============================================================================

# Generation time distribution: Gamma(shape=2.5, mean=15.4)
Tg_shape <- 2.5
Tg_rate  <- Tg_shape / 15.4

# Natural history means (from fiber defaults)
T_incub        <- 8.5
T_onset_hosp   <- 1.0
T_onset_death  <- 9.3
T_onset_recov  <- 13.0
T_hosp_death   <- 4.5
T_hosp_recov   <- 8.0

# Phase 1 end time (infection → hospitalisation or community outcome)
# For hospitalised: infection + incubation + onset_to_hosp
# For non-hospitalised: infection + incubation + onset_to_outcome
t_phase1_hosp    <- T_incub + T_onset_hosp
t_phase1_no_hosp <- T_incub + (prob_death_comm * T_onset_death +
                                 (1 - prob_death_comm) * T_onset_recov)

# P(Tg falls in Phase 1 window)
p_tg_phase1 <- prob_hospitalised_genPop * pgamma(t_phase1_hosp,    Tg_shape, Tg_rate) +
  (1 - prob_hospitalised_genPop) * pgamma(t_phase1_no_hosp, Tg_shape, Tg_rate)

# Phase 2 window: hospitalisation → outcome
# P(Tg falls in phase2) = P(t_phase1_hosp < Tg < t_phase1_hosp + T_hosp_outcome)
t_hosp_outcome <- T_hosp_death * prob_death_hosp + T_hosp_recov * (1 - prob_death_hosp)
p_tg_phase2 <- pgamma(t_phase1_hosp + t_hosp_outcome, Tg_shape, Tg_rate) -
  pgamma(t_phase1_hosp,                   Tg_shape, Tg_rate)

# Funeral window: after death (community or hospital)
t_death_comm <- T_incub + T_onset_death
t_death_hosp <- T_incub + T_onset_hosp + T_hosp_death
t_death_mean <- (1 - prob_hospitalised_genPop) * t_death_comm +
  prob_hospitalised_genPop * t_death_hosp
# Funeral Tg is short (shape=20, rate=10 → mean=2 days after death)
# P(Tg ~ Gamma(20,10) < some window) ≈ 1 since it's very tight
# So funeral fraction ≈ 1 (all funeral contacts exposed at time of death)
p_tg_funeral <- 1.0

p_eff_funeral <- p_unsafe_funeral * p_inf_funeral_unsafe +
  (1 - p_unsafe_funeral) * p_inf_funeral_safe

# ── Frequency dependent: HH ──────────────────────────────────────────────────
# p per contact = p_inf_household / hh_size
# R_hh = p_inf_household / hh_size × (hh_size - 1)
#       = p_inf_household × (hh_size - 1) / hh_size
# Need E[(n-1)/n] at individual level (not mean(n-1)/mean(n))
hh_sizes      <- deg_df$hh_deg + 1L          # hh_size = degree + self
hh_sizes      <- hh_sizes[hh_sizes > 1]      # exclude isolates
mean_freq_hh  <- mean((hh_sizes - 1L) / hh_sizes)  # E[(n-1)/n]

R_hh <- p_inf_household * mean_freq_hh * p_tg_phase1

# ── Density dependent: community ─────────────────────────────────────────────
R_comm <- p_inf_community * mean_comm_deg * p_tg_phase1

# ── Frequency dependent: hospital ────────────────────────────────────────────
# p per HCW = p_inf_patient_to_hcw / n_hcw
# R_hosp = prob_hosp × Σ(p/n_hcw) over n_hcw contacts = prob_hosp × p_inf_patient_to_hcw
# (frequency dependent: total contribution = p regardless of n_hcw)
R_hosp <- prob_hospitalised_genPop * p_inf_patient_to_hcw * p_tg_phase2

# ── Funeral ──────────────────────────────────────────────────────────────────
R_fun <- prob_death_comm * p_eff_funeral * (mean_hh_deg + mean_comm_deg) * p_tg_funeral

R0 <- R_hh + R_comm + R_hosp + R_fun

cat(sprintf("=== R0 approximation: %s ===\n", case_tag))
cat(sprintf("  N                    : %d\n", N))
cat(sprintf("  mean HH size         : %.2f\n", mean(hh_sizes)))
cat(sprintf("  E[(n-1)/n] HH        : %.4f\n", mean_freq_hh))
cat(sprintf("  mean comm degree     : %.2f\n", mean_comm_deg))
cat(sprintf("  mean adm degree      : %.2f\n", mean_adm_deg))
cat(sprintf("  p_eff_funeral        : %.4f\n", p_eff_funeral))
cat(sprintf("  P(Tg in phase1)      : %.4f\n", p_tg_phase1))
cat(sprintf("  P(Tg in phase2)      : %.4f\n", p_tg_phase2))
cat(sprintf("\n  R_household (freq) : %.3f\n", R_hh))
cat(sprintf("  R_community (dens) : %.3f\n", R_comm))
cat(sprintf("  R_hospital  (freq) : %.3f\n", R_hosp))
cat(sprintf("  R_funeral          : %.3f\n", R_fun))
cat(sprintf("\n  R0 ≈ %.3f\n", R0))
