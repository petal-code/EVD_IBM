# ==============================================================================
# C2_sim_p1_testrun.R
# Purpose:
#   Single test run of IBM simulation
#   Requires sim_prep to be built first via C2_sim_p0_simprep.R
# ==============================================================================

source("function/COD_IBM_ebola_sim.R")
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# ==============================================================================
# [Configuration]
# ==============================================================================

case_tag      <- "case1_1M"
index_case_id <- 42L          # NULL = random
sim_seed      <- 42L
max_time      <- 600
max_infected  <- Inf

network_dir <- "output/network"
fig_dir     <- "figure/C2_sim"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# [Section 1] Load data
# ==============================================================================

message(sprintf("=== Loading: %s ===", case_tag))

nodes    <- readRDS(file.path(network_dir, sprintf("%s_nodes.rds",      case_tag)))
sim_prep <- readRDS(file.path(network_dir, sprintf("%s_sim_prep.rds",   case_tag)))

# Household contact matrices from p6
mats_path <- "output/MPMmat/DRC_network_input_matrices.rds"
mats_p6   <- readRDS(mats_path)
close_only_home <- mats_p6$close_only_home  # 16x16
phys_only_home  <- mats_p6$phys_only_home   # 16x16

message(sprintf("  N = %d | HCWs = %d", nrow(nodes), sum(nodes$is_hcw)))

if (!is.null(index_case_id)) {
  if (!index_case_id %in% nodes$person_id)
    stop(sprintf("index_case_id %d not found", index_case_id))
  idx_info <- nodes[nodes$person_id == index_case_id, ]
  message(sprintf("  Index case: person %d | age=%d | HCW=%s | cell=%d",
                  index_case_id, idx_info$age,
                  ifelse(idx_info$is_hcw, "YES", "no"), idx_info$cell_id))
}

# ==============================================================================
# [Section 2] Natural history parameters
# ==============================================================================

make_gamma_fn <- function(mean, sd) {
  shape <- (mean / sd)^2; rate <- mean / sd^2
  function(n) if (n <= 0L) numeric(0) else rgamma(n, shape = shape, rate = rate)
}

# Natural history — matched to fiber DEFAULT_SCALAR_INPUTS
# Source: setup_model_parameters.R
incubation_period_fn           <- make_gamma_fn(8.5,  4.5)  # fiber: mean=8.5  sd=4.5
onset_to_death_fn              <- make_gamma_fn(9.3,  3.0)  # fiber: mean=9.3  sd=3.0
onset_to_recovery_fn           <- make_gamma_fn(13.0, 4.0)  # fiber: mean=13.0 sd=4.0
hospitalisation_to_death_fn    <- make_gamma_fn(4.5,  2.0)  # fiber: mean=4.5  sd=2.0
hospitalisation_to_recovery_fn <- make_gamma_fn(8.0,  2.5)  # fiber: mean=8.0  sd=2.5

# IMPORTANT: raw onset_to_hosp mean = 1.0d (fiber pattern).
# Actual hospitalisation delay = onset_to_hospitalisation_fn(1) x delay_hosp(t)
# where delay_hosp(t) is from the scenario CSV (units: days).
# Do NOT set to expected calendar delay — hospitalisation_delay_factor handles scaling.
onset_to_hospitalisation_fn    <- make_gamma_fn(1.0,  0.35) # fiber: raw=1.0 sd=0.35

# Community/hospital generation time — fiber: shape=2.5, rate=2.5/15.4 (mean=15.4d)
generation_time_fn <- local({
  shape <- 2.5; rate <- 2.5 / 15.4   # mean=15.4d, sd~9.76d
  function(n) if (n <= 0L) numeric(0) else rgamma(n, shape = shape, rate = rate)
})

# Funeral generation time — separate from community (fiber: shape=20, rate=10, mean=2.0d)
# Infections tightly clustered around burial event; very low variance
funeral_generation_time_fn <- local({
  shape <- 20; rate <- 10   # mean=2.0d, sd=0.447d
  function(n) if (n <= 0L) numeric(0) else rgamma(n, shape = shape, rate = rate)
})

# ==============================================================================
# [Section 3] Transmission parameters
# ==============================================================================

p_inf_household_close            <- 0.00
p_inf_household_physical         <- 0.1

# Per single contact-event probability (3-week effective computed internally)
p_inf_community_close_daily      <- 0.00
p_inf_community_close_weekly     <- 0.00
p_inf_community_close_monthly    <- 0.00

p_inf_community_physical_daily   <- 0.0015*21
p_inf_community_physical_weekly  <- 0.00152*3
p_inf_community_physical_monthly <- 0.0015*3/4

p_inf_hcw_to_hcw                 <- 0.02
p_inf_patient_to_hcw             <- 0.02

funeral_avg                      <- 15
funeral_k                        <- 0.3   # fiber: overdisp_offspring_funeral=0.30
p_unsafe_funeral                 <- 0.50   # fallback scalar (overridden by TV comm/hosp fns)
p_inf_funeral_household          <- 0.25
p_inf_funeral_community          <- 0.25
funeral_unsafe_multiplier        <- 1.0
funeral_safe_multiplier          <- 0.20   # fiber: 1 - safe_funeral_efficacy = 1 - 0.80 = 0.20

prob_hospitalised_genPop         <- 0.30   # baseline (overridden by TV below)
prob_hospitalised_hcw            <- 0.60
prob_death_comm                  <- 0.70   # fiber: 0.70
prob_death_hosp                  <- 0.50   # fiber: 0.50

# ==============================================================================
# [Section 4b] Time-varying response parameters
# Loaded directly from fiber scenario CSV used in ABC calibration.
# Source: final_six_scenario_values_revised_methodology.csv
# Scenario: middle_drc_conflict (DRC-like archetype with conflict disruption)
# ==============================================================================

# Helper: piecewise-linear interpolation clamped at endpoints (fiber-identical)
make_time_varying <- function(times, values) {
  stopifnot(length(times) == length(values), length(times) >= 2L,
            !is.unsorted(times, strictly = FALSE))
  fn <- stats::approxfun(times, values, rule = 2)
  class(fn) <- c("time_varying_fn", "function")
  fn
}

clip01 <- function(x) pmin(pmax(x, 0), 1)

# ── Load scenario matrix ──────────────────────────────────────────────────
tv_scenario_id  <- "middle_drc_conflict"   # change to switch scenario
tv_csv_path     <- "data_processed/final_six_scenario_values_revised_methodology.csv"

tv_matrix_full  <- read.csv(tv_csv_path, stringsAsFactors = FALSE)
tv_matrix       <- tv_matrix_full[tv_matrix_full$scenario == tv_scenario_id, ]
tv_matrix       <- tv_matrix[order(tv_matrix$relative_day), ]

if (nrow(tv_matrix) == 0L)
  stop(sprintf("Scenario '%s' not found in %s", tv_scenario_id, tv_csv_path))

message(sprintf("  Loaded TV params: scenario='%s', %d time points (day 0-%.0f)",
                tv_scenario_id, nrow(tv_matrix), max(tv_matrix$relative_day)))

times_tv <- tv_matrix$relative_day

# ── hospital unsafe funeral: ETU-weighted blend ───────────────────────────
# Mirrors fiber build_time_varying_args():
#   p_unsafe_hosp(t) = (1-prop_etu(t)) * prob_unsafe_funeral_hosp
#                    + prop_etu(t)      * prob_unsafe_funeral_etu
p_unsafe_hosp_values <- clip01(
  (1 - tv_matrix$prop_etu) * tv_matrix$prob_unsafe_funeral_hosp +
    tv_matrix$prop_etu     * tv_matrix$prob_unsafe_funeral_etu
)

# ── Build time-varying functions ──────────────────────────────────────────
prob_hospitalised_genPop_fn <- make_time_varying(
  times_tv, clip01(tv_matrix$prob_hosp)
)

# delay_hosp in CSV is in days; used directly as hospitalisation_delay_factor
# (fiber raw onset_to_hosp mean = 1d, so delay_hosp IS the mean day count)
hospitalisation_delay_factor <- make_time_varying(
  times_tv, pmax(tv_matrix$delay_hosp, 0.01)
)

p_unsafe_funeral_comm_fn <- make_time_varying(
  times_tv, clip01(tv_matrix$prob_unsafe_funeral_comm)
)

p_unsafe_funeral_hosp_fn <- make_time_varying(
  times_tv, p_unsafe_hosp_values
)

prop_etu_fn <- make_time_varying(
  times_tv, clip01(tv_matrix$prop_etu)
)

# ipc_helper drives both ETU efficacy uplift and PPE coverage (fiber pattern)
ipc_index_fn <- make_time_varying(
  times_tv, clip01(tv_matrix$ipc_helper)
)

# PPE coverage = ipc_helper curve (fiber: ppe_coverage_hcw = ipc_helper)
ppe_coverage_fn <- make_time_varying(
  times_tv, clip01(tv_matrix$ipc_helper)
)

# ── Fixed efficacy scalars ────────────────────────────────────────────────
etu_efficacy_baseline     <- 0.90   # fiber DEFAULT_SCALAR_INPUTS$etu_efficacy
non_etu_hospital_efficacy <- 0.30   # fiber DEFAULT_SCALAR_INPUTS$general_hospital_quarantine_efficacy
ppe_efficacy_hcw          <- 0.70   # fiber DEFAULT_SCALAR_INPUTS$ppe_efficacy

# ==============================================================================
# [Section 4c] Diagnostic plot — time-varying response parameters
# ==============================================================================

t_seq <- seq(0, 365, by = 1)

tv_diag <- data.frame(
  day                   = t_seq,
  prob_hosp_genPop      = prob_hospitalised_genPop_fn(t_seq),
  hosp_delay_factor     = hospitalisation_delay_factor(t_seq),
  p_unsafe_funeral_comm = p_unsafe_funeral_comm_fn(t_seq),
  p_unsafe_funeral_hosp = p_unsafe_funeral_hosp_fn(t_seq),
  prop_etu              = prop_etu_fn(t_seq),
  ipc_index             = ipc_index_fn(t_seq),
  ppe_coverage          = ppe_coverage_fn(t_seq)
)

# Derived efficacy curves for reference
tv_diag$etu_eff       <- etu_efficacy_baseline +
  (1 - etu_efficacy_baseline) * tv_diag$ipc_index
tv_diag$hosp_quar_eff <- tv_diag$prop_etu * tv_diag$etu_eff +
  (1 - tv_diag$prop_etu) * tv_diag$ipc_index
tv_diag$ppe_eff       <- tv_diag$ppe_coverage * ppe_efficacy_hcw

tv_long <- tidyr::pivot_longer(
  tv_diag, cols = -day, names_to = "parameter", values_to = "value"
)

param_labels <- c(
  prob_hosp_genPop      = "P(hospitalised | genPop)",
  hosp_delay_factor     = "Hospitalisation delay factor",
  p_unsafe_funeral_comm = "P(unsafe funeral | comm death)",
  p_unsafe_funeral_hosp = "P(unsafe funeral | hosp death)",
  prop_etu              = "Prop. cases in ETU",
  ipc_index             = "IPC maturity index",
  ppe_coverage          = "PPE coverage",
  etu_eff               = "ETU efficacy [derived]",
  hosp_quar_eff         = "Hospital quarantine eff. [derived]",
  ppe_eff               = "Effective PPE eff. [derived]"
)
tv_long$param_label <- param_labels[tv_long$parameter]
tv_long$group <- ifelse(
  tv_long$parameter %in% c("etu_eff", "hosp_quar_eff", "ppe_eff"),
  "Derived", "Input"
)

p_tv <- ggplot(tv_long, aes(x = day, y = value, color = group)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ param_label, scales = "free_y", ncol = 3) +
  scale_color_manual(values = c("Input" = "#185FA5", "Derived" = "#993C1D"),
                     name = NULL) +
  labs(title    = sprintf("%s — Time-varying response parameters (DRC-like)", case_tag),
       subtitle = "Input parameters (blue) | derived hospital/PPE efficacy (red)",
       x = "Day since outbreak start", y = "Value") +
  theme_bw() +
  theme(plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "grey40"),
        strip.text    = element_text(size = 8,  face = "bold"),
        legend.position = "top")

out_tv <- file.path(fig_dir, sprintf("%s_tv_params.png", case_tag))
ggsave(out_tv, plot = p_tv, width = 14, height = 10, dpi = 150)
message(sprintf("Saved: %s", out_tv))

# ==============================================================================
# [Section 4] Run simulation
# ==============================================================================

message("\n=== Running simulation ===")
t_start <- proc.time()[["elapsed"]]

result <- ebola_network_sim(
  sim_prep  = sim_prep,
  nodes     = nodes,

  seeding_cases     = 25L,
  seeding_ids       = NULL,

  incubation_period_fn           = incubation_period_fn,
  onset_to_hospitalisation_fn    = onset_to_hospitalisation_fn,
  onset_to_death_fn              = onset_to_death_fn,
  onset_to_recovery_fn           = onset_to_recovery_fn,
  hospitalisation_to_death_fn    = hospitalisation_to_death_fn,
  hospitalisation_to_recovery_fn = hospitalisation_to_recovery_fn,
  generation_time_fn             = generation_time_fn,
  funeral_generation_time_fn     = funeral_generation_time_fn,

  prob_hospitalised_genPop  = prob_hospitalised_genPop,
  prob_hospitalised_hcw     = prob_hospitalised_hcw,
  prob_death_comm           = prob_death_comm,
  prob_death_hosp           = prob_death_hosp,

  p_inf_household_close            = p_inf_household_close,
  p_inf_household_physical         = p_inf_household_physical,
  close_only_home                  = close_only_home,
  phys_only_home                   = phys_only_home,
  p_inf_community_close_daily      = p_inf_community_close_daily,
  p_inf_community_close_weekly     = p_inf_community_close_weekly,
  p_inf_community_close_monthly    = p_inf_community_close_monthly,
  p_inf_community_physical_daily   = p_inf_community_physical_daily,
  p_inf_community_physical_weekly  = p_inf_community_physical_weekly,
  p_inf_community_physical_monthly = p_inf_community_physical_monthly,
  p_inf_hcw_to_hcw                 = p_inf_hcw_to_hcw,
  p_inf_patient_to_hcw             = p_inf_patient_to_hcw,

  funeral_avg               = funeral_avg,
  funeral_k                 = funeral_k,
  p_unsafe_funeral          = p_unsafe_funeral,
  p_inf_funeral_household   = p_inf_funeral_household,
  p_inf_funeral_community   = p_inf_funeral_community,
  funeral_unsafe_multiplier = funeral_unsafe_multiplier,
  funeral_safe_multiplier   = funeral_safe_multiplier,

  # Time-varying response parameters
  prob_hospitalised_genPop_fn   = prob_hospitalised_genPop_fn,
  hospitalisation_delay_factor  = hospitalisation_delay_factor,
  p_unsafe_funeral_comm_fn      = p_unsafe_funeral_comm_fn,
  p_unsafe_funeral_hosp_fn      = p_unsafe_funeral_hosp_fn,
  prop_etu_fn                   = prop_etu_fn,
  ipc_index_fn                  = ipc_index_fn,
  etu_efficacy_baseline         = etu_efficacy_baseline,
  non_etu_hospital_efficacy     = non_etu_hospital_efficacy,
  ppe_coverage_fn               = ppe_coverage_fn,
  ppe_efficacy_hcw              = ppe_efficacy_hcw,

  max_time           = max_time,
  max_infected       = max_infected,
  seed               = sim_seed,
  monitoring_console = TRUE
)

elapsed <- round(proc.time()[["elapsed"]] - t_start, 1)
message(sprintf("\n=== Done: %.1f sec ===", elapsed))

# ==============================================================================
# [Section 5] Summary
# ==============================================================================

inf_df <- result$infected

ctype_labels <- c("0"="Index", "1"="Household", "2"="Comm close",
                  "3"="Comm physical", "4"="Hospital",
                  "5"="Funeral HH", "6"="Funeral comm")

cat(sprintf("\n── Simulation summary ──\n"))
cat(sprintf("  Case          : %s\n", case_tag))
cat(sprintf("  Stop reason   : %s\n", result$stop_reason))
cat(sprintf("  Total infected: %d\n", result$n_cumul_infected))
cat(sprintf("  Deaths        : %d (%.1f%%)\n",
            sum(inf_df$outcome_death, na.rm=TRUE),
            100*mean(inf_df$outcome_death, na.rm=TRUE)))
cat(sprintf("  Hospitalised  : %d (%.1f%%)\n",
            sum(inf_df$hospitalised, na.rm=TRUE),
            100*mean(inf_df$hospitalised, na.rm=TRUE)))
cat(sprintf("  HCWs infected : %d\n", sum(inf_df$is_hcw, na.rm=TRUE)))
cat(sprintf("  Generations   : %d\n", max(inf_df$generation, na.rm=TRUE)))
cat(sprintf("  Duration      : %.0f days\n",
            diff(range(inf_df$time_infection, na.rm=TRUE))))
cat(sprintf("  Funeral attendees infected: %d\n",
            sum(!is.na(inf_df$funeral_attended_for))))

ct <- table(inf_df$contact_type)
cat("\n  Transmission by route:\n")
for (i in names(ct))
  cat(sprintf("    %-16s: %d (%.1f%%)\n",
              ctype_labels[i], ct[i], 100*ct[i]/sum(ct)))

# ==============================================================================
# [Section 6] Plots
# ==============================================================================

inf_df$day_inf     <- floor(inf_df$time_infection)
inf_df$ctype_label <- ctype_labels[as.character(inf_df$contact_type)]

epi_curve <- inf_df %>%
  group_by(day_inf) %>%
  summarise(n_new = n(), .groups = "drop")

p_epi <- ggplot(epi_curve, aes(x = day_inf, y = n_new)) +
  geom_col(fill = "tomato", alpha = 0.8) +
  labs(title    = sprintf("%s — Epidemic curve", case_tag),
       subtitle = sprintf("Total: %d | Stop: %s",
                          result$n_cumul_infected, result$stop_reason),
       x = "Day", y = "New infections") +
  theme_bw() +
  theme(plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "grey40"))

p_route <- ggplot(inf_df %>% filter(!is.na(day_inf), contact_type > 0),
                  aes(x = day_inf, fill = ctype_label)) +
  geom_bar(alpha = 0.8) +
  scale_fill_brewer(palette = "Set2", name = "Route") +
  labs(title = "Transmission route over time", x = "Day", y = "Count") +
  theme_bw() + theme(legend.position = "top")

p_gen <- ggplot(inf_df %>% filter(!is.na(generation)),
                aes(x = factor(generation))) +
  geom_bar(fill = "steelblue", alpha = 0.8) +
  labs(title = "By generation", x = "Generation", y = "Count") +
  theme_bw()

p_funeral <- ggplot(inf_df %>% filter(!is.na(funeral_role)),
                    aes(x = funeral_role, fill = funeral_role)) +
  geom_bar(alpha = 0.8) +
  scale_fill_manual(values = c("household"="#993C1D","community"="#185FA5"),
                    guide = "none") +
  labs(title = "Funeral attendees by role", x = NULL, y = "Count") +
  theme_bw()

p_combined <- (p_epi + p_route) / (p_gen + p_funeral)

out_fig <- file.path(fig_dir,
                     sprintf("%s_seed%s_epi.png", case_tag,
                             ifelse(is.null(index_case_id),"rand",index_case_id)))
ggsave(out_fig, plot = p_combined, width = 12, height = 10, dpi = 150)
cat(sprintf("\nSaved: %s\n", out_fig))

# ==============================================================================
# [Section 7] Weekly incident deaths
# ==============================================================================

death_df <- inf_df %>%
  filter(isTRUE(outcome_death) | outcome_death == TRUE) %>%
  filter(!is.na(time_outcome)) %>%
  mutate(
    week      = floor(time_outcome / 7),
    is_hcw    = isTRUE(is_hcw) | is_hcw == TRUE,
    death_loc = outcome_location
  )

weekly_deaths <- death_df %>%
  group_by(week) %>%
  summarise(
    n_total   = n(),
    n_hcw     = sum(is_hcw, na.rm = TRUE),
    n_genpop  = n_total - n_hcw,
    n_comm    = sum(death_loc == "community", na.rm = TRUE),
    n_hosp    = sum(death_loc == "hospital",  na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(week_start = week * 7)

# Stacked bar: genPop community / genPop hospital / HCW
weekly_long <- weekly_deaths %>%
  select(week_start, n_comm, n_hosp, n_hcw) %>%
  tidyr::pivot_longer(
    cols      = c(n_comm, n_hosp, n_hcw),
    names_to  = "group",
    values_to = "deaths"
  ) %>%
  mutate(group = factor(group,
                        levels = c("n_hcw", "n_hosp", "n_comm"),
                        labels = c("HCW", "Hospital (genPop)", "Community (genPop)")))

p_weekly_death <- ggplot(weekly_long,
                         aes(x = week_start, y = deaths, fill = group)) +
  geom_col(alpha = 0.85, width = 6) +
  scale_fill_manual(
    values = c("HCW"                = "#993C1D",
               "Hospital (genPop)"  = "#E07B39",
               "Community (genPop)" = "#185FA5"),
    name = NULL
  ) +
  labs(
    title    = sprintf("%s — Weekly incident deaths", case_tag),
    subtitle = sprintf("Total deaths: %d | HCW: %d (%.1f%%)",
                       sum(weekly_deaths$n_total),
                       sum(weekly_deaths$n_hcw),
                       100 * sum(weekly_deaths$n_hcw) / max(sum(weekly_deaths$n_total), 1)),
    x = "Day (week start)", y = "Deaths per week"
  ) +
  theme_bw() +
  theme(
    plot.title      = element_text(size = 12, face = "bold"),
    plot.subtitle   = element_text(size = 9,  color = "grey40"),
    legend.position = "top"
  )

out_death <- file.path(fig_dir,
                       sprintf("%s_seed%s_weekly_deaths.png", case_tag,
                               ifelse(is.null(index_case_id), "rand", index_case_id)))
ggsave(out_death, plot = p_weekly_death, width = 10, height = 5, dpi = 150)
cat(sprintf("Saved: %s\n", out_death))

invisible(result)
