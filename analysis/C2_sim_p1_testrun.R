# ==============================================================================
# C2_sim_p1_testrun.R
# Purpose:
#   Single test run of IBM simulation
#   Requires sim_prep to be built first via C2_sim_p0_simprep.R
# ==============================================================================

source("function/COD_IBM_ebola_sim.R")
library(dplyr)
library(ggplot2)
library(patchwork)

# ==============================================================================
# [Configuration]
# ==============================================================================

case_tag      <- "case1_1M"
index_case_id <- 42L          # NULL = random
sim_seed      <- 42L
max_time      <- 365
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

incubation_period_fn           <- make_gamma_fn(9,  4)
onset_to_hospitalisation_fn    <- make_gamma_fn(5,  2)
onset_to_death_fn              <- make_gamma_fn(12, 4)
onset_to_recovery_fn           <- make_gamma_fn(18, 5)
hospitalisation_to_death_fn    <- make_gamma_fn(8,  3)
hospitalisation_to_recovery_fn <- make_gamma_fn(13, 4)
generation_time_fn             <- make_gamma_fn(12, 4)

# ==============================================================================
# [Section 3] Transmission parameters
# ==============================================================================

p_inf_household_close            <- 0.05
p_inf_household_physical         <- 0.15

# Per single contact-event probability (3-week effective computed internally)
p_inf_community_close_daily      <- 0.00
p_inf_community_close_weekly     <- 0.00
p_inf_community_close_monthly    <- 0.00

p_inf_community_physical_daily   <- 0.05
p_inf_community_physical_weekly  <- 0.08
p_inf_community_physical_monthly <- 0.10

p_inf_hcw_to_hcw                 <- 0.03
p_inf_patient_to_hcw             <- 0.02

funeral_avg                      <- 20
funeral_k                        <- 2
p_unsafe_funeral                 <- 0.50
p_inf_funeral_household          <- 0.15
p_inf_funeral_community          <- 0.05
funeral_unsafe_multiplier        <- 1.0
funeral_safe_multiplier          <- 0.1

prob_hospitalised_genPop         <- 0.30
prob_hospitalised_hcw            <- 0.60
prob_death_comm                  <- 0.05
prob_death_hosp                  <- 0.02

# ==============================================================================
# [Section 4] Run simulation
# ==============================================================================

message("\n=== Running simulation ===")
t_start <- proc.time()[["elapsed"]]

result <- ebola_network_sim(
  sim_prep  = sim_prep,
  nodes     = nodes,

  seeding_cases     = 1L,
  seeding_ids       = index_case_id,

  incubation_period_fn           = incubation_period_fn,
  onset_to_hospitalisation_fn    = onset_to_hospitalisation_fn,
  onset_to_death_fn              = onset_to_death_fn,
  onset_to_recovery_fn           = onset_to_recovery_fn,
  hospitalisation_to_death_fn    = hospitalisation_to_death_fn,
  hospitalisation_to_recovery_fn = hospitalisation_to_recovery_fn,
  generation_time_fn             = generation_time_fn,

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

invisible(result)
