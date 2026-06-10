# ==============================================================================
# C1_ebola_sim_run.R
# Purpose:
#   Single run of Ebola IBM simulation for one of the 3 cases
#   Uses pre-built network (p6) and sim_prep (p7)
#
# Change [Configuration] section to select case and index case
# ==============================================================================

source("function/COD_IBM_ebola_sim.R")
library(dplyr)
library(ggplot2)
library(patchwork)
# ==============================================================================
# [Configuration] — Change ONLY these
# ==============================================================================

# case_tag      <- "case1_1M"   # "case1_1M", "case2_Ituri", "case3_Kivu"
# index_case_id <- 500000L + 42L           # person_id of index case (NULL = random)
case_tag      <- "case2_Ituri"   # "case1_1M", "case2_Ituri", "case3_Kivu"
index_case_id <- 3200000L + 42L           # person_id of index case (NULL = random)
sim_seed      <- 42L
max_time      <- 150           # max simulation days
max_infected  <- Inf           # stop if cumulative infections exceed this

network_dir  <- "output/network"
sim_prep_dir <- file.path(network_dir, "sim_prep")
fig_dir      <- "figure/C1_sim"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# [Section 1] Load network data
# ==============================================================================

message(sprintf("=== Loading: %s ===", case_tag))

nodes    <- readRDS(file.path(network_dir, sprintf("%s_nodes.rds", case_tag)))
sim_prep <- readRDS(file.path(sim_prep_dir, sprintf("%s_sim_prep.rds", case_tag)))

message(sprintf("  N = %d | HCWs = %d", nrow(nodes), sum(nodes$is_hcw)))

# Validate index case
if (!is.null(index_case_id)) {
  if (!index_case_id %in% nodes$person_id)
    stop(sprintf("index_case_id %d not found in nodes", index_case_id))
  idx_info <- nodes[nodes$person_id == index_case_id, ]
  message(sprintf("  Index case: person %d | age=%d | HCW=%s | cell=%d",
                  index_case_id, idx_info$age,
                  ifelse(idx_info$is_hcw, "YES", "no"),
                  idx_info$cell_id))
}

# ==============================================================================
# [Section 2] Natural history parameters — from fiber DEFAULT_SCALAR_INPUTS
# ==============================================================================

# Gamma distribution helper
make_gamma_sampler <- function(mean, sd) {
  shape <- (mean / sd)^2
  rate  <- mean / (sd^2)
  function(n) if (n <= 0L) numeric(0) else rgamma(n=n, shape=shape, rate=rate)
}

incubation_period_fn           <- make_gamma_sampler(mean=8.5,  sd=4.5)
onset_to_hospitalisation_fn    <- make_gamma_sampler(mean=1.0,  sd=0.35)
onset_to_death_fn              <- make_gamma_sampler(mean=9.3,  sd=3.0)
onset_to_recovery_fn           <- make_gamma_sampler(mean=13.0, sd=4.0)
hospitalisation_to_death_fn    <- make_gamma_sampler(mean=4.5,  sd=2.0)
hospitalisation_to_recovery_fn <- make_gamma_sampler(mean=8.0,  sd=2.5)
generation_time_fn             <- make_gamma_sampler(mean=15.4, sd=sqrt(15.4 / 2.5))
# generation time: shape=2.5, rate=shape/mean → sd=sqrt(mean/shape)

# ==============================================================================
# [Section 3] Transmission parameters — from fiber DEFAULT_SCALAR_INPUTS
# ==============================================================================

prob_symptomatic          <- 1.0   # Ebola: always symptomatic
prob_hospitalised_genPop  <- 0.50  # baseline (time-varying in fiber scenarios)
prob_hospitalised_hcw     <- 0.50
prob_death_comm           <- 0.70
prob_death_hosp           <- 0.50

# Per-contact transmission probabilities
# fiber uses offspring counts (mn_offspring_genPop=1.25); here we map to
# per-contact probabilities consistent with network structure
p_inf_household       <- 0.20
p_inf_community       <- 0.02
p_inf_hcw_to_hcw      <- 0.02
# p_inf_patient_to_hcw  <- 0.02
p_inf_patient_to_hcw  <- 0.02

# Funeral
p_inf_funeral_unsafe  <- 0.02
p_unsafe_funeral      <- 0.70   # baseline prob funeral is unsafe
p_inf_funeral_safe    <- p_inf_funeral_unsafe * (1 - 0.80)  # safe_funeral_efficacy=0.80

ppe_efficacy_hcw              <- 0.70
prob_hospital_cond_hcw_preAdm <- 0.50

# ==============================================================================
# [Section 4] Run simulation
# ==============================================================================

message("\n=== Running simulation ===")
t_start <- proc.time()[["elapsed"]]

result <- ebola_network_sim(
  # Network
  sim_prep          = sim_prep,
  nodes             = nodes,

  # Seeding
  seeding_cases     = 1L,
  seeding_ids       = index_case_id,

  # Natural history
  incubation_period_fn           = incubation_period_fn,
  onset_to_hospitalisation_fn    = onset_to_hospitalisation_fn,
  onset_to_death_fn              = onset_to_death_fn,
  onset_to_recovery_fn           = onset_to_recovery_fn,
  hospitalisation_to_death_fn    = hospitalisation_to_death_fn,
  hospitalisation_to_recovery_fn = hospitalisation_to_recovery_fn,
  generation_time_fn             = generation_time_fn,

  # Severity
  prob_hospitalised_genPop = prob_hospitalised_genPop,
  prob_hospitalised_hcw    = prob_hospitalised_hcw,
  prob_death_comm          = prob_death_comm,
  prob_death_hosp          = prob_death_hosp,

  # Transmission
  p_inf_household      = p_inf_household,
  p_inf_community      = p_inf_community,
  p_inf_hcw_to_hcw     = p_inf_hcw_to_hcw,
  p_inf_patient_to_hcw = p_inf_patient_to_hcw,
  p_unsafe_funeral     = p_unsafe_funeral,
  p_inf_funeral_unsafe = p_inf_funeral_unsafe,
  p_inf_funeral_safe   = p_inf_funeral_safe,

  # Stopping
  max_time      = max_time,
  max_infected  = max_infected,

  # Misc
  seed               = sim_seed,
  monitoring_console = TRUE
)

elapsed <- round(proc.time()[["elapsed"]] - t_start, 1)
message(sprintf("\n=== Done: %.1f sec ===", elapsed))

# ==============================================================================
# [Section 5] Summary
# ==============================================================================

inf_df <- result$infected

cat(sprintf("\n── Simulation summary ──\n"))
cat(sprintf("  Case          : %s\n", case_tag))
cat(sprintf("  Index case    : person %s\n",
            ifelse(is.null(index_case_id), "random", index_case_id)))
cat(sprintf("  Stop reason   : %s\n", result$stop_reason))
cat(sprintf("  Total infected: %d\n", result$n_cumul_infected))
cat(sprintf("  Deaths        : %d (%.0f%%)\n",
            sum(inf_df$outcome_death, na.rm=TRUE),
            100*mean(inf_df$outcome_death, na.rm=TRUE)))
cat(sprintf("  Hospitalised  : %d (%.0f%%)\n",
            sum(inf_df$hospitalised, na.rm=TRUE),
            100*mean(inf_df$hospitalised, na.rm=TRUE)))
cat(sprintf("  HCWs infected : %d\n", sum(inf_df$is_hcw, na.rm=TRUE)))
cat(sprintf("  Generations   : %d\n", max(inf_df$generation, na.rm=TRUE)))
cat(sprintf("  Duration      : %.0f days\n",
            diff(range(inf_df$time_infection, na.rm=TRUE))))

# Contact type breakdown
ctype_labels <- c("0"="Index","1"="Household","2"="Community",
                  "3"="Hospital","4"="Funeral")
ct <- table(inf_df$contact_type)
cat("\n  Transmission by route:\n")
for (i in names(ct))
  cat(sprintf("    %-12s: %d (%.0f%%)\n",
              ctype_labels[i], ct[i], 100*ct[i]/sum(ct)))

# ==============================================================================
# [Section 6] Epidemic curve
# ==============================================================================

# Daily incidence
inf_df$day_inf <- floor(inf_df$time_infection)
epi_curve <- inf_df %>%
  group_by(day_inf) %>%
  summarise(n_new = n(), .groups="drop") %>%
  arrange(day_inf)

p_epi <- ggplot(epi_curve, aes(x=day_inf, y=n_new)) +
  geom_col(fill="tomato", alpha=0.8) +
  labs(title    = sprintf("%s — Epidemic curve (index: person %s)",
                          case_tag,
                          ifelse(is.null(index_case_id),"random",index_case_id)),
       subtitle = sprintf("Total: %d infected | Stop: %s",
                          result$n_cumul_infected, result$stop_reason),
       x="Day", y="New infections") +
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"),
        plot.subtitle=element_text(size=9, color="grey40"))

# Generation distribution
p_gen <- ggplot(inf_df %>% filter(!is.na(generation)),
                aes(x=factor(generation))) +
  geom_bar(fill="steelblue", alpha=0.8) +
  labs(title="Infections by generation",
       x="Generation", y="Count") +
  theme_bw()

# Contact type over time
inf_df$ctype_label <- ctype_labels[as.character(inf_df$contact_type)]
p_route <- ggplot(inf_df %>% filter(!is.na(day_inf), contact_type > 0),
                  aes(x=day_inf, fill=ctype_label)) +
  geom_bar(alpha=0.8) +
  scale_fill_brewer(palette="Set2", name="Route") +
  labs(title="Transmission route over time",
       x="Day", y="New infections") +
  theme_bw() +
  theme(legend.position="top")

# Save
out_fig <- file.path(fig_dir,
                     sprintf("%s_seed%s_epi.png", case_tag,
                             ifelse(is.null(index_case_id),"rand",index_case_id)))
p_combined <- (p_epi / (p_gen + p_route))
ggsave(out_fig, plot=p_combined, width=12, height=10, dpi=150)
cat(sprintf("\nSaved: %s\n", out_fig))

# Return result for further analysis
invisible(result)

dir.create("output/sim", showWarnings = FALSE, recursive = TRUE)
saveRDS(result, sprintf("output/sim/%s_sim_result.rds", case_tag))
