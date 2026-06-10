# ==============================================================================
# COD_IBM_ebola_test.R
# Purpose:
#   Build a minimal synthetic network and run the Ebola IBM simulation
#   to verify the model runs correctly end-to-end.
#
# Synthetic population:
#   500 individuals, 120 households, 3 hospitals, ~15 HCWs
#   Layer 1: household cliques
#   Layer 2: random community edges
#   Layer 3: HCW-HCW cliques + admission lookup
# ==============================================================================

library(dplyr)

source("function/COD_IBM_ebola_sim.R")

set.seed(42)

# ==============================================================================
# [1] Synthetic population
# ==============================================================================
n_people    <- 500L
n_hospitals <- 3L
n_hh        <- 120L

# Assign individuals to households (random sizes 2-6)
hh_ids <- sample(seq_len(n_hh), n_people, replace = TRUE)

# Age groups (Prem 16 groups)
age_groups <- sample(1:16, n_people, replace = TRUE,
                     prob = c(0.18, 0.16, 0.15, 0.12, 0.10,
                              0.08, 0.06, 0.05, 0.04, 0.02,
                              0.01, 0.01, 0.01, 0.005, 0.005, 0.005))
ages <- age_groups * 5 - 2  # Approximate midpoint age

# HCWs: ~3% of adults
is_adult  <- ages >= 18
is_hcw    <- logical(n_people)
adult_idx <- which(is_adult)
n_hcw     <- round(n_people * 0.03)
hcw_idx   <- sample(adult_idx, n_hcw)
is_hcw[hcw_idx] <- TRUE

# Assign hospital (nearest = random for synthetic data)
hospital_ids <- sample(seq_len(n_hospitals), n_people, replace = TRUE)

# Cell IDs (simplified: one cell per 10 people)
cell_ids <- ceiling(seq_len(n_people) / 10)

nodes <- data.frame(
  person_id   = seq_len(n_people),
  hh_id       = hh_ids,
  cell_id     = cell_ids,
  hospital_id = hospital_ids,
  age_group   = age_groups,
  is_hcw      = is_hcw,
  is_adult    = is_adult,
  stringsAsFactors = FALSE
)

cat(sprintf("Population  : %d individuals\n", n_people))
cat(sprintf("Households  : %d\n", n_hh))
cat(sprintf("HCWs        : %d (%.1f%%)\n", sum(is_hcw), 100*mean(is_hcw)))
cat(sprintf("Hospitals   : %d\n", n_hospitals))

# ==============================================================================
# [2] Layer 1 — Household edges (full clique within each household)
# ==============================================================================
hh_edges <- lapply(split(nodes$person_id, nodes$hh_id), function(members) {
  if (length(members) < 2) return(NULL)
  pairs <- combn(members, 2)
  data.frame(from = pairs[1,], to = pairs[2,], stringsAsFactors = FALSE)
}) %>% Filter(Negate(is.null), .) %>% bind_rows()

cat(sprintf("Layer 1 HH edges   : %d\n", nrow(hh_edges)))

# ==============================================================================
# [3] Layer 2 — Community edges (random, ~8 per person)
# ==============================================================================
n_comm_edges <- n_people * 4L  # ~8 per person (both directions → /2)
from_comm    <- sample(seq_len(n_people), n_comm_edges, replace = TRUE)
to_comm      <- sample(seq_len(n_people), n_comm_edges, replace = TRUE)

# Remove self-loops and same-household edges, keep i < j
same_hh <- nodes$hh_id[from_comm] == nodes$hh_id[to_comm]
comm_df <- data.frame(from = from_comm, to = to_comm) %>%
  filter(from != to, !same_hh, from < to) %>%
  distinct(from, to) %>%
  mutate(from_cell = nodes$cell_id[from],
         to_cell   = nodes$cell_id[to])

cat(sprintf("Layer 2 comm edges : %d\n", nrow(comm_df)))

# ==============================================================================
# [4] Layer 3 — Healthcare edges
# ==============================================================================

# (a) HCW-HCW edges within same hospital
hcw_nodes     <- nodes[nodes$is_hcw, ]
hcw_hcw_edges <- hcw_nodes %>%
  group_by(hospital_id) %>%
  group_map(~ {
    members <- .x$person_id
    if (length(members) < 2) return(NULL)
    pairs <- combn(members, 2)
    data.frame(from        = pairs[1,],
               to          = pairs[2,],
               hospital_id = .y$hospital_id,
               stringsAsFactors = FALSE)
  }) %>% bind_rows()

cat(sprintf("Layer 3 HCW edges  : %d\n", nrow(hcw_hcw_edges)))

# (b) Admission lookup: non-HCW → hospital + HCW list
hospital_hcw_lookup <- hcw_nodes %>%
  group_by(hospital_id) %>%
  summarise(hcw_list = list(person_id), .groups = "drop")

non_hcw_nodes     <- nodes[!nodes$is_hcw, ]
admission_lookup  <- non_hcw_nodes %>%
  select(person_id, hh_id, hospital_id) %>%
  left_join(hospital_hcw_lookup, by = "hospital_id")

cat(sprintf("Layer 3 admission  : %d individuals\n", nrow(admission_lookup)))

# ==============================================================================
# [5] Ebola natural history parameter functions
# ==============================================================================

# Incubation period: Gamma(shape=1.5, rate=0.15) → mean ~10 days
incubation_period_fn <- function(n)
  rgamma(n, shape = 1.5, rate = 0.15)

# Onset to hospitalisation: Gamma(shape=2, rate=0.4) → mean ~5 days
onset_to_hospitalisation_fn <- function(n)
  rgamma(n, shape = 2, rate = 0.4)

# Onset to death (community): Gamma(shape=2, rate=0.18) → mean ~11 days
onset_to_death_fn <- function(n)
  rgamma(n, shape = 2, rate = 0.18)

# Onset to recovery (community): Gamma(shape=2, rate=0.13) → mean ~15 days
onset_to_recovery_fn <- function(n)
  rgamma(n, shape = 2, rate = 0.13)

# Hospitalisation to death: Gamma(shape=2, rate=0.25) → mean ~8 days
hospitalisation_to_death_fn <- function(n)
  rgamma(n, shape = 2, rate = 0.25)

# Hospitalisation to recovery: Gamma(shape=2, rate=0.17) → mean ~12 days
hospitalisation_to_recovery_fn <- function(n)
  rgamma(n, shape = 2, rate = 0.17)

# Generation time: Gamma(shape=2, rate=0.22) → mean ~9 days
generation_time_fn <- function(n)
  rgamma(n, shape = 2, rate = 0.22)

# ==============================================================================
# [6] Run simulation
# ==============================================================================
cat("\n=== Running Ebola IBM simulation ===\n")
t0_sim <- proc.time()

result <- ebola_network_sim(

  # Network
  nodes              = nodes,
  layer1_hh          = hh_edges,
  layer2_comm        = comm_df,
  layer3_hcw_edges   = hcw_hcw_edges,
  layer3_admission   = admission_lookup,
  cell_dist          = NULL,  # Not used in transmission logic directly

  # Seeding
  seeding_cases      = 1L,
  t0                 = 0,

  # Natural history
  incubation_period_fn           = incubation_period_fn,
  onset_to_hospitalisation_fn    = onset_to_hospitalisation_fn,
  onset_to_death_fn              = onset_to_death_fn,
  onset_to_recovery_fn           = onset_to_recovery_fn,
  hospitalisation_to_death_fn    = hospitalisation_to_death_fn,
  hospitalisation_to_recovery_fn = hospitalisation_to_recovery_fn,
  generation_time_fn             = generation_time_fn,

  # Disease severity
  prob_symptomatic          = 1.0,
  prob_hospitalised_genPop  = 0.5,
  prob_hospitalised_hcw     = 0.8,   # HCWs more likely to be hospitalised
  prob_death_comm           = 0.5,
  prob_death_hosp           = 0.3,   # Better outcomes in hospital

  # Transmission probabilities
  p_inf_household           = 0.15,
  p_inf_community           = 0.03,
  p_inf_hcw_to_hcw          = 0.10,
  p_inf_patient_to_hcw      = 0.08,
  p_unsafe_funeral          = 0.7,
  p_inf_funeral_unsafe      = 0.20,
  p_inf_funeral_safe        = 0.01,

  # HCW-specific
  ppe_efficacy_hcw              = 0.5,
  prob_hospital_cond_hcw_preAdm = 0.6,

  # No interventions for baseline test
  antiviral_start    = Inf,
  quarantine_start   = Inf,
  vax_start          = Inf,
  prob_vax           = 0,
  prob_treat_self    = 0,
  prob_treat_household   = 0,
  prob_treat_community   = 0,
  quarantine_efficacy    = 0,
  prob_quarantine_self       = 0,
  prob_quarantine_household  = 0,
  prob_quarantine_community  = 0,
  logistical_delay_fn    = function(n) rep(0, n),
  time_to_quarantine_fn  = function(n) rep(0, n),

  # Stopping
  max_infected       = 5000L,
  max_time           = 365,

  seed               = 42,
  monitoring_console = TRUE
)

elapsed <- round((proc.time() - t0_sim)[["elapsed"]], 1)

# ==============================================================================
# [7] Summary output
# ==============================================================================
inf  <- result$infected
full <- result$full

cat(sprintf("\n=== Simulation Results ===\n"))
cat(sprintf("  Stop reason        : %s\n", result$stop_reason))
cat(sprintf("  Total infected     : %d / %d (%.1f%%)\n",
            result$n_cumul_infected, n_people,
            100 * result$n_cumul_infected / n_people))
cat(sprintf("  Elapsed time       : %.1f sec\n", elapsed))

if (result$n_cumul_infected > 0) {
  cat(sprintf("\n  --- Case breakdown ---\n"))
  cat(sprintf("  HCW infected       : %d / %d (%.1f%%)\n",
              sum(inf$is_hcw), sum(nodes$is_hcw),
              100 * sum(inf$is_hcw) / max(sum(nodes$is_hcw), 1)))
  cat(sprintf("  genPop infected    : %d\n", sum(!inf$is_hcw)))

  cat(sprintf("\n  --- Transmission routes ---\n"))
  ct <- table(inf$contact_type)
  route_labels <- c("0" = "seed", "1" = "household", "2" = "community",
                    "3" = "hospital", "4" = "funeral")
  for (nm in names(ct))
    cat(sprintf("  %-12s : %d\n", route_labels[nm], ct[[nm]]))

  cat(sprintf("\n  --- Outcomes ---\n"))
  cat(sprintf("  Deaths             : %d (CFR %.1f%%)\n",
              sum(inf$outcome_death, na.rm = TRUE),
              100 * mean(inf$outcome_death, na.rm = TRUE)))
  cat(sprintf("  Hospitalised       : %d (%.1f%%)\n",
              sum(inf$hospitalised, na.rm = TRUE),
              100 * mean(inf$hospitalised, na.rm = TRUE)))
  cat(sprintf("  Community deaths   : %d\n",
              sum(inf$outcome_death & inf$outcome_location == "community",
                  na.rm = TRUE)))
  cat(sprintf("  Hospital deaths    : %d\n",
              sum(inf$outcome_death & inf$outcome_location == "hospital",
                  na.rm = TRUE)))
  cat(sprintf("  Unsafe funerals    : %d\n",
              sum(inf$funeral_unsafe, na.rm = TRUE)))

  cat(sprintf("\n  --- Epidemic timing ---\n"))
  cat(sprintf("  First infection    : day %.1f\n",
              min(inf$time_infection, na.rm = TRUE)))
  cat(sprintf("  Last outcome       : day %.1f\n",
              max(inf$time_outcome, na.rm = TRUE)))
  cat(sprintf("  Peak active        : %d cases\n",
              max(tabulate(floor(inf$time_infection[!is.na(inf$time_infection)])
                           + 1L), na.rm = TRUE)))
}
