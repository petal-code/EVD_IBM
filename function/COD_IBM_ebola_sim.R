# ==============================================================================
# COD_IBM_ebola_sim.R
# Purpose:
#   Individual-based Ebola transmission model for Kinshasa
#   Network structure: pre-built Kinshasa contact network (Layer 1/2/3)
#   Natural history: fiber-based Ebola disease progression
#   Transmission: network_bp_sim_boost-style gate mechanism
#
# Disease progression (fiber):
#   Infection → Incubation → Symptoms (prob_symptomatic = 1 for Ebola)
#     → Hospitalisation (prob varies by genPop/HCW)
#       → Community: death (prob_death_comm) or recovery
#       → Hospital:  death (prob_death_hosp) or recovery
#         → Death: unsafe/safe funeral (p_unsafe_funeral)
#
# Transmission phases:
#   Phase 1 (infection ~ hospitalisation):
#     Layer 1 (household): p_inf_household
#     Layer 2 (community): p_inf_community
#     HCW only: Layer 3 HCW-HCW: p_inf_hcw_to_hcw × (1 - ppe_efficacy_hcw)
#   Phase 2 (hospitalisation ~ outcome):
#     genPop: Layer 3 admission (patient→HCW): p_inf_patient_to_hcw
#     HCW:    Layer 3 HCW-HCW: p_inf_hcw_to_hcw
#   Phase 3 (death → funeral):
#     unsafe: Layer 1 + Layer 2: p_inf_funeral_unsafe
#     safe:   Layer 1 + Layer 2: p_inf_funeral_safe
#
# Each transmission attempt passes 4 gates:
#   Gate 1: transmission probability
#   Gate 2: source-side antiviral/vaccine efficacy
#   Gate 3: source quarantine
#   Gate 4: target-side antiviral/vaccine efficacy
# ==============================================================================

ebola_network_sim <- function(

  # ── Pre-built sim prep (from p7) ─────────────────────────────────────────────
  sim_prep,                       # output of C1_network_p7_sim_prep.R
  # list(N, person_ids, pid_to_idx,
  #      hh_nbrs, comm_nbrs, hcw_nbrs, adm_nbrs,
  #      hh_ring, ring1, ring2)
  nodes,                          # node dataframe from network build

  # ── Seeding ───────────────────────────────────────────────────────────────────
  seeding_cases    = 1L,          # number of index cases
  seeding_ids      = NULL,        # specific person_ids to seed (NULL = random)
  t0               = 0,           # simulation start time

  # ── Natural history (fiber-style, all functions of n) ────────────────────────
  incubation_period_fn,           # function(n): incubation period (days)
  onset_to_hospitalisation_fn,    # function(n): symptom onset to hospitalisation
  onset_to_death_fn,              # function(n): symptom onset to death (community)
  onset_to_recovery_fn,           # function(n): symptom onset to recovery (community)
  hospitalisation_to_death_fn,    # function(n): hospitalisation to death
  hospitalisation_to_recovery_fn, # function(n): hospitalisation to recovery
  generation_time_fn,             # function(n): generation time distribution

  # ── Disease severity ─────────────────────────────────────────────────────────
  prob_symptomatic          = 1.0,  # Ebola: always symptomatic
  prob_hospitalised_genPop,         # scalar or function(t): genPop hospitalisation prob
  prob_hospitalised_hcw,            # scalar or function(t): HCW hospitalisation prob
  prob_death_comm,                  # probability of death in community
  prob_death_hosp,                  # probability of death in hospital

  # ── Transmission probabilities ───────────────────────────────────────────────
  p_inf_household,                  # per-contact transmission prob, household
  p_inf_community,                  # per-contact transmission prob, community
  p_inf_hcw_to_hcw,                 # HCW → HCW transmission prob (hospital)
  p_inf_patient_to_hcw,             # patient → HCW transmission prob (hospital)
  p_unsafe_funeral,                 # scalar or function(t): prob funeral is unsafe
  p_inf_funeral_unsafe,             # transmission prob at unsafe funeral
  p_inf_funeral_safe,               # transmission prob at safe funeral

  # ── HCW-specific parameters ──────────────────────────────────────────────────
  ppe_efficacy_hcw              = 0,    # PPE efficacy reducing HCW pre-admission hospital transmission
  prob_hospital_cond_hcw_preAdm = 0.5, # prob HCW works in hospital while infectious pre-admission

  # ── Interventions ─────────────────────────────────────────────────────────────
  antiviral_start               = Inf,  # time antivirals become available
  quarantine_start              = Inf,  # time quarantine starts
  vax_start                     = Inf,  # time vaccination starts
  prob_vax                      = 0,    # probability of vaccination
  vax_eff_inf                   = 0,    # vaccine efficacy against infection
  vax_eff_trans                 = 0,    # vaccine efficacy against transmission
  vax_target_ages               = NULL, # age groups to target (NULL = all)
  prob_treat_self               = 0,
  prob_treat_household          = 0,
  prob_treat_community          = 0,
  quarantine_efficacy           = 0,
  prob_quarantine_self          = 0,
  prob_quarantine_household     = 0,
  prob_quarantine_community     = 0,
  logistical_delay_fn           = function(n) rep(0, n),
  time_to_quarantine_fn         = function(n) rep(0, n),
  pk_params                     = NULL,
  pd_params_inf                 = NULL,
  pd_params_trans               = NULL,

  # ── Stopping criteria ────────────────────────────────────────────────────────
  max_infected                  = Inf,
  max_time                      = Inf,

  # ── Misc ─────────────────────────────────────────────────────────────────────
  seed                          = 42,
  monitoring_console            = FALSE

) {

  set.seed(seed)

  # ── Helper: resolve scalar or function(t) ──────────────────────────────────
  resolve_tv <- function(param, t) {
    if (is.function(param)) param(t) else rep(param, length(t))
  }

  # ── Load pre-built neighbor lookups from sim_prep ─────────────────────────
  message("Loading pre-built neighbor lookup tables...")
  N              <- sim_prep$N
  pid_to_idx     <- sim_prep$pid_to_idx
  hh_nbrs        <- sim_prep$hh_nbrs
  comm_nbrs      <- sim_prep$comm_nbrs
  hcw_nbrs       <- sim_prep$hcw_nbrs
  admission_nbrs <- sim_prep$adm_nbrs
  message(sprintf("  Neighbor lookup loaded for %d individuals", N))

  # Ring helpers — computed on-demand from adjacency lists
  get_ring1 <- function(idx)
    unique(c(hh_nbrs[[idx]], comm_nbrs[[idx]]))

  get_ring2 <- function(idx) {
    r1 <- get_ring1(idx)
    if (length(r1) == 0) return(integer(0))
    unique(unlist(lapply(r1, get_ring1), use.names = FALSE))
  }

  # ── Initialize simulation state dataframe ──────────────────────────────────
  # status: 1L = susceptible, 2L = infected/active, 3L = recovered/dead
  tdf <- data.frame(
    person_id                      = nodes$person_id,
    hh_id                          = nodes$hh_id,
    cell_id                        = nodes$cell_id,
    hospital_id                    = nodes$hospital_id,
    age_group                      = nodes$age_group,
    is_hcw                         = nodes$is_hcw,
    status                         = rep(1L, N),
    # Natural history timings (absolute)
    time_infection                 = NA_real_,
    time_onset                     = NA_real_,
    time_hospitalisation           = NA_real_,
    time_outcome                   = NA_real_,
    # Natural history flags
    symptomatic                    = NA,
    hospitalised                   = NA,
    outcome_death                  = NA,        # TRUE = death, FALSE = recovery
    outcome_location               = NA_character_,  # "community" or "hospital"
    funeral_unsafe                 = NA,
    # Transmission tracking
    generation                     = NA_integer_,
    ancestor_id                    = NA_integer_,
    contact_type                   = NA_integer_,  # 1=HH, 2=comm, 3=hospital, 4=funeral
    # Intervention state
    treated                        = 0L,
    time_treated                   = NA_real_,
    quarantined                    = 0L,
    time_quarantined               = NA_real_,
    vaccinated                     = 0L,
    time_vaccinated                = NA_real_,
    stringsAsFactors               = FALSE
  )

  # ── Seed cases ─────────────────────────────────────────────────────────────
  seed_idx <- if (!is.null(seeding_ids)) {
    pid_to_idx[seeding_ids]
  } else {
    sample(seq_len(N), seeding_cases, replace = FALSE)
  }

  for (sid in seed_idx) {
    inc_p  <- incubation_period_fn(1)
    t_inf  <- t0
    t_onset <- t_inf + inc_p

    # Natural history for seed case
    p_hosp_t <- resolve_tv(
      if (tdf$is_hcw[sid]) prob_hospitalised_hcw else prob_hospitalised_genPop,
      t_onset
    )
    is_hosp   <- rbinom(1, 1, p_hosp_t) == 1L
    t_hosp    <- if (is_hosp) t_onset + onset_to_hospitalisation_fn(1) else NA_real_

    # Community outcome first, then check if hospitalisation precedes it
    is_death_comm <- rbinom(1, 1, prob_death_comm) == 1L
    t_comm_outcome <- if (is_death_comm) t_onset + onset_to_death_fn(1) else
      t_onset + onset_to_recovery_fn(1)

    if (is_hosp && !is.na(t_hosp) && t_hosp < t_comm_outcome) {
      # Successfully hospitalised
      second_chance_death <- (prob_death_hosp / prob_death_comm)
      is_death <- if (is_death_comm) rbinom(1, 1, second_chance_death) == 1L else FALSE
      t_outcome <- t_hosp + if (is_death) hospitalisation_to_death_fn(1) else
        hospitalisation_to_recovery_fn(1)
      outcome_loc <- "hospital"
    } else {
      is_hosp   <- FALSE
      t_hosp    <- NA_real_
      is_death  <- is_death_comm
      t_outcome <- t_comm_outcome
      outcome_loc <- "community"
    }

    # Funeral safety for deaths
    funeral_uns <- NA
    if (is_death) {
      p_unsafe_t <- resolve_tv(p_unsafe_funeral, t_outcome)
      funeral_uns <- rbinom(1, 1, p_unsafe_t) == 1L
    }

    tdf$status[sid]             <- 2L
    tdf$generation[sid]         <- 1L
    tdf$contact_type[sid]       <- 0L
    tdf$time_infection[sid]     <- t_inf
    tdf$time_onset[sid]         <- t_onset
    tdf$symptomatic[sid]        <- TRUE
    tdf$hospitalised[sid]       <- is_hosp
    tdf$time_hospitalisation[sid] <- t_hosp
    tdf$outcome_death[sid]      <- is_death
    tdf$outcome_location[sid]   <- outcome_loc
    tdf$time_outcome[sid]       <- t_outcome
    tdf$funeral_unsafe[sid]     <- funeral_uns
  }

  active_queue     <- seed_idx
  n_cumul_infected <- length(seed_idx)
  stop_reason      <- "extinction"
  vax_deployed     <- FALSE

  # ── Pre-loop vaccination ───────────────────────────────────────────────────
  if (prob_vax > 0 && is.finite(vax_start) && vax_start <= t0) {
    eligible <- which(tdf$status == 1L)
    if (!is.null(vax_target_ages))
      eligible <- eligible[tdf$age_group[eligible] %in% vax_target_ages]
    if (length(eligible) > 0) {
      vax_ids <- eligible[rbinom(length(eligible), 1L, prob_vax) == 1L]
      tdf$vaccinated[vax_ids]    <- 1L
      tdf$time_vaccinated[vax_ids] <- t0
    }
    vax_deployed <- TRUE
  }

  # ── Main simulation loop ───────────────────────────────────────────────────
  event_count <- 0L

  while (length(active_queue) > 0) {

    if (monitoring_console) {
      event_count <- event_count + 1L
      if (event_count %% 100L == 0L)
        message(sprintf("  [event %6d] t = %6.1f | active = %4d | total = %6d",
                        event_count,
                        min(tdf$time_infection[active_queue], na.rm = TRUE),
                        length(active_queue),
                        n_cumul_infected))
    }

    # Stopping criteria
    if (n_cumul_infected >= max_infected) { stop_reason <- "max_infected"; break }
    t_earliest <- min(tdf$time_infection[active_queue], na.rm = TRUE)
    if (t_earliest >= max_time)           { stop_reason <- "max_time";     break }

    # One-shot vaccination trigger
    if (prob_vax > 0 && is.finite(vax_start) && !vax_deployed &&
        t_earliest >= vax_start) {
      eligible <- which(tdf$status == 1L)
      if (!is.null(vax_target_ages))
        eligible <- eligible[tdf$age_group[eligible] %in% vax_target_ages]
      if (length(eligible) > 0) {
        vax_ids <- eligible[rbinom(length(eligible), 1L, prob_vax) == 1L]
        tdf$vaccinated[vax_ids]    <- 1L
        tdf$time_vaccinated[vax_ids] <- vax_start
      }
      vax_deployed <- TRUE
    }

    # Pick earliest case
    active_queue <- active_queue[order(tdf$time_infection[active_queue])]
    idx          <- active_queue[1]
    active_queue <- active_queue[-1]

    t_inf_idx  <- tdf$time_infection[idx]
    t_onset_idx <- tdf$time_onset[idx]
    t_hosp_idx  <- tdf$time_hospitalisation[idx]
    t_out_idx   <- tdf$time_outcome[idx]
    is_hcw_idx  <- tdf$is_hcw[idx]
    gen_idx     <- tdf$generation[idx]

    # Self-treatment (at symptom onset)
    logistical_delay <- logistical_delay_fn(1)
    if (!is.na(t_onset_idx) && t_onset_idx >= antiviral_start &&
        tdf$treated[idx] == 0L &&
        rbinom(1, 1, prob_treat_self) == 1L) {
      tdf$treated[idx]      <- 1L
      tdf$time_treated[idx] <- t_onset_idx + logistical_delay
    }

    # Self-quarantine
    if (!is.na(t_onset_idx) && t_onset_idx >= quarantine_start &&
        tdf$quarantined[idx] == 0L &&
        rbinom(1, 1, prob_quarantine_self) == 1L) {
      tdf$quarantined[idx]      <- 1L
      tdf$time_quarantined[idx] <- t_onset_idx + time_to_quarantine_fn(1)
    }

    # Source vaccine efficacy on transmission
    vax_eff_trans_idx <- if (tdf$vaccinated[idx] == 1L) vax_eff_trans else 0

    # Ring trigger time
    t_ring <- if (!is.na(t_onset_idx) &&
                  (t_onset_idx >= antiviral_start || t_onset_idx >= quarantine_start))
      t_onset_idx + logistical_delay else Inf

    # ===========================================================================
    # Phase 1: infection → hospitalisation (or outcome if not hospitalised)
    # Layer 1 (household) + Layer 2 (community)
    # HCW also: Layer 3 HCW-HCW with PPE reduction
    # ===========================================================================
    t_phase1_end <- if (!is.na(t_hosp_idx)) t_hosp_idx else t_out_idx

    phase1_contacts <- list(
      list(nbrs = hh_nbrs[[idx]],
           p_inf = if (length(hh_nbrs[[idx]]) > 0)
             p_inf_household / length(hh_nbrs[[idx]] + 1) else 0,
           ctype = 1L,
           p_treat = prob_treat_household, p_quar = prob_quarantine_household,
           treat_code = 2L, quar_code = 2L),
      list(nbrs = comm_nbrs[[idx]], p_inf = p_inf_community, ctype = 2L,
           p_treat = prob_treat_community, p_quar = prob_quarantine_community,
           treat_code = 3L, quar_code = 3L)
    )

    # HCW phase 1 hospital contacts (while working, before admission)
    if (is_hcw_idx && length(hcw_nbrs[[idx]]) > 0) {
      p_inf_hcw_preAdm <- p_inf_hcw_to_hcw * (1 - ppe_efficacy_hcw)
      phase1_contacts[[3]] <- list(
        nbrs       = hcw_nbrs[[idx]],
        p_inf      = p_inf_hcw_preAdm,
        ctype      = 3L,
        p_treat    = 0, p_quar = 0,
        treat_code = 3L, quar_code = 3L
      )
    }

    for (grp in phase1_contacts) {
      if (length(grp$nbrs) == 0) next

      # Ring 1 treatment / quarantine (direct contacts)
      if (is.finite(t_ring)) {
        for (nbr_id in grp$nbrs) {
          if (tdf$treated[nbr_id] == 0L && grp$p_treat > 0 &&
              t_ring >= antiviral_start &&
              rbinom(1, 1, grp$p_treat) == 1L)
            tdf$treated[nbr_id]      <- 1L; tdf$time_treated[nbr_id] <- t_ring
            if (tdf$quarantined[nbr_id] == 0L && grp$p_quar > 0 &&
                t_ring >= quarantine_start &&
                rbinom(1, 1, grp$p_quar) == 1L) {
              tdf$quarantined[nbr_id]      <- 1L
              tdf$time_quarantined[nbr_id] <- t_ring
            }
        }

        # Ring 2 treatment / quarantine (on-demand from adjacency lists)
        if (prob_treat_community > 0 || prob_quarantine_community > 0) {
          ring2_ids <- get_ring2(idx)
          for (nbr_id in ring2_ids) {
            if (tdf$treated[nbr_id] == 0L && prob_treat_community > 0 &&
                t_ring >= antiviral_start &&
                rbinom(1, 1, prob_treat_community) == 1L) {
              tdf$treated[nbr_id]      <- 1L
              tdf$time_treated[nbr_id] <- t_ring
            }
            if (tdf$quarantined[nbr_id] == 0L && prob_quarantine_community > 0 &&
                t_ring >= quarantine_start &&
                rbinom(1, 1, prob_quarantine_community) == 1L) {
              tdf$quarantined[nbr_id]      <- 1L
              tdf$time_quarantined[nbr_id] <- t_ring
            }
          }
        }
      }

      for (nbr_id in grp$nbrs) {
        if (tdf$status[nbr_id] != 1L) next

        # Transmission time within phase 1 window
        t_inf_nbr <- t_inf_idx + generation_time_fn(1)
        if (t_inf_nbr > t_phase1_end) next  # Outside phase 1 window

        # Gate 1: transmission probability
        if (rbinom(1, 1, grp$p_inf) == 0L) next

        # Gate 2: source antiviral + vaccine (max rule)
        eff_trans_src <- max(
          if (!is.null(pk_params) && tdf$treated[idx] == 1L)
            compute_drug_eff_source(tdf, idx, pd_params_trans, pk_params) else 0,
          vax_eff_trans_idx
        )
        if (eff_trans_src > 0 && rbinom(1, 1, eff_trans_src) == 1L) next

        # Gate 3: source quarantine
        if (tdf$quarantined[idx] == 1L &&
            !is.na(tdf$time_quarantined[idx]) &&
            t_inf_nbr > tdf$time_quarantined[idx] &&
            rbinom(1, 1, quarantine_efficacy) == 1L) next

        # Gate 4: target antiviral + vaccine (max rule)
        eff_inf_tgt <- max(
          if (!is.null(pk_params) && tdf$treated[nbr_id] == 1L)
            compute_drug_eff_target(tdf, nbr_id, t_inf_nbr,
                                    pd_params_inf, pk_params) else 0,
          if (tdf$vaccinated[nbr_id] == 1L) vax_eff_inf else 0
        )
        if (eff_inf_tgt > 0 && rbinom(1, 1, eff_inf_tgt) == 1L) next

        # Infection confirmed → assign natural history
        tdf <- infect_individual(
          tdf, nbr_id, idx, t_inf_nbr, gen_idx + 1L, grp$ctype,
          incubation_period_fn, onset_to_hospitalisation_fn,
          onset_to_death_fn, onset_to_recovery_fn,
          hospitalisation_to_death_fn, hospitalisation_to_recovery_fn,
          prob_hospitalised_genPop, prob_hospitalised_hcw,
          prob_death_comm, prob_death_hosp, p_unsafe_funeral,
          resolve_tv
        )
        active_queue     <- c(active_queue, nbr_id)
        n_cumul_infected <- n_cumul_infected + 1L
      }
    }

    # ===========================================================================
    # Phase 2: hospitalisation → outcome
    # genPop: patient → HCW (admission_nbrs)
    # HCW:    HCW-HCW (hcw_nbrs, no PPE after admission)
    # ===========================================================================
    if (!is.na(t_hosp_idx)) {

      if (is_hcw_idx) {
        phase2_nbrs <- hcw_nbrs[[idx]]
        p_inf_p2    <- p_inf_hcw_to_hcw
      } else {
        phase2_nbrs <- admission_nbrs[[idx]]
        p_inf_p2    <- if (length(phase2_nbrs) > 0)
          p_inf_patient_to_hcw / length(phase2_nbrs) else 0
      }

      for (nbr_id in phase2_nbrs) {
        if (tdf$status[nbr_id] != 1L) next

        t_inf_nbr <- t_hosp_idx + generation_time_fn(1)
        if (t_inf_nbr > t_out_idx) next  # Outside phase 2 window

        if (rbinom(1, 1, p_inf_p2) == 0L) next

        eff_trans_src <- max(
          if (!is.null(pk_params) && tdf$treated[idx] == 1L)
            compute_drug_eff_source(tdf, idx, pd_params_trans, pk_params) else 0,
          vax_eff_trans_idx
        )
        if (eff_trans_src > 0 && rbinom(1, 1, eff_trans_src) == 1L) next

        if (tdf$quarantined[idx] == 1L &&
            !is.na(tdf$time_quarantined[idx]) &&
            t_inf_nbr > tdf$time_quarantined[idx] &&
            rbinom(1, 1, quarantine_efficacy) == 1L) next

        eff_inf_tgt <- max(
          if (!is.null(pk_params) && tdf$treated[nbr_id] == 1L)
            compute_drug_eff_target(tdf, nbr_id, t_inf_nbr,
                                    pd_params_inf, pk_params) else 0,
          if (tdf$vaccinated[nbr_id] == 1L) vax_eff_inf else 0
        )
        if (eff_inf_tgt > 0 && rbinom(1, 1, eff_inf_tgt) == 1L) next

        tdf <- infect_individual(
          tdf, nbr_id, idx, t_inf_nbr, gen_idx + 1L, 3L,
          incubation_period_fn, onset_to_hospitalisation_fn,
          onset_to_death_fn, onset_to_recovery_fn,
          hospitalisation_to_death_fn, hospitalisation_to_recovery_fn,
          prob_hospitalised_genPop, prob_hospitalised_hcw,
          prob_death_comm, prob_death_hosp, p_unsafe_funeral,
          resolve_tv
        )
        active_queue     <- c(active_queue, nbr_id)
        n_cumul_infected <- n_cumul_infected + 1L
      }
    }

    # ===========================================================================
    # Phase 3: funeral (death only)
    # Layer 1 + Layer 2 with unsafe/safe transmission probability
    # ===========================================================================
    if (isTRUE(tdf$outcome_death[idx])) {

      p_inf_funeral <- if (isTRUE(tdf$funeral_unsafe[idx]))
        p_inf_funeral_unsafe else p_inf_funeral_safe

      funeral_nbrs <- c(hh_nbrs[[idx]], comm_nbrs[[idx]])

      for (nbr_id in funeral_nbrs) {
        if (tdf$status[nbr_id] != 1L) next

        if (rbinom(1, 1, p_inf_funeral) == 0L) next

        # Funeral transmission time: shortly after death
        t_inf_nbr <- t_out_idx + generation_time_fn(1)

        eff_inf_tgt <- if (tdf$vaccinated[nbr_id] == 1L) vax_eff_inf else 0
        if (eff_inf_tgt > 0 && rbinom(1, 1, eff_inf_tgt) == 1L) next

        tdf <- infect_individual(
          tdf, nbr_id, idx, t_inf_nbr, gen_idx + 1L, 4L,
          incubation_period_fn, onset_to_hospitalisation_fn,
          onset_to_death_fn, onset_to_recovery_fn,
          hospitalisation_to_death_fn, hospitalisation_to_recovery_fn,
          prob_hospitalised_genPop, prob_hospitalised_hcw,
          prob_death_comm, prob_death_hosp, p_unsafe_funeral,
          resolve_tv
        )
        active_queue     <- c(active_queue, nbr_id)
        n_cumul_infected <- n_cumul_infected + 1L
      }
    }

    # Mark individual as resolved
    tdf$status[idx] <- 3L

  } # end main loop

  # ── Output ─────────────────────────────────────────────────────────────────
  infected_df <- tdf[!is.na(tdf$time_infection), ]
  infected_df <- infected_df[order(infected_df$time_infection,
                                   infected_df$person_id), ]

  list(
    full             = tdf,
    infected         = infected_df,
    stop_reason      = stop_reason,
    n_cumul_infected = n_cumul_infected,
    n_active_at_stop = length(active_queue)
  )
}

# ==============================================================================
# Internal helper: assign natural history to a newly infected individual
# ==============================================================================
infect_individual <- function(
    tdf, idx, ancestor_idx, t_inf, generation, contact_type,
    incubation_period_fn, onset_to_hospitalisation_fn,
    onset_to_death_fn, onset_to_recovery_fn,
    hospitalisation_to_death_fn, hospitalisation_to_recovery_fn,
    prob_hospitalised_genPop, prob_hospitalised_hcw,
    prob_death_comm, prob_death_hosp, p_unsafe_funeral,
    resolve_tv
) {
  # Ebola: always symptomatic
  inc_p   <- incubation_period_fn(1)
  t_onset <- t_inf + inc_p

  # Hospitalisation
  p_hosp_t <- resolve_tv(
    if (tdf$is_hcw[idx]) prob_hospitalised_hcw else prob_hospitalised_genPop,
    t_onset
  )
  is_hosp <- rbinom(1, 1, p_hosp_t) == 1L
  t_hosp  <- if (is_hosp) t_onset + onset_to_hospitalisation_fn(1) else NA_real_

  # Community outcome
  is_death_comm  <- rbinom(1, 1, prob_death_comm) == 1L
  t_comm_outcome <- if (is_death_comm) t_onset + onset_to_death_fn(1) else
    t_onset + onset_to_recovery_fn(1)

  # Final outcome
  if (is_hosp && !is.na(t_hosp) && t_hosp < t_comm_outcome) {
    second_chance <- prob_death_hosp / prob_death_comm
    is_death      <- if (is_death_comm) rbinom(1, 1, second_chance) == 1L else FALSE
    t_outcome     <- t_hosp + if (is_death) hospitalisation_to_death_fn(1) else
      hospitalisation_to_recovery_fn(1)
    outcome_loc   <- "hospital"
  } else {
    is_hosp     <- FALSE
    t_hosp      <- NA_real_
    is_death    <- is_death_comm
    t_outcome   <- t_comm_outcome
    outcome_loc <- "community"
  }

  # Funeral safety
  funeral_uns <- NA
  if (is_death) {
    p_unsafe_t  <- resolve_tv(p_unsafe_funeral, t_outcome)
    funeral_uns <- rbinom(1, 1, p_unsafe_t) == 1L
  }

  tdf$status[idx]               <- 2L
  tdf$generation[idx]           <- generation
  tdf$ancestor_id[idx]          <- tdf$person_id[ancestor_idx]
  tdf$contact_type[idx]         <- contact_type
  tdf$time_infection[idx]       <- t_inf
  tdf$time_onset[idx]           <- t_onset
  tdf$symptomatic[idx]          <- TRUE
  tdf$hospitalised[idx]         <- is_hosp
  tdf$time_hospitalisation[idx] <- t_hosp
  tdf$outcome_death[idx]        <- is_death
  tdf$outcome_location[idx]     <- outcome_loc
  tdf$time_outcome[idx]         <- t_outcome
  tdf$funeral_unsafe[idx]       <- funeral_uns

  tdf
}
