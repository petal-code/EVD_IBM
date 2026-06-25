# ==============================================================================
# COD_IBM_ebola_sim.R
# Purpose:
#   Individual-based transmission model for DRC
#   Network structure: pre-built contact network (Layer 1/2/3)
#
# Architecture:
#   - sim_prep: pre-built adj lists loaded once, O(1) neighbor lookup
#   - Individual state stored as flat vectors (not data.frame) for speed
#   - Pre-allocated queue with grow-on-demand
#   - infect_individual() defined inside main function → <<- direct vector writes
#   - p_inf varies by frequency stratum (daily/weekly/monthly)
#     via 1-(1-p)^n over 3-week window
#
# contact_type codes:
#   0L = index
#   1L = household
#   2L = community close  (is_physical=0)
#   3L = community physical (is_physical=1)
#   4L = hospital
#   5L = funeral household
#   6L = funeral community
# ==============================================================================

ebola_network_sim <- function(

  # ── sim_prep (from C2_sim_p0_simprep.R) ──────────────────────────────────
  sim_prep,   # list(N, pid_to_idx,
  #      hh_nbrs,
  #      comm_close_daily_nbrs, comm_close_weekly_nbrs, comm_close_monthly_nbrs,
  #      comm_phys_daily_nbrs,  comm_phys_weekly_nbrs,  comm_phys_monthly_nbrs,
  #      hcw_nbrs, adm_nbrs)
  nodes,

  # ── Seeding ───────────────────────────────────────────────────────────────
  seeding_cases    = 1L,
  seeding_ids      = NULL,
  t0               = 0,

  # ── Natural history ───────────────────────────────────────────────────────
  incubation_period_fn,
  onset_to_hospitalisation_fn,
  onset_to_death_fn,
  onset_to_recovery_fn,
  hospitalisation_to_death_fn,
  hospitalisation_to_recovery_fn,
  generation_time_fn,

  # ── Disease severity ──────────────────────────────────────────────────────
  prob_symptomatic             = 1.0,
  prob_hospitalised_genPop,
  prob_hospitalised_hcw,
  prob_death_comm,
  prob_death_hosp,

  # ── Transmission probabilities ────────────────────────────────────────────
  # Per single contact-event probability.
  # 3-week effective prob = 1-(1-p)^n computed internally:
  #   daily n=21, weekly n=3, monthly n=0.75
  p_inf_household_close,
  p_inf_household_physical,
  close_only_home,    # 16x16: Prem home × (1 - blended_ratio$home), from p6
  phys_only_home,     # 16x16: Prem home × blended_ratio$home, from p6

  p_inf_community_close_daily    = 0,
  p_inf_community_close_weekly   = 0,
  p_inf_community_close_monthly  = 0,

  p_inf_community_physical_daily,
  p_inf_community_physical_weekly,
  p_inf_community_physical_monthly,

  p_inf_hcw_to_hcw,
  p_inf_patient_to_hcw,

  # ── Funeral ───────────────────────────────────────────────────────────────
  funeral_avg,
  funeral_k,
  p_unsafe_funeral,
  p_inf_funeral_household,
  p_inf_funeral_community,
  funeral_unsafe_multiplier    = 1.0,
  funeral_safe_multiplier      = 0.1,

  # ── HCW-specific ──────────────────────────────────────────────────────────
  ppe_efficacy_hcw              = 0,
  prob_hospital_cond_hcw_preAdm = 0.5,

  # ── Treatment ─────────────────────────────────────────────────────────────
  antiviral_start               = Inf,
  drug_eff_trans                = 0,
  drug_eff_inf                  = 0,
  drug_eff_death                = 0,
  pk_params                     = NULL,
  pd_params_trans               = NULL,
  pd_params_inf                 = NULL,

  prob_treat_self               = 0,
  hcw_prep_start                = Inf,
  prob_treat_hcw_prep           = 0,
  prob_treat_hcw_pep            = 0,

  prob_treat_given_trace_household = 0,
  prob_treat_given_trace_close     = 0,
  prob_treat_given_trace_physical  = 0,
  prob_treat_given_trace_funeral   = 0,

  # ── Quarantine ────────────────────────────────────────────────────────────
  quarantine_start              = Inf,
  quarantine_efficacy           = 0,
  prob_quarantine_self          = 0,
  prob_quarantine_given_trace_household = 0,
  prob_quarantine_given_trace_close     = 0,
  prob_quarantine_given_trace_physical  = 0,
  prob_quarantine_given_trace_funeral   = 0,
  time_to_quarantine_fn         = function(n) rep(0, n),
  logistical_delay_fn           = function(n) rep(0, n),

  # ── Ring tracing ──────────────────────────────────────────────────────────
  prob_trace_household          = 0,
  prob_trace_close              = 0,
  prob_trace_physical           = 0,
  prob_trace_funeral            = 0,

  # ── Stopping ──────────────────────────────────────────────────────────────
  max_infected                  = Inf,
  max_time                      = Inf,

  # ── Misc ──────────────────────────────────────────────────────────────────
  seed                          = 42L,
  monitoring_console            = FALSE

) {

  set.seed(seed)

  # ── Helpers ───────────────────────────────────────────────────────────────
  resolve_tv <- function(param, t) {
    if (is.function(param)) param(t) else param
  }

  # 3-week cumulative p from per-contact p and frequency n
  p_cumul <- function(p, n) {
    if (p <= 0) return(0)
    if (p >= 1) return(1)
    1 - (1 - p)^n
  }

  # Pre-compute 3-week effective p_inf per stratum (once at start)
  p_eff_close_daily    <- p_cumul(p_inf_community_close_daily,    21)
  p_eff_close_weekly   <- p_cumul(p_inf_community_close_weekly,   3)
  p_eff_close_monthly  <- p_cumul(p_inf_community_close_monthly,  0.75)
  p_eff_phys_daily     <- p_cumul(p_inf_community_physical_daily,   21)
  p_eff_phys_weekly    <- p_cumul(p_inf_community_physical_weekly,  3)
  p_eff_phys_monthly   <- p_cumul(p_inf_community_physical_monthly, 0.75)

  # PK/PD drug efficacy (one-compartment PK + Emax PD)
  compute_drug_eff <- function(t_current, t_treated, pd_params, pk_params) {
    if (is.null(pk_params) || is.null(pd_params)) return(0)
    dt <- t_current - t_treated
    if (is.na(dt) || dt < 0) return(0)
    conc <- pk_params$dose * exp(-pk_params$ke * dt)
    pd_params$emax * conc / (pd_params$ec50 + conc)
  }

  # ── Load sim_prep (O(1) neighbor lookup) ──────────────────────────────────
  message("Loading sim_prep...")
  N                       <- sim_prep$N
  pid_to_idx              <- sim_prep$pid_to_idx
  hh_nbrs                 <- sim_prep$hh_nbrs
  comm_close_daily_nbrs   <- sim_prep$comm_close_daily_nbrs
  comm_close_weekly_nbrs  <- sim_prep$comm_close_weekly_nbrs
  comm_close_monthly_nbrs <- sim_prep$comm_close_monthly_nbrs
  comm_phys_daily_nbrs    <- sim_prep$comm_phys_daily_nbrs
  comm_phys_weekly_nbrs   <- sim_prep$comm_phys_weekly_nbrs
  comm_phys_monthly_nbrs  <- sim_prep$comm_phys_monthly_nbrs
  hcw_nbrs                <- sim_prep$hcw_nbrs
  adm_nbrs                <- sim_prep$adm_nbrs
  message(sprintf("  Loaded for %d individuals", N))

  # ── Ring helpers ──────────────────────────────────────────────────────────
  get_all_comm_nbrs <- function(idx) {
    unique(c(comm_close_daily_nbrs[[idx]],  comm_phys_daily_nbrs[[idx]],
             comm_close_weekly_nbrs[[idx]], comm_phys_weekly_nbrs[[idx]],
             comm_close_monthly_nbrs[[idx]],comm_phys_monthly_nbrs[[idx]]))
  }

  get_ring1_all <- function(idx)
    unique(c(hh_nbrs[[idx]], get_all_comm_nbrs(idx)))

  get_ring2_all <- function(idx) {
    r1 <- get_ring1_all(idx)
    if (length(r1) == 0L) return(integer(0))
    unique(unlist(lapply(r1, get_ring1_all), use.names = FALSE))
  }

  # ── State vectors ─────────────────────────────────────────────────────────
  status               <- rep(1L,           N)
  time_infection       <- rep(NA_real_,     N)
  time_onset           <- rep(NA_real_,     N)
  time_hospitalisation <- rep(NA_real_,     N)
  time_outcome         <- rep(NA_real_,     N)
  hospitalised         <- rep(FALSE,        N)
  outcome_death        <- rep(NA,           N)
  outcome_location     <- rep(NA_character_,N)
  funeral_unsafe_vec   <- rep(NA,           N)
  funeral_attended_for <- rep(NA_integer_,  N)
  funeral_role_vec     <- rep(NA_character_,N)
  generation           <- rep(NA_integer_,  N)
  ancestor_id          <- rep(NA_integer_,  N)
  contact_type_vec     <- rep(NA_integer_,  N)
  treated              <- rep(0L,           N)
  time_treated         <- rep(NA_real_,     N)
  quarantined          <- rep(0L,           N)
  time_quarantined     <- rep(NA_real_,     N)
  traced_via           <- rep(NA_character_,N)

  v_person_id   <- nodes$person_id
  v_hh_id       <- nodes$hh_id
  v_cell_id     <- nodes$cell_id
  v_hospital_id <- nodes$hospital_id
  v_age_group   <- nodes$age_group
  v_is_hcw      <- nodes$is_hcw

  # ── infect_individual (inner — writes via <<-) ────────────────────────────
  infect_individual <- function(idx, ancestor_idx, t_inf,
                                generation_val, contact_type_val) {
    t_onset <- t_inf + incubation_period_fn(1)

    p_hosp_t <- resolve_tv(
      if (isTRUE(v_is_hcw[idx])) prob_hospitalised_hcw else prob_hospitalised_genPop,
      t_onset)
    is_hosp <- rbinom(1, 1, p_hosp_t) == 1L
    t_hosp  <- if (is_hosp) t_onset + onset_to_hospitalisation_fn(1) else NA_real_

    eff_death        <- drug_eff_death * (treated[idx] == 1L)
    p_death_comm_eff <- prob_death_comm * (1 - eff_death)
    p_death_hosp_eff <- prob_death_hosp * (1 - eff_death)

    is_death_comm  <- rbinom(1, 1, p_death_comm_eff) == 1L
    t_comm_outcome <- if (is_death_comm) t_onset + onset_to_death_fn(1) else
      t_onset + onset_to_recovery_fn(1)

    if (is_hosp && !is.na(t_hosp) && t_hosp < t_comm_outcome) {
      second_chance <- p_death_hosp_eff / max(p_death_comm_eff, 1e-9)
      is_death      <- if (is_death_comm)
        rbinom(1, 1, min(second_chance, 1)) == 1L else FALSE
      t_outcome   <- t_hosp + if (is_death) hospitalisation_to_death_fn(1) else
        hospitalisation_to_recovery_fn(1)
      outcome_loc <- "hospital"
    } else {
      is_hosp     <- FALSE
      t_hosp      <- NA_real_
      is_death    <- is_death_comm
      t_outcome   <- t_comm_outcome
      outcome_loc <- "community"
    }

    funeral_uns <- NA
    if (is_death) {
      p_unsafe_t  <- resolve_tv(p_unsafe_funeral, t_outcome)
      funeral_uns <- rbinom(1, 1, p_unsafe_t) == 1L
    }

    status[idx]               <<- 2L
    generation[idx]           <<- generation_val
    ancestor_id[idx]          <<- v_person_id[ancestor_idx]
    contact_type_vec[idx]     <<- contact_type_val
    time_infection[idx]       <<- t_inf
    time_onset[idx]           <<- t_onset
    hospitalised[idx]         <<- is_hosp
    time_hospitalisation[idx] <<- t_hosp
    outcome_death[idx]        <<- is_death
    outcome_location[idx]     <<- outcome_loc
    time_outcome[idx]         <<- t_outcome
    funeral_unsafe_vec[idx]   <<- funeral_uns

    invisible(NULL)
  }

  # ── Pre-allocated queue ───────────────────────────────────────────────────
  queue_cap   <- max(1000L, N %/% 100L)
  queue_buf   <- integer(queue_cap)
  queue_times <- numeric(queue_cap)
  queue_head  <- 1L
  queue_tail  <- 0L

  enqueue <- function(idx, t_inf) {
    queue_tail <<- queue_tail + 1L
    if (queue_tail > queue_cap) {
      queue_cap   <<- queue_cap * 2L
      queue_buf   <<- c(queue_buf,   integer(queue_cap %/% 2L))
      queue_times <<- c(queue_times, numeric(queue_cap %/% 2L))
    }
    queue_buf[queue_tail]   <<- idx
    queue_times[queue_tail] <<- t_inf
  }

  dequeue_earliest <- function() {
    ar      <- queue_head:queue_tail
    min_pos <- ar[which.min(queue_times[ar])]
    idx     <- queue_buf[min_pos]
    queue_buf[min_pos]   <<- queue_buf[queue_head]
    queue_times[min_pos] <<- queue_times[queue_head]
    queue_head <<- queue_head + 1L
    idx
  }

  queue_size <- function() max(0L, queue_tail - queue_head + 1L)

  # ── HCW PrEP at or before t0 ─────────────────────────────────────────────
  hcw_prep_deployed <- FALSE
  if (prob_treat_hcw_prep > 0 && is.finite(hcw_prep_start) &&
      hcw_prep_start <= t0) {
    hcw_idx  <- which(v_is_hcw)
    prep_idx <- hcw_idx[rbinom(length(hcw_idx), 1L, prob_treat_hcw_prep) == 1L]
    treated[prep_idx]      <- 1L
    time_treated[prep_idx] <- hcw_prep_start
    hcw_prep_deployed      <- TRUE
    message(sprintf("  HCW PrEP at t0: %d HCWs treated", length(prep_idx)))
  }

  # ── Seed cases ────────────────────────────────────────────────────────────
  seed_idx <- if (!is.null(seeding_ids)) {
    pid_to_idx[seeding_ids + 1L]
  } else {
    sample(seq_len(N), seeding_cases, replace = FALSE)
  }

  for (sid in seed_idx) {
    infect_individual(sid, sid, t0, 1L, 0L)
    enqueue(sid, t0)
  }
  n_cumul_infected <- length(seed_idx)
  stop_reason      <- "extinction"

  # ── Gate check ────────────────────────────────────────────────────────────
  passes_gates <- function(idx, nbr_id, t_inf_nbr, p_inf, eff_trans_src) {
    if (rbinom(1, 1, p_inf) == 0L) return(FALSE)
    if (eff_trans_src > 0 && rbinom(1, 1, eff_trans_src) == 1L) return(FALSE)
    if (quarantined[idx] == 1L && !is.na(time_quarantined[idx]) &&
        t_inf_nbr > time_quarantined[idx] &&
        rbinom(1, 1, quarantine_efficacy) == 1L) return(FALSE)
    eff_inf_tgt <- if (treated[nbr_id] == 1L) {
      if (!is.null(pk_params))
        compute_drug_eff(t_inf_nbr, time_treated[nbr_id], pd_params_inf, pk_params)
      else drug_eff_inf
    } else 0
    if (eff_inf_tgt > 0 && rbinom(1, 1, eff_inf_tgt) == 1L) return(FALSE)
    TRUE
  }

  # ── Ring intervention helper ──────────────────────────────────────────────
  apply_ring <- function(nbr_ids, trace_prob, treat_prob, quar_prob,
                         trace_label, t_ring, delay) {
    for (nbr_id in nbr_ids) {
      if (!is.na(traced_via[nbr_id])) next
      if (rbinom(1, 1, trace_prob) == 0L) next
      traced_via[nbr_id] <<- trace_label
      if (treat_prob > 0 && treated[nbr_id] == 0L &&
          is.finite(t_ring) && t_ring >= antiviral_start &&
          rbinom(1, 1, treat_prob) == 1L) {
        treated[nbr_id]      <<- 1L
        time_treated[nbr_id] <<- t_ring + delay
      }
      if (quar_prob > 0 && quarantined[nbr_id] == 0L &&
          is.finite(t_ring) && t_ring >= quarantine_start &&
          rbinom(1, 1, quar_prob) == 1L) {
        quarantined[nbr_id]      <<- 1L
        time_quarantined[nbr_id] <<- t_ring + time_to_quarantine_fn(1)
      }
    }
  }

  # ── Main simulation loop ──────────────────────────────────────────────────
  event_count <- 0L

  while (queue_size() > 0L) {

    if (monitoring_console) {
      event_count <- event_count + 1L
      if (event_count %% 100L == 0L)
        message(sprintf("  [event %6d] active=%4d | total=%6d",
                        event_count, queue_size(), n_cumul_infected))
    }

    if (n_cumul_infected >= max_infected) { stop_reason <- "max_infected"; break }
    t_earliest <- queue_times[queue_head]
    if (t_earliest >= max_time)           { stop_reason <- "max_time";     break }

    # HCW PrEP mid-simulation trigger
    if (prob_treat_hcw_prep > 0 && is.finite(hcw_prep_start) &&
        !hcw_prep_deployed && t_earliest >= hcw_prep_start) {
      hcw_idx  <- which(v_is_hcw & status == 1L)
      prep_idx <- hcw_idx[rbinom(length(hcw_idx), 1L, prob_treat_hcw_prep) == 1L]
      treated[prep_idx]      <- 1L
      time_treated[prep_idx] <- hcw_prep_start
      hcw_prep_deployed      <- TRUE
      message(sprintf("  HCW PrEP deployed at t=%.1f: %d HCWs",
                      hcw_prep_start, length(prep_idx)))
    }

    idx <- dequeue_earliest()

    t_inf_idx   <- time_infection[idx]
    t_onset_idx <- time_onset[idx]
    t_hosp_idx  <- time_hospitalisation[idx]
    t_out_idx   <- time_outcome[idx]
    is_hcw_idx  <- isTRUE(v_is_hcw[idx])
    gen_idx     <- generation[idx]

    eff_trans_src <- if (treated[idx] == 1L) {
      if (!is.null(pk_params))
        compute_drug_eff(t_inf_idx, time_treated[idx], pd_params_trans, pk_params)
      else drug_eff_trans
    } else 0

    logistical_delay <- logistical_delay_fn(1)

    # Self-treatment
    if (!is.na(t_onset_idx) && t_onset_idx >= antiviral_start &&
        treated[idx] == 0L && rbinom(1, 1, prob_treat_self) == 1L) {
      treated[idx]      <- 1L
      time_treated[idx] <- t_onset_idx + logistical_delay
    }

    # Self-quarantine
    if (!is.na(t_onset_idx) && t_onset_idx >= quarantine_start &&
        quarantined[idx] == 0L && rbinom(1, 1, prob_quarantine_self) == 1L) {
      quarantined[idx]      <- 1L
      time_quarantined[idx] <- t_onset_idx + time_to_quarantine_fn(1)
    }

    t_ring <- if (!is.na(t_onset_idx) &&
                  (t_onset_idx >= antiviral_start ||
                   t_onset_idx >= quarantine_start))
      t_onset_idx + logistical_delay else Inf

    # Ring interventions (O(1) lookup from sim_prep)
    if (is.finite(t_ring)) {
      apply_ring(hh_nbrs[[idx]],
                 prob_trace_household,
                 prob_treat_given_trace_household,
                 prob_quarantine_given_trace_household,
                 "household", t_ring, logistical_delay)

      apply_ring(unique(c(comm_close_daily_nbrs[[idx]],
                          comm_close_weekly_nbrs[[idx]],
                          comm_close_monthly_nbrs[[idx]])),
                 prob_trace_close,
                 prob_treat_given_trace_close,
                 prob_quarantine_given_trace_close,
                 "close", t_ring, logistical_delay)

      apply_ring(unique(c(comm_phys_daily_nbrs[[idx]],
                          comm_phys_weekly_nbrs[[idx]],
                          comm_phys_monthly_nbrs[[idx]])),
                 prob_trace_physical,
                 prob_treat_given_trace_physical,
                 prob_quarantine_given_trace_physical,
                 "physical", t_ring, logistical_delay)
    }

    # =========================================================================
    # Phase 1: infection → hospitalisation (or outcome)
    # =========================================================================
    t_phase1_end <- if (!is.na(t_hosp_idx)) t_hosp_idx else t_out_idx

    # Household: p_eff per edge based on age-group combination
    # p_eff = p_inf_household_close    * close_only_home[age_i, age_j]
    #       + p_inf_household_physical * phys_only_home[age_i, age_j]
    age_i_idx <- v_age_group[idx]
    hh_nbrs_now <- hh_nbrs[[idx]]
    for (nbr_id in hh_nbrs_now) {
      if (status[nbr_id] != 1L) next
      age_j_idx <- v_age_group[nbr_id]
      p_eff_hh  <- p_inf_household_close    * close_only_home[age_i_idx, age_j_idx] +
        p_inf_household_physical * phys_only_home[age_i_idx,  age_j_idx]
      if (p_eff_hh <= 0) next
      t_inf_nbr <- t_inf_idx + generation_time_fn(1)
      if (t_inf_nbr > t_phase1_end) next
      if (!passes_gates(idx, nbr_id, t_inf_nbr, p_eff_hh, eff_trans_src)) next
      infect_individual(nbr_id, idx, t_inf_nbr, gen_idx + 1L, 1L)
      enqueue(nbr_id, t_inf_nbr)
      n_cumul_infected <- n_cumul_infected + 1L
    }

    phase1_groups <- list(
      list(nbrs = comm_close_daily_nbrs[[idx]],
           p_inf = p_eff_close_daily,    ctype = 2L),
      list(nbrs = comm_close_weekly_nbrs[[idx]],
           p_inf = p_eff_close_weekly,   ctype = 2L),
      list(nbrs = comm_close_monthly_nbrs[[idx]],
           p_inf = p_eff_close_monthly,  ctype = 2L),
      list(nbrs = comm_phys_daily_nbrs[[idx]],
           p_inf = p_eff_phys_daily,     ctype = 3L),
      list(nbrs = comm_phys_weekly_nbrs[[idx]],
           p_inf = p_eff_phys_weekly,    ctype = 3L),
      list(nbrs = comm_phys_monthly_nbrs[[idx]],
           p_inf = p_eff_phys_monthly,   ctype = 3L)
    )

    # HCW pre-admission hospital contacts
    if (is_hcw_idx && length(hcw_nbrs[[idx]]) > 0 &&
        rbinom(1, 1, prob_hospital_cond_hcw_preAdm) == 1L) {
      phase1_groups[[length(phase1_groups) + 1L]] <- list(
        nbrs  = hcw_nbrs[[idx]],
        p_inf = p_inf_hcw_to_hcw * (1 - ppe_efficacy_hcw),
        ctype = 4L)
    }

    for (grp in phase1_groups) {
      if (length(grp$nbrs) == 0L || grp$p_inf <= 0) next
      for (nbr_id in grp$nbrs) {
        if (status[nbr_id] != 1L) next
        t_inf_nbr <- t_inf_idx + generation_time_fn(1)
        if (t_inf_nbr > t_phase1_end) next
        if (!passes_gates(idx, nbr_id, t_inf_nbr, grp$p_inf, eff_trans_src)) next
        infect_individual(nbr_id, idx, t_inf_nbr, gen_idx + 1L, grp$ctype)
        enqueue(nbr_id, t_inf_nbr)
        n_cumul_infected <- n_cumul_infected + 1L
      }
    }

    # =========================================================================
    # Phase 2: hospitalisation → outcome
    # =========================================================================
    if (!is.na(t_hosp_idx)) {
      if (is_hcw_idx) {
        phase2_nbrs <- hcw_nbrs[[idx]]
        p_inf_p2    <- if (length(phase2_nbrs) > 0)
          p_inf_hcw_to_hcw / length(phase2_nbrs) else 0
      } else {
        phase2_nbrs <- adm_nbrs[[idx]]
        p_inf_p2    <- if (length(phase2_nbrs) > 0)
          p_inf_patient_to_hcw / length(phase2_nbrs) else 0
      }

      for (nbr_id in phase2_nbrs) {
        if (status[nbr_id] != 1L) next
        t_inf_nbr <- t_hosp_idx + generation_time_fn(1)
        if (t_inf_nbr > t_out_idx) next
        if (!passes_gates(idx, nbr_id, t_inf_nbr, p_inf_p2, eff_trans_src)) next

        # HCW PEP on exposure
        if (isTRUE(v_is_hcw[nbr_id]) && prob_treat_hcw_pep > 0 &&
            treated[nbr_id] == 0L &&
            rbinom(1, 1, prob_treat_hcw_pep) == 1L) {
          treated[nbr_id]      <- 1L
          time_treated[nbr_id] <- t_inf_nbr + logistical_delay
        }

        infect_individual(nbr_id, idx, t_inf_nbr, gen_idx + 1L, 4L)
        enqueue(nbr_id, t_inf_nbr)
        n_cumul_infected <- n_cumul_infected + 1L
      }
    }

    # =========================================================================
    # Phase 3: funeral (death only)
    # NegBin attendee count; recruit household → 1-ring → 2-ring
    # =========================================================================
    if (isTRUE(outcome_death[idx])) {

      funedgenum   <- rnbinom(1, size = funeral_k, mu = funeral_avg)
      is_unsafe    <- isTRUE(funeral_unsafe_vec[idx])
      funeral_mult <- if (is_unsafe) funeral_unsafe_multiplier else funeral_safe_multiplier

      # Pool 1: household (all attend)
      hh_att    <- hh_nbrs[[idx]]
      remaining <- max(0L, funedgenum - length(hh_att))

      # Pool 2: 1-ring community
      comm_1ring <- setdiff(get_all_comm_nbrs(idx), hh_att)
      comm_att   <- integer(0)
      if (remaining > 0L && length(comm_1ring) > 0L) {
        n_draw    <- min(remaining, length(comm_1ring))
        comm_att  <- sample(comm_1ring, n_draw, replace = FALSE)
        remaining <- remaining - n_draw
      }

      # Pool 3: 2-ring
      ring2_att <- integer(0)
      if (remaining > 0L) {
        ring2_all <- setdiff(get_ring2_all(idx), c(hh_att, comm_att, idx))
        if (length(ring2_all) > 0L) {
          n_draw    <- min(remaining, length(ring2_all))
          ring2_att <- sample(ring2_all, n_draw, replace = FALSE)
        }
      }

      all_comm_att <- unique(c(comm_att, ring2_att))

      # Record funeral attendance
      # Use local() with explicit parent environment assignment to avoid
      # <<- scope issues caused by get_ring2_all() creating nested environments
      fa  <- funeral_attended_for
      frv <- funeral_role_vec
      for (att_id in c(hh_att, all_comm_att)) {
        if (is.na(fa[att_id])) {
          fa[att_id]  <- v_person_id[idx]
          frv[att_id] <- if (att_id %in% hh_att) "household" else "community"
        }
      }
      funeral_attended_for <<- fa
      funeral_role_vec     <<- frv
      rm(fa, frv)

      # Funeral ring tracing
      if (is.finite(t_ring)) {
        apply_ring(c(hh_att, all_comm_att),
                   prob_trace_funeral,
                   prob_treat_given_trace_funeral,
                   prob_quarantine_given_trace_funeral,
                   "funeral", t_ring, logistical_delay)
      }

      # Transmission
      p_inf_f_hh   <- p_inf_funeral_household * funeral_mult
      p_inf_f_comm <- p_inf_funeral_community * funeral_mult

      for (nbr_id in hh_att) {
        if (status[nbr_id] != 1L) next
        t_inf_nbr <- t_out_idx + generation_time_fn(1)
        if (!passes_gates(idx, nbr_id, t_out_idx, p_inf_f_hh, eff_trans_src)) next
        infect_individual(nbr_id, idx, t_inf_nbr, gen_idx + 1L, 5L)
        enqueue(nbr_id, t_inf_nbr)
        n_cumul_infected <- n_cumul_infected + 1L
      }

      for (nbr_id in all_comm_att) {
        if (status[nbr_id] != 1L) next
        t_inf_nbr <- t_out_idx + generation_time_fn(1)
        if (!passes_gates(idx, nbr_id, t_out_idx, p_inf_f_comm, eff_trans_src)) next
        infect_individual(nbr_id, idx, t_inf_nbr, gen_idx + 1L, 6L)
        enqueue(nbr_id, t_inf_nbr)
        n_cumul_infected <- n_cumul_infected + 1L
      }
    }

    status[idx] <- 3L

  } # end main loop

  # ── Assemble output ───────────────────────────────────────────────────────
  infected_mask <- !is.na(time_infection)

  infected_df <- data.frame(
    person_id            = v_person_id[infected_mask],
    hh_id                = v_hh_id[infected_mask],
    cell_id              = v_cell_id[infected_mask],
    hospital_id          = v_hospital_id[infected_mask],
    age_group            = v_age_group[infected_mask],
    is_hcw               = v_is_hcw[infected_mask],
    time_infection       = time_infection[infected_mask],
    time_onset           = time_onset[infected_mask],
    time_hospitalisation = time_hospitalisation[infected_mask],
    time_outcome         = time_outcome[infected_mask],
    hospitalised         = hospitalised[infected_mask],
    outcome_death        = outcome_death[infected_mask],
    outcome_location     = outcome_location[infected_mask],
    funeral_unsafe       = funeral_unsafe_vec[infected_mask],
    funeral_attended_for = funeral_attended_for[infected_mask],
    funeral_role         = funeral_role_vec[infected_mask],
    generation           = generation[infected_mask],
    ancestor_id          = ancestor_id[infected_mask],
    contact_type         = contact_type_vec[infected_mask],
    treated              = treated[infected_mask],
    time_treated         = time_treated[infected_mask],
    quarantined          = quarantined[infected_mask],
    time_quarantined     = time_quarantined[infected_mask],
    traced_via           = traced_via[infected_mask],
    stringsAsFactors     = FALSE
  )
  infected_df <- infected_df[order(infected_df$time_infection,
                                   infected_df$person_id), ]

  list(
    infected         = infected_df,
    stop_reason      = stop_reason,
    n_cumul_infected = n_cumul_infected,
    n_active_at_stop = queue_size()
  )
}
