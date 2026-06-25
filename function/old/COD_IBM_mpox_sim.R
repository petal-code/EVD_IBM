# ==============================================================================
# COD_IBM_mpox_sim.R
# Purpose:
#   Individual-based mpox transmission model for DRC
#
# Architecture:
#   - Individual state stored as flat vectors for speed
#   - Pre-allocated queue with grow-on-demand
#   - infect_individual() defined inside main function → <<- direct vector writes
#   - Neighbor lookup: data.table with key on 'from' (bidirectional edge list)
#     → only infected individuals' neighbors ever queried (memory efficient)
#   - p_inf varies by contact frequency stratum (daily/weekly/monthly)
#     via 1-(1-p)^n accumulation over 3-week window
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

mpox_network_sim <- function(

  # ── Network (edge lists) ──────────────────────────────────────────────────
  nodes,    # person_id, hh_id, cell_id, hospital_id, age_group, is_hcw
  layer1,   # from, to
  layer2d,  # from, to, is_physical  (daily)
  layer2w,  # from, to, is_physical  (weekly)
  layer2m,  # from, to, is_physical  (monthly)
  layer3h,  # from, to, hospital_id  (HCW-HCW)
  layer3a,  # person_id, hospital_id, hcw_list

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
  p_inf_household,

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
  library(data.table)

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

  # Pre-compute 3-week effective p_inf per stratum
  p_eff_close_daily    <- p_cumul(p_inf_community_close_daily,    21)
  p_eff_close_weekly   <- p_cumul(p_inf_community_close_weekly,   3)
  p_eff_close_monthly  <- p_cumul(p_inf_community_close_monthly,  0.75)
  p_eff_phys_daily     <- p_cumul(p_inf_community_physical_daily,   21)
  p_eff_phys_weekly    <- p_cumul(p_inf_community_physical_weekly,  3)
  p_eff_phys_monthly   <- p_cumul(p_inf_community_physical_monthly, 0.75)

  # PK/PD drug efficacy
  compute_drug_eff <- function(t_current, t_treated, pd_params, pk_params) {
    if (is.null(pk_params) || is.null(pd_params)) return(0)
    dt <- t_current - t_treated
    if (is.na(dt) || dt < 0) return(0)
    conc <- pk_params$dose * exp(-pk_params$ke * dt)
    pd_params$emax * conc / (pd_params$ec50 + conc)
  }

  # ── Person ID lookup ──────────────────────────────────────────────────────
  N          <- nrow(nodes)
  person_ids <- nodes$person_id
  max_pid    <- max(person_ids)

  pid_to_idx <- integer(max_pid + 1L)
  pid_to_idx[person_ids + 1L] <- seq_len(N)
  pid2idx <- function(pid) pid_to_idx[pid + 1L]

  # ── Build bidirectional data.tables with key on 'from' ───────────────────
  # Duplicate each edge in both directions so querying by 'from' covers both ends.
  # This avoids a full scan on 'to' while keeping memory at 2x original edge list.
  message("Preparing edge lookup tables...")

  make_bidir_dt <- function(edges, extra_cols = NULL) {
    cols <- c("from", "to", extra_cols)
    e    <- as.data.table(edges)[, ..cols]
    rev  <- copy(e)
    rev[, c("from","to") := .(to, from)]
    dt <- rbind(e, rev)
    setkey(dt, from)
    dt
  }

  dt_layer1  <- make_bidir_dt(layer1)
  dt_layer2d <- make_bidir_dt(layer2d, "is_physical")
  dt_layer2w <- make_bidir_dt(layer2w, "is_physical")
  dt_layer2m <- make_bidir_dt(layer2m, "is_physical")
  dt_layer3h <- make_bidir_dt(layer3h)

  # Admission lookup: person_id -> HCW indices
  adm_lookup <- list()
  for (k in seq_len(nrow(layer3a))) {
    pid  <- as.character(layer3a$person_id[k])
    hcws <- layer3a$hcw_list[[k]]
    if (!is.null(hcws) && length(hcws) > 0)
      adm_lookup[[pid]] <- pid2idx(hcws)
  }
  message("  Edge tables ready")

  # ── Neighbor accessor functions ───────────────────────────────────────────
  get_nbrs_l1 <- function(pid) {
    e <- dt_layer1[.(pid), to]
    idx <- pid2idx(e[!is.na(e)])
    idx[idx > 0L]
  }

  get_nbrs_l2 <- function(pid, dt, physical_only) {
    e <- dt[.(pid)]
    if (nrow(e) == 0L) return(integer(0))
    e <- if (physical_only) e[is_physical == 1L] else e[is_physical == 0L]
    idx <- pid2idx(e$to)
    idx[!is.na(idx) & idx > 0L]
  }

  get_nbrs_hcw <- function(pid) {
    e <- dt_layer3h[.(pid), to]
    idx <- pid2idx(e[!is.na(e)])
    idx[idx > 0L]
  }

  get_adm_nbrs <- function(pid) {
    v <- adm_lookup[[as.character(pid)]]
    if (is.null(v)) integer(0) else v
  }

  get_all_comm_nbrs <- function(pid) {
    unique(c(
      get_nbrs_l2(pid, dt_layer2d, FALSE), get_nbrs_l2(pid, dt_layer2d, TRUE),
      get_nbrs_l2(pid, dt_layer2w, FALSE), get_nbrs_l2(pid, dt_layer2w, TRUE),
      get_nbrs_l2(pid, dt_layer2m, FALSE), get_nbrs_l2(pid, dt_layer2m, TRUE)
    ))
  }

  get_ring1_all <- function(pid)
    unique(c(get_nbrs_l1(pid), get_all_comm_nbrs(pid)))

  get_ring2_all <- function(pid) {
    r1 <- get_ring1_all(pid)
    if (length(r1) == 0L) return(integer(0))
    r1_pids <- person_ids[r1]
    unique(unlist(lapply(r1_pids, get_ring1_all), use.names = FALSE))
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
    pid2idx(seeding_ids)
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
    pid_idx     <- v_person_id[idx]

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

    # Ring interventions
    if (is.finite(t_ring)) {
      apply_ring(get_nbrs_l1(pid_idx),
                 prob_trace_household,
                 prob_treat_given_trace_household,
                 prob_quarantine_given_trace_household,
                 "household", t_ring, logistical_delay)

      apply_ring(unique(c(
        get_nbrs_l2(pid_idx, dt_layer2d, FALSE),
        get_nbrs_l2(pid_idx, dt_layer2w, FALSE),
        get_nbrs_l2(pid_idx, dt_layer2m, FALSE))),
        prob_trace_close,
        prob_treat_given_trace_close,
        prob_quarantine_given_trace_close,
        "close", t_ring, logistical_delay)

      apply_ring(unique(c(
        get_nbrs_l2(pid_idx, dt_layer2d, TRUE),
        get_nbrs_l2(pid_idx, dt_layer2w, TRUE),
        get_nbrs_l2(pid_idx, dt_layer2m, TRUE))),
        prob_trace_physical,
        prob_treat_given_trace_physical,
        prob_quarantine_given_trace_physical,
        "physical", t_ring, logistical_delay)
    }

    # =========================================================================
    # Phase 1: infection → hospitalisation (or outcome)
    # =========================================================================
    t_phase1_end <- if (!is.na(t_hosp_idx)) t_hosp_idx else t_out_idx

    phase1_groups <- list(
      list(nbrs = get_nbrs_l1(pid_idx),
           p_inf = p_inf_household,       ctype = 1L),
      list(nbrs = get_nbrs_l2(pid_idx, dt_layer2d, FALSE),
           p_inf = p_eff_close_daily,     ctype = 2L),
      list(nbrs = get_nbrs_l2(pid_idx, dt_layer2w, FALSE),
           p_inf = p_eff_close_weekly,    ctype = 2L),
      list(nbrs = get_nbrs_l2(pid_idx, dt_layer2m, FALSE),
           p_inf = p_eff_close_monthly,   ctype = 2L),
      list(nbrs = get_nbrs_l2(pid_idx, dt_layer2d, TRUE),
           p_inf = p_eff_phys_daily,      ctype = 3L),
      list(nbrs = get_nbrs_l2(pid_idx, dt_layer2w, TRUE),
           p_inf = p_eff_phys_weekly,     ctype = 3L),
      list(nbrs = get_nbrs_l2(pid_idx, dt_layer2m, TRUE),
           p_inf = p_eff_phys_monthly,    ctype = 3L)
    )

    if (is_hcw_idx && rbinom(1, 1, prob_hospital_cond_hcw_preAdm) == 1L) {
      hcw_now <- get_nbrs_hcw(pid_idx)
      if (length(hcw_now) > 0)
        phase1_groups[[length(phase1_groups) + 1L]] <- list(
          nbrs  = hcw_now,
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
        phase2_nbrs <- get_nbrs_hcw(pid_idx)
        p_inf_p2    <- if (length(phase2_nbrs) > 0)
          p_inf_hcw_to_hcw / length(phase2_nbrs) else 0
      } else {
        phase2_nbrs <- get_adm_nbrs(pid_idx)
        p_inf_p2    <- if (length(phase2_nbrs) > 0)
          p_inf_patient_to_hcw / length(phase2_nbrs) else 0
      }

      for (nbr_id in phase2_nbrs) {
        if (status[nbr_id] != 1L) next
        t_inf_nbr <- t_hosp_idx + generation_time_fn(1)
        if (t_inf_nbr > t_out_idx) next
        if (!passes_gates(idx, nbr_id, t_inf_nbr, p_inf_p2, eff_trans_src)) next

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
    # =========================================================================
    if (isTRUE(outcome_death[idx])) {

      funedgenum   <- rnbinom(1, size = funeral_k, mu = funeral_avg)
      is_unsafe    <- isTRUE(funeral_unsafe_vec[idx])
      funeral_mult <- if (is_unsafe) funeral_unsafe_multiplier else funeral_safe_multiplier

      hh_att    <- get_nbrs_l1(pid_idx)
      remaining <- max(0L, funedgenum - length(hh_att))

      comm_1ring <- setdiff(get_all_comm_nbrs(pid_idx), hh_att)
      comm_att   <- integer(0)
      if (remaining > 0L && length(comm_1ring) > 0L) {
        n_draw    <- min(remaining, length(comm_1ring))
        comm_att  <- sample(comm_1ring, n_draw, replace = FALSE)
        remaining <- remaining - n_draw
      }

      ring2_att <- integer(0)
      if (remaining > 0L) {
        ring2_all <- setdiff(get_ring2_all(pid_idx), c(hh_att, comm_att, idx))
        if (length(ring2_all) > 0L) {
          n_draw    <- min(remaining, length(ring2_all))
          ring2_att <- sample(ring2_all, n_draw, replace = FALSE)
        }
      }

      all_comm_att <- unique(c(comm_att, ring2_att))

      for (att_id in c(hh_att, all_comm_att)) {
        if (is.na(funeral_attended_for[att_id])) {
          funeral_attended_for[att_id] <<- v_person_id[idx]
          funeral_role_vec[att_id]     <<- if (att_id %in% hh_att)
            "household" else "community"
        }
      }

      if (is.finite(t_ring)) {
        apply_ring(c(hh_att, all_comm_att),
                   prob_trace_funeral,
                   prob_treat_given_trace_funeral,
                   prob_quarantine_given_trace_funeral,
                   "funeral", t_ring, logistical_delay)
      }

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
