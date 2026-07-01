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
#   - p_inf supplied directly as per-edge probability (no frequency accumulation)
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
  funeral_generation_time_fn    = NULL,  # separate Gamma for funeral infections (fiber: shape=20, rate=10, mean=2d)

  # ── Disease severity ──────────────────────────────────────────────────────
  prob_symptomatic             = 1.0,
  prob_hospitalised_genPop,
  prob_hospitalised_hcw,
  prob_death_comm,
  prob_death_hosp,

  # ── Transmission probabilities ────────────────────────────────────────────
  # Per-edge probability supplied directly (frequency-dependent, not accumulated).
  # daily/weekly/monthly edges are separate network layers; p_inf is applied
  # once per edge per transmission attempt.
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

  # ── Time-varying response parameters (scalar or function(t)) ─────────────
  # Following fiber framework: each can be a fixed scalar or function(t)
  prob_hospitalised_genPop_fn   = NULL,  # overrides prob_hospitalised_genPop if supplied
  hospitalisation_delay_factor  = 1.0,   # scalar or function(t_onset): multiplier on onset_to_hosp draws
  p_unsafe_funeral_comm_fn      = NULL,  # overrides p_unsafe_funeral for community deaths
  p_unsafe_funeral_hosp_fn      = NULL,  # probability of unsafe funeral after hospital death
  prop_etu_fn                   = NULL,  # scalar or function(t): proportion of hosp cases in ETU
  ipc_index_fn                  = NULL,  # scalar or function(t): latent IPC/PPE response maturity [0,1]
  etu_efficacy_baseline         = 0.80,  # fixed: baseline ETU transmission-blocking efficacy
  non_etu_hospital_efficacy     = 0.30,  # fixed: non-ETU hospital transmission reduction

  # ── HCW-specific ──────────────────────────────────────────────────────────
  ppe_coverage_fn               = 0,     # scalar or function(t): PPE coverage among HCWs
  ppe_efficacy_hcw              = 0,     # fixed: PPE efficacy when worn correctly
  prob_hospital_cond_hcw_preAdm = 0.5,

  # ── Treatment ─────────────────────────────────────────────────────────────
  antiviral_start               = Inf,
  drug_eff_trans                = 0,
  drug_eff_inf                  = 0,
  drug_eff_death                = 0,
  pk_params                     = NULL,
  pd_params_trans               = NULL,
  pd_params_inf                 = NULL,
  drug_eff_inf_data             = NULL,  # data.frame(t_since_treatment, eff):
  # 투약 후 경과시간 기준 감염예방효과 곡선
  # (상승→피크→하강). 공급 시 pk_params보다 우선.
  dpc_efficacy_data             = NULL,  # DPC_fixed_efficacy_varied_d50.rds 형태:
  # data.frame(dpc, efficacy, efficacy_lo, efficacy_hi, ...)
  # dpc(onset→treatment 총 지연일수) 기준 효과곡선.
  # 공급 시 drug_eff_inf_data/pk_params보다 우선.
  dpc_efficacy_col              = "efficacy",  # 위 data.frame에서 사용할 컬럼명
  # ("efficacy","efficacy_lo","efficacy_hi",
  #  "eighty_efficacy_lo","eighty_efficacy_hi" 등)

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
  notify_delay_fn               = NULL,  # 상수/function(n)/time_varying_fn:
  # 발현 → 보고/인지(notify)까지 지연.
  # NULL이면 임시로 입원지연과 동일하게 계산
  # (onset_to_hospitalisation_fn(1) * hospitalisation_delay_factor(t))
  logistical_delay_fn           = function(n) rep(0, n),
  # = distribute_delay: 보고 → 실제 투약/격리까지 지연
  # (DPC 인풋 데이터는 여기로 들어감)

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
  # fiber-style resolve: scalar or function(t), vectorised over t
  resolve_tv <- function(param, t) {
    if (is.function(param)) param(t) else rep(param, length(t))
  }

  # Resolve and clamp to [0, 1]
  resolve_prob <- function(param, t) {
    val <- resolve_tv(param, t)
    pmin(pmax(val, 0), 1)
  }

  # Compute effective hospital quarantine efficacy at time t (fiber formula):
  #   etu_eff(t) = etu_efficacy_baseline + (1 - etu_efficacy_baseline) * ipc_index(t)
  #   hosp_quar_eff(t) = prop_etu(t) * etu_eff(t) + (1 - prop_etu(t)) * ipc_index(t)
  # Falls back to non_etu_hospital_efficacy if prop_etu_fn/ipc_index_fn not supplied.
  compute_hospital_quarantine_eff <- function(t) {
    if (!is.null(prop_etu_fn) && !is.null(ipc_index_fn)) {
      p_etu_t   <- resolve_prob(prop_etu_fn,   t)
      ipc_t     <- resolve_prob(ipc_index_fn,  t)
      etu_eff_t <- etu_efficacy_baseline + (1 - etu_efficacy_baseline) * ipc_t
      p_etu_t * etu_eff_t + (1 - p_etu_t) * ipc_t
    } else {
      # Fallback: weighted average of ETU and non-ETU efficacy with fixed prop_etu
      p_etu_t <- if (!is.null(prop_etu_fn)) resolve_prob(prop_etu_fn, t) else 0
      p_etu_t * etu_efficacy_baseline + (1 - p_etu_t) * non_etu_hospital_efficacy
    }
  }

  # Effective PPE efficacy at time t:
  #   eff_ppe(t) = ppe_coverage(t) * ppe_efficacy_hcw
  compute_ppe_eff <- function(t) {
    cov_t <- resolve_prob(ppe_coverage_fn, t)
    cov_t * ppe_efficacy_hcw
  }

  # Community p_inf used as-is (no frequency accumulation)
  # User supplies the intended per-edge transmission probability directly
  p_eff_close_daily    <- p_inf_community_close_daily
  p_eff_close_weekly   <- p_inf_community_close_weekly
  p_eff_close_monthly  <- p_inf_community_close_monthly
  p_eff_phys_daily     <- p_inf_community_physical_daily
  p_eff_phys_weekly    <- p_inf_community_physical_weekly
  p_eff_phys_monthly   <- p_inf_community_physical_monthly

  # PK/PD 드러그 효율 (1구획 PK + Emax PD)
  compute_drug_eff <- function(t_current, t_treated, pd_params, pk_params) {
    if (is.null(pk_params) || is.null(pd_params)) return(0)
    dt <- t_current - t_treated
    if (is.na(dt) || dt < 0) return(0)
    conc <- pk_params$dose * exp(-pk_params$ke * dt)
    pd_params$emax * conc / (pd_params$ec50 + conc)
  }

  # 데이터 기반 감염예방효과 곡선 (투약 후 경과시간 dt 기준, 상승→피크→하강 가능)
  # dt < 0 (투약 전 노출)이면 효과 0. dt가 데이터 범위를 넘으면 마지막 값으로 고정.
  drug_eff_inf_fn <- NULL
  if (!is.null(drug_eff_inf_data)) {
    drug_eff_inf_fn <- approxfun(
      x    = drug_eff_inf_data$t_since_treatment,
      y    = drug_eff_inf_data$eff,
      rule = 2
    )
  }
  compute_drug_eff_data <- function(eff_fn, t_current, t_treated) {
    if (is.null(eff_fn) || is.na(t_treated)) return(0)
    dt <- t_current - t_treated
    if (is.na(dt) || dt < 0) return(0)
    eff_fn(dt)
  }

  # DPC(onset → treatment 총 지연일수) 기준 효과 곡선
  # 예: DPC_fixed_efficacy_varied_d50.rds — dpc_efficacy_col로 컬럼 선택
  dpc_eff_fn <- NULL
  dpc_at_peak <- 0  # PrEP(사전투약)에 쓸 "최댓값 효과" 지점 (기본값: dpc=0)
  if (!is.null(dpc_efficacy_data)) {
    if (!dpc_efficacy_col %in% names(dpc_efficacy_data))
      stop(sprintf("dpc_efficacy_col '%s' not found in dpc_efficacy_data", dpc_efficacy_col))
    dpc_eff_fn <- approxfun(
      x    = dpc_efficacy_data$dpc,
      y    = dpc_efficacy_data[[dpc_efficacy_col]],
      rule = 2
    )
    # 곡선이 비단조(non-monotonic)일 수 있으므로, dpc=0이 아니라 실제 최댓값이 있는
    # dpc 지점을 찾아서 PrEP에 사용 (항상 곡선의 진짜 최댓값을 보장)
    dpc_at_peak <- dpc_efficacy_data$dpc[which.max(dpc_efficacy_data[[dpc_efficacy_col]])]
  }

  # logistical_delay_fn이 make_time_varying()으로 만든 시간가변 함수(class "time_varying_fn")면
  # 캘린더 시간 t에서 결정론적 지연값을 반환, 아니면 기존처럼 확률적 draw(n=1)
  # fn: time_varying_fn(클래스 태그) -> 캘린더시간 t에서 결정론적 값
  #     일반 function(n)         -> 기존처럼 확률적 draw(n=1)
  #     순수 숫자(상수)           -> 그대로 상수 반환 (DPC를 상수로 고정하고 싶을 때)
  resolve_delay <- function(fn, t) {
    if (inherits(fn, "time_varying_fn")) fn(t)
    else if (is.function(fn)) fn(1)
    else fn
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
  dpc_value            <- rep(NA_real_,     N)  # 각 개인이 실제로 받은 onset→treatment 지연(dpc)

  # ── 약물로 막힌(prevented) 노출 이벤트 로그 (queue와 동일한 grow-on-demand 버퍼) ──
  prevented_cap        <- max(1000L, N %/% 100L)
  prevented_n           <- 0L
  prevented_target_id   <- integer(prevented_cap)   # 노출당했지만 감염 안 된 사람
  prevented_source_id   <- integer(prevented_cap)   # 감염원
  prevented_t_exposure  <- numeric(prevented_cap)
  prevented_dpc         <- numeric(prevented_cap)
  prevented_eff         <- numeric(prevented_cap)
  prevented_ctype       <- integer(prevented_cap)

  log_prevented <- function(target_id, source_id, t_exposure, dpc, eff, ctype) {
    prevented_n <<- prevented_n + 1L
    if (prevented_n > prevented_cap) {
      prevented_cap         <<- prevented_cap * 2L
      prevented_target_id   <<- c(prevented_target_id,  integer(prevented_cap %/% 2L))
      prevented_source_id   <<- c(prevented_source_id,  integer(prevented_cap %/% 2L))
      prevented_t_exposure  <<- c(prevented_t_exposure, numeric(prevented_cap %/% 2L))
      prevented_dpc         <<- c(prevented_dpc,        numeric(prevented_cap %/% 2L))
      prevented_eff         <<- c(prevented_eff,        numeric(prevented_cap %/% 2L))
      prevented_ctype       <<- c(prevented_ctype,      integer(prevented_cap %/% 2L))
    }
    prevented_target_id[prevented_n]  <<- target_id
    prevented_source_id[prevented_n]  <<- source_id
    prevented_t_exposure[prevented_n] <<- t_exposure
    prevented_dpc[prevented_n]        <<- dpc
    prevented_eff[prevented_n]        <<- eff
    prevented_ctype[prevented_n]      <<- ctype
  }
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

    # Hospitalisation probability: time-varying override takes precedence
    p_hosp_base <- if (isTRUE(v_is_hcw[idx])) prob_hospitalised_hcw else prob_hospitalised_genPop
    p_hosp_tv   <- if (!isTRUE(v_is_hcw[idx]) && !is.null(prob_hospitalised_genPop_fn))
      prob_hospitalised_genPop_fn else p_hosp_base
    p_hosp_t    <- resolve_prob(p_hosp_tv, t_onset)
    is_hosp     <- rbinom(1, 1, p_hosp_t) == 1L

    # Hospitalisation delay: raw draw × time-varying delay factor (fiber pattern)
    t_hosp <- if (is_hosp) {
      delay_factor_t <- resolve_tv(hospitalisation_delay_factor, t_onset)
      t_onset + onset_to_hospitalisation_fn(1) * delay_factor_t
    } else NA_real_

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

    # Funeral safety: community vs hospital death, time-varying (fiber pattern)
    funeral_uns <- NA
    if (is_death) {
      if (outcome_loc == "hospital" && !is.null(p_unsafe_funeral_hosp_fn)) {
        p_unsafe_t <- resolve_prob(p_unsafe_funeral_hosp_fn, t_outcome)
      } else if (outcome_loc == "community" && !is.null(p_unsafe_funeral_comm_fn)) {
        p_unsafe_t <- resolve_prob(p_unsafe_funeral_comm_fn, t_outcome)
      } else {
        p_unsafe_t <- resolve_prob(p_unsafe_funeral, t_outcome)
      }
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
    dpc_value[prep_idx]    <- dpc_at_peak  # 사전(prophylactic) 투약: 항상 최댓값 효과 사용
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

  # Resolve funeral generation time function (fallback to community Tg if not supplied)
  funeral_Tg_fn <- if (!is.null(funeral_generation_time_fn))
    funeral_generation_time_fn else generation_time_fn
  passes_gates <- function(idx, nbr_id, t_inf_nbr, p_inf, eff_trans_src, ctype = NA_integer_) {
    if (rbinom(1, 1, p_inf) == 0L) return(FALSE)
    if (eff_trans_src > 0 && rbinom(1, 1, eff_trans_src) == 1L) return(FALSE)
    if (quarantined[idx] == 1L && !is.na(time_quarantined[idx]) &&
        t_inf_nbr > time_quarantined[idx] &&
        rbinom(1, 1, quarantine_efficacy) == 1L) return(FALSE)
    eff_inf_tgt <- if (treated[nbr_id] == 1L) {
      if (!is.null(dpc_eff_fn) && !is.na(dpc_value[nbr_id]))
        dpc_eff_fn(dpc_value[nbr_id])
      else if (!is.null(drug_eff_inf_fn))
        compute_drug_eff_data(drug_eff_inf_fn, t_inf_nbr, time_treated[nbr_id])
      else if (!is.null(pk_params))
        compute_drug_eff(t_inf_nbr, time_treated[nbr_id], pd_params_inf, pk_params)
      else drug_eff_inf
    } else 0
    if (eff_inf_tgt > 0 && rbinom(1, 1, eff_inf_tgt) == 1L) {
      log_prevented(nbr_id, idx, t_inf_nbr, dpc_value[nbr_id], eff_inf_tgt, ctype)
      return(FALSE)
    }
    TRUE
  }

  # ── Ring intervention helper ──────────────────────────────────────────────
  apply_ring <- function(nbr_ids, trace_prob, treat_prob, quar_prob,
                         trace_label, t_ring, delay) {
    # treat_prob: 상수 또는 function(t) (antiviral coverage가 시간에 따라 변하는 경우)
    treat_prob_t <- resolve_prob(treat_prob, t_ring)
    for (nbr_id in nbr_ids) {
      if (!is.na(traced_via[nbr_id])) next
      if (rbinom(1, 1, trace_prob) == 0L) next
      traced_via[nbr_id] <<- trace_label
      if (treat_prob_t > 0 && treated[nbr_id] == 0L &&
          is.finite(t_ring) && t_ring >= antiviral_start &&
          rbinom(1, 1, treat_prob_t) == 1L) {
        treated[nbr_id]      <<- 1L
        time_treated[nbr_id] <<- t_ring + delay
        dpc_value[nbr_id]    <<- delay  # 보고→투약(distribute_delay) 그 자체가 DPC
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
        message(sprintf("T = %4f  [event %6d] active=%4d | total=%6d",
                        t_earliest, event_count, queue_size(), n_cumul_infected))
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
      dpc_value[prep_idx]    <- dpc_at_peak  # 사전(prophylactic) 투약: 항상 최댓값 효과 사용
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

    # notify_delay: 발현 → 보고/인지. notify_delay_fn 없으면 임시로 입원지연과 동일하게.
    notify_delay <- if (!is.null(notify_delay_fn)) {
      resolve_delay(notify_delay_fn, t_onset_idx)
    } else {
      onset_to_hospitalisation_fn(1) * resolve_tv(hospitalisation_delay_factor, t_onset_idx)
    }
    # distribute_delay: 보고 → 실제 투약/격리. DPC 인풋 데이터가 여기로 들어옴.
    distribute_delay <- resolve_delay(logistical_delay_fn, t_onset_idx)

    # Self-treatment (본인이 이미 인지한 상태이므로 notify 단계 없이 distribute_delay만 적용)
    prob_treat_self_t <- resolve_prob(prob_treat_self, t_onset_idx)
    if (!is.na(t_onset_idx) && t_onset_idx >= antiviral_start &&
        treated[idx] == 0L && rbinom(1, 1, prob_treat_self_t) == 1L) {
      treated[idx]      <- 1L
      time_treated[idx] <- t_onset_idx + distribute_delay
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
      t_onset_idx + notify_delay else Inf

    # Ring interventions (O(1) lookup from sim_prep)
    if (is.finite(t_ring)) {
      apply_ring(hh_nbrs[[idx]],
                 prob_trace_household,
                 prob_treat_given_trace_household,
                 prob_quarantine_given_trace_household,
                 "household", t_ring, distribute_delay)

      apply_ring(unique(c(comm_close_daily_nbrs[[idx]],
                          comm_close_weekly_nbrs[[idx]],
                          comm_close_monthly_nbrs[[idx]])),
                 prob_trace_close,
                 prob_treat_given_trace_close,
                 prob_quarantine_given_trace_close,
                 "close", t_ring, distribute_delay)

      apply_ring(unique(c(comm_phys_daily_nbrs[[idx]],
                          comm_phys_weekly_nbrs[[idx]],
                          comm_phys_monthly_nbrs[[idx]])),
                 prob_trace_physical,
                 prob_treat_given_trace_physical,
                 prob_quarantine_given_trace_physical,
                 "physical", t_ring, distribute_delay)
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
      if (!passes_gates(idx, nbr_id, t_inf_nbr, p_eff_hh, eff_trans_src, ctype = 1L)) next
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
    # PPE efficacy resolved at transmission time: ppe_coverage(t) * ppe_efficacy_hcw
    if (is_hcw_idx && length(hcw_nbrs[[idx]]) > 0 &&
        rbinom(1, 1, prob_hospital_cond_hcw_preAdm) == 1L) {
      ppe_eff_t <- compute_ppe_eff(t_inf_idx)
      phase1_groups[[length(phase1_groups) + 1L]] <- list(
        nbrs  = hcw_nbrs[[idx]],
        p_inf = p_inf_hcw_to_hcw * (1 - ppe_eff_t),
        ctype = 4L)
    }

    for (grp in phase1_groups) {
      if (length(grp$nbrs) == 0L || grp$p_inf <= 0) next
      for (nbr_id in grp$nbrs) {
        if (status[nbr_id] != 1L) next
        t_inf_nbr <- t_inf_idx + generation_time_fn(1)
        if (t_inf_nbr > t_phase1_end) next
        if (!passes_gates(idx, nbr_id, t_inf_nbr, grp$p_inf, eff_trans_src, ctype = grp$ctype)) next
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
        # Hospital quarantine thinning (fiber pattern):
        # ETU isolation + IPC maturity reduce post-admission transmission
        hosp_quar_eff_t <- compute_hospital_quarantine_eff(t_inf_nbr)
        if (hosp_quar_eff_t > 0 && rbinom(1, 1, hosp_quar_eff_t) == 1L) next
        if (!passes_gates(idx, nbr_id, t_inf_nbr, p_inf_p2, eff_trans_src, ctype = 4L)) next

        # HCW PEP on exposure
        if (isTRUE(v_is_hcw[nbr_id]) && prob_treat_hcw_pep > 0 &&
            treated[nbr_id] == 0L &&
            rbinom(1, 1, prob_treat_hcw_pep) == 1L) {
          treated[nbr_id]      <- 1L
          time_treated[nbr_id] <- t_inf_nbr + distribute_delay
          dpc_value[nbr_id]    <- distribute_delay  # 노출→PEP 투약까지 지연
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

      # Pool 1: household (sampled, not mandatory)
      hh_pool   <- hh_nbrs[[idx]]
      n_hh_draw <- min(funedgenum, length(hh_pool))
      hh_att    <- if (n_hh_draw > 0L) sample(hh_pool, n_hh_draw, replace = FALSE) else integer(0)
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
                   "funeral", t_ring, distribute_delay)
      }

      # Transmission
      p_inf_f_hh   <- p_inf_funeral_household * funeral_mult
      p_inf_f_comm <- p_inf_funeral_community * funeral_mult

      for (nbr_id in hh_att) {
        if (status[nbr_id] != 1L) next
        t_inf_nbr <- t_out_idx + funeral_Tg_fn(1)
        if (!passes_gates(idx, nbr_id, t_out_idx, p_inf_f_hh, eff_trans_src, ctype = 5L)) next
        infect_individual(nbr_id, idx, t_inf_nbr, gen_idx + 1L, 5L)
        enqueue(nbr_id, t_inf_nbr)
        n_cumul_infected <- n_cumul_infected + 1L
      }

      for (nbr_id in all_comm_att) {
        if (status[nbr_id] != 1L) next
        t_inf_nbr <- t_out_idx + funeral_Tg_fn(1)
        if (!passes_gates(idx, nbr_id, t_out_idx, p_inf_f_comm, eff_trans_src, ctype = 6L)) next
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
    dpc_value            = dpc_value[infected_mask],
    quarantined          = quarantined[infected_mask],
    time_quarantined     = time_quarantined[infected_mask],
    traced_via           = traced_via[infected_mask],
    stringsAsFactors     = FALSE
  )
  infected_df <- infected_df[order(infected_df$time_infection,
                                   infected_df$person_id), ]

  # ── prevented_df: 약물로 막힌 노출 이벤트 (1건당 1행, 같은 사람이 여러 번 나올 수 있음) ──
  prevented_idx <- seq_len(prevented_n)
  prevented_df <- data.frame(
    target_person_id   = v_person_id[prevented_target_id[prevented_idx]],
    source_person_id   = v_person_id[prevented_source_id[prevented_idx]],
    target_is_hcw       = v_is_hcw[prevented_target_id[prevented_idx]],
    contact_type        = prevented_ctype[prevented_idx],
    time_exposure       = prevented_t_exposure[prevented_idx],
    dpc_value            = prevented_dpc[prevented_idx],
    eff_inf_used         = prevented_eff[prevented_idx],
    time_treated         = time_treated[prevented_target_id[prevented_idx]],
    stringsAsFactors     = FALSE
  )
  prevented_df <- prevented_df[order(prevented_df$time_exposure,
                                     prevented_df$target_person_id), ]

  list(
    infected         = infected_df,
    prevented        = prevented_df,
    stop_reason      = stop_reason,
    n_cumul_infected = n_cumul_infected,
    n_active_at_stop = queue_size()
  )
}
