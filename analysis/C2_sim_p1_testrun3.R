# ==============================================================================
# C2_sim_p1_testrun_dpc.R
# Purpose:
#   Single test run of IBM simulation with DPC-based antiviral efficacy +
#   time-varying antiviral coverage / DPC delay (conflict disruption scenario).
#   Requires sim_prep to be built first via C2_sim_p0_simprep.R
#
# New vs. C2_sim_p1_testrun.R:
#   - dpc_efficacy_data / dpc_efficacy_col: efficacy(dpc) curve replaces the
#     old PK/PD time-since-treatment curve for infection-prevention effect.
#   - PrEP always uses the curve's peak efficacy (dpc_at_peak), regardless of
#     calendar time.
#   - logistical_delay_fn can now be a time-varying function(t) (class
#     "time_varying_fn", built with make_time_varying()) representing the
#     onset->treatment delay (DPC) at outbreak day t. Falls back to the
#     original stochastic function(n) behaviour if not time-varying.
#   - prob_treat_given_trace_household / _close / _physical and
#     prob_treat_self can now also be scalar OR function(t) — i.e. antiviral
#     coverage (분배율) can fall over the outbreak.
# ==============================================================================

source("function/COD_IBM_ebola_sim_obv.R")
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(here)

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

# DPC 효과곡선 컬럼 선택: "efficacy" / "efficacy_lo" / "efficacy_hi" /
#                          "eighty_efficacy_lo" / "eighty_efficacy_hi"
dpc_efficacy_col <- "efficacy"

# ── ① 약물 배포 범위 (셋 중 하나 선택) ──────────────────────────────────────
#   "household"           : household만
#   "household_close"     : household + 1-ring close (넓은 레이어)
#   "household_physical"  : household + 1-ring physical (좁은 레이어)
antiviral_scope <- "household_close"

# ── ② 커버리지(분배율): 상수 또는 time-varying ─────────────────────────────
coverage_mode     <- "tv"      # "constant" 또는 "tv"
coverage_constant <- 0.5       # coverage_mode == "constant"일 때 쓸 값 [0,1]

# ── ③ DPC(distribute_delay): 상수 또는 time-varying ────────────────────────
dpc_mode     <- "tv"           # "constant" 또는 "tv"
dpc_constant <- 3              # dpc_mode == "constant"일 때 쓸 값 (days)

OUT_BASE <- here("outputs", "simulation", "conflict_dpc_max7")
dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# [Section 1] Load data
# ==============================================================================

message(sprintf("=== Loading: %s ===", case_tag))

nodes    <- readRDS(file.path(network_dir, sprintf("%s_nodes.rds",      case_tag)))
sim_prep <- readRDS(file.path(network_dir, sprintf("%s_sim_prep.rds",   case_tag)))

mats_path <- "output/MPMmat/DRC_network_input_matrices.rds"
mats_p6   <- readRDS(mats_path)
close_only_home <- mats_p6$close_only_home  # 16x16
phys_only_home  <- mats_p6$phys_only_home   # 16x16

# DPC -> efficacy 곡선 (예: DPC_fixed_efficacy_varied_d50.rds)
dpc_efficacy_data <- readRDS("data_processed/DPC_fixed_efficacy_varied_d50.rds")

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
# [Section 2] Natural history parameters (unchanged)
# ==============================================================================

make_gamma_fn <- function(mean, sd) {
  shape <- (mean / sd)^2; rate <- mean / sd^2
  function(n) if (n <= 0L) numeric(0) else rgamma(n, shape = shape, rate = rate)
}

incubation_period_fn           <- make_gamma_fn(8.5,  4.5)
onset_to_death_fn              <- make_gamma_fn(9.3,  3.0)
onset_to_recovery_fn           <- make_gamma_fn(13.0, 4.0)
hospitalisation_to_death_fn    <- make_gamma_fn(4.5,  2.0)
hospitalisation_to_recovery_fn <- make_gamma_fn(8.0,  2.5)
onset_to_hospitalisation_fn    <- make_gamma_fn(1.0,  0.35)

generation_time_fn <- local({
  shape <- 2.5; rate <- 2.5 / 15.4
  function(n) if (n <= 0L) numeric(0) else rgamma(n, shape = shape, rate = rate)
})

funeral_generation_time_fn <- local({
  shape <- 20; rate <- 10
  function(n) if (n <= 0L) numeric(0) else rgamma(n, shape = shape, rate = rate)
})

# ==============================================================================
# [Section 3] Transmission parameters (unchanged)
# ==============================================================================

p_inf_household_close            <- 0.00
p_inf_household_physical         <- 0.1

p_inf_community_close_daily      <- 0.00
p_inf_community_close_weekly     <- 0.00
p_inf_community_close_monthly    <- 0.00

p_inf_community_physical_daily   <- 0.0015*21
p_inf_community_physical_weekly  <- 0.00152*3
p_inf_community_physical_monthly <- 0.0015*3/4

p_inf_hcw_to_hcw                 <- 0.02
p_inf_patient_to_hcw             <- 0.02

funeral_avg                      <- 15
funeral_k                        <- 0.3
p_unsafe_funeral                 <- 0.50
p_inf_funeral_household          <- 0.25
p_inf_funeral_community          <- 0.25
funeral_unsafe_multiplier        <- 1.0
funeral_safe_multiplier          <- 0.20

prob_hospitalised_genPop         <- 0.30
prob_hospitalised_hcw            <- 0.60
prob_death_comm                  <- 0.70
prob_death_hosp                  <- 0.50

# ==============================================================================
# [Section 4b] Time-varying response parameters (epi response — unchanged)
# ==============================================================================

make_time_varying <- function(times, values) {
  stopifnot(length(times) == length(values), length(times) >= 2L,
            !is.unsorted(times, strictly = FALSE))
  fn <- stats::approxfun(times, values, rule = 2)
  class(fn) <- c("time_varying_fn", "function")
  fn
}

clip01 <- function(x) pmin(pmax(x, 0), 1)

tv_scenario_id  <- "middle_drc_conflict"
tv_csv_path     <- "data_processed/final_six_scenario_values_revised_methodology.csv"

tv_matrix_full  <- read.csv(tv_csv_path, stringsAsFactors = FALSE)
tv_matrix       <- tv_matrix_full[tv_matrix_full$scenario == tv_scenario_id, ]
tv_matrix       <- tv_matrix[order(tv_matrix$relative_day), ]

if (nrow(tv_matrix) == 0L)
  stop(sprintf("Scenario '%s' not found in %s", tv_scenario_id, tv_csv_path))

times_tv <- tv_matrix$relative_day

p_unsafe_hosp_values <- clip01(
  (1 - tv_matrix$prop_etu) * tv_matrix$prob_unsafe_funeral_hosp +
    tv_matrix$prop_etu     * tv_matrix$prob_unsafe_funeral_etu
)

prob_hospitalised_genPop_fn <- make_time_varying(times_tv, clip01(tv_matrix$prob_hosp))
hospitalisation_delay_factor <- make_time_varying(times_tv, pmax(tv_matrix$delay_hosp, 0.01))
p_unsafe_funeral_comm_fn <- make_time_varying(times_tv, clip01(tv_matrix$prob_unsafe_funeral_comm))
p_unsafe_funeral_hosp_fn <- make_time_varying(times_tv, p_unsafe_hosp_values)
prop_etu_fn <- make_time_varying(times_tv, clip01(tv_matrix$prop_etu))
ipc_index_fn <- make_time_varying(times_tv, clip01(tv_matrix$ipc_helper))
ppe_coverage_fn <- make_time_varying(times_tv, clip01(tv_matrix$ipc_helper))

etu_efficacy_baseline     <- 0.90
non_etu_hospital_efficacy <- 0.30
ppe_efficacy_hcw          <- 0.70

# ==============================================================================
# [Section 4c] Conflict-disruption antiviral coverage & DPC curves
# "외부 원인(분쟁 등)으로 분배율이 떨어지고 지연이 생기면?" 시나리오
# ==============================================================================

sdb <- readRDS(here("data_processed", "SDB_communityDeath_blended.rds"))

# ---- Tweak sdb$value in the 150-400 day window only ----
# Shape-preserving x-axis rescale: curve form maintained, trough shifted
# from day 200 to day 325.
idx_150_325 <- sdb$day >= 150 & sdb$day <= 325
idx_325_400 <- sdb$day >  325 & sdb$day <= 400

rescale_sdb_segment <- function(day_out, orig_from, orig_to) {
  t        <- (day_out - day_out[1]) / (day_out[length(day_out)] - day_out[1])
  day_orig <- orig_from + t * (orig_to - orig_from)
  approx(sdb$day, sdb$value, xout = day_orig, rule = 2)$y
}

sdb_tweaked <- sdb$value
sdb_tweaked[idx_150_325] <- rescale_sdb_segment(sdb$day[idx_150_325], 150, 200)
sdb_tweaked[idx_325_400] <- rescale_sdb_segment(sdb$day[idx_325_400], 200, 350)
sdb$value_tweaked <- sdb_tweaked

# ---- Coverage and DPC curves derived from tweaked sdb ----
sdb$coverage_conflict <- sdb$value_tweaked * 80 / max(sdb$value_tweaked)
sdb$dpc_conflict      <- 1 + 6 * (1 - (sdb$value_tweaked / max(sdb$value_tweaked)))

# Find peak coverage day restricted to day < 200
sub      <- sdb[sdb$day < 200, ]
peak_row <- sub[which.max(sub$coverage_conflict), ]
peak_day <- peak_row$day

# Hold DPC at 1 up until peak_day (before conflict disruption begins)
sdb$dpc_conflict[sdb$day <= peak_day] <- 1

message(sprintf("peak_day = %d, peak_coverage = %.2f", peak_day, peak_row$coverage_conflict))

# ---- Flat coverage comparator: same ramp up to peak_day, held constant after ----
sdb$coverage_flat <- sdb$coverage_conflict
sdb$coverage_flat[sdb$day > peak_day] <- peak_row$coverage_conflict

# ---- Build the TV functions actually fed into the simulator ----
# coverage_conflict is on a 0-100 scale -> rescale to probability [0,1]
antiviral_coverage_fn <- make_time_varying(sdb$day, sdb$coverage_conflict / 100)
dpc_delay_fn           <- make_time_varying(sdb$day, sdb$dpc_conflict)

# ==============================================================================
# [Section 4e] ①②③ 토글 -> 실제 ebola_network_sim() 인자로 변환
# ==============================================================================

# ② 커버리지: 상수 또는 time-varying 함수 중 실제로 쓸 값
coverage_arg <- if (coverage_mode == "constant") coverage_constant else antiviral_coverage_fn

# ③ DPC: 상수 또는 time-varying 함수 중 실제로 쓸 값
#    (상수는 함수로 안 감싸도 됨 — resolve_delay()가 순수 숫자도 그대로 받음)
dpc_arg <- if (dpc_mode == "constant") dpc_constant else dpc_delay_fn

# ① 배포 범위: 선택 안 된 레이어는 coverage 0(=치료 없음)
zero_or_coverage <- function(layer_on) if (layer_on) coverage_arg else 0

cov_household <- zero_or_coverage(TRUE)  # household는 세 옵션 모두 공통으로 포함
cov_close     <- zero_or_coverage(antiviral_scope == "household_close")
cov_physical  <- zero_or_coverage(antiviral_scope == "household_physical")

message(sprintf("  Antiviral scope = %s | coverage_mode = %s | dpc_mode = %s",
                antiviral_scope, coverage_mode, dpc_mode))

# ==============================================================================
# [Section 4d] Diagnostic plot — conflict scenario coverage / DPC / response TVs
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
  ppe_coverage          = ppe_coverage_fn(t_seq),
  # 아래 둘은 ②③ 토글에서 실제로 시뮬에 들어간 값(coverage_arg/dpc_arg) 기준
  # — constant 모드면 평평한 선, tv 모드면 곡선이 그대로 보여야 함
  antiviral_coverage    = if (coverage_mode == "constant") rep(coverage_constant, length(t_seq))
  else antiviral_coverage_fn(t_seq),
  dpc_delay             = if (dpc_mode == "constant") rep(dpc_constant, length(t_seq))
  else dpc_delay_fn(t_seq)
)

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
  antiviral_coverage    = "Antiviral coverage (분배율)",
  dpc_delay             = "DPC delay (days)",
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
  labs(title    = sprintf("%s — Time-varying response parameters (conflict scenario)", case_tag),
       subtitle = "Input parameters (blue) | derived hospital/PPE efficacy (red)",
       x = "Day since outbreak start", y = "Value") +
  theme_bw() +
  theme(plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "grey40"),
        strip.text    = element_text(size = 8,  face = "bold"),
        legend.position = "top")

out_tv <- file.path(OUT_BASE, sprintf("%s_tv_params_conflict.png", case_tag))
ggsave(out_tv, plot = p_tv, width = 14, height = 10, dpi = 150)
message(sprintf("Saved: %s", out_tv))

# DPC efficacy curve itself (sanity check)
p_dpc_curve <- ggplot(dpc_efficacy_data, aes(x = dpc, y = .data[[dpc_efficacy_col]])) +
  geom_line(linewidth = 0.9, color = "#993C1D") +
  labs(title = sprintf("DPC efficacy curve (column: %s)", dpc_efficacy_col),
       x = "DPC (days, onset → treatment)", y = "Efficacy") +
  theme_bw()
ggsave(file.path(OUT_BASE, sprintf("%s_dpc_efficacy_curve.png", case_tag)),
       plot = p_dpc_curve, width = 7, height = 5, dpi = 150)

# ==============================================================================
# [Section 5] Run simulation
# ==============================================================================

message("\n=== Running simulation (conflict-disrupted antiviral coverage/DPC) ===")
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

  # Time-varying epi response parameters
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

  # ── Antiviral treatment: DPC-driven efficacy + TV coverage/delay ──────────
  antiviral_start               = 0,                   # antiviral 프로그램이 켜져 있는 시점
  dpc_efficacy_data             = dpc_efficacy_data,    # DPC_fixed_efficacy_varied_d50.rds
  dpc_efficacy_col              = dpc_efficacy_col,

  logistical_delay_fn           = dpc_arg,              # ★ distribute_delay = ③에서 고른 DPC(상수/TV)
  # notify_delay_fn은 지정 안 함 → 함수 내부에서 임시로 입원지연(onset_to_hospitalisation_fn × hospitalisation_delay_factor)과 동일하게 자동 적용됨

  prob_treat_self                  = coverage_arg,  # ②에서 고른 커버리지(상수/TV)
  prob_treat_given_trace_household = cov_household,  # ① 항상 포함
  prob_treat_given_trace_close     = cov_close,      # ① scope=="household_close"일 때만 coverage, 아니면 0
  prob_treat_given_trace_physical  = cov_physical,   # ① scope=="household_physical"일 때만 coverage, 아니면 0
  prob_treat_given_trace_funeral   = 0,                       # 장례는 일단 비활성 (필요시 켜기)

  prob_trace_household           = 1.0,   # 가구 내는 항상 추적된다고 가정 (필요시 조정)
  prob_trace_close               = 0.5,
  prob_trace_physical            = 0.5,
  prob_trace_funeral             = 0,

  # PrEP: 항상 곡선의 최댓값 효과 사용 (dpc_at_peak, 시뮬함수 내부에서 자동 처리)
  hcw_prep_start                = Inf,    # 필요시 활성화
  prob_treat_hcw_prep           = 0,
  prob_treat_hcw_pep            = 0,

  max_time           = max_time,
  max_infected       = max_infected,
  seed               = sim_seed,
  monitoring_console = TRUE
)

elapsed <- round(proc.time()[["elapsed"]] - t_start, 1)
message(sprintf("\n=== Done: %.1f sec ===", elapsed))

# ==============================================================================
# [Section 6] Summary
# ==============================================================================

inf_df  <- result$infected
prev_df <- result$prevented

ctype_labels <- c("0"="Index", "1"="Household", "2"="Comm close",
                  "3"="Comm physical", "4"="Hospital",
                  "5"="Funeral HH", "6"="Funeral comm")

cat(sprintf("\n── Simulation summary (conflict DPC scenario) ──\n"))
cat(sprintf("  Case          : %s\n", case_tag))
cat(sprintf("  Stop reason   : %s\n", result$stop_reason))
cat(sprintf("  Total infected: %d\n", result$n_cumul_infected))
cat(sprintf("  Deaths        : %d (%.1f%%)\n",
            sum(inf_df$outcome_death, na.rm=TRUE),
            100*mean(inf_df$outcome_death, na.rm=TRUE)))
cat(sprintf("  Hospitalised  : %d (%.1f%%)\n",
            sum(inf_df$hospitalised, na.rm=TRUE),
            100*mean(inf_df$hospitalised, na.rm=TRUE)))
cat(sprintf("  Treated       : %d (%.1f%%)\n",
            sum(inf_df$treated, na.rm=TRUE),
            100*mean(inf_df$treated, na.rm=TRUE)))
cat(sprintf("  Median DPC among treated: %.1f days\n",
            median(inf_df$dpc_value[inf_df$treated == 1L], na.rm = TRUE)))
cat(sprintf("  HCWs infected : %d\n", sum(inf_df$is_hcw, na.rm=TRUE)))
cat(sprintf("  Generations   : %d\n", max(inf_df$generation, na.rm=TRUE)))
cat(sprintf("  Duration      : %.0f days\n",
            diff(range(inf_df$time_infection, na.rm=TRUE))))

ct <- table(inf_df$contact_type)
cat("\n  Transmission by route:\n")
for (i in names(ct))
  cat(sprintf("    %-16s: %d (%.1f%%)\n",
              ctype_labels[i], ct[i], 100*ct[i]/sum(ct)))

# ── 약물로 막은(prevented) 노출 이벤트 ──────────────────────────────────────
n_prevented <- nrow(prev_df)
n_attempts  <- n_prevented + nrow(inf_df %>% filter(contact_type > 0))  # index case 제외
cat(sprintf("\n  Prevented exposures (약으로 막음): %d\n", n_prevented))
cat(sprintf("  Prevention rate (전체 노출시도 중): %.1f%%\n",
            100 * n_prevented / max(n_attempts, 1)))
cat(sprintf("  Median DPC among prevented: %.1f days\n",
            median(prev_df$dpc_value, na.rm = TRUE)))
cat(sprintf("  Mean eff_inf used (prevented): %.2f\n",
            mean(prev_df$eff_inf_used, na.rm = TRUE)))

pt <- table(prev_df$contact_type)
cat("\n  Prevented by route:\n")
for (i in names(pt))
  cat(sprintf("    %-16s: %d (%.1f%%)\n",
              ctype_labels[i], pt[i], 100*pt[i]/sum(pt)))

# ==============================================================================
# [Section 7] Plots
# ==============================================================================

inf_df$day_inf     <- floor(inf_df$time_infection)
inf_df$ctype_label <- ctype_labels[as.character(inf_df$contact_type)]

epi_curve <- inf_df %>%
  group_by(day_inf) %>%
  summarise(n_new = n(), .groups = "drop")

p_epi <- ggplot(epi_curve, aes(x = day_inf, y = n_new)) +
  geom_col(fill = "tomato", alpha = 0.8) +
  labs(title    = sprintf("%s — Epidemic curve (conflict DPC scenario)", case_tag),
       subtitle = sprintf("Total: %d | Stop: %s",
                          result$n_cumul_infected, result$stop_reason),
       x = "Day", y = "New infections") +
  theme_bw() +
  theme(plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9,  color = "grey40"))

prev_df$day_exposure <- floor(prev_df$time_exposure)
prev_df$ctype_label  <- ctype_labels[as.character(prev_df$contact_type)]

# DPC 인풋 곡선(라인) + 실제 개인들이 받은 DPC 값(점, 감염/예방 모두 포함)
dpc_observed_pts <- bind_rows(
  inf_df %>% filter(treated == 1L, !is.na(dpc_value)) %>%
    transmute(day = day_inf, dpc_value, outcome = "Infected (despite treatment)"),
  prev_df %>% filter(!is.na(dpc_value)) %>%
    transmute(day = day_exposure, dpc_value, outcome = "Prevented")
)

p_dpc_observed <- ggplot() +
  geom_jitter(data = dpc_observed_pts,
              aes(x = day, y = dpc_value, color = outcome),
              width = 0.3, height = 0, alpha = 0.35, size = 1.2) +
  geom_line(data = tv_diag, aes(x = day, y = dpc_delay),
            linewidth = 1, color = "black") +
  scale_color_manual(values = c("Infected (despite treatment)" = "tomato",
                                "Prevented" = "#185FA5"), name = NULL) +
  labs(title = "DPC delay: input curve (black) vs. individuals' realized DPC (points)",
       subtitle = sprintf("dpc_mode = %s", dpc_mode),
       x = "Day", y = "DPC (days)") +
  coord_cartesian(xlim = c(0, max(dpc_observed_pts$day, na.rm = TRUE))) +
  theme_bw() + theme(legend.position = "top")

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

# ── Prevented (약으로 막은 노출) 플롯 ───────────────────────────────────────

# 경로별 prevented 건수
p_prev_route <- ggplot(prev_df, aes(x = ctype_label, fill = ctype_label)) +
  geom_bar(alpha = 0.8) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "Prevented exposures by route", x = NULL, y = "Count") +
  theme_bw()

# DPC 값에 따른 prevented 분포 (DPC가 짧을수록 막을 확률 높아야 정상)
p_prev_dpc <- ggplot(prev_df, aes(x = dpc_value)) +
  geom_histogram(binwidth = 1, fill = "#185FA5", alpha = 0.8) +
  labs(title = "Prevented exposures by DPC", x = "DPC (days)", y = "Count") +
  theme_bw()

# 시간에 따른 infected vs prevented 비교 (분배율/DPC 악화 시 prevented가 줄어드는지 확인용)
daily_compare <- bind_rows(
  inf_df  %>% filter(contact_type > 0) %>%
    mutate(day = floor(time_infection), outcome = "Infected"),
  prev_df %>% mutate(day = floor(time_exposure), outcome = "Prevented")
) %>%
  group_by(day, outcome) %>%
  summarise(n = n(), .groups = "drop")

p_prev_vs_inf <- ggplot(daily_compare, aes(x = day, y = n, fill = outcome)) +
  geom_col(position = "stack", alpha = 0.85) +
  scale_fill_manual(values = c("Infected" = "tomato", "Prevented" = "#185FA5"), name = NULL) +
  labs(title = "Exposure outcomes over time: infected vs. drug-prevented",
       x = "Day", y = "Count") +
  theme_bw() + theme(legend.position = "top")

p_combined <- (p_epi + p_dpc_observed) / (p_route + p_gen) /
  (p_prev_vs_inf) / (p_prev_route + p_prev_dpc)

out_fig <- file.path(OUT_BASE,
                     sprintf("%s_seed%s_epi_dpc.png", case_tag,
                             ifelse(is.null(index_case_id),"rand",index_case_id)))
ggsave(out_fig, plot = p_combined, width = 12, height = 18, dpi = 150)
cat(sprintf("\nSaved: %s\n", out_fig))

# prevented_df도 별도로 저장해두면 추후 분석에 편함
out_prev <- file.path(OUT_BASE,
                      sprintf("%s_seed%s_prevented.rds", case_tag,
                              ifelse(is.null(index_case_id),"rand",index_case_id)))
saveRDS(prev_df, out_prev)
cat(sprintf("Saved: %s\n", out_prev))

invisible(result)
