# ============================================================
# Build frequency-stratified 3-week contact matrices for DRC (Congo),
# now generalized to THREE unique-partner models (linear / Bernoulli s=1 /
# Bernoulli s=0.5), and produce visualization figures for each.
#
# v3 changes (this merge):
#   - Merges the matrix-construction script (p6) with the rarefaction /
#     scale-comparison diagnostic code into a single script.
#   - appear_factor is now computed for THREE models instead of one:
#       "linear" : unique = events (appear_factor = 1 for every stratum;
#                  this is what the "fully independent, brand-new partner
#                  every event" assumption reduces to at t = L = 21)
#       "s1"     : Bernoulli, T_eff = T_avg (this is the ORIGINAL v2 model —
#                  kept as the default/backward-compatible close_unique_3wk /
#                  physical_unique_3wk used by p7)
#       "s05"    : Bernoulli, T_eff = max(T_avg * 0.5, 1) — assumes true
#                  revisit interval is half the nominal label (people meet
#                  ~2x more often than "weekly"/"monthly" suggests). Daily
#                  is clamped at T_eff=1 since it's already the fastest
#                  representable rate in this discrete daily model.
#   - close_unique_3wk / physical_unique_3wk (top-level, used by p7) remain
#     the "s1" (current) model, unchanged in name/shape for compatibility.
#   - All three models' full unique matrices are ALSO saved under
#     close_unique_3wk_by_model / physical_unique_3wk_by_model.
#   - Visualization: age-stratified cumulative bar charts (both the
#     close-only 3-segment stack, and the close/physical 6-segment stack)
#     are now produced once PER MODEL (3 models x 2 chart types = 6 figures),
#     all on a SHARED y-axis scale so the three models are visually
#     comparable. The events row (model-independent) is a separate figure.
#
# Earlier v2 design choices (still in effect):
#   1. Population (WorldPop) weighting is NOT used in matrix construction —
#      everything is a per-capita rate matrix. Population belongs to a
#      downstream network/node generation script.
#   2. Home is NOT frequency-stratified (close/physical split only).
#   3. events = raw 3-week contact counts (repeats included); unique =
#      events * appear_factor[stratum], now parametrized by model.
#   4. Physical-contact versions use the LIC-blended physical/close ratio.
#   5. Community (work+school+other) assumed additive/independent.
# ============================================================

library(socialmixr)
library(data.table)
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)

dir.create("figure/MPMmat", recursive = TRUE, showWarnings = FALSE)
dir.create("output/MPMmat", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Step 0: Config
# ------------------------------------------------------------
prem_dir     <- "data/Prem_contact"
prem_country <- "Congo"

age_labels <- c("0-4","5-9","10-14","15-19","20-24","25-29",
                "30-34","35-39","40-44","45-49","50-54","55-59",
                "60-64","65-69","70-74","75+")
age_breaks <- c(seq(0, 75, by = 5), Inf)

locations       <- c("home", "work", "school", "other")   # all Prem settings
split_locations <- c("work", "school", "other")            # frequency-stratified
loc_labels       <- c("Home", "Work", "School", "Other")

freq_strata <- list(
  daily   = 1L,
  weekly  = 2L,
  monthly = c(3L, 4L, 5L)
)

freq_category_labels <- c(
  `1` = "daily", `2` = "weekly", `3` = "monthly",
  `4` = "less_often", `5` = "first_time"
)

avg_revisit_days <- c(daily = 1, weekly = 7, monthly = 28)
window_L <- 21

# ------------------------------------------------------------
# Step 0b: three unique-partner models
#   linear : appear_factor = 1 for every stratum (unique == events)
#   s1     : Bernoulli, T_eff = T_avg (current/default — used by p7)
#   s05    : Bernoulli, T_eff = max(T_avg * 0.5, 1)
# ------------------------------------------------------------
model_defs <- list(
  linear = list(type = "linear"),
  s1     = list(type = "bernoulli", s = 1.0),
  s05    = list(type = "bernoulli", s = 0.5)
)
model_labels <- c(
  linear = "Linear (unique = events)",
  s1     = "Bernoulli, s = 1 (current)",
  s05    = "Bernoulli, s = 0.5"
)

compute_appear_factor <- function(model) {
  if (model$type == "linear") {
    return(c(daily = 1, weekly = 1, monthly = 1))
  }
  sapply(avg_revisit_days, function(T_avg) {
    T_eff <- max(T_avg * model$s, 1)
    p     <- 1 / T_eff
    (1 - (1 - p)^window_L) / (p * window_L)
  })
}

appear_factor_by_model <- lapply(model_defs, compute_appear_factor)

message("== appear_factor by model ==")
for (mn in names(appear_factor_by_model)) {
  af <- appear_factor_by_model[[mn]]
  message(sprintf("  [%s] daily=%.4f weekly=%.4f monthly=%.4f",
                  mn, af["daily"], af["weekly"], af["monthly"]))
}

mousa_lic_physical <- c(
  home      = 0.690,
  school    = 0.755,
  community = 0.620,
  work      = 0.565
)

# ------------------------------------------------------------
# Step 1: Load Prem (Congo) close-contact matrices
# ------------------------------------------------------------
load_prem_matrices <- function(data_dir, country) {
  file_map <- list(
    home   = "MUestimates_home_1.xlsx",
    work   = "MUestimates_work_1.xlsx",
    school = "MUestimates_school_1.xlsx",
    other  = "MUestimates_other_locations_1.xlsx"
  )
  load_one <- function(setting) {
    fpath <- file.path(data_dir, file_map[[setting]])
    if (!file.exists(fpath)) return(NULL)
    raw <- read_excel(fpath, sheet = country, col_names = FALSE)
    raw <- raw[-1, ]
    mat <- matrix(as.numeric(as.matrix(raw)), nrow = 16, ncol = 16)
    (mat + t(mat)) / 2
  }
  mats <- lapply(setNames(names(file_map), names(file_map)), load_one)
  comm_mat <- matrix(0, 16, 16)
  for (s in c("work", "school", "other"))
    if (!is.null(mats[[s]])) comm_mat <- comm_mat + mats[[s]]
  mats$community <- comm_mat
  mats
}

prem_mats <- load_prem_matrices(prem_dir, prem_country)

# ------------------------------------------------------------
# Step 2: Compute frequency proportions from POLYMOD (Mossong et al.)
# ------------------------------------------------------------
data(polymod)
pm_contacts     <- as.data.table(polymod$contacts)
pm_participants <- as.data.table(polymod$participants)

pm_merged <- merge(
  pm_contacts,
  pm_participants[, .(part_id, part_age_exact)],
  by = "part_id", all.x = TRUE
)
pm_merged[, age_group := cut(part_age_exact, breaks = age_breaks,
                             labels = age_labels, right = FALSE)]

prem_to_polymod_loc <- list(
  home   = "cnt_home",
  work   = "cnt_work",
  school = "cnt_school",
  other  = c("cnt_transport", "cnt_leisure", "cnt_otherplace")
)

compute_frequency_proportions_by_age <- function(loc_cols, loc_name) {
  is_loc <- rowSums(pm_merged[, ..loc_cols, drop = FALSE]) > 0
  sub <- pm_merged[is_loc]

  prop_table <- sapply(age_labels, function(ag) {
    rows <- sub[age_group == ag & !is.na(frequency_multi) & frequency_multi %in% 1:5]
    if (nrow(rows) == 0) return(setNames(rep(NA_real_, 5), as.character(1:5)))
    tab <- table(factor(rows$frequency_multi, levels = 1:5))
    as.numeric(tab) / sum(tab)
  })
  rownames(prop_table) <- freq_category_labels[as.character(1:5)]
  prop_df <- as.data.frame(t(prop_table))
  prop_df$age_group <- age_labels
  prop_df$location  <- loc_name
  prop_df[, c("age_group", "location", freq_category_labels[as.character(1:5)])]
}

frequency_proportion_diagnostic <- do.call(
  rbind,
  lapply(names(prem_to_polymod_loc), function(loc) {
    compute_frequency_proportions_by_age(prem_to_polymod_loc[[loc]], loc)
  })
)

message("\n== Frequency category proportions by age group and location ==")
for (loc in locations) {
  message(sprintf("\n-- %s --", loc))
  sub_diag <- frequency_proportion_diagnostic[
    frequency_proportion_diagnostic$location == loc,
    c("age_group", "daily", "weekly", "monthly", "less_often", "first_time")
  ]
  sub_diag[, -1] <- round(sub_diag[, -1], 3)
  print(sub_diag, row.names = FALSE)
}

compute_stratum_proportions_by_age <- function(loc_cols) {
  is_loc <- rowSums(pm_merged[, ..loc_cols, drop = FALSE]) > 0
  sub <- pm_merged[is_loc]

  lapply(freq_strata, function(codes) {
    sapply(age_labels, function(ag) {
      rows <- sub[age_group == ag & !is.na(frequency_multi) & frequency_multi %in% 1:5]
      if (nrow(rows) == 0) return(NA_real_)
      mean(rows$frequency_multi %in% codes)
    })
  })
}

stratum_props <- lapply(prem_to_polymod_loc[split_locations], compute_stratum_proportions_by_age)

for (loc in names(stratum_props)) {
  for (st in names(freq_strata)) {
    v <- stratum_props[[loc]][[st]]
    na_idx <- is.na(v)
    if (any(na_idx)) stratum_props[[loc]][[st]][na_idx] <- mean(v, na.rm = TRUE)
  }
}

# ------------------------------------------------------------
# Step 3: Build 3-week CONTACT EVENT matrices (raw, repeats included),
#         then split by frequency stratum. Population-free; model-independent.
# ------------------------------------------------------------
total_events_3wk <- setNames(
  lapply(split_locations, function(loc) prem_mats[[loc]] * 21),
  split_locations
)

message("\n== Total 3-week event rowSums (= prem_mats rowSums * 21) ==")
for (loc in split_locations) {
  message(sprintf("  %-6s: %.2f - %.2f", loc,
                  min(rowSums(total_events_3wk[[loc]])),
                  max(rowSums(total_events_3wk[[loc]]))))
}

apply_stratum_split <- function(total_events_mat, prop_vec) {
  M_row <- total_events_mat * matrix(prop_vec, nrow = 16, ncol = 16, byrow = FALSE)
  M_col <- total_events_mat * matrix(prop_vec, nrow = 16, ncol = 16, byrow = TRUE)
  (M_row + M_col) / 2
}

close_events_3wk <- setNames(vector("list", length(freq_strata)), names(freq_strata))
for (st in names(freq_strata)) {
  close_events_3wk[[st]] <- setNames(
    lapply(split_locations, function(loc)
      apply_stratum_split(total_events_3wk[[loc]], stratum_props[[loc]][[st]])),
    split_locations
  )
  close_events_3wk[[st]]$community <- Reduce(`+`, close_events_3wk[[st]][split_locations])
}

message("\n== Stratum-split event recovery check (should equal total_events_3wk) ==")
for (loc in split_locations) {
  recovered <- close_events_3wk$daily[[loc]] + close_events_3wk$weekly[[loc]] + close_events_3wk$monthly[[loc]]
  diff_range <- range(recovered - total_events_3wk[[loc]])
  message(sprintf("  %-6s max abs diff: %.2e", loc, max(abs(diff_range))))
}

# ------------------------------------------------------------
# Step 4: Compute POLYMOD physical/total ratio (16x16, symmetric)
# ------------------------------------------------------------
pm_for_ratio <- copy(polymod$contacts)
pm_for_ratio$cnt_other_combined <- as.integer(
  rowSums(pm_for_ratio[, c("cnt_transport", "cnt_leisure", "cnt_otherplace")],
          na.rm = TRUE) > 0
)
polymod_for_ratio          <- polymod
polymod_for_ratio$contacts <- pm_for_ratio

compute_physical_ratio_single <- function(loc_col) {
  cm_total_raw <- contact_matrix(
    polymod_for_ratio,
    age.limits = age_breaks[-length(age_breaks)],
    filter     = setNames(list(1), loc_col),
    symmetric  = FALSE
  )$matrix
  cm_total <- (cm_total_raw + t(cm_total_raw)) / 2

  cm_phys_raw <- contact_matrix(
    polymod_for_ratio,
    age.limits = age_breaks[-length(age_breaks)],
    filter     = setNames(list(1, 1), c(loc_col, "phys_contact")),
    symmetric  = FALSE
  )$matrix
  cm_phys <- (cm_phys_raw + t(cm_phys_raw)) / 2

  ratio <- cm_phys / cm_total
  ratio[!is.finite(ratio)] <- NA
  list(total = cm_total, physical = cm_phys, ratio = ratio)
}

polymod_loc_cols <- list(
  home   = "cnt_home",
  work   = "cnt_work",
  school = "cnt_school",
  other  = "cnt_other_combined"
)

polymod_ratios <- lapply(polymod_loc_cols, compute_physical_ratio_single)

get_ratio <- function(loc_name) {
  r <- polymod_ratios[[loc_name]]$ratio
  r[is.na(r)] <- 0
  r
}

ratio_home   <- get_ratio("home")
ratio_work   <- get_ratio("work")
ratio_school <- get_ratio("school")
ratio_other  <- get_ratio("other")

# ------------------------------------------------------------
# Step 5: Blend POLYMOD ratio shape with Mousa LIC/LMIC level
# ------------------------------------------------------------
logit     <- function(p) log(p / (1 - p))
inv_logit <- function(x) exp(x) / (1 + exp(x))

rescale_to_lic_logit <- function(ratio_mat, lic_target_mean, loc_name = "",
                                 tol = 1e-7, max_iter = 200) {
  nonzero <- ratio_mat > 0
  p  <- pmin(pmax(ratio_mat[nonzero], 1e-4), 1 - 1e-4)
  lp <- logit(p)

  lo <- -20; hi <- 20
  for (i in seq_len(max_iter)) {
    mid      <- (lo + hi) / 2
    achieved <- mean(inv_logit(lp + mid))
    if (abs(achieved - lic_target_mean) < tol) break
    if (achieved < lic_target_mean) lo <- mid else hi <- mid
  }

  scaled          <- ratio_mat
  scaled[nonzero] <- inv_logit(lp + mid)

  message(sprintf("  [%s] offset=%.4f | target=%.4f | achieved=%.4f (iter=%d)",
                  loc_name, mid, lic_target_mean, achieved, i))
  scaled
}

blended_ratio <- list(
  home   = rescale_to_lic_logit(ratio_home,   mousa_lic_physical["home"],      loc_name = "home"),
  work   = rescale_to_lic_logit(ratio_work,   mousa_lic_physical["work"],      loc_name = "work"),
  school = rescale_to_lic_logit(ratio_school, mousa_lic_physical["school"],    loc_name = "school"),
  other  = rescale_to_lic_logit(ratio_other,  mousa_lic_physical["community"], loc_name = "other")
)

message("\n== Blended physical ratio — achieved vs target ==")
for (loc in names(blended_ratio)) {
  message(sprintf("  %s: target=%.4f | achieved=%.4f",
                  loc, mousa_lic_physical[ifelse(loc == "other", "community", loc)],
                  mean(blended_ratio[[loc]][blended_ratio[[loc]] > 0])))
}

blended_ratio_comm <- (blended_ratio$work + blended_ratio$school + blended_ratio$other) / 3

# ------------------------------------------------------------
# Step 6: Apply blended ratio -> physical EVENT matrices (per stratum)
#         model-independent, same as events
# ------------------------------------------------------------
physical_events_3wk <- setNames(vector("list", length(freq_strata)), names(freq_strata))
for (st in names(freq_strata)) {
  physical_events_3wk[[st]] <- list()
  for (loc in split_locations) {
    raw <- close_events_3wk[[st]][[loc]] * blended_ratio[[loc]]
    physical_events_3wk[[st]][[loc]] <- (raw + t(raw)) / 2
  }
  physical_events_3wk[[st]]$community <- Reduce(`+`, physical_events_3wk[[st]][split_locations])
}

# ------------------------------------------------------------
# Step 7: Convert EVENTS -> UNIQUE PARTNERS, for EACH of the 3 models
# ------------------------------------------------------------
close_unique_3wk_by_model    <- list()
physical_unique_3wk_by_model <- list()

for (mn in names(model_defs)) {
  af <- appear_factor_by_model[[mn]]
  close_unique_3wk_by_model[[mn]] <- setNames(vector("list", length(freq_strata)), names(freq_strata))
  physical_unique_3wk_by_model[[mn]] <- setNames(vector("list", length(freq_strata)), names(freq_strata))
  for (st in names(freq_strata)) {
    close_unique_3wk_by_model[[mn]][[st]]    <- lapply(close_events_3wk[[st]],    function(m) m * af[st])
    physical_unique_3wk_by_model[[mn]][[st]] <- lapply(physical_events_3wk[[st]], function(m) m * af[st])
  }
}

message("\n== 3-week UNIQUE partner rowSums by model x stratum (close, community) ==")
for (mn in names(model_defs)) {
  for (st in names(freq_strata)) {
    mat <- close_unique_3wk_by_model[[mn]][[st]]$community
    message(sprintf("  [%-5s] %-8s: %.2f - %.2f", mn, st,
                    min(rowSums(mat)), max(rowSums(mat))))
  }
}

# Backward-compatible top-level aliases (used by p7 / network build):
# these remain the "s1" (current/default) model, unchanged in shape.
close_unique_3wk    <- close_unique_3wk_by_model$s1
physical_unique_3wk <- physical_unique_3wk_by_model$s1

# ------------------------------------------------------------
# Step 8: Home — close/physical split ONLY, no frequency stratification
#         (model-independent; home was never frequency-stratified)
# ------------------------------------------------------------
close_only_home <- prem_mats$home * (1 - blended_ratio$home)
phys_only_home  <- prem_mats$home * blended_ratio$home

close_only_home <- (close_only_home + t(close_only_home)) / 2
phys_only_home  <- (phys_only_home  + t(phys_only_home))  / 2

message("\n== Household contact matrices (close/physical only, no frequency split) ==")
message(sprintf("  close_only_home rowSum range: %.3f - %.3f",
                min(rowSums(close_only_home)), max(rowSums(close_only_home))))
message(sprintf("  phys_only_home  rowSum range: %.3f - %.3f",
                min(rowSums(phys_only_home)),  max(rowSums(phys_only_home))))
message(sprintf("  sum check (should equal prem_mats$home rowSums): %.3f - %.3f",
                min(rowSums(close_only_home + phys_only_home)),
                max(rowSums(close_only_home + phys_only_home))))

# ------------------------------------------------------------
# Step 9: Community-level stratum allocation probabilities
#         (uses the s1/current model, since this feeds p7's network build)
# ------------------------------------------------------------
prem_unique_community <- close_unique_3wk$daily$community +
  close_unique_3wk$weekly$community +
  close_unique_3wk$monthly$community

message("\n== prem_unique_community rowSum range (s1 model) ==")
message(sprintf("  %.2f - %.2f",
                min(rowSums(prem_unique_community)),
                max(rowSums(prem_unique_community))))

stratum_prob_mat <- sapply(names(freq_strata), function(st) {
  rowSums(close_unique_3wk[[st]]$community) / (rowSums(prem_unique_community) + 1e-12)
})
rownames(stratum_prob_mat) <- age_labels

stratum_check <- rowSums(stratum_prob_mat)
message(sprintf("  Stratum prob sum range: [%.4f, %.4f]",
                min(stratum_check), max(stratum_check)))

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------
saveRDS(
  list(
    # events = raw 3-week contact counts (repeats included), by stratum
    # (model-independent)
    close_events_3wk    = close_events_3wk,
    physical_events_3wk = physical_events_3wk,

    # unique = estimated distinct partners over 3 weeks, by stratum.
    # Top-level close_unique_3wk / physical_unique_3wk = "s1" model,
    # kept for backward compatibility with p7.
    close_unique_3wk    = close_unique_3wk,
    physical_unique_3wk = physical_unique_3wk,

    # ALL THREE models' unique matrices (linear / s1 / s05)
    close_unique_3wk_by_model    = close_unique_3wk_by_model,
    physical_unique_3wk_by_model = physical_unique_3wk_by_model,
    appear_factor_by_model       = appear_factor_by_model,

    # home: close/physical split only, no frequency stratification
    close_only_home = close_only_home,
    phys_only_home  = phys_only_home,

    # supporting outputs
    prem_unique_community = prem_unique_community,
    stratum_prob_mat       = stratum_prob_mat,
    blended_ratio_comm     = blended_ratio_comm,
    blended_ratio           = blended_ratio,
    stratum_props                   = stratum_props,
    frequency_proportion_diagnostic = frequency_proportion_diagnostic,
    prem_mats = prem_mats
  ),
  "output/MPMmat/DRC_network_input_matrices.rds"
)
cat("\nSaved: output/MPMmat/DRC_network_input_matrices.rds\n")
cat("NOTE: population (WorldPop) is intentionally NOT used anywhere in this\n")
cat("      script. All matrices above are per-capita rates. Population\n")
cat("      weighting belongs to the downstream node/network generation step.\n")
cat("NOTE: close_unique_3wk / physical_unique_3wk (top-level) = 's1' model,\n")
cat("      unchanged for p7 compatibility. All 3 models are also saved\n")
cat("      under *_by_model for comparison / diagnostic use.\n")

# ============================================================
# Step 10: Visualization
# ============================================================

mat_to_long <- function(mat, age_labs) {
  df <- as.data.frame(mat)
  colnames(df) <- age_labs
  df$participant <- age_labs
  pivot_longer(df, cols = -participant,
               names_to = "contact", values_to = "value")
}

make_heatmap_tile <- function(mat,
                              col_title,
                              low_col, high_col,
                              fill_label, scale_limits,
                              show_y = FALSE) {
  df <- mat_to_long(mat, age_labels)
  df$participant <- factor(df$participant, levels = age_labels)
  df$contact     <- factor(df$contact,     levels = age_labels)

  ggplot(df, aes(x = contact, y = participant, fill = value)) +
    geom_tile(color = NA) +
    scale_fill_gradient(low = low_col, high = high_col,
                        limits = scale_limits, na.value = "grey90",
                        name = fill_label) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(
      title = col_title,
      x     = "Contact age",
      y     = if (show_y) "Participant age" else NULL
    ) +
    theme_minimal(base_size = 8) +
    theme(
      plot.title        = element_text(size = 8, face = "bold", hjust = 0.5),
      axis.text         = element_text(size = 7),
      axis.title        = element_text(size = 7),
      legend.position   = "right",
      legend.title      = element_text(size = 7),
      legend.text       = element_text(size = 7),
      legend.key.height = unit(0.45, "cm"),
      panel.grid        = element_blank(),
      plot.margin       = margin(3, 4, 3, 4)
    )
}

# ----------------------------------------------------------
# Plot (a): 1x3 heatmap — Prem close contact  [Panel A]  (unchanged)
# ----------------------------------------------------------
plot_locs_a       <- c("work", "school", "other")
plot_loc_labels_a <- c("Work", "School", "Other")

shared_max_a <- max(sapply(plot_locs_a,
                           function(l) max(prem_mats[[l]], na.rm = TRUE)))

panels_a <- mapply(function(loc, llab, i) {
  make_heatmap_tile(
    prem_mats[[loc]],
    col_title    = llab,
    low_col      = "#FAEEDA", high_col = "#BA7517",
    fill_label   = "contacts/day",
    scale_limits = c(0, shared_max_a),
    show_y       = (i == 1)
  )
}, plot_locs_a, plot_loc_labels_a, seq_along(plot_locs_a), SIMPLIFY = FALSE)

grid_a <- plot_grid(plotlist = panels_a, nrow = 1, ncol = 3)

ggsave("figure/MPMmat/DRC_prem_close_heatmap_1x3.png",
       grid_a, width = 11, height = 3.3, dpi = 300)
message("Saved: figure/MPMmat/DRC_prem_close_heatmap_1x3.png")

# ----------------------------------------------------------
# Plot (b): 1x3 heatmap — blended physical ratio  [Panel B]  (unchanged)
# ----------------------------------------------------------
plot_locs_b       <- c("work", "school", "other")
plot_loc_labels_b <- c("Work", "School", "Other")

panels_b <- mapply(function(loc, llab, i) {
  make_heatmap_tile(
    blended_ratio[[loc]],
    col_title    = llab,
    low_col      = "#E6F1FB", high_col = "#0C447C",
    fill_label   = "proportion",
    scale_limits = c(0, 1),
    show_y       = (i == 1)
  )
}, plot_locs_b, plot_loc_labels_b, seq_along(plot_locs_b), SIMPLIFY = FALSE)

grid_b <- plot_grid(plotlist = panels_b, nrow = 1, ncol = 3)

ggsave("figure/MPMmat/DRC_blended_ratio_heatmap_1x3.png",
       grid_b,
       width = 10, height = 4.5, dpi = 300)
message("Saved: figure/MPMmat/DRC_blended_ratio_heatmap_1x3.png")

# ----------------------------------------------------------
# Plot (c): 1x3 heatmap — household contact structure  [Panel C]  (unchanged)
# ----------------------------------------------------------
p_c1 <- make_heatmap_tile(
  prem_mats$home,
  col_title    = "Close contact (total)",
  low_col      = "#FAEEDA", high_col = "#BA7517",
  fill_label   = "contacts/day",
  scale_limits = c(0, max(prem_mats$home, na.rm = TRUE)),
  show_y       = TRUE
)
p_c2 <- make_heatmap_tile(
  blended_ratio$home,
  col_title    = "Physical contact ratio",
  low_col      = "#E6F1FB", high_col = "#0C447C",
  fill_label   = "proportion",
  scale_limits = c(0, 1),
  show_y       = FALSE
)
p_c3 <- make_heatmap_tile(
  phys_only_home,
  col_title    = "Physical contact",
  low_col      = "#FAEEDA", high_col = "#BA7517",
  fill_label   = "contacts/day",
  scale_limits = c(0, max(prem_mats$home, na.rm = TRUE)),
  show_y       = FALSE
)

grid_c <- plot_grid(p_c1, p_c2, p_c3, nrow = 1)

ggsave("figure/MPMmat/DRC_home_heatmap_1x3.png",
       grid_c, width = 11, height = 3.3, dpi = 300)
message("Saved: figure/MPMmat/DRC_home_heatmap_1x3.png")

# ----------------------------------------------------------
# Plot (d0): 1x3 age-stratified bar chart — 3-week EVENTS
#            (model-independent; unchanged across all 3 models)
# ----------------------------------------------------------
plot_locs_d       <- c("work", "school", "other")
plot_loc_labels_d <- c("Work", "School", "Other")

freq_colors <- c("Daily" = "#993C1D", "Weekly" = "#0F6E56", "Monthly+" = "#185FA5")
stratum_display <- c(daily = "Daily", weekly = "Weekly", monthly = "Monthly+")

build_events_df <- function(loc) {
  rows <- lapply(seq_along(age_labels), function(i) {
    ag   <- age_labels[i]
    vals <- sapply(names(freq_strata), function(st) rowSums(close_events_3wk[[st]][[loc]])[i])
    data.frame(age_group = ag,
               stratum   = unname(stratum_display[names(freq_strata)]),
               contacts  = as.numeric(vals))
  })
  do.call(rbind, rows)
}

y_max_events <- max(unlist(lapply(plot_locs_d, function(loc) {
  df <- build_events_df(loc)
  tapply(df$contacts, df$age_group, sum)
}))) * 1.08

make_bar_generic <- function(df, fill_col, color_map, level_order, llab, y_max,
                             show_col_title = FALSE, show_x = FALSE,
                             show_y_title = FALSE, y_title_text = "") {
  df[[fill_col]] <- factor(df[[fill_col]], levels = level_order)
  df$age_group   <- factor(df$age_group, levels = age_labels)

  ggplot(df, aes(x = age_group, y = contacts, fill = .data[[fill_col]])) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = color_map, name = fill_col) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)),
                       limits = c(0, y_max)) +
    labs(
      title = if (show_col_title) llab else NULL,
      x     = if (show_x) "Age group" else NULL,
      y     = if (show_y_title) y_title_text else NULL
    ) +
    theme_minimal(base_size = 8) +
    theme(
      plot.title         = element_text(size = 9, face = "bold", hjust = 0.5),
      plot.margin        = margin(3, 4, 3, 4),
      axis.text.x        = if (show_x) element_text(size = 5, angle = 45, hjust = 1)
      else element_blank(),
      axis.ticks.x       = element_blank(),
      axis.text.y        = element_text(size = 6),
      axis.title         = element_text(size = 7),
      legend.position    = "none",
      panel.grid.major.x = element_blank()
    )
}

events_levels <- c("Monthly+", "Weekly", "Daily")

panels_events <- mapply(function(loc, llab, i)
  make_bar_generic(build_events_df(loc), "stratum", freq_colors, events_levels,
                   llab, y_max_events,
                   show_col_title = TRUE, show_x = TRUE,
                   show_y_title = (i == 1), y_title_text = "3-week contact events"),
  plot_locs_d, plot_loc_labels_d, seq_along(plot_locs_d), SIMPLIFY = FALSE)

legend_events <- get_legend(
  make_bar_generic(build_events_df("work"), "stratum", freq_colors, events_levels,
                   "Work", y_max_events) +
    theme(legend.position = "right")
)

grid_events <- plot_grid(plotlist = panels_events, nrow = 1, ncol = 3)
grid_events_full <- plot_grid(grid_events, legend_events, nrow = 1, rel_widths = c(1, 0.12))

ggsave("figure/MPMmat/DRC_rarefaction_events_1x3.png",
       grid_events_full, width = 12, height = 4, dpi = 300)
message("Saved: figure/MPMmat/DRC_rarefaction_events_1x3.png")

# ----------------------------------------------------------
# Plot (d1-d3): per-model UNIQUE partner charts (2x3 each):
#   Row 1 = close-only, 3-segment stack (Daily/Weekly/Monthly+)
#   Row 2 = close/physical split, 6-segment stack
# All three models share ONE y-axis scale for direct comparability.
# ----------------------------------------------------------
freq_colors_6 <- c(
  "Monthly+ (non-physical)" = "#B9D4EC",  # light blue
  "Monthly+ (physical)"     = "#185FA5",  # dark blue
  "Weekly (non-physical)"   = "#A9D9CC",  # light green
  "Weekly (physical)"       = "#0F6E56",  # dark green
  "Daily (non-physical)"    = "#E4B6A2",  # light brown
  "Daily (physical)"        = "#993C1D"   # dark brown
)
stack_levels_6 <- names(freq_colors_6)

# Row1 (close, 3-seg) uses levels=c("Monthly+","Weekly","Daily") directly,
# which stacks top-to-bottom as Monthly+(blue) / Weekly(green) / Daily(brown).
# To make Row2 (physical-split, 6-seg) match that SAME top-to-bottom color
# order (blue on top, brown on bottom), the macro group order must also go
# Monthly -> Weekly -> Daily from top to bottom, with physical stacked above
# non-physical within each color pair (unchanged sub-order):
level_order_6 <- c(
  "Monthly+ (physical)", "Monthly+ (non-physical)",
  "Weekly (physical)",   "Weekly (non-physical)",
  "Daily (physical)",    "Daily (non-physical)"
)

build_unique_close_df <- function(loc, model_name) {
  src <- close_unique_3wk_by_model[[model_name]]
  rows <- lapply(seq_along(age_labels), function(i) {
    ag   <- age_labels[i]
    vals <- sapply(names(freq_strata), function(st) rowSums(src[[st]][[loc]])[i])
    data.frame(age_group = ag,
               stratum   = unname(stratum_display[names(freq_strata)]),
               contacts  = as.numeric(vals))
  })
  do.call(rbind, rows)
}

build_unique_split_df <- function(loc, model_name) {
  src_close <- close_unique_3wk_by_model[[model_name]]
  src_phys  <- physical_unique_3wk_by_model[[model_name]]
  rows <- lapply(seq_along(age_labels), function(i) {
    ag <- age_labels[i]
    do.call(rbind, lapply(names(freq_strata), function(st) {
      total_i   <- rowSums(src_close[[st]][[loc]])[i]
      phys_i    <- rowSums(src_phys[[st]][[loc]])[i]
      nonphys_i <- total_i - phys_i
      lab <- stratum_display[[st]]
      data.frame(
        age_group = ag,
        segment   = c(paste0(lab, " (non-physical)"), paste0(lab, " (physical)")),
        contacts  = c(nonphys_i, phys_i)
      )
    }))
  })
  do.call(rbind, rows)
}

# Shared y-axis across ALL models and BOTH chart types
y_max_unique <- max(unlist(lapply(names(model_defs), function(mn) {
  sapply(plot_locs_d, function(loc) {
    df <- build_unique_close_df(loc, mn)
    max(tapply(df$contacts, df$age_group, sum))
  })
}))) * 1.08

message(sprintf("\nShared y_max for unique-partner charts (all models): %.2f", y_max_unique))

render_model_unique_figure <- function(model_name) {
  panels_close <- mapply(function(loc, llab, i)
    make_bar_generic(build_unique_close_df(loc, model_name), "stratum",
                     freq_colors, events_levels, llab, y_max_unique,
                     show_col_title = TRUE, show_x = FALSE,
                     show_y_title = (i == 1), y_title_text = "3-week unique partners"),
    plot_locs_d, plot_loc_labels_d, seq_along(plot_locs_d), SIMPLIFY = FALSE)

  panels_split <- mapply(function(loc, llab, i)
    make_bar_generic(build_unique_split_df(loc, model_name), "segment",
                     freq_colors_6, level_order_6, llab, y_max_unique,
                     show_col_title = FALSE, show_x = TRUE,
                     show_y_title = (i == 1), y_title_text = "3-week unique partners"),
    plot_locs_d, plot_loc_labels_d, seq_along(plot_locs_d), SIMPLIFY = FALSE)

  legend_close <- get_legend(
    make_bar_generic(build_unique_close_df("work", model_name), "stratum",
                     freq_colors, events_levels, "Work", y_max_unique) +
      theme(legend.position = "right")
  )
  legend_split <- get_legend(
    make_bar_generic(build_unique_split_df("work", model_name), "segment",
                     freq_colors_6, level_order_6, "Work", y_max_unique) +
      theme(legend.position = "right")
  )

  grid_body <- plot_grid(
    panels_close[[1]], panels_close[[2]], panels_close[[3]],
    panels_split[[1]], panels_split[[2]], panels_split[[3]],
    nrow = 2, ncol = 3, rel_heights = c(1, 1.15)
  )
  legends_stacked <- plot_grid(legend_close, legend_split, ncol = 1, rel_heights = c(1, 1))
  grid_with_legend <- plot_grid(grid_body, legends_stacked, nrow = 1, rel_widths = c(1, 0.18))

  title_bar <- ggdraw() +
    draw_label(model_labels[model_name], fontface = "bold", size = 12)

  grid_full <- plot_grid(title_bar, grid_with_legend, ncol = 1, rel_heights = c(0.07, 1))

  out_path <- sprintf("figure/MPMmat/DRC_unique_edges_2x3_%s.png", model_name)
  ggsave(out_path, grid_full, width = 12, height = 7.5, dpi = 300)
  message(sprintf("Saved: %s", out_path))
}

for (mn in names(model_defs)) {
  render_model_unique_figure(mn)
}

# ----------------------------------------------------------
# Plot (g1-g3): per-model UNIQUE contact matrix heatmaps (1x3: work/school/
#               other), summed across daily+weekly+monthly. Same orange
#               scale as Plot A, shared across all three models for direct
#               comparison. Cell values shown as integers, no colorbar
#               (legend removed — values are printed directly on tiles).
# ----------------------------------------------------------
plot_locs_g       <- c("work", "school", "other")
plot_loc_labels_g <- c("Work", "School", "Other")

total_unique_mat_by_model <- lapply(names(model_defs), function(mn) {
  src <- close_unique_3wk_by_model[[mn]]
  setNames(
    lapply(plot_locs_g, function(loc)
      src$daily[[loc]] + src$weekly[[loc]] + src$monthly[[loc]]),
    plot_locs_g
  )
})
names(total_unique_mat_by_model) <- names(model_defs)

shared_max_g <- max(sapply(names(model_defs), function(mn)
  max(sapply(plot_locs_g, function(loc)
    max(total_unique_mat_by_model[[mn]][[loc]], na.rm = TRUE)))))

message(sprintf("\nShared color scale max for unique-matrix heatmaps (all models): %.2f",
                shared_max_g))

# Numbered, no-legend heatmap tile (used only for the unique-edge matrices)
make_heatmap_tile_numbered <- function(mat,
                                       col_title,
                                       low_col, high_col,
                                       scale_limits,
                                       show_y = FALSE) {
  df <- mat_to_long(mat, age_labels)
  df$participant <- factor(df$participant, levels = age_labels)
  df$contact     <- factor(df$contact,     levels = age_labels)
  df$label       <- as.character(round(df$value))

  # text color flips to white on dark tiles for readability
  mid_val   <- mean(scale_limits)
  df$txt_col <- ifelse(df$value > mid_val, "white", "black")

  ggplot(df, aes(x = contact, y = participant, fill = value)) +
    geom_tile(color = NA) +
    geom_text(aes(label = label, color = txt_col), size = 1.9) +
    scale_color_identity() +
    scale_fill_gradient(low = low_col, high = high_col,
                        limits = scale_limits, na.value = "grey90",
                        guide = "none") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(
      title = col_title,
      x     = "Contact age",
      y     = if (show_y) "Participant age" else NULL
    ) +
    theme_minimal(base_size = 8) +
    theme(
      plot.title  = element_text(size = 8, face = "bold", hjust = 0.5),
      axis.text   = element_text(size = 7),
      axis.title  = element_text(size = 7),
      panel.grid  = element_blank(),
      plot.margin = margin(3, 4, 3, 4)
    )
}

render_model_matrix_heatmap <- function(model_name) {
  mats_this <- total_unique_mat_by_model[[model_name]]

  panels_g <- mapply(function(loc, llab, i) {
    make_heatmap_tile_numbered(
      mats_this[[loc]],
      col_title    = llab,
      low_col      = "#FAEEDA", high_col = "#BA7517",
      scale_limits = c(0, shared_max_g),
      show_y       = (i == 1)
    )
  }, plot_locs_g, plot_loc_labels_g, seq_along(plot_locs_g), SIMPLIFY = FALSE)

  grid_g <- plot_grid(plotlist = panels_g, nrow = 1, ncol = 3)

  title_bar <- ggdraw() +
    draw_label(model_labels[model_name], fontface = "bold", size = 12)
  grid_full <- plot_grid(title_bar, grid_g, ncol = 1, rel_heights = c(0.1, 1))

  out_path <- sprintf("figure/MPMmat/DRC_unique_matrix_heatmap_1x3_%s.png", model_name)
  ggsave(out_path, grid_full, width = 11, height = 3.7, dpi = 300)
  message(sprintf("Saved: %s", out_path))
}

for (mn in names(model_defs)) {
  render_model_matrix_heatmap(mn)
}

# ----------------------------------------------------------
# Plot (b-2): 1x3 heatmap — POLYMOD original physical ratio
#             + blended ratio + difference (blended - original)
#             (unchanged)
# ----------------------------------------------------------
plot_locs_b2 <- c("work", "school", "other")
plot_loc_labels_b2 <- c("Work", "School", "Other")

raw_ratio <- list(
  work   = get_ratio("work"),
  school = get_ratio("school"),
  other  = get_ratio("other")
)

panels_b2_raw <- mapply(function(loc, llab, i) {
  make_heatmap_tile(
    raw_ratio[[loc]],
    col_title    = llab,
    low_col      = "#E6F1FB", high_col = "#0C447C",
    fill_label   = "proportion",
    scale_limits = c(0, 1),
    show_y       = (i == 1)
  )
}, plot_locs_b2, plot_loc_labels_b2, seq_along(plot_locs_b2), SIMPLIFY = FALSE)

grid_b2_raw <- plot_grid(plotlist = panels_b2_raw, nrow = 1, ncol = 3)

panels_b2_blended <- mapply(function(loc, llab, i) {
  make_heatmap_tile(
    blended_ratio[[loc]],
    col_title    = llab,
    low_col      = "#E6F1FB", high_col = "#0C447C",
    fill_label   = "proportion",
    scale_limits = c(0, 1),
    show_y       = (i == 1)
  )
}, plot_locs_b2, plot_loc_labels_b2, seq_along(plot_locs_b2), SIMPLIFY = FALSE)

grid_b2_blended <- plot_grid(plotlist = panels_b2_blended, nrow = 1, ncol = 3)

make_diff_tile <- function(mat_diff, col_title, show_y = FALSE) {
  df <- mat_to_long(mat_diff, age_labels)
  df$participant <- factor(df$participant, levels = age_labels)
  df$contact     <- factor(df$contact,     levels = age_labels)

  abs_max <- max(abs(df$value), na.rm = TRUE)

  ggplot(df, aes(x = contact, y = participant, fill = value)) +
    geom_tile(color = NA) +
    scale_fill_gradient2(
      low      = "#2166AC",
      mid      = "white",
      high     = "#B2182B",
      midpoint = 0,
      limits   = c(-abs_max, abs_max),
      name     = "difference"
    ) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(
      title = col_title,
      x     = "Contact age",
      y     = if (show_y) "Participant age" else NULL
    ) +
    theme_minimal(base_size = 8) +
    theme(
      plot.title        = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.text         = element_text(size = 5),
      axis.title        = element_text(size = 7),
      legend.position   = "right",
      legend.title      = element_text(size = 6),
      legend.text       = element_text(size = 5),
      legend.key.height = unit(0.45, "cm"),
      panel.grid        = element_blank(),
      plot.margin       = margin(3, 4, 3, 4)
    )
}

panels_b2_diff <- mapply(function(loc, llab, i) {
  diff_mat <- blended_ratio[[loc]] - raw_ratio[[loc]]
  make_diff_tile(diff_mat, col_title = llab, show_y = (i == 1))
}, plot_locs_b2, plot_loc_labels_b2, seq_along(plot_locs_b2), SIMPLIFY = FALSE)

grid_b2_diff <- plot_grid(plotlist = panels_b2_diff, nrow = 1, ncol = 3)

label_raw     <- ggdraw() + draw_label("POLYMOD original", angle = 90, size = 8, fontface = "bold")
label_blended <- ggdraw() + draw_label("Blended (LIC/LMIC)", angle = 90, size = 8, fontface = "bold")
label_diff    <- ggdraw() + draw_label("Difference\n(blended - original)", angle = 90, size = 8, fontface = "bold")

row_width <- c(0.04, 1)

row_raw     <- plot_grid(label_raw,     grid_b2_raw,     nrow = 1, rel_widths = row_width)
row_blended <- plot_grid(label_blended, grid_b2_blended, nrow = 1, rel_widths = row_width)
row_diff    <- plot_grid(label_diff,    grid_b2_diff,    nrow = 1, rel_widths = row_width)

grid_b2_full <- plot_grid(row_raw, row_blended, row_diff,
                          ncol = 1, rel_heights = c(1, 1, 1))

ggsave("figure/MPMmat/DRC_physical_ratio_comparison_3x3.png",
       grid_b2_full,
       width = 11, height = 10, dpi = 300)
message("Saved: figure/MPMmat/DRC_physical_ratio_comparison_3x3.png")

# ----------------------------------------------------------
# Plot (f): POLYMOD raw frequency proportions — stacked bar by age
#           work / school / other (1x3)  (unchanged)
# ----------------------------------------------------------
plot_locs_f       <- c("work", "school", "other")
plot_loc_labels_f <- c("Work", "School", "Other")

freq_colors_3 <- c("Daily" = "#993C1D", "Weekly" = "#0F6E56", "Monthly+" = "#185FA5")

freq_long <- frequency_proportion_diagnostic %>%
  mutate(
    `Monthly+` = monthly + less_often + first_time,
    Daily      = daily,
    Weekly     = weekly
  ) %>%
  select(age_group, location, Daily, Weekly, `Monthly+`) %>%
  pivot_longer(cols = c("Daily", "Weekly", "Monthly+"),
               names_to = "stratum", values_to = "proportion") %>%
  mutate(
    age_group = factor(age_group, levels = age_labels),
    stratum   = factor(stratum, levels = c("Monthly+", "Weekly", "Daily")),
    location  = case_when(
      location == "work"   ~ "Work",
      location == "school" ~ "School",
      location == "other"  ~ "Other",
      location == "home"   ~ "Home"
    )
  ) %>%
  filter(location %in% plot_loc_labels_f)

make_freq_bar <- function(loc_label, show_y = FALSE) {
  df <- freq_long %>% filter(location == loc_label)

  ggplot(df, aes(x = age_group, y = proportion, fill = stratum)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = freq_colors_3, name = "Frequency") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)),
                       limits = c(0, 1),
                       labels = scales::percent) +
    labs(
      title = loc_label,
      x     = "Age group",
      y     = if (show_y) "Proportion" else NULL
    ) +
    theme_minimal(base_size = 8) +
    theme(
      plot.title         = element_text(size = 9, face = "bold", hjust = 0.5),
      plot.margin        = margin(3, 4, 3, 4),
      axis.text.x        = element_text(size = 5, angle = 45, hjust = 1),
      axis.text.y        = element_text(size = 6),
      axis.title         = element_text(size = 7),
      legend.position    = "none",
      panel.grid.major.x = element_blank()
    )
}

panels_f <- mapply(function(llab, i)
  make_freq_bar(llab, show_y = (i == 1)),
  plot_loc_labels_f, seq_along(plot_loc_labels_f), SIMPLIFY = FALSE)

grid_f <- plot_grid(plotlist = panels_f, nrow = 1, ncol = 3)

ggsave("figure/MPMmat/DRC_polymod_freq_proportion_by_age.png",
       grid_f, width = 11, height = 3.5, dpi = 300)
message("Saved: figure/MPMmat/DRC_polymod_freq_proportion_by_age.png")
