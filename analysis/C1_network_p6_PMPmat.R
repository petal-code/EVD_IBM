# ============================================================
# Build frequency-stratified 3-week contact matrices for DRC (Congo)
# and produce visualization figures.
#
# Pipeline:
#   1. Load Prem (Congo) close-contact matrices by location
#      (home / work / school / other), already symmetrized.
#   2. Compute age- and location-specific frequency proportions
#      from POLYMOD frequency_multi (daily / weekly / monthly+).
#   3. Apply frequency-specific multipliers to Prem matrices ->
#      three separate 3-week close contact matrices per location:
#        - prem_3wk_daily    : contacts occurring daily
#        - prem_3wk_weekly   : contacts occurring weekly
#        - prem_3wk_monthly  : contacts occurring monthly or less
#      Each is also re-aggregated into "community".
#   4. Compute location-specific physical/total contact ratio
#      (16x16, from POLYMOD raw data, manually symmetrized).
#   5. Rescale POLYMOD ratio to Mousa et al. LIC/LMIC targets
#      via logit-scale binary search (arithmetic mean matched)
#      -> blended physical ratio.
#   6. Apply blended physical ratio to each frequency-stratified
#      3-week close contact matrix -> physical contact matrices.
#   7. Save all matrices to output/MPMmat/.
#   8. Visualize:
#      (a) 1x3 heatmap: Prem close contact (home / work / other)
#      (b) 1x3 heatmap: blended physical ratio (work / school / other)
#      (c) 1x3 heatmap: household contact structure
#          (Prem close | physical ratio | Prem x ratio)
#      (d) 2x3 bar chart: age-stratified contact-slots and unique
#          contacts by frequency stratum (work / school / other)
#      (e) Population-weighted Panel A/B summary bar chart
#          (work / school / other)
#
# Locations mapping: POLYMOD (cnt_home, cnt_work, cnt_school,
# cnt_transport, cnt_leisure, cnt_otherplace) -> Prem (home, work,
# school, other), where Prem "other" = POLYMOD transport+leisure+otherplace
# ============================================================

# ------------------------------------------------------------
# Libraries
# ------------------------------------------------------------
library(socialmixr)
library(data.table)
library(readxl)
library(ggplot2)
library(tidyr)
library(terra)
library(cowplot)

# ------------------------------------------------------------
# Output directories
# ------------------------------------------------------------
dir.create("figure/MPMmat", recursive = TRUE, showWarnings = FALSE)
dir.create("output/MPMmat", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# Step 0: Config
# ------------------------------------------------------------
prem_dir     <- "data/Prem_contact"
prem_country <- "Congo"
worldpop_dir <- "data/worldpop/DRC_1km"

age_labels <- c("0-4","5-9","10-14","15-19","20-24","25-29",
                "30-34","35-39","40-44","45-49","50-54","55-59",
                "60-64","65-69","70-74","75+")
age_breaks <- c(seq(0, 75, by = 5), Inf)

locations  <- c("home", "work", "school", "other")
loc_labels <- c("Home", "Work", "School", "Other")

# Frequency strata: POLYMOD frequency_multi codes
#   1=daily, 2=weekly, 3=monthly, 4=less_often, 5=first_time
# "monthly" absorbs codes 3/4/5
freq_strata <- list(
  daily   = 1L,
  weekly  = 2L,
  monthly = c(3L, 4L, 5L)
)

freq_category_labels <- c(
  `1` = "daily", `2` = "weekly", `3` = "monthly",
  `4` = "less_often", `5` = "first_time"
)

# Per-stratum 3-week multipliers
# daily  : 21 occurrences in 21 days -> weight = 1 (per-day rate unchanged)
# weekly : 3 occurrences in 21 days  -> weight = 3
# monthly: <1 occurrence             -> weight = 3/4
freq_stratum_weights <- c(daily = 1, weekly = 3, monthly = 3/4)

# Mousa et al. (eLife 2021), Appendix 2-figure 6(A)
# LIC/LMIC physical contact proportions by location
mousa_lic_physical <- c(
  home      = 0.690,
  school    = 0.755,
  community = 0.620,  # maps to Prem "other"
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
# Step 2: Compute frequency proportions from POLYMOD
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

# Mapping: Prem location -> POLYMOD location column(s)
prem_to_polymod_loc <- list(
  home   = "cnt_home",
  work   = "cnt_work",
  school = "cnt_school",
  other  = c("cnt_transport", "cnt_leisure", "cnt_otherplace")
)

# Diagnostic: raw frequency proportions by age group and location
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

message("== Frequency category proportions by age group and location ==")
for (loc in locations) {
  message(sprintf("\n-- %s --", loc))
  sub_diag <- frequency_proportion_diagnostic[
    frequency_proportion_diagnostic$location == loc,
    c("age_group", "daily", "weekly", "monthly", "less_often", "first_time")
  ]
  sub_diag[, -1] <- round(sub_diag[, -1], 3)
  print(sub_diag, row.names = FALSE)
}

# Stratum proportions per age x location
# stratum_props[[location]][[stratum]] = 16-vector of proportions
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

stratum_props <- lapply(prem_to_polymod_loc, compute_stratum_proportions_by_age)

# Impute NA with location-level stratum mean
for (loc in names(stratum_props)) {
  for (st in names(freq_strata)) {
    v <- stratum_props[[loc]][[st]]
    na_idx <- is.na(v)
    if (any(na_idx)) stratum_props[[loc]][[st]][na_idx] <- mean(v, na.rm = TRUE)
  }
}

# ------------------------------------------------------------
# Step 3: Build frequency-stratified 3-week close contact matrices
#
# prem_3wk_<s>[[l]] = Prem[[l]] * proportion_in_stratum(s) * weight(s)
# Symmetrized via (row-perspective + col-perspective) / 2.
# ------------------------------------------------------------
apply_stratum_matrix <- function(mat, prop_vec, stratum_weight) {
  scale_row <- prop_vec * stratum_weight
  scale_col <- prop_vec * stratum_weight
  M_row <- mat * matrix(scale_row, nrow = 16, ncol = 16, byrow = FALSE)
  M_col <- mat * matrix(scale_col, nrow = 16, ncol = 16, byrow = TRUE)
  (M_row + M_col) / 2
}

prem_3wk_daily   <- list()
prem_3wk_weekly  <- list()
prem_3wk_monthly <- list()

for (loc in locations) {
  prem_3wk_daily[[loc]]   <- apply_stratum_matrix(
    prem_mats[[loc]], stratum_props[[loc]]$daily,   freq_stratum_weights["daily"])
  prem_3wk_weekly[[loc]]  <- apply_stratum_matrix(
    prem_mats[[loc]], stratum_props[[loc]]$weekly,  freq_stratum_weights["weekly"])
  prem_3wk_monthly[[loc]] <- apply_stratum_matrix(
    prem_mats[[loc]], stratum_props[[loc]]$monthly, freq_stratum_weights["monthly"])
}

prem_3wk_daily$community   <- prem_3wk_daily$work   + prem_3wk_daily$school   + prem_3wk_daily$other
prem_3wk_weekly$community  <- prem_3wk_weekly$work  + prem_3wk_weekly$school  + prem_3wk_weekly$other
prem_3wk_monthly$community <- prem_3wk_monthly$work + prem_3wk_monthly$school + prem_3wk_monthly$other

message("\n== 3-week close contact rowSums by stratum ==")
for (loc in c("home", "community")) {
  for (st in c("daily", "weekly", "monthly")) {
    mat <- get(paste0("prem_3wk_", st))[[loc]]
    message(sprintf("  %-8s %s: %.2f - %.2f", st, loc,
                    min(rowSums(mat)), max(rowSums(mat))))
  }
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
  # symmetric = FALSE: avoid socialmixr's population-weighted symmetrization,
  # which distorts sparse age cells (e.g. 75+ in work/school).
  # Manual (matrix + transpose) / 2 symmetrization used instead.
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
#
# Logit-scale binary search: finds offset such that the arithmetic
# mean of nonzero cells exactly matches the Mousa LIC/LMIC target.
# Guarantees all values remain in (0, 1).
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

# ------------------------------------------------------------
# Step 6: Apply blended ratio -> frequency-stratified physical matrices
# blended_ratio is frequency-agnostic; same ratio applied per stratum.
# ------------------------------------------------------------
apply_physical_ratio <- function(close_list) {
  phys_list <- list()
  for (loc in locations) {
    raw             <- close_list[[loc]] * blended_ratio[[loc]]
    phys_list[[loc]] <- (raw + t(raw)) / 2
  }
  phys_list$community <- phys_list$work + phys_list$school + phys_list$other
  phys_list
}

prem_3wk_daily_physical   <- apply_physical_ratio(prem_3wk_daily)
prem_3wk_weekly_physical  <- apply_physical_ratio(prem_3wk_weekly)
prem_3wk_monthly_physical <- apply_physical_ratio(prem_3wk_monthly)

# ------------------------------------------------------------
# Step 7: Compute prem_unique and stratum allocation probabilities
#
# prem_unique = sum of three frequency-stratified unique contact matrices.
# Used in p7 network build: sample all contacts once from prem_unique,
# then allocate each edge to daily/weekly/monthly via Multinomial,
# and assign is_physical flag via Bernoulli(blended_ratio).
# This eliminates cross-stratum duplicate edges by construction.
#
# stratum_probs_comm: 16x3 matrix of [p_daily, p_weekly, p_monthly]
# per participant age group, for Multinomial allocation in C++.
# ------------------------------------------------------------

# Community unique contact matrix (sum of three strata)
prem_unique_community <- prem_3wk_daily$community +
  prem_3wk_weekly$community +
  prem_3wk_monthly$community

message("\n== prem_unique_community rowSum range ==")
message(sprintf("  %.2f - %.2f",
                min(rowSums(prem_unique_community)),
                max(rowSums(prem_unique_community))))

# Community stratum allocation probabilities per age group (16x3 matrix)
# Row i: [p_daily_i, p_weekly_i, p_monthly_i] — sums to 1 per row
# Derived from stratum_props (community = mean of work/school/other,
# weighted by contact volume in prem_unique_community)
stratum_props_comm <- list(
  daily   = (stratum_props$work$daily   * rowSums(prem_3wk_daily$work)   +
               stratum_props$school$daily * rowSums(prem_3wk_daily$school) +
               stratum_props$other$daily  * rowSums(prem_3wk_daily$other)) /
    (rowSums(prem_unique_community) + 1e-12),
  weekly  = (stratum_props$work$weekly   * rowSums(prem_3wk_weekly$work)   +
               stratum_props$school$weekly * rowSums(prem_3wk_weekly$school) +
               stratum_props$other$weekly  * rowSums(prem_3wk_weekly$other)) /
    (rowSums(prem_unique_community) + 1e-12),
  monthly = (stratum_props$work$monthly   * rowSums(prem_3wk_monthly$work)   +
               stratum_props$school$monthly * rowSums(prem_3wk_monthly$school) +
               stratum_props$other$monthly  * rowSums(prem_3wk_monthly$other)) /
    (rowSums(prem_unique_community) + 1e-12)
)

# Verify: should sum to ~1 per age group
stratum_check <- stratum_props_comm$daily + stratum_props_comm$weekly +
  stratum_props_comm$monthly
message(sprintf("  Stratum prob sum range: [%.4f, %.4f]",
                min(stratum_check), max(stratum_check)))

# 16x3 matrix for C++ (rows = age groups, cols = daily/weekly/monthly)
stratum_prob_mat <- cbind(
  daily   = stratum_props_comm$daily,
  weekly  = stratum_props_comm$weekly,
  monthly = stratum_props_comm$monthly
)
rownames(stratum_prob_mat) <- age_labels

# Community blended physical ratio (contact-volume weighted mean)
blended_ratio_comm <- (blended_ratio$work   * rowSums(prem_unique_community) +
                         blended_ratio$school * rowSums(prem_unique_community) +
                         blended_ratio$other  * rowSums(prem_unique_community)) /
  (3 * rowSums(prem_unique_community) + 1e-12)
# Simpler: unweighted mean across three locations
blended_ratio_comm <- (blended_ratio$work + blended_ratio$school + blended_ratio$other) / 3

message(sprintf("  blended_ratio_comm mean (nonzero): %.4f",
                mean(blended_ratio_comm[blended_ratio_comm > 0])))

# ------------------------------------------------------------
# Household close-only and physical-only matrices
# Used in simulation for age-group-specific household transmission:
#   p_eff = p_inf_household_close    * close_only_home[age_i, age_j]
#         + p_inf_household_physical * phys_only_home[age_i, age_j]
#
# close_only_home: Prem home daily rate × (1 - blended_ratio$home)
#   → contacts that are close but NOT physical
# phys_only_home:  Prem home daily rate × blended_ratio$home
#   → contacts that are physical (subset of close)
# Both use raw Prem daily rate (no rarefaction) since household
# contacts are effectively permanent relationships.
# ------------------------------------------------------------
close_only_home <- prem_mats$home * (1 - blended_ratio$home)
phys_only_home  <- prem_mats$home * blended_ratio$home

# Symmetrize
close_only_home <- (close_only_home + t(close_only_home)) / 2
phys_only_home  <- (phys_only_home  + t(phys_only_home))  / 2

message("\n== Household contact matrices ==")
message(sprintf("  close_only_home rowSum range: %.3f - %.3f",
                min(rowSums(close_only_home)), max(rowSums(close_only_home))))
message(sprintf("  phys_only_home  rowSum range: %.3f - %.3f",
                min(rowSums(phys_only_home)),  max(rowSums(phys_only_home))))
message(sprintf("  sum check (should equal prem_mats$home rowSums): %.3f - %.3f",
                min(rowSums(close_only_home + phys_only_home)),
                max(rowSums(close_only_home + phys_only_home))))

# ------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------
saveRDS(
  list(
    # Frequency-stratified matrices (retained for diagnostics/visualization)
    close_3wk_daily      = prem_3wk_daily,
    close_3wk_weekly     = prem_3wk_weekly,
    close_3wk_monthly    = prem_3wk_monthly,
    physical_3wk_daily   = prem_3wk_daily_physical,
    physical_3wk_weekly  = prem_3wk_weekly_physical,
    physical_3wk_monthly = prem_3wk_monthly_physical,
    # Network build inputs (used in p7 C++ edge builder)
    prem_unique_community = prem_unique_community,  # 16x16 unique contact matrix
    stratum_prob_mat      = stratum_prob_mat,        # 16x3 stratum allocation probs
    blended_ratio_comm    = blended_ratio_comm,      # 16x16 physical ratio (community)
    blended_ratio         = blended_ratio,           # per-location blended ratios
    # Household simulation inputs
    close_only_home       = close_only_home,  # 16x16: Prem home × (1-phys_ratio)
    phys_only_home        = phys_only_home,   # 16x16: Prem home × phys_ratio
    # Supporting objects
    stratum_props                   = stratum_props,
    frequency_proportion_diagnostic = frequency_proportion_diagnostic
  ),
  "output/MPMmat/DRC_network_input_matrices.rds"
)
cat("\nSaved: output/MPMmat/DRC_network_input_matrices.rds\n")

# ============================================================
# Step 8: Visualization
# ============================================================

# --- Shared helper: 16x16 matrix -> long data.frame ---
mat_to_long <- function(mat, age_labs) {
  df <- as.data.frame(mat)
  colnames(df) <- age_labs
  df$participant <- age_labs
  pivot_longer(df, cols = -participant,
               names_to = "contact", values_to = "value")
}

# --- Shared helper: build one heatmap tile ---
make_heatmap_tile <- function(mat, title, low_col, high_col,
                              fill_label, scale_limits,
                              show_y = FALSE, subtitle = NULL) {
  df <- mat_to_long(mat, age_labels)
  df$participant <- factor(df$participant, levels = age_labels)
  df$contact     <- factor(df$contact,     levels = age_labels)

  ggplot(df, aes(x = contact, y = participant, fill = value)) +
    geom_tile(color = NA) +
    scale_fill_gradient(low = low_col, high = high_col,
                        limits = scale_limits, na.value = "grey90",
                        name = fill_label) +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(title    = title,
         subtitle = subtitle,
         x = "Contact age",
         y = if (show_y) "Participant age" else NULL) +
    theme_minimal(base_size = 8) +
    theme(
      plot.title        = element_text(size = 9, face = "bold", hjust = 0.5),
      plot.subtitle     = element_text(size = 6, hjust = 0.5, color = "grey40"),
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

# ----------------------------------------------------------
# Plot (a): 1x3 heatmap — Prem close contact
#           home / work / other
# ----------------------------------------------------------
plot_locs_a       <- c("home", "work", "other")
plot_loc_labels_a <- c("Home", "Work", "Other")

shared_max_a <- max(sapply(plot_locs_a,
                           function(l) max(prem_mats[[l]], na.rm = TRUE)))

panels_a <- mapply(function(loc, llab, i) {
  make_heatmap_tile(
    prem_mats[[loc]],
    title        = llab,
    low_col      = "#FAEEDA", high_col = "#BA7517",
    fill_label   = "contacts/day",
    scale_limits = c(0, shared_max_a),
    show_y       = (i == 1)
  )
}, plot_locs_a, plot_loc_labels_a, seq_along(plot_locs_a), SIMPLIFY = FALSE)

grid_a <- plot_grid(plotlist = panels_a, nrow = 1, ncol = 3)
title_a <- ggdraw() +
  draw_label("Prem close contact matrices — DRC (home / work / other)",
             fontface = "bold", size = 11, x = 0.5, hjust = 0.5)

ggsave("figure/MPMmat/DRC_prem_close_heatmap_1x3.png",
       plot_grid(title_a, grid_a, ncol = 1, rel_heights = c(0.08, 1)),
       width = 10, height = 4.5, dpi = 300)
message("Saved: figure/MPMmat/DRC_prem_close_heatmap_1x3.png")

# ----------------------------------------------------------
# Plot (b): 1x3 heatmap — blended physical ratio
#           work / school / other
# ----------------------------------------------------------
plot_locs_b       <- c("work", "school", "other")
plot_loc_labels_b <- c("Work", "School", "Other")
mousa_targets     <- c(work = 0.565, school = 0.755, other = 0.620)

panels_b <- mapply(function(loc, llab, i) {
  achieved <- mean(blended_ratio[[loc]][blended_ratio[[loc]] > 0])
  make_heatmap_tile(
    blended_ratio[[loc]],
    title        = llab,
    subtitle     = sprintf("target %.3f  |  achieved %.3f",
                           mousa_targets[loc], achieved),
    low_col      = "#E6F1FB", high_col = "#0C447C",
    fill_label   = "proportion",
    scale_limits = c(0, 1),
    show_y       = (i == 1)
  )
}, plot_locs_b, plot_loc_labels_b, seq_along(plot_locs_b), SIMPLIFY = FALSE)

grid_b <- plot_grid(plotlist = panels_b, nrow = 1, ncol = 3)
title_b <- ggdraw() +
  draw_label("Blended physical contact ratio — DRC (work / school / other)",
             fontface = "bold", size = 11, x = 0.5, hjust = 0.5)
subtitle_b <- ggdraw() +
  draw_label("POLYMOD shape rescaled to Mousa et al. LIC/LMIC targets via logit-scale binary search",
             size = 8, color = "grey40", x = 0.5, hjust = 0.5)

ggsave("figure/MPMmat/DRC_blended_ratio_heatmap_1x3.png",
       plot_grid(title_b, subtitle_b, grid_b,
                 ncol = 1, rel_heights = c(0.07, 0.05, 1)),
       width = 10, height = 4.5, dpi = 300)
message("Saved: figure/MPMmat/DRC_blended_ratio_heatmap_1x3.png")

# ----------------------------------------------------------
# Plot (c): 1x3 heatmap — household contact structure
#   Col 1: Prem close  |  Col 2: physical ratio  |  Col 3: Prem x ratio
# ----------------------------------------------------------
home_physical <- {
  raw <- prem_mats$home * blended_ratio$home
  (raw + t(raw)) / 2
}

p_c1 <- make_heatmap_tile(
  prem_mats$home,
  title        = "Prem close contact (daily rate)",
  low_col      = "#FAEEDA", high_col = "#BA7517",
  fill_label   = "contacts/day",
  scale_limits = c(0, max(prem_mats$home, na.rm = TRUE)),
  show_y       = TRUE
)
p_c2 <- make_heatmap_tile(
  blended_ratio$home,
  title        = "Physical contact ratio",
  low_col      = "#E6F1FB", high_col = "#0C447C",
  fill_label   = "proportion",
  scale_limits = c(0, 1),
  show_y       = FALSE
)
p_c3 <- make_heatmap_tile(
  home_physical,
  title        = "Physical contact (Prem \u00d7 ratio)",
  low_col      = "#FAEEDA", high_col = "#BA7517",
  fill_label   = "contacts/day",
  scale_limits = c(0, max(home_physical, na.rm = TRUE)),
  show_y       = FALSE
)

title_c <- ggdraw() +
  draw_label("Household contact structure — DRC",
             fontface = "bold", size = 11, x = 0.5, hjust = 0.5)
subtitle_c <- ggdraw() +
  draw_label("Left: Prem close contact  |  Centre: blended physical ratio  |  Right: physical contact",
             size = 8, color = "grey40", x = 0.5, hjust = 0.5)

ggsave("figure/MPMmat/DRC_home_heatmap_1x3.png",
       plot_grid(title_c, subtitle_c,
                 plot_grid(p_c1, p_c2, p_c3, nrow = 1),
                 ncol = 1, rel_heights = c(0.07, 0.05, 1)),
       width = 11, height = 4.5, dpi = 300)
message("Saved: figure/MPMmat/DRC_home_heatmap_1x3.png")

# ----------------------------------------------------------
# Plot (d): 2x3 age-stratified bar chart
#   Row 1: contact-slots (x21)  |  Row 2: unique contacts (rarefied)
#   Cols : work / school / other
# ----------------------------------------------------------
plot_locs_d       <- c("work", "school", "other")
plot_loc_labels_d <- c("Work", "School", "Other")

freq_colors <- c("Daily" = "#993C1D", "Weekly" = "#0F6E56", "Monthly+" = "#185FA5")

# Load WorldPop population weights
worldpop_age_map <- list(
  "0-4" = c("00","01"), "5-9" = c("05"), "10-14" = c("10"),
  "15-19" = c("15"), "20-24" = c("20"), "25-29" = c("25"),
  "30-34" = c("30"), "35-39" = c("35"), "40-44" = c("40"),
  "45-49" = c("45"), "50-54" = c("50"), "55-59" = c("55"),
  "60-64" = c("60"), "65-69" = c("65"), "70-74" = c("70"),
  "75+" = c("75","80","85","90")
)

message("Loading WorldPop cod_t_ files...")
pop_by_age <- sapply(age_labels, function(ag) {
  codes   <- worldpop_age_map[[ag]]
  pop_sum <- 0
  for (code in codes) {
    fname   <- file.path(worldpop_dir,
                         sprintf("cod_t_%s_2026_CN_1km_R2025A_UA_v1.tif", code))
    r       <- rast(fname)
    pop_sum <- pop_sum + global(r, "sum", na.rm = TRUE)[[1]]
  }
  pop_sum
})
names(pop_by_age) <- age_labels
pop_weights <- pop_by_age / sum(pop_by_age)
message(sprintf("  Total DRC population: %.0f", sum(pop_by_age)))

get_freq_props_age <- function(loc_name, ag) {
  row <- frequency_proportion_diagnostic[
    frequency_proportion_diagnostic$location  == loc_name &
      frequency_proportion_diagnostic$age_group == ag, ]
  if (nrow(row) == 0) return(list(p_daily = NA, p_weekly = NA, p_monthly_plus = NA))
  list(p_daily        = row$daily,
       p_weekly       = row$weekly,
       p_monthly_plus = row$monthly + row$less_often + row$first_time)
}

build_age_df <- function(loc, panel_label) {
  rows <- lapply(age_labels, function(ag) {
    w           <- pop_weights[ag]
    daily_rate  <- rowSums(prem_mats[[loc]])[which(age_labels == ag)]
    total_slots <- daily_rate * 21
    fp          <- get_freq_props_age(loc, ag)
    if (is.na(fp$p_daily)) return(NULL)

    if (panel_label == "A") {
      vals <- c(w * total_slots * fp$p_daily,
                w * total_slots * fp$p_weekly,
                w * total_slots * fp$p_monthly_plus)
    } else {
      vals <- c(w * total_slots * fp$p_daily        * (1 / 21),
                w * total_slots * fp$p_weekly       * (1 / 7),
                w * total_slots * fp$p_monthly_plus * (3 / 4))
    }
    data.frame(age_group = ag,
               stratum   = c("Daily", "Weekly", "Monthly+"),
               contacts  = vals)
  })
  do.call(rbind, rows)
}

compute_y_max <- function(panel_label) {
  vals <- unlist(lapply(plot_locs_d, function(loc) {
    df <- build_age_df(loc, panel_label)
    tapply(df$contacts, df$age_group, sum)
  }))
  max(vals, na.rm = TRUE) * 1.08
}

y_max_A <- compute_y_max("A")
y_max_B <- compute_y_max("B")

make_bar <- function(loc, llab, panel_label, y_max,
                     show_title = FALSE, show_x = FALSE, show_y_title = FALSE) {
  df <- build_age_df(loc, panel_label)
  df$age_group <- factor(df$age_group, levels = age_labels)
  df$stratum   <- factor(df$stratum, levels = c("Monthly+", "Weekly", "Daily"))

  ggplot(df, aes(x = age_group, y = contacts, fill = stratum)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = freq_colors, name = "Frequency") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)),
                       limits = c(0, y_max)) +
    labs(title = if (show_title) llab else NULL,
         x     = if (show_x) "Age group" else NULL,
         y     = if (show_y_title)
           ifelse(panel_label == "A", "Contact-slots (×21)", "Unique contacts")
         else NULL) +
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

plots_A <- mapply(function(loc, llab, i)
  make_bar(loc, llab, "A", y_max_A,
           show_title = TRUE, show_x = FALSE, show_y_title = (i == 1)),
  plot_locs_d, plot_loc_labels_d, seq_along(plot_locs_d), SIMPLIFY = FALSE)

plots_B <- mapply(function(loc, llab, i)
  make_bar(loc, llab, "B", y_max_B,
           show_title = FALSE, show_x = TRUE, show_y_title = (i == 1)),
  plot_locs_d, plot_loc_labels_d, seq_along(plot_locs_d), SIMPLIFY = FALSE)

shared_legend <- get_legend(
  make_bar("work", "Work", "A", y_max_A) + theme(legend.position = "right")
)

grid_d <- plot_grid(
  plots_A[[1]], plots_A[[2]], plots_A[[3]],
  plots_B[[1]], plots_B[[2]], plots_B[[3]],
  nrow = 2, ncol = 3, rel_heights = c(1, 1.15)
)

title_d <- ggdraw() +
  draw_label("Contact structure by age group — DRC population-weighted (work / school / other)",
             fontface = "bold", size = 11, x = 0.5, hjust = 0.5)
subtitle_d <- ggdraw() +
  draw_label("A: raw contact-slots (daily rate × 21)  |  B: unique contacts (rarefied)",
             size = 8, color = "grey40", x = 0.5, hjust = 0.5)

ggsave("figure/MPMmat/DRC_rarefaction_AB_by_age_2x3.png",
       plot_grid(title_d, subtitle_d,
                 plot_grid(grid_d, shared_legend, nrow = 1, rel_widths = c(1, 0.1)),
                 ncol = 1, rel_heights = c(0.06, 0.04, 1)),
       width = 12, height = 7, dpi = 300)
message("Saved: figure/MPMmat/DRC_rarefaction_AB_by_age_2x3.png")

# ----------------------------------------------------------
# Plot (e): Population-weighted Panel A/B summary bar chart
#           work / school / other
# ----------------------------------------------------------
plot_locs_e       <- c("work", "school", "other")
plot_loc_labels_e <- c("Work", "School", "Other")

panel_df <- do.call(rbind, lapply(seq_along(plot_locs_e), function(i) {
  loc  <- plot_locs_e[i]
  llab <- plot_loc_labels_e[i]

  slot_daily_w <- slot_weekly_w <- slot_monthly_plus_w <- 0

  for (ag in age_labels) {
    w           <- pop_weights[ag]
    daily_rate  <- rowSums(prem_mats[[loc]])[which(age_labels == ag)]
    total_slots <- daily_rate * 21
    fp          <- get_freq_props_age(loc, ag)
    if (is.na(fp$p_daily)) next
    slot_daily_w        <- slot_daily_w        + w * total_slots * fp$p_daily
    slot_weekly_w       <- slot_weekly_w       + w * total_slots * fp$p_weekly
    slot_monthly_plus_w <- slot_monthly_plus_w + w * total_slots * fp$p_monthly_plus
  }

  data.frame(
    location = llab,
    panel    = rep(c("A: Total contact-slots (x21)", "B: Unique contacts (rarefied)"), each = 3),
    category = rep(c("Daily", "Weekly", "Monthly+"), 2),
    contacts = c(slot_daily_w, slot_weekly_w, slot_monthly_plus_w,
                 slot_daily_w * (1/21), slot_weekly_w * (1/7),
                 slot_monthly_plus_w * (3/4))
  )
}))

panel_df$location <- factor(panel_df$location, levels = plot_loc_labels_e)
panel_df$category <- factor(panel_df$category, levels = c("Monthly+", "Weekly", "Daily"))
panel_df$panel    <- factor(panel_df$panel,
                            levels = c("A: Total contact-slots (x21)",
                                       "B: Unique contacts (rarefied)"))

p_e <- ggplot(panel_df, aes(x = location, y = contacts, fill = category)) +
  geom_col(width = 0.65) +
  facet_wrap(~ panel, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = freq_colors, name = "Frequency") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    title    = "Contact structure — DRC population-weighted (work / school / other)",
    subtitle = "A: raw contact-slots (daily rate x 21)  |  B: unique contacts (rarefied)",
    x = "Location", y = "Number of contacts"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title         = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.subtitle      = element_text(size = 8,  hjust = 0.5, color = "grey40"),
    strip.text         = element_text(size = 9,  face = "bold"),
    axis.text.x        = element_text(size = 9),
    axis.text.y        = element_text(size = 8),
    axis.title         = element_text(size = 9),
    legend.position    = "right",
    panel.grid.major.x = element_blank()
  )

ggsave("figure/MPMmat/DRC_rarefaction_panel_AB_no_home.png",
       p_e, width = 9, height = 5, dpi = 300)
message("Saved: figure/MPMmat/DRC_rarefaction_panel_AB_no_home.png")
