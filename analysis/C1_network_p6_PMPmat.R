# ============================================================
# Build frequency-stratified 3-week contact matrices for DRC (Congo)
# and produce visualization figures.
# ============================================================

library(socialmixr)
library(data.table)
library(readxl)
library(ggplot2)
library(tidyr)
library(terra)
library(cowplot)

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

freq_strata <- list(
  daily   = 1L,
  weekly  = 2L,
  monthly = c(3L, 4L, 5L)
)

freq_category_labels <- c(
  `1` = "daily", `2` = "weekly", `3` = "monthly",
  `4` = "less_often", `5` = "first_time"
)

freq_stratum_weights <- c(daily = 1, weekly = 3, monthly = 3/4)

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

for (loc in names(stratum_props)) {
  for (st in names(freq_strata)) {
    v <- stratum_props[[loc]][[st]]
    na_idx <- is.na(v)
    if (any(na_idx)) stratum_props[[loc]][[st]][na_idx] <- mean(v, na.rm = TRUE)
  }
}

# ------------------------------------------------------------
# Step 3: Build frequency-stratified 3-week close contact matrices
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

# ------------------------------------------------------------
# Step 6: Apply blended ratio -> frequency-stratified physical matrices
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
# ------------------------------------------------------------
prem_unique_community <- prem_3wk_daily$community +
  prem_3wk_weekly$community +
  prem_3wk_monthly$community

message("\n== prem_unique_community rowSum range ==")
message(sprintf("  %.2f - %.2f",
                min(rowSums(prem_unique_community)),
                max(rowSums(prem_unique_community))))

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

stratum_check <- stratum_props_comm$daily + stratum_props_comm$weekly +
  stratum_props_comm$monthly
message(sprintf("  Stratum prob sum range: [%.4f, %.4f]",
                min(stratum_check), max(stratum_check)))

stratum_prob_mat <- cbind(
  daily   = stratum_props_comm$daily,
  weekly  = stratum_props_comm$weekly,
  monthly = stratum_props_comm$monthly
)
rownames(stratum_prob_mat) <- age_labels

blended_ratio_comm <- (blended_ratio$work + blended_ratio$school + blended_ratio$other) / 3

message(sprintf("  blended_ratio_comm mean (nonzero): %.4f",
                mean(blended_ratio_comm[blended_ratio_comm > 0])))

close_only_home <- prem_mats$home * (1 - blended_ratio$home)
phys_only_home  <- prem_mats$home * blended_ratio$home

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
    close_3wk_daily      = prem_3wk_daily,
    close_3wk_weekly     = prem_3wk_weekly,
    close_3wk_monthly    = prem_3wk_monthly,
    physical_3wk_daily   = prem_3wk_daily_physical,
    physical_3wk_weekly  = prem_3wk_weekly_physical,
    physical_3wk_monthly = prem_3wk_monthly_physical,
    prem_unique_community = prem_unique_community,
    stratum_prob_mat      = stratum_prob_mat,
    blended_ratio_comm    = blended_ratio_comm,
    blended_ratio         = blended_ratio,
    close_only_home       = close_only_home,
    phys_only_home        = phys_only_home,
    stratum_props                   = stratum_props,
    frequency_proportion_diagnostic = frequency_proportion_diagnostic
  ),
  "output/MPMmat/DRC_network_input_matrices.rds"
)
cat("\nSaved: output/MPMmat/DRC_network_input_matrices.rds\n")

# ============================================================
# Step 8: Visualization
# ============================================================

mat_to_long <- function(mat, age_labs) {
  df <- as.data.frame(mat)
  colnames(df) <- age_labs
  df$participant <- age_labs
  pivot_longer(df, cols = -participant,
               names_to = "contact", values_to = "value")
}

# Helper: add panel label outside the plot area using cowplot
# Wraps a ggplot (or grob) in an ggdraw canvas with a label drawn
# at the top-left corner, outside the plot region.
add_panel_label <- function(plot_obj, label,
                            x = 0.01, y = 0.99,
                            size = 12, fontface = "bold") {
  ggdraw(plot_obj) +
    draw_label(label, x = x, y = y,
               hjust = 0, vjust = 1,
               fontface = fontface, size = size)
}

# Shared heatmap builder — no titles, no internal annotations
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
# Plot (a): 1×3 heatmap — Prem close contact  [Panel A]
# ----------------------------------------------------------
plot_locs_a       <- c("home", "work", "other")
plot_loc_labels_a <- c("Home", "Work", "Other")

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
       # add_panel_label(grid_a, "A"),
       width = 11, height = 3.3, dpi = 300)
message("Saved: figure/MPMmat/DRC_prem_close_heatmap_1x3.png")

# ----------------------------------------------------------
# Plot (b): 1×3 heatmap — blended physical ratio  [Panel B]
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
       add_panel_label(grid_b, "B"),
       width = 10, height = 4.5, dpi = 300)
message("Saved: figure/MPMmat/DRC_blended_ratio_heatmap_1x3.png")

# ----------------------------------------------------------
# Plot (c): 1×3 heatmap — household contact structure  [Panel C]
# ----------------------------------------------------------
home_physical <- {
  raw <- prem_mats$home * blended_ratio$home
  (raw + t(raw)) / 2
}

p_c1 <- make_heatmap_tile(
  prem_mats$home,
  col_title    = "Close contact",
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
  home_physical,
  col_title    = "Physical contact",
  low_col      = "#FAEEDA", high_col = "#BA7517",
  fill_label   = "contacts/day",
  scale_limits = c(0, max(prem_mats$home, na.rm = TRUE)),
  show_y       = FALSE
)

grid_c <- plot_grid(p_c1, p_c2, p_c3, nrow = 1)

ggsave("figure/MPMmat/DRC_home_heatmap_1x3.png",
       # add_panel_label(grid_c, "C"),
       width = 11, height = 3.3, dpi = 300)
message("Saved: figure/MPMmat/DRC_home_heatmap_1x3.png")

# ----------------------------------------------------------
# Plot (d): 2×3 age-stratified bar chart  [Panel D]
# ----------------------------------------------------------
plot_locs_d       <- c("work", "school", "other")
plot_loc_labels_d <- c("Work", "School", "Other")

freq_colors <- c("Daily" = "#993C1D", "Weekly" = "#0F6E56", "Monthly+" = "#185FA5")

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

build_age_df <- function(loc, row_label) {
  rows <- lapply(age_labels, function(ag) {
    w           <- pop_weights[ag]
    daily_rate  <- rowSums(prem_mats[[loc]])[which(age_labels == ag)]
    total_slots <- daily_rate * 21
    fp          <- get_freq_props_age(loc, ag)
    if (is.na(fp$p_daily)) return(NULL)

    if (row_label == "A") {
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

compute_y_max <- function(row_label) {
  vals <- unlist(lapply(plot_locs_d, function(loc) {
    df <- build_age_df(loc, row_label)
    tapply(df$contacts, df$age_group, sum)
  }))
  max(vals, na.rm = TRUE) * 1.08
}

y_max_A <- compute_y_max("A")
y_max_B <- compute_y_max("B")

make_bar <- function(loc, llab, row_label,
                     y_max,
                     show_col_title = FALSE,
                     show_x         = FALSE,
                     show_y_title   = FALSE) {
  df <- build_age_df(loc, row_label)
  df$age_group <- factor(df$age_group, levels = age_labels)
  df$stratum   <- factor(df$stratum, levels = c("Monthly+", "Weekly", "Daily"))

  y_title <- if (show_y_title) {
    if (row_label == "A") "Contact-slots (×21)" else "Unique contacts"
  } else NULL

  ggplot(df, aes(x = age_group, y = contacts, fill = stratum)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = freq_colors, name = "Frequency") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)),
                       limits = c(0, y_max)) +
    labs(
      title = if (show_col_title) llab else NULL,
      x     = if (show_x) "Age group" else NULL,
      y     = y_title
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

plots_A <- mapply(function(loc, llab, i)
  make_bar(loc, llab, row_label = "A", y_max = y_max_A,
           show_col_title = TRUE,
           show_x         = FALSE,
           show_y_title   = (i == 1)),
  plot_locs_d, plot_loc_labels_d, seq_along(plot_locs_d), SIMPLIFY = FALSE)

plots_B <- mapply(function(loc, llab, i)
  make_bar(loc, llab, row_label = "B", y_max = y_max_B,
           show_col_title = FALSE,
           show_x         = TRUE,
           show_y_title   = (i == 1)),
  plot_locs_d, plot_loc_labels_d, seq_along(plot_locs_d), SIMPLIFY = FALSE)

shared_legend <- get_legend(
  make_bar("work", "Work", "A", y_max_A) + theme(legend.position = "right")
)

grid_d <- plot_grid(
  plots_A[[1]], plots_A[[2]], plots_A[[3]],
  plots_B[[1]], plots_B[[2]], plots_B[[3]],
  nrow = 2, ncol = 3, rel_heights = c(1, 1.15)
)

# Attach legend, then add panel label outside
grid_d_with_legend <- plot_grid(grid_d, shared_legend,
                                nrow = 1, rel_widths = c(1, 0.1))

ggsave("figure/MPMmat/DRC_rarefaction_AB_by_age_2x3.png",
       add_panel_label(grid_d_with_legend, "D"),
       width = 12, height = 7, dpi = 300)
message("Saved: figure/MPMmat/DRC_rarefaction_AB_by_age_2x3.png")

# ----------------------------------------------------------
# Plot (e): Population-weighted summary bar chart  [Panel E]
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
    panel    = rep(c("A: Contact-slots (×21)", "B: Unique contacts"), each = 3),
    category = rep(c("Daily", "Weekly", "Monthly+"), 2),
    contacts = c(slot_daily_w, slot_weekly_w, slot_monthly_plus_w,
                 slot_daily_w * (1/21), slot_weekly_w * (1/7),
                 slot_monthly_plus_w * (3/4))
  )
}))

panel_df$location <- factor(panel_df$location, levels = plot_loc_labels_e)
panel_df$category <- factor(panel_df$category, levels = c("Monthly+", "Weekly", "Daily"))
panel_df$panel    <- factor(panel_df$panel,
                            levels = c("A: Contact-slots (×21)",
                                       "B: Unique contacts"))

p_e <- ggplot(panel_df, aes(x = location, y = contacts, fill = category)) +
  geom_col(width = 0.65) +
  facet_wrap(~ panel, scales = "fixed", nrow = 1) +
  scale_fill_manual(values = freq_colors, name = "Frequency") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Location", y = "Number of contacts") +
  theme_minimal(base_size = 10) +
  theme(
    strip.text         = element_text(size = 9, face = "bold"),
    axis.text.x        = element_text(size = 9),
    axis.text.y        = element_text(size = 8),
    axis.title         = element_text(size = 9),
    legend.position    = "right",
    panel.grid.major.x = element_blank()
  )

ggsave("figure/MPMmat/DRC_rarefaction_panel_AB_no_home.png",
       add_panel_label(p_e, "E"),
       width = 9, height = 5, dpi = 300)
message("Saved: figure/MPMmat/DRC_rarefaction_panel_AB_no_home.png")


# ----------------------------------------------------------
# Plot (b-2): 1×3 heatmap — POLYMOD 원본 physical ratio
#             + blended ratio + 차이 (blended - original)
#             Cols: work / school / other
# ----------------------------------------------------------
plot_locs_b2 <- c("work", "school", "other")
plot_loc_labels_b2 <- c("Work", "School", "Other")

# 원본 POLYMOD ratio (NA -> 0)
raw_ratio <- list(
  work   = get_ratio("work"),
  school = get_ratio("school"),
  other  = get_ratio("other")
)

# (1) 원본 POLYMOD ratio 패널
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

# (2) blended ratio 패널 (기존 Plot b 재활용)
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

# (3) 차이 패널 (blended - original), 대칭 색상 스케일
make_diff_tile <- function(mat_diff, col_title, show_y = FALSE) {
  df <- mat_to_long(mat_diff, age_labels)
  df$participant <- factor(df$participant, levels = age_labels)
  df$contact     <- factor(df$contact,     levels = age_labels)

  abs_max <- max(abs(df$value), na.rm = TRUE)

  ggplot(df, aes(x = contact, y = participant, fill = value)) +
    geom_tile(color = NA) +
    scale_fill_gradient2(
      low      = "#2166AC",   # 파란색: blended < original
      mid      = "white",
      high     = "#B2182B",   # 빨간색: blended > original
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

# Row 레이블
label_raw     <- ggdraw() + draw_label("POLYMOD original", angle = 90, size = 8, fontface = "bold")
label_blended <- ggdraw() + draw_label("Blended (LIC/LMIC)", angle = 90, size = 8, fontface = "bold")
label_diff    <- ggdraw() + draw_label("Difference\n(blended − original)", angle = 90, size = 8, fontface = "bold")

row_width <- c(0.04, 1)

row_raw     <- plot_grid(label_raw,     grid_b2_raw,     nrow = 1, rel_widths = row_width)
row_blended <- plot_grid(label_blended, grid_b2_blended, nrow = 1, rel_widths = row_width)
row_diff    <- plot_grid(label_diff,    grid_b2_diff,    nrow = 1, rel_widths = row_width)

grid_b2_full <- plot_grid(row_raw, row_blended, row_diff,
                          ncol = 1, rel_heights = c(1, 1, 1))

ggsave("figure/MPMmat/DRC_physical_ratio_comparison_3x3.png",
       add_panel_label(grid_b2_full, "B"),
       width = 11, height = 10, dpi = 300)
message("Saved: figure/MPMmat/DRC_physical_ratio_comparison_3x3.png")


# ----------------------------------------------------------
# Plot (f): POLYMOD 원본 빈도 비율 — 연령별 stacked bar
#           work / school / other  (3×1)
# ----------------------------------------------------------
plot_locs_f       <- c("work", "school", "other")
plot_loc_labels_f <- c("Work", "School", "Other")

freq_colors_3 <- c("Daily" = "#993C1D", "Weekly" = "#0F6E56", "Monthly+" = "#185FA5")

# frequency_proportion_diagnostic에서 long format으로 변환
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

shared_legend_f <- get_legend(
  make_freq_bar("Work") + theme(legend.position = "right")
)

grid_f <- plot_grid(plotlist = panels_f, nrow = 1, ncol = 3)

ggsave("figure/MPMmat/DRC_polymod_freq_proportion_by_age.png",
       grid_f, width = 11, height = 3.5, dpi = 300)
message("Saved: figure/MPMmat/DRC_polymod_freq_proportion_by_age.png")
