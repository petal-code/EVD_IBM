# ==============================================================================
# DRC_scale_comparison_combined.R
#
# Merges DRC_total_scale_comparison_1x3.R and DRC_weekly_scale_comparison_1x3.R
# into one script (shared data loading + model functions), and fixes the
# population-weighting scope:
#
#   - "work"   representative individual: population-weighted using ONLY
#              ages 20+ (adults — no one under 20 is meaningfully "at work")
#   - "school" representative individual: population-weighted using ONLY
#              ages 5-19 (school-age population)
#   - "other"  representative individual: population-weighted using the
#              FULL population (everyone can plausibly have "other" contacts)
#
# Using the full population for work/school would dilute the representative
# rate with age groups that essentially never participate in that setting,
# understating E_total for a person who is actually a worker/student.
#
# Produces two figures:
#   (A) DRC_total_scale_comparison_1x3.png  — daily+weekly+monthly COMBINED,
#       3 models (Linear / Bernoulli s=1 / Bernoulli s=0.5), 28 days,
#       vertical line at day 21. This is the corrected/primary comparison.
#   (B) DRC_weekly_scale_comparison_1x3.png — same 3 models, WEEKLY stratum
#       only, kept as a secondary breakdown view.
#
# Model definitions (unchanged from before, with the T_eff >= 1 day clamp
# that fixes the s=0.5 daily-stratum NaN bug):
#   Linear:            cumulative(t) = E_relevant * t / L
#   Bernoulli(s):       T_eff = max(T_avg * s, 1); p_eff = 1/T_eff
#                       N_eff = E * T_eff / L
#                       unique(t) = N_eff * [1 - (1-p_eff)^t]
#   Combined Bernoulli: sum of the per-stratum Bernoulli(s) curves
# ==============================================================================

library(ggplot2)
library(cowplot)
library(terra)

matrices_path <- "output/MPMmat/DRC_network_input_matrices.rds"
worldpop_dir  <- "data/worldpop/DRC_1km"
dir.create("figure/MPMmat", recursive = TRUE, showWarnings = FALSE)

mats <- readRDS(matrices_path)

split_locations <- c("work", "school", "other")
loc_labels       <- c("Work", "School", "Other")
strata           <- c("daily", "weekly", "monthly")

T_avg    <- c(daily = 1, weekly = 7, monthly = 28)   # nominal revisit period (days)
window_L <- 21                                        # p6's 3-week observation window
t_max    <- 28                                        # show 4 weeks

age_labels <- c("0-4","5-9","10-14","15-19","20-24","25-29",
                "30-34","35-39","40-44","45-49","50-54","55-59",
                "60-64","65-69","70-74","75+")

# ------------------------------------------------------------
# Step 0: WorldPop population weights, then per-location relevant subsets
# ------------------------------------------------------------
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
pop_weights_full <- pop_by_age / sum(pop_by_age)

# Location -> relevant age subset (indices into age_labels)
loc_age_subset <- list(
  work   = which(age_labels %in% c("20-24","25-29","30-34","35-39","40-44",
                                   "45-49","50-54","55-59","60-64","65-69",
                                   "70-74","75+")),                       # 20+
  school = which(age_labels %in% c("5-9","10-14","15-19")),              # 5-19
  other  = seq_along(age_labels)                                          # full population
)

message("\n== Population-relevant age subset per location ==")
for (loc in split_locations) {
  message(sprintf("  %-6s: %s", loc,
                  paste(age_labels[loc_age_subset[[loc]]], collapse = ", ")))
}

# Population-weighted mean restricted to a location's relevant age subset
# (weights renormalized to sum to 1 within the subset)
pop_weighted_mean_loc <- function(x, loc) {
  idx <- loc_age_subset[[loc]]
  w   <- pop_weights_full[idx]
  w   <- w / sum(w)
  sum(x[idx] * w)
}

# ------------------------------------------------------------
# Step 1: representative E per stratum x location, using location-specific
#         population-weighted mean
# ------------------------------------------------------------
E_rep <- sapply(split_locations, function(loc) {
  sapply(strata, function(st) {
    pop_weighted_mean_loc(rowSums(mats$close_events_3wk[[st]][[loc]]), loc)
  })
})
rownames(E_rep) <- strata
colnames(E_rep) <- split_locations

E_total <- colSums(E_rep)   # daily+weekly+monthly combined, per location
E_weekly <- E_rep["weekly", ]

message("\n== E_rep (location-relevant pop-weighted 21-day events) ==")
print(round(E_rep, 3))
message("\n== E_total (daily+weekly+monthly combined) ==")
print(round(E_total, 3))

# ------------------------------------------------------------
# Step 2: model functions
# ------------------------------------------------------------
linear_curve <- function(E_relevant, L, t_seq) {
  (E_relevant / L) * t_seq
}

bernoulli_scaled_curve <- function(E, T_avg_days, s, L, t_seq) {
  # T_eff clamped to >= 1 day: a discrete daily-Bernoulli model can't
  # represent "more than once per day." Without this, daily (T_avg=1)
  # with s<1 pushes p_eff above 1, producing NaN for non-integer t.
  T_eff <- pmax(T_avg_days * s, 1)
  p_eff <- 1 / T_eff
  N_eff <- E * T_eff / L
  N_eff * (1 - (1 - p_eff)^t_seq)
}

bernoulli_combined <- function(loc, s, t_seq) {
  Reduce(`+`, lapply(strata, function(st)
    bernoulli_scaled_curve(E_rep[st, loc], T_avg[st], s, window_L, t_seq)))
}

t_seq <- seq(0, t_max, by = 0.25)

model_colors <- c(
  "Linear (fully independent)" = "#B5651D",
  "Bernoulli, s = 1" = "#0F6E56",
  "Bernoulli, s = 0.5"         = "#185FA5"
)

# ------------------------------------------------------------
# Step 3: sanity checks at t = 21
# ------------------------------------------------------------
message("\n== Sanity check @t=21: combined s=1 vs p6 saved close_unique_3wk (summed) ==")
for (loc in split_locations) {
  at21  <- bernoulli_combined(loc, 1.0, window_L)
  saved <- sum(sapply(strata, function(st)
    pop_weighted_mean_loc(rowSums(mats$close_unique_3wk[[st]][[loc]]), loc)))
  message(sprintf("  %-6s: curve@21=%.4f | saved=%.4f | diff=%.2e",
                  loc, at21, saved, abs(at21 - saved)))
}

message("\n== Sanity check @t=21: weekly-only s=1 vs p6 saved close_unique_3wk$weekly ==")
for (loc in split_locations) {
  at21  <- bernoulli_scaled_curve(E_weekly[loc], T_avg["weekly"], 1.0, window_L, window_L)
  saved <- pop_weighted_mean_loc(rowSums(mats$close_unique_3wk$weekly[[loc]]), loc)
  message(sprintf("  %-6s: curve@21=%.4f | saved=%.4f | diff=%.2e",
                  loc, at21, saved, abs(at21 - saved)))
}

# ------------------------------------------------------------
# Step 4: shared panel builder
# ------------------------------------------------------------
make_panel <- function(loc, llab, df_builder, y_label, show_y = FALSE, show_legend = FALSE) {
  df <- df_builder(loc)
  df$model <- factor(df$model, levels = names(model_colors))

  ggplot(df, aes(x = day, y = unique, color = model)) +
    geom_line(linewidth = 0.9) +
    geom_vline(xintercept = window_L, linetype = "dashed",
               color = "grey40", linewidth = 0.5) +
    scale_color_manual(values = model_colors, name = "Model") +
    scale_x_continuous(breaks = seq(0, t_max, by = 7)) +
    labs(
      title = llab,
      x     = "Days",
      y     = if (show_y) y_label else NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title       = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title       = element_text(size = 8),
      axis.text        = element_text(size = 7),
      legend.position  = if (show_legend) "right" else "none",
      panel.grid.minor = element_blank()
    )
}

render_grid <- function(df_builder, y_label, out_path, width = 12, height = 3) {
  panels <- lapply(seq_along(split_locations), function(i)
    make_panel(split_locations[i], loc_labels[i], df_builder, y_label, show_y = (i == 1)))
  legend_shared <- get_legend(
    make_panel(split_locations[1], loc_labels[1], df_builder, y_label, show_legend = TRUE)
  )
  grid_main <- plot_grid(plotlist = panels, nrow = 1, ncol = 3)
  grid_full <- plot_grid(grid_main, legend_shared, nrow = 1, rel_widths = c(1, 0.28))
  ggsave(out_path, grid_full, width = width, height = height, dpi = 300)
  message(sprintf("Saved: %s", out_path))
}

# ------------------------------------------------------------
# Figure A: daily+weekly+monthly COMBINED (primary/corrected figure)
# ------------------------------------------------------------
build_df_total <- function(loc) {
  rbind(
    data.frame(day = t_seq, unique = linear_curve(E_total[loc], window_L, t_seq),
               model = "Linear (fully independent)"),
    data.frame(day = t_seq, unique = bernoulli_combined(loc, 1.0, t_seq),
               model = "Bernoulli, s = 1"),
    data.frame(day = t_seq, unique = bernoulli_combined(loc, 0.5, t_seq),
               model = "Bernoulli, s = 0.5")
  )
}

render_grid(
  build_df_total,
  y_label  = "Cumulative unique edges (location-relevant pop-weighted avg individual)",
  out_path = "figure/MPMmat/DRC_total_scale_comparison_1x3.png"
)

# ------------------------------------------------------------
# Figure B: WEEKLY-only breakdown (secondary figure)
# ------------------------------------------------------------
build_df_weekly <- function(loc) {
  rbind(
    data.frame(day = t_seq, unique = linear_curve(E_weekly[loc], window_L, t_seq),
               model = "Linear (fully independent)"),
    data.frame(day = t_seq,
               unique = bernoulli_scaled_curve(E_weekly[loc], T_avg["weekly"], 1.0, window_L, t_seq),
               model = "Bernoulli, s = 1"),
    data.frame(day = t_seq,
               unique = bernoulli_scaled_curve(E_weekly[loc], T_avg["weekly"], 0.5, window_L, t_seq),
               model = "Bernoulli, s = 0.5")
  )
}

render_grid(
  build_df_weekly,
  y_label  = "Cumulative unique weekly edges\n(location-relevant pop-weighted avg individual)",
  out_path = "figure/MPMmat/DRC_weekly_scale_comparison_1x3.png"
)
