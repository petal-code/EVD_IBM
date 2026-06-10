# ==============================================================================
# C1_network_p4_validation.R
# Purpose:
#   Validate synthetic population output from p4:
#     1. Map of household density per cell
#     2. Map of mean household size per cell
#     3. Age distribution per cell (sample)
#     4. Household-level age distribution
#     5. Random household viewer
# ==============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# ==============================================================================
# [Configuration] — Change this to check different territories
# ==============================================================================

synpop_dir <- "output/household"
fig_dir    <- "figure/C1_p4_validation"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# List available files
synpop_files <- list.files(synpop_dir,
                           pattern = "_synthetic_population\\.rds$",
                           full.names = TRUE)
cat("Available territories:\n")
for (i in seq_along(synpop_files))
  cat(sprintf("  [%d] %s\n", i, sub("_synthetic_population\\.rds$", "",
                                    basename(synpop_files[i]))))

# ── Load one territory ─────────────────────────────────────────
# Change index to check different territory
territory_idx <- 1L
pop_file <- synpop_files[territory_idx]
tag      <- sub("_synthetic_population\\.rds$", "", basename(pop_file))

cat(sprintf("\nLoading: %s\n", tag))
res  <- readRDS(pop_file)
pers <- res$individuals   # person_id, hh_id, age
hh   <- res$households    # hh_id, cell_id
cel  <- res$cells         # cell_id, x, y, cell_pop

cat(sprintf("  Individuals : %d\n", nrow(pers)))
cat(sprintf("  Households  : %d\n", nrow(hh)))
cat(sprintf("  Cells       : %d\n", nrow(cel)))

# ==============================================================================
# [1] Join all tables
# ==============================================================================

pop_full <- pers %>%
  left_join(hh,  by = "hh_id") %>%
  left_join(cel, by = "cell_id")

# ==============================================================================
# [2] Cell-level summary
# ==============================================================================

cell_summary <- pop_full %>%
  group_by(cell_id, x, y) %>%
  summarise(
    n_indiv = n(),
    n_hh    = n_distinct(hh_id),
    mean_age = mean(age, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mean_hh_size = n_indiv / n_hh)

# Map 1: Household count per cell
p1 <- ggplot(cell_summary, aes(x = x, y = y, fill = log1p(n_hh))) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma",
                       name   = "log(HH+1)",
                       labels = function(x) round(expm1(x))) +
  coord_fixed() +
  labs(title    = sprintf("%s", tag),
       subtitle = sprintf("Households per cell | total HH: %d", nrow(hh)),
       x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position  = "right",
        plot.title       = element_text(size = 11, face = "bold"),
        plot.subtitle    = element_text(size = 9, color = "grey40"))

# Map 2: Mean household size per cell
p2 <- ggplot(cell_summary, aes(x = x, y = y, fill = mean_hh_size)) +
  geom_raster() +
  scale_fill_viridis_c(option = "plasma", name = "Mean HH size") +
  coord_fixed() +
  labs(title    = "Mean household size per cell",
       x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "right",
        plot.title      = element_text(size = 11, face = "bold"))

# Map 3: Mean age per cell
p3 <- ggplot(cell_summary, aes(x = x, y = y, fill = mean_age)) +
  geom_raster() +
  scale_fill_viridis_c(option = "cividis", name = "Mean age") +
  coord_fixed() +
  labs(title = "Mean age per cell",
       x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "right",
        plot.title      = element_text(size = 11, face = "bold"))

# Save maps
p_maps <- p1 + p2 + p3 + plot_layout(ncol = 3)
ggsave(file.path(fig_dir, sprintf("%s_cell_maps.png", tag)),
       plot = p_maps, width = 18, height = 6, dpi = 150)
cat(sprintf("  Saved: %s_cell_maps.png\n", tag))

# ==============================================================================
# [3] Household size distribution
# ==============================================================================

hh_sizes <- pop_full %>%
  group_by(hh_id) %>%
  summarise(hh_size = n(), .groups = "drop")

p_hh <- ggplot(hh_sizes, aes(x = hh_size)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white", alpha = 0.8) +
  scale_x_continuous(breaks = seq(1, max(hh_sizes$hh_size), by = 2)) +
  labs(title    = sprintf("%s — Household size distribution", tag),
       subtitle = sprintf("Mean=%.1f | Median=%.0f | Max=%d",
                          mean(hh_sizes$hh_size),
                          median(hh_sizes$hh_size),
                          max(hh_sizes$hh_size)),
       x = "Household size", y = "Count") +
  theme_bw()

ggsave(file.path(fig_dir, sprintf("%s_hh_size.png", tag)),
       plot = p_hh, width = 8, height = 5, dpi = 150)
cat(sprintf("  Saved: %s_hh_size.png\n", tag))

# ==============================================================================
# [4] Age distribution overall + by Prem age group
# ==============================================================================

age_grp <- cut(pop_full$age,
               breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40,
                          45, 50, 55, 60, 65, 70, 75, Inf),
               labels = c("0-4","5-9","10-14","15-19","20-24","25-29",
                          "30-34","35-39","40-44","45-49","50-54",
                          "55-59","60-64","65-69","70-74","75+"),
               right  = FALSE)

age_dist <- data.frame(age_grp = age_grp) %>%
  group_by(age_grp) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(pct = 100 * n / sum(n))

p_age <- ggplot(age_dist, aes(x = age_grp, y = pct)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  labs(title    = sprintf("%s — Age distribution (Prem 16 groups)", tag),
       subtitle = sprintf("Mean age=%.1f | Median=%.0f",
                          mean(pop_full$age, na.rm=TRUE),
                          median(pop_full$age, na.rm=TRUE)),
       x = "Age group", y = "% of population") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, sprintf("%s_age_dist.png", tag)),
       plot = p_age, width = 8, height = 5, dpi = 150)
cat(sprintf("  Saved: %s_age_dist.png\n", tag))

# ==============================================================================
# [5] Random household viewer
# ==============================================================================

view_random_households <- function(pop_full, n = 10L, seed_val = 42L) {
  set.seed(seed_val)
  sample_hh_ids <- sample(unique(pop_full$hh_id), min(n, n_distinct(pop_full$hh_id)))

  cat(sprintf("\n=== Random Household Viewer (n=%d) ===\n", length(sample_hh_ids)))
  for (hid in sample_hh_ids) {
    members <- pop_full %>% filter(hh_id == hid)
    cell_xy <- members %>% distinct(x, y)
    cat(sprintf(
      "  HH %-6s | cell (%7.4f, %7.4f) | size=%d | ages: %s\n",
      hid,
      cell_xy$x[1], cell_xy$y[1],
      nrow(members),
      paste(sort(members$age), collapse = ", ")
    ))
  }
}

view_random_households(pop_full, n = 15L, seed_val = 42L)

# ==============================================================================
# [6] Cell-level household age diversity check
# Sample 5 cells and show their household compositions
# ==============================================================================

cat("\n=== Cell-level household composition (5 sample cells) ===\n")
set.seed(42)
sample_cells <- sample(unique(pop_full$cell_id), 5)

for (cid in sample_cells) {
  cell_data <- pop_full %>% filter(cell_id == cid)
  cell_xy   <- cell_data %>% distinct(x, y)
  hh_in_cell <- cell_data %>%
    group_by(hh_id) %>%
    summarise(size = n(),
              ages = paste(sort(age), collapse=","),
              .groups = "drop")

  cat(sprintf("\n  Cell %d (%7.4f, %7.4f) | %d indiv | %d households:\n",
              cid, cell_xy$x[1], cell_xy$y[1],
              nrow(cell_data), nrow(hh_in_cell)))
  for (i in seq_len(min(5, nrow(hh_in_cell)))) {
    cat(sprintf("    HH %s: size=%d ages=[%s]\n",
                hh_in_cell$hh_id[i],
                hh_in_cell$size[i],
                hh_in_cell$ages[i]))
  }
}

cat(sprintf("\nAll figures saved to: %s/\n", fig_dir))
