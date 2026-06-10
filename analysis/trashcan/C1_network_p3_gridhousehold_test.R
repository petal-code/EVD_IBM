library(dplyr)

# ── Load pre-built household data (from DHS) ──────────────────
province_hh_list <- readRDS("output/household/province_household_list.rds")

# ── Province-to-DHS mapping (old 11 → new 26 provinces) ───────
# Kinshasa maps directly
province_to_dhs <- list(
  "Kinshasa" = "kinshasa"
)

# ── Step 1: Convert fractional population to integer counts ───
# Using Poisson rounding to preserve total population distribution
set.seed(42)
pop_int <- pop_df
age_sex_cols <- names(pop_df)[-(1:2)]

for (col in age_sex_cols) {
  # rpois handles fractional expected values probabilistically
  pop_int[[col]] <- rpois(nrow(pop_df), lambda = pmax(pop_df[[col]], 0))
}

cat(sprintf("Total population after rounding: %.0f\n",
            sum(pop_int[, age_sex_cols])))

# ── Step 2: Build age-sex lookup table per grid cell ──────────
# Reshape to long format: each row = one person's (cell, sex, age_group)
pop_long <- pop_int %>%
  tidyr::pivot_longer(cols = all_of(age_sex_cols),
                      names_to  = "agesex",
                      values_to = "count") %>%
  filter(count > 0) %>%
  mutate(
    sex       = substr(agesex, 1, 1),                        # "f" or "m"
    age_lower = as.integer(gsub("[fm]", "", agesex))         # numeric age lower bound
  )

# Expand: repeat each row by count (one row per individual)
pop_individuals <- pop_long %>%
  filter(count > 0) %>%
  tidyr::uncount(count) %>%                                  # expand by count
  mutate(
    # Assign random age within 5-year bin
    age = age_lower + sample(0:4, n(), replace = TRUE)
  ) %>%
  select(x, y, sex, age)

cat(sprintf("Total individuals generated: %d\n", nrow(pop_individuals)))
cat("Age-sex distribution sample:\n")
print(head(pop_individuals, 10))

# ── Step 3: Assign individuals to households per grid cell ─────
# Strategy:
# 1. For each cell, get total population count
# 2. Sample households from DHS until population is filled
# 3. Assign individuals to households, matching age-sex distribution

dhs_province <- province_to_dhs[[province_name]]
dhs_hh       <- province_hh_list[[dhs_province]]

# Helper: sample households from DHS until target population is reached
sample_households <- function(target_n, dhs_hh) {
  households <- list()
  total      <- 0
  hh_id      <- 1

  while (total < target_n) {
    # Sample one household from DHS (weighted by survey weight)
    idx    <- sample(nrow(dhs_hh), 1, prob = dhs_hh$hv005 / sum(dhs_hh$hv005))
    hh     <- dhs_hh[idx, ]
    ages   <- unlist(hh$hh_ages)

    if (length(ages) == 0) next

    households[[hh_id]] <- data.frame(
      hh_id = hh_id,
      age   = ages
    )
    total  <- total + length(ages)
    hh_id  <- hh_id + 1
  }

  # Trim excess individuals to match target_n exactly
  result <- bind_rows(households)
  result[1:min(nrow(result), target_n), ]
}

# ── Step 4: Process each grid cell ────────────────────────────
cat("Assigning households to grid cells...\n")

cells <- pop_int %>%
  mutate(total_pop = rowSums(across(all_of(age_sex_cols)))) %>%
  filter(total_pop > 0) %>%
  select(x, y, total_pop)

cat(sprintf("Total cells to process: %d\n", nrow(cells)))
cat(sprintf("Total population: %.0f\n", sum(cells$total_pop)))

# Process (test with first 100 cells first)
set.seed(42)
test_cells <- head(cells, 100)

household_list <- lapply(1:nrow(test_cells), function(i) {
  cell     <- test_cells[i, ]
  hh_data  <- sample_households(cell$total_pop, dhs_hh)

  # Attach cell coordinates
  hh_data$cell_x   <- cell$x
  hh_data$cell_y   <- cell$y
  # Random location within 100m cell
  hh_data$indiv_x  <- cell$x + runif(nrow(hh_data), -0.0004, 0.0004)
  hh_data$indiv_y  <- cell$y + runif(nrow(hh_data), -0.0004, 0.0004)

  # Make hh_id unique across cells
  hh_data$hh_id    <- paste0("cell", i, "_hh", hh_data$hh_id)
  hh_data

  if (i %% 20 == 0) cat(sprintf("  Processed %d / %d cells\n", i, nrow(test_cells)))

  hh_data
})

result_test <- bind_rows(household_list)
cat(sprintf("\nTest result: %d individuals in %d unique households\n",
            nrow(result_test),
            length(unique(result_test$hh_id))))

head(result_test, 15)


library(dplyr)
library(tidyr)

# ── Full run with progress tracking ───────────────────────────
set.seed(42)

cat(sprintf("Processing %d cells, total population: %.0f\n",
            nrow(cells), sum(cells$total_pop)))

start_time <- Sys.time()

household_list_full <- lapply(1:nrow(cells), function(i) {
  cell    <- cells[i, ]
  hh_data <- sample_households(cell$total_pop, dhs_hh)

  # Attach coordinates
  hh_data$cell_x  <- cell$x
  hh_data$cell_y  <- cell$y

  # Random scatter within ~100m cell boundary
  hh_data$indiv_x <- cell$x + runif(nrow(hh_data), -0.0004, 0.0004)
  hh_data$indiv_y <- cell$y + runif(nrow(hh_data), -0.0004, 0.0004)

  # Unique household ID across all cells
  hh_data$hh_id   <- paste0("cell", i, "_hh", hh_data$hh_id)

  # Progress log every 5000 cells
  if (i %% 1000 == 0) {
    elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
    cat(sprintf("  [%s min] %d / %d cells (%.1f%%)\n",
                elapsed, i, nrow(cells), 100 * i / nrow(cells)))
  }

  hh_data
})

result_full <- bind_rows(household_list_full)

elapsed_total <- round(difftime(Sys.time(), start_time, units = "mins"), 1)
cat(sprintf("\nDone in %s min\n", elapsed_total))
cat(sprintf("Total individuals : %d\n", nrow(result_full)))
cat(sprintf("Total households  : %d\n", length(unique(result_full$hh_id))))

# ── Save result ────────────────────────────────────────────────
saveRDS(result_full,
        file.path(output_dir, sprintf("%s_synthetic_population.rds", province_name)))

cat(sprintf("Saved to output/household/%s_synthetic_population.rds\n", province_name))



library(ggplot2)
library(dplyr)
library(tidyr)

# ── Load data ──────────────────────────────────────────────────
pop      <- readRDS("output/household/Kinshasa_synthetic_population.rds")
dhs_kin  <- province_hh_list[["kinshasa"]]

# ── Plot 1: Household size distribution (synthetic vs DHS) ────
hh_sizes_synthetic <- pop %>%
  group_by(hh_id) %>%
  summarise(hh_size = n(), .groups = "drop") %>%
  mutate(source = "Synthetic")

hh_sizes_dhs <- dhs_kin %>%
  mutate(hh_size = hh_size_from_ages,
         source  = "DHS 2013") %>%
  select(hh_size, source)

hh_sizes_combined <- bind_rows(hh_sizes_synthetic, hh_sizes_dhs)

p1 <- ggplot(hh_sizes_combined, aes(x = hh_size, fill = source)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 30, color = "black", alpha = 0.7,
                 position = "identity") +
  scale_fill_manual(values = c("Synthetic" = "steelblue", "DHS 2013" = "coral")) +
  coord_cartesian(xlim = c(0, 25)) +
  labs(title = "Household size: Synthetic vs DHS",
       x = "Household size", y = "Density", fill = "") +
  theme_bw()

# ── Plot 2: Age distribution (synthetic vs DHS) ───────────────
ages_synthetic <- pop %>%
  mutate(source = "Synthetic") %>%
  select(age, source)

ages_dhs <- dhs_kin %>%
  rowwise() %>%
  mutate(ages = list(unlist(hh_ages))) %>%
  unnest(ages) %>%
  mutate(source = "DHS 2013",
         age    = as.integer(ages)) %>%
  select(age, source)

ages_combined <- bind_rows(ages_synthetic, ages_dhs)

p2 <- ggplot(ages_combined, aes(x = age, fill = source)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40, color = "black", alpha = 0.7,
                 position = "identity") +
  scale_fill_manual(values = c("Synthetic" = "steelblue", "DHS 2013" = "coral")) +
  labs(title = "Age distribution: Synthetic vs DHS",
       x = "Age", y = "Density", fill = "") +
  theme_bw()

# ── Plot 3: Spatial density heatmap ───────────────────────────
# Sample to avoid overplotting
pop_sample <- pop %>% slice_sample(n = 100000)

p3 <- ggplot(pop_sample, aes(x = indiv_x, y = indiv_y)) +
  stat_density_2d(aes(fill = after_stat(density)),
                  geom = "raster", contour = FALSE) +
  scale_fill_viridis_c(option = "magma", name = "Density") +
  coord_fixed() +
  labs(title = "Population density - Kinshasa",
       x = "Longitude", y = "Latitude") +
  theme_bw()

# ── Plot 4: Household size spatial variation ──────────────────
# Average household size per grid cell
hh_cell_size <- pop %>%
  group_by(cell_x, cell_y, hh_id) %>%
  summarise(hh_size = n(), .groups = "drop") %>%
  group_by(cell_x, cell_y) %>%
  summarise(mean_hh_size = mean(hh_size), .groups = "drop")

p4 <- ggplot(hh_cell_size, aes(x = cell_x, y = cell_y, fill = mean_hh_size)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", name = "Mean HH size") +
  coord_fixed() +
  labs(title = "Mean household size per grid cell - Kinshasa",
       x = "Longitude", y = "Latitude") +
  theme_bw()

# ── Combine and save all plots ─────────────────────────────────
library(patchwork)

p_combined <- (p1 + p2) / (p3 + p4) +
  plot_annotation(title    = "Kinshasa Synthetic Population Summary",
                  subtitle = sprintf("%.0f individuals | %.0f households",
                                     nrow(pop),
                                     length(unique(pop$hh_id))))

ggsave("output/household/Kinshasa_population_summary.png",
       plot   = p_combined,
       width  = 16, height = 12, dpi = 150)

print(p_combined)



