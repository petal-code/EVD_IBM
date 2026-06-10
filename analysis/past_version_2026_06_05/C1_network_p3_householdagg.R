library(rdhs); library(haven); library(ggplot2); library(dplyr)

# ── Load DHS data ──────────────────────────────────────────────
set_rdhs_config(email = "charles.whittaker16@imperial.ac.uk",
                project = "mpox Outbreak Response - Understanding Household Size Distributions in DRC")
survs    <- dhs_surveys(countryIds = c("CD"), surveyYear = 2013)
datasets <- dhs_datasets(surveyIds = survs$SurveyId, fileType = "HR", fileFormat = "flat")
datasets$path <- unlist(get_datasets(datasets$FileName))

x <- readRDS(datasets$path)
y <- x[, c("hv001", "hv004", "hv024", "hv002", "hv005", "hv009")]
y$province <- haven::as_factor(y$hv024)

# ── Extract household age compositions by province ────────────
hh_age_labels <- c(paste0("hv105_0", 1:9),
                   paste0("hv105_1", 0:9),
                   paste0("hv105_2", 0:4))

cat("Extracting household age data...\n")
for (i in 1:nrow(y)) {
  temp_ages <- x[i, hh_age_labels]
  temp_ages <- temp_ages[which(!is.na(temp_ages))]
  y$hh_size_from_ages[i] <- length(temp_ages)
  y$hh_ages[i]           <- list(unname(temp_ages))
  if (i %% 1000 == 0) cat(sprintf("  %d / %d\n", i, nrow(y)))
}

# ── Summarise household structure per province ─────────────────
province_summary <- y %>%
  group_by(province) %>%
  summarise(
    n_hh          = n(),
    mean_hh_size  = mean(hh_size_from_ages),
    median_hh_size = median(hh_size_from_ages),
    sd_hh_size    = sd(hh_size_from_ages),
    .groups = "drop"
  )

print(province_summary)

# ── Plot 1: Household size distribution by province ───────────
ggplot(y, aes(x = hh_size_from_ages, fill = province)) +
  geom_histogram(bins = 20, col = "black", alpha = 0.8) +
  facet_wrap(. ~ province, scales = "free_y") +
  labs(title = "Household size distribution by province (DRC DHS 2013)",
       x = "Household size", y = "Count") +
  theme_bw() +
  theme(legend.position = "none")

# ── Plot 2: Age distribution by province ──────────────────────
age_long <- y %>%
  select(province, hh_ages) %>%
  rowwise() %>%
  mutate(ages = list(unlist(hh_ages))) %>%
  tidyr::unnest(ages)

ggplot(age_long, aes(x = ages, fill = province)) +
  geom_histogram(bins = 30, col = "black", alpha = 0.8) +
  facet_wrap(. ~ province, scales = "free_y") +
  labs(title = "Age distribution by province (DRC DHS 2013)",
       x = "Age", y = "Count") +
  theme_bw() +
  theme(legend.position = "none")


# ── Save household structure data by province ──────────────────
output_dir <- "output/household"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save 1: summary table as CSV
write.csv(province_summary,
          file.path(output_dir, "province_household_summary.csv"),
          row.names = FALSE)

# Save 2: full household-level data (with ages) as RDS
# Keep only essential columns
y_save <- y %>%
  select(province, hv001, hv002, hv005, hv009, hh_size_from_ages, hh_ages)

saveRDS(y_save,
        file.path(output_dir, "province_household_data.rds"))

# Save 3: per-province household list (for easy lookup later)
# Format: named list where each element = one province's households
province_hh_list <- y_save %>%
  group_by(province) %>%
  group_split() %>%
  setNames(levels(y_save$province))

saveRDS(province_hh_list,
        file.path(output_dir, "province_household_list.rds"))

# Save 4: plots
ggsave(file.path(output_dir, "plot_hh_size_by_province.png"),
       plot = last_plot(),   # age distribution plot
       width = 14, height = 10, dpi = 150)

cat("Saved files:\n")
list.files(output_dir)
