# ==============================================================================
# C1_network_p2_lv2dist_map.R
# Purpose:
#   1. Load WorldPop R2025A 2025 population rasters for DRC
#   2. For each Level 2 territory:
#      - Crop and mask raster to territory boundary
#      - Save population dataframe as RDS (used by p4)
#      - Save population density map as PNG
#
# Inputs:
#   data/worldpop/DRC_1km   : WorldPop R2025A GeoTIFF files (_t_ total files only)
#   data/shpmap/     : GADM Level 1 and Level 2 shapefiles
#
# Outputs:
#   output/popdata/{NAME_1}_{NAME_2}_pop.rds  — population dataframe per territory
#   figure/C1_p2/{NAME_1}_{NAME_2}.png        — population density map per territory
# ==============================================================================

library(terra)
library(sf)
library(dplyr)
library(ggplot2)

# ==============================================================================
# [Section 1] Load shapefiles
# ==============================================================================

message("=== Section 1: Loading shapefiles ===")

cod2_sf <- st_as_sf(vect("data/shpmap/gadm41_COD_2.shp"))
message(sprintf("  Territories (Level 2) : %d", nrow(cod2_sf)))

# ==============================================================================
# [Section 2] Load and sum WorldPop R2025A population rasters
# Use _t_ files (age-stratified totals = female + male already combined)
# Avoids loading 40 f/m files — only 20 _t_ files needed
# ==============================================================================

message("\n=== Section 2: Loading population rasters ===")

files_all <- list.files("data/worldpop/DRC_1km", pattern = "\\.tif$", full.names = TRUE)
files_t   <- files_all[grepl("/cod_t_", files_all)]  # lowercase _t_ only (excludes _T_F_, _T_M_)
message(sprintf("  Using %d files (age-stratified totals)", length(files_t)))

# Sum all age layers to get total population raster
pop_rast <- NULL
for (f in files_t) {
  r <- rast(f)
  if (is.null(pop_rast)) pop_rast <- r else pop_rast <- pop_rast + r
}
message(sprintf("  Total DRC population : %.0f",
                global(pop_rast, "sum", na.rm = TRUE)[[1]]))

# Parse age lower bound from each _t_ file for age breakdown storage
# Format: cod_t_{age}_2025_CN_100m_R2025A_v1.tif
age_layers <- lapply(files_t, function(f) {
  bn  <- basename(f)
  age <- as.integer(sub("cod_t_(\\d+)_.*", "\\1", bn))
  list(f = f, age = age)
})

# ==============================================================================
# [Section 3] Create output directories
# ==============================================================================

message("\n=== Section 3: Setting up output directories ===")

dir.create("output/popdata", showWarnings = FALSE, recursive = TRUE)
dir.create("figure/C1_p2",   showWarnings = FALSE, recursive = TRUE)
message("  output/popdata/  — population RDS files")
message("  figure/C1_p2/    — population density maps")

# ==============================================================================
# [Section 4] Process each Level 2 territory
# ==============================================================================

message("\n=== Section 4: Processing territories ===")

n_total   <- nrow(cod2_sf)
t_start   <- proc.time()[["elapsed"]]
n_saved   <- 0L
n_skipped <- 0L

for (i in seq_len(n_total)) {

  terr  <- cod2_sf[i, ]
  name1 <- terr$NAME_1  # Province
  name2 <- terr$NAME_2  # Territory

  # Safe base name for file outputs
  safe_name <- sprintf("%s_%s",
                       gsub("[^A-Za-z0-9]", "_", name1),
                       gsub("[^A-Za-z0-9]", "_", name2))

  rds_path <- file.path("output/popdata", paste0(safe_name, "_pop.rds"))
  png_path <- file.path("figure/C1_p2",   paste0(safe_name, ".png"))

  # Skip if both outputs already exist
  if (file.exists(rds_path) && file.exists(png_path)) {
    n_skipped <- n_skipped + 1L
    next
  }

  # Crop total population raster to territory bounding box
  bbox     <- st_bbox(terr)
  terr_ext <- ext(bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"])

  pop_crop <- tryCatch(crop(pop_rast, terr_ext), error = function(e) NULL)
  if (is.null(pop_crop)) {
    message(sprintf("  [%d/%d] SKIP (crop failed): %s - %s", i, n_total, name1, name2))
    n_skipped <- n_skipped + 1L; next
  }

  # Mask to exact territory boundary
  pop_mask <- tryCatch(mask(pop_crop, vect(terr)), error = function(e) NULL)
  if (is.null(pop_mask)) {
    message(sprintf("  [%d/%d] SKIP (mask failed): %s - %s", i, n_total, name1, name2))
    n_skipped <- n_skipped + 1L; next
  }

  # Convert to dataframe (total population)
  pop_df_total <- as.data.frame(pop_mask, xy = TRUE) |>
    setNames(c("x", "y", "pop_total")) |>
    filter(!is.na(pop_total), pop_total > 0)

  if (nrow(pop_df_total) == 0) {
    n_skipped <- n_skipped + 1L; next
  }

  # Build age breakdown per cell for this territory
  # Each _t_ age layer cropped + masked, stored as age_lower + count
  age_list <- lapply(age_layers, function(lyr) {
    r_crop <- tryCatch({
      r <- rast(lyr$f)
      crop(mask(crop(r, terr_ext), vect(terr)), terr_ext)
    }, error = function(e) NULL)
    if (is.null(r_crop)) return(NULL)

    df <- as.data.frame(r_crop, xy = TRUE) |>
      setNames(c("x", "y", "count")) |>
      filter(!is.na(count)) |>
      mutate(age_lower = lyr$age)
    df
  })

  # Combine all age layers
  pop_df_age <- bind_rows(Filter(Negate(is.null), age_list))

  # Save RDS: list with total + age breakdown + territory metadata
  pop_save <- list(
    name1      = name1,
    name2      = name2,
    safe_name  = safe_name,
    total_pop  = sum(pop_df_total$pop_total),
    pop_total  = pop_df_total,   # x, y, pop_total
    pop_age    = pop_df_age      # x, y, count, age_lower
  )

  if (!file.exists(rds_path)) {
    saveRDS(pop_save, rds_path)
  }

  # Save PNG map
  if (!file.exists(png_path)) {
    p <- ggplot() +
      geom_raster(data = pop_df_total,
                  aes(x = x, y = y, fill = log1p(pop_total))) +
      scale_fill_viridis_c(option = "magma") +
      geom_sf(data = terr, fill = NA,
              color = "black", linewidth = 0.8, alpha = 0.7) +
      coord_sf() +
      labs(title    = sprintf("%s — %s", name1, name2),
           subtitle = sprintf("Population: %.0f", sum(pop_df_total$pop_total)),
           x = NULL, y = NULL) +
      theme_void() +
      theme(legend.position = "none",
            plot.title      = element_text(size = 11, face = "bold"),
            plot.subtitle   = element_text(size = 9, color = "grey40"),
            plot.margin     = margin(5, 5, 5, 5))

    ggsave(png_path, plot = p, width = 6, height = 5, dpi = 120)
  }

  n_saved <- n_saved + 1L

  # Progress every 10 territories
  if (i %% 10 == 0 || i == n_total) {
    elapsed <- proc.time()[["elapsed"]] - t_start
    rate    <- n_saved / max(elapsed, 0.1)
    eta     <- round((n_total - i) / max(rate, 0.01) / 60, 1)
    message(sprintf("  [%d/%d] saved: %d | skipped: %d | %.1f terr/sec | ETA %.1f min",
                    i, n_total, n_saved, n_skipped, rate, eta))
  }
}

# ==============================================================================
# [Done] Summary
# ==============================================================================

elapsed_total <- round(proc.time()[["elapsed"]] - t_start, 1)
message("\n=== Done ===")
message(sprintf("  Territories processed : %d", n_total))
message(sprintf("  Saved                 : %d", n_saved))
message(sprintf("  Skipped               : %d", n_skipped))
message(sprintf("  Total time            : %.1f sec (%.1f min)",
                elapsed_total, elapsed_total / 60))
message("  RDS → output/popdata/")
message("  PNG → figure/C1_p2/")

# ==============================================================================
# [Section 5] Metro regions — Ituri + Kivu (Nord + Sud)
# ==============================================================================

message("\n=== Section 5: Metro regions ===")

dir.create("output/popdata/metro", showWarnings = FALSE, recursive = TRUE)
dir.create("figure/C1_p2/metro",  showWarnings = FALSE, recursive = TRUE)

pop_files_lv2 <- list.files("output/popdata",
                            pattern = "_pop\\.rds$",
                            full.names = TRUE)

library(ggspatial)
# Metro definitions
metros <- list(
  Ituri = list(
    tag      = "metro_Ituri",
    name1    = "Ituri",
    pattern  = "^Ituri_"
  ),
  Kivu = list(
    tag      = "metro_Kivu",
    name1    = "Nord-Kivu + Sud-Kivu",
    pattern  = "^Nord_Kivu_|^Sud_Kivu_"
  )
)

for (metro in metros) {
  rds_path <- file.path("output/popdata/metro",
                        paste0(metro$tag, "_pop.rds"))
  png_path <- file.path("figure/C1_p2/metro",
                        paste0(metro$tag, ".png"))

  # if (file.exists(rds_path) && file.exists(png_path)) {
  #   message(sprintf("  SKIP (exists): %s", metro$tag)); next
  # }

  message(sprintf("  Building: %s", metro$tag))

  # Load matching territory files
  m_files <- pop_files_lv2[grepl(metro$pattern, basename(pop_files_lv2))]
  message(sprintf("    Territories: %d", length(m_files)))

  # Merge pop_total and pop_age across territories
  pop_total_list <- lapply(m_files, function(f) readRDS(f)$pop_total)
  pop_age_list   <- lapply(m_files, function(f) readRDS(f)$pop_age)

  pop_df_total <- bind_rows(pop_total_list)
  pop_df_age   <- bind_rows(pop_age_list)

  total_pop <- sum(pop_df_total$pop_total)
  n_cells   <- nrow(pop_df_total)

  message(sprintf("    Population: %.0f | Cells: %d", total_pop, n_cells))

  # Save RDS
  if (!file.exists(rds_path)) {
    saveRDS(list(
      name1     = metro$name1,
      name2     = metro$tag,
      safe_name = metro$tag,
      total_pop = total_pop,
      pop_total = pop_df_total,
      pop_age   = pop_df_age
    ), rds_path)
    message(sprintf("    Saved RDS: %s", basename(rds_path)))
  }

  # Save PNG map
  # if (!file.exists(png_path)) {
    # Load province boundaries for overlay
    prov_names <- unique(unlist(lapply(m_files, function(f) readRDS(f)$name1)))
    prov_sf    <- cod2_sf[cod2_sf$NAME_1 %in% prov_names, ]


    p <- ggplot() +
      geom_raster(data = pop_df_total,
                  aes(x = x, y = y, fill = log1p(pop_total))) +
      scale_fill_viridis_c(option = "magma") +
      geom_sf(data = prov_sf, fill = NA,
              color = "black", linewidth = 0.3, alpha = 0.4) +
      annotation_scale(location = "bl", width_hint = 0.3,
                       style = "ticks", text_col = "black",
                       line_col = "black") +
      annotation_north_arrow(location = "tl",
                             which_north = "true",
                             style = north_arrow_fancy_orienteering(
                               fill = c("grey", "black"),
                               text_col = "black"
                             )) +
      coord_sf() +
      labs(title    = metro$name1,
           subtitle = sprintf("Population: %.0f | Cells: %d",
                              total_pop, n_cells),
           x = NULL, y = NULL) +
      theme_void() +
      theme(legend.position = "none",
            plot.title      = element_text(size = 13, face = "bold"),
            plot.subtitle   = element_text(size = 10, color = "grey40"),
            plot.margin     = margin(5, 5, 5, 5))

    ggsave(png_path, plot = p, width = 5, height = 5, dpi = 500)
    message(sprintf("    Saved PNG: %s", basename(png_path)))
  # }
}

message("\n=== Metro Done ===")
