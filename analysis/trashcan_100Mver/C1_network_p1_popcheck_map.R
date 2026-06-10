# ==============================================================================
# C1_network_p1_popcheck_map.R
# Purpose:
#   Load WorldPop R2025A 2025 population data for DRC
#   Plot full DRC map + Nord-Kivu zoom with Level 2 boundaries
#
# Inputs:
#   data/worldpop/   : WorldPop R2025A GeoTIFF files (_t_ total files only)
#   data/shpmap/     : GADM Level 1 and Level 2 shapefiles
#
# Outputs:
#   figure/C1_p1/drc_map.png
#   figure/C1_p1/nordkivu_map.png
# ==============================================================================

library(terra)
library(sf)
library(dplyr)
library(ggplot2)

# ==============================================================================
# [Section 1] Load shapefiles
# ==============================================================================

message("=== Section 1: Loading shapefiles ===")

cod1_sf <- st_as_sf(vect("data/shpmap/gadm41_COD_1.shp"))
cod2_sf <- st_as_sf(vect("data/shpmap/gadm41_COD_2.shp"))

# Extract Nord-Kivu boundaries
nk_l1 <- cod1_sf[cod1_sf$NAME_1 == "Nord-Kivu", ]
nk_l2 <- cod2_sf[cod2_sf$NAME_1 == "Nord-Kivu", ]

message(sprintf("  Provinces (Level 1)   : %d", nrow(cod1_sf)))
message(sprintf("  Territories (Level 2) : %d", nrow(cod2_sf)))
message(sprintf("  Nord-Kivu territories : %d", nrow(nk_l2)))

# ==============================================================================
# [Section 2] Load and sum WorldPop R2025A population rasters
# Use _t_ files (age-stratified totals = female + male already combined)
# Avoids loading 40 f/m files — only 20 _t_ files needed
# ==============================================================================

message("\n=== Section 2: Loading population rasters ===")

files_all <- list.files("data/worldpop", pattern = "\\.tif$", full.names = TRUE)
files_t   <- files_all[grepl("/cod_t_", files_all)]  # lowercase _t_ only
message(sprintf("  Using %d files (age-stratified totals)", length(files_t)))

# Sum all age layers to get total population raster
pop_rast <- NULL
for (f in files_t) {
  r <- rast(f)
  if (is.null(pop_rast)) pop_rast <- r else pop_rast <- pop_rast + r
}
message(sprintf("  Total DRC population : %.0f",
                global(pop_rast, "sum", na.rm = TRUE)[[1]]))

# ==============================================================================
# [Section 3] Convert raster to dataframe
# ==============================================================================

message("\n=== Section 3: Converting raster to dataframe ===")

pop_df <- as.data.frame(pop_rast, xy = TRUE) |>
  setNames(c("x", "y", "pop")) |>
  filter(!is.na(pop), pop > 0)

message(sprintf("  Populated cells : %d", nrow(pop_df)))
message(sprintf("  Total pop       : %.0f", sum(pop_df$pop)))

# Subset to Nord-Kivu bounding box
nk_bbox   <- st_bbox(nk_l1)
pop_nk_df <- pop_df |>
  filter(x >= nk_bbox["xmin"], x <= nk_bbox["xmax"],
         y >= nk_bbox["ymin"], y <= nk_bbox["ymax"])

message(sprintf("  Nord-Kivu pop   : %.0f", sum(pop_nk_df$pop)))

# ==============================================================================
# [Section 4] Plot
# ==============================================================================

message("\n=== Section 4: Plotting ===")

# Full DRC map with Nord-Kivu highlighted
p_drc <- ggplot() +
  geom_raster(data = pop_df,
              aes(x = x, y = y, fill = log1p(pop))) +
  scale_fill_viridis_c(option = "magma") +
  geom_sf(data = cod2_sf, fill = NA,
          color = "black", linewidth = 0.5, alpha = 0.4) +
  geom_sf(data = nk_l1, fill = NA,
          color = "black", linewidth = 1.0, alpha = 0.8) +
  coord_sf() +
  labs(
       x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "none")

# Nord-Kivu zoom
p_nk <- ggplot() +
  geom_raster(data = pop_nk_df,
              aes(x = x, y = y, fill = log1p(pop))) +
  scale_fill_viridis_c(option = "magma") +
  geom_sf(data = nk_l2, fill = NA,
          color = "black", linewidth = 0.5, alpha = 0.4) +
  geom_sf(data = nk_l1, fill = NA,
          color = "black", linewidth = 1.2, alpha = 0.8) +
  coord_sf() +
  labs(
       x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "none")

# ==============================================================================
# [Section 5] Save separately
# ==============================================================================

message("\n=== Section 5: Saving ===")

dir.create("figure/C1_p1", showWarnings = FALSE, recursive = TRUE)

ggsave("figure/C1_p1/drc_map.png",
       plot = p_drc, width = 10, height = 10, dpi = 500)
message("  Saved: figure/C1_p1/drc_map.png")

ggsave("figure/C1_p1/nordkivu_map.png",
       plot = p_nk, width = 8, height = 10, dpi = 500)
message("  Saved: figure/C1_p1/nordkivu_map.png")

message("\nDone!")
