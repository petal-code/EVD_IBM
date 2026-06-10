library(terra)
library(ggplot2)

# Load level 2 shapefile
cod2 <- vect("data/shpmap/gadm41_COD_2.shp")

# Use total population raster (cod_t files summed)
# Already have files_fm — use _t_ files for cleaner total
files_t <- files[grepl("/cod_t_", files)]

cat("Loading total population raster...\n")
pop_total_drc <- NULL
for (f in files_t) {
  r <- rast(f)
  if (is.null(pop_total_drc)) pop_total_drc <- r else pop_total_drc <- pop_total_drc + r
}
cat(sprintf("Total DRC pop: %.0f\n", global(pop_total_drc, "sum", na.rm=TRUE)[[1]]))

# Convert to dataframe for plotting
cat("Converting to dataframe (this may take a moment)...\n")
pop_df_drc <- as.data.frame(pop_total_drc, xy = TRUE) |>
  setNames(c("x", "y", "pop")) |>
  filter(!is.na(pop), pop > 0)

# Convert shapefile to dataframe for ggplot
cod2_df <- as.data.frame(geom(cod2))

# Plot
cat("Plotting...\n")
library(sf)

# Convert terra SpatVector to sf for ggplot2
cod2_sf <- st_as_sf(cod2)

# Plot
ggplot() +
  geom_raster(data = pop_df_drc,
              aes(x = x, y = y, fill = log1p(pop))) +
  geom_sf(data = cod2_sf, fill = NA,
          color = "white", linewidth = 0.1) +
  scale_fill_viridis_c(option = "magma",
                       name   = "log(pop+1)",
                       labels = function(x) round(expm1(x), 1)) +
  coord_sf() +
  labs(title    = "DRC Population density (WorldPop R2025A 2025)",
       subtitle = sprintf("Total: %.0f | Level 2 boundaries overlay",
                          sum(pop_df_drc$pop)),
       x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(legend.position = "right")

ggsave("output/drc_pop_level2_map.png", width = 14, height = 10, dpi = 150)

# ggsave("output/drc_pop_level2_map.png", width = 14, height = 10, dpi = 150)
cat("Saved to output/drc_pop_level2_map.png\n")


ggplot() +
  geom_raster(data = pop_df_drc,
              aes(x = x, y = y, fill = log1p(pop))) +
  scale_fill_viridis_c(option = "magma",
                       name   = "log(pop+1)",
                       labels = function(x) round(expm1(x), 1)) +
  geom_sf(data = cod2_sf, fill = NA,
          color = "black", linewidth = 0.2, alpha = 0.3) +
  coord_sf() +
  labs(title    = "DRC Population density (WorldPop R2025A 2025)",
       subtitle = sprintf("Total: %.0f | Level 2 boundaries overlay",
                          sum(pop_df_drc$pop)),
       x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position  = "none",
        plot.title       = element_text(size = 14, face = "bold"),
        plot.subtitle    = element_text(size = 10))

ggsave("output/drc_pop_level2_map.png", width = 14, height = 10, dpi = 150)



p_nk <- ggplot() +
  geom_raster(data = pop_nk_df,
              aes(x = x, y = y, fill = log1p(pop))) +
  scale_fill_viridis_c(option = "magma") +
  geom_sf(data = nk_l2, fill = NA,
          color = "black", linewidth = 0.3, alpha = 0.4) +
  geom_sf(data = nk_l1, fill = NA,
          color = "white", linewidth = 1.0) +
  coord_sf() +
  labs(title    = "Nord-Kivu (Level 2 boundaries)",
       subtitle = sprintf("Population: %.0f", sum(pop_nk_df$pop)),
       x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "none",
        plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9))

# Combine
p_combined <- p_drc + p_nk +
  plot_annotation(
    title    = "DRC Population density (WorldPop R2025A 2025)",
    subtitle = sprintf("Total DRC: %.0f", sum(pop_df_drc$pop))
  )

ggsave("output/drc_nordkivu_map.png",
       plot = p_combined, width = 16, height = 9, dpi = 150)
cat("Saved to output/drc_nordkivu_map.png\n")


library(sf)
library(ggplot2)
library(patchwork)

# Filter Nord-Kivu from level 1 and level 2
cod1_sf <- st_as_sf(cod1)
cod2_sf <- st_as_sf(cod2)

nk_l1 <- cod1_sf[cod1_sf$NAME_1 == "Nord-Kivu", ]
nk_l2 <- cod2_sf[cod2_sf$NAME_1 == "Nord-Kivu", ]

# Crop population raster to Nord-Kivu extent
nk_ext  <- ext(vect(nk_l1))
pop_nk_df <- pop_df_drc |>
  filter(x >= nk_ext[1], x <= nk_ext[2],
         y >= nk_ext[3], y <= nk_ext[4])

# Full DRC map
p_drc <- ggplot() +
  geom_raster(data = pop_df_drc,
              aes(x = x, y = y, fill = log1p(pop))) +
  scale_fill_viridis_c(option = "magma") +
  geom_sf(data = cod2_sf, fill = NA,
          color = "black", linewidth = 0.2, alpha = 0.3) +
  # Highlight Nord-Kivu boundary
  geom_sf(data = nk_l1, fill = NA,
          color = "white", linewidth = 0.8) +
  coord_sf() +
  labs(title = "DRC", x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, face = "bold"))

# Nord-Kivu zoom
p_nk <- ggplot() +
  geom_raster(data = pop_nk_df,
              aes(x = x, y = y, fill = log1p(pop))) +
  scale_fill_viridis_c(option = "magma") +
  geom_sf(data = nk_l2, fill = NA,
          color = "black", linewidth = 0.3, alpha = 0.4) +
  geom_sf(data = nk_l1, fill = NA,
          color = "white", linewidth = 1.0) +
  coord_sf() +
  labs(title    = "Nord-Kivu (Level 2 boundaries)",
       subtitle = sprintf("Population: %.0f", sum(pop_nk_df$pop)),
       x = NULL, y = NULL) +
  theme_void() +
  theme(legend.position = "none",
        plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9))

# # Combine
# p_combined <- p_drc + p_nk +
#   plot_annotation(
#     title    = "DRC Population density (WorldPop R2025A 2025)",
#     subtitle = sprintf("Total DRC: %.0f", sum(pop_df_drc$pop))
#   )

ggsave("output/drc_nordkivu_map.png",
       plot = p_nk, width = 16, height = 9, dpi = 150)
cat("Saved to output/drc_nordkivu_map.png\n")
