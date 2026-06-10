library(terra)

# ── Load travel time rasters ───────────────────────────────────
motorized <- rast("data/2020_motorized_travel_time_to_healthcare.geotiff")
walking   <- rast("data/2020_walking_only_travel_time_to_healthcare.geotiff")

# ── Check metadata ─────────────────────────────────────────────
cat("=== Motorized ===\n")
print(motorized)

cat("\n=== Walking ===\n")
print(walking)

# ── Crop to Kinshasa extent ────────────────────────────────────
# Get Kinshasa bounding box from population data
kin_ext <- ext(
  min(pop$indiv_x) - 0.05,
  max(pop$indiv_x) + 0.05,
  min(pop$indiv_y) - 0.05,
  max(pop$indiv_y) + 0.05
)
cat("Kinshasa extent:\n")
print(kin_ext)

# Crop both rasters
motorized_kin <- crop(motorized, kin_ext)
walking_kin   <- crop(walking,   kin_ext)

cat("\n=== Motorized (Kinshasa) ===\n")
print(motorized_kin)
cat(sprintf("Resolution: ~%.0fm x ~%.0fm\n",
            res(motorized_kin)[1] * 111000,
            res(motorized_kin)[2] * 111000))

# ── Quick summary of travel times ─────────────────────────────
cat("\nMotorized travel time (min):\n")
print(summary(values(motorized_kin, na.rm = TRUE)))

cat("\nWalking travel time (min):\n")
print(summary(values(walking_kin, na.rm = TRUE)))

# ── Plot both side by side ─────────────────────────────────────
par(mfrow = c(1, 2))
plot(motorized_kin,
     main   = "Motorized travel time (min)\nKinshasa",
     col    = hcl.colors(100, "RdYlGn", rev = TRUE))

plot(walking_kin,
     main   = "Walking travel time (min)\nKinshasa",
     col    = hcl.colors(100, "RdYlGn", rev = TRUE))
par(mfrow = c(1, 1))




# ── Filter hospital-level facilities only ─────────────────────
hf_hospital <- hf_kin %>%
  filter(facility_group == "Hospital")

cat(sprintf("Total Kinshasa facilities : %d\n", nrow(hf_kin)))
cat(sprintf("Hospital-level only       : %d\n", nrow(hf_hospital)))
cat("\nBreakdown:\n")
print(table(hf_hospital$esstype))

# ── Plot hospitals only on population density ──────────────────
ggplot() +
  stat_density_2d(data    = pop_sample,
                  aes(x = indiv_x, y = indiv_y, fill = after_stat(density)),
                  geom    = "raster",
                  contour = FALSE,
                  alpha   = 0.8) +
  scale_fill_viridis_c(option = "magma", name = "Pop density") +

  geom_point(data  = hf_hospital,
             aes(x = lon, y = lat, color = esstype),
             size  = 0.5,
             alpha = 0.4) +
  scale_color_manual(
    name   = "Hospital type",
    values = c("Hôpital"                      = "#e41a1c",
               "Hôpital Général de Référence" = "#ff7f00",
               "Centre Hopitalier"            = "#4daf4a")
  ) +
  coord_fixed() +
  labs(title    = "Kinshasa: Population density + Hospitals",
       subtitle = sprintf("%d hospitals | %.0fk individuals sampled",
                          nrow(hf_hospital), nrow(pop_sample) / 1000),
       x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(legend.position = "right")

ggsave("output/household/Kinshasa_pop_hospitals.png",
       width = 14, height = 10, dpi = 150)
