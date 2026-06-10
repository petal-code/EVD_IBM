# ==============================================================================
# C2_sim_p1_gifgen.R
# Purpose:
#   Animated GIF of Ebola IBM simulation result
#   - WorldPop raster background
#   - Epidemic curve by route
#   - SIR curve
#   - Age distribution of infections
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(gifski)
library(terra)

# ==============================================================================
# [Configuration]
# ==============================================================================

case_tag   <- "case1_1M"
sim_result <- readRDS(sprintf("output/sim/%s_sim_result.rds", case_tag))
nodes      <- readRDS(sprintf("output/network/%s_nodes.rds",  case_tag))

fig_dir    <- "figure/C2_gif"
frames_dir <- file.path(fig_dir, "frames")
dir.create(frames_dir, showWarnings = FALSE, recursive = TRUE)

t_max               <- 200    # max days to animate
t_step              <- 5      # days per frame
fps                 <- 10
infectious_duration <- 14     # approx infectious window for map coloring
frame_width         <- 8
frame_height        <- 4

# ==============================================================================
# [Load + prepare data]
# ==============================================================================

inf_df <- sim_result$infected
N      <- nrow(nodes)

# Route labels
ctype_labels <- c("0"="Seed", "1"="Household", "2"="Community",
                  "3"="Hospital", "4"="Funeral")
route_colors <- c("Seed"      = "#888888",
                  "Household" = "#4472C4",
                  "Community" = "#ED7D31",
                  "Hospital"  = "#9B59B6",
                  "Funeral"   = "#E74C3C")

age_labels_16 <- c(paste0(seq(0, 70, 5), "-", seq(4, 74, 5)), "75+")

# Add coordinates + labels
set.seed(42)
inf_df <- inf_df %>%
  left_join(nodes %>% select(person_id, x, y, cell_id), by="person_id") %>%
  mutate(
    route     = factor(ctype_labels[as.character(contact_type)],
                       levels=names(route_colors)),
    day       = floor(time_infection),
    age_label = factor(age_labels_16[pmin(age_group, 16)],
                       levels=age_labels_16),
    # Jitter within cell (~500m)
    lon_j = x + rnorm(n(), 0, 0.003),
    lat_j = y + rnorm(n(), 0, 0.003)
  ) %>%
  filter(time_infection <= t_max)

# Node coordinates for background dots
node_coords <- nodes %>%
  group_by(x, y) %>%
  summarise(n = n(), .groups = "drop") %>%
  slice_sample(n = min(5000, nrow(.)))  # Sample for speed

# WorldPop raster background
pop_files <- list.files("data/worldpop/DRC_1km", pattern="\\.tif$", full.names=TRUE)
pop_files <- pop_files[grepl("/cod_t_", pop_files)]
pop_rast  <- NULL
for (f in pop_files) {
  r <- rast(f)
  if (is.null(pop_rast)) pop_rast <- r else pop_rast <- pop_rast + r
}

# Crop to case bbox
x_range <- range(nodes$x)
y_range <- range(nodes$y)
pad <- 0.02
pop_crop <- as.data.frame(
  crop(pop_rast, ext(x_range[1]-pad, x_range[2]+pad,
                     y_range[1]-pad, y_range[2]+pad)),
  xy = TRUE
) %>% setNames(c("x","y","pop")) %>% filter(!is.na(pop), pop > 0)

# Fixed y limits
epi_all  <- inf_df %>% count(day, route)
ylim_epi <- max(tapply(epi_all$n, epi_all$day, sum), na.rm=TRUE) * 1.15

age_all  <- inf_df %>% count(age_label, .drop=FALSE)
ylim_age <- max(age_all$n, na.rm=TRUE) * 1.15

x_breaks <- pretty(0:t_max, n=8)
x_breaks <- x_breaks[x_breaks >= 0 & x_breaks <= t_max]

# ==============================================================================
# [Theme]
# ==============================================================================

theme_pub <- function(base_size=10) {
  theme_minimal(base_size=base_size) +
    theme(panel.grid.minor  = element_blank(),
          axis.text         = element_text(size=base_size),
          axis.title        = element_text(size=base_size),
          plot.title        = element_text(size=base_size, face="bold", hjust=0.5),
          legend.text       = element_text(size=base_size))
}

# ==============================================================================
# [Render frames]
# ==============================================================================

t_seq <- seq(0, t_max, by=t_step)
message(sprintf("Rendering %d frames (single core)...", length(t_seq)))
t0 <- proc.time()[["elapsed"]]

frame_paths <- character(length(t_seq))

for (fi in seq_along(t_seq)) {

  t_val <- t_seq[fi]
  sub   <- inf_df %>% filter(time_infection <= t_val)

  # ── Map status ──────────────────────────────────────────────────
  map_df <- sub %>%
    filter(!is.na(lon_j), !is.na(lat_j))

  p_map <- ggplot() +
    geom_raster(data=pop_crop, aes(x=x, y=y, fill=log1p(pop)), alpha=0.6) +
    scale_fill_viridis_c(option="magma", guide="none") +
    geom_point(data=node_coords, aes(x=x, y=y),
               color="grey70", size=0.15, alpha=0.3) +
    {if (nrow(map_df) > 0)
      geom_point(data=map_df,
                 aes(x=lon_j, y=lat_j, color=route),
                 size=0.5, alpha=0.8)} +
    scale_color_manual(values=route_colors, name=NULL, drop=FALSE) +
    coord_fixed(
      xlim=c(x_range[1]-pad, x_range[2]+pad),
      ylim=c(y_range[1]-pad, y_range[2]+pad)
    ) +
    labs(title=sprintf("%s | Day %.0f | N infected: %d",
                       case_tag, t_val, nrow(sub))) +
    theme_void() +
    theme(plot.title        = element_text(size=10, face="bold", hjust=0.5),
          legend.position   = c(0.97, 0.97),
          legend.justification = c("right","top"),
          legend.background = element_rect(fill=alpha("white",0.7), color=NA),
          legend.text       = element_text(size=9))

  # ── Epidemic curve ───────────────────────────────────────────────
  epi <- sub %>% count(day, route)

  p_epi <- ggplot() +
    {if (nrow(epi) > 0)
      geom_col(data=epi, aes(x=day, y=n, fill=route), width=0.8)} +
    scale_fill_manual(values=route_colors, drop=FALSE) +
    scale_x_continuous(limits=c(-0.5, t_max+0.5), breaks=x_breaks) +
    scale_y_continuous(limits=c(0, ylim_epi),
                       expand=expansion(mult=c(0,0))) +
    labs(x=NULL, y="New cases", title="Epidemic curve") +
    theme_pub() +
    theme(legend.position="none", panel.grid.major.x=element_blank())

  # ── Age distribution ─────────────────────────────────────────────
  age_cum <- sub %>%
    count(age_label, route, .drop=FALSE) %>%
    complete(
      age_label = factor(age_labels_16, levels=age_labels_16),
      route     = factor(names(route_colors), levels=names(route_colors)),
      fill = list(n=0L)
    )

  p_age <- ggplot() +
    {if (nrow(age_cum) > 0)
      geom_col(data=age_cum, aes(x=age_label, y=n, fill=route), width=0.75)} +
    scale_fill_manual(values=route_colors, drop=FALSE) +
    scale_x_discrete(limits=age_labels_16) +
    scale_y_continuous(limits=c(0, ylim_age),
                       expand=expansion(mult=c(0,0))) +
    labs(x=NULL, y="Cumulative cases", title="Age distribution") +
    theme_pub() +
    theme(axis.text.x     = element_text(angle=45, hjust=1, size=8),
          legend.position = "none",
          panel.grid.major.x = element_blank())

  # ── Assemble ─────────────────────────────────────────────────────
  p_full <- p_map + p_epi + p_age + plot_layout(widths=c(2, 1.2, 1))

  frame_path <- file.path(frames_dir, sprintf("frame_%04d.png", fi))
  ggsave(frame_path, p_full,
         width=frame_width, height=frame_height, dpi=120)
  frame_paths[fi] <- frame_path

  if (fi %% 20 == 0)
    message(sprintf("  Frame %d/%d (%.1f sec)", fi, length(t_seq),
                    proc.time()[["elapsed"]] - t0))
}

# ==============================================================================
# [Assemble GIF]
# ==============================================================================

message("Assembling GIF...")
gif_path <- file.path(fig_dir, sprintf("%s_epidemic.gif", case_tag))
gifski(png_files=frame_paths, gif_file=gif_path,
       width=1440, height=720, delay=1/fps)
message(sprintf("Saved: %s", gif_path))
