# ==============================================================================
# C1_network_p6_validation.R
# Purpose:
#   Visualize network for one territory:
#     1. Population density map with hospitals
#     2. Household locations (jittered within cell)
#     3. Layer 1 (household) edges sample
#     4. Layer 2 (community) edges sample
#     5. Layer 3 HCW locations + edges sample
# ==============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# ==============================================================================
# [Configuration]
# ==============================================================================

network_dir <- "output/network"
synpop_dir  <- "output/household"
hf_path     <- "data/COD_GRID3_health_facilities_v8.csv"
fig_dir     <- "figure/C1_p6_validation"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Change index to check different territory
territory_idx <- 1L

# Jitter amount within cell (~100m cell = 0.0009 deg)
jitter_deg <- 0.0004

# Max edges to show per layer (for readability)
max_edges_show <- 500L

# ==============================================================================
# [Section 1] Load data
# ==============================================================================

# Find matching files
net_files <- list.files(network_dir, pattern = "_nodes\\.rds$",
                        full.names = TRUE)
cat(sprintf("Available territories:\n"))
for (i in seq_along(net_files))
  cat(sprintf("  [%d] %s\n", i,
              sub("_nodes\\.rds$", "", basename(net_files[i]))))

tag <- sub("_nodes\\.rds$", "", basename(net_files[territory_idx]))
cat(sprintf("\nLoading: %s\n", tag))

nodes    <- readRDS(file.path(network_dir, sprintf("%s_nodes.rds",            tag)))
layer1   <- readRDS(file.path(network_dir, sprintf("%s_layer1_household.rds", tag)))
layer2   <- readRDS(file.path(network_dir, sprintf("%s_layer2_community.rds", tag)))
layer3h  <- readRDS(file.path(network_dir, sprintf("%s_layer3_hcw_edges.rds", tag)))
layer3a  <- readRDS(file.path(network_dir, sprintf("%s_layer3_admission.rds", tag)))

cat(sprintf("  Nodes       : %d\n", nrow(nodes)))
cat(sprintf("  HCWs        : %d\n", sum(nodes$is_hcw)))
cat(sprintf("  Layer 1 HH  : %d edges\n", nrow(layer1)))
cat(sprintf("  Layer 2 comm: %d edges\n", nrow(layer2)))
cat(sprintf("  Layer 3 HCW : %d edges\n", nrow(layer3h)))

# Load hospitals for this territory
hf_all <- read.csv(hf_path)
hf_terr <- hf_all %>%
  filter(
    esstype %in% c("Hôpital", "Hôpital Général de Référence", "Centre Hopitalier"),
    !is.na(lon), !is.na(lat)
  ) %>%
  filter(OBJECTID %in% unique(nodes$hospital_id))

cat(sprintf("  Hospitals   : %d\n", nrow(hf_terr)))

# ==============================================================================
# [Section 2] Add jitter to node coordinates
# ==============================================================================

set.seed(42)
nodes_j <- nodes %>%
  mutate(
    xj = x + runif(n(), -jitter_deg, jitter_deg),
    yj = y + runif(n(), -jitter_deg, jitter_deg)
  )

# Quick lookup: person_id → jittered coords
coord_lookup <- nodes_j %>% select(person_id, xj, yj)

# ==============================================================================
# [Section 3] Map 1 — Population density + hospitals
# ==============================================================================

cell_pop <- nodes %>%
  group_by(x, y) %>%
  summarise(n = n(), .groups = "drop")

p_density <- ggplot() +
  geom_raster(data = cell_pop,
              aes(x = x, y = y, fill = log1p(n))) +
  scale_fill_viridis_c(option = "magma", name = "log(pop+1)") +
  geom_point(data = hf_terr,
             aes(x = lon, y = lat),
             color = "cyan", size = 3, shape = 17) +
  geom_point(data = hf_terr,
             aes(x = lon, y = lat),
             color = "white", size = 1.2, shape = 17) +
  coord_fixed() +
  labs(title    = sprintf("%s", tag),
       subtitle = sprintf("Population: %d | Hospitals: %d (▲)",
                          nrow(nodes), nrow(hf_terr)),
       x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        legend.position = "right")

# ==============================================================================
# [Section 4] Map 2 — Layer 1 household edges (sample)
# ==============================================================================

set.seed(42)
l1_sample <- layer1 %>%
  slice_sample(n = min(max_edges_show, nrow(layer1)))

l1_edges <- l1_sample %>%
  left_join(coord_lookup, by = c("from" = "person_id")) %>%
  rename(x1 = xj, y1 = yj) %>%
  left_join(coord_lookup, by = c("to" = "person_id")) %>%
  rename(x2 = xj, y2 = yj)

p_layer1 <- ggplot() +
  geom_raster(data = cell_pop,
              aes(x = x, y = y, fill = log1p(n)), alpha = 0.4) +
  scale_fill_viridis_c(option = "magma", guide = "none") +
  geom_segment(data = l1_edges,
               aes(x=x1, y=y1, xend=x2, yend=y2),
               color = "steelblue", alpha = 0.4, linewidth = 0.3) +
  geom_point(data = nodes_j,
             aes(x = xj, y = yj),
             color = "white", size = 0.3, alpha = 0.5) +
  coord_fixed() +
  labs(title    = "Layer 1: Household edges",
       subtitle = sprintf("%d edges shown (of %d total)",
                          nrow(l1_sample), nrow(layer1)),
       x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"))

# ==============================================================================
# [Section 5] Map 3 — Layer 2 community edges (sample)
# ==============================================================================

set.seed(42)
l2_sample <- layer2 %>%
  slice_sample(n = min(max_edges_show, nrow(layer2)))

l2_edges <- l2_sample %>%
  left_join(coord_lookup, by = c("from" = "person_id")) %>%
  rename(x1 = xj, y1 = yj) %>%
  left_join(coord_lookup, by = c("to" = "person_id")) %>%
  rename(x2 = xj, y2 = yj)

p_layer2 <- ggplot() +
  geom_raster(data = cell_pop,
              aes(x = x, y = y, fill = log1p(n)), alpha = 0.4) +
  scale_fill_viridis_c(option = "magma", guide = "none") +
  geom_segment(data = l2_edges,
               aes(x=x1, y=y1, xend=x2, yend=y2),
               color = "coral", alpha = 0.3, linewidth = 0.3) +
  geom_point(data = nodes_j,
             aes(x = xj, y = yj),
             color = "white", size = 0.3, alpha = 0.5) +
  coord_fixed() +
  labs(title    = "Layer 2: Community edges",
       subtitle = sprintf("%d edges shown (of %d total)",
                          nrow(l2_sample), nrow(layer2)),
       x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"))

# ==============================================================================
# [Section 6] Map 4 — Layer 3 HCW locations + edges
# ==============================================================================

hcw_nodes_j <- nodes_j %>% filter(is_hcw)

set.seed(42)
l3_sample <- layer3h %>%
  slice_sample(n = min(max_edges_show, nrow(layer3h)))

l3_edges <- l3_sample %>%
  left_join(coord_lookup, by = c("from" = "person_id")) %>%
  rename(x1 = xj, y1 = yj) %>%
  left_join(coord_lookup, by = c("to" = "person_id")) %>%
  rename(x2 = xj, y2 = yj)

p_layer3 <- ggplot() +
  geom_raster(data = cell_pop,
              aes(x = x, y = y, fill = log1p(n)), alpha = 0.4) +
  scale_fill_viridis_c(option = "magma", guide = "none") +
  geom_segment(data = l3_edges,
               aes(x=x1, y=y1, xend=x2, yend=y2),
               color = "yellow", alpha = 0.5, linewidth = 0.4) +
  geom_point(data = hcw_nodes_j,
             aes(x = xj, y = yj),
             color = "yellow", size = 1.0, alpha = 0.8) +
  geom_point(data = hf_terr,
             aes(x = lon, y = lat),
             color = "cyan", size = 3, shape = 17) +
  coord_fixed() +
  labs(title    = "Layer 3: HCW network",
       subtitle = sprintf("HCWs: %d | HCW-HCW edges: %d | Hospitals: %d",
                          nrow(hcw_nodes_j), nrow(layer3h), nrow(hf_terr)),
       x = NULL, y = NULL) +
  theme_void() +
  theme(plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"))

# ==============================================================================
# [Section 7] Combine and save
# ==============================================================================

p_combined <- (p_density + p_layer1) / (p_layer2 + p_layer3) +
  plot_annotation(
    title    = sprintf("Network layers — %s", tag),
    subtitle = sprintf("Pop: %d | HH: %d | L1: %d | L2: %d | L3: %d edges",
                       nrow(nodes),
                       n_distinct(nodes$hh_id),
                       nrow(layer1), nrow(layer2), nrow(layer3h)),
    theme = theme(plot.title    = element_text(size = 13, face = "bold"),
                  plot.subtitle = element_text(size = 10, color = "grey40"))
  )

out_path <- file.path(fig_dir, sprintf("%s_network_layers.png", tag))
ggsave(out_path, plot = p_combined, width = 16, height = 12, dpi = 150)
cat(sprintf("\nSaved: %s\n", out_path))

# ==============================================================================
# [Section 8] Network summary stats
# ==============================================================================

cat("\n=== Network Summary ===\n")

# Degree distribution (Layer 2)
l2_degree <- c(layer2$from, layer2$to) %>%
  table() %>%
  as.data.frame() %>%
  setNames(c("person_id", "degree"))

cat(sprintf("Layer 2 degree: mean=%.1f | median=%.0f | max=%d\n",
            mean(l2_degree$degree),
            median(l2_degree$degree),
            max(l2_degree$degree)))

# HCW admission coverage
cat(sprintf("Admission lookup: %d non-HCW individuals\n", nrow(layer3a)))
no_hcw <- sum(sapply(layer3a$hcw_list, is.null))
cat(sprintf("  No HCW assigned: %d (%.1f%%)\n",
            no_hcw, 100 * no_hcw / nrow(layer3a)))

# Hospital coverage
cat(sprintf("\nHospitals in territory: %d\n", nrow(hf_terr)))
cat(sprintf("Unique hospitals assigned to nodes: %d\n",
            n_distinct(nodes$hospital_id)))


# ==============================================================================
# C1_network_p6_single_person.R
# Purpose:
#   Pick one individual and visualize all their network connections:
#     - Layer 1: household members
#     - Layer 2: community contacts (kernel-weighted distance effect visible)
#     - Layer 3: HCW connections (if applicable)
#   Shows kernel distance decay effect on community contact distances
# ==============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# ==============================================================================
# [Configuration] — Change these
# ==============================================================================

network_dir   <- "output/network"
hf_path       <- "data/COD_GRID3_health_facilities_v8.csv"
fig_dir       <- "figure/C1_p6_validation"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

for (mmm in 1:30)
{
territory_idx <- 1L    # Which territory to look at
person_id_sel <- mmm   # ← Change this (1 to 100) to pick a different individual
jitter_deg    <- 0.0004  # ~44m jitter within cell

# ==============================================================================
# [Section 1] Load data
# ==============================================================================

net_files <- list.files(network_dir, pattern = "_nodes\\.rds$",
                        full.names = TRUE)
tag <- sub("_nodes\\.rds$", "", basename(net_files[territory_idx]))
cat(sprintf("Territory: %s\n", tag))

nodes   <- readRDS(file.path(network_dir, sprintf("%s_nodes.rds",            tag)))
layer1  <- readRDS(file.path(network_dir, sprintf("%s_layer1_household.rds", tag)))
layer2  <- readRDS(file.path(network_dir, sprintf("%s_layer2_community.rds", tag)))
layer3h <- readRDS(file.path(network_dir, sprintf("%s_layer3_hcw_edges.rds", tag)))

# Clamp person_id to valid range
pid <- max(1L, min(as.integer(person_id_sel), 100L))
cat(sprintf("Selected person_id: %d\n", pid))

# Show selected person info
sel <- nodes[nodes$person_id == pid, ]
cat(sprintf("  Age: %d | Age group: %d | HCW: %s | HH: %s | Cell: %d\n",
            sel$age, sel$age_group,
            ifelse(sel$is_hcw, "YES", "no"),
            sel$hh_id, sel$cell_id))

# ==============================================================================
# [Section 2] Find all connections of selected person
# ==============================================================================

# Layer 1: household members
l1_contacts <- c(
  layer1$to[layer1$from == pid],
  layer1$from[layer1$to == pid]
)

# Layer 2: community contacts
l2_contacts <- c(
  layer2$to[layer2$from == pid],
  layer2$from[layer2$to == pid]
)

# Layer 3: HCW contacts
l3_contacts <- c(
  layer3h$to[layer3h$from == pid],
  layer3h$from[layer3h$to == pid]
)

cat(sprintf("  Layer 1 (HH)       : %d contacts\n", length(l1_contacts)))
cat(sprintf("  Layer 2 (community): %d contacts\n", length(l2_contacts)))
cat(sprintf("  Layer 3 (HCW)      : %d contacts\n", length(l3_contacts)))

# ==============================================================================
# [Section 3] Build plotting dataframe with jitter
# ==============================================================================

set.seed(42)
nodes_j <- nodes %>%
  mutate(
    xj = x + runif(n(), -jitter_deg, jitter_deg),
    yj = y + runif(n(), -jitter_deg, jitter_deg)
  )

# Selected person jittered coords
sel_j <- nodes_j[nodes_j$person_id == pid, ]

# All contacts combined
all_contacts <- unique(c(l1_contacts, l2_contacts, l3_contacts))
contact_nodes <- nodes_j %>%
  filter(person_id %in% all_contacts) %>%
  mutate(layer = case_when(
    person_id %in% l1_contacts & person_id %in% l3_contacts ~ "HH + HCW",
    person_id %in% l1_contacts ~ "Layer 1 (HH)",
    person_id %in% l3_contacts ~ "Layer 3 (HCW)",
    TRUE                        ~ "Layer 2 (Community)"
  ))

# Edge dataframes
coord_lkp <- nodes_j %>% select(person_id, xj, yj)

make_edges <- function(contacts, pid, coord_lkp) {
  if (length(contacts) == 0) return(NULL)
  data.frame(to = contacts) %>%
    left_join(coord_lkp, by = c("to" = "person_id")) %>%
    rename(x2 = xj, y2 = yj) %>%
    mutate(x1 = sel_j$xj, y1 = sel_j$yj)
}

edges_l1 <- make_edges(l1_contacts, pid, coord_lkp)
edges_l2 <- make_edges(l2_contacts, pid, coord_lkp)
edges_l3 <- make_edges(l3_contacts, pid, coord_lkp)

# ==============================================================================
# [Section 4] Compute distance to each contact (for kernel check)
# ==============================================================================

# Degrees to km (approximate)
lat0_rad <- mean(nodes$y) * pi / 180
deg_to_km_x <- 111.32 * cos(lat0_rad)
deg_to_km_y <- 110.54

compute_dist_km <- function(contacts, sel, nodes) {
  if (length(contacts) == 0) return(numeric(0))
  ctct_nodes <- nodes[nodes$person_id %in% contacts, ]
  dx <- (ctct_nodes$x - sel$x) * deg_to_km_x
  dy <- (ctct_nodes$y - sel$y) * deg_to_km_y
  sqrt(dx^2 + dy^2)
}

dist_l1 <- compute_dist_km(l1_contacts, sel, nodes)
dist_l2 <- compute_dist_km(l2_contacts, sel, nodes)

# ==============================================================================
# [Section 5] Background population density
# ==============================================================================

cell_pop <- nodes %>%
  group_by(x, y) %>%
  summarise(n = n(), .groups = "drop")

# Zoom extent around selected person (±0.15 deg ~ ±15km)
zoom <- 0.15
xlim <- c(sel$x - zoom, sel$x + zoom)
ylim <- c(sel$y - zoom, sel$y + zoom)

# ==============================================================================
# [Section 6] Plot 1 — Full network map (zoomed)
# ==============================================================================

layer_colors <- c(
  "Layer 1 (HH)"       = "#4EC9FF",
  "Layer 2 (Community)" = "#FF8C42",
  "Layer 3 (HCW)"      = "#FFD700",
  "HH + HCW"           = "#FF4EFF"
)

p_network <- ggplot() +
  geom_raster(data = cell_pop %>%
                filter(x >= xlim[1], x <= xlim[2],
                       y >= ylim[1], y <= ylim[2]),
              aes(x = x, y = y, fill = log1p(n)), alpha = 0.5) +
  scale_fill_viridis_c(option = "magma", guide = "none") +

  # Layer 2 edges (draw first — behind)
  { if (!is.null(edges_l2))
    geom_segment(data = edges_l2,
                 aes(x=x1, y=y1, xend=x2, yend=y2),
                 color = "#FF8C42", alpha = 0.5, linewidth = 0.4) } +

  # Layer 3 edges
  { if (!is.null(edges_l3))
    geom_segment(data = edges_l3,
                 aes(x=x1, y=y1, xend=x2, yend=y2),
                 color = "#FFD700", alpha = 0.7, linewidth = 0.6) } +

  # Layer 1 edges
  { if (!is.null(edges_l1))
    geom_segment(data = edges_l1,
                 aes(x=x1, y=y1, xend=x2, yend=y2),
                 color = "#4EC9FF", alpha = 0.8, linewidth = 0.6) } +

  # Contact nodes
  geom_point(data = contact_nodes,
             aes(x = xj, y = yj, color = layer),
             size = 2.0, alpha = 0.9) +
  scale_color_manual(values = layer_colors, name = "Layer") +

  # Selected person (highlighted)
  geom_point(data = sel_j,
             aes(x = xj, y = yj),
             color = "white", size = 5, shape = 21, fill = "red",
             stroke = 1.5) +

  coord_fixed(xlim = xlim, ylim = ylim) +
  labs(title    = sprintf("Person %d — All network connections", pid),
       subtitle = sprintf("Age: %d | HCW: %s | HH: %d L1 | %d L2 | %d L3 contacts",
                          sel$age,
                          ifelse(sel$is_hcw, "YES", "no"),
                          length(l1_contacts),
                          length(l2_contacts),
                          length(l3_contacts)),
       x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        legend.position = "right")

# ==============================================================================
# [Section 7] Plot 2 — Distance distribution of Layer 2 contacts
#   Shows kernel decay effect
# ==============================================================================

# Kernel function for overlay
kernel_path <- "output/kernel/community_distance_kernel.rds"
kernel <- readRDS(kernel_path)
p_hat  <- kernel$p_hat
a1_hat <- kernel$a1_hat
a2_hat <- kernel$a2_hat

d_seq <- seq(0, max(c(dist_l2, 0.1)) * 1.1, length.out = 300)
kernel_df <- data.frame(
  d = d_seq,
  w = p_hat * a1_hat * exp(-a1_hat * d_seq) +
    (1 - p_hat) * a2_hat * exp(-a2_hat * d_seq)
)
# Normalize kernel to overlay on histogram
kernel_df$w_scaled <- kernel_df$w / max(kernel_df$w)

p_dist <- ggplot() +
  # Observed community contact distances
  { if (length(dist_l2) > 0)
    geom_histogram(data = data.frame(d = dist_l2),
                   aes(x = d, y = after_stat(density)),
                   bins = 30, fill = "#FF8C42", alpha = 0.7, color = "white") } +
  # Kernel curve
  geom_line(data = kernel_df,
            aes(x = d, y = w_scaled * max(
              if (length(dist_l2) > 1) density(dist_l2)$y else 1)),
            color = "white", linewidth = 1.2, linetype = "dashed") +
  labs(title    = "Layer 2: Contact distance distribution",
       subtitle = sprintf("Observed contacts: %d | Kernel overlay (dashed)",
                          length(dist_l2)),
       x = "Distance (km)", y = "Density") +
  theme_bw() +
  theme(plot.title    = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey40"))

# ==============================================================================
# [Section 8] Plot 3 — Layer 1 household members age distribution
# ==============================================================================

if (length(l1_contacts) > 0) {
  hh_members <- nodes %>%
    filter(person_id %in% c(pid, l1_contacts)) %>%
    mutate(role = ifelse(person_id == pid, "Selected", "HH member"))

  p_hh_age <- ggplot(hh_members, aes(x = age, fill = role)) +
    geom_histogram(bins = 20, color = "white", alpha = 0.8, position = "stack") +
    scale_fill_manual(values = c("Selected" = "red", "HH member" = "#4EC9FF")) +
    labs(title    = sprintf("Household composition (HH size: %d)",
                            length(l1_contacts) + 1L),
         subtitle = sprintf("Ages: %s",
                            paste(sort(hh_members$age), collapse = ", ")),
         x = "Age", y = "Count", fill = "") +
    theme_bw() +
    theme(plot.title    = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 9, color = "grey40"),
          legend.position = "top")
} else {
  p_hh_age <- ggplot() +
    labs(title = "No household members") +
    theme_void()
}

# ==============================================================================
# [Section 9] Combine and save
# ==============================================================================

p_combined <- p_network / (p_dist + p_hh_age) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(
    title = sprintf("%s — Person %d network", tag, pid),
    theme = theme(plot.title = element_text(size = 14, face = "bold"))
  )

out_path <- file.path(fig_dir, sprintf("%s_person%d_network.png", tag, pid))
ggsave(out_path, plot = p_combined, width = 14, height = 12, dpi = 150)
cat(sprintf("\nSaved: %s\n", out_path))

# ==============================================================================
# [Section 10] Console summary
# ==============================================================================

cat("\n=== Contact summary ===\n")
cat(sprintf("  HH contacts     : %s\n",
            paste(l1_contacts, collapse=", ")))
if (length(dist_l2) > 0) {
  cat(sprintf("  Community contacts: %d\n", length(dist_l2)))
  cat(sprintf("    Distance: mean=%.2f km | max=%.2f km\n",
              mean(dist_l2), max(dist_l2)))
  cat(sprintf("    Within 1km : %d (%.0f%%)\n",
              sum(dist_l2 <= 1), 100*mean(dist_l2 <= 1)))
  cat(sprintf("    1-5km      : %d (%.0f%%)\n",
              sum(dist_l2 > 1 & dist_l2 <= 5),
              100*mean(dist_l2 > 1 & dist_l2 <= 5)))
  cat(sprintf("    5-10km     : %d (%.0f%%)\n",
              sum(dist_l2 > 5), 100*mean(dist_l2 > 5)))
}
if (length(l3_contacts) > 0)
  cat(sprintf("  HCW contacts    : %d\n", length(l3_contacts)))

}
