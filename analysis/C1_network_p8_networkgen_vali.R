# ==============================================================================
# C1_network_p7_networkgen_vali.R
# Purpose:
#   Validate network output from p7 for a selected simulation case.
#
#   Part 1 : Network layer maps (density + L1 + L2 daily/weekly/monthly + L3)
#   Part 1b: L2 edge distance — observed vs kernel expected (per stratum)
#   Part 1c: Prem matrix vs observed contact age structure (per stratum)
#   Part 1d: Physical contact ratio — observed vs blended_ratio target (per stratum)
#   Part 1e: Household contact matrices — close_only_home and phys_only_home
#            (expected values from p6, no is_physical flag in Layer 1)
#   Part 2 : Single person network viewer
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(cowplot)

# ==============================================================================
# [Configuration]
# ==============================================================================

network_dir   <- "output/network"
matrices_path <- "output/MPMmat/DRC_network_input_matrices.rds"
hf_path       <- "data/COD_GRID3_health_facilities_v8.csv"
kernel_path   <- "output/kernel/community_distance_kernel.rds"
fig_dir       <- "figure/C1_p7_validation"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

case_tag      <- "case1_1M"
person_id_sel <- 42L
max_edges_show <- 500L
jitter_deg    <- 0.004

# ==============================================================================
# [Section 1] Load data
# ==============================================================================

message("=== Section 1: Loading data ===")

tag <- case_tag

nodes   <- readRDS(file.path(network_dir, sprintf("%s_nodes.rds",            tag)))
layer1  <- readRDS(file.path(network_dir, sprintf("%s_layer1_household.rds", tag)))
layer2d <- readRDS(file.path(network_dir, sprintf("%s_layer2_daily.rds",     tag)))
layer2w <- readRDS(file.path(network_dir, sprintf("%s_layer2_weekly.rds",    tag)))
layer2m <- readRDS(file.path(network_dir, sprintf("%s_layer2_monthly.rds",   tag)))
layer3h <- readRDS(file.path(network_dir, sprintf("%s_layer3_hcw_edges.rds", tag)))
layer3a <- readRDS(file.path(network_dir, sprintf("%s_layer3_admission.rds", tag)))

mats <- readRDS(matrices_path)

blended_ratio_comm <- (mats$blended_ratio$work +
                         mats$blended_ratio$school +
                         mats$blended_ratio$other) / 3
blended_ratio_home <- mats$blended_ratio$home
close_only_home    <- mats$close_only_home
phys_only_home     <- mats$phys_only_home

hf_all  <- read.csv(hf_path)
hf_terr <- hf_all %>%
  filter(esstype %in% c("Hôpital", "Hôpital Général de Référence",
                        "Centre Hopitalier"),
         !is.na(lon), !is.na(lat)) %>%
  filter(OBJECTID %in% unique(nodes$hospital_id))

message(sprintf("  Nodes         : %d", nrow(nodes)))
message(sprintf("  HCWs          : %d", sum(nodes$is_hcw)))
message(sprintf("  Layer 1 HH    : %d edges", nrow(layer1)))
message(sprintf("  Layer 2 daily : %d edges", nrow(layer2d)))
message(sprintf("  Layer 2 weekly: %d edges", nrow(layer2w)))
message(sprintf("  Layer 2 monthly:%d edges", nrow(layer2m)))
message(sprintf("  Layer 3 HCW   : %d edges", nrow(layer3h)))
message(sprintf("  Hospitals     : %d", nrow(hf_terr)))

kernel          <- readRDS(kernel_path)
kernel_integral <- function(a, b)
  kernel$p_hat     * (exp(-kernel$a1_hat * a) - exp(-kernel$a1_hat * b)) +
  (1-kernel$p_hat) * (exp(-kernel$a2_hat * a) - exp(-kernel$a2_hat * b))

age_grp_labels <- c("0-4","5-9","10-14","15-19","20-24","25-29",
                    "30-34","35-39","40-44","45-49","50-54",
                    "55-59","60-64","65-69","70-74","75+")

# ==============================================================================
# [Section 2] Shared objects
# ==============================================================================

set.seed(42)
nodes_j <- nodes %>%
  mutate(xj = x + runif(n(), -jitter_deg, jitter_deg),
         yj = y + runif(n(), -jitter_deg, jitter_deg))
coord_lookup <- nodes_j %>% select(person_id, xj, yj)

cell_pop <- nodes %>%
  group_by(x, y) %>%
  summarise(n = n(), .groups = "drop")

lat0_rad <- mean(nodes$y) * pi / 180
deg_km_x <- 111.32 * cos(lat0_rad)
deg_km_y <- 110.54

age_lookup <- nodes %>% select(person_id, age_group)

# ==============================================================================
# [Part 1] Network layer maps
# ==============================================================================

message("\n=== Part 1: Network layer maps ===")

sample_edges <- function(layer, n = max_edges_show, seed = 42) {
  set.seed(seed)
  layer %>%
    slice_sample(n = min(n, nrow(layer))) %>%
    left_join(coord_lookup, by = c("from" = "person_id")) %>%
    rename(x1 = xj, y1 = yj) %>%
    left_join(coord_lookup, by = c("to" = "person_id")) %>%
    rename(x2 = xj, y2 = yj)
}

base_map <- function() {
  list(
    geom_raster(data = cell_pop, aes(x = x, y = y, fill = log1p(n)), alpha = 0.5),
    scale_fill_viridis_c(option = "magma", guide = "none"),
    coord_fixed(),
    theme_void(),
    theme(plot.title    = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8,  color = "grey40"))
  )
}

p_density <- ggplot() +
  base_map() +
  geom_raster(data = cell_pop, aes(x = x, y = y, fill = log1p(n))) +
  scale_fill_viridis_c(option = "magma", name = "log(pop+1)") +
  geom_point(data = hf_terr, aes(x = lon, y = lat),
             color = "cyan", size = 3, shape = 17) +
  labs(title    = tag,
       subtitle = sprintf("Pop: %d | Hospitals: %d (▲)", nrow(nodes), nrow(hf_terr)),
       x = NULL, y = NULL) +
  theme(legend.position = "right")

l1_edges <- sample_edges(layer1)
p_layer1 <- ggplot() + base_map() +
  geom_segment(data = l1_edges,
               aes(x = x1, y = y1, xend = x2, yend = y2),
               color = "steelblue", alpha = 0.4, linewidth = 0.3) +
  labs(title    = "Layer 1: Household",
       subtitle = sprintf("%d shown / %d total", nrow(l1_edges), nrow(layer1)))

make_l2_map <- function(layer, label, color) {
  edges <- sample_edges(layer)
  ggplot() + base_map() +
    geom_segment(data = edges,
                 aes(x = x1, y = y1, xend = x2, yend = y2),
                 color = color, alpha = 0.35, linewidth = 0.3) +
    labs(title    = sprintf("Layer 2: %s", label),
         subtitle = sprintf("%d shown / %d total", nrow(edges), nrow(layer)))
}

p_l2d <- make_l2_map(layer2d, "Community daily",   "#993C1D")
p_l2w <- make_l2_map(layer2w, "Community weekly",  "#0F6E56")
p_l2m <- make_l2_map(layer2m, "Community monthly", "#185FA5")

hcw_nodes_j <- nodes_j %>% filter(is_hcw)
l3_edges    <- sample_edges(layer3h)
p_layer3 <- ggplot() + base_map() +
  geom_segment(data = l3_edges,
               aes(x = x1, y = y1, xend = x2, yend = y2),
               color = "yellow", alpha = 0.5, linewidth = 0.4) +
  geom_point(data = hcw_nodes_j, aes(x = xj, y = yj),
             color = "yellow", size = 0.8, alpha = 0.7) +
  geom_point(data = hf_terr, aes(x = lon, y = lat),
             color = "cyan", size = 3, shape = 17) +
  labs(title    = "Layer 3: HCW",
       subtitle = sprintf("HCWs: %d | edges: %d | hospitals: %d",
                          nrow(hcw_nodes_j), nrow(layer3h), nrow(hf_terr)))

p_layers <- (p_density + p_layer1 + p_l2d) /
  (p_l2w    + p_l2m    + p_layer3) +
  plot_annotation(
    title    = sprintf("Network layers — %s", tag),
    subtitle = sprintf(
      "Pop: %d | HH: %d | L1: %d | L2 daily: %d | weekly: %d | monthly: %d | L3: %d",
      nrow(nodes), n_distinct(nodes$hh_id),
      nrow(layer1), nrow(layer2d), nrow(layer2w), nrow(layer2m), nrow(layer3h)),
    theme = theme(plot.title    = element_text(size = 13, face = "bold"),
                  plot.subtitle = element_text(size = 9,  color = "grey40"))
  )

out_p1 <- file.path(fig_dir, sprintf("%s_network_layers.png", tag))
ggsave(out_p1, plot = p_layers, width = 18, height = 12, dpi = 150)
message(sprintf("  Saved: %s", out_p1))

# ==============================================================================
# [Part 1b] L2 edge distance — observed vs kernel expected, per stratum
# ==============================================================================

message("\n=== Part 1b: Edge distance validation ===")

bucket_breaks <- c(0, 1.5, 10.5, 100.5, Inf)
bucket_labels <- c("0-1.5km", "1.5-10.5km", "10.5-100.5km", "100.5km+")

bw <- c(kernel_integral(0, 1.5), kernel_integral(1.5, 10.5),
        kernel_integral(10.5, 100.5), kernel_integral(100.5, Inf))
exp_df <- data.frame(
  bucket  = factor(bucket_labels, levels = bucket_labels),
  pct     = 100 * bw / sum(bw),
  stratum = "Expected (kernel)"
)

coords <- nodes %>% select(person_id, x, y)

compute_bucket_obs <- function(layer, stratum_label) {
  layer %>%
    left_join(coords, by = c("from" = "person_id")) %>% rename(x1=x, y1=y) %>%
    left_join(coords, by = c("to"   = "person_id")) %>% rename(x2=x, y2=y) %>%
    mutate(dist_km = sqrt(((x2-x1)*deg_km_x)^2 + ((y2-y1)*deg_km_y)^2),
           bucket  = cut(dist_km, breaks = bucket_breaks,
                         labels = bucket_labels, right = FALSE)) %>%
    group_by(bucket) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(pct = 100 * n / sum(n), stratum = stratum_label)
}

bucket_comp <- bind_rows(
  compute_bucket_obs(layer2d, "Daily")   %>% select(bucket, pct, stratum),
  compute_bucket_obs(layer2w, "Weekly")  %>% select(bucket, pct, stratum),
  compute_bucket_obs(layer2m, "Monthly") %>% select(bucket, pct, stratum),
  exp_df
) %>%
  mutate(bucket  = factor(bucket,  levels = bucket_labels),
         stratum = factor(stratum, levels = c("Daily","Weekly","Monthly",
                                              "Expected (kernel)")))

p_bucket <- ggplot(bucket_comp, aes(x = bucket, y = pct, fill = stratum)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 3) +
  scale_fill_manual(
    values = c("Daily"="#993C1D","Weekly"="#0F6E56",
               "Monthly"="#185FA5","Expected (kernel)"="grey50"),
    name = NULL) +
  labs(title    = sprintf("%s — L2 edge distance: observed vs kernel", tag),
       subtitle = "Observed per frequency stratum vs kernel expected proportions",
       x = "Distance bucket", y = "% of edges") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        legend.position = "top")

out_bucket <- file.path(fig_dir, sprintf("%s_l2_bucket_dist.png", tag))
ggsave(out_bucket, plot = p_bucket, width = 9, height = 5, dpi = 150)
message(sprintf("  Saved: %s", out_bucket))

# ==============================================================================
# [Part 1c] Prem matrix vs observed contact age structure, per stratum
# ==============================================================================

message("\n=== Part 1c: Prem vs observed age structure ===")

build_obs_mat <- function(layer) {
  obs <- layer %>%
    left_join(age_lookup, by = c("from" = "person_id")) %>% rename(ag_from = age_group) %>%
    left_join(age_lookup, by = c("to"   = "person_id")) %>% rename(ag_to   = age_group) %>%
    filter(!is.na(ag_from), !is.na(ag_to))
  mat <- matrix(0L, 16, 16)
  for (k in seq_len(nrow(obs))) {
    mat[obs$ag_from[k], obs$ag_to[k]] <- mat[obs$ag_from[k], obs$ag_to[k]] + 1L
    mat[obs$ag_to[k], obs$ag_from[k]] <- mat[obs$ag_to[k], obs$ag_from[k]] + 1L
  }
  mat
}

normalize_rows <- function(mat) {
  rs <- rowSums(mat); rs[rs == 0] <- 1; mat / rs
}

mat_to_long_hm <- function(mat, label) {
  df <- as.data.frame(mat)
  colnames(df) <- age_grp_labels
  df$from_age  <- age_grp_labels
  df$type      <- label
  pivot_longer(df, cols = all_of(age_grp_labels),
               names_to = "to_age", values_to = "value")
}

prem_list  <- list(Daily   = mats$close_3wk_daily$community,
                   Weekly  = mats$close_3wk_weekly$community,
                   Monthly = mats$close_3wk_monthly$community)
layer2_list <- list(Daily = layer2d, Weekly = layer2w, Monthly = layer2m)

prem_plots <- lapply(names(prem_list), function(st) {
  prem_norm <- normalize_rows((prem_list[[st]] + t(prem_list[[st]])) / 2)
  obs_norm  <- normalize_rows({m <- build_obs_mat(layer2_list[[st]]); (m + t(m)) / 2})
  r <- cor(as.vector(prem_norm), as.vector(obs_norm))
  message(sprintf("  [%s] Pearson r = %.4f", st, r))

  bind_rows(mat_to_long_hm(prem_norm, "Prem (expected)"),
            mat_to_long_hm(obs_norm,  "Observed")) %>%
    mutate(from_age = factor(from_age, levels = age_grp_labels),
           to_age   = factor(to_age,   levels = age_grp_labels)) %>%
    ggplot(aes(x = to_age, y = from_age, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", name = "Proportion") +
    facet_wrap(~ type) +
    labs(title    = sprintf("%s contacts", st),
         subtitle = sprintf("Pearson r = %.4f", r),
         x = "Contact age", y = "Participant age") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          axis.text.y = element_text(size = 6),
          plot.title  = element_text(size = 10, face = "bold"),
          strip.text  = element_text(face = "bold"))
})

out_prem <- file.path(fig_dir, sprintf("%s_prem_comparison.png", tag))
ggsave(out_prem,
       wrap_plots(prem_plots, ncol = 1) +
         plot_annotation(title = sprintf("%s — Prem vs Observed", tag),
                         theme = theme(plot.title = element_text(size = 13, face = "bold"))),
       width = 12, height = 15, dpi = 150)
message(sprintf("  Saved: %s", out_prem))

# ==============================================================================
# [Part 1d] Physical contact ratio — observed vs blended_ratio target (L2)
# ==============================================================================

message("\n=== Part 1d: Physical ratio validation (community) ===")

compute_phys_ratio_obs <- function(layer, stratum_label) {
  layer %>%
    left_join(age_lookup, by = c("from" = "person_id")) %>% rename(ag_from = age_group) %>%
    left_join(age_lookup, by = c("to"   = "person_id")) %>% rename(ag_to   = age_group) %>%
    filter(!is.na(ag_from), !is.na(ag_to)) %>%
    group_by(ag_from, ag_to) %>%
    summarise(obs_ratio = mean(is_physical), n = n(), .groups = "drop") %>%
    mutate(exp_ratio = mapply(function(i,j) blended_ratio_comm[i,j], ag_from, ag_to),
           stratum   = stratum_label)
}

phys_all <- bind_rows(
  compute_phys_ratio_obs(layer2d, "Daily"),
  compute_phys_ratio_obs(layer2w, "Weekly"),
  compute_phys_ratio_obs(layer2m, "Monthly")
) %>% mutate(stratum = factor(stratum, levels = c("Daily","Weekly","Monthly")))

for (st in c("Daily","Weekly","Monthly")) {
  sub <- phys_all %>% filter(stratum == st)
  message(sprintf("  [%s] mean obs=%.4f | mean exp=%.4f",
                  st, weighted.mean(sub$obs_ratio, sub$n),
                  weighted.mean(sub$exp_ratio, sub$n)))
}

p_phys <- ggplot(phys_all, aes(x = exp_ratio, y = obs_ratio, size = n)) +
  geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
  geom_point(alpha = 0.5, color = "#185FA5") +
  facet_wrap(~ stratum) +
  scale_size_continuous(range = c(0.5, 4), name = "Edge count") +
  labs(title    = sprintf("%s — Community physical ratio: observed vs expected", tag),
       subtitle = "Each point = one (age_i, age_j) cell | dashed = perfect agreement",
       x = "Expected (blended_ratio_comm)", y = "Observed (is_physical mean)") +
  theme_bw() +
  theme(plot.title = element_text(size = 12, face = "bold"),
        strip.text = element_text(face = "bold"))

out_phys <- file.path(fig_dir, sprintf("%s_physical_ratio.png", tag))
ggsave(out_phys, plot = p_phys, width = 12, height = 5, dpi = 150)
message(sprintf("  Saved: %s", out_phys))

# ==============================================================================
# [Part 1e] Household contact matrices — close_only_home and phys_only_home
# Layer 1 has no is_physical flag; physical structure handled in simulation
# via age-group-specific matrices from p6.
# Here we visualize the expected matrices and check row sums.
# ==============================================================================

message("\n=== Part 1e: Household contact matrices (from p6) ===")

message(sprintf("  close_only_home rowSum: %.3f - %.3f",
                min(rowSums(close_only_home)), max(rowSums(close_only_home))))
message(sprintf("  phys_only_home  rowSum: %.3f - %.3f",
                min(rowSums(phys_only_home)),  max(rowSums(phys_only_home))))
message(sprintf("  sum check vs prem_mats$home: %.3f - %.3f",
                min(rowSums(close_only_home + phys_only_home)),
                max(rowSums(close_only_home + phys_only_home))))
message(sprintf("  phys proportion (mean): %.3f",
                mean(rowSums(phys_only_home) /
                       (rowSums(close_only_home) + rowSums(phys_only_home)))))

make_hh_heatmap <- function(mat, title, high_col) {
  df <- as.data.frame(mat)
  colnames(df) <- age_grp_labels
  df$from_age  <- age_grp_labels
  df_long <- pivot_longer(df, cols = all_of(age_grp_labels),
                          names_to = "to_age", values_to = "value") %>%
    mutate(from_age = factor(from_age, levels = age_grp_labels),
           to_age   = factor(to_age,   levels = age_grp_labels))

  ggplot(df_long, aes(x = to_age, y = from_age, fill = value)) +
    geom_tile(color = NA) +
    scale_fill_gradient(low = "#FAEEDA", high = high_col,
                        name = "contacts/day") +
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    labs(title = title, x = "Contact age", y = "Participant age") +
    theme_minimal(base_size = 8) +
    theme(plot.title  = element_text(size = 9, face = "bold", hjust = 0.5),
          axis.text   = element_text(size = 5),
          axis.title  = element_text(size = 7),
          panel.grid  = element_blank(),
          legend.title = element_text(size = 6),
          legend.text  = element_text(size = 5),
          legend.key.height = unit(0.4, "cm"))
}

p_close_hm <- make_hh_heatmap(close_only_home,
                              "Close-only (Prem × (1 - phys_ratio))",
                              "#BA7517")
p_phys_hm  <- make_hh_heatmap(phys_only_home,
                              "Physical-only (Prem × phys_ratio)",
                              "#0C447C")
p_ratio_hm <- make_hh_heatmap(blended_ratio_home,
                              "Physical ratio (blended_ratio$home)",
                              "#0C447C")

# Age-specific physical proportion bar chart
phys_prop_df <- data.frame(
  age_group  = factor(age_grp_labels, levels = age_grp_labels),
  close_only = rowSums(close_only_home),
  phys_only  = rowSums(phys_only_home)
) %>%
  mutate(phys_pct = 100 * phys_only / (close_only + phys_only))

p_phys_pct <- ggplot(phys_prop_df, aes(x = age_group, y = phys_pct)) +
  geom_col(fill = "#0C447C", alpha = 0.8, width = 0.7) +
  geom_hline(yintercept = mean(phys_prop_df$phys_pct),
             color = "tomato", linetype = "dashed", linewidth = 0.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title    = sprintf("Physical proportion by age group (mean = %.1f%%)",
                          mean(phys_prop_df$phys_pct)),
       x = "Age group", y = "% physical of household contacts") +
  theme_bw() +
  theme(plot.title  = element_text(size = 9, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

title_1e <- ggdraw() +
  draw_label(sprintf("%s — Household contact matrices (p6 output)", tag),
             fontface = "bold", size = 12, x = 0.5, hjust = 0.5)

p_hh_mats <- plot_grid(p_close_hm, p_ratio_hm, p_phys_hm,
                       nrow = 1, ncol = 3)

p_1e <- plot_grid(title_1e,
                  p_hh_mats,
                  p_phys_pct,
                  ncol = 1,
                  rel_heights = c(0.06, 1, 0.5))

out_1e <- file.path(fig_dir, sprintf("%s_household_matrices.png", tag))
ggsave(out_1e, plot = p_1e, width = 13, height = 10, dpi = 150)
message(sprintf("  Saved: %s", out_1e))

# ==============================================================================
# [Part 2] Single person network viewer
# ==============================================================================

message("\n=== Part 2: Single person viewer ===")

pid   <- max(1L, min(as.integer(person_id_sel), nrow(nodes)))
sel   <- nodes[nodes$person_id == pid, ]
sel_j <- nodes_j[nodes_j$person_id == pid, ]

message(sprintf("  Person %d | Age: %d | HCW: %s | Cell: %d",
                pid, sel$age, ifelse(sel$is_hcw, "YES", "no"), sel$cell_id))

l1_contacts  <- c(layer1$to[layer1$from==pid],   layer1$from[layer1$to==pid])
l2d_contacts <- c(layer2d$to[layer2d$from==pid],  layer2d$from[layer2d$to==pid])
l2w_contacts <- c(layer2w$to[layer2w$from==pid],  layer2w$from[layer2w$to==pid])
l2m_contacts <- c(layer2m$to[layer2m$from==pid],  layer2m$from[layer2m$to==pid])
l3_contacts  <- c(layer3h$to[layer3h$from==pid],  layer3h$from[layer3h$to==pid])

l2_contacts_all <- unique(c(l2d_contacts, l2w_contacts, l2m_contacts))

message(sprintf("  L1: %d | L2 daily: %d | weekly: %d | monthly: %d | L3: %d",
                length(l1_contacts), length(l2d_contacts),
                length(l2w_contacts), length(l2m_contacts), length(l3_contacts)))

all_contacts  <- unique(c(l1_contacts, l2_contacts_all, l3_contacts))
contact_nodes <- nodes_j %>%
  filter(person_id %in% all_contacts) %>%
  mutate(layer = case_when(
    person_id %in% l1_contacts  ~ "Layer 1 (HH)",
    person_id %in% l3_contacts  ~ "Layer 3 (HCW)",
    person_id %in% l2d_contacts ~ "L2 Daily",
    person_id %in% l2w_contacts ~ "L2 Weekly",
    TRUE                        ~ "L2 Monthly"
  ))

make_edges <- function(contacts, coord_lkp, sel_j) {
  if (length(contacts) == 0) return(NULL)
  data.frame(to = contacts) %>%
    left_join(coord_lkp, by = c("to" = "person_id")) %>%
    rename(x2 = xj, y2 = yj) %>%
    mutate(x1 = sel_j$xj, y1 = sel_j$yj)
}

edges_l1  <- make_edges(l1_contacts,  coord_lookup, sel_j)
edges_l2d <- make_edges(l2d_contacts, coord_lookup, sel_j)
edges_l2w <- make_edges(l2w_contacts, coord_lookup, sel_j)
edges_l2m <- make_edges(l2m_contacts, coord_lookup, sel_j)
edges_l3  <- make_edges(l3_contacts,  coord_lookup, sel_j)

all_x <- c(sel$x, if (length(all_contacts) > 0) nodes$x[nodes$person_id %in% all_contacts])
all_y <- c(sel$y, if (length(all_contacts) > 0) nodes$y[nodes$person_id %in% all_contacts])
pad   <- 0.005
xlim  <- c(min(all_x) - pad, max(all_x) + pad)
ylim  <- c(min(all_y) - pad, max(all_y) + pad)

layer_colors <- c("Layer 1 (HH)"="#4EC9FF","L2 Daily"="#993C1D",
                  "L2 Weekly"="#0F6E56","L2 Monthly"="#185FA5",
                  "Layer 3 (HCW)"="#FFD700")

p_network <- ggplot() +
  geom_raster(data = cell_pop %>%
                filter(x >= xlim[1], x <= xlim[2], y >= ylim[1], y <= ylim[2]),
              aes(x = x, y = y, fill = log1p(n)), alpha = 0.5) +
  scale_fill_viridis_c(option = "magma", guide = "none") +
  { if (!is.null(edges_l2m)) geom_segment(data=edges_l2m, aes(x=x1,y=y1,xend=x2,yend=y2), color="#185FA5", alpha=0.4, linewidth=0.4) } +
  { if (!is.null(edges_l2w)) geom_segment(data=edges_l2w, aes(x=x1,y=y1,xend=x2,yend=y2), color="#0F6E56", alpha=0.4, linewidth=0.4) } +
  { if (!is.null(edges_l2d)) geom_segment(data=edges_l2d, aes(x=x1,y=y1,xend=x2,yend=y2), color="#993C1D", alpha=0.5, linewidth=0.5) } +
  { if (!is.null(edges_l3))  geom_segment(data=edges_l3,  aes(x=x1,y=y1,xend=x2,yend=y2), color="#FFD700", alpha=0.7, linewidth=0.6) } +
  { if (!is.null(edges_l1))  geom_segment(data=edges_l1,  aes(x=x1,y=y1,xend=x2,yend=y2), color="#4EC9FF", alpha=0.8, linewidth=0.6) } +
  geom_point(data = contact_nodes, aes(x=xj, y=yj, color=layer), size=2.0, alpha=0.9) +
  scale_color_manual(values = layer_colors, name = "Layer") +
  geom_point(data = sel_j, aes(x=xj, y=yj),
             color="white", size=5, shape=21, fill="red", stroke=1.5) +
  coord_fixed(xlim = xlim, ylim = ylim) +
  labs(title    = sprintf("Person %d — network connections", pid),
       subtitle = sprintf("Age: %d | HCW: %s | L1: %d | L2 d/w/m: %d/%d/%d | L3: %d",
                          sel$age, ifelse(sel$is_hcw,"YES","no"),
                          length(l1_contacts), length(l2d_contacts),
                          length(l2w_contacts), length(l2m_contacts), length(l3_contacts)),
       x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(plot.title = element_text(size=12, face="bold"),
        legend.position = "right")

dist_km_fn <- function(contacts) {
  if (length(contacts) == 0) return(numeric(0))
  cn <- nodes[nodes$person_id %in% contacts, ]
  sqrt(((cn$x - sel$x) * deg_km_x)^2 + ((cn$y - sel$y) * deg_km_y)^2)
}

dist_df <- bind_rows(
  data.frame(d = dist_km_fn(l2d_contacts), stratum = "Daily"),
  data.frame(d = dist_km_fn(l2w_contacts), stratum = "Weekly"),
  data.frame(d = dist_km_fn(l2m_contacts), stratum = "Monthly")
) %>% filter(is.finite(d))

d_all   <- dist_df$d
d_seq   <- seq(0, max(c(d_all, 0.1)) * 1.1, length.out = 300)
kern_df <- data.frame(
  d = d_seq,
  w = kernel$p_hat * kernel$a1_hat * exp(-kernel$a1_hat * d_seq) +
    (1-kernel$p_hat) * kernel$a2_hat * exp(-kernel$a2_hat * d_seq)
)
kern_df$w_scaled <- ifelse(length(d_all) > 1,
                           kern_df$w / max(kern_df$w) * max(density(d_all)$y), kern_df$w)

p_dist <- ggplot() +
  { if (nrow(dist_df) > 0)
    geom_histogram(data=dist_df, aes(x=d, y=after_stat(density), fill=stratum),
                   bins=30, alpha=0.6, position="identity") } +
  scale_fill_manual(values=c("Daily"="#993C1D","Weekly"="#0F6E56","Monthly"="#185FA5"),
                    name="Stratum") +
  geom_line(data=kern_df, aes(x=d, y=w_scaled),
            color="white", linewidth=1.2, linetype="dashed") +
  labs(title    = "L2: Contact distance distribution",
       subtitle = sprintf("%d contacts | kernel (dashed)", nrow(dist_df)),
       x = "Distance (km)", y = "Density") +
  theme_bw() +
  theme(plot.title = element_text(size=11, face="bold"))

if (length(l1_contacts) > 0) {
  hh_members <- nodes %>%
    filter(person_id %in% c(pid, l1_contacts)) %>%
    mutate(role = ifelse(person_id == pid, "Selected", "HH member"))
  p_hh_age <- ggplot(hh_members, aes(x = age, fill = role)) +
    geom_histogram(bins=20, color="white", alpha=0.8, position="stack") +
    scale_fill_manual(values=c("Selected"="red","HH member"="#4EC9FF")) +
    labs(title    = sprintf("Household (size: %d)", length(l1_contacts)+1L),
         subtitle = sprintf("Ages: %s", paste(sort(hh_members$age), collapse=", ")),
         x="Age", y="Count", fill="") +
    theme_bw() +
    theme(plot.title = element_text(size=11, face="bold"), legend.position="top")
} else {
  p_hh_age <- ggplot() + labs(title="No household members") + theme_void()
}

p_person <- p_network / (p_dist + p_hh_age) +
  plot_layout(heights = c(2, 1)) +
  plot_annotation(title = sprintf("%s — Person %d", tag, pid),
                  theme = theme(plot.title = element_text(size=14, face="bold")))

out_p2 <- file.path(fig_dir, sprintf("%s_person%d.png", tag, pid))
ggsave(out_p2, plot = p_person, width=14, height=12, dpi=150)
message(sprintf("  Saved: %s", out_p2))

message("\n=== Contact summary ===")
message(sprintf("  HH contacts : %s", paste(l1_contacts, collapse=", ")))
if (nrow(dist_df) > 0) {
  message(sprintf("  L2 daily    : %d | weekly: %d | monthly: %d",
                  length(l2d_contacts), length(l2w_contacts), length(l2m_contacts)))
  message(sprintf("  L2 dist     : mean=%.2f km | max=%.2f km", mean(d_all), max(d_all)))
}
if (length(l3_contacts) > 0)
  message(sprintf("  HCW contacts: %d", length(l3_contacts)))

message("\n=== Done ===")
message(sprintf("  Figures saved to: %s/", fig_dir))
