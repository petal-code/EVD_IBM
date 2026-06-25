# ==============================================================================
# C1_network_p6_validation.R
# Purpose:
#   Visualize and validate network for simulation cases:
#     Part 1: Network layer maps (density + L1 + L2 + L3)
#     Part 2: Single person network + kernel distance check
# ==============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# ==============================================================================
# [Configuration] — Change ONLY these
# ==============================================================================

network_dir <- "output/network"
hf_path     <- "data/COD_GRID3_health_facilities_v8.csv"
kernel_path <- "output/kernel/community_distance_kernel.rds"
fig_dir     <- "figure/C1_p6_validation"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Select case: "case1_1M", "case2_Ituri", "case3_Kivu"
case_tag <- "case1_1M"

# Person ID to inspect (Part 2)
person_id_sel <- 42L

# Max edges to show per layer
max_edges_show <- 500L

# Jitter within cell (1km cell)
jitter_deg <- 0.004

# ==============================================================================
# [Section 1] Load data
# ==============================================================================

net_files  <- list.files(network_dir, pattern="_nodes\\.rds$", full.names=TRUE)
avail_tags <- sub("_nodes\\.rds$", "", basename(net_files))
cat("Available cases:\n")
for (t in avail_tags) cat(sprintf("  %s\n", t))

if (!case_tag %in% avail_tags)
  stop(sprintf("case_tag '%s' not found. Choose from: %s",
               case_tag, paste(avail_tags, collapse=", ")))

tag <- case_tag
cat(sprintf("\nLoading: %s\n", tag))

nodes   <- readRDS(file.path(network_dir, sprintf("%s_nodes.rds",            tag)))
layer1  <- readRDS(file.path(network_dir, sprintf("%s_layer1_household.rds", tag)))
layer2  <- readRDS(file.path(network_dir, sprintf("%s_layer2_community.rds", tag)))
layer3h <- readRDS(file.path(network_dir, sprintf("%s_layer3_hcw_edges.rds", tag)))
layer3a <- readRDS(file.path(network_dir, sprintf("%s_layer3_admission.rds", tag)))

cat(sprintf("  Nodes       : %d\n", nrow(nodes)))
cat(sprintf("  HCWs        : %d\n", sum(nodes$is_hcw)))
cat(sprintf("  Layer 1 HH  : %d edges\n", nrow(layer1)))
cat(sprintf("  Layer 2 comm: %d edges\n", nrow(layer2)))
cat(sprintf("  Layer 3 HCW : %d edges\n", nrow(layer3h)))

hf_all  <- read.csv(hf_path)
hf_terr <- hf_all %>%
  filter(esstype %in% c("Hôpital", "Hôpital Général de Référence",
                        "Centre Hopitalier"),
         !is.na(lon), !is.na(lat)) %>%
  filter(OBJECTID %in% unique(nodes$hospital_id))
cat(sprintf("  Hospitals   : %d\n", nrow(hf_terr)))

# ==============================================================================
# [Section 2] Jitter + shared objects
# ==============================================================================

set.seed(42)
nodes_j <- nodes %>%
  mutate(xj = x + runif(n(), -jitter_deg, jitter_deg),
         yj = y + runif(n(), -jitter_deg, jitter_deg))
coord_lookup <- nodes_j %>% select(person_id, xj, yj)

cell_pop <- nodes %>%
  group_by(x, y) %>%
  summarise(n = n(), .groups = "drop")

# ==============================================================================
# [Part 1] Network layer maps
# ==============================================================================

p_density <- ggplot() +
  geom_raster(data=cell_pop, aes(x=x,y=y,fill=log1p(n))) +
  scale_fill_viridis_c(option="magma", name="log(pop+1)") +
  geom_point(data=hf_terr, aes(x=lon,y=lat),
             color="cyan", size=3, shape=17) +
  geom_point(data=hf_terr, aes(x=lon,y=lat),
             color="white", size=1.2, shape=17) +
  coord_fixed() +
  labs(title=tag,
       subtitle=sprintf("Population: %d | Hospitals: %d (▲)",
                        nrow(nodes), nrow(hf_terr)),
       x=NULL, y=NULL) +
  theme_void() +
  theme(plot.title=element_text(size=11, face="bold"),
        plot.subtitle=element_text(size=9, color="grey40"),
        legend.position="right")

set.seed(42)
l1_sample <- layer1 %>% slice_sample(n=min(max_edges_show, nrow(layer1)))
l1_edges  <- l1_sample %>%
  left_join(coord_lookup, by=c("from"="person_id")) %>% rename(x1=xj,y1=yj) %>%
  left_join(coord_lookup, by=c("to"="person_id"))   %>% rename(x2=xj,y2=yj)

p_layer1 <- ggplot() +
  geom_raster(data=cell_pop, aes(x=x,y=y,fill=log1p(n)), alpha=0.4) +
  scale_fill_viridis_c(option="magma", guide="none") +
  geom_segment(data=l1_edges, aes(x=x1,y=y1,xend=x2,yend=y2),
               color="steelblue", alpha=0.4, linewidth=0.3) +
  geom_point(data=nodes_j, aes(x=xj,y=yj), color="white", size=0.3, alpha=0.5) +
  coord_fixed() +
  labs(title="Layer 1: Household edges",
       subtitle=sprintf("%d shown (of %d total)", nrow(l1_sample), nrow(layer1)),
       x=NULL, y=NULL) +
  theme_void() +
  theme(plot.title=element_text(size=11, face="bold"),
        plot.subtitle=element_text(size=9, color="grey40"))

set.seed(42)
l2_sample <- layer2 %>% slice_sample(n=min(max_edges_show, nrow(layer2)))
l2_edges  <- l2_sample %>%
  left_join(coord_lookup, by=c("from"="person_id")) %>% rename(x1=xj,y1=yj) %>%
  left_join(coord_lookup, by=c("to"="person_id"))   %>% rename(x2=xj,y2=yj)

p_layer2 <- ggplot() +
  geom_raster(data=cell_pop, aes(x=x,y=y,fill=log1p(n)), alpha=0.4) +
  scale_fill_viridis_c(option="magma", guide="none") +
  geom_segment(data=l2_edges, aes(x=x1,y=y1,xend=x2,yend=y2),
               color="coral", alpha=0.3, linewidth=0.3) +
  geom_point(data=nodes_j, aes(x=xj,y=yj), color="white", size=0.3, alpha=0.5) +
  coord_fixed() +
  labs(title="Layer 2: Community edges",
       subtitle=sprintf("%d shown (of %d total)", nrow(l2_sample), nrow(layer2)),
       x=NULL, y=NULL) +
  theme_void() +
  theme(plot.title=element_text(size=11, face="bold"),
        plot.subtitle=element_text(size=9, color="grey40"))

hcw_nodes_j <- nodes_j %>% filter(is_hcw)
set.seed(42)
l3_sample <- layer3h %>% slice_sample(n=min(max_edges_show, nrow(layer3h)))
l3_edges  <- l3_sample %>%
  left_join(coord_lookup, by=c("from"="person_id")) %>% rename(x1=xj,y1=yj) %>%
  left_join(coord_lookup, by=c("to"="person_id"))   %>% rename(x2=xj,y2=yj)

p_layer3 <- ggplot() +
  geom_raster(data=cell_pop, aes(x=x,y=y,fill=log1p(n)), alpha=0.4) +
  scale_fill_viridis_c(option="magma", guide="none") +
  geom_segment(data=l3_edges, aes(x=x1,y=y1,xend=x2,yend=y2),
               color="yellow", alpha=0.5, linewidth=0.4) +
  geom_point(data=hcw_nodes_j, aes(x=xj,y=yj),
             color="yellow", size=1.0, alpha=0.8) +
  geom_point(data=hf_terr, aes(x=lon,y=lat), color="cyan", size=3, shape=17) +
  coord_fixed() +
  labs(title="Layer 3: HCW network",
       subtitle=sprintf("HCWs: %d | edges: %d | Hospitals: %d",
                        nrow(hcw_nodes_j), nrow(layer3h), nrow(hf_terr)),
       x=NULL, y=NULL) +
  theme_void() +
  theme(plot.title=element_text(size=11, face="bold"),
        plot.subtitle=element_text(size=9, color="grey40"))

p_layers <- (p_density + p_layer1) / (p_layer2 + p_layer3) +
  plot_annotation(
    title    = sprintf("Network layers — %s", tag),
    subtitle = sprintf("Pop: %d | HH: %d | L1: %d | L2: %d | L3: %d edges",
                       nrow(nodes), n_distinct(nodes$hh_id),
                       nrow(layer1), nrow(layer2), nrow(layer3h)),
    theme = theme(plot.title=element_text(size=13, face="bold"),
                  plot.subtitle=element_text(size=10, color="grey40"))
  )

out_p1 <- file.path(fig_dir, sprintf("%s_network_layers.png", tag))
ggsave(out_p1, plot=p_layers, width=16, height=12, dpi=150)
cat(sprintf("\nSaved: %s\n", out_p1))

# ==============================================================================
# [Part 1b] Layer 2 edge distance — bucket observed vs expected
# ==============================================================================

lat0_rad <- mean(nodes$y) * pi / 180
deg_km_x <- 111.32 * cos(lat0_rad)
deg_km_y <- 110.54

coords <- nodes %>% select(person_id, x, y)
edge_dist_all <- layer2 %>%
  left_join(coords, by=c("from"="person_id")) %>% rename(x1=x, y1=y) %>%
  left_join(coords, by=c("to"="person_id"))   %>% rename(x2=x, y2=y) %>%
  mutate(dist_km = sqrt(((x2-x1)*deg_km_x)^2 + ((y2-y1)*deg_km_y)^2))

# Assign to buckets
bucket_breaks <- c(0, 1.5, 10.5, 100.5, Inf)
bucket_labels <- c("0-1.5km", "1.5-10.5km", "10.5-100.5km", "100.5km+")
edge_dist_all <- edge_dist_all %>%
  mutate(bucket = cut(dist_km, breaks=bucket_breaks,
                      labels=bucket_labels, right=FALSE))

obs_df <- edge_dist_all %>%
  group_by(bucket) %>%
  summarise(n=n(), .groups="drop") %>%
  mutate(pct = 100 * n / sum(n),
         type = "Observed")

# Expected from kernel integral weights (buckets 0+1 merged as 0-1.5km)
kernel  <- readRDS(kernel_path)
kernel_integral <- function(a, b)
  kernel$p_hat * (exp(-kernel$a1_hat*a) - exp(-kernel$a1_hat*b)) +
  (1-kernel$p_hat) * (exp(-kernel$a2_hat*a) - exp(-kernel$a2_hat*b))

bw <- c(
  kernel_integral(0,     1.5),    # Bucket 0+1: 0-1.5km
  kernel_integral(1.5,  10.5),    # Bucket 2
  kernel_integral(10.5, 100.5),   # Bucket 3
  kernel_integral(100.5, Inf)     # Bucket 4
)

exp_df <- data.frame(
  bucket = factor(bucket_labels, levels=bucket_labels),
  pct    = 100 * bw / sum(bw),
  type   = "Expected (kernel)"
)

comp_df <- bind_rows(obs_df %>% select(bucket, pct, type), exp_df)

p_bucket <- ggplot(comp_df, aes(x=bucket, y=pct, fill=type)) +
  geom_col(position="dodge", alpha=0.85, width=0.6) +
  scale_fill_manual(values=c("Observed"="steelblue",
                             "Expected (kernel)"="tomato"),
                    name=NULL) +
  geom_text(aes(label=sprintf("%.1f%%", pct)),
            position=position_dodge(width=0.6),
            vjust=-0.5, size=3.5) +
  labs(title    = sprintf("%s — L2 edges: observed vs expected by bucket", tag),
       subtitle = sprintf("Total edges: %d", nrow(layer2)),
       x="Distance bucket", y="% of edges") +
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"),
        plot.subtitle=element_text(size=9, color="grey40"),
        legend.position="top")

out_bucket <- file.path(fig_dir, sprintf("%s_l2_bucket_dist.png", tag))
ggsave(out_bucket, plot=p_bucket, width=8, height=5, dpi=150)
cat(sprintf("Saved: %s\n", out_bucket))

# ==============================================================================
# [Part 2] Single person network
# ==============================================================================

pid   <- max(1L, min(as.integer(person_id_sel), nrow(nodes)))
sel   <- nodes[nodes$person_id == pid, ]
sel_j <- nodes_j[nodes_j$person_id == pid, ]
cat(sprintf("\nPerson %d | Age: %d | HCW: %s | Cell: %d\n",
            pid, sel$age, ifelse(sel$is_hcw,"YES","no"), sel$cell_id))

l1_contacts <- c(layer1$to[layer1$from==pid],  layer1$from[layer1$to==pid])
l2_contacts <- c(layer2$to[layer2$from==pid],  layer2$from[layer2$to==pid])
l3_contacts <- c(layer3h$to[layer3h$from==pid], layer3h$from[layer3h$to==pid])

cat(sprintf("  Layer 1: %d | Layer 2: %d | Layer 3: %d contacts\n",
            length(l1_contacts), length(l2_contacts), length(l3_contacts)))

all_contacts  <- unique(c(l1_contacts, l2_contacts, l3_contacts))
contact_nodes <- nodes_j %>%
  filter(person_id %in% all_contacts) %>%
  mutate(layer = case_when(
    person_id %in% l1_contacts & person_id %in% l3_contacts ~ "HH + HCW",
    person_id %in% l1_contacts ~ "Layer 1 (HH)",
    person_id %in% l3_contacts ~ "Layer 3 (HCW)",
    TRUE                        ~ "Layer 2 (Community)"
  ))

make_edges <- function(contacts, coord_lkp, sel_j) {
  if (length(contacts) == 0) return(NULL)
  data.frame(to=contacts) %>%
    left_join(coord_lkp, by=c("to"="person_id")) %>%
    rename(x2=xj, y2=yj) %>%
    mutate(x1=sel_j$xj, y1=sel_j$yj)
}
edges_l1 <- make_edges(l1_contacts, coord_lookup, sel_j)
edges_l2 <- make_edges(l2_contacts, coord_lookup, sel_j)
edges_l3 <- make_edges(l3_contacts, coord_lookup, sel_j)

lat0_rad <- mean(nodes$y) * pi / 180
deg_km_x <- 111.32 * cos(lat0_rad)
deg_km_y <- 110.54
dist_km  <- function(contacts, sel, nodes) {
  if (length(contacts) == 0) return(numeric(0))
  cn <- nodes[nodes$person_id %in% contacts, ]
  sqrt(((cn$x - sel$x)*deg_km_x)^2 + ((cn$y - sel$y)*deg_km_y)^2)
}
dist_l2 <- dist_km(l2_contacts, sel, nodes)

# Set extent to fit all contacts + selected person
all_x <- c(sel$x, if(length(all_contacts)>0) nodes$x[nodes$person_id %in% all_contacts])
all_y <- c(sel$y, if(length(all_contacts)>0) nodes$y[nodes$person_id %in% all_contacts])
pad   <- 0.005
xlim  <- c(min(all_x) - pad, max(all_x) + pad)
ylim  <- c(min(all_y) - pad, max(all_y) + pad)

layer_colors <- c("Layer 1 (HH)"        = "#4EC9FF",
                  "Layer 2 (Community)" = "#FF8C42",
                  "Layer 3 (HCW)"       = "#FFD700",
                  "HH + HCW"            = "#FF4EFF")

p_network <- ggplot() +
  geom_raster(data=cell_pop %>% filter(x>=xlim[1],x<=xlim[2],y>=ylim[1],y<=ylim[2]),
              aes(x=x,y=y,fill=log1p(n)), alpha=0.5) +
  scale_fill_viridis_c(option="magma", guide="none") +
  { if (!is.null(edges_l2))
    geom_segment(data=edges_l2, aes(x=x1,y=y1,xend=x2,yend=y2),
                 color="#FF8C42", alpha=0.5, linewidth=0.4) } +
  { if (!is.null(edges_l3))
    geom_segment(data=edges_l3, aes(x=x1,y=y1,xend=x2,yend=y2),
                 color="#FFD700", alpha=0.7, linewidth=0.6) } +
  { if (!is.null(edges_l1))
    geom_segment(data=edges_l1, aes(x=x1,y=y1,xend=x2,yend=y2),
                 color="#4EC9FF", alpha=0.8, linewidth=0.6) } +
  geom_point(data=contact_nodes, aes(x=xj,y=yj,color=layer),
             size=2.0, alpha=0.9) +
  scale_color_manual(values=layer_colors, name="Layer") +
  geom_point(data=sel_j, aes(x=xj,y=yj),
             color="white", size=5, shape=21, fill="red", stroke=1.5) +
  coord_fixed(xlim=xlim, ylim=ylim) +
  labs(title    = sprintf("Person %d — network connections", pid),
       subtitle = sprintf("Age: %d | HCW: %s | L1: %d | L2: %d | L3: %d",
                          sel$age, ifelse(sel$is_hcw,"YES","no"),
                          length(l1_contacts), length(l2_contacts),
                          length(l3_contacts)),
       x="Longitude", y="Latitude") +
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"),
        plot.subtitle=element_text(size=9, color="grey40"),
        legend.position="right")

kernel  <- readRDS(kernel_path)
d_seq   <- seq(0, max(c(dist_l2, 0.1)) * 1.1, length.out=300)
kern_df <- data.frame(
  d = d_seq,
  w = kernel$p_hat * kernel$a1_hat * exp(-kernel$a1_hat * d_seq) +
    (1-kernel$p_hat) * kernel$a2_hat * exp(-kernel$a2_hat * d_seq)
)
kern_df$w_scaled <- kern_df$w / max(kern_df$w)

p_dist <- ggplot() +
  { if (length(dist_l2) > 0)
    geom_histogram(data=data.frame(d=dist_l2),
                   aes(x=d, y=after_stat(density)),
                   bins=30, fill="#FF8C42", alpha=0.7, color="white") } +
  geom_line(data=kern_df,
            aes(x=d, y=w_scaled * max(if(length(dist_l2)>1) density(dist_l2)$y else 1)),
            color="white", linewidth=1.2, linetype="dashed") +
  labs(title    = "Layer 2: Contact distance distribution",
       subtitle = sprintf("%d contacts | kernel (dashed)", length(dist_l2)),
       x="Distance (km)", y="Density") +
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"),
        plot.subtitle=element_text(size=9, color="grey40"))

if (length(l1_contacts) > 0) {
  hh_members <- nodes %>%
    filter(person_id %in% c(pid, l1_contacts)) %>%
    mutate(role=ifelse(person_id==pid,"Selected","HH member"))
  p_hh_age <- ggplot(hh_members, aes(x=age, fill=role)) +
    geom_histogram(bins=20, color="white", alpha=0.8, position="stack") +
    scale_fill_manual(values=c("Selected"="red","HH member"="#4EC9FF")) +
    labs(title    = sprintf("Household (size: %d)", length(l1_contacts)+1L),
         subtitle = sprintf("Ages: %s", paste(sort(hh_members$age), collapse=", ")),
         x="Age", y="Count", fill="") +
    theme_bw() +
    theme(plot.title=element_text(size=12, face="bold"),
          legend.position="top")
} else {
  p_hh_age <- ggplot() + labs(title="No household members") + theme_void()
}

p_person <- p_network / (p_dist + p_hh_age) +
  plot_layout(heights=c(2,1)) +
  plot_annotation(
    title = sprintf("%s — Person %d", tag, pid),
    theme = theme(plot.title=element_text(size=14, face="bold"))
  )

out_p2 <- file.path(fig_dir, sprintf("%s_person%d.png", tag, pid))
ggsave(out_p2, plot=p_person, width=14, height=12, dpi=150)
cat(sprintf("Saved: %s\n", out_p2))

# Console summary
cat("\n=== Contact summary ===\n")
cat(sprintf("  HH contacts: %s\n", paste(l1_contacts, collapse=", ")))
if (length(dist_l2) > 0) {
  cat(sprintf("  Community  : %d | mean=%.2f km | max=%.2f km\n",
              length(dist_l2), mean(dist_l2), max(dist_l2)))
  cat(sprintf("    ≤1km  : %d (%.0f%%)\n",
              sum(dist_l2<=1), 100*mean(dist_l2<=1)))
  cat(sprintf("    1-10km: %d (%.0f%%)\n",
              sum(dist_l2>1 & dist_l2<=10), 100*mean(dist_l2>1 & dist_l2<=10)))
  cat(sprintf("    >10km : %d (%.0f%%)\n",
              sum(dist_l2>10), 100*mean(dist_l2>10)))
}
if (length(l3_contacts) > 0)
  cat(sprintf("  HCW contacts: %d\n", length(l3_contacts)))

# ==============================================================================
# [Part 1c] Prem matrix vs observed contact age structure
# ==============================================================================

prem_dir     <- "data/Prem_contact"
prem_country <- "Congo"

load_prem_community <- function(data_dir, country) {
  file_map <- list(
    work   = "MUestimates_work_1.xlsx",
    school = "MUestimates_school_1.xlsx",
    other  = "MUestimates_other_locations_1.xlsx"
  )
  comm_mat <- matrix(0, 16, 16)
  for (s in names(file_map)) {
    fpath <- file.path(data_dir, file_map[[s]])
    if (!file.exists(fpath)) next
    raw <- readxl::read_excel(fpath, sheet=country, col_names=FALSE)
    raw <- raw[-1, ]
    mat <- matrix(as.numeric(as.matrix(raw)), nrow=16, ncol=16)
    mat <- (mat + t(mat)) / 2
    comm_mat <- comm_mat + mat
  }
  comm_mat
}

prem_comm <- load_prem_community(prem_dir, prem_country)

# Observed contact age structure from layer2
age_grp_labels <- c("0-4","5-9","10-14","15-19","20-24","25-29",
                    "30-34","35-39","40-44","45-49","50-54",
                    "55-59","60-64","65-69","70-74","75+")

# Build age_group per person
age_lookup <- nodes %>% select(person_id, age_group)

obs_contacts <- layer2 %>%
  left_join(age_lookup, by=c("from"="person_id")) %>% rename(ag_from=age_group) %>%
  left_join(age_lookup, by=c("to"="person_id"))   %>% rename(ag_to=age_group) %>%
  filter(!is.na(ag_from), !is.na(ag_to))

# Count observed contacts: symmetric (add both directions)
obs_mat <- matrix(0L, 16, 16)
for (k in seq_len(nrow(obs_contacts))) {
  i <- obs_contacts$ag_from[k]
  j <- obs_contacts$ag_to[k]
  obs_mat[i, j] <- obs_mat[i, j] + 1L
  obs_mat[j, i] <- obs_mat[j, i] + 1L
}

# Normalize both matrices row-wise (proportion of contacts per age group)
normalize_rows <- function(mat) {
  row_sums <- rowSums(mat)
  row_sums[row_sums == 0] <- 1
  mat / row_sums
}

prem_norm <- normalize_rows((prem_comm + t(prem_comm)) / 2)
obs_norm  <- normalize_rows((obs_mat   + t(obs_mat))   / 2)

# Convert to long format for ggplot
mat_to_long <- function(mat, type_label) {
  df <- as.data.frame(mat)
  colnames(df) <- age_grp_labels
  df$from_age  <- age_grp_labels
  df$type      <- type_label
  tidyr::pivot_longer(df, cols=all_of(age_grp_labels),
                      names_to="to_age", values_to="value")
}

library(tidyr)
prem_long <- mat_to_long(prem_norm, "Prem (expected)")
obs_long  <- mat_to_long(obs_norm,  "Observed")

comp_long <- bind_rows(prem_long, obs_long) %>%
  mutate(from_age = factor(from_age, levels=age_grp_labels),
         to_age   = factor(to_age,   levels=age_grp_labels))

# Heatmap comparison
p_prem_comp <- ggplot(comp_long, aes(x=to_age, y=from_age, fill=value)) +
  geom_tile() +
  scale_fill_viridis_c(option="magma", name="Proportion\nof contacts") +
  facet_wrap(~type) +
  labs(title    = sprintf("%s — Contact age structure: Prem vs Observed", tag),
       subtitle = "Row-normalized: proportion of contacts with each age group",
       x="Contact age group", y="Index age group") +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, hjust=1, size=7),
        axis.text.y  = element_text(size=7),
        plot.title   = element_text(size=12, face="bold"),
        plot.subtitle= element_text(size=9, color="grey40"),
        strip.text   = element_text(face="bold"))

out_prem <- file.path(fig_dir, sprintf("%s_prem_comparison.png", tag))
ggsave(out_prem, plot=p_prem_comp, width=12, height=6, dpi=150)
cat(sprintf("Saved: %s\n", out_prem))

# Correlation check
cat(sprintf("\n=== Prem vs Observed correlation ===\n"))
cat(sprintf("  Pearson r: %.4f\n",
            cor(as.vector(prem_norm), as.vector(obs_norm))))
cat(sprintf("  Max abs diff: %.4f\n",
            max(abs(prem_norm - obs_norm))))
