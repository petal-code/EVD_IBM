# ==============================================================================
# COD_IBM_network_build.R
# Purpose:
#   Build all 4 contact network layers for Kinshasa EVD IBM
#   Layer 1: Household (full clique)
#   Layer 2: Community (Prem NB, 10km radius, cell-level distance stored)
#   Layer 3: Healthcare (HCW-HCW clique + admission lookup)
#   Layer 4: Funeral extra (Prem NB, 10km radius, no decay)
#             — duplicates with Layer 1 and Layer 2 pre-removed at build time
#
# Distance note:
#   Bounding box + haversine filtering uses individual coords (indiv_x/y)
#   Stored distance_km uses cell-center coords (cell_x/y) to save memory
#
# Usage:
#   sample_fraction = 1.0  → full network
#   sample_fraction = 0.1  → 10% test run (cell-level sampling)
#
# Outputs (saved to output/network/):
#   {prefix}_nodes.rds            : node dataframe
#   {prefix}_layer1_household.rds : household edges
#   {prefix}_layer2_community.rds : community edges (distance_km stored)
#   {prefix}_layer3_hcw_edges.rds : HCW-HCW edges
#   {prefix}_layer3_admission.rds : admission lookup
#   {prefix}_layer4_funeral.rds   : funeral extra edges (HH + comm duplicates removed)
# ==============================================================================

library(dplyr)
library(tidyr)
library(readxl)
library(parallel)

# ==============================================================================
# [Configuration] — Modify only this section
# ==============================================================================

province_name   <- "Kinshasa"
prem_dir        <- "data/Prem_contact"
prem_country    <- "Congo"
hf_path         <- "data/COD_GRID3_health_facilities_v8.csv"
synpop_path     <- "output/household/Kinshasa_synthetic_population.rds"
output_dir      <- "output/network"

# ── Sampling fraction ─────────────────────────────────────────
# 1.0 = full network | 0.1 = 10% test run
# Sampling is done at the CELL level to preserve spatial structure
sample_fraction <- 1

# Network parameters
hcw_rate        <- 13.78 / 10000  # HCWs per total population (DRC)
comm_nb_size    <- 0.1             # NB dispersion for community contacts
funeral_nb_size <- 0.1             # NB dispersion for funeral extra contacts
radius_km       <- 10              # Contact radius in km

n_cores         <- min(15L, max(1L, detectCores() - 1L))
network_seed    <- 42

# Output file prefix
is_test    <- sample_fraction < 1.0
file_tag   <- if (is_test) sprintf("_test%02d", round(sample_fraction * 100)) else ""
out_prefix <- sprintf("%s%s", province_name, file_tag)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
t_start <- proc.time()

message("=== COD IBM Network Builder ===")
message(sprintf("  Province        : %s", province_name))
message(sprintf("  Sample fraction : %.0f%%", sample_fraction * 100))
message(sprintf("  Output prefix   : %s", out_prefix))
message(sprintf("  Cores           : %d", n_cores))

# ==============================================================================
# [Section 1] Load and prepare data
# ==============================================================================

message("\n=== Section 1: Loading data ===")

pop_raw <- readRDS(synpop_path)

pop_full <- pop_raw %>%
  mutate(
    person_id = row_number(),
    age_group = pmin(floor(age / 5L) + 1L, 16L),  # Prem 16 age groups
    is_adult  = age >= 18L
  )

message(sprintf("  Full population   : %d individuals", nrow(pop_full)))

# ── Subsample at cell level if sample_fraction < 1 ────────────
if (is_test) {
  set.seed(network_seed)
  sampled_cells <- pop_full %>%
    distinct(cell_x, cell_y) %>%
    slice_sample(prop = sample_fraction)

  pop_nodes <- pop_full %>%
    semi_join(sampled_cells, by = c("cell_x", "cell_y")) %>%
    mutate(person_id = row_number())  # Reassign sequential person_id

  message(sprintf("  Sampled cells     : %d / %d",
                  nrow(sampled_cells),
                  n_distinct(pop_full$cell_x, pop_full$cell_y)))
  message(sprintf("  Sampled pop       : %d (%.1f%%)",
                  nrow(pop_nodes), 100 * nrow(pop_nodes) / nrow(pop_full)))
} else {
  pop_nodes <- pop_full
  message(sprintf("  Using full population: %d individuals", nrow(pop_nodes)))
}

rm(pop_raw, pop_full); gc()
message(sprintf("  Adults (>= 18)    : %d", sum(pop_nodes$is_adult)))

# Load Prem matrices
message("  Loading Prem contact matrices...")
load_prem_matrices <- function(data_dir, country) {
  file_map <- list(
    home   = "MUestimates_home_1.xlsx",
    work   = "MUestimates_work_1.xlsx",
    school = "MUestimates_school_1.xlsx",
    other  = "MUestimates_other_locations_1.xlsx"
  )
  load_one <- function(setting) {
    fpath <- file.path(data_dir, file_map[[setting]])
    if (!file.exists(fpath)) { warning(sprintf("Missing: %s", fpath)); return(NULL) }
    raw <- read_excel(fpath, sheet = country, col_names = FALSE)
    raw <- raw[-1, ]  # Drop header row (X1, X2, ...)
    mat <- matrix(as.numeric(as.matrix(raw)), nrow = 16, ncol = 16)
    (mat + t(mat)) / 2  # Symmetrize
  }
  mats <- lapply(setNames(names(file_map), names(file_map)), load_one)
  comm_mat <- matrix(0, 16, 16)
  for (s in c("work", "school", "other"))
    if (!is.null(mats[[s]])) comm_mat <- comm_mat + mats[[s]]
  mats$community <- comm_mat
  mats
}

prem_mats    <- load_prem_matrices(prem_dir, prem_country)
prem_comm_mu <- rowSums(prem_mats$community)

message(sprintf("  Prem community mu range: %.2f - %.2f",
                min(prem_comm_mu), max(prem_comm_mu)))

# Load hospital data
message("  Loading hospital data...")
hf_all <- read.csv(hf_path)
hf_hospital <- hf_all %>%
  filter(province == province_name, !is.na(lon), !is.na(lat)) %>%
  mutate(facility_group = case_when(
    esstype %in% c("Hôpital", "Hôpital Général de Référence",
                   "Centre Hopitalier") ~ "Hospital",
    TRUE ~ "Other"
  )) %>%
  filter(facility_group == "Hospital")

message(sprintf("  Hospitals: %d", nrow(hf_hospital)))

# ==============================================================================
# [Section 2] Assign HCWs and nearest hospitals
# ==============================================================================

message("\n=== Section 2: HCW assignment ===")

set.seed(network_seed)
n_hcw     <- round(nrow(pop_nodes) * hcw_rate)
adult_ids <- pop_nodes$person_id[pop_nodes$is_adult]
hcw_ids   <- sample(adult_ids, size = n_hcw, replace = FALSE)

pop_nodes <- pop_nodes %>%
  mutate(is_hcw = person_id %in% hcw_ids)

message(sprintf("  HCWs assigned: %d (%.4f%%)",
                sum(pop_nodes$is_hcw), 100 * mean(pop_nodes$is_hcw)))

haversine_km_vec <- function(lon1, lat1, lon2_vec, lat2_vec) {
  R   <- 6371
  phi <- (lat2_vec - lat1) * pi / 180
  lam <- (lon2_vec - lon1) * pi / 180
  a   <- sin(phi/2)^2 + cos(lat1*pi/180) * cos(lat2_vec*pi/180) * sin(lam/2)^2
  2 * R * asin(sqrt(a))
}

hosp_coords <- as.matrix(hf_hospital[, c("lon", "lat")])

message("  Assigning HCWs to nearest hospitals...")
hcw_nodes <- pop_nodes %>% filter(is_hcw)
nearest_hosp_hcw <- apply(
  as.matrix(hcw_nodes[, c("indiv_x", "indiv_y")]), 1, function(coord) {
    which.min(haversine_km_vec(coord[1], coord[2],
                               hosp_coords[, 1], hosp_coords[, 2]))
  })
hcw_nodes <- hcw_nodes %>%
  mutate(hospital_id = hf_hospital$OBJECTID[nearest_hosp_hcw])

message("  Assigning non-HCWs to nearest hospitals...")
non_hcw_nodes <- pop_nodes %>% filter(!is_hcw)
nearest_hosp_non <- apply(
  as.matrix(non_hcw_nodes[, c("indiv_x", "indiv_y")]), 1, function(coord) {
    which.min(haversine_km_vec(coord[1], coord[2],
                               hosp_coords[, 1], hosp_coords[, 2]))
  })
non_hcw_nodes <- non_hcw_nodes %>%
  mutate(hospital_id = hf_hospital$OBJECTID[nearest_hosp_non])

pop_nodes <- pop_nodes %>%
  left_join(
    bind_rows(
      hcw_nodes     %>% select(person_id, hospital_id),
      non_hcw_nodes %>% select(person_id, hospital_id)
    ),
    by = "person_id"
  )

message(sprintf("  Hospitals with HCWs: %d / %d",
                length(unique(hcw_nodes$hospital_id)), nrow(hf_hospital)))

# ==============================================================================
# [Section 3] Save node dataframe
# ==============================================================================

message("\n=== Section 3: Saving nodes ===")

nodes_save <- pop_nodes %>%
  select(person_id, hh_id, hospital_id, age, age_group,
         is_hcw, is_adult, cell_x, cell_y, indiv_x, indiv_y)

out_path <- file.path(output_dir, sprintf("%s_nodes.rds", out_prefix))
saveRDS(nodes_save, out_path)
message(sprintf("  Saved: %s (%.1f MB)", basename(out_path),
                file.info(out_path)$size / 1e6))

# ==============================================================================
# [Section 4] Layer 1 — Household edges
# ==============================================================================

message("\n=== Section 4: Layer 1 — Household edges ===")
t0 <- proc.time()

hh_edges <- lapply(
  split(pop_nodes$person_id, pop_nodes$hh_id),
  function(members) {
    if (length(members) < 2) return(NULL)
    pairs <- combn(members, 2)
    data.frame(from = pairs[1, ], to = pairs[2, ],
               stringsAsFactors = FALSE)
  }) %>%
  Filter(Negate(is.null), .) %>%
  bind_rows()

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
message(sprintf("  Household edges: %d (%.1f sec)", nrow(hh_edges), elapsed))

out_path <- file.path(output_dir, sprintf("%s_layer1_household.rds", out_prefix))
saveRDS(hh_edges, out_path)
message(sprintf("  Saved: %s (%.1f MB)", basename(out_path),
                file.info(out_path)$size / 1e6))

# Keep hh_edges in memory for funeral deduplication — rm() after Layer 4

# ==============================================================================
# [Section 5] Layer 2 — Community edges
# Filtering  : individual coords (indiv_x/y) for accurate 10km radius
# Stored dist: cell coords (cell_x/y) to reduce memory
# ==============================================================================

message("\n=== Section 5: Layer 2 — Community edges ===")

create_comm_edges_parallel <- function(node_df, prem_mu, size = 0.1,
                                       radius_km = 10, n_cores = 15L, seed = 42) {
  N            <- nrow(node_df)
  coords_indiv <- as.matrix(node_df[, c("indiv_x", "indiv_y")])
  coords_cell  <- as.matrix(node_df[, c("cell_x",  "cell_y")])
  hh_ids       <- node_df$hh_id
  age_grps     <- node_df$age_group
  chunks       <- split(seq_len(N), ceiling(seq_len(N) / ceiling(N / n_cores)))

  worker <- function(chunk_ids, coords_indiv, coords_cell, hh_ids, age_grps,
                     prem_mu, size, radius_km, seed) {
    hav_km <- function(lon1, lat1, lon2, lat2) {
      R <- 6371
      a <- sin(((lat2-lat1)*pi/180)/2)^2 +
        cos(lat1*pi/180)*cos(lat2*pi/180)*sin(((lon2-lon1)*pi/180)/2)^2
      2*R*asin(sqrt(a))
    }
    set.seed(seed)
    lat_rad <- radius_km / 111.0
    from_v  <- integer(0); to_v <- integer(0); dist_v <- numeric(0)

    for (i in chunk_ids) {
      mu_i <- prem_mu[age_grps[i]]
      n_c  <- rnbinom(1, size = size, mu = mu_i)
      if (n_c == 0) next

      lon_i   <- coords_indiv[i, 1]; lat_i <- coords_indiv[i, 2]
      lon_rad <- radius_km / (111.0 * cos(lat_i * pi / 180))

      cands <- which(
        abs(coords_indiv[, 2] - lat_i) <= lat_rad &
          abs(coords_indiv[, 1] - lon_i) <= lon_rad &
          hh_ids != hh_ids[i] &
          seq_len(nrow(coords_indiv)) != i
      )
      if (length(cands) == 0) next

      dists_filter <- hav_km(lon_i, lat_i,
                             coords_indiv[cands, 1], coords_indiv[cands, 2])
      cands <- cands[dists_filter <= radius_km]
      if (length(cands) == 0) next

      n_s     <- min(n_c, length(cands))
      idx_s   <- sample(length(cands), n_s, replace = FALSE)
      sampled <- cands[idx_s]

      # Keep i < j to avoid within-layer duplicates
      keep <- sampled > i
      if (any(keep)) {
        sampled_keep <- sampled[keep]
        from_v <- c(from_v, rep(i, sum(keep)))
        to_v   <- c(to_v, sampled_keep)
        # Store cell-level distance (reduces unique values → better compression)
        dist_v <- c(dist_v, hav_km(
          coords_cell[i, 1],            coords_cell[i, 2],
          coords_cell[sampled_keep, 1], coords_cell[sampled_keep, 2]
        ))
      }
    }
    if (length(from_v) == 0) return(NULL)
    data.frame(from = from_v, to = to_v, distance_km = dist_v,
               stringsAsFactors = FALSE)
  }

  cl <- makeCluster(n_cores, type = "PSOCK")
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl,
                varlist = c("coords_indiv", "coords_cell", "hh_ids",
                            "age_grps", "prem_mu", "size", "radius_km"),
                envir = environment())
  chunk_args <- lapply(seq_along(chunks), function(ci)
    list(chunk_ids = chunks[[ci]], seed = seed + ci))

  results <- parLapply(cl, chunk_args, function(args)
    worker(args$chunk_ids, coords_indiv, coords_cell, hh_ids, age_grps,
           prem_mu, size, radius_km, args$seed))

  out <- bind_rows(Filter(Negate(is.null), results))
  out[!duplicated(paste(out$from, out$to)), ]
}

t0         <- proc.time()
comm_edges <- create_comm_edges_parallel(
  pop_nodes, prem_comm_mu,
  size = comm_nb_size, radius_km = radius_km,
  n_cores = n_cores, seed = network_seed
)
elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
message(sprintf("  Community edges: %d (%.1f sec / %.1f min)",
                nrow(comm_edges), elapsed, elapsed / 60))
message(sprintf("  Mean distance  : %.2f km", mean(comm_edges$distance_km)))

out_path <- file.path(output_dir, sprintf("%s_layer2_community.rds", out_prefix))
saveRDS(comm_edges, out_path)
message(sprintf("  Saved: %s (%.1f MB)", basename(out_path),
                file.info(out_path)$size / 1e6))

# Keep comm_edges in memory for funeral deduplication — rm() after Layer 4

# ==============================================================================
# [Section 6] Layer 3 — Healthcare edges
# ==============================================================================

message("\n=== Section 6: Layer 3 — Healthcare edges ===")

# (a) HCW-HCW edges within same hospital (always active)
hcw_hcw_edges <- hcw_nodes %>%
  group_by(hospital_id) %>%
  group_map(~ {
    members <- .x$person_id
    if (length(members) < 2) return(NULL)
    pairs <- combn(members, 2)
    data.frame(from        = pairs[1, ],
               to          = pairs[2, ],
               hospital_id = .y$hospital_id,
               stringsAsFactors = FALSE)
  }) %>%
  bind_rows()

message(sprintf("  HCW-HCW edges     : %d", nrow(hcw_hcw_edges)))
message(sprintf("  Hospitals covered : %d / %d",
                length(unique(hcw_hcw_edges$hospital_id)), nrow(hf_hospital)))

out_path <- file.path(output_dir, sprintf("%s_layer3_hcw_edges.rds", out_prefix))
saveRDS(hcw_hcw_edges, out_path)
message(sprintf("  Saved: %s (%.1f MB)", basename(out_path),
                file.info(out_path)$size / 1e6))
rm(hcw_hcw_edges); gc()

# (b) Admission lookup: person_id + hh_id + hospital_id + HCW list
hospital_hcw_lookup <- hcw_nodes %>%
  group_by(hospital_id) %>%
  summarise(hcw_list = list(person_id), .groups = "drop")

admission_lookup <- non_hcw_nodes %>%
  select(person_id, hh_id, hospital_id) %>%
  left_join(hospital_hcw_lookup, by = "hospital_id")

message(sprintf("  Admission lookup  : %d individuals", nrow(admission_lookup)))
message(sprintf("  Mean HCWs/hospital: %.1f",
                mean(sapply(hospital_hcw_lookup$hcw_list, length))))

out_path <- file.path(output_dir, sprintf("%s_layer3_admission.rds", out_prefix))
saveRDS(admission_lookup, out_path)
message(sprintf("  Saved: %s (%.1f MB)", basename(out_path),
                file.info(out_path)$size / 1e6))
rm(admission_lookup, hospital_hcw_lookup); gc()

# ==============================================================================
# [Section 7] Layer 4 — Funeral extra edges
# Prem community NB, 10km radius, no distance decay
# Activated on death: added on top of existing HH + community edges
# Pre-deduplicated against Layer 1 (household) and Layer 2 (community)
# ==============================================================================

message("\n=== Section 7: Layer 4 — Funeral extra edges ===")

create_funeral_edges_parallel <- function(node_df, prem_mu, size = 0.1,
                                          radius_km = 10, n_cores = 15L, seed = 77) {
  N            <- nrow(node_df)
  coords_indiv <- as.matrix(node_df[, c("indiv_x", "indiv_y")])
  hh_ids       <- node_df$hh_id
  age_grps     <- node_df$age_group
  chunks       <- split(seq_len(N), ceiling(seq_len(N) / ceiling(N / n_cores)))

  worker <- function(chunk_ids, coords_indiv, hh_ids, age_grps,
                     prem_mu, size, radius_km, seed) {
    hav_km <- function(lon1, lat1, lon2, lat2) {
      R <- 6371
      a <- sin(((lat2-lat1)*pi/180)/2)^2 +
        cos(lat1*pi/180)*cos(lat2*pi/180)*sin(((lon2-lon1)*pi/180)/2)^2
      2*R*asin(sqrt(a))
    }
    set.seed(seed)
    lat_rad <- radius_km / 111.0
    from_v  <- integer(0)
    to_v    <- integer(0)

    for (i in chunk_ids) {
      mu_i <- prem_mu[age_grps[i]]
      n_c  <- rnbinom(1, size = size, mu = mu_i)
      if (n_c == 0) next

      lon_i   <- coords_indiv[i, 1]; lat_i <- coords_indiv[i, 2]
      lon_rad <- radius_km / (111.0 * cos(lat_i * pi / 180))

      cands <- which(
        abs(coords_indiv[, 2] - lat_i) <= lat_rad &
          abs(coords_indiv[, 1] - lon_i) <= lon_rad &
          hh_ids != hh_ids[i] &
          seq_len(nrow(coords_indiv)) != i
      )
      if (length(cands) == 0) next

      dists <- hav_km(lon_i, lat_i,
                      coords_indiv[cands, 1], coords_indiv[cands, 2])
      cands <- cands[dists <= radius_km]
      if (length(cands) == 0) next

      n_s     <- min(n_c, length(cands))
      sampled <- sample(cands, n_s, replace = FALSE)

      # Keep i < j to avoid within-layer duplicates
      keep <- sampled > i
      if (any(keep)) {
        from_v <- c(from_v, rep(i, sum(keep)))
        to_v   <- c(to_v, sampled[keep])
      }
    }
    if (length(from_v) == 0) return(NULL)
    data.frame(from = from_v, to = to_v, stringsAsFactors = FALSE)
  }

  cl <- makeCluster(n_cores, type = "PSOCK")
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl,
                varlist = c("coords_indiv", "hh_ids", "age_grps",
                            "prem_mu", "size", "radius_km"),
                envir = environment())
  chunk_args <- lapply(seq_along(chunks), function(ci)
    list(chunk_ids = chunks[[ci]], seed = seed + ci))

  results <- parLapply(cl, chunk_args, function(args)
    worker(args$chunk_ids, coords_indiv, hh_ids, age_grps,
           prem_mu, size, radius_km, args$seed))

  out <- bind_rows(Filter(Negate(is.null), results))
  out[!duplicated(paste(out$from, out$to)), ]
}

t0            <- proc.time()
funeral_edges <- create_funeral_edges_parallel(
  pop_nodes, prem_comm_mu,
  size = funeral_nb_size, radius_km = radius_km,
  n_cores = n_cores, seed = network_seed + 35
)
elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
message(sprintf("  Funeral extra edges (raw)  : %d (%.1f sec / %.1f min)",
                nrow(funeral_edges), elapsed, elapsed / 60))

# ── Remove duplicates with Layer 1 (household) and Layer 2 (community) ────────
# Build lookup keys from existing layers: paste(min, max) for undirected edges
message("  Removing duplicates with Layer 1 and Layer 2...")

# Canonical edge key: always smaller id first
edge_key <- function(df) paste(pmin(df$from, df$to), pmax(df$from, df$to), sep = "_")

existing_keys <- c(edge_key(hh_edges),
                   edge_key(comm_edges))

funeral_keys  <- edge_key(funeral_edges)
funeral_edges <- funeral_edges[!funeral_keys %in% existing_keys, ]

message(sprintf("  Funeral extra edges (dedup): %d (removed %d duplicates)",
                nrow(funeral_edges),
                length(funeral_keys) - nrow(funeral_edges)))

out_path <- file.path(output_dir, sprintf("%s_layer4_funeral.rds", out_prefix))
saveRDS(funeral_edges, out_path)
message(sprintf("  Saved: %s (%.1f MB)", basename(out_path),
                file.info(out_path)$size / 1e6))

# Free all remaining large objects
rm(hh_edges, comm_edges, funeral_edges); gc()

# ==============================================================================
# [Done] Summary
# ==============================================================================

elapsed_total <- round((proc.time() - t_start)[["elapsed"]], 1)

message("\n=== Network Build Complete ===")
message(sprintf("  Province        : %s", province_name))
message(sprintf("  Sample fraction : %.0f%%", sample_fraction * 100))
message(sprintf("  Nodes           : %d", nrow(pop_nodes)))
message(sprintf("  Total time      : %.1f sec (%.1f min)",
                elapsed_total, elapsed_total / 60))
message("\n  Output files:")
for (f in list.files(output_dir, pattern = out_prefix)) {
  size_mb <- round(file.info(file.path(output_dir, f))$size / 1e6, 1)
  message(sprintf("    %-55s %6.1f MB", f, size_mb))
}
