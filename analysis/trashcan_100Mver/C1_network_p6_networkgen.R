# ==============================================================================
# C1_network_p6_network_build.R
# Purpose:
#   For each Level 2 territory:
#     1. Load synthetic population (from p4)
#     2. Build cell-cell sparse distance matrix
#     3. Build Layer 1: household edges (full clique)
#     4. Build Layer 2: community edges (Prem NB + kernel-weighted sampling)
#     5. Build Layer 3: healthcare edges (HCW-HCW + admission lookup)
#     6. Save all network layers
#
# Inputs:
#   output/household/{NAME1}_{NAME2}_synthetic_population.rds  (from p4)
#   output/kernel/community_distance_kernel.rds                (from p5_kernel)
#   data/COD_GRID3_health_facilities_v8.csv
#   data/Prem_contact/
#
# Outputs (per territory):
#   output/network/{NAME1}_{NAME2}_nodes.rds
#   output/network/{NAME1}_{NAME2}_cell_dist.rds
#   output/network/{NAME1}_{NAME2}_layer1_household.rds
#   output/network/{NAME1}_{NAME2}_layer2_community.rds
#   output/network/{NAME1}_{NAME2}_layer3_hcw_edges.rds
#   output/network/{NAME1}_{NAME2}_layer3_admission.rds
# ==============================================================================

library(dplyr)
library(readxl)
library(dbscan)
library(Matrix)
library(Rcpp)
library(data.table)
library(sf)

# ==============================================================================
# [Configuration]
# ==============================================================================

prem_dir   <- "data/Prem_contact"
prem_country <- "Congo"
hf_path    <- "data/COD_GRID3_health_facilities_v8.csv"
synpop_dir <- "output/household"
kernel_path <- "output/kernel/community_distance_kernel.rds"
output_dir <- "output/network"
shp_path   <- "data/shpmap/gadm41_COD_2.shp"

hcw_rate     <- 13.78 / 10000  # HCWs per total population (DRC)
comm_nb_size <- 0.1             # NB dispersion
radius_km    <- 10              # Contact radius in km
radius_m     <- radius_km * 1000
network_seed <- 42
chunk_size   <- 100000L

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# [Section 1] Load shared data (once, reused across territories)
# ==============================================================================

message("=== Section 1: Loading shared data ===")

# Prem community contact matrix
load_prem_matrices <- function(data_dir, country) {
  file_map <- list(
    home   = "MUestimates_home_1.xlsx",
    work   = "MUestimates_work_1.xlsx",
    school = "MUestimates_school_1.xlsx",
    other  = "MUestimates_other_locations_1.xlsx"
  )
  load_one <- function(setting) {
    fpath <- file.path(data_dir, file_map[[setting]])
    if (!file.exists(fpath)) return(NULL)
    raw <- read_excel(fpath, sheet = country, col_names = FALSE)
    raw <- raw[-1, ]
    mat <- matrix(as.numeric(as.matrix(raw)), nrow = 16, ncol = 16)
    (mat + t(mat)) / 2
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

# Distance kernel parameters (from p5_kernel)
kernel      <- readRDS(kernel_path)
p_hat       <- kernel$p_hat
a1_hat      <- kernel$a1_hat
a2_hat      <- kernel$a2_hat
message(sprintf("  Kernel: p=%.4f a1=%.4f a2=%.4f", p_hat, a1_hat, a2_hat))

# Hospital data (all DRC)
hf_all <- read.csv(hf_path)
hf_hospital_all <- hf_all %>%
  filter(
    esstype %in% c("Hôpital", "Hôpital Général de Référence", "Centre Hopitalier"),
    !is.na(lon), !is.na(lat)
  )
message(sprintf("  Hospitals (DRC): %d", nrow(hf_hospital_all)))

# Level 2 shapefile for territory-hospital spatial join
cod2_sf  <- st_as_sf(vect(shp_path))
hf_sf    <- st_as_sf(hf_hospital_all, coords = c("lon", "lat"), crs = 4326)
hf_joined <- st_join(hf_sf, cod2_sf[, c("NAME_1", "NAME_2")]) %>%
  st_drop_geometry()

# Compile C++ edge builder
message("  Compiling C++ edge builder...")
sourceCpp("function/build_edges.cpp")
message("  Compiled successfully")

# ==============================================================================
# [Section 2] List synthetic population files
# ==============================================================================

message("\n=== Section 2: Finding synthetic population files ===")

pop_files <- list.files(synpop_dir,
                        pattern = "_synthetic_population\\.rds$",
                        full.names = TRUE)
message(sprintf("  Territories found: %d", length(pop_files)))

# ==============================================================================
# [Helper functions]
# ==============================================================================

# Build cell-cell sparse distance matrix for a territory
build_cell_dist <- function(cell_lookup, radius_m) {
  cell_coords <- as.matrix(cell_lookup[, c("cell_x", "cell_y")])
  lat0_rad    <- mean(cell_coords[, 2]) * pi / 180
  cell_mx     <- cell_coords[, 1] * 111320 * cos(lat0_rad)
  cell_my     <- cell_coords[, 2] * 110540
  cell_m      <- cbind(cell_mx, cell_my)

  nn       <- frNN(cell_m, eps = radius_m)
  nn_id    <- nn$id
  n_cells  <- nrow(cell_lookup)

  from_l <- vector("list", n_cells)
  to_l   <- vector("list", n_cells)
  dist_l <- vector("list", n_cells)
  k      <- 0L

  for (ci in seq_len(n_cells)) {
    nbrs <- nn_id[[ci]]
    nbrs <- nbrs[nbrs > ci]  # Upper triangle only
    if (length(nbrs) == 0) next

    dx    <- cell_mx[nbrs] - cell_mx[ci]
    dy    <- cell_my[nbrs] - cell_my[ci]
    dists <- sqrt(dx^2 + dy^2) / 1000  # km

    keep <- dists <= radius_km
    if (!any(keep)) next

    k           <- k + 1L
    from_l[[k]] <- rep(ci, sum(keep))
    to_l[[k]]   <- nbrs[keep]
    dist_l[[k]] <- dists[keep]
  }

  if (k == 0L) return(list(
    dist_mat    = sparseMatrix(i=integer(0), j=integer(0), x=numeric(0),
                               dims=c(n_cells, n_cells)),
    cell_lookup = cell_lookup,
    nn_cell_id  = nn_id,
    cell_mx     = cell_mx,
    cell_my     = cell_my,
    from_vec    = integer(0),
    to_vec      = integer(0),
    dist_vec    = numeric(0)
  ))

  from_c <- unlist(from_l[seq_len(k)], use.names = FALSE)
  to_c   <- unlist(to_l[seq_len(k)],   use.names = FALSE)
  dist_c <- unlist(dist_l[seq_len(k)], use.names = FALSE)

  dist_mat <- sparseMatrix(i = from_c, j = to_c, x = dist_c,
                           dims = c(n_cells, n_cells))

  list(dist_mat    = dist_mat,
       cell_lookup = cell_lookup,
       nn_cell_id  = nn_id,
       cell_mx     = cell_mx,
       cell_my     = cell_my,
       from_vec    = from_c,   # Pass to C++ directly
       to_vec      = to_c,
       dist_vec    = dist_c)
}

# ==============================================================================
# [Section 3] Process each territory
# ==============================================================================

message("\n=== Section 3: Building networks ===")

n_total   <- length(pop_files)
t_start   <- proc.time()[["elapsed"]]
n_saved   <- 0L
n_skipped <- 0L

for (fi in seq_along(pop_files)) {

  pop_raw <- readRDS(pop_files[fi])

  # Parse territory name from filename
  fname  <- basename(pop_files[fi])
  tag    <- sub("_synthetic_population\\.rds$", "", fname)

  # Check if all outputs already exist
  out_nodes   <- file.path(output_dir, sprintf("%s_nodes.rds",            tag))
  out_cdist   <- file.path(output_dir, sprintf("%s_cell_dist.rds",        tag))
  out_layer1  <- file.path(output_dir, sprintf("%s_layer1_household.rds", tag))
  out_layer2  <- file.path(output_dir, sprintf("%s_layer2_community.rds", tag))
  out_layer3h <- file.path(output_dir, sprintf("%s_layer3_hcw_edges.rds", tag))
  out_layer3a <- file.path(output_dir, sprintf("%s_layer3_admission.rds", tag))

  all_exist <- all(file.exists(c(out_nodes, out_cdist, out_layer1,
                                 out_layer2, out_layer3h, out_layer3a)))
  if (all_exist) {
    n_skipped <- n_skipped + 1L; next
  }

  message(sprintf("\n  [%d/%d] %s", fi, n_total, tag))

  # ── Prepare node dataframe ─────────────────────────────────────────────────
  pop_nodes <- pop_raw %>%
    mutate(
      person_id = row_number(),
      age_group = pmin(floor(age / 5L) + 1L, 16L),
      is_adult  = age >= 18L
    )

  # Assign cell_id
  cell_lookup <- pop_nodes %>%
    distinct(cell_x, cell_y) %>%
    arrange(cell_x, cell_y) %>%
    mutate(cell_id = row_number())

  pop_nodes <- pop_nodes %>%
    left_join(cell_lookup, by = c("cell_x", "cell_y"))

  n_cells <- nrow(cell_lookup)

  # ── Build cell distance matrix ─────────────────────────────────────────────
  cdist <- build_cell_dist(cell_lookup, radius_m)

  # ── Assign HCWs ────────────────────────────────────────────────────────────
  set.seed(network_seed + fi)
  n_hcw     <- max(1L, round(nrow(pop_nodes) * hcw_rate))
  adult_ids <- pop_nodes$person_id[pop_nodes$is_adult]

  if (length(adult_ids) == 0) {
    message(sprintf("    SKIP: no adults in %s", tag))
    n_skipped <- n_skipped + 1L; next
  }

  hcw_ids   <- sample(adult_ids, size = min(n_hcw, length(adult_ids)),
                      replace = FALSE)
  pop_nodes <- pop_nodes %>% mutate(is_hcw = person_id %in% hcw_ids)

  # ── Assign nearest hospital at cell level ──────────────────────────────────
  # Use only hospitals in this territory; fallback to nearest DRC hospital
  name1 <- pop_raw$cell_x[1]  # Placeholder; parse from tag below
  # Get hospitals for this territory
  terr_hospitals <- hf_joined %>%
    filter(paste0(gsub("[^A-Za-z0-9]", "_", NAME_1), "_",
                  gsub("[^A-Za-z0-9]", "_", NAME_2)) == tag)

  if (nrow(terr_hospitals) == 0) {
    # No hospitals in territory — use nearest DRC hospital
    terr_hospitals <- hf_hospital_all
  }

  lat0_rad <- mean(cell_lookup$cell_y) * pi / 180
  hosp_mx  <- terr_hospitals$lon * 111320 * cos(lat0_rad)
  hosp_my  <- terr_hospitals$lat * 110540
  cell_mx  <- cdist$cell_mx
  cell_my  <- cdist$cell_my

  cell_hosp_idx <- max.col(
    -(outer(cell_mx, hosp_mx, "-")^2 +
        outer(cell_my, hosp_my, "-")^2)
  )

  cell_lookup <- cell_lookup %>%
    mutate(hospital_id = terr_hospitals$OBJECTID[cell_hosp_idx])

  pop_nodes <- pop_nodes %>%
    left_join(cell_lookup %>% select(cell_id, hospital_id), by = "cell_id")

  hcw_nodes     <- pop_nodes %>% filter(is_hcw)
  non_hcw_nodes <- pop_nodes %>% filter(!is_hcw)

  # ── Save nodes ─────────────────────────────────────────────────────────────
  nodes_save <- pop_nodes %>%
    select(person_id, hh_id, cell_id, hospital_id,
           age, age_group, is_hcw, is_adult, cell_x, cell_y)
  saveRDS(nodes_save, out_nodes)

  # ── Save cell distance matrix ───────────────────────────────────────────────
  saveRDS(list(dist_mat    = cdist$dist_mat,
               cell_lookup = cell_lookup),
          out_cdist)

  # ── Layer 1: Household edges ───────────────────────────────────────────────
  hh_edges <- lapply(
    split(pop_nodes$person_id, pop_nodes$hh_id),
    function(members) {
      if (length(members) < 2) return(NULL)
      pairs <- combn(members, 2)
      data.frame(from = pairs[1,], to = pairs[2,], stringsAsFactors = FALSE)
    }) %>% Filter(Negate(is.null), .) %>% bind_rows()

  saveRDS(hh_edges, out_layer1)
  rm(hh_edges); gc()

  # ── Layer 2: Community edges (Rcpp + kernel-weighted sampling) ─────────────
  cell_sizes   <- tabulate(pop_nodes$cell_id, nbins = n_cells)
  max_per_cell <- max(cell_sizes)
  cell_mat     <- matrix(0L, nrow = n_cells, ncol = max_per_cell)
  cell_fill    <- integer(n_cells)
  for (idx in seq_len(nrow(pop_nodes))) {
    ci                          <- pop_nodes$cell_id[idx]
    cell_fill[ci]               <- cell_fill[ci] + 1L
    cell_mat[ci, cell_fill[ci]] <- idx
  }
  rm(cell_fill)

  hh_ids_int <- as.integer(as.factor(pop_nodes$hh_id))
  cell_ids_int <- as.integer(pop_nodes$cell_id)
  age_grps     <- pop_nodes$age_group
  N            <- nrow(pop_nodes)
  n_chunks_run <- ceiling(N / chunk_size)
  chunk_results <- vector("list", n_chunks_run)

  t_layer2 <- proc.time()[["elapsed"]]
  for (ci in seq_len(n_chunks_run)) {
    idx_start <- (ci - 1L) * chunk_size + 1L
    idx_end   <- min(ci * chunk_size, N)
    chunk_ids <- idx_start:idx_end

    set.seed(network_seed + ci)
    chunk_mu <- prem_comm_mu[age_grps[chunk_ids]]
    n_draws  <- as.integer(rnbinom(length(chunk_ids), size = comm_nb_size,
                                   mu = chunk_mu))

    active_mask  <- n_draws > 0L
    if (!any(active_mask)) next

    active_ids   <- chunk_ids[active_mask]
    active_draws <- n_draws[active_mask]

    chunk_edges <- build_edges_cpp(
      active_ids    = active_ids,
      active_draws  = active_draws,
      hh_ids        = hh_ids_int,
      cell_ids      = cell_ids_int,
      cell_mat      = cell_mat,
      cell_sizes    = cell_sizes,
      nn_cell_id    = cdist$nn_cell_id,
      cell_dist_vec = cdist$dist_vec,
      cell_dist_from = as.integer(cdist$from_vec),
      cell_dist_to   = as.integer(cdist$to_vec),
      p_hat         = p_hat,
      a1_hat        = a1_hat,
      a2_hat        = a2_hat,
      seed          = network_seed + ci * 999
    )

    if (nrow(chunk_edges) > 0)
      chunk_results[[ci]] <- as.data.table(chunk_edges)
  }

  comm_edges <- rbindlist(chunk_results[!sapply(chunk_results, is.null)])
  if (nrow(comm_edges) > 0) {
    comm_edges <- unique(comm_edges, by = c("from", "to"))
    setDF(comm_edges)
  }

  elapsed_l2 <- round(proc.time()[["elapsed"]] - t_layer2, 1)
  message(sprintf("    Layer 2: %d edges (%.1f sec)", nrow(comm_edges), elapsed_l2))
  saveRDS(comm_edges, out_layer2)
  rm(comm_edges, cell_mat, cell_sizes, chunk_results); gc()

  # ── Layer 3a: HCW-HCW edges ────────────────────────────────────────────────
  hcw_hcw_edges <- hcw_nodes %>%
    group_by(hospital_id) %>%
    group_map(~ {
      members <- .x$person_id
      if (length(members) < 2) return(NULL)
      pairs <- combn(members, 2)
      data.frame(from = pairs[1,], to = pairs[2,],
                 hospital_id = .y$hospital_id,
                 stringsAsFactors = FALSE)
    }) %>% bind_rows()

  saveRDS(hcw_hcw_edges, out_layer3h)
  rm(hcw_hcw_edges); gc()

  # ── Layer 3b: Admission lookup ─────────────────────────────────────────────
  hospital_hcw_lookup <- hcw_nodes %>%
    group_by(hospital_id) %>%
    summarise(hcw_list = list(person_id), .groups = "drop")

  admission_lookup <- non_hcw_nodes %>%
    select(person_id, hh_id, hospital_id) %>%
    left_join(hospital_hcw_lookup, by = "hospital_id")

  saveRDS(admission_lookup, out_layer3a)
  rm(admission_lookup, hospital_hcw_lookup); gc()

  n_saved <- n_saved + 1L

  elapsed <- proc.time()[["elapsed"]] - t_start
  rate    <- n_saved / max(elapsed, 0.1)
  eta     <- round((n_total - fi) / max(rate, 0.01) / 60, 1)
  message(sprintf("    Done | pop: %d | cells: %d | hcw: %d | ETA %.1f min",
                  nrow(pop_nodes), n_cells, sum(pop_nodes$is_hcw), eta))
}

# ==============================================================================
# [Done] Summary
# ==============================================================================

elapsed_total <- round(proc.time()[["elapsed"]] - t_start, 1)
message("\n=== Network Build Complete ===")
message(sprintf("  Territories : %d", n_total))
message(sprintf("  Saved       : %d", n_saved))
message(sprintf("  Skipped     : %d", n_skipped))
message(sprintf("  Total time  : %.1f sec (%.1f min)",
                elapsed_total, elapsed_total / 60))
