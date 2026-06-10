# ==============================================================================
# C1_network_p6_network_build.R
# Purpose:
#   For each Level 2 territory (+ metro regions):
#     1. Load synthetic population (from p4)
#     2. Build full cell-cell distance matrix → bucket matrix (uint8)
#     3. Build cell-age member index
#     4. Build Layer 1: household edges (full clique)
#     5. Build Layer 2: community edges (NB → Multinomial → kernel-weighted)
#     6. Build Layer 3: healthcare edges
#     7. Save all network layers
# ==============================================================================

library(dplyr)
library(readxl)
library(Rcpp)
library(data.table)
library(sf)

# ==============================================================================
# [Configuration]
# ==============================================================================

prem_dir     <- "data/Prem_contact"
prem_country <- "Congo"
hf_path      <- "data/COD_GRID3_health_facilities_v8.csv"
synpop_dir   <- "output/household"
kernel_path  <- "output/kernel/community_distance_kernel.rds"
output_dir   <- "output/network"
shp_path     <- "data/shpmap/gadm41_COD_2.shp"

hcw_rate     <- 13.78 / 10000  # HCWs per total population (DRC)
comm_nb_size <- 0.1             # NB dispersion
network_seed <- 42
chunk_size   <- 100000L

# Bucket boundaries (km) — based on Meta Data for Good mobility bins
# Bucket 0: same cell (d=0.5), 1: (0,1], 2: (1,10], 3: (10,100], 4: (100,1000]
BUCKET_BREAKS <- c(1.0, 10.0, 100.0, 1000.0)

# ==============================================================================
# [Batch control] — Change ONLY these numbers
# ==============================================================================
batch_id  <- 1L
n_batches <- 1L
n_test    <- NULL  # Set to e.g. 3L for test run

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# [Section 1] Load shared data
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
prem_comm    <- prem_mats$community  # 16 × 16
prem_comm_mu <- rowSums(prem_comm)
message(sprintf("  Prem community mu range: %.2f - %.2f",
                min(prem_comm_mu), max(prem_comm_mu)))

# Distance kernel parameters
kernel <- readRDS(kernel_path)
p_hat  <- kernel$p_hat
a1_hat <- kernel$a1_hat
a2_hat <- kernel$a2_hat
message(sprintf("  Kernel: p=%.4f a1=%.4f a2=%.4f", p_hat, a1_hat, a2_hat))

# Bucket kernel weights
# Bucket 0: same cell d=0.5km
# Bucket 1: mid=0.5km (0,1]  — same midpoint, different from self
# Bucket 2: mid=5.5km (1,10]
# Bucket 3: mid=55km  (10,100]
# Bucket 4: survival at 100km = p*exp(-a1*100) + (1-p)*exp(-a2*100)
kernel_fn <- function(d)
  p_hat * a1_hat * exp(-a1_hat * d) + (1-p_hat) * a2_hat * exp(-a2_hat * d)

bucket_weights <- c(
  kernel_fn(0.5),   # Bucket 0: same cell
  kernel_fn(0.5),   # Bucket 1: 0-1km (mid 0.5km)
  kernel_fn(5.5),   # Bucket 2: 1-10km
  kernel_fn(55.0),  # Bucket 3: 10-100km
  p_hat * exp(-a1_hat * 100) + (1-p_hat) * exp(-a2_hat * 100)  # Bucket 4: survival
)
message(sprintf("  Bucket weights: %.4f %.4f %.4f %.4f %.4f",
                bucket_weights[1], bucket_weights[2], bucket_weights[3],
                bucket_weights[4], bucket_weights[5]))

# Hospital data
hf_all <- read.csv(hf_path)
hf_hospital_all <- hf_all %>%
  filter(esstype %in% c("Hôpital", "Hôpital Général de Référence",
                        "Centre Hopitalier"),
         !is.na(lon), !is.na(lat))
message(sprintf("  Hospitals (DRC): %d", nrow(hf_hospital_all)))

# Level 2 shapefile
cod2_sf   <- st_read(shp_path, quiet = TRUE)
hf_sf     <- st_as_sf(hf_hospital_all, coords = c("lon", "lat"), crs = 4326)
hf_joined <- st_join(hf_sf, cod2_sf[, c("NAME_1", "NAME_2")]) %>%
  st_drop_geometry()

# Compile C++ edge builder
message("  Compiling C++ edge builder...")
sourceCpp("function/build_edges.cpp")
message("  Compiled successfully")

# ==============================================================================
# [Section 2] List and batch files
# ==============================================================================

message("\n=== Section 2: Finding synthetic population files ===")

pop_files_all <- list.files(synpop_dir,
                            pattern = "_synthetic_population\\.rds$",
                            full.names = TRUE)
message(sprintf("  Total: %d", length(pop_files_all)))

# batch_idx <- seq(batch_id, length(pop_files_all), by = n_batches)
pop_files <- pop_files_all[batch_idx]
pop_files <- pop_files_all[c(162, 163)]

if (!is.null(n_test)) {
  pop_files <- head(pop_files, n_test)
  message(sprintf("  TEST MODE: %d territories", n_test))
}
message(sprintf("  Batch %d/%d: %d territories", batch_id, n_batches, length(pop_files)))

# ==============================================================================
# [Helper] Build full distance bucket matrix
# ==============================================================================

build_bucket_matrix <- function(cell_tbl) {
  n_cells  <- nrow(cell_tbl)
  lat0_rad <- mean(cell_tbl$y) * pi / 180
  cell_mx  <- cell_tbl$x * 111320 * cos(lat0_rad)
  cell_my  <- cell_tbl$y * 110540

  # Full n_cells × n_cells bucket matrix (integer 0-4)
  bucket_mat <- matrix(0L, nrow = n_cells, ncol = n_cells)

  t0 <- proc.time()[["elapsed"]]
  for (ci in seq_len(n_cells)) {
    if (ci %% 1000L == 0L)
      message(sprintf("    bucket_mat: %d / %d cells", ci, n_cells))
    dx   <- cell_mx - cell_mx[ci]
    dy   <- cell_my - cell_my[ci]
    dist <- sqrt(dx^2 + dy^2) / 1000  # km

    # Assign buckets
    b <- integer(n_cells)
    b[dist <= 1.0]                        <- 1L
    b[dist > 1.0  & dist <= 10.0]        <- 2L
    b[dist > 10.0 & dist <= 100.0]       <- 3L
    b[dist > 100.0]                       <- 4L
    b[ci] <- 0L  # Same cell: bucket 0 (d=0.5km)

    bucket_mat[ci, ] <- b
  }
  message(sprintf("    bucket_mat done: %.1f sec", proc.time()[["elapsed"]] - t0))

  list(bucket_mat = bucket_mat, cell_mx = cell_mx, cell_my = cell_my)
}

# ==============================================================================
# [Helper] Build cell × age_group member index
# ==============================================================================

build_cell_age_members <- function(pop_nodes, n_cells) {
  # Returns flat list of length n_cells * 16
  # Index: (cell_id - 1) * 16 + (age_group - 1)
  result <- vector("list", n_cells * 16L)
  for (idx in seq_along(result)) result[[idx]] <- integer(0)

  for (row_i in seq_len(nrow(pop_nodes))) {
    cid <- pop_nodes$cell_id[row_i]
    ag  <- pop_nodes$age_group[row_i]
    flat_idx <- (cid - 1L) * 16L + ag  # 1-based flat index
    result[[flat_idx]] <- c(result[[flat_idx]], pop_nodes$person_id[row_i])
  }
  result
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

  fname <- basename(pop_files[fi])
  tag   <- sub("_synthetic_population\\.rds$", "", fname)

  out_nodes   <- file.path(output_dir, sprintf("%s_nodes.rds",            tag))
  out_cdist   <- file.path(output_dir, sprintf("%s_cell_dist.rds",        tag))
  out_layer1  <- file.path(output_dir, sprintf("%s_layer1_household.rds", tag))
  out_layer2  <- file.path(output_dir, sprintf("%s_layer2_community.rds", tag))
  out_layer3h <- file.path(output_dir, sprintf("%s_layer3_hcw_edges.rds", tag))
  out_layer3a <- file.path(output_dir, sprintf("%s_layer3_admission.rds", tag))

  if (all(file.exists(c(out_nodes, out_cdist, out_layer1,
                        out_layer2, out_layer3h, out_layer3a)))) {
    n_skipped <- n_skipped + 1L; next
  }

  message(sprintf("\n  [%d/%d] %s", fi, n_total, tag))

  # ── Load p4 output ─────────────────────────────────────────────────────────
  pop_raw  <- readRDS(pop_files[fi])
  pers     <- pop_raw$individuals   # person_id, hh_id, age
  hh_tbl   <- pop_raw$households    # hh_id, cell_id
  cell_tbl <- pop_raw$cells         # cell_id, x, y, cell_pop

  # ── Build flat node dataframe ──────────────────────────────────────────────
  pop_nodes <- pers %>%
    left_join(hh_tbl,  by = "hh_id") %>%
    left_join(cell_tbl %>% select(cell_id, x, y), by = "cell_id") %>%
    mutate(
      age_group = pmin(floor(age / 5L) + 1L, 16L),
      is_adult  = age >= 18L
    )

  n_cells <- nrow(cell_tbl)
  message(sprintf("    pop: %d | cells: %d", nrow(pop_nodes), n_cells))

  # ── Build full distance bucket matrix ──────────────────────────────────────
  cdist <- build_bucket_matrix(cell_tbl)

  # ── Assign HCWs ────────────────────────────────────────────────────────────
  set.seed(network_seed + fi)
  n_hcw     <- max(1L, round(nrow(pop_nodes) * hcw_rate))
  adult_ids <- pop_nodes$person_id[pop_nodes$is_adult]

  if (length(adult_ids) == 0) {
    message(sprintf("    SKIP: no adults")); n_skipped <- n_skipped + 1L; next
  }

  hcw_ids   <- sample(adult_ids, size = min(n_hcw, length(adult_ids)),
                      replace = FALSE)
  pop_nodes <- pop_nodes %>% mutate(is_hcw = person_id %in% hcw_ids)

  # ── Assign nearest hospital ─────────────────────────────────────────────────
  terr_hospitals <- hf_joined %>%
    filter(paste0(gsub("[^A-Za-z0-9]", "_", NAME_1), "_",
                  gsub("[^A-Za-z0-9]", "_", NAME_2)) == tag)

  if (nrow(terr_hospitals) == 0)
    terr_hospitals <- hf_hospital_all

  lat0_rad      <- mean(cell_tbl$y) * pi / 180
  hosp_mx       <- terr_hospitals$lon * 111320 * cos(lat0_rad)
  hosp_my       <- terr_hospitals$lat * 110540
  cell_hosp_idx <- max.col(
    -(outer(cdist$cell_mx, hosp_mx, "-")^2 +
        outer(cdist$cell_my, hosp_my, "-")^2)
  )

  cell_tbl  <- cell_tbl %>%
    mutate(hospital_id = terr_hospitals$OBJECTID[cell_hosp_idx])
  pop_nodes <- pop_nodes %>%
    left_join(cell_tbl %>% select(cell_id, hospital_id), by = "cell_id")

  hcw_nodes     <- pop_nodes %>% filter(is_hcw)
  non_hcw_nodes <- pop_nodes %>% filter(!is_hcw)

  # ── Save nodes + cell table ────────────────────────────────────────────────
  saveRDS(pop_nodes %>%
            select(person_id, hh_id, cell_id, hospital_id,
                   age, age_group, is_hcw, is_adult, x, y),
          out_nodes)
  saveRDS(cell_tbl, out_cdist)

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

  # ── Layer 2: Community edges ───────────────────────────────────────────────
  # Build cell-age member index (flat list, length = n_cells * 16)
  t0 <- proc.time()[["elapsed"]]
  cell_age_members <- build_cell_age_members(pop_nodes, n_cells)
  message(sprintf("    cell_age_members: %.1f sec", proc.time()[["elapsed"]] - t0))

  # Initialize C++ global state once per territory
  init_edge_builder(
    bucket_mat_r     = cdist$bucket_mat,
    cell_age_members = cell_age_members,
    prem_matrix      = prem_comm,
    bucket_weights   = bucket_weights,
    nb_size          = comm_nb_size
  )
  rm(cell_age_members); gc()

  # Process in chunks
  N            <- nrow(pop_nodes)
  n_chunks_run <- ceiling(N / chunk_size)
  chunk_results <- vector("list", n_chunks_run)

  t_layer2 <- proc.time()[["elapsed"]]
  for (ci in seq_len(n_chunks_run)) {
    idx_start <- (ci - 1L) * chunk_size + 1L
    idx_end   <- min(ci * chunk_size, N)
    chunk_ids <- idx_start:idx_end

    set.seed(network_seed + ci)
    chunk_edges <- build_edges_cpp(
      active_ids  = as.integer(chunk_ids),
      cell_ids    = as.integer(pop_nodes$cell_id),
      hh_ids      = as.integer(pop_nodes$hh_id),
      age_groups  = as.integer(pop_nodes$age_group),
      seed        = network_seed + ci * 999L
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
  rm(comm_edges, chunk_results, cdist); gc()

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
  message(sprintf("    Done | hcw: %d | ETA %.1f min",
                  sum(pop_nodes$is_hcw), eta))
}

# ==============================================================================
# [Done]
# ==============================================================================

elapsed_total <- round(proc.time()[["elapsed"]] - t_start, 1)
message("\n=== Network Build Complete ===")
message(sprintf("  Batch   : %d / %d", batch_id, n_batches))
message(sprintf("  Saved   : %d / %d", n_saved, n_total))
message(sprintf("  Skipped : %d", n_skipped))
message(sprintf("  Time    : %.1f sec (%.1f min)",
                elapsed_total, elapsed_total / 60))
