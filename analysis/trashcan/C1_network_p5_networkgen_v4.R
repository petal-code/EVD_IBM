# ==============================================================================
# COD_IBM_network_build.R
# Purpose:
#   Build all contact network layers for Kinshasa EVD IBM
#   Layer 1: Household (full clique)
#   Layer 2: Community (Prem NB, 10km radius)
#   Layer 3: Healthcare (HCW-HCW clique + admission lookup)
#   Funeral : Layer 1 + Layer 2 at simulation time
#
# Key optimizations:
#   1. Cell-level distance sparse matrix (89k x 89k upper triangle)
#      → precomputed once via frNN, stored as dgCMatrix
#      → community edges store cell_id pairs, not float distances
#      → kernel applied instantly at simulation time via sparse lookup
#   2. frNN spatial index for O(log N) candidate lookup
#   3. Cell-level distance: ~89k unique cells vs 7.6M individuals
#
# Usage:
#   sample_fraction = 1.0  → full network
#   sample_fraction = 0.1  → 10% test run (cell-level sampling)
#
# Outputs (saved to output/network/):
#   {prefix}_nodes.rds            : node dataframe (person_id, hh_id, cell_id, ...)
#   {prefix}_cell_dist.rds        : sparse cell-cell distance matrix (dgCMatrix)
#   {prefix}_layer1_household.rds : household edges
#   {prefix}_layer2_community.rds : community edges (from_cell, to_cell stored)
#   {prefix}_layer3_hcw_edges.rds : HCW-HCW edges
#   {prefix}_layer3_admission.rds : admission lookup
# ==============================================================================

library(dplyr)
library(tidyr)
library(readxl)
library(parallel)
library(dbscan)
library(Matrix)
library(pbapply)  # Progress bar for parallel tasks

# ==============================================================================
# [Configuration] — Modify only this section
# ==============================================================================

province_name   <- "Kinshasa"
prem_dir        <- "data/Prem_contact"
prem_country    <- "Congo"
hf_path         <- "data/COD_GRID3_health_facilities_v8.csv"
synpop_path     <- "output/household/Kinshasa_synthetic_population.rds"
output_dir      <- "output/network"

# 1.0 = full network | 0.1 = 10% test run (cell-level sampling)
sample_fraction <- 1.0

# Network parameters
hcw_rate     <- 13.78 / 10000  # HCWs per total population (DRC)
comm_nb_size <- 0.1             # NB dispersion for community contacts
radius_km    <- 10              # Contact radius in km

n_cores      <- min(15L, max(1L, detectCores() - 1L))
network_seed <- 42

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

pop_raw  <- readRDS(synpop_path)
pop_full <- pop_raw %>%
  mutate(
    person_id = row_number(),
    age_group = pmin(floor(age / 5L) + 1L, 16L),  # Prem 16 age groups
    is_adult  = age >= 18L
  )

message(sprintf("  Full population : %d individuals", nrow(pop_full)))

# Subsample at cell level if sample_fraction < 1
if (is_test) {
  set.seed(network_seed)
  sampled_cells <- pop_full %>%
    distinct(cell_x, cell_y) %>%
    slice_sample(prop = sample_fraction)

  pop_nodes <- pop_full %>%
    semi_join(sampled_cells, by = c("cell_x", "cell_y")) %>%
    mutate(person_id = row_number())

  message(sprintf("  Sampled cells   : %d / %d",
                  nrow(sampled_cells),
                  n_distinct(pop_full$cell_x, pop_full$cell_y)))
  message(sprintf("  Sampled pop     : %d (%.1f%%)",
                  nrow(pop_nodes), 100 * nrow(pop_nodes) / nrow(pop_full)))
} else {
  pop_nodes <- pop_full
  message(sprintf("  Using full population: %d individuals", nrow(pop_nodes)))
}

rm(pop_raw, pop_full); gc()

# ── Assign cell_id (integer index for each unique 100m grid cell) ──────────
cell_lookup <- pop_nodes %>%
  distinct(cell_x, cell_y) %>%
  arrange(cell_x, cell_y) %>%
  mutate(cell_id = row_number())

pop_nodes <- pop_nodes %>%
  left_join(cell_lookup, by = c("cell_x", "cell_y"))

n_cells <- nrow(cell_lookup)
message(sprintf("  Adults (>= 18)  : %d", sum(pop_nodes$is_adult)))
message(sprintf("  Unique cells    : %d", n_cells))

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
# [Section 2] Build cell-cell sparse distance matrix
# Only store upper triangle (i < j) for cells within radius_km
# Used at simulation time for instant kernel weight lookup
# ==============================================================================

message("\n=== Section 2: Cell-cell distance matrix ===")
t0 <- proc.time()

# Cell center coordinates (lon/lat degrees)
cell_coords <- as.matrix(cell_lookup[, c("cell_x", "cell_y")])

# Pre-convert cell lon/lat to meters (planar approximation)
# Kinshasa is near equator (~lat -4) so error is negligible (<0.1%)
lat0_rad <- mean(cell_coords[, 2]) * pi / 180
cell_mx  <- cell_coords[, 1] * 111320 * cos(lat0_rad)  # longitude → meters
cell_my  <- cell_coords[, 2] * 110540                  # latitude  → meters
cell_m   <- cbind(cell_mx, cell_my)                    # meter coordinate matrix

radius_m <- radius_km * 1000  # 10km → 10,000m

# Build spatial index on meter coordinates — eps in meters, no degree conversion needed
message("  Building cell spatial index (frNN, meter coords)...")
nn_cells   <- frNN(cell_m, eps = radius_m)
nn_cell_id <- nn_cells$id
rm(nn_cells); gc()

# Build sparse upper-triangle distance matrix using Euclidean distance (meters)
message("  Computing cell-cell distances...")

# Use list accumulation — avoid repeated c() memory reallocation
from_l <- vector("list", n_cells)
to_l   <- vector("list", n_cells)
dist_l <- vector("list", n_cells)
k      <- 0L

for (ci in seq_len(n_cells)) {
  nbrs <- nn_cell_id[[ci]]
  nbrs <- nbrs[nbrs > ci]  # Upper triangle only
  if (length(nbrs) == 0) next

  # Euclidean distance in meters — no trig functions needed
  dx    <- cell_mx[nbrs] - cell_mx[ci]
  dy    <- cell_my[nbrs] - cell_my[ci]
  dists <- sqrt(dx^2 + dy^2) / 1000  # Convert to km for storage

  keep <- dists <= radius_km
  if (!any(keep)) next

  k           <- k + 1L
  from_l[[k]] <- rep(ci, sum(keep))
  to_l[[k]]   <- nbrs[keep]
  dist_l[[k]] <- dists[keep]

  if (ci %% 10000 == 0)
    message(sprintf("  Cell distances: %d / %d (%.1f%%)",
                    ci, n_cells, 100 * ci / n_cells))
}

# Collapse lists in one shot — much faster than repeated c()
from_c <- unlist(from_l[seq_len(k)], use.names = FALSE)
to_c   <- unlist(to_l[seq_len(k)],   use.names = FALSE)
dist_c <- unlist(dist_l[seq_len(k)], use.names = FALSE)
rm(from_l, to_l, dist_l); gc()

# Store as sparse matrix (dgCMatrix) — upper triangle only
cell_dist_mat <- sparseMatrix(
  i    = from_c, j = to_c, x = dist_c,
  dims = c(n_cells, n_cells),
  symmetric = FALSE  # Upper triangle only; access via get_cell_dist()
)

elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
message(sprintf("  Cell-cell pairs within %dkm: %d (%.1f sec)",
                radius_km, length(dist_c), elapsed))
message(sprintf("  Sparse matrix size: %.1f MB",
                object.size(cell_dist_mat) / 1e6))

out_path <- file.path(output_dir, sprintf("%s_cell_dist.rds", out_prefix))
saveRDS(list(dist_mat   = cell_dist_mat,
             cell_lookup = cell_lookup),
        out_path)
message(sprintf("  Saved: %s (%.1f MB)", basename(out_path),
                file.info(out_path)$size / 1e6))

rm(from_c, to_c, dist_c); gc()  # Keep nn_cell_id for Section 6

# ==============================================================================
# [Section 3] Assign HCWs and nearest hospitals
# ==============================================================================

message("\n=== Section 3: HCW assignment ===")

set.seed(network_seed)
n_hcw     <- round(nrow(pop_nodes) * hcw_rate)
adult_ids <- pop_nodes$person_id[pop_nodes$is_adult]
hcw_ids   <- sample(adult_ids, size = n_hcw, replace = FALSE)

pop_nodes <- pop_nodes %>%
  mutate(is_hcw = person_id %in% hcw_ids)

message(sprintf("  HCWs assigned: %d (%.4f%%)",
                sum(pop_nodes$is_hcw), 100 * mean(pop_nodes$is_hcw)))

# Convert hospital coords to meters (same planar projection as cell_mx/my)
hosp_mx <- hf_hospital$lon * 111320 * cos(lat0_rad)
hosp_my <- hf_hospital$lat * 110540

# Assign nearest hospital at CELL level using vectorized Euclidean distance
# 89k cells × 251 hospitals — fully vectorized, no apply loop
message("  Computing nearest hospital per cell...")
cell_hosp_idx <- max.col(
  -(outer(cell_mx, hosp_mx, "-")^2 +
      outer(cell_my, hosp_my, "-")^2)   # Negative squared distance → max.col finds minimum
)

cell_lookup <- cell_lookup %>%
  mutate(hospital_id = hf_hospital$OBJECTID[cell_hosp_idx])

# Merge hospital_id into pop_nodes via cell_id
pop_nodes <- pop_nodes %>%
  left_join(cell_lookup %>% select(cell_id, hospital_id),
            by = "cell_id")

hcw_nodes     <- pop_nodes %>% filter(is_hcw)
non_hcw_nodes <- pop_nodes %>% filter(!is_hcw)

message(sprintf("  Hospitals with HCWs: %d / %d",
                length(unique(hcw_nodes$hospital_id)), nrow(hf_hospital)))

# ==============================================================================
# [Section 4] Save node dataframe
# ==============================================================================

message("\n=== Section 4: Saving nodes ===")

nodes_save <- pop_nodes %>%
  select(person_id, hh_id, cell_id, hospital_id,
         age, age_group, is_hcw, is_adult,
         cell_x, cell_y, indiv_x, indiv_y)

out_path <- file.path(output_dir, sprintf("%s_nodes.rds", out_prefix))
saveRDS(nodes_save, out_path)
message(sprintf("  Saved: %s (%.1f MB)", basename(out_path),
                file.info(out_path)$size / 1e6))

# ==============================================================================
# [Section 5] Layer 1 — Household edges
# ==============================================================================

message("\n=== Section 5: Layer 1 — Household edges ===")
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
rm(hh_edges); gc()

# save.image(file = "my_workspace_2026_05_29.RData", compress = FALSE)
load("my_workspace_2026_05_29.RData")

# ==============================================================================
# [Section 6] Layer 2 — Community edges
# Each edge stores from_cell / to_cell instead of float distance
# Kernel weight at simulation time: cell_dist_mat[from_cell, to_cell]
#
# Candidate lookup strategy (no individual-level frNN needed):
#   cell frNN (already built in Section 2) → neighbour cell list
#   cell_members[[cell_id]] → individuals in that cell
#   → combine neighbour cell members as candidate pool
# ==============================================================================
# message("\n=== Section 6: Layer 2 — Community edges ===")
#
# # Build cell membership as integer matrix (cache-friendly vs list)
# # Each row = one cell, columns = person indices (0-padded)
# # Contiguous memory → much better cache hit rate than list of vectors
# message("  Building cell membership matrix...")
# cell_sizes   <- tabulate(pop_nodes$cell_id, nbins = n_cells)
# max_per_cell <- max(cell_sizes)
# cell_mat     <- matrix(0L, nrow = n_cells, ncol = max_per_cell)
#
# # Fill matrix: for each cell, store its person row-indices
# cell_fill <- integer(n_cells)
# for (idx in seq_len(nrow(pop_nodes))) {
#   ci                          <- pop_nodes$cell_id[idx]
#   cell_fill[ci]               <- cell_fill[ci] + 1L
#   cell_mat[ci, cell_fill[ci]] <- idx
# }
# rm(cell_fill); gc()
#
# message(sprintf("  Cell matrix: %d x %d (%.1f MB)",
#                 n_cells, max_per_cell,
#                 object.size(cell_mat) / 1e6))
#
# # Optimized single-core community edge builder
# #
# # Key improvements over the original:
# #   1. rnbinom() vectorized per chunk — eliminates per-individual RNG call overhead
# #   2. cand_buf pre-allocated once per chunk — no lapply/unlist allocation inside loop
# #   3. edge_buf pre-allocated matrix — replaces list-of-lists accumulation,
# #      removes bind_rows / unlist at flush time
# #   4. gc() removed from inner loop — avoids per-chunk stop-the-world GC pause;
# #      called once via on.exit() instead
# #   5. chunk_size enlarged (default 50k) — reduces loop overhead proportionally
# #   6. Deduplication via data.table::unique() — keyed hash, in-place, no extra copy
# #   7. Per-chunk disk flush — prevents GC slowdown from in-memory accumulation;
# #      only one chunk lives in RAM at a time; combined at the end via rbindlist
# create_comm_edges_chunked <- function(node_df, cell_mat, cell_sizes,
#                                       nn_cell_id, prem_mu, size = 0.1,
#                                       chunk_size = 50000L, seed = 42) {
#
#   # Ensure data.table is available for fast deduplication at the end
#   if (!requireNamespace("data.table", quietly = TRUE)) {
#     stop("data.table is required — install.packages('data.table')")
#   }
#   library(data.table)
#
#   # Single gc() at function exit — never inside the hot loop
#   on.exit(gc(), add = TRUE)
#
#   set.seed(seed)
#   N        <- nrow(node_df)
#   hh_ids   <- node_df$hh_id
#   cell_ids <- node_df$cell_id
#   age_grps <- node_df$age_group
#
#   n_chunks    <- ceiling(N / chunk_size)
#   tmp_dir     <- file.path(tempdir(), "comm_edges_chunks")
#   dir.create(tmp_dir, showWarnings = FALSE)
#   chunk_files <- character(n_chunks)   # paths of flushed chunk files; "" = empty chunk
#
#   # Clean up temp directory on exit regardless of success or error
#   on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
#
#   t_start <- proc.time()[["elapsed"]]
#
#   # Upper bound on candidate pool size for cand_buf pre-allocation.
#   # max neighbors per cell * max cells per neighborhood.
#   max_nbr_cells <- max(lengths(nn_cell_id))
#   cand_buf_size <- max_per_cell * max_nbr_cells
#   cand_buf      <- integer(cand_buf_size)   # reused every individual
#
#   for (ci in seq_len(n_chunks)) {
#     idx_start <- (ci - 1L) * chunk_size + 1L
#     idx_end   <- min(ci * chunk_size, N)
#     chunk_ids <- idx_start:idx_end
#     n_chunk   <- length(chunk_ids)
#
#     # --- Improvement 1: vectorized rnbinom for the entire chunk ---
#     # One call instead of n_chunk separate calls; avoids RNG per-iter overhead
#     chunk_mu  <- prem_mu[age_grps[chunk_ids]]
#     n_draws   <- rnbinom(n_chunk, size = size, mu = chunk_mu)
#
#     # Skip zero-draw individuals before entering the loop
#     active_pos <- which(n_draws > 0L)   # positions within chunk_ids
#     if (length(active_pos) == 0L) next
#
#     # --- Improvement 3: pre-allocated flat edge buffer ---
#     # Worst-case edges = active individuals * max possible degree.
#     # Multiplier 429 = ceil(qnbinom(0.99, size, mu=max(prem_mu)) * 1.2)
#     max_edges      <- length(active_pos) * 429L
#     edge_from      <- integer(max_edges)
#     edge_to        <- integer(max_edges)
#     edge_from_cell <- integer(max_edges)
#     edge_to_cell   <- integer(max_edges)
#     eptr <- 0L   # write pointer into edge buffer
#
#     for (ii in active_pos) {
#       i   <- chunk_ids[ii]
#       n_c <- n_draws[ii]
#
#       nbr_cells <- nn_cell_id[[cell_ids[i]]]
#       if (length(nbr_cells) == 0L) next
#
#       # --- Improvement 2: cand_buf pre-allocated, refilled in-place ---
#       # No lapply, no unlist, no intermediate list allocation
#       ptr <- 0L
#       for (cj in nbr_cells) {
#         sz <- cell_sizes[cj]
#         if (sz == 0L) next
#         cand_buf[ptr + seq_len(sz)] <- cell_mat[cj, seq_len(sz)]
#         ptr <- ptr + sz
#       }
#       if (ptr == 0L) next
#       cands <- cand_buf[seq_len(ptr)]
#
#       # Exclude self and same-household candidates
#       keep_mask <- (cands != i) & (hh_ids[cands] != hh_ids[i])
#       cands     <- cands[keep_mask]
#       if (length(cands) == 0L) next
#
#       n_s     <- min(n_c, length(cands))
#       sampled <- cands[sample(length(cands), n_s, replace = FALSE)]
#
#       # Keep i < j only — avoids within-layer duplicates at source
#       keep   <- sampled > i
#       n_keep <- sum(keep)
#       if (n_keep == 0L) next
#
#       # Grow buffer if needed (rare; only when multiplier underestimates)
#       if (eptr + n_keep > length(edge_from)) {
#         extra          <- max(n_keep, max_edges)   # double-or-more strategy
#         edge_from      <- c(edge_from,      integer(extra))
#         edge_to        <- c(edge_to,        integer(extra))
#         edge_from_cell <- c(edge_from_cell, integer(extra))
#         edge_to_cell   <- c(edge_to_cell,   integer(extra))
#       }
#
#       rows         <- eptr + seq_len(n_keep)
#       sampled_keep <- sampled[keep]
#       edge_from[rows]      <- i
#       edge_to[rows]        <- sampled_keep
#       edge_from_cell[rows] <- cell_ids[i]
#       edge_to_cell[rows]   <- cell_ids[sampled_keep]
#       eptr <- eptr + n_keep
#     }
#
#     # --- Improvement 7: flush chunk to disk immediately ---
#     # Only the current chunk's edges live in RAM; no accumulation across chunks.
#     # compress=FALSE → fast sequential write, read speed matches in-memory rbindlist
#     if (eptr > 0L) {
#       tmp_path    <- file.path(tmp_dir, sprintf("chunk_%05d.rds", ci))
#       saveRDS(
#         data.table(
#           from      = edge_from[seq_len(eptr)],
#           to        = edge_to[seq_len(eptr)],
#           from_cell = edge_from_cell[seq_len(eptr)],
#           to_cell   = edge_to_cell[seq_len(eptr)]
#         ),
#         tmp_path, compress = FALSE
#       )
#       chunk_files[ci] <- tmp_path
#     }
#
#     # Progress report
#     elapsed <- proc.time()[["elapsed"]] - t_start
#     rate    <- idx_end / elapsed
#     eta     <- round((N - idx_end) / rate / 60, 1)
#     message(sprintf(
#       "  [Layer 2] chunk %d / %d | %d / %d (%.1f%%) | %.0f indiv/sec | ETA %.1f min",
#       ci, n_chunks, idx_end, N, 100 * idx_end / N, rate, eta
#     ))
#
#     # No gc() here — let R schedule GC naturally
#   }
#
#   # Read all flushed chunks and combine — only non-empty paths
#   message("  Combining chunks...")
#   written <- chunk_files[nzchar(chunk_files)]
#   out     <- rbindlist(lapply(written, readRDS), use.names = FALSE)
#
#   if (nrow(out) == 0L) {
#     return(data.frame(
#       from = integer(0), to = integer(0),
#       from_cell = integer(0), to_cell = integer(0)
#     ))
#   }
#
#   # --- Improvement 6: data.table keyed deduplication ---
#   # Avoids integer overflow risk of hash multiplication; in-place, no extra copy
#   out <- unique(out, by = c("from", "to"))
#   setDF(out)   # return as plain data.frame to match downstream expectations
#
#   out
# }
#
# t0         <- proc.time()
# comm_edges <- create_comm_edges_chunked(
#   pop_nodes, cell_mat, cell_sizes, nn_cell_id, prem_comm_mu,
#   size       = comm_nb_size,
#   chunk_size = 1000L,
#   seed       = network_seed
# )
# elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
# message(sprintf("  Community edges: %d (%.1f sec / %.1f min)",
#                 nrow(comm_edges), elapsed, elapsed / 60))
#
# out_path <- file.path(output_dir, sprintf("%s_layer2_community.rds", out_prefix))
# saveRDS(comm_edges, out_path)
# message(sprintf("  Saved: %s (%.1f MB)", basename(out_path),
#                 file.info(out_path)$size / 1e6))
# rm(comm_edges, nn_cell_id, cell_mat, cell_sizes); gc()

message("\n=== Section 6: Layer 2 — Community edges ===")

# Build cell membership as integer matrix (cache-friendly vs list)
# Each row = one cell, columns = person indices (0-padded)
# Contiguous memory → much better cache hit rate than list of vectors
message("  Building cell membership matrix...")
cell_sizes   <- tabulate(pop_nodes$cell_id, nbins = n_cells)
max_per_cell <- max(cell_sizes)
cell_mat     <- matrix(0L, nrow = n_cells, ncol = max_per_cell)

# Fill matrix: for each cell, store its person row-indices
cell_fill <- integer(n_cells)
for (idx in seq_len(nrow(pop_nodes))) {
  ci                          <- pop_nodes$cell_id[idx]
  cell_fill[ci]               <- cell_fill[ci] + 1L
  cell_mat[ci, cell_fill[ci]] <- idx
}
rm(cell_fill); gc()

message(sprintf("  Cell matrix: %d x %d (%.1f MB)",
                n_cells, max_per_cell,
                object.size(cell_mat) / 1e6))

# Parallel community edge builder — PSOCK cluster (Windows-compatible)
#
# Key improvements over the original:
#   1. rnbinom() vectorized per chunk
#   2. cand_buf pre-allocated once per chunk
#   3. edge_buf pre-allocated flat matrix
#   4. gc() removed from inner loop
#   5. chunk_size enlarged
#   6. Deduplication via data.table::unique()
#   7. Per-chunk disk flush — only one chunk lives in RAM per worker
#   8. PSOCK parallel cluster — 10 workers; large objects exported once at startup
#   9. Batch-based progress monitoring — main process reports after every batch
create_comm_edges_chunked <- function(node_df, cell_mat, cell_sizes,
                                      nn_cell_id, prem_mu, size = 0.1,
                                      chunk_size = 50000L, n_cores = 10L,
                                      report_every = 50L,   # batches between progress lines
                                      seed = 42) {

  if (!requireNamespace("data.table", quietly = TRUE))
    stop("data.table is required")
  if (!requireNamespace("parallel",   quietly = TRUE))
    stop("parallel is required")

  library(data.table)
  library(parallel)

  N        <- nrow(node_df)
  hh_ids   <- node_df$hh_id
  cell_ids <- node_df$cell_id
  age_grps <- node_df$age_group

  n_chunks    <- ceiling(N / chunk_size)
  tmp_dir     <- file.path(tempdir(), "comm_edges_chunks")
  dir.create(tmp_dir, showWarnings = FALSE)
  on.exit({ unlink(tmp_dir, recursive = TRUE); gc() }, add = TRUE)

  max_nbr_cells <- max(lengths(nn_cell_id))
  cand_buf_size <- max_per_cell * max_nbr_cells

  # Worker function — all needed objects exported once via clusterExport
  process_chunk <- function(ci) {
    idx_start <- (ci - 1L) * chunk_size + 1L
    idx_end   <- min(ci * chunk_size, N)
    chunk_ids <- idx_start:idx_end
    n_chunk   <- length(chunk_ids)

    set.seed(seed + ci)

    chunk_mu   <- prem_mu[age_grps[chunk_ids]]
    n_draws    <- rnbinom(n_chunk, size = size, mu = chunk_mu)
    active_pos <- which(n_draws > 0L)
    if (length(active_pos) == 0L) return(list(path = NULL, edges = 0L))

    max_edges      <- length(active_pos) * 429L
    edge_from      <- integer(max_edges)
    edge_to        <- integer(max_edges)
    edge_from_cell <- integer(max_edges)
    edge_to_cell   <- integer(max_edges)
    eptr           <- 0L
    cand_buf       <- integer(cand_buf_size)

    for (ii in active_pos) {
      i   <- chunk_ids[ii]
      n_c <- n_draws[ii]

      nbr_cells <- nn_cell_id[[cell_ids[i]]]
      if (length(nbr_cells) == 0L) next

      ptr <- 0L
      for (cj in nbr_cells) {
        sz <- cell_sizes[cj]
        if (sz == 0L) next
        cand_buf[ptr + seq_len(sz)] <- cell_mat[cj, seq_len(sz)]
        ptr <- ptr + sz
      }
      if (ptr == 0L) next
      cands <- cand_buf[seq_len(ptr)]

      keep_mask <- (cands != i) & (hh_ids[cands] != hh_ids[i])
      cands     <- cands[keep_mask]
      if (length(cands) == 0L) next

      n_s     <- min(n_c, length(cands))
      sampled <- cands[sample(length(cands), n_s, replace = FALSE)]

      keep   <- sampled > i
      n_keep <- sum(keep)
      if (n_keep == 0L) next

      if (eptr + n_keep > length(edge_from)) {
        extra          <- max(n_keep, max_edges)
        edge_from      <- c(edge_from,      integer(extra))
        edge_to        <- c(edge_to,        integer(extra))
        edge_from_cell <- c(edge_from_cell, integer(extra))
        edge_to_cell   <- c(edge_to_cell,   integer(extra))
      }

      rows         <- eptr + seq_len(n_keep)
      sampled_keep <- sampled[keep]
      edge_from[rows]      <- i
      edge_to[rows]        <- sampled_keep
      edge_from_cell[rows] <- cell_ids[i]
      edge_to_cell[rows]   <- cell_ids[sampled_keep]
      eptr <- eptr + n_keep
    }

    if (eptr == 0L) return(list(path = NULL, edges = 0L))

    tmp_path <- file.path(tmp_dir, sprintf("chunk_%05d.rds", ci))
    saveRDS(
      data.table(
        from      = edge_from[seq_len(eptr)],
        to        = edge_to[seq_len(eptr)],
        from_cell = edge_from_cell[seq_len(eptr)],
        to_cell   = edge_to_cell[seq_len(eptr)]
      ),
      tmp_path, compress = FALSE
    )

    list(path = tmp_path, edges = eptr)
  }

  # Start PSOCK cluster
  message(sprintf("  Starting PSOCK cluster with %d workers...", n_cores))
  cl <- makeCluster(n_cores, type = "PSOCK")
  on.exit(stopCluster(cl), add = TRUE)

  message("  Exporting data to workers (may take a moment)...")
  t_export <- proc.time()[["elapsed"]]
  clusterExport(cl, varlist = c(
    "N", "hh_ids", "cell_ids", "age_grps",
    "chunk_size", "seed", "size", "prem_mu",
    "cell_mat", "cell_sizes", "nn_cell_id",
    "cand_buf_size", "max_per_cell", "tmp_dir"
  ), envir = environment())
  clusterEvalQ(cl, library(data.table))
  message(sprintf("  Export done (%.1f sec)", proc.time()[["elapsed"]] - t_export))

  # Split chunk indices into batches of n_cores
  # Each batch is submitted together; progress reported after each batch completes
  batch_size  <- n_cores
  batch_ids   <- split(seq_len(n_chunks), ceiling(seq_len(n_chunks) / batch_size))
  n_batches   <- length(batch_ids)

  chunk_files   <- character(n_chunks)
  total_edges   <- 0L
  chunks_done   <- 0L
  t_start       <- proc.time()[["elapsed"]]

  message(sprintf("  Processing %d chunks in %d batches of %d...",
                  n_chunks, n_batches, batch_size))

  for (bi in seq_len(n_batches)) {
    batch   <- batch_ids[[bi]]
    results <- clusterApplyLB(cl, batch, process_chunk)

    # Collect results from this batch
    for (ri in seq_along(results)) {
      r <- results[[ri]]
      if (!is.null(r$path)) {
        chunk_files[batch[ri]] <- r$path
        total_edges <- total_edges + r$edges
      }
    }
    chunks_done <- chunks_done + length(batch)

    # Progress report every report_every batches (or on the last batch)
    if (bi %% report_every == 0L || bi == n_batches) {
      elapsed <- proc.time()[["elapsed"]] - t_start
      indiv_done <- min(chunks_done * chunk_size, N)
      rate       <- indiv_done / elapsed
      eta        <- round((N - indiv_done) / rate / 60, 1)
      message(sprintf(
        "  [Layer 2] batch %d / %d | %d / %d indiv (%.1f%%) | %.0f indiv/sec | edges so far: %d | ETA %.1f min",
        bi, n_batches, indiv_done, N,
        100 * indiv_done / N, rate, total_edges, eta
      ))
    }
  }

  elapsed_total <- round(proc.time()[["elapsed"]] - t_start, 1)
  message(sprintf("  All chunks complete: %d edges | %.1f sec (%.1f min)",
                  total_edges, elapsed_total, elapsed_total / 60))

  written <- chunk_files[nzchar(chunk_files)]
  if (length(written) == 0L) {
    return(data.frame(from = integer(0), to = integer(0),
                      from_cell = integer(0), to_cell = integer(0)))
  }

  message(sprintf("  Combining %d chunk files...", length(written)))
  out <- rbindlist(lapply(written, readRDS), use.names = FALSE)

  message("  Deduplicating...")
  out <- unique(out, by = c("from", "to"))
  setDF(out)

  out
}

t0         <- proc.time()
comm_edges <- create_comm_edges_chunked(
  pop_nodes, cell_mat, cell_sizes, nn_cell_id, prem_comm_mu,
  size         = comm_nb_size,
  chunk_size   = 10000L,
  n_cores      = 10L,
  report_every = 1L,   # print progress every 10 batches
  seed         = network_seed
)
elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
message(sprintf("  Community edges: %d (%.1f sec / %.1f min)",
                nrow(comm_edges), elapsed, elapsed / 60))

out_path <- file.path(output_dir, sprintf("%s_layer2_community.rds", out_prefix))
saveRDS(comm_edges, out_path)
message(sprintf("  Saved: %s (%.1f MB)", basename(out_path),
                file.info(out_path)$size / 1e6))
rm(comm_edges, nn_cell_id, cell_mat, cell_sizes); gc()


# ===========================================================================
# [Section 7] Layer 3 — Healthcare edges
# ==============================================================================

message("\n=== Section 7: Layer 3 — Healthcare edges ===")

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
# [Done] Summary
# ==============================================================================

elapsed_total <- round((proc.time() - t_start)[["elapsed"]], 1)

message("\n=== Network Build Complete ===")
message(sprintf("  Province        : %s", province_name))
message(sprintf("  Sample fraction : %.0f%%", sample_fraction * 100))
message(sprintf("  Nodes           : %d", nrow(pop_nodes)))
message(sprintf("  Cells           : %d", n_cells))
message(sprintf("  Total time      : %.1f sec (%.1f min)",
                elapsed_total, elapsed_total / 60))
message(sprintf("\n  Funeral network : Layer 1 + Layer 2 (activated at simulation time)"))
message(sprintf("  Kernel lookup   : cell_dist_mat[from_cell, to_cell]"))
message("\n  Output files:")
for (f in list.files(output_dir, pattern = out_prefix)) {
  size_mb <- round(file.info(file.path(output_dir, f))$size / 1e6, 1)
  message(sprintf("    %-55s %6.1f MB", f, size_mb))
}
