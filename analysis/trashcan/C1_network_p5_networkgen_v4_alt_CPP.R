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
library(Rcpp)
library(data.table)

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
# cell_mat 다시 만들기
cell_sizes   <- tabulate(pop_nodes$cell_id, nbins = n_cells)
max_per_cell <- max(cell_sizes)
cell_mat     <- matrix(0L, nrow = n_cells, ncol = max_per_cell)

cell_fill <- integer(n_cells)
for (idx in seq_len(nrow(pop_nodes))) {
  ci                          <- pop_nodes$cell_id[idx]
  cell_fill[ci]               <- cell_fill[ci] + 1L
  cell_mat[ci, cell_fill[ci]] <- idx
}
rm(cell_fill)
message(sprintf("  Cell matrix: %d x %d (%.1f MB)",
                n_cells, max_per_cell,
                object.size(cell_mat) / 1e6))

message("\n=== Section 6: Layer 2 — Community edges (Rcpp Version) ===")

# 1. C++ 코드를 R 함수로 컴파일 (처음 1회 실행 시 몇 초 소요됨)
# - C++의 동적 배열(std::vector)을 사용하여 R의 메모리 사전 할당/복사 문제 해결
# - Partial Fisher-Yates Shuffle 알고리즘으로 필요한 만큼만 초고속 샘플링 (O(k) 시간 복잡도)
# cppFunction 대신 sourceCpp를 사용합니다.

# 저장한 C++ 파일을 불러옵니다. 여기서 에러가 나는지 콘솔을 잘 확인하세요.
sourceCpp("function/build_edges.cpp")

# ------------------------------------------------------------------------------
# 3. 메인 함수 (자료형 에러 수정본)
# ------------------------------------------------------------------------------
create_comm_edges_rcpp <- function(node_df, cell_mat, cell_sizes,
                                   nn_cell_id, prem_mu, size = 0.1,
                                   chunk_size = 50000L, seed = 42) {

  on.exit(gc(), add = TRUE)

  N <- nrow(node_df)

  # ====================================================================
  # ★ 핵심 수정: C++이 에러를 뱉지 않도록 문자형을 정수형(Integer)으로 안전하게 변환
  # ====================================================================
  hh_ids   <- as.integer(as.factor(node_df$hh_id))  # 문자형 가구 ID -> 정수 변환
  cell_ids <- as.integer(node_df$cell_id)           # 셀 ID도 확실하게 정수 변환
  age_grps <- node_df$age_group

  n_chunks <- ceiling(N / chunk_size)
  chunk_results <- vector("list", n_chunks)
  t_start <- proc.time()[["elapsed"]]

  for (ci in seq_len(n_chunks)) {
    idx_start <- (ci - 1L) * chunk_size + 1L
    idx_end   <- min(ci * chunk_size, N)
    chunk_ids <- idx_start:idx_end
    n_chunk   <- length(chunk_ids)

    set.seed(seed + ci)
    chunk_mu  <- prem_mu[age_grps[chunk_ids]]

    # 난수 결과도 확실하게 정수로 변환해서 넘깁니다.
    n_draws   <- as.integer(rnbinom(n_chunk, size = size, mu = chunk_mu))

    active_mask <- n_draws > 0L
    if (!any(active_mask)) next

    active_ids   <- chunk_ids[active_mask]
    active_draws <- n_draws[active_mask]

    # C++ 함수 호출
    chunk_edges <- build_edges_cpp(
      active_ids = active_ids,
      active_draws = active_draws,
      hh_ids = hh_ids,
      cell_ids = cell_ids,
      cell_mat = cell_mat,
      cell_sizes = cell_sizes,
      nn_cell_id = nn_cell_id,
      seed = seed + ci * 999
    )

    if (nrow(chunk_edges) > 0) {
      chunk_results[[ci]] <- as.data.table(chunk_edges)
    }

    elapsed <- proc.time()[["elapsed"]] - t_start
    rate    <- idx_end / elapsed
    eta     <- round((N - idx_end) / rate / 60, 1)
    message(sprintf(
      "  [Layer 2 Rcpp] chunk %d/%d | %.1f%% | %.0f indiv/sec | ETA %.1f min",
      ci, n_chunks, 100 * idx_end / N, rate, eta
    ))
  }

  message("  Combining chunks and deduplicating...")
  out <- rbindlist(chunk_results[!sapply(chunk_results, is.null)])

  if (nrow(out) == 0L) {
    return(data.frame(from=integer(0), to=integer(0), from_cell=integer(0), to_cell=integer(0)))
  }

  out <- unique(out, by = c("from", "to"))
  setDF(out)

  out
}

# 3. 실행 (chunk_size를 1,000에서 100,000으로 대폭 올렸습니다!)
t0 <- proc.time()
comm_edges <- create_comm_edges_rcpp(
  pop_nodes, cell_mat, cell_sizes, nn_cell_id, prem_comm_mu,
  size       = comm_nb_size,
  chunk_size = 100000L,
  seed       = network_seed
)
elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
message(sprintf("  Community edges: %d (%.1f sec / %.1f min)",
                nrow(comm_edges), elapsed, elapsed / 60))

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
