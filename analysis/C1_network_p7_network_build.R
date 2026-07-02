# ==============================================================================
# C1_network_p7_network_build.R  (v2)
# Purpose:
#   For each simulation case (case1_1M, case2_Ituri, case3_Kivu):
#     1. Load synthetic population (from p4)
#     2. Build full cell-cell distance bucket matrix
#     3. Build cell-age member index
#     4. Build Layer 1: household edges (full clique)
#     5. Build Layer 2: community edges — single C++ pass.
#        v2 CHANGE: physical and non-physical contacts are now drawn as two
#        INDEPENDENT processes (not close-then-Bernoulli-thin):
#          - non-physical degree ~ Poisson(mu_nonphys_i)
#          - physical     degree ~ NegBin(mu_phys_i, size = phys_nb_size)
#        Each process also uses its OWN row-normalized age-partner
#        distribution (physical and non-physical contacts can have
#        different age patterns). Stratum (daily/weekly/monthly) is still
#        allocated via Multinomial(stratum_prob_mat) per participant age,
#        independent of physical status.
#        Collision handling: for a given participant, physical partners are
#        drawn first; if the non-physical pass later picks the same
#        partner, that pick is dropped ("physical wins" on collision).
#        Output: three separate edge files (daily/weekly/monthly),
#        each with columns: from, to, is_physical
#     6. Build Layer 3: healthcare edges
#        - HCW-HCW edges per hospital
#        - Non-HCW admission lookup
#     7. Save all network layers
#
# Input matrices from p6 (build_DRC_network_matrices.R):
#   prem_unique_community  : 16x16 unique contact matrix (sum of three strata)
#   stratum_prob_mat       : 16x3 stratum allocation probs [p_daily, p_weekly, p_monthly]
#   physical_unique_3wk    : list(daily/weekly/monthly)$community — used to
#                             build physical_unique_community
#   (nonphysical_unique_community is derived here as
#    prem_unique_community - physical_unique_community)
#
# Output per case (output/network/):
#   {tag}_nodes.rds
#   {tag}_cell_dist.rds
#   {tag}_layer1_household.rds
#   {tag}_layer2_daily.rds       — from, to, is_physical
#   {tag}_layer2_weekly.rds      — from, to, is_physical
#   {tag}_layer2_monthly.rds     — from, to, is_physical
#   {tag}_layer3_hcw_edges.rds
#   {tag}_layer3_admission.rds
# ==============================================================================

library(dplyr)
library(Rcpp)
library(sf)

# ==============================================================================
# [Configuration]
# ==============================================================================

matrices_path <- "output/MPMmat/DRC_network_input_matrices.rds"
hf_path       <- "data/COD_GRID3_health_facilities_v8.csv"
synpop_dir    <- "output/household"
kernel_path   <- "output/kernel/community_distance_kernel.rds"
output_dir    <- "output/network"
shp_path      <- "data/shpmap/gadm41_COD_2.shp"

hcw_rate      <- 13.78 / 10000  # HCWs per total population (DRC)
phys_nb_size  <- 0.25           # NB dispersion parameter for PHYSICAL contact draws
network_seed  <- 42L

# Cases to run — set to NULL to run all three
# Available: "case1_1M", "case2_Ituri", "case3_Kivu"
run_cases <- c("case1_1M")

# Bucket boundaries (km) — based on Meta Data for Good mobility bins
# Bucket 0: same cell (d=0.5km)
# Bucket 1: (0, 1.5]km
# Bucket 2: (1.5, 10.5]km
# Bucket 3: (10.5, 100.5]km
# Bucket 4: >100.5km
BUCKET_BREAKS <- c(1.5, 10.5, 100.5)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# [Section 1] Load shared data
# ==============================================================================

message("=== Section 1: Loading shared data ===")

# Load p6 output
mats <- readRDS(matrices_path)

# Total unique community contact matrix (daily+weekly+monthly combined)
prem_unique      <- mats$prem_unique_community   # 16x16
stratum_prob_mat <- mats$stratum_prob_mat        # 16x3 [p_daily, p_weekly, p_monthly]

# v2: physical / non-physical are split into two independent matrices,
# each row-normalized separately inside the C++ builder.
physical_unique_community <- mats$physical_unique_3wk$daily$community +
  mats$physical_unique_3wk$weekly$community +
  mats$physical_unique_3wk$monthly$community

nonphysical_unique_community <- prem_unique - physical_unique_community
# guard against tiny negative floating-point residue
nonphysical_unique_community[nonphysical_unique_community < 0] <- 0

# Also load household blended_ratio for Layer 1 (not used in C++ but kept for reference)
blended_ratio_home <- mats$blended_ratio$home

message(sprintf("  prem_unique rowSum range: %.2f - %.2f",
                min(rowSums(prem_unique)), max(rowSums(prem_unique))))
message(sprintf("  physical_unique_community rowSum range: %.2f - %.2f",
                min(rowSums(physical_unique_community)), max(rowSums(physical_unique_community))))
message(sprintf("  nonphysical_unique_community rowSum range: %.2f - %.2f",
                min(rowSums(nonphysical_unique_community)), max(rowSums(nonphysical_unique_community))))
message(sprintf("  recovery check (physical+nonphysical should equal prem_unique): max abs diff = %.2e",
                max(abs((physical_unique_community + nonphysical_unique_community) - prem_unique))))
message(sprintf("  stratum_prob_mat check (should sum to ~1): [%.4f, %.4f]",
                min(rowSums(stratum_prob_mat)), max(rowSums(stratum_prob_mat))))
message(sprintf("  phys_nb_size (physical draw dispersion): %.3f", phys_nb_size))

# Distance kernel parameters
kernel <- readRDS(kernel_path)
p_hat  <- kernel$p_hat
a1_hat <- kernel$a1_hat
a2_hat <- kernel$a2_hat
message(sprintf("  Kernel: p=%.4f a1=%.4f a2=%.4f", p_hat, a1_hat, a2_hat))

# Bucket kernel weights: integral of mixture exponential over each bucket
kernel_integral <- function(a, b) {
  p_hat     * (exp(-a1_hat * a) - exp(-a1_hat * b)) +
    (1-p_hat) * (exp(-a2_hat * a) - exp(-a2_hat * b))
}

bucket_weights <- c(
  kernel_integral(0,      0.5),    # Bucket 0: same cell (~0.5km)
  kernel_integral(0.5,    1.5),    # Bucket 1: 0.5-1.5km
  kernel_integral(1.5,   10.5),    # Bucket 2: 1.5-10.5km
  kernel_integral(10.5, 100.5),    # Bucket 3: 10.5-100.5km
  kernel_integral(100.5,  Inf)     # Bucket 4: 100.5km+
)
message(sprintf("  Bucket weights: %.4f %.4f %.4f %.4f %.4f",
                bucket_weights[1], bucket_weights[2], bucket_weights[3],
                bucket_weights[4], bucket_weights[5]))

# Hospital data
hf_all          <- read.csv(hf_path)
hf_hospital_all <- hf_all %>%
  filter(esstype %in% c("Hôpital", "Hôpital Général de Référence",
                        "Centre Hopitalier"),
         !is.na(lon), !is.na(lat))
message(sprintf("  Hospitals (DRC): %d", nrow(hf_hospital_all)))

# Level 2 shapefile — for hospital-territory assignment
cod2_sf   <- st_read(shp_path, quiet = TRUE)
hf_sf     <- st_as_sf(hf_hospital_all, coords = c("lon", "lat"), crs = 4326)
hf_joined <- st_join(hf_sf, cod2_sf[, c("NAME_1", "NAME_2")]) %>%
  st_drop_geometry()

# Compile C++ edge builder
message("  Compiling C++ edge builder...")
sourceCpp("function/build_edges.cpp")
message("  Done")

# ==============================================================================
# [Section 2] List simulation case files
# ==============================================================================

message("\n=== Section 2: Finding synthetic population files ===")

pop_files <- c(
  file.path(synpop_dir, "case1_1M_synthetic_population.rds"),
  file.path(synpop_dir, "case2_Ituri_synthetic_population.rds"),
  file.path(synpop_dir, "case3_Kivu_synthetic_population.rds")
)
pop_files <- pop_files[file.exists(pop_files)]

# Filter to selected cases if run_cases is specified
if (!is.null(run_cases)) {
  pop_files <- pop_files[grepl(paste(run_cases, collapse = "|"),
                               basename(pop_files))]
  message(sprintf("  Running cases: %s", paste(run_cases, collapse = ", ")))
}
message(sprintf("  Case files to process: %d", length(pop_files)))

# ==============================================================================
# [Helper] Build cell-cell distance bucket matrix
# ==============================================================================

build_bucket_matrix <- function(cell_tbl) {
  n_cells  <- nrow(cell_tbl)
  lat0_rad <- mean(cell_tbl$y) * pi / 180
  cell_mx  <- cell_tbl$x * 111320 * cos(lat0_rad)
  cell_my  <- cell_tbl$y * 110540

  bucket_mat <- matrix(0L, nrow = n_cells, ncol = n_cells)
  t0 <- proc.time()[["elapsed"]]

  for (ci in seq_len(n_cells)) {
    if (ci %% 1000L == 0L)
      message(sprintf("    bucket_mat: %d / %d", ci, n_cells))
    dx   <- cell_mx - cell_mx[ci]
    dy   <- cell_my - cell_my[ci]
    dist <- sqrt(dx^2 + dy^2) / 1000  # km

    b        <- integer(n_cells)
    b[dist <= 1.5]                   <- 1L
    b[dist > 1.5  & dist <= 10.5]   <- 2L
    b[dist > 10.5 & dist <= 100.5]  <- 3L
    b[dist > 100.5]                  <- 4L
    b[ci]    <- 0L  # same cell

    bucket_mat[ci, ] <- b
  }
  message(sprintf("    bucket_mat done: %.1f sec",
                  proc.time()[["elapsed"]] - t0))

  list(bucket_mat = bucket_mat, cell_mx = cell_mx, cell_my = cell_my)
}

# ==============================================================================
# [Helper] Build cell × age_group member index (flat list, length = n_cells*16)
# ==============================================================================

build_cell_age_members <- function(pop_nodes, n_cells) {
  result <- vector("list", n_cells * 16L)
  for (idx in seq_along(result)) result[[idx]] <- integer(0)

  for (row_i in seq_len(nrow(pop_nodes))) {
    cid      <- pop_nodes$cell_id[row_i]
    ag       <- pop_nodes$age_group[row_i]
    flat_idx <- (cid - 1L) * 16L + ag  # 1-based flat index
    result[[flat_idx]] <- c(result[[flat_idx]], pop_nodes$person_id[row_i])
  }
  result
}

# ==============================================================================
# [Helper] Build all community layers in one C++ pass
# v2: physical and non-physical are drawn as independent processes
# (Poisson vs NegBin, each with its own age-partner distribution); on
# collision within a participant's own draws, physical wins. Each edge is
# also allocated to daily/weekly/monthly via Multinomial(stratum_prob_mat).
# Returns list of three data.frames: daily, weekly, monthly
# each with columns: from, to, is_physical
# ==============================================================================

build_community_layers <- function(pop_nodes, cdist, seed) {
  t0 <- proc.time()[["elapsed"]]

  n_cells          <- nrow(pop_nodes %>% distinct(cell_id))
  cell_age_members <- build_cell_age_members(pop_nodes, n_cells)
  message(sprintf("    cell_age_members built: %.1f sec",
                  proc.time()[["elapsed"]] - t0))

  init_edge_builder(
    bucket_mat_r      = cdist$bucket_mat,
    cell_age_members  = cell_age_members,
    nonphys_unique    = nonphysical_unique_community,
    phys_unique       = physical_unique_community,
    stratum_prob_mat  = stratum_prob_mat,
    bucket_weights    = bucket_weights,
    phys_nb_size      = phys_nb_size
  )
  rm(cell_age_members); gc()

  t1    <- proc.time()[["elapsed"]]
  N     <- nrow(pop_nodes)

  edges <- as.data.frame(build_edges_cpp(
    active_ids = as.integer(seq_len(N)),
    cell_ids   = as.integer(pop_nodes$cell_id),
    hh_ids     = as.integer(pop_nodes$hh_id),
    age_groups = as.integer(pop_nodes$age_group),
    seed       = seed
  ))
  # stratum: 0=daily, 1=weekly, 2=monthly
  # is_physical is now determined by which pass (physical/non-physical)
  # generated the edge, with physical-wins collision handling already
  # applied per participant inside the C++ builder.

  message(sprintf("    C++ edges built: %d total (%.1f sec)",
                  nrow(edges), proc.time()[["elapsed"]] - t1))
  message(sprintf("    physical=%d | non-physical=%d",
                  sum(edges$is_physical == 1L), sum(edges$is_physical == 0L)))

  # Split by stratum and deduplicate within each.
  # NOTE: physical/non-physical collisions are already resolved inside the
  # C++ builder (physical wins) for a given participant's own draws. This
  # dedup step still guards against the pre-existing residual case where
  # two DIFFERENT participants independently sample each other within the
  # same stratum (kept with is_physical=max, i.e. physical wins there too).
  split_and_dedup <- function(s) {
    edges[edges$stratum == s, c("from","to","is_physical")] %>%
      group_by(from, to) %>%
      summarise(is_physical = max(is_physical), .groups = "drop")
  }

  list(
    daily   = split_and_dedup(0L),
    weekly  = split_and_dedup(1L),
    monthly = split_and_dedup(2L)
  )
}

# ==============================================================================
# [Section 3] Process each simulation case
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
  out_layer2d <- file.path(output_dir, sprintf("%s_layer2_daily.rds",     tag))
  out_layer2w <- file.path(output_dir, sprintf("%s_layer2_weekly.rds",    tag))
  out_layer2m <- file.path(output_dir, sprintf("%s_layer2_monthly.rds",   tag))
  out_layer3h <- file.path(output_dir, sprintf("%s_layer3_hcw_edges.rds", tag))
  out_layer3a <- file.path(output_dir, sprintf("%s_layer3_admission.rds", tag))

  all_out <- c(out_nodes, out_cdist, out_layer1,
               out_layer2d, out_layer2w, out_layer2m,
               out_layer3h, out_layer3a)

  if (all(file.exists(all_out))) {
    message(sprintf("  SKIP (all outputs exist): %s", tag))
    n_skipped <- n_skipped + 1L; next
  }

  message(sprintf("\n  [%d/%d] %s", fi, n_total, tag))

  # ── Load p4 output ──────────────────────────────────────────────────────────
  pop_raw  <- readRDS(pop_files[fi])
  pers     <- pop_raw$individuals   # person_id, hh_id, age
  hh_tbl   <- pop_raw$households    # hh_id, cell_id
  cell_tbl <- pop_raw$cells         # cell_id, x, y, cell_pop

  # ── Build flat node dataframe ────────────────────────────────────────────────
  pop_nodes <- pers %>%
    left_join(hh_tbl,  by = "hh_id") %>%
    left_join(cell_tbl %>% select(cell_id, x, y), by = "cell_id") %>%
    mutate(
      age_group = pmin(floor(age / 5L) + 1L, 16L),
      is_adult  = age >= 18L
    )

  n_cells <- nrow(cell_tbl)
  message(sprintf("    pop: %d | cells: %d", nrow(pop_nodes), n_cells))

  # ── Build full distance bucket matrix ────────────────────────────────────────
  cdist <- build_bucket_matrix(cell_tbl)

  # ── Assign HCWs (random sample of adults) ────────────────────────────────────
  set.seed(network_seed + fi)
  n_hcw     <- max(1L, round(nrow(pop_nodes) * hcw_rate))
  adult_ids <- pop_nodes$person_id[pop_nodes$is_adult]

  if (length(adult_ids) == 0) {
    message("    SKIP: no adults"); n_skipped <- n_skipped + 1L; next
  }

  hcw_ids   <- sample(adult_ids, size = min(n_hcw, length(adult_ids)),
                      replace = FALSE)
  pop_nodes <- pop_nodes %>% mutate(is_hcw = person_id %in% hcw_ids)
  message(sprintf("    HCWs assigned: %d (%.2f%%)",
                  sum(pop_nodes$is_hcw), mean(pop_nodes$is_hcw) * 100))

  # ── Assign nearest hospital ──────────────────────────────────────────────────
  terr_hospitals <- hf_joined %>%
    filter(paste0(gsub("[^A-Za-z0-9]", "_", NAME_1), "_",
                  gsub("[^A-Za-z0-9]", "_", NAME_2)) == tag)

  if (nrow(terr_hospitals) == 0) {
    x_range <- range(cell_tbl$x)
    y_range <- range(cell_tbl$y)
    terr_hospitals <- hf_hospital_all %>%
      filter(lon >= x_range[1], lon <= x_range[2],
             lat >= y_range[1], lat <= y_range[2])
  }
  if (nrow(terr_hospitals) == 0)
    terr_hospitals <- hf_hospital_all  # Last resort: nearest in DRC

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
  message(sprintf("    Hospitals assigned: %d unique",
                  n_distinct(pop_nodes$hospital_id)))

  # ── Save nodes + cell table ──────────────────────────────────────────────────
  saveRDS(pop_nodes %>%
            select(person_id, hh_id, cell_id, hospital_id,
                   age, age_group, is_hcw, is_adult, x, y),
          out_nodes)
  saveRDS(cell_tbl, out_cdist)
  message("    Saved: nodes + cell_dist")

  # ── Layer 1: Household edges (full clique per household) ─────────────────────
  # All household members connected; physical contact structure handled in
  # simulation via close_only_home and phys_only_home matrices (from p6)
  hh_edges <- lapply(
    split(pop_nodes$person_id, pop_nodes$hh_id),
    function(members) {
      if (length(members) < 2L) return(NULL)
      pairs <- combn(members, 2L)
      data.frame(from = pairs[1L,], to = pairs[2L,],
                 stringsAsFactors = FALSE)
    }) %>% Filter(Negate(is.null), .) %>% bind_rows()

  saveRDS(hh_edges, out_layer1)
  message(sprintf("    Layer 1 (household): %d edges", nrow(hh_edges)))
  rm(hh_edges); gc()

  # ── Layer 2: Community edges — single C++ pass, split into daily/weekly/monthly
  # v2: physical (NegBin, k=0.25) and non-physical (Poisson) sampled
  # independently per participant; physical wins on within-person collision.

  if (all(file.exists(c(out_layer2d, out_layer2w, out_layer2m)))) {
    message("    Layer 2: all strata exist, skip")
  } else {
    layer2_all <- build_community_layers(
      pop_nodes = pop_nodes,
      cdist     = cdist,
      seed      = network_seed + fi * 10L
    )
    saveRDS(layer2_all$daily,   out_layer2d)
    saveRDS(layer2_all$weekly,  out_layer2w)
    saveRDS(layer2_all$monthly, out_layer2m)
    message(sprintf("    Layer 2 saved: daily=%d | weekly=%d | monthly=%d edges",
                    nrow(layer2_all$daily), nrow(layer2_all$weekly),
                    nrow(layer2_all$monthly)))
    rm(layer2_all); gc()
  }

  # ── Layer 3a: HCW-HCW edges (within hospital) ────────────────────────────────
  hcw_nodes <- pop_nodes %>% filter(is_hcw)

  hcw_hcw_edges <- hcw_nodes %>%
    group_by(hospital_id) %>%
    group_map(~ {
      members <- .x$person_id
      if (length(members) < 2L) return(NULL)
      pairs <- combn(members, 2L)
      data.frame(from        = pairs[1L,],
                 to          = pairs[2L,],
                 hospital_id = .y$hospital_id,
                 stringsAsFactors = FALSE)
    }) %>% bind_rows()

  saveRDS(hcw_hcw_edges, out_layer3h)
  message(sprintf("    Layer 3 HCW-HCW: %d edges", nrow(hcw_hcw_edges)))
  rm(hcw_hcw_edges); gc()

  # ── Layer 3b: Admission lookup ────────────────────────────────────────────────
  # For each non-HCW: hospital_id + list of HCWs at that hospital
  # Used at simulation time to draw patient-HCW contacts upon admission
  hospital_hcw_lookup <- hcw_nodes %>%
    group_by(hospital_id) %>%
    summarise(hcw_list = list(person_id), .groups = "drop")

  admission_lookup <- pop_nodes %>%
    filter(!is_hcw) %>%
    select(person_id, hh_id, hospital_id) %>%
    left_join(hospital_hcw_lookup, by = "hospital_id")

  saveRDS(admission_lookup, out_layer3a)
  message(sprintf("    Layer 3 admission lookup: %d persons", nrow(admission_lookup)))
  rm(admission_lookup, hospital_hcw_lookup, hcw_nodes); gc()

  # ── Done ─────────────────────────────────────────────────────────────────────
  n_saved <- n_saved + 1L
  elapsed <- proc.time()[["elapsed"]] - t_start
  message(sprintf("    [%s] Complete | elapsed %.1f min",
                  tag, elapsed / 60))
}

# ==============================================================================
# [Done]
# ==============================================================================

elapsed_total <- round(proc.time()[["elapsed"]] - t_start, 1)
message("\n=== Network Build Complete ===")
message(sprintf("  Saved   : %d / %d", n_saved,   n_total))
message(sprintf("  Skipped : %d",       n_skipped))
message(sprintf("  Time    : %.1f sec (%.1f min)",
                elapsed_total, elapsed_total / 60))
