library(readxl)
library(dplyr)

prem_dir     <- "data/Prem_contact"
prem_country <- "Congo"

# ── Fixed loader: skip header row, force numeric ───────────────
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

    # First row is label row (X1, X2, ...) — drop it
    raw <- raw[-1, ]

    # Convert to numeric matrix (16 x 16)
    mat <- matrix(as.numeric(as.matrix(raw)), nrow = 16, ncol = 16)

    # Symmetrize
    (mat + t(mat)) / 2
  }

  mats <- lapply(setNames(names(file_map), names(file_map)), load_one)

  # Community matrix = work + school + other
  comm_mat <- matrix(0, 16, 16)
  for (s in c("work", "school", "other")) {
    if (!is.null(mats[[s]])) comm_mat <- comm_mat + mats[[s]]
  }
  mats$community <- comm_mat
  mats
}

prem_mats <- load_prem_matrices(prem_dir, prem_country)

# ── Verify ─────────────────────────────────────────────────────
cat("Prem matrices loaded:\n")
for (nm in names(prem_mats)) {
  cat(sprintf("  %-12s: %dx%d | sum=%.2f\n",
              nm, nrow(prem_mats[[nm]]), ncol(prem_mats[[nm]]),
              sum(prem_mats[[nm]])))
}



library(dplyr)

# ── Step 1: Prepare node dataframe ────────────────────────────
# Add person_id and Prem age group (16 groups: 0-4, 5-9, ..., 75+)
pop_nodes <- pop %>%
  mutate(
    person_id = row_number(),
    # Prem age groups: 1=0-4, 2=5-9, ..., 16=75+
    age_group = pmin(floor(age / 5) + 1L, 16L),
    # Adult flag for HCW assignment (age >= 18)
    is_adult  = age >= 18
  )

cat(sprintf("Total individuals : %d\n", nrow(pop_nodes)))
cat(sprintf("Adults (>= 18)    : %d\n", sum(pop_nodes$is_adult)))
cat(sprintf("Age group range   : %d - %d\n",
            min(pop_nodes$age_group), max(pop_nodes$age_group)))
cat("\nAge group distribution:\n")
print(table(pop_nodes$age_group))


# ── Step 2: Assign HCWs from adult population ─────────────────
set.seed(42)

# HCW rate: 13.78 per 10,000 total population
hcw_rate    <- 13.78 / 10000
n_hcw       <- round(nrow(pop_nodes) * hcw_rate)
cat(sprintf("Target HCWs: %d (%.4f%% of population)\n",
            n_hcw, 100 * hcw_rate))

# Sample HCWs from adults only
adult_ids   <- pop_nodes$person_id[pop_nodes$is_adult]
hcw_ids     <- sample(adult_ids, size = n_hcw, replace = FALSE)

# Flag HCWs in node dataframe
pop_nodes <- pop_nodes %>%
  mutate(is_hcw = person_id %in% hcw_ids)

cat(sprintf("HCWs assigned     : %d\n", sum(pop_nodes$is_hcw)))
cat(sprintf("HCW age range     : %d - %d\n",
            min(pop_nodes$age[pop_nodes$is_hcw]),
            max(pop_nodes$age[pop_nodes$is_hcw])))

# ── Step 3: Assign each HCW to nearest hospital ───────────────
# Haversine distance (vectorized)
haversine_km <- function(lon1, lat1, lon2, lat2) {
  R   <- 6371
  phi <- (lat2 - lat1) * pi / 180
  lam <- (lon2 - lon1) * pi / 180
  a   <- sin(phi/2)^2 + cos(lat1*pi/180) * cos(lat2*pi/180) * sin(lam/2)^2
  2 * R * asin(sqrt(a))
}

hcw_nodes <- pop_nodes %>% filter(is_hcw)

# For each HCW find nearest hospital
cat("Assigning HCWs to nearest hospital...\n")
hosp_coords <- as.matrix(hf_hospital[, c("lon", "lat")])

nearest_hospital <- apply(as.matrix(hcw_nodes[, c("indiv_x", "indiv_y")]),
                          1, function(coord) {
                            dists <- haversine_km(coord[1], coord[2],
                                                  hosp_coords[, 1], hosp_coords[, 2])
                            which.min(dists)
                          })

hcw_nodes$hospital_idx  <- nearest_hospital
hcw_nodes$hospital_id   <- hf_hospital$OBJECTID[nearest_hospital]
hcw_nodes$hospital_dist <- apply(
  cbind(hcw_nodes$indiv_x, hcw_nodes$indiv_y,
        hosp_coords[nearest_hospital, 1],
        hosp_coords[nearest_hospital, 2]),
  1, function(r) haversine_km(r[1], r[2], r[3], r[4])
)

cat(sprintf("HCWs assigned to hospitals\n"))
cat(sprintf("  Mean dist to hospital : %.2f km\n", mean(hcw_nodes$hospital_dist)))
cat(sprintf("  Max dist to hospital  : %.2f km\n", max(hcw_nodes$hospital_dist)))
cat(sprintf("  Hospitals with HCWs   : %d / %d\n",
            length(unique(hcw_nodes$hospital_id)), nrow(hf_hospital)))
cat("HCW per hospital (top 10):\n")
print(sort(table(hcw_nodes$hospital_id), decreasing = TRUE)[1:10])




library(dplyr)
library(igraph)
library(parallel)

# ── Configuration ──────────────────────────────────────────────
set.seed(42)

# Use test subset (1000 cells) for now
test_cells <- pop_nodes %>%
  distinct(cell_x, cell_y) %>%
  slice_sample(n = 1000)

nodes <- pop_nodes %>%
  semi_join(test_cells, by = c("cell_x", "cell_y"))

cat(sprintf("Test nodes: %d individuals | %d households\n",
            nrow(nodes), length(unique(nodes$hh_id))))

# ── Layer 1: Household edges ───────────────────────────────────
create_hh_edges <- function(node_df) {
  edge_list <- lapply(split(node_df$person_id, node_df$hh_id), function(members) {
    if (length(members) < 2) return(NULL)
    pairs <- combn(members, 2)
    data.frame(from  = pairs[1, ],
               to    = pairs[2, ],
               layer = "household",
               stringsAsFactors = FALSE)
  })
  bind_rows(Filter(Negate(is.null), edge_list))
}

cat("Building Layer 1: household edges...\n")
hh_edges <- create_hh_edges(nodes)
cat(sprintf("  Household edges: %d\n", nrow(hh_edges)))

# ── Layer 2: Community edges (NB, 10km radius) ────────────────
# Distance stored for later kernel weighting — NOT applied here
create_comm_edges_parallel <- function(node_df, mu = 10, size = 0.1,
                                       radius_km = 10, n_cores = 15L, seed = 42) {
  N      <- nrow(node_df)
  coords <- as.matrix(node_df[, c("indiv_x", "indiv_y")])
  hh_ids <- node_df$hh_id
  chunks <- split(seq_len(N), ceiling(seq_len(N) / ceiling(N / n_cores)))

  worker <- function(chunk_ids, coords, hh_ids, mu, size, radius_km, seed) {
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
      n_c <- rnbinom(1, size = size, mu = mu)
      if (n_c == 0) next
      lon_i   <- coords[i, 1]; lat_i <- coords[i, 2]
      lon_rad <- radius_km / (111.0 * cos(lat_i * pi / 180))

      cands <- which(abs(coords[, 2] - lat_i) <= lat_rad &
                       abs(coords[, 1] - lon_i) <= lon_rad &
                       hh_ids != hh_ids[i] & seq_len(nrow(coords)) != i)
      if (length(cands) == 0) next

      dists <- hav_km(lon_i, lat_i, coords[cands, 1], coords[cands, 2])
      cands <- cands[dists <= radius_km]
      dists <- dists[dists <= radius_km]
      if (length(cands) == 0) next

      n_s     <- min(n_c, length(cands))
      idx_s   <- sample(length(cands), n_s, replace = FALSE)
      sampled <- cands[idx_s]
      d_s     <- dists[idx_s]

      # Keep i < j only to avoid duplicate edges
      keep <- sampled > i
      if (any(keep)) {
        from_v <- c(from_v, rep(i, sum(keep)))
        to_v   <- c(to_v, sampled[keep])
        dist_v <- c(dist_v, d_s[keep])  # Store distance for later
      }
    }
    if (length(from_v) == 0) return(NULL)
    data.frame(from        = from_v,
               to          = to_v,
               distance_km = dist_v,  # Distance stored, weight applied later
               layer       = "community",
               stringsAsFactors = FALSE)
  }

  cl <- makeCluster(n_cores, type = "PSOCK")
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, varlist = c("coords", "hh_ids", "mu", "size", "radius_km"),
                envir = environment())
  chunk_args <- lapply(seq_along(chunks), function(ci)
    list(chunk_ids = chunks[[ci]], seed = seed + ci))

  results <- parLapply(cl, chunk_args, function(args)
    worker(args$chunk_ids, coords, hh_ids, mu, size, radius_km, args$seed))

  out <- bind_rows(Filter(Negate(is.null), results))
  out[!duplicated(paste(out$from, out$to)), ]
}

cat("Building Layer 2: community edges (parallel)...\n")
t0         <- proc.time()
comm_edges <- create_comm_edges_parallel(nodes, n_cores = 15L)
elapsed    <- round((proc.time() - t0)[["elapsed"]], 1)
cat(sprintf("  Community edges: %d (%.1f sec)\n", nrow(comm_edges), elapsed))
cat(sprintf("  Mean distance  : %.2f km\n", mean(comm_edges$distance_km)))
# ── Check Layer 3 HCW coverage in test subset ─────────────────
cat(sprintf("HCWs in test subset    : %d\n", nrow(test_hcw)))
cat(sprintf("Hospitals with HCWs    : %d\n", length(unique(test_hcw$hospital_id))))
cat(sprintf("Mean HCWs per hospital : %.1f\n",
            mean(table(test_hcw$hospital_id))))

# HCW-HCW edges: n*(n-1)/2 per hospital
# With ~124 HCWs spread across hospitals, edges would be sparse
hcw_per_hosp <- table(test_hcw$hospital_id)
cat(sprintf("\nHCWs per hospital distribution:\n"))
print(summary(as.numeric(hcw_per_hosp)))

# Expected edges if evenly distributed
cat(sprintf("\nExpected HCW-HCW edges : %d\n",
            sum(sapply(as.numeric(hcw_per_hosp),
                       function(n) n*(n-1)/2))))

# ── Layer 3: Healthcare edges ──────────────────────────────────
cat("Building Layer 3: healthcare edges...\n")

# (a) HCW-HCW edges within same hospital (full clique)
test_hcw <- hcw_nodes %>%
  semi_join(test_cells, by = c("cell_x", "cell_y"))

hcw_hcw_edges <- test_hcw %>%
  group_by(hospital_id) %>%
  group_map(~ {
    members <- .x$person_id
    if (length(members) < 2) return(NULL)
    pairs <- combn(members, 2)
    data.frame(from        = pairs[1, ],
               to          = pairs[2, ],
               hospital_id = .y$hospital_id,
               layer       = "healthcare_hcw",
               stringsAsFactors = FALSE)
  }) %>%
  bind_rows()

cat(sprintf("  HCW-HCW edges     : %d\n", nrow(hcw_hcw_edges)))
cat(sprintf("  Hospitals covered : %d\n", length(unique(hcw_hcw_edges$hospital_id))))

# (b) Non-HCW admission lookup: person_id → nearest hospital
# Activated during simulation when individual is hospitalized
hosp_coords   <- as.matrix(hf_hospital[, c("lon", "lat")])
non_hcw_nodes <- nodes %>% filter(!is_hcw)

nearest_hosp <- apply(as.matrix(non_hcw_nodes[, c("indiv_x", "indiv_y")]),
                      1, function(coord) {
                        dists <- sqrt((hosp_coords[, 1] - coord[1])^2 +
                                        (hosp_coords[, 2] - coord[2])^2)
                        which.min(dists)
                      })

admission_lookup <- data.frame(
  person_id   = non_hcw_nodes$person_id,
  hospital_id = hf_hospital$OBJECTID[nearest_hosp]
)

cat(sprintf("  Admission lookup  : %d individuals\n", nrow(admission_lookup)))

# ── Layer 4: Funeral extra community edges ────────────────────
# Additional contacts activated on death
# NB(mu = rowSums(prem_community)[age_group], size = 0.1), 10km radius
# No distance decay

# Precompute expected contacts per age group from Prem community matrix
prem_comm_mu <- rowSums(prem_mats$community)  # length 16, one per age group
cat("Prem community mu by age group:\n")
print(round(prem_comm_mu, 2))

create_funeral_extra_edges <- function(node_df, prem_mu, size = 0.1,
                                       radius_km = 10, n_cores = 15L, seed = 77) {
  N         <- nrow(node_df)
  coords    <- as.matrix(node_df[, c("indiv_x", "indiv_y")])
  hh_ids    <- node_df$hh_id
  age_grps  <- node_df$age_group
  pids      <- node_df$person_id
  chunks    <- split(seq_len(N), ceiling(seq_len(N) / ceiling(N / n_cores)))

  worker <- function(chunk_ids, coords, hh_ids, age_grps, pids,
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
      # Use Prem community mu for this person's age group
      mu_i <- prem_mu[age_grps[i]]
      n_c  <- rnbinom(1, size = size, mu = mu_i)
      if (n_c == 0) next

      lon_i   <- coords[i, 1]; lat_i <- coords[i, 2]
      lon_rad <- radius_km / (111.0 * cos(lat_i * pi / 180))

      cands <- which(abs(coords[, 2] - lat_i) <= lat_rad &
                       abs(coords[, 1] - lon_i) <= lon_rad &
                       hh_ids != hh_ids[i] & seq_len(nrow(coords)) != i)
      if (length(cands) == 0) next

      dists <- hav_km(lon_i, lat_i, coords[cands, 1], coords[cands, 2])
      cands <- cands[dists <= radius_km]
      if (length(cands) == 0) next

      n_s     <- min(n_c, length(cands))
      sampled <- sample(cands, n_s, replace = FALSE)

      # Keep i < j only to avoid duplicates
      keep <- sampled > i
      if (any(keep)) {
        from_v <- c(from_v, rep(i, sum(keep)))
        to_v   <- c(to_v, sampled[keep])
      }
    }
    if (length(from_v) == 0) return(NULL)
    data.frame(from  = from_v,
               to    = to_v,
               layer = "funeral_extra",
               stringsAsFactors = FALSE)
  }

  cl <- makeCluster(n_cores, type = "PSOCK")
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl,
                varlist = c("coords", "hh_ids", "age_grps", "pids",
                            "prem_mu", "size", "radius_km"),
                envir = environment())
  chunk_args <- lapply(seq_along(chunks), function(ci)
    list(chunk_ids = chunks[[ci]], seed = seed + ci))

  results <- parLapply(cl, chunk_args, function(args)
    worker(args$chunk_ids, coords, hh_ids, age_grps, pids,
           prem_mu, size, radius_km, args$seed))

  out <- bind_rows(Filter(Negate(is.null), results))
  out[!duplicated(paste(out$from, out$to)), ]
}

cat("Building Layer 4: funeral extra edges (parallel)...\n")
t0             <- proc.time()
funeral_edges  <- create_funeral_extra_edges(nodes, prem_comm_mu, n_cores = 15L)
elapsed        <- round((proc.time() - t0)[["elapsed"]], 1)
cat(sprintf("  Funeral extra edges: %d (%.1f sec)\n", nrow(funeral_edges), elapsed))

# ── Summary ────────────────────────────────────────────────────
cat("\n=== Network Layer Summary (test subset) ===\n")
cat(sprintf("  Nodes              : %d\n",   nrow(nodes)))
cat(sprintf("  Layer 1 HH edges   : %d\n",   nrow(hh_edges)))
cat(sprintf("  Layer 2 comm edges : %d\n",   nrow(comm_edges)))
cat(sprintf("  Layer 3 HCW edges  : %d\n",   nrow(hcw_hcw_edges)))
cat(sprintf("  Layer 3 admission  : %d\n",   nrow(admission_lookup)))
cat(sprintf("  Layer 4 funeral    : %d\n",   length(funeral_pool)))

# ── Layer 3 Extended: hospital → HCW list lookup ──────────────
# At simulation time:
#   - person gets hospitalized
#   - look up their hospital_id from admission_lookup
#   - look up HCW list from hospital_hcw_lookup
#   - activate edges: person ↔ all HCWs in that hospital

# Build hospital → HCW list (from full HCW data, not just test subset)
hospital_hcw_lookup <- hcw_nodes %>%
  group_by(hospital_id) %>%
  summarise(hcw_list = list(person_id), .groups = "drop")

cat(sprintf("Hospitals in lookup : %d\n", nrow(hospital_hcw_lookup)))
cat(sprintf("Total HCWs indexed  : %d\n",
            sum(sapply(hospital_hcw_lookup$hcw_list, length))))
cat(sprintf("Mean HCWs/hospital  : %.1f\n",
            mean(sapply(hospital_hcw_lookup$hcw_list, length))))
cat(sprintf("Max HCWs/hospital   : %d\n",
            max(sapply(hospital_hcw_lookup$hcw_list, length))))

# Merge into admission_lookup for easy access at simulation time
admission_lookup <- admission_lookup %>%
  left_join(hospital_hcw_lookup, by = "hospital_id")

cat(sprintf("\nAdmission lookup with HCW lists: %d individuals\n",
            nrow(admission_lookup)))
cat("Sample:\n")
print(head(admission_lookup %>% select(person_id, hospital_id) %>%
             mutate(n_hcw = sapply(admission_lookup$hcw_list, length))))

# ── Save all network layers ────────────────────────────────────
output_dir <- "output/network"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Save node dataframe (includes HCW flag, age_group, hospital assignment)
saveRDS(pop_nodes,
        file.path(output_dir, "Kinshasa_nodes.rds"))

# Save Layer 1: household edges
saveRDS(hh_edges,
        file.path(output_dir, "Kinshasa_layer1_household.rds"))

# Save Layer 2: community edges (distance stored, kernel applied later)
saveRDS(comm_edges,
        file.path(output_dir, "Kinshasa_layer2_community.rds"))

# Save Layer 3: HCW-HCW edges (always active)
saveRDS(hcw_hcw_edges,
        file.path(output_dir, "Kinshasa_layer3_hcw_edges.rds"))

# Save Layer 3: admission lookup (person → hospital + HCW list, activated on hospitalization)
saveRDS(admission_lookup,
        file.path(output_dir, "Kinshasa_layer3_admission_lookup.rds"))

# Save Layer 4: funeral pool (activated on death)
saveRDS(funeral_pool,
        file.path(output_dir, "Kinshasa_layer4_funeral_pool.rds"))

# Save HCW node info separately for reference
saveRDS(hcw_nodes,
        file.path(output_dir, "Kinshasa_hcw_nodes.rds"))

cat("Saved files:\n")
for (f in list.files(output_dir)) {
  size_mb <- round(file.info(file.path(output_dir, f))$size / 1e6, 1)
  cat(sprintf("  %-50s %.1f MB\n", f, size_mb))
}

