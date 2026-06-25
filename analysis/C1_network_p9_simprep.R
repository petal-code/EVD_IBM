# ==============================================================================
# C1_network_p9_simprep.R
# Purpose:
#   Build neighbor adjacency lists from p7 network output for simulation.
#   Splits Layer 2 into 6 sub-lists by stratum (daily/weekly/monthly)
#   and contact type (close/physical via is_physical flag).
#   Saves result as {tag}_sim_prep.rds for fast O(1) lookup in simulation.
#
# Output: output/network/{tag}_sim_prep.rds
#   list(N, pid_to_idx,
#        hh_nbrs,
#        comm_close_daily_nbrs, comm_close_weekly_nbrs, comm_close_monthly_nbrs,
#        comm_phys_daily_nbrs,  comm_phys_weekly_nbrs,  comm_phys_monthly_nbrs,
#        hcw_nbrs, adm_nbrs)
# ==============================================================================

library(dplyr)

# ==============================================================================
# [Configuration]
# ==============================================================================

network_dir <- "output/network"
case_tag    <- "case1_1M"   # "case1_1M", "case2_Ituri", "case3_Kivu"

# ==============================================================================
# [Section 1] Load network outputs
# ==============================================================================

message(sprintf("=== Building sim_prep: %s ===", case_tag))

nodes   <- readRDS(file.path(network_dir, sprintf("%s_nodes.rds",            case_tag)))
layer1  <- readRDS(file.path(network_dir, sprintf("%s_layer1_household.rds", case_tag)))
layer2d <- readRDS(file.path(network_dir, sprintf("%s_layer2_daily.rds",     case_tag)))
layer2w <- readRDS(file.path(network_dir, sprintf("%s_layer2_weekly.rds",    case_tag)))
layer2m <- readRDS(file.path(network_dir, sprintf("%s_layer2_monthly.rds",   case_tag)))
layer3h <- readRDS(file.path(network_dir, sprintf("%s_layer3_hcw_edges.rds", case_tag)))
layer3a <- readRDS(file.path(network_dir, sprintf("%s_layer3_admission.rds", case_tag)))

N <- nrow(nodes)
message(sprintf("  N = %d", N))

# person_id -> row index (O(1) lookup)
pid_to_idx <- integer(max(nodes$person_id) + 1L)
pid_to_idx[nodes$person_id + 1L] <- seq_len(N)

# ==============================================================================
# [Helper] Build adjacency list from edge dataframe
# Undirected: each edge (from, to) adds both directions
# Optional filter on a column
# ==============================================================================

build_adj_list <- function(edges, n, pid_to_idx,
                           filter_col = NULL, filter_val = NULL) {
  result <- vector("list", n)
  for (i in seq_len(n)) result[[i]] <- integer(0)

  if (!is.null(filter_col) && !is.null(filter_val))
    edges <- edges[edges[[filter_col]] == filter_val, ]

  if (nrow(edges) == 0L) return(result)

  fi <- pid_to_idx[edges$from + 1L]
  ti <- pid_to_idx[edges$to   + 1L]

  for (k in seq_len(nrow(edges))) {
    if (is.na(fi[k]) || fi[k] == 0L || is.na(ti[k]) || ti[k] == 0L) next
    result[[fi[k]]] <- c(result[[fi[k]]], ti[k])
    result[[ti[k]]] <- c(result[[ti[k]]], fi[k])
  }

  lapply(result, unique)
}

# ==============================================================================
# [Section 2] Build neighbor lists
# ==============================================================================

message("  Building Layer 1 (household)...")
hh_nbrs <- build_adj_list(layer1, N, pid_to_idx)

message("  Building Layer 2 daily close...")
comm_close_daily_nbrs <- build_adj_list(layer2d, N, pid_to_idx,
                                        "is_physical", 0L)
message("  Building Layer 2 daily physical...")
comm_phys_daily_nbrs  <- build_adj_list(layer2d, N, pid_to_idx,
                                        "is_physical", 1L)
message("  Building Layer 2 weekly close...")
comm_close_weekly_nbrs <- build_adj_list(layer2w, N, pid_to_idx,
                                         "is_physical", 0L)
message("  Building Layer 2 weekly physical...")
comm_phys_weekly_nbrs  <- build_adj_list(layer2w, N, pid_to_idx,
                                         "is_physical", 1L)
message("  Building Layer 2 monthly close...")
comm_close_monthly_nbrs <- build_adj_list(layer2m, N, pid_to_idx,
                                          "is_physical", 0L)
message("  Building Layer 2 monthly physical...")
comm_phys_monthly_nbrs  <- build_adj_list(layer2m, N, pid_to_idx,
                                          "is_physical", 1L)

message("  Building Layer 3 HCW-HCW...")
hcw_nbrs <- build_adj_list(layer3h, N, pid_to_idx)

message("  Building Layer 3 admission lookup...")
adm_nbrs <- vector("list", N)
for (i in seq_len(N)) adm_nbrs[[i]] <- integer(0)

for (k in seq_len(nrow(layer3a))) {
  pid  <- layer3a$person_id[k]
  hcws <- layer3a$hcw_list[[k]]
  if (is.null(hcws) || length(hcws) == 0L) next
  idx      <- pid_to_idx[pid + 1L]
  hcw_idxs <- pid_to_idx[hcws + 1L]
  hcw_idxs <- hcw_idxs[!is.na(hcw_idxs) & hcw_idxs > 0L]
  if (!is.na(idx) && idx > 0L && length(hcw_idxs) > 0L)
    adm_nbrs[[idx]] <- hcw_idxs
}

# ==============================================================================
# [Section 3] Diagnostics
# ==============================================================================

mean_deg <- function(lst) mean(sapply(lst, length))
message("\n  === Neighbor list summary ===")
message(sprintf("  HH                 : %.2f", mean_deg(hh_nbrs)))
message(sprintf("  Comm close daily   : %.2f", mean_deg(comm_close_daily_nbrs)))
message(sprintf("  Comm close weekly  : %.2f", mean_deg(comm_close_weekly_nbrs)))
message(sprintf("  Comm close monthly : %.2f", mean_deg(comm_close_monthly_nbrs)))
message(sprintf("  Comm phys daily    : %.2f", mean_deg(comm_phys_daily_nbrs)))
message(sprintf("  Comm phys weekly   : %.2f", mean_deg(comm_phys_weekly_nbrs)))
message(sprintf("  Comm phys monthly  : %.2f", mean_deg(comm_phys_monthly_nbrs)))
message(sprintf("  HCW-HCW            : %.2f", mean_deg(hcw_nbrs)))
message(sprintf("  Admission HCWs     : %.2f", mean_deg(adm_nbrs)))

# ==============================================================================
# [Section 4] Save
# ==============================================================================

sim_prep <- list(
  N                       = N,
  pid_to_idx              = pid_to_idx,
  hh_nbrs                 = hh_nbrs,
  comm_close_daily_nbrs   = comm_close_daily_nbrs,
  comm_close_weekly_nbrs  = comm_close_weekly_nbrs,
  comm_close_monthly_nbrs = comm_close_monthly_nbrs,
  comm_phys_daily_nbrs    = comm_phys_daily_nbrs,
  comm_phys_weekly_nbrs   = comm_phys_weekly_nbrs,
  comm_phys_monthly_nbrs  = comm_phys_monthly_nbrs,
  hcw_nbrs                = hcw_nbrs,
  adm_nbrs                = adm_nbrs
)

out_path <- file.path(network_dir, sprintf("%s_sim_prep.rds", case_tag))
saveRDS(sim_prep, out_path)
message(sprintf("\n  Saved: %s", out_path))
message("=== Done ===")
