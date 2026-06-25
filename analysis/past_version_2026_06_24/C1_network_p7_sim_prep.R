# ==============================================================================
# C1_network_p7_sim_prep.R
# Purpose:
#   Pre-process network outputs for simulation for all 3 cases:
#     1. Build neighbor lookup lists (hh, comm, hcw, admission)
#     2. Build ring networks (hh_ring, ring1, ring2)
#
# Inputs:
#   output/network/{tag}_nodes.rds
#   output/network/{tag}_layer1_household.rds
#   output/network/{tag}_layer2_community.rds
#   output/network/{tag}_layer3_hcw_edges.rds
#   output/network/{tag}_layer3_admission.rds
#
# Outputs:
#   output/network/sim_prep/{tag}_nbrs.rds
#     list(
#       N          : number of individuals
#       person_ids : integer vector (row index → person_id)
#       pid_to_idx : integer vector (person_id → row index)
#       hh_nbrs    : list[N] — household neighbors (row indices)
#       comm_nbrs  : list[N] — community neighbors (row indices)
#       hcw_nbrs   : list[N] — HCW-HCW neighbors (row indices)
#       adm_nbrs   : list[N] — admission HCW neighbors (row indices)
#       hh_ring    : list[N] — hh_nbrs only
#       ring1      : list[N] — hh + comm direct contacts
#       ring2      : list[N] — 2-hop expansion of ring1
#     )
# ==============================================================================

library(dplyr)
library(data.table)

# ==============================================================================
# [Configuration]
# ==============================================================================

network_dir  <- "output/network"
sim_prep_dir <- file.path(network_dir, "sim_prep")
dir.create(sim_prep_dir, showWarnings = FALSE, recursive = TRUE)

case_tags <- c("case1_1M", "case2_Ituri", "case3_Kivu")

# ==============================================================================
# [Helper] Build adjacency list from edge dataframe (vectorized)
# ==============================================================================

build_adj <- function(edges, pid_to_idx, N, col_from = "from", col_to = "to") {
  if (is.null(edges) || nrow(edges) == 0)
    return(vector("list", N))

  from_idx <- pid_to_idx[edges[[col_from]]]
  to_idx   <- pid_to_idx[edges[[col_to]]]

  # Both directions — use data.table for speed
  dt <- data.table(
    src  = c(from_idx, to_idx),
    dest = c(to_idx,   from_idx)
  )
  dt <- dt[!is.na(src) & !is.na(dest)]
  setkey(dt, src)

  adj <- vector("list", N)
  for (i in seq_len(N)) {
    rows <- dt[.(i), dest, nomatch = NULL]
    if (length(rows) > 0) adj[[i]] <- rows
  }
  adj
}

# ==============================================================================
# [Section] Process each case
# ==============================================================================

for (tag in case_tags) {

  out_path <- file.path(sim_prep_dir, sprintf("%s_sim_prep.rds", tag))
  if (file.exists(out_path)) {
    message(sprintf("SKIP (exists): %s", tag)); next
  }

  message(sprintf("\n=== %s ===", tag))
  t0 <- proc.time()[["elapsed"]]

  # ── Load network files ──────────────────────────────────────────────────────
  nodes  <- readRDS(file.path(network_dir, sprintf("%s_nodes.rds",            tag)))
  layer1 <- readRDS(file.path(network_dir, sprintf("%s_layer1_household.rds", tag)))
  layer2 <- readRDS(file.path(network_dir, sprintf("%s_layer2_community.rds", tag)))
  layer3h <- readRDS(file.path(network_dir, sprintf("%s_layer3_hcw_edges.rds", tag)))
  layer3a <- readRDS(file.path(network_dir, sprintf("%s_layer3_admission.rds", tag)))

  N          <- nrow(nodes)
  person_ids <- nodes$person_id
  pid_to_idx <- integer(max(person_ids))
  pid_to_idx[person_ids] <- seq_len(N)

  message(sprintf("  N = %d", N))

  # ── Build neighbor lookup lists ─────────────────────────────────────────────
  message("  Building hh_nbrs...")
  t1 <- proc.time()[["elapsed"]]
  hh_nbrs <- build_adj(layer1, pid_to_idx, N)
  message(sprintf("    Done: %.1f sec", proc.time()[["elapsed"]] - t1))

  message("  Building comm_nbrs...")
  t1 <- proc.time()[["elapsed"]]
  comm_nbrs <- build_adj(layer2, pid_to_idx, N)
  message(sprintf("    Done: %.1f sec", proc.time()[["elapsed"]] - t1))

  message("  Building hcw_nbrs...")
  t1 <- proc.time()[["elapsed"]]
  hcw_nbrs <- build_adj(layer3h, pid_to_idx, N)
  message(sprintf("    Done: %.1f sec", proc.time()[["elapsed"]] - t1))

  # Admission: person → HCW row indices (asymmetric)
  message("  Building adm_nbrs...")
  t1 <- proc.time()[["elapsed"]]
  adm_nbrs <- vector("list", N)
  for (k in seq_len(nrow(layer3a))) {
    i <- pid_to_idx[layer3a$person_id[k]]
    if (is.na(i) || i == 0L) next
    hcw_ids <- layer3a$hcw_list[[k]]
    if (length(hcw_ids) > 0)
      adm_nbrs[[i]] <- pid_to_idx[hcw_ids]
  }
  message(sprintf("    Done: %.1f sec", proc.time()[["elapsed"]] - t1))

  # ── Build ring networks ─────────────────────────────────────────────────────

  # ── Summary stats ───────────────────────────────────────────────────────────
  message(sprintf("  hh_nbrs  : mean=%.1f | max=%d",
                  mean(lengths(hh_nbrs)), max(lengths(hh_nbrs))))
  message(sprintf("  comm_nbrs: mean=%.1f | max=%d",
                  mean(lengths(comm_nbrs)), max(lengths(comm_nbrs))))
  message(sprintf("  hcw_nbrs : %d HCWs with connections",
                  sum(lengths(hcw_nbrs) > 0)))
  message(sprintf("  adm_nbrs : %d patients with HCW assignment",
                  sum(lengths(adm_nbrs) > 0)))

  # ── Save ────────────────────────────────────────────────────────────────────
  # ring1 and ring2 are computed on-demand in simulation:
  #   ring1[[i]] = unique(c(hh_nbrs[[i]], comm_nbrs[[i]]))
  #   ring2      = unique(unlist(ring1[ring1[[idx]]]))
  saveRDS(list(
    N          = N,
    person_ids = person_ids,
    pid_to_idx = pid_to_idx,
    hh_nbrs    = hh_nbrs,
    comm_nbrs  = comm_nbrs,
    hcw_nbrs   = hcw_nbrs,
    adm_nbrs   = adm_nbrs
  ), out_path)

  message(sprintf("  Saved: %s (%.1f sec total)",
                  basename(out_path),
                  proc.time()[["elapsed"]] - t0))
}

message("\n=== Done ===")
