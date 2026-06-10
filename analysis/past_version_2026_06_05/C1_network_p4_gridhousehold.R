# ==============================================================================
# C1_network_p4_gridhousehold.R
# Purpose:
#   For each Level 2 territory:
#     1. Load age-stratified population from p2 (pop_age)
#     2. Poisson rounding per age group per cell
#     3. Expand to individual level with WorldPop age structure
#     4. Assign individuals to DHS households (similarity-weighted by age group)
#     5. Save three tables: individuals, households, cells
#
# Output structure:
#   individuals : person_id, hh_id, age
#   households  : hh_id, cell_id
#   cells       : cell_id, x, y
# ==============================================================================

library(dplyr)
library(tidyr)

# ==============================================================================
# [Configuration]
# ==============================================================================

province_to_dhs <- list(
  "Bas-Uele"       = "orientale",
  "Équateur"       = "equateur",
  "Haut-Katanga"   = "katanga",
  "Haut-Lomami"    = "katanga",
  "Haut-Uele"      = "orientale",
  "Ituri"          = "orientale",
  "Kasaï"          = "kasai-occidental",
  "Kasaï-Central"  = "kasai-occidental",
  "Kasaï-Oriental" = "kasai-oriental",
  "Kinshasa"       = "kinshasa",
  "Kongo-Central"  = "bas-congo",
  "Kwango"         = "bandundu",
  "Kwilu"          = "bandundu",
  "Lomami"         = "kasai-oriental",
  "Lualaba"        = "katanga",
  "Mai-Ndombe"     = "bandundu",
  "Maniema"        = "maniema",
  "Mongala"        = "equateur",
  "Nord-Kivu"      = "nord-kivu",
  "Nord-Ubangi"    = "equateur",
  "Sankuru"        = "kasai-oriental",
  "Sud-Kivu"       = "sud-kivu",
  "Sud-Ubangi"     = "equateur",
  "Tanganyika"     = "katanga",
  "Tshopo"         = "orientale",
  "Tshuapa"        = "equateur",
  # Metro regions
  "Nord-Kivu + Sud-Kivu" = "nord-kivu"
)

output_dir <- "output/household"
seed       <- 42L
n_test     <- NULL    # Set to NULL for full run

# ==============================================================================
# [Batch control] — Change ONLY this number (1 to n_batches)
# ==============================================================================
batch_id  <- 1L
n_batches <- 1L

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# [Section 1] Load shared data
# ==============================================================================

message("=== Section 1: Loading shared data ===")

province_hh_list <- readRDS("output/household/province_household_list.rds")

pop_files <- c(
  list.files("output/popdata",       pattern = "_pop\\.rds$", full.names = TRUE),
  list.files("output/popdata/metro", pattern = "_pop\\.rds$", full.names = TRUE)
)

make_out_path <- function(pop_file, out_dir) {
  tag <- sub("_pop\\.rds$", "", basename(pop_file))
  file.path(out_dir, paste0(tag, "_synthetic_population.rds"))
}

pop_files_todo <- pop_files[!file.exists(
  sapply(pop_files, make_out_path, out_dir = output_dir)
)]

# Batch split
all_todo       <- pop_files_todo
batch_idx      <- seq(batch_id, length(all_todo), by = n_batches)
pop_files_todo <- all_todo[batch_idx]

if (!is.null(n_test)) {
  pop_files_todo <- head(pop_files_todo, n_test)
  message(sprintf("  TEST MODE: %d territories", n_test))
}

message(sprintf("  Batch %d/%d : %d territories",
                batch_id, n_batches, length(pop_files_todo)))

if (length(pop_files_todo) == 0) {
  message("All done!"); stop("Done", call. = FALSE)
}

# ==============================================================================
# [Section 2] Process each territory
# ==============================================================================

message("\n=== Section 2: Building synthetic populations ===")

n_total   <- length(pop_files_todo)
t_start   <- proc.time()[["elapsed"]]
n_saved   <- 0L
n_skipped <- 0L

for (fi in seq_along(pop_files_todo)) {

  pop_file <- pop_files_todo[fi]
  tag      <- sub("_pop\\.rds$", "", basename(pop_file))
  out_path <- make_out_path(pop_file, output_dir)

  if (file.exists(out_path)) { n_skipped <- n_skipped + 1L; next }

  pop_data <- readRDS(pop_file)
  name1    <- pop_data$name1
  name2    <- pop_data$name2

  dhs_key <- province_to_dhs[[name1]]
  if (is.null(dhs_key) || !dhs_key %in% names(province_hh_list)) {
    message(sprintf("  [%d/%d] SKIP (no DHS): %s - %s", fi, n_total, name1, name2))
    n_skipped <- n_skipped + 1L; next
  }

  dhs_hh  <- province_hh_list[[dhs_key]]
  pop_age <- pop_data$pop_age  # x, y, count, age_lower

  if (is.null(pop_age) || nrow(pop_age) == 0) {
    n_skipped <- n_skipped + 1L; next
  }

  # ── Step 1: Poisson rounding per age group per cell ──────────
  set.seed(seed + sum(utf8ToInt(tag)))
  pop_age$count_int <- rpois(nrow(pop_age),
                             lambda = pmax(pop_age$count, 0))

  # ── Step 2: Build cell table ──────────────────────────────────
  cell_tbl <- pop_age %>%
    group_by(x, y) %>%
    summarise(cell_pop = sum(count_int), .groups = "drop") %>%
    filter(cell_pop > 0) %>%
    mutate(cell_id = row_number())

  if (nrow(cell_tbl) == 0) { n_skipped <- n_skipped + 1L; next }

  # ── Step 3: Expand to individual level ───────────────────────
  # Add cell_id to pop_age, expand by count, assign age within bin
  # Use cell_idx (integer) instead of paste(x,y) for fast split
  pop_age_cells <- pop_age %>%
    filter(count_int > 0) %>%
    left_join(cell_tbl %>% select(x, y, cell_id), by = c("x", "y")) %>%
    mutate(row_id = row_number())

  pop_indiv <- pop_age_cells %>%
    uncount(count_int) %>%
    mutate(age = age_lower + sample(0:4, n(), replace = TRUE)) %>%
    select(cell_id, age, age_lower)

  # Fast split by integer cell_id — much faster than paste(x,y)
  indiv_by_cell <- split(seq_len(nrow(pop_indiv)), pop_indiv$cell_id)

  # ── Step 4: Pre-compute DHS household age-group profiles ──────
  # Done once per territory — reused for all cells
  dhs_hh_valid <- dhs_hh[sapply(dhs_hh$hh_ages,
                                function(x) ncol(x) > 0), ]

  # 16 x n_hh matrix: age-group profile per household
  dhs_age_profiles <- vapply(seq_len(nrow(dhs_hh_valid)), function(i) {
    ages <- as.numeric(dhs_hh_valid$hh_ages[[i]])
    ages <- ages[!is.na(ages)]
    if (length(ages) == 0L) return(rep(0, 16))
    grp  <- pmin(floor(ages / 5L), 15L)
    tabulate(grp + 1L, nbins = 16L) / length(ages)
  }, numeric(16))  # vapply faster than sapply

  base_weights  <- dhs_hh_valid$hv005 / sum(dhs_hh_valid$hv005)
  n_hh_dhs      <- nrow(dhs_hh_valid)

  # Pre-compute mean household size for pre-sampling
  mean_hh_size <- mean(sapply(seq_len(n_hh_dhs), function(i) {
    length(as.numeric(dhs_hh_valid$hh_ages[[i]]))
  }))

  # ── Step 5: Assign households to each cell ────────────────────
  all_persons    <- vector("list", nrow(cell_tbl))
  # Pre-allocate based on expected household count
  expected_n_hh  <- ceiling(sum(cell_tbl$cell_pop) / mean_hh_size) + 100L
  hh_id_vec      <- integer(expected_n_hh)
  hh_cell_vec    <- integer(expected_n_hh)
  hh_ptr         <- 0L  # Write pointer
  global_hh_id   <- 0L

  for (ci in seq_len(nrow(cell_tbl))) {
    cell_id    <- cell_tbl$cell_id[ci]
    indiv_idx  <- indiv_by_cell[[as.character(cell_id)]]
    if (is.null(indiv_idx) || length(indiv_idx) == 0L) next

    target_n   <- length(indiv_idx)
    ages_cell  <- pop_indiv$age[indiv_idx]

    # Cell age-group profile
    grp_cell     <- pmin(floor(ages_cell / 5L), 15L)
    cell_profile <- tabulate(grp_cell + 1L, nbins = 16L) / target_n

    # Similarity-weighted sampling weights
    sim <- colSums(cell_profile * dhs_age_profiles)
    sim <- pmax(sim, 0)
    w   <- base_weights * (sim + 1e-6)
    w   <- w / sum(w)

    # Pre-sample enough households at once
    set.seed(seed + ci * 1000L)
    n_presample <- ceiling(target_n / mean_hh_size) * 2L + 10L
    hh_idx_pool <- sample(n_hh_dhs, n_presample, replace = TRUE, prob = w)

    # Assign individuals to households
    # WorldPop ages are substituted per age group — preserving DHS household
    # age-group structure while using WorldPop age distribution within each group
    person_hh  <- integer(target_n)
    person_age <- ages_cell  # Start with WorldPop ages (will be rearranged)
    pos        <- 1L

    # Pre-build age-group pools from WorldPop individuals in this cell
    # Pool for each Prem age group (0-15): shuffled indices into ages_cell
    grp_cell  <- pmin(floor(ages_cell / 5L), 15L)
    age_pools <- lapply(0:15, function(g) {
      idx <- which(grp_cell == g)
      if (length(idx) > 0) sample(idx) else integer(0)  # Shuffle within group
    })
    age_pool_ptr <- integer(16)  # Write pointer per group

    for (hh_pool_i in hh_idx_pool) {
      if (pos > target_n) break

      dhs_ages <- as.numeric(dhs_hh_valid$hh_ages[[hh_pool_i]])
      dhs_ages <- dhs_ages[!is.na(dhs_ages)]
      hh_size  <- length(dhs_ages)
      if (hh_size == 0L) next

      global_hh_id <- global_hh_id + 1L
      end_pos      <- min(pos + hh_size - 1L, target_n)
      actual_size  <- end_pos - pos + 1L

      # Replace DHS ages with WorldPop ages from matching age group
      # If pool for that group is exhausted, fall back to any available age
      for (slot in seq_len(actual_size)) {
        dhs_grp  <- pmin(floor(dhs_ages[slot] / 5L), 15L)
        pool     <- age_pools[[dhs_grp + 1L]]
        ptr      <- age_pool_ptr[dhs_grp + 1L]

        if (ptr < length(pool)) {
          # Use WorldPop age from matching group
          age_pool_ptr[dhs_grp + 1L] <- ptr + 1L
          person_age[pos + slot - 1L] <- ages_cell[pool[ptr + 1L]]
        }
        # If pool exhausted: keep original WorldPop age already in person_age
      }

      person_hh[pos:end_pos] <- global_hh_id

      # Write to pre-allocated vectors; expand if needed
      hh_ptr <- hh_ptr + 1L
      if (hh_ptr > length(hh_id_vec)) {
        extra       <- ceiling(length(hh_id_vec) * 0.2)
        hh_id_vec   <- c(hh_id_vec,   integer(extra))
        hh_cell_vec <- c(hh_cell_vec, integer(extra))
      }
      hh_id_vec[hh_ptr]   <- global_hh_id
      hh_cell_vec[hh_ptr] <- cell_id

      pos <- end_pos + 1L
    }

    all_persons[[ci]] <- data.frame(
      hh_id = person_hh,
      age   = person_age
    )
  }

  # ── Step 6: Combine and save ──────────────────────────────────
  persons_df <- bind_rows(all_persons) %>%
    mutate(person_id = row_number()) %>%
    select(person_id, hh_id, age)

  # Trim pre-allocated vectors to actual size
  households_df <- data.frame(
    hh_id   = hh_id_vec[seq_len(hh_ptr)],
    cell_id = hh_cell_vec[seq_len(hh_ptr)]
  )

  result <- list(
    individuals = persons_df,    # person_id, hh_id, age
    households  = households_df, # hh_id, cell_id
    cells       = cell_tbl       # cell_id, x, y, cell_pop
  )

  saveRDS(result, out_path)
  n_saved <- n_saved + 1L

  elapsed <- proc.time()[["elapsed"]] - t_start
  rate    <- n_saved / max(elapsed, 0.1)
  eta     <- round((n_total - fi) / max(rate, 0.01) / 60, 1)
  message(sprintf(
    "  [%d/%d] %s - %s | indiv: %d | hh: %d | cells: %d | ETA %.1f min",
    fi, n_total, name1, name2,
    nrow(persons_df), nrow(households_df), nrow(cell_tbl), eta
  ))
}

# ==============================================================================
# [Section 3] Validation (test mode only)
# ==============================================================================

if (!is.null(n_test) && n_saved > 0) {
  message("\n=== Section 3: Validation ===")

  saved_files <- list.files(output_dir,
                            pattern = "_synthetic_population\\.rds$",
                            full.names = TRUE)
  check_files <- tail(saved_files, n_saved)

  for (f in check_files) {
    res  <- readRDS(f)
    tag  <- sub("_synthetic_population\\.rds$", "", basename(f))
    pers <- res$individuals
    hh   <- res$households
    cel  <- res$cells

    cat(sprintf("\n── %s ──\n", tag))
    cat(sprintf("  Individuals : %d\n", nrow(pers)))
    cat(sprintf("  Households  : %d\n", nrow(hh)))
    cat(sprintf("  Cells       : %d\n", nrow(cel)))

    # Household size distribution
    hh_sz <- pers %>% group_by(hh_id) %>% summarise(n=n(), .groups="drop")
    cat(sprintf("  HH size     : mean=%.1f | median=%.0f | max=%d\n",
                mean(hh_sz$n), median(hh_sz$n), max(hh_sz$n)))

    # Age distribution
    cat(sprintf("  Age         : mean=%.1f | median=%.0f\n",
                mean(pers$age, na.rm=TRUE), median(pers$age, na.rm=TRUE)))
    age_grp <- cut(pers$age,
                   breaks = c(0,5,15,30,45,60,Inf),
                   labels = c("0-4","5-14","15-29","30-44","45-59","60+"),
                   right  = FALSE)
    cat("  Age groups (%):\n")
    print(round(prop.table(table(age_grp)) * 100, 1))

    # Cell population
    cell_pop <- pers %>%
      left_join(hh, by = "hh_id") %>%
      group_by(cell_id) %>%
      summarise(n = n(), .groups = "drop")
    cat(sprintf("  Indiv/cell  : mean=%.1f | max=%d\n",
                mean(cell_pop$n), max(cell_pop$n)))
  }
}

# ==============================================================================
# [Done]
# ==============================================================================

elapsed_total <- round(proc.time()[["elapsed"]] - t_start, 1)
message("\n=== Done ===")
message(sprintf("  Saved   : %d / %d", n_saved, n_total))
message(sprintf("  Skipped : %d", n_skipped))
message(sprintf("  Time    : %.1f sec (%.1f min)",
                elapsed_total, elapsed_total / 60))
