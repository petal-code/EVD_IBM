# ==============================================================================
# C1_network_p4_gridhousehold.R
# Purpose:
#   For each Level 2 territory:
#     1. Load age-stratified population data from p2 (pop_age)
#     2. Poisson rounding per age group per cell
#     3. Expand to individual level (one row per person, age assigned within bin)
#     4. Assign individuals to DHS households per cell
#     5. Save synthetic population
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
  "Tshuapa"        = "equateur"
)

output_dir  <- "output/household"
seed        <- 42L
n_test      <- NULL   # Set to 3L for testing, NULL for full run

# ==============================================================================
# [Batch control] — Change ONLY this number (1 to n_batches)
# ==============================================================================
batch_id <- 10L   # ← Change this per RStudio instance (1, 2, 3, ..., 10)
n_batches <- 10L  # Total number of RStudio instances

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# [Section 1] Load shared data
# ==============================================================================

message("=== Section 1: Loading shared data ===")

province_hh_list <- readRDS("output/household/province_household_list.rds")

pop_files <- list.files("output/popdata", pattern = "_pop\\.rds$",
                        full.names = TRUE)
message(sprintf("  Territory files: %d", length(pop_files)))

make_out_path <- function(pop_file, out_dir) {
  tag <- sub("_pop\\.rds$", "", basename(pop_file))
  file.path(out_dir, paste0(tag, "_synthetic_population.rds"))
}

pop_files_todo <- pop_files[!file.exists(
  sapply(pop_files, make_out_path, out_dir = output_dir)
)]

# Split into batches — each RStudio instance handles one batch
all_todo   <- pop_files_todo
batch_idx  <- seq(batch_id, length(all_todo), by = n_batches)
pop_files_todo <- all_todo[batch_idx]

if (!is.null(n_test)) {
  pop_files_todo <- head(pop_files_todo, n_test)
  message(sprintf("  TEST MODE: %d territories", n_test))
}

message(sprintf("  Batch %d / %d : %d territories assigned",
                batch_id, n_batches, length(pop_files_todo)))

if (length(pop_files_todo) == 0) {
  message("All done!"); stop("Done", call. = FALSE)
}

# ==============================================================================
# [Helper] Sample DHS households to fill a cell's population
# ==============================================================================

sample_households <- function(target_n, dhs_hh, seed_val) {
  set.seed(seed_val)

  # Filter households with valid age data
  valid  <- sapply(dhs_hh$hh_ages, function(x) ncol(x) > 0)
  dhs_hh <- dhs_hh[valid, ]
  if (nrow(dhs_hh) == 0 || target_n == 0)
    return(data.frame(hh_id = integer(0), age = integer(0)))

  households <- list()
  total <- 0L; hh_id <- 1L

  while (total < target_n) {
    idx  <- sample(nrow(dhs_hh), 1,
                   prob = dhs_hh$hv005 / sum(dhs_hh$hv005))
    # hh_ages[[idx]]: 1-row data.frame, one column per member
    ages <- as.numeric(dhs_hh$hh_ages[[idx]])
    ages <- ages[!is.na(ages)]
    if (length(ages) == 0) next

    households[[hh_id]] <- data.frame(hh_id = hh_id, age = ages)
    total <- total + length(ages)
    hh_id <- hh_id + 1L
  }

  result <- bind_rows(households)
  result[seq_len(min(nrow(result), target_n)), ]
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
  pop_age$count_int <- rpois(nrow(pop_age), lambda = pmax(pop_age$count, 0))

  # ── Step 2: Expand to individual level ───────────────────────
  # One row per person, age randomly assigned within 5-year bin
  pop_indiv <- pop_age %>%
    filter(count_int > 0) %>%
    uncount(count_int) %>%
    mutate(
      age = age_lower + sample(0:4, n(), replace = TRUE)
    ) %>%
    select(x, y, age)

  if (nrow(pop_indiv) == 0) { n_skipped <- n_skipped + 1L; next }

  # ── Step 3: Assign individuals to households per cell ─────────
  # For each cell: take individuals already assigned there,
  # then assign them into DHS-sampled households
  cells <- pop_indiv %>%
    group_by(x, y) %>%
    summarise(cell_pop = n(), .groups = "drop")

  household_list <- lapply(seq_len(nrow(cells)), function(ci) {
    cell <- cells[ci, ]

    # Get individuals in this cell (already have correct age structure)
    indiv_cell <- pop_indiv %>%
      filter(x == cell$x, y == cell$y)

    # Sample household structure from DHS
    hh_dat <- sample_households(nrow(indiv_cell), dhs_hh,
                                seed_val = seed + ci * 1000L)
    if (nrow(hh_dat) == 0) return(NULL)

    # Replace DHS ages with WorldPop ages (preserve age structure from raster)
    # hh_dat gives household membership; use WorldPop individual ages
    hh_dat$age   <- indiv_cell$age[seq_len(nrow(hh_dat))]
    hh_dat$cell_x <- cell$x
    hh_dat$cell_y <- cell$y
    hh_dat$hh_id  <- paste0(tag, "_c", ci, "_hh", hh_dat$hh_id)
    hh_dat
  })

  result <- bind_rows(Filter(Negate(is.null), household_list))
  saveRDS(result, out_path)
  n_saved <- n_saved + 1L

  elapsed <- proc.time()[["elapsed"]] - t_start
  rate    <- n_saved / max(elapsed, 0.1)
  eta     <- round((n_total - fi) / max(rate, 0.01) / 60, 1)
  message(sprintf(
    "  [%d/%d] %s - %s | pop: %d | hh: %d | cells: %d | ETA %.1f min",
    fi, n_total, name1, name2,
    nrow(result), length(unique(result$hh_id)), nrow(cells), eta
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
    pop <- readRDS(f)
    tag <- sub("_synthetic_population\\.rds$", "", basename(f))

    cat(sprintf("\n── %s ──\n", tag))
    cat(sprintf("  Individuals : %d\n", nrow(pop)))
    cat(sprintf("  Households  : %d\n", length(unique(pop$hh_id))))
    cat(sprintf("  Cells       : %d\n",
                length(unique(paste(pop$cell_x, pop$cell_y)))))

    # Household size
    hh_sz <- pop %>% group_by(hh_id) %>% summarise(n=n(), .groups="drop")
    cat(sprintf("  HH size     : mean=%.1f | median=%.0f | max=%d\n",
                mean(hh_sz$n), median(hh_sz$n), max(hh_sz$n)))

    # Age distribution
    cat(sprintf("  Age         : mean=%.1f | median=%.0f\n",
                mean(pop$age, na.rm=TRUE), median(pop$age, na.rm=TRUE)))
    age_grp <- cut(pop$age,
                   breaks = c(0, 5, 15, 30, 45, 60, Inf),
                   labels = c("0-4","5-14","15-29","30-44","45-59","60+"),
                   right  = FALSE)
    cat("  Age groups (%):\n")
    print(round(prop.table(table(age_grp)) * 100, 1))

    # Cell population
    cell_pop <- pop %>%
      group_by(cell_x, cell_y) %>%
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
