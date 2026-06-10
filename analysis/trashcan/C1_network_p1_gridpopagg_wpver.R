library(httr)
library(rvest)

# ── Configuration ──────────────────────────────────────────────
base_url <- "https://data.worldpop.org/repo/wopr/COD/population/v4.4/Province/"
save_dir <- "data/gridpop"
dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

# ── Helper: fetch all href links from a URL ────────────────────
fetch_links <- function(url, pattern = NULL) {
  res  <- GET(url, user_agent("Mozilla/5.0"), timeout(30))
  page <- content(res, as = "text", encoding = "UTF-8")
  links <- read_html(page) |>
    html_elements("a") |>
    html_attr("href")
  if (!is.null(pattern)) links <- grep(pattern, links, value = TRUE)
  links
}

# ── Step 1: Get province list ──────────────────────────────────
provinces <- fetch_links(base_url)
provinces <- provinces[grepl("^[A-Z]", provinces)]  # keep only province directories
cat(sprintf("Found %d provinces\n", length(provinces)))

# ── Step 2: Build all .tif download URLs ───────────────────────
# File structure: Province/Kinshasa/COD_Kinshasa_province_population_v4.4_agesex/xxx.tif
get_agesex_urls <- function(province) {
  prov_name   <- gsub("/$", "", province)  # remove trailing slash
  agesex_url  <- sprintf("%s%sCOD_%s_province_population_v4.4_agesex/",
                         base_url, province, prov_name)

  files <- tryCatch(
    fetch_links(agesex_url, pattern = "\\.tif$"),
    error = function(e) {
      cat(sprintf("  [FAIL] %s: %s\n", province, e$message))
      character(0)
    }
  )

  paste0(agesex_url, files)
}

cat("Scanning agesex folders...\n")
all_urls <- unlist(lapply(provinces, function(p) {
  cat(sprintf("  Scanning: %s\n", p))
  get_agesex_urls(p)
}))

cat(sprintf("\nTotal .tif files found: %d\n", length(all_urls)))

# ── Step 3: Download all files ─────────────────────────────────
download_tif <- function(url, save_dir) {
  province <- url |>
    gsub(base_url, "", x = _) |>   # strip base URL
    strsplit("/") |>
    unlist() |>
    head(1)                          # extract province name

  filename <- basename(url)
  prov_dir <- file.path(save_dir, province)
  dir.create(prov_dir, showWarnings = FALSE)
  dest     <- file.path(prov_dir, filename)

  # Skip if already downloaded successfully
  if (file.exists(dest) && file.info(dest)$size > 10000) {
    cat(sprintf("  [SKIP] %s\n", filename))
    return(invisible(NULL))
  }

  res <- tryCatch(
    GET(url,
        user_agent("Mozilla/5.0"),
        write_disk(dest, overwrite = TRUE),
        timeout(300)),   # 5 min timeout per file
    error = function(e) {
      cat(sprintf("  [FAIL] %s: %s\n", filename, e$message))
      NULL
    }
  )

  # Validate file size
  if (!is.null(res) && file.info(dest)$size < 10000) {
    cat(sprintf("  [WARN] %s too small, may be corrupted\n", filename))
  }
}

cat("\nStarting downloads...\n")
for (i in seq_along(all_urls)) {
  cat(sprintf("[%d/%d] %s\n", i, length(all_urls), basename(all_urls[i])))
  download_tif(all_urls[i], save_dir)
}

cat("\nAll downloads complete!\n")

