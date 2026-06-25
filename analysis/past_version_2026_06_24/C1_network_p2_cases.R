# ==============================================================================
# C1_network_p2_cases.R
# Purpose:
#   Build population data for 3 simulation cases
#
# Cases:
#   Case 1 (~1M)      : Bunia-centered square, half=0.19 deg (~42km)
#   Case 2 (Ituri)    : Ituri province
#   Case 3 (Ituri+Kivu): Ituri + Nord-Kivu + Sud-Kivu
#
# Outputs:
#   output/popdata/{tag}_pop.rds
#   figure/C1_p2/{tag}.png
#   figure/C1_p2/overview_all_cases.png
# ==============================================================================

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(ggspatial)

# ==============================================================================
# [Configuration]
# ==============================================================================

cases <- list(
  list(tag        = "case1_1M",
       name1      = "Ituri_Aru_centered",
       use_square = TRUE,
       half       = 0.19),
  list(tag        = "case2_Ituri",
       name1      = "Ituri",
       use_square = FALSE,
       provinces  = "Ituri"),
  list(tag        = "case3_Kivu",
       name1      = "Nord-Kivu + Sud-Kivu",
       use_square = FALSE,
       provinces  = c("Ituri", "Nord-Kivu", "Sud-Kivu"))
)

dir.create("output/popdata", showWarnings = FALSE, recursive = TRUE)
dir.create("figure/C1_p2",   showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# [Section 1] Load WorldPop rasters
# ==============================================================================

message("=== Section 1: Loading WorldPop rasters ===")

files_all <- list.files("data/worldpop/DRC_1km", pattern = "\\.tif$",
                        full.names = TRUE)
files_t   <- files_all[grepl("/cod_t_", files_all)]
message(sprintf("  Age files (_t_): %d", length(files_t)))

pop_rast <- NULL
for (f in files_t) {
  r <- rast(f)
  if (is.null(pop_rast)) pop_rast <- r else pop_rast <- pop_rast + r
}
message(sprintf("  Total DRC pop: %.0f",
                global(pop_rast, "sum", na.rm = TRUE)[[1]]))

age_layers <- lapply(files_t, function(f) {
  age <- as.integer(sub("cod_t_(\\d+)_.*", "\\1", basename(f)))
  list(f = f, age = age)
})

# ==============================================================================
# [Section 2] Load shapefiles + find Bunia center
# ==============================================================================

message("\n=== Section 2: Finding Bunia center ===")

cod1_sf <- st_read("data/shpmap/gadm41_COD_1.shp", quiet = TRUE)
cod2_sf <- st_read("data/shpmap/gadm41_COD_2.shp", quiet = TRUE)

pop_df <- as.data.frame(pop_rast, xy = TRUE) |>
  setNames(c("x", "y", "pop")) |>
  filter(!is.na(pop), pop > 0)

bunia_sf   <- cod2_sf[cod2_sf$NAME_1 == "Ituri" &
                        grepl("Bunia", cod2_sf$NAME_2, ignore.case = TRUE), ]
bbox_bunia <- st_bbox(bunia_sf)
pop_bunia  <- pop_df |>
  filter(x >= bbox_bunia["xmin"], x <= bbox_bunia["xmax"],
         y >= bbox_bunia["ymin"], y <= bbox_bunia["ymax"])

cx <- pop_bunia$x[which.max(pop_bunia$pop)]
cy <- pop_bunia$y[which.max(pop_bunia$pop)]
message(sprintf("  Bunia center: lon=%.4f, lat=%.4f", cx, cy))

# ==============================================================================
# [Section 3] Process each case
# ==============================================================================

message("\n=== Section 3: Processing cases ===")

for (case in cases) {
  tag      <- case$tag
  rds_path <- file.path("output/popdata", paste0(tag, "_pop.rds"))
  png_path <- file.path("figure/C1_p2",   paste0(tag, ".png"))

  if (file.exists(rds_path) && file.exists(png_path)) {
    message(sprintf("  SKIP (exists): %s", tag)); next
  }

  message(sprintf("\n  Building: %s", tag))

  if (case$use_square) {
    # ── Case 1: Bunia-centered square ───────────────────────────
    half <- case$half
    xmin <- cx - half; xmax <- cx + half
    ymin <- cy - half; ymax <- cy + half
    pop_sq  <- pop_df |> filter(x >= xmin, x <= xmax, y >= ymin, y <= ymax)
    km_side <- round(half * 2 * 111, 0)

    age_list <- lapply(age_layers, function(lyr) {
      r      <- rast(lyr$f)
      r_crop <- tryCatch(crop(r, ext(xmin, xmax, ymin, ymax)),
                         error = function(e) NULL)
      if (is.null(r_crop)) return(NULL)
      as.data.frame(r_crop, xy = TRUE) |>
        setNames(c("x","y","count")) |>
        filter(!is.na(count), count > 0,
               x >= xmin, x <= xmax, y >= ymin, y <= ymax) |>
        mutate(age_lower = lyr$age)
    })
    pop_age   <- bind_rows(Filter(Negate(is.null), age_list))
    total_pop <- sum(pop_sq$pop)
    n_cells   <- nrow(pop_sq)
    message(sprintf("    Pop: %.0f | Cells: %d | ~%dx%d km",
                    total_pop, n_cells, km_side, km_side))

    if (!file.exists(rds_path))
      saveRDS(list(name1=case$name1, name2=tag, safe_name=tag,
                   total_pop=total_pop, cx=cx, cy=cy, half_deg=half,
                   pop_total=pop_sq, pop_age=pop_age), rds_path)

    if (!file.exists(png_path)) {
      adm_sf <- st_make_valid(cod2_sf) |>
        st_intersection(st_as_sfc(st_bbox(
          c(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), crs=4326)))
      p <- ggplot() +
        geom_raster(data=pop_sq, aes(x=x,y=y,fill=log1p(pop))) +
        scale_fill_viridis_c(option="magma") +
        geom_sf(data=adm_sf, fill=NA, color="white",
                linewidth=0.2, alpha=0.5) +
        annotate("point", x=cx, y=cy,
                 color="cyan", size=3, shape=4, stroke=2) +
        annotation_scale(location="bl", style="ticks",
                         text_col="white", line_col="white") +
        annotation_north_arrow(location="tl", which_north="true",
                               style=north_arrow_fancy_orienteering(
                                 fill=c("white","grey40"), text_col="white")) +
        coord_sf(xlim=c(xmin,xmax), ylim=c(ymin,ymax)) +
        labs(title    = sprintf("%s — Bunia-centered square", tag),
             subtitle = sprintf("Pop: %.1fM | Cells: %d | ~%dx%d km",
                                total_pop/1e6, n_cells, km_side, km_side),
             x=NULL, y=NULL) +
        theme_void() +
        theme(legend.position="none",
              plot.title=element_text(size=13,face="bold",color="white"),
              plot.subtitle=element_text(size=10,color="grey70"),
              plot.background=element_rect(fill="grey10",color=NA),
              plot.margin=margin(8,8,8,8))
      ggsave(png_path, plot=p, width=8, height=8, dpi=300)
    }

  } else {
    # ── Case 2/3: province boundaries ───────────────────────────
    provs      <- case$provinces
    prov2_sf   <- st_make_valid(cod2_sf[cod2_sf$NAME_1 %in% provs, ])
    prov1_sf   <- st_make_valid(cod1_sf[cod1_sf$NAME_1 %in% provs, ])
    prov_union <- st_union(prov2_sf)
    bbox_p     <- st_bbox(prov_union)

    xmin <- as.numeric(bbox_p["xmin"]); xmax <- as.numeric(bbox_p["xmax"])
    ymin <- as.numeric(bbox_p["ymin"]); ymax <- as.numeric(bbox_p["ymax"])

    # Mask cells to province boundaries
    pop_bbox <- pop_df |> filter(x >= xmin, x <= xmax, y >= ymin, y <= ymax)
    pop_sf   <- st_as_sf(pop_bbox, coords=c("x","y"), crs=4326)
    in_prov  <- st_intersects(pop_sf, prov_union, sparse=FALSE)[,1]
    pop_sq   <- pop_bbox[in_prov, ]

    # Age breakdown — mask to province
    age_list <- lapply(age_layers, function(lyr) {
      r      <- rast(lyr$f)
      r_crop <- tryCatch(crop(r, ext(xmin, xmax, ymin, ymax)),
                         error = function(e) NULL)
      if (is.null(r_crop)) return(NULL)
      df <- as.data.frame(r_crop, xy=TRUE) |>
        setNames(c("x","y","count")) |>
        filter(!is.na(count), count > 0)
      if (nrow(df) == 0) return(NULL)
      df_sf  <- st_as_sf(df, coords=c("x","y"), crs=4326)
      in_p   <- st_intersects(df_sf, prov_union, sparse=FALSE)[,1]
      df[in_p, ] |> mutate(age_lower = lyr$age)
    })
    pop_age   <- bind_rows(Filter(Negate(is.null), age_list))
    total_pop <- sum(pop_sq$pop)
    n_cells   <- nrow(pop_sq)
    message(sprintf("    Provinces: %s | Pop: %.0f | Cells: %d",
                    paste(provs, collapse="+"), total_pop, n_cells))

    if (!file.exists(rds_path))
      saveRDS(list(name1=case$name1, name2=tag, safe_name=tag,
                   total_pop=total_pop,
                   pop_total=pop_sq, pop_age=pop_age), rds_path)

    if (!file.exists(png_path)) {
      p <- ggplot() +
        geom_raster(data=pop_sq, aes(x=x,y=y,fill=log1p(pop))) +
        scale_fill_viridis_c(option="magma") +
        geom_sf(data=prov2_sf, fill=NA, color="white",
                linewidth=0.2, alpha=0.4) +
        geom_sf(data=prov1_sf, fill=NA, color="white", linewidth=0.8) +
        annotation_scale(location="bl", style="ticks",
                         text_col="white", line_col="white") +
        annotation_north_arrow(location="tl", which_north="true",
                               style=north_arrow_fancy_orienteering(
                                 fill=c("white","grey40"), text_col="white")) +
        coord_sf(xlim=c(xmin,xmax), ylim=c(ymin,ymax)) +
        labs(title    = sprintf("%s — %s", tag, paste(provs, collapse="+")),
             subtitle = sprintf("Pop: %.1fM | Cells: %d",
                                total_pop/1e6, n_cells),
             x=NULL, y=NULL) +
        theme_void() +
        theme(legend.position="none",
              plot.title=element_text(size=13,face="bold",color="white"),
              plot.subtitle=element_text(size=10,color="grey70"),
              plot.background=element_rect(fill="grey10",color=NA),
              plot.margin=margin(8,8,8,8))
      ggsave(png_path, plot=p, width=8, height=10, dpi=300)
    }
  }
  message(sprintf("    Saved: %s", tag))
}

# ==============================================================================
# [Section 4] Overview map
# ==============================================================================

message("\n=== Section 4: Overview map ===")

ituri_l1_sf <- st_make_valid(cod1_sf[cod1_sf$NAME_1 == "Ituri", ])
kivu_l1_sf  <- st_make_valid(cod1_sf[cod1_sf$NAME_1 %in%
                                       c("Ituri","Nord-Kivu","Sud-Kivu"), ])

bbox_kv <- st_bbox(kivu_l1_sf)
pad     <- 0.5
xmin_ov <- as.numeric(bbox_kv["xmin"]) - pad
xmax_ov <- as.numeric(bbox_kv["xmax"]) + pad
ymin_ov <- as.numeric(bbox_kv["ymin"]) - pad
ymax_ov <- as.numeric(bbox_kv["ymax"]) + pad

pop_overview <- pop_df |>
  filter(x >= xmin_ov, x <= xmax_ov, y >= ymin_ov, y <= ymax_ov)

adm_ov <- st_make_valid(cod2_sf) |>
  st_intersection(st_as_sfc(st_bbox(
    c(xmin=xmin_ov,xmax=xmax_ov,ymin=ymin_ov,ymax=ymax_ov), crs=4326)))

case1_rect <- data.frame(
  xmin=cx-0.19, xmax=cx+0.19,
  ymin=cy-0.19, ymax=cy+0.19
)

p_overview <- ggplot() +
  geom_raster(data = pop_overview,
              aes(x = x, y = y, fill = log1p(pop))) +
  scale_fill_viridis_c(option = "magma") +
  # Level 2 admin background lines
  geom_sf(data = adm_ov, fill = NA,
          color = "black", linewidth = 0.12, alpha = 0.25) +

  # Case 3: Ituri+Kivu province outline (green dashed, top layer)
  geom_sf(data = kivu_l1_sf, fill = NA,
          color = "green", linewidth = 1, linetype = "solid") +
  # Case 2: Ituri province outline (gold solid)
  geom_sf(data = ituri_l1_sf, fill = NA,
          color = "black", linewidth = 1, linetype = "dashed") +
  # Case 1: 1M square (cyan solid)
  geom_rect(data = case1_rect,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            color = "cyan", fill = NA, linewidth = 1.2) +
  # Centroid marker
  # annotate("point", x=cx, y=cy,
  #          color="black", size=3, shape=4, stroke=2) +
  annotation_scale(location="bl", width_hint=0.2,
                   style="ticks", text_col="black", line_col="black") +
  annotation_north_arrow(location="tr", which_north="true",
                         style=north_arrow_fancy_orienteering(
                           fill=c("white","black"), text_col="black")) +
  coord_sf(xlim=c(xmin_ov, xmax_ov), ylim=c(ymin_ov, ymax_ov)) +
  labs(
    x=NULL, y=NULL) +
  theme_void() +
  theme(plot.title      = element_text(size=14, face="bold", color="black"),
        plot.subtitle   = element_text(size=9, color="grey40"),
        plot.background = element_rect(fill="white", color=NA),
        legend.position = "none",
        plot.margin     = margin(8,8,8,8))

ggsave("figure/C1_p2/overview_all_cases.png",
       plot=p_overview, width=4, height=6, dpi=300)
message("  Saved: figure/C1_p2/overview_all_cases.png")

# ==============================================================================
# [Done]
# ==============================================================================

message("\n=== Done ===")
message("  RDS → output/popdata/cases/")
message("  PNG → figure/C1_p2/cases/")
