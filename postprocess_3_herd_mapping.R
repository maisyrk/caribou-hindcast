# postprocess_3_herd_mapping.R ----------

# Description:
# Run this script to create the publication-style single-panel herd boundary map.
# It uses a raster template only to define the study-area footprint and map extent.
# The raster itself is not plotted.
#
# Required inputs:
# - HSI_output/HSI_2020_masked.tif (or another raster used as the footprint template)
# - data/post_processing/herds/Aire_repartition_populations_caribouForestier.shp
# - data/post_processing/herds/Caribou_range_boundary.shp
# - map_data/forest_limit/LIM_FOR_ATT_2018.shp
# - map_data/forest_limit/northern_extent_managed_forest.shp
# - map_data/natural_earth/ne_50m_ocean/ne_50m_ocean.shp
# - map_data/natural_earth/ne_50m_lakes/ne_50m_lakes.shp
#
# Expected outputs:
# - map_layouts/herd_boundaries_single_map_labeled_timber_limit.png
#
# Notes:
# - Herd polygons are filled by herd and labelled directly with herd codes.
# - The legend includes the study area and the timber allocation northern limit.
# - The herd-code lookup table is drawn below the map.

## 1. SETUP ======================================================================

### 1.1 Define project root ------------------------------------------------------

# This makes the script portable.
# Open or set the working directory to the packaged project folder.

get_project_root <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(dirname(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))))
  }
  
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    editor_path <- tryCatch(
      rstudioapi::getSourceEditorContext()$path,
      error = function(e) ""
    )
    
    if (nzchar(editor_path)) {
      return(dirname(dirname(normalizePath(editor_path, winslash = "/", mustWork = TRUE))))
    }
  }
  
  script_path <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  
  if (!is.null(script_path)) {
    return(dirname(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))))
  }
  
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  
  if (dir.exists(file.path(wd, "scripts")) && dir.exists(file.path(wd, "data"))) {
    return(wd)
  }
  
  stop(
    "Could not determine the project root automatically.\n",
    "Open the 'caribou_hindcast_packaged' folder as the working directory, ",
    "or run this script from that project root."
  )
}

project_root <- get_project_root()
setwd(project_root)

### 1.2 Load required packages ---------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggspatial)
  library(scales)
  library(RColorBrewer)
  library(ggrepel)
  library(patchwork)
  library(tibble)
})

### 1.3 Check optional package ---------------------------------------------------

if (!requireNamespace("nngeo", quietly = TRUE)) {
  stop(
    "Package 'nngeo' is required. Install it before running this script:\n",
    "install.packages('nngeo')"
  )
}

### 1.4 User inputs --------------------------------------------------------------

# Raster used only to define the study-area footprint and map extent.
# The raster is not plotted in the final figure.
r_template_path <- "HSI_output/HSI_2020_masked.tif"

# Herd shapefiles
qc_shp <- "data/post_processing/herds/Aire_repartition_populations_caribouForestier.shp"
on_shp <- "data/post_processing/herds/Caribou_range_boundary.shp"

# Timber allocation limits (QC + ON)
timber_limit_qc_path <- "map_data/forest_limit/LIM_FOR_ATT_2018.shp"
timber_limit_on_path <- "map_data/forest_limit/northern_extent_managed_forest.shp"

# Local Natural Earth layers
ocean_path <- "map_data/natural_earth/ne_50m_ocean/ne_50m_ocean.shp"
lakes_path <- "map_data/natural_earth/ne_50m_lakes/ne_50m_lakes.shp"

# Projection used for the main map
crs_map <- "EPSG:3979"

# Output
out_png <- "map_layouts/herd_boundaries_single_map_labeled_timber_limit.png"

# Map framing buffer around raster footprint
frame_buffer_m <- 200000
y_pad_m <- 25000

# Optional speed-up when polygonizing the raster footprint
foot_fact <- 2

# Transparency
herd_fill_alpha <- 0.60
herd_line_alpha <- 1.00

# Line widths
lw_country   <- 0.20
lw_prov_over <- 0.20
lw_qcon      <- 0.20
lw_outline   <- 0.30
lw_herd      <- 0.30
lw_timber    <- 0.40

# Label aesthetics
label_size <- 3.2

# Herd-key text aesthetics
key_text_size <- 2.9
key_lineheight <- 1.04

### 1.5 Create output folder -----------------------------------------------------

dir.create(dirname(out_png), recursive = TRUE, showWarnings = FALSE)

## 2. DEFINE LOOKUP TABLES =======================================================

### 2.1 Herd label to herd code map ----------------------------------------------

herd_code_map <- c(
  "Assinica" = "ASN",
  "Basse Cote Nord" = "BCN",
  "Caniapiscau" = "CAN",
  "Charlevoix" = "CVX",
  "Detour" = "DET",
  "Joir River" = "JR",
  "Lac Joseph" = "LJ",
  "Manicouagan" = "MAN",
  "Nottaway" = "NOT",
  "Outardes" = "OUT",
  "Pipmuacan" = "PIP",
  "Secteur Baie James" = "SBJ",
  "Secteur Matamec" = "SM",
  "Temiscamie" = "TEM",
  "Val Dor" = "VdO",
  "Berens" = "BER",
  "Brightsand" = "BRS",
  "Churchill" = "CHU",
  "Discontinuous Distribution" = "DD",
  "James Bay" = "JB",
  "Kesagami" = "KES",
  "Kinloch" = "KIN",
  "Lake Superior Coast" = "LSC",
  "Missisa" = "MIS",
  "Nipigon" = "NIP",
  "Ozhiski" = "OZH",
  "Pagwachuan" = "PAG",
  "Spirit" = "SPR",
  "Swan" = "SL",
  "Sydney" = "SYD"
)

### 2.2 Optional pretty names for the key ----------------------------------------

pretty_name_map <- c(
  "Val Dor" = "Val d'Or",
  "Basse Cote Nord" = "Basse-Cote-Nord"
)

## 3. DEFINE HELPER FUNCTIONS ====================================================

### 3.1 Safely read an sf file ---------------------------------------------------

safe_read_sf <- function(path, crs_out) {
  st_read(path, quiet = TRUE) |>
    st_make_valid() |>
    st_transform(crs_out)
}

### 3.2 Standardize herd-name field ----------------------------------------------

# Uses the first likely herd-name field it finds.
# If none is found, it creates fallback names.

standardize_herd_name <- function(x, province_tag) {
  candidate_fields <- c(
    "herd", "HERD",
    "herd_name", "HERD_NAME",
    "name", "NAME",
    "nom", "NOM",
    "population", "POPULATION",
    "range_name", "RANGE_NAME",
    "caribou_pop", "NOM_POP",
    "area_name", "AREA_NAME"
  )
  
  hit <- candidate_fields[candidate_fields %in% names(x)][1]
  
  if (!is.na(hit)) {
    x <- x |>
      mutate(herd_name = as.character(.data[[hit]]))
  } else {
    x <- x |>
      mutate(herd_name = paste0(province_tag, "_Herd_", row_number()))
  }
  
  x |>
    mutate(
      herd_name = trimws(herd_name),
      herd_name = if_else(
        is.na(herd_name) | herd_name == "",
        paste0(province_tag, "_Herd_", row_number()),
        herd_name
      )
    )
}

### 3.3 Clean herd names for display ---------------------------------------------

make_herd_label <- function(x) {
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x <- gsub("['`]", "", x, perl = TRUE)
  x <- gsub("&", " And ", x, perl = TRUE)
  x <- gsub("[^A-Za-z0-9]+", " ", x, perl = TRUE)
  x <- trimws(x)
  x <- gsub("\\s+", " ", x, perl = TRUE)
  x <- tools::toTitleCase(x)
  x
}

### 3.4 Normalize labels for robust matching -------------------------------------

normalize_key <- function(x) {
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "", x, perl = TRUE)
  x
}

### 3.5 Dissolve herd polygons by herd -------------------------------------------

dissolve_by_herd <- function(x) {
  x |>
    group_by(province, herd_name, herd_label) |>
    summarise(geometry = st_union(geometry), .groups = "drop") |>
    st_make_valid()
}

### 3.6 Clip a line layer to a province polygon ----------------------------------

clip_line_to_poly <- function(x, poly) {
  if (is.null(x) || nrow(x) == 0) return(NULL)
  
  x <- x |>
    st_make_valid() |>
    st_zm(drop = TRUE, what = "ZM")
  
  poly <- poly |>
    st_make_valid() |>
    st_zm(drop = TRUE, what = "ZM")
  
  out <- suppressWarnings(st_intersection(x, poly))
  if (nrow(out) == 0) return(NULL)
  
  out_ls <- suppressWarnings(st_collection_extract(out, "LINESTRING", warn = FALSE))
  if (nrow(out_ls) > 0) {
    return(st_sf(geometry = st_geometry(out_ls), crs = st_crs(x)))
  }
  
  out_ml <- suppressWarnings(st_collection_extract(out, "MULTILINESTRING", warn = FALSE))
  if (nrow(out_ml) == 0) return(NULL)
  
  out_ls <- suppressWarnings(st_cast(out_ml, "LINESTRING", warn = FALSE))
  if (nrow(out_ls) == 0) return(NULL)
  
  st_sf(geometry = st_geometry(out_ls), crs = st_crs(x))
}

### 3.7 Dissolve linework --------------------------------------------------------

dissolve_lines <- function(x) {
  if (is.null(x) || nrow(x) == 0) return(NULL)
  
  x <- x |>
    st_make_valid() |>
    st_zm(drop = TRUE, what = "ZM")
  
  g <- suppressWarnings(st_union(st_geometry(x)))
  
  if (!inherits(g, "sfc")) {
    g <- st_sfc(g, crs = st_crs(x))
  }
  
  st_sf(geometry = g, crs = st_crs(x))
}

### 3.8 Build herd-code key table ------------------------------------------------

make_key_table <- function(x, code_lookup, pretty_name_map) {
  x |>
    st_drop_geometry() |>
    distinct(province, herd_label) |>
    mutate(label_key = normalize_key(herd_label)) |>
    left_join(code_lookup |> select(label_key, herd_code), by = "label_key") |>
    mutate(
      herd_code = if_else(is.na(herd_code), herd_label, herd_code),
      full_name = recode(herd_label, !!!pretty_name_map, .default = herd_label),
      key_line  = paste0(herd_code, " - ", full_name)
    ) |>
    arrange(herd_code) |>
    select(province, herd_code, full_name, key_line)
}

### 3.9 Split a character vector into two columns --------------------------------

split_two_cols <- function(x) {
  n <- length(x)
  n1 <- ceiling(n / 2)
  
  list(
    col1 = x[seq_len(n1)],
    col2 = if (n1 < n) x[(n1 + 1):n] else character(0)
  )
}

## 4. LOAD RASTER TEMPLATE =======================================================

### 4.1 Read the raster footprint template ---------------------------------------

r_template <- rast(r_template_path)

### 4.2 Reproject template if needed ---------------------------------------------

r_crs <- sf::st_crs(crs(r_template))
target_crs <- sf::st_crs(crs_map)

if (is.na(r_crs) || r_crs != target_crs) {
  r_template <- project(r_template, crs_map)
}

## 5. LOAD BASEMAP DATA ==========================================================

### 5.1 Load local Natural Earth layers ------------------------------------------

ocean <- st_read(ocean_path, quiet = TRUE) |>
  st_transform(crs_map)

lakes <- st_read(lakes_path, quiet = TRUE) |>
  st_transform(crs_map)

### 5.2 Load country and province boundaries -------------------------------------

countries <- ne_countries(scale = "medium", returnclass = "sf") |>
  filter(admin %in% c("Canada", "United States of America")) |>
  st_transform(crs_map)

provinces <- ne_states(country = "canada", returnclass = "sf") |>
  filter(name_en %in% c("Ontario", "Quebec")) |>
  st_transform(crs_map)

qc <- provinces |> filter(name_en == "Quebec")
on <- provinces |> filter(name_en == "Ontario")

### 5.3 Build QC-ON border line --------------------------------------------------

qc_on_border <- st_intersection(st_boundary(qc), st_boundary(on)) |>
  st_collection_extract("LINESTRING")

## 6. LOAD AND PREPARE TIMBER LIMIT ==============================================

### 6.1 Read raw timber limit layers ---------------------------------------------

timber_qc_raw <- safe_read_sf(timber_limit_qc_path, crs_map)
timber_on_raw <- safe_read_sf(timber_limit_on_path, crs_map)

### 6.2 Clip each timber line to its province ------------------------------------

timber_qc <- dissolve_lines(clip_line_to_poly(timber_qc_raw, qc))
timber_on <- dissolve_lines(clip_line_to_poly(timber_on_raw, on))

### 6.3 Combine QC and ON timber lines -------------------------------------------

timber_parts <- Filter(Negate(is.null), list(timber_qc, timber_on))

if (length(timber_parts) == 0) {
  stop("Timber allocation line could not be created after clipping to Quebec and Ontario.")
}

timber_limit <- do.call(rbind, timber_parts)
timber_limit <- dissolve_lines(timber_limit)

## 7. LOAD AND PREPARE HERD POLYGONS =============================================

### 7.1 Read Quebec herd polygons ------------------------------------------------

qc_herds <- safe_read_sf(qc_shp, crs_map) |>
  standardize_herd_name("QC") |>
  mutate(
    province   = "QC",
    herd_label = make_herd_label(herd_name)
  ) |>
  select(province, herd_name, herd_label, geometry)

### 7.2 Read Ontario herd polygons -----------------------------------------------

on_herds <- safe_read_sf(on_shp, crs_map) |>
  standardize_herd_name("ON") |>
  mutate(
    province   = "ON",
    herd_label = make_herd_label(herd_name)
  ) |>
  select(province, herd_name, herd_label, geometry)

### 7.3 Combine and dissolve herd polygons ---------------------------------------

herds <- bind_rows(qc_herds, on_herds) |>
  st_make_valid() |>
  dissolve_by_herd()

### 7.4 Build robust herd-code lookup --------------------------------------------

code_lookup <- tibble(
  map_label = names(herd_code_map),
  herd_code = unname(herd_code_map),
  label_key = normalize_key(names(herd_code_map))
) |>
  distinct(label_key, .keep_all = TRUE)

### 7.5 Add herd codes to the herd polygons --------------------------------------

herds <- herds |>
  mutate(label_key = normalize_key(herd_label)) |>
  left_join(code_lookup |> select(label_key, herd_code), by = "label_key") |>
  mutate(herd_code = if_else(is.na(herd_code), herd_label, herd_code)) |>
  select(-label_key)

## 8. BUILD FOOTPRINT FROM RASTER TEMPLATE =======================================

### 8.1 Convert raster valid area to a polygon footprint -------------------------

foot_mask <- terra::ifel(!is.na(r_template[[1]]), 1, NA)

if (foot_fact > 1) {
  foot_mask <- terra::aggregate(foot_mask, fact = foot_fact, fun = "max", na.rm = TRUE)
}

foot_poly <- terra::as.polygons(foot_mask, dissolve = TRUE, na.rm = TRUE)
foot_poly_sf <- sf::st_as_sf(foot_poly) |>
  sf::st_make_valid()

foot_poly_sf <- nngeo::st_remove_holes(foot_poly_sf)

footprint_poly <- sf::st_union(foot_poly_sf)

footprint_sf <- sf::st_sf(
  id = 1,
  geometry = sf::st_sfc(footprint_poly, crs = sf::st_crs(foot_poly_sf))
) |>
  sf::st_make_valid()

r_outline <- sf::st_boundary(footprint_sf)

### 8.2 Build the final map extent -----------------------------------------------

roi_geom <- sf::st_buffer(footprint_sf, dist = frame_buffer_m)
bb <- sf::st_bbox(roi_geom)

## 9. BUILD LABEL POINTS AND HERD COLOURS ========================================

### 9.1 Build herd label points --------------------------------------------------

label_pts_sf <- suppressWarnings(
  st_centroid(herds, of_largest_polygon = TRUE)
)

label_xy <- st_coordinates(label_pts_sf)

label_df <- label_pts_sf |>
  st_drop_geometry() |>
  mutate(
    X = label_xy[, 1],
    Y = label_xy[, 2]
  )

### 9.2 Build herd colour palette ------------------------------------------------

herd_levels <- herds$herd_label |> unique() |> sort()

dark2_base <- brewer.pal(8, "Dark2")
herd_cols <- colorRampPalette(dark2_base)(length(herd_levels))
names(herd_cols) <- herd_levels

## 10. BUILD HERD-CODE KEY TABLES ================================================

### 10.1 Build ON and QC key tables ----------------------------------------------

qc_key_tbl <- make_key_table(qc_herds, code_lookup, pretty_name_map)
on_key_tbl <- make_key_table(on_herds, code_lookup, pretty_name_map)

### 10.2 Split each key into two columns -----------------------------------------

on_split <- split_two_cols(on_key_tbl$key_line)
qc_split <- split_two_cols(qc_key_tbl$key_line)

on_key_text_1 <- paste(on_split$col1, collapse = "\n")
on_key_text_2 <- paste(on_split$col2, collapse = "\n")
qc_key_text_1 <- paste(qc_split$col1, collapse = "\n")
qc_key_text_2 <- paste(qc_split$col2, collapse = "\n")

## 11. BUILD MAIN MAP ============================================================

### 11.1 Build the single-panel herd map -----------------------------------------

p_main <- ggplot() +
  
  # Basemap
  geom_sf(data = ocean,     fill = "grey85", color = NA) +
  geom_sf(data = countries, fill = "grey92", color = NA) +
  geom_sf(data = lakes,     fill = "grey85", color = NA) +
  geom_sf(data = countries, fill = NA, color = "grey35", linewidth = lw_country) +
  
  # Herd polygons
  geom_sf(
    data = herds,
    aes(fill = herd_label, color = herd_label),
    linewidth = lw_herd,
    show.legend = FALSE
  ) +
  
  # Province lines over polygons
  geom_sf(data = provinces, fill = NA, color = "grey25", linewidth = lw_prov_over) +
  geom_sf(data = qc_on_border, color = "grey10", linewidth = lw_qcon, show.legend = FALSE) +
  
  # Study area outline
  geom_sf(
    data = r_outline,
    aes(linetype = "Study Area"),
    color = "black",
    linewidth = lw_outline,
    fill = NA,
    show.legend = TRUE,
    key_glyph = "path"
  ) +
  
  # Timber allocation northern limit
  geom_sf(
    data = timber_limit,
    aes(linetype = "Timber Allocation Northern Limit"),
    color = "black",
    linewidth = lw_timber,
    fill = NA,
    show.legend = TRUE,
    key_glyph = "path"
  ) +
  
  # Herd labels
  ggrepel::geom_text_repel(
    data = label_df,
    aes(x = X, y = Y, label = herd_code),
    seed = 123,
    size = label_size,
    fontface = "plain",
    colour = "black",
    force = 0,
    force_pull = 0,
    box.padding = 0.02,
    point.padding = 0,
    min.segment.length = Inf,
    segment.color = NA,
    max.overlaps = Inf,
    bg.color = alpha("white", 0.90),
    bg.r = 0.10
  ) +
  
  coord_sf(
    crs = st_crs(crs_map),
    datum = st_crs(4326),
    xlim = c(bb["xmin"], bb["xmax"]),
    ylim = c(bb["ymin"], bb["ymax"] + y_pad_m),
    expand = FALSE
  ) +
  
  ggspatial::annotation_scale(
    location = "bl",
    width_hint = 0.16,
    pad_x = unit(0.20, "cm"),
    pad_y = unit(0.20, "cm"),
    text_cex = 0.7
  ) +
  
  scale_fill_manual(
    values = alpha(herd_cols, herd_fill_alpha),
    guide = "none"
  ) +
  scale_colour_manual(
    values = alpha(herd_cols, herd_line_alpha),
    guide = "none"
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c(
      "Study Area" = "dashed",
      "Timber Allocation Northern Limit" = "solid"
    ),
    breaks = c("Study Area", "Timber Allocation Northern Limit")
  ) +
  
  guides(
    linetype = guide_legend(
      order = 1,
      override.aes = list(
        colour = "black",
        linewidth = c(lw_outline, lw_timber)
      ),
      keywidth = unit(1.8, "cm"),
      keyheight = unit(0.35, "cm")
    )
  ) +
  
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major = element_line(color = "grey75", linewidth = 0.20),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    
    legend.position = "bottom",
    legend.box = "vertical",
    
    legend.key = element_rect(fill = NA, colour = NA),
    legend.background = element_rect(fill = NA, colour = NA),
    legend.box.background = element_blank(),
    
    legend.spacing.y = unit(0.01, "cm"),
    legend.margin = margin(t = -8, r = 2, b = 0, l = 2),
    legend.box.margin = margin(t = -10, r = 0, b = 0, l = 0),
    
    legend.text = element_text(size = 9),
    legend.title = element_blank(),
    
    plot.margin = margin(t = 4, r = 8, b = -8, l = 8)
  )

## 12. BUILD HERD-KEY PANEL ======================================================

### 12.1 Build the text-only herd-code panel -------------------------------------

p_key <- ggplot() +
  annotate(
    "text", x = 0.02, y = 0.99,
    label = on_key_text_1,
    hjust = 0, vjust = 1,
    size = key_text_size,
    fontface = "plain",
    lineheight = key_lineheight
  ) +
  annotate(
    "text", x = 0.27, y = 0.99,
    label = on_key_text_2,
    hjust = 0, vjust = 1,
    size = key_text_size,
    fontface = "plain",
    lineheight = key_lineheight
  ) +
  annotate(
    "text", x = 0.55, y = 0.99,
    label = qc_key_text_1,
    hjust = 0, vjust = 1,
    size = key_text_size,
    fontface = "plain",
    lineheight = key_lineheight
  ) +
  annotate(
    "text", x = 0.79, y = 0.99,
    label = qc_key_text_2,
    hjust = 0, vjust = 1,
    size = key_text_size,
    fontface = "plain",
    lineheight = key_lineheight
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(t = 4, r = 8, b = 2, l = 8)
  )

## 13. COMBINE AND EXPORT ========================================================

### 13.1 Combine map and key panel -----------------------------------------------

final_plot <- p_main / p_key +
  plot_layout(heights = unit(c(1, 7), c("null", "cm")))

### 13.2 Export final figure -----------------------------------------------------

ggsave(
  filename = out_png,
  plot = final_plot,
  width = 8.5,
  height = 8.2,
  dpi = 600
)


