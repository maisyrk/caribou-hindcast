# postprocess_5_herd_spatial_classification.R ----------

# Description:
# Run this script to classify each herd by its position relative to the
# Northern Limit of Timber Allocation using the heavily simplified split method.
#
# Method:
# - split Quebec using a heavily simplified QC timber line
# - split Ontario using a heavily simplified ON timber line
# - merge QC + ON south polygons and QC + ON north polygons
# - clip herds to the same QC + ON mainland polygons used for the split
# - calculate % south and % north for each herd
#
# Required inputs:
# - data/post_processing/herds/Aire_repartition_populations_caribouForestier.shp
# - data/post_processing/herds/Caribou_range_boundary.shp
# - map_data/forest_limit/LIM_FOR_ATT_2018.shp
# - map_data/forest_limit/northern_extent_managed_forest.shp
#
# Expected outputs:
# - tables/herd_timber_line_position_by_herd_combined_ON_QC.csv
#

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
    "Open the packaged project folder as the working directory, ",
    "or run this script from that project root."
  )
}

project_root <- get_project_root()
setwd(project_root)

### 1.2 Load required packages ---------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tibble)
  library(rnaturalearth)
  library(rnaturalearthdata)
})

### 1.3 User inputs --------------------------------------------------------------

# Herd shapefiles
qc_shp <- file.path(
  project_root,
  "data", "post_processing", "herds",
  "Aire_repartition_populations_caribouForestier.shp"
)

on_shp <- file.path(
  project_root,
  "data", "post_processing", "herds",
  "Caribou_range_boundary.shp"
)

# Timber allocation limits
timber_limit_qc_path <- file.path(
  project_root,
  "map_data", "forest_limit",
  "LIM_FOR_ATT_2018.shp"
)

timber_limit_on_path <- file.path(
  project_root,
  "map_data", "forest_limit",
  "northern_extent_managed_forest.shp"
)

# Projection
crs_map <- "EPSG:3979"

# Output
out_dir <- file.path(project_root, "tables")
out_csv_main <- file.path(out_dir, "herd_timber_line_position_by_herd_combined_ON_QC.csv")

# Tuning
simplify_tol_m   <- 15000
snap_tol_m       <- 10000
min_face_area_m2 <- 1e8

# Rounding
digits_area <- 2
digits_pct  <- 2

### 1.4 Create output folder -----------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## 2. DEFINE LOOKUPS =============================================================

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

## 3. DEFINE HELPER FUNCTIONS ====================================================

### 3.1 Drop Z and M safely ------------------------------------------------------

drop_zm_safe <- function(x) {
  tryCatch(st_zm(x, drop = TRUE, what = "ZM"), error = function(e) x)
}

### 3.2 Read and clean an sf file ------------------------------------------------

read_clean_sf <- function(path, crs_out) {
  st_read(path, quiet = TRUE) |>
    st_make_valid() |>
    drop_zm_safe() |>
    st_transform(crs_out)
}

### 3.3 Standardize herd-name field ----------------------------------------------

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
    x$herd_name <- as.character(x[[hit]])
  } else {
    x$herd_name <- paste0(province_tag, "_Herd_", seq_len(nrow(x)))
  }
  
  x$herd_name <- trimws(x$herd_name)
  bad <- is.na(x$herd_name) | x$herd_name == ""
  x$herd_name[bad] <- paste0(province_tag, "_Herd_", which(bad))
  
  x
}

### 3.4 Clean herd names for matching --------------------------------------------

make_herd_label <- function(x) {
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x <- gsub("['`]", "", x, perl = TRUE)
  x <- gsub("&", " And ", x, perl = TRUE)
  x <- gsub("[^A-Za-z0-9]+", " ", x, perl = TRUE)
  x <- trimws(x)
  x <- gsub("\\s+", " ", x, perl = TRUE)
  tools::toTitleCase(x)
}

### 3.5 Normalize a text key -----------------------------------------------------

normalize_key <- function(x) {
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "", x, perl = TRUE)
  x
}

### 3.6 Split one province with a simplified timber line -------------------------

split_one_province <- function(prov_poly,
                               timber_raw,
                               prov_name,
                               simplify_tol_m,
                               snap_tol_m,
                               min_face_area_m2) {
  
  prov_parts <- suppressWarnings(st_cast(st_make_valid(prov_poly), "POLYGON", warn = FALSE))
  prov_parts <- st_sf(geometry = st_geometry(prov_parts), crs = st_crs(prov_poly))
  prov_parts$area_m2 <- as.numeric(st_area(prov_parts))
  prov_main <- prov_parts[which.max(prov_parts$area_m2), "geometry"]
  
  timber_clip <- suppressWarnings(st_intersection(timber_raw, prov_main)) |>
    st_make_valid()
  
  timber_ls <- suppressWarnings(st_collection_extract(timber_clip, "LINESTRING", warn = FALSE))
  if (nrow(timber_ls) == 0) {
    timber_ml <- suppressWarnings(st_collection_extract(timber_clip, "MULTILINESTRING", warn = FALSE))
    timber_ls <- suppressWarnings(st_cast(timber_ml, "LINESTRING", warn = FALSE))
  }
  
  if (nrow(timber_ls) == 0) {
    stop("No timber line remained after clipping for ", prov_name)
  }
  
  timber_ls <- st_sf(geometry = st_geometry(timber_ls), crs = st_crs(prov_main))
  
  timber_union <- st_union(timber_ls)
  timber_simple <- st_simplify(
    timber_union,
    dTolerance = simplify_tol_m,
    preserveTopology = TRUE
  )
  
  timber_simple_sf <- st_sf(geometry = timber_simple, crs = st_crs(prov_main)) |>
    st_make_valid()
  
  timber_simple_ls <- suppressWarnings(st_collection_extract(timber_simple_sf, "LINESTRING", warn = FALSE))
  if (nrow(timber_simple_ls) == 0) {
    timber_simple_ml <- suppressWarnings(st_collection_extract(timber_simple_sf, "MULTILINESTRING", warn = FALSE))
    timber_simple_ls <- suppressWarnings(st_cast(timber_simple_ml, "LINESTRING", warn = FALSE))
  }
  
  if (nrow(timber_simple_ls) == 0) {
    stop("No usable simplified timber line for ", prov_name)
  }
  
  timber_simple_ls$len_m <- as.numeric(st_length(timber_simple_ls))
  timber_main <- timber_simple_ls[which.max(timber_simple_ls$len_m), "geometry"]
  
  prov_boundary <- st_boundary(prov_main)
  timber_main$geometry <- suppressWarnings(
    st_snap(st_geometry(timber_main), st_geometry(prov_boundary), tolerance = snap_tol_m)
  )
  timber_main <- st_make_valid(timber_main)
  
  endpoints <- st_boundary(timber_main)
  endpoints <- suppressWarnings(st_cast(endpoints, "POINT", warn = FALSE))
  endpoints <- st_sf(geometry = endpoints, crs = st_crs(prov_main))
  
  if (nrow(endpoints) < 2) {
    stop("Could not extract timber-line endpoints for ", prov_name)
  }
  
  connectors <- st_nearest_points(endpoints, prov_boundary)
  connectors <- st_sf(geometry = connectors, crs = st_crs(prov_main))
  
  linework <- c(
    st_geometry(prov_boundary),
    st_geometry(timber_main),
    st_geometry(connectors)
  )
  st_crs(linework) <- st_crs(prov_main)
  
  linework <- st_union(linework)
  
  if ("st_node" %in% getNamespaceExports("sf")) {
    linework <- sf::st_node(linework)
  }
  
  faces <- st_polygonize(linework)
  faces <- st_collection_extract(faces, "POLYGON", warn = FALSE)
  faces <- st_sf(geometry = faces, crs = st_crs(prov_main))
  faces <- suppressWarnings(st_intersection(faces, prov_main)) |>
    st_make_valid()
  
  if (nrow(faces) == 0) {
    stop("Polygonization produced no faces for ", prov_name)
  }
  
  faces$area_m2 <- as.numeric(st_area(faces))
  faces <- faces[faces$area_m2 > min_face_area_m2, ]
  
  if (nrow(faces) < 2) {
    stop("Failed to split ", prov_name, " into at least 2 usable faces")
  }
  
  faces <- faces[order(faces$area_m2, decreasing = TRUE), ]
  faces <- faces[1:2, ]
  faces$face_id <- 1:2
  
  face_pts <- st_point_on_surface(faces)
  y <- st_coordinates(face_pts)[, "Y"]
  
  south_poly <- faces[which.min(y), ] |>
    summarise(do_union = TRUE) |>
    st_make_valid()
  south_poly$side <- "south"
  south_poly$province_split <- prov_name
  
  north_poly <- faces[which.max(y), ] |>
    summarise(do_union = TRUE) |>
    st_make_valid()
  north_poly$side <- "north"
  north_poly$province_split <- prov_name
  
  list(
    province_main = prov_main,
    south_poly = south_poly,
    north_poly = north_poly
  )
}

## 4. LOAD PROVINCES =============================================================

### 4.1 Read Quebec and Ontario polygons -----------------------------------------

provinces <- ne_states(country = "canada", returnclass = "sf") |>
  filter(name_en %in% c("Ontario", "Quebec")) |>
  st_transform(crs_map)

qc <- provinces |> filter(name_en == "Quebec")
on <- provinces |> filter(name_en == "Ontario")

## 5. LOAD RAW TIMBER LINES ======================================================

### 5.1 Read QC and ON timber lines ----------------------------------------------

timber_qc_raw <- read_clean_sf(timber_limit_qc_path, crs_map)
timber_on_raw <- read_clean_sf(timber_limit_on_path, crs_map)

## 6. SPLIT QC AND ON SEPARATELY =================================================

### 6.1 Split Quebec -------------------------------------------------------------

qc_split <- split_one_province(
  prov_poly = qc,
  timber_raw = timber_qc_raw,
  prov_name = "QC",
  simplify_tol_m = simplify_tol_m,
  snap_tol_m = snap_tol_m,
  min_face_area_m2 = min_face_area_m2
)

### 6.2 Split Ontario ------------------------------------------------------------

on_split <- split_one_province(
  prov_poly = on,
  timber_raw = timber_on_raw,
  prov_name = "ON",
  simplify_tol_m = simplify_tol_m,
  snap_tol_m = snap_tol_m,
  min_face_area_m2 = min_face_area_m2
)

### 6.3 Build final north and south zones ----------------------------------------

south_zone <- rbind(qc_split$south_poly, on_split$south_poly) |>
  summarise(do_union = TRUE) |>
  st_make_valid()
south_zone$side <- "south"

north_zone <- rbind(qc_split$north_poly, on_split$north_poly) |>
  summarise(do_union = TRUE) |>
  st_make_valid()
north_zone$side <- "north"

### 6.4 Build combined mainland polygon ------------------------------------------

combined_main <- rbind(qc_split$province_main, on_split$province_main) |>
  summarise(do_union = TRUE) |>
  st_make_valid()

## 7. LOAD AND CLIP HERDS ========================================================

### 7.1 Read and prepare Quebec herds --------------------------------------------

qc_herds <- read_clean_sf(qc_shp, crs_map)
qc_herds <- standardize_herd_name(qc_herds, "QC")
qc_herds$province <- "QC"
qc_herds$herd_label <- make_herd_label(qc_herds$herd_name)
qc_herds <- qc_herds[, c("province", "herd_name", "herd_label", attr(qc_herds, "sf_column"))]
qc_herds <- qc_herds |>
  group_by(province, herd_name, herd_label) |>
  summarise(do_union = TRUE)

### 7.2 Read and prepare Ontario herds -------------------------------------------

on_herds <- read_clean_sf(on_shp, crs_map)
on_herds <- standardize_herd_name(on_herds, "ON")
on_herds$province <- "ON"
on_herds$herd_label <- make_herd_label(on_herds$herd_name)
on_herds <- on_herds[, c("province", "herd_name", "herd_label", attr(on_herds, "sf_column"))]
on_herds <- on_herds |>
  group_by(province, herd_name, herd_label) |>
  summarise(do_union = TRUE)

### 7.3 Combine and clip herds ---------------------------------------------------

herds <- rbind(qc_herds, on_herds) |>
  st_make_valid()

herds_clipped <- suppressWarnings(st_intersection(herds, combined_main)) |>
  st_make_valid() |>
  group_by(province, herd_name, herd_label) |>
  summarise(do_union = TRUE)

### 7.4 Add herd codes -----------------------------------------------------------

code_lookup <- tibble(
  herd_code = unname(herd_code_map),
  label_key = normalize_key(names(herd_code_map))
)

herds_clipped$label_key <- normalize_key(herds_clipped$herd_label)

herds_clipped <- left_join(herds_clipped, code_lookup, by = "label_key")
herds_clipped$herd_code[is.na(herds_clipped$herd_code)] <-
  herds_clipped$herd_label[is.na(herds_clipped$herd_code)]

herds_clipped$label_key <- NULL

## 8. CALCULATE NORTH/SOUTH AREAS ================================================

### 8.1 Build base results table -------------------------------------------------

results <- st_drop_geometry(herds_clipped)[, c("province", "herd_name", "herd_label", "herd_code")]
results$total_area_m2 <- as.numeric(st_area(herds_clipped))

### 8.2 Intersect with south zone ------------------------------------------------

south_int <- suppressWarnings(st_intersection(herds_clipped, south_zone))

if (nrow(south_int) > 0) {
  south_int$area_south_m2 <- as.numeric(st_area(south_int))
  south_tbl <- st_drop_geometry(south_int) |>
    group_by(province, herd_name, herd_label, herd_code) |>
    summarise(area_south_m2 = sum(area_south_m2), .groups = "drop")
} else {
  south_tbl <- results |>
    transmute(province, herd_name, herd_label, herd_code, area_south_m2 = 0)
}

### 8.3 Intersect with north zone ------------------------------------------------

north_int <- suppressWarnings(st_intersection(herds_clipped, north_zone))

if (nrow(north_int) > 0) {
  north_int$area_north_m2 <- as.numeric(st_area(north_int))
  north_tbl <- st_drop_geometry(north_int) |>
    group_by(province, herd_name, herd_label, herd_code) |>
    summarise(area_north_m2 = sum(area_north_m2), .groups = "drop")
} else {
  north_tbl <- results |>
    transmute(province, herd_name, herd_label, herd_code, area_north_m2 = 0)
}

### 8.4 Combine area summaries ---------------------------------------------------

results <- results |>
  left_join(south_tbl, by = c("province", "herd_name", "herd_label", "herd_code")) |>
  left_join(north_tbl, by = c("province", "herd_name", "herd_label", "herd_code"))

results$area_south_m2[is.na(results$area_south_m2)] <- 0
results$area_north_m2[is.na(results$area_north_m2)] <- 0

### 8.5 Calculate percentages and classes ----------------------------------------

results$pct_south <- 100 * results$area_south_m2 / results$total_area_m2
results$pct_north <- 100 * results$area_north_m2 / results$total_area_m2
results$pct_total_check <- results$pct_south + results$pct_north

results$dominant_side <- ifelse(
  results$pct_south > results$pct_north, "south",
  ifelse(results$pct_north > results$pct_south, "north", "tie")
)

results$position_class <- ifelse(
  results$pct_south > 75, "South (>75%)",
  ifelse(results$pct_north > 75, "North (>75%)", "Line straddler")
)

### 8.6 Convert to km² and round -------------------------------------------------

results$total_area_km2 <- results$total_area_m2 / 1e6
results$area_south_km2 <- results$area_south_m2 / 1e6
results$area_north_km2 <- results$area_north_m2 / 1e6

results <- results[, c(
  "province",
  "herd_code",
  "herd_name",
  "herd_label",
  "position_class",
  "dominant_side",
  "total_area_km2",
  "area_south_km2",
  "area_north_km2",
  "pct_south",
  "pct_north",
  "pct_total_check"
)]

results$position_class <- factor(
  results$position_class,
  levels = c("South (>75%)", "Line straddler", "North (>75%)")
)

results <- results |>
  arrange(province, position_class, herd_code, herd_label)

results$total_area_km2  <- round(results$total_area_km2, digits_area)
results$area_south_km2  <- round(results$area_south_km2, digits_area)
results$area_north_km2  <- round(results$area_north_km2, digits_area)
results$pct_south       <- round(results$pct_south, digits_pct)
results$pct_north       <- round(results$pct_north, digits_pct)
results$pct_total_check <- round(results$pct_total_check, digits_pct)

## 9. EXPORT OUTPUT ==============================================================
write.csv(results, out_csv_main, row.names = FALSE, na = "")

