# postprocess_2_map_visualization.R ----------

# Description:
# Run this script to create publication-style 3-panel map figures for:
# - masked HSI rasters
# - dHSI rasters
#
# The script uses the same map layout for both figure types, including:
# - QC + ON basemap
# - raster-derived study area footprint
# - combined QC + ON timber allocation northern limit
# - 3 stacked panels
# - no inset
#
# Required inputs:
# - HSI_output/HSI_1985_masked.tif
# - HSI_output/HSI_2000_masked.tif
# - HSI_output/HSI_2020_masked.tif
# - HSI_output/dHSI_85-20.tif
# - HSI_output/dHSI_85-00.tif
# - HSI_output/dHSI_00-20.tif
# - map_data/forest_limit/LIM_FOR_ATT_2018.shp
# - map_data/forest_limit/northern_extent_managed_forest.shp
# - map_data/natural_earth/ne_50m_ocean/ne_50m_ocean.shp
# - map_data/natural_earth/ne_50m_lakes/ne_50m_lakes.shp
#
# Expected outputs:
# - map_layouts/hsi_3panel_final.png
# - map_layouts/dhsi_3panel_no_inset.png

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
  library(tidyterra)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggspatial)
})

if (!requireNamespace("nngeo", quietly = TRUE)) {
  stop(
    "Package 'nngeo' is required. Install it before running this script:\n",
    "install.packages('nngeo')"
  )
}

### 1.3 User inputs --------------------------------------------------------------

# Choose which figure type(s) to create.
# Use any combination of: "HSI", "dHSI"
run_map_types <- c("HSI", "dHSI")

# Masked HSI rasters
r_path_1985 <- "HSI_output/HSI_1985_masked.tif"
r_path_2000 <- "HSI_output/HSI_2000_masked.tif"
r_path_2020 <- "HSI_output/HSI_2020_masked.tif"

# dHSI rasters
r_path_85_20 <- "HSI_output/dHSI_85-20.tif"
r_path_85_00 <- "HSI_output/dHSI_85-00.tif"
r_path_00_20 <- "HSI_output/dHSI_00-20.tif"

# Timber allocation limits (QC + ON)
timber_limit_qc_path <- "map_data/forest_limit/LIM_FOR_ATT_2018.shp"
timber_limit_on_path <- "map_data/forest_limit/northern_extent_managed_forest.shp"

# Local Natural Earth layers
ocean_path <- "map_data/natural_earth/ne_50m_ocean/ne_50m_ocean.shp"
lakes_path <- "map_data/natural_earth/ne_50m_lakes/ne_50m_lakes.shp"

# Output folder
output_dir <- "map_layouts"

# Target plotting density per raster layer.
# Increase for more detail, decrease for faster plotting.
target_cells_per_layer <- 1e6

# Projection used for the main maps.
# Keep as EPSG:3979 to match the current raster workflow.
crs_map <- "EPSG:3979"

# Raster color ramps
hsi_raster_cols  <- c("darkred", "lightgoldenrodyellow", "dodgerblue4")
dhsi_raster_cols <- c("darkred", "grey92", "dodgerblue4")

### 1.4 Validate map type inputs -------------------------------------------------

allowed_map_types <- c("HSI", "dHSI")

if (!all(run_map_types %in% allowed_map_types)) {
  stop(
    "run_map_types must only contain: ",
    paste(allowed_map_types, collapse = ", ")
  )
}

### 1.5 Create output folder -----------------------------------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

## 2. DEFINE HELPER FUNCTIONS ====================================================

### 2.1 Downsample raster for plotting -------------------------------------------

downsample_for_plot <- function(r, target_cells = 1e6) {
  n <- ncell(r)
  
  if (is.na(n) || n <= target_cells) {
    return(r)
  }
  
  fact <- ceiling(sqrt(n / target_cells))
  terra::aggregate(r, fact = fact, fun = mean, na.rm = TRUE)
}

### 2.2 Safely read sf file ------------------------------------------------------

safe_read_sf <- function(path, crs_out) {
  if (is.null(path) || !nzchar(path)) {
    return(NULL)
  }
  
  st_read(path, quiet = TRUE) |>
    st_transform(crs_out)
}

### 2.3 Crop a line layer to a polygon -------------------------------------------

crop_line_to_poly <- function(x, poly) {
  if (is.null(x)) {
    return(NULL)
  }
  
  x <- x |>
    st_make_valid() |>
    st_zm(drop = TRUE, what = "ZM")
  
  x <- suppressWarnings(st_crop(x, st_bbox(poly)))
  x <- suppressWarnings(st_intersection(x, poly))
  x <- suppressWarnings(st_collection_extract(x, "LINESTRING", warn = FALSE))
  
  if (nrow(x) == 0) {
    return(NULL)
  }
  
  x
}

### 2.4 Dissolve line segments safely --------------------------------------------

dissolve_lines <- function(x) {
  if (is.null(x) || nrow(x) == 0) {
    return(NULL)
  }
  
  x <- x |>
    st_make_valid() |>
    st_zm(drop = TRUE, what = "ZM")
  
  g <- st_collection_extract(st_geometry(x), "LINESTRING", warn = FALSE)
  
  if (length(g) == 0) {
    return(NULL)
  }
  
  g <- st_union(g)
  g <- st_sfc(g, crs = st_crs(x))
  g <- st_collection_extract(g, "LINESTRING", warn = FALSE)
  
  if (length(g) == 0) {
    return(NULL)
  }
  
  g_ml <- st_cast(st_combine(g), "MULTILINESTRING")
  g_m  <- suppressWarnings(st_line_merge(g_ml))
  
  g_out <- st_collection_extract(g_m, "LINESTRING", warn = FALSE)
  
  if (length(g_out) == 0) {
    g_out <- st_collection_extract(g_ml, "LINESTRING", warn = FALSE)
  }
  
  st_sf(geometry = g_out)
}

### 2.5 Build raster footprint and extent ----------------------------------------

# This creates the irregular raster footprint used for the study-area outline
# and the buffered plotting extent.

build_footprint_objects <- function(r_plot, buffer_dist = 200000, foot_fact = 2) {
  foot_src <- r_plot[[1]]
  foot_mask <- terra::ifel(!is.na(foot_src), 1, NA)
  
  if (foot_fact > 1) {
    foot_mask <- terra::aggregate(foot_mask, fact = foot_fact, fun = "max", na.rm = TRUE)
  }
  
  foot_poly <- terra::as.polygons(foot_mask, dissolve = TRUE, values = FALSE)
  
  foot_poly_sf <- sf::st_as_sf(foot_poly) |>
    sf::st_make_valid()
  
  foot_poly_sf <- nngeo::st_remove_holes(foot_poly_sf)
  
  footprint_poly <- sf::st_union(foot_poly_sf)
  
  footprint_sf <- sf::st_sf(
    id = 1,
    geometry = sf::st_sfc(footprint_poly),
    crs = sf::st_crs(foot_poly_sf)
  ) |>
    sf::st_make_valid()
  
  r_outline <- sf::st_boundary(footprint_sf)
  roi_geom  <- sf::st_buffer(footprint_sf, dist = buffer_dist)
  bb        <- sf::st_bbox(roi_geom)
  
  list(
    footprint_sf = footprint_sf,
    r_outline    = r_outline,
    roi_geom     = roi_geom,
    bb           = bb
  )
}

### 2.6 Build timber allocation limit overlay ------------------------------------

build_timber_limit <- function(timber_qc_raw, timber_on_raw, qc, on) {
  timber_qc <- dissolve_lines(crop_line_to_poly(timber_qc_raw, qc))
  timber_on <- dissolve_lines(crop_line_to_poly(timber_on_raw, on))
  
  timber_limit <- bind_rows(timber_qc, timber_on)
  timber_limit <- dissolve_lines(timber_limit)
  
  timber_limit
}

### 2.7 Build HSI fill scale -----------------------------------------------------

build_hsi_fill_scale <- function(cols) {
  ggplot2::scale_fill_gradientn(
    colours  = scales::alpha(cols, 0.95),
    values   = scales::rescale(c(0, 0.5, 1)),
    limits   = c(0, 1),
    oob      = scales::squish,
    name     = "HSI",
    na.value = "transparent"
  )
}

### 2.8 Build dHSI fill scale ----------------------------------------------------

build_dhsi_fill_scale <- function(cols) {
  ggplot2::scale_fill_gradientn(
    colours  = scales::alpha(cols, 0.95),
    values   = scales::rescale(c(-1, 0, 1)),
    limits   = c(-1, 1),
    oob      = scales::squish,
    name     = "dHSI",
    na.value = "transparent"
  )
}

### 2.9 Build the final 3-panel map ----------------------------------------------

build_three_panel_map <- function(r_plot,
                                  provinces,
                                  countries,
                                  ocean,
                                  lakes,
                                  qc_on_border,
                                  r_outline,
                                  timber_limit,
                                  bb,
                                  crs_map,
                                  fill_scale) {
  
  lw_country   <- 0.20
  lw_prov_over <- 0.20
  lw_qcon      <- 0.20
  lw_timber    <- 0.22
  lw_outline   <- 0.25
  y_pad_m      <- 25000
  
  lt_scale <- scale_linetype_manual(
    name = NULL,
    breaks = c("Study Area", "Northern Limit for Timber Allocation"),
    values = c(
      "Study Area"                           = "dashed",
      "Northern Limit for Timber Allocation" = "solid"
    ),
    labels = c(
      "Study Area"                           = "Study Area",
      "Northern Limit for Timber Allocation" = "Timber Allocation Northern Limit"
    )
  )
  
  ggplot() +
    geom_sf(data = ocean, fill = "grey85", color = NA) +
    geom_sf(data = countries, fill = "grey92", color = NA) +
    geom_sf(data = lakes, fill = "grey85", color = NA) +
    geom_sf(data = countries, fill = NA, color = "grey35", linewidth = lw_country) +
    geom_spatraster(data = r_plot) +
    geom_sf(data = provinces, fill = NA, color = "grey25", linewidth = lw_prov_over) +
    geom_sf(data = qc_on_border, color = "grey10", linewidth = lw_qcon, show.legend = FALSE) +
    geom_sf(
      data = r_outline,
      fill = NA,
      aes(linetype = "Study Area"),
      color = "black",
      linewidth = lw_outline,
      key_glyph = "path"
    ) +
    {
      if (!is.null(timber_limit)) {
        geom_sf(
          data = timber_limit,
          aes(linetype = "Northern Limit for Timber Allocation"),
          color = "black",
          linewidth = lw_timber,
          key_glyph = "path"
        )
      }
    } +
    facet_wrap(~lyr, ncol = 1) +
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
      text_cex = 0.45
    ) +
    fill_scale +
    lt_scale +
    guides(
      fill = guide_colorbar(
        order = 1,
        direction = "horizontal",
        barwidth = unit(6.0, "cm"),
        barheight = unit(0.35, "cm"),
        title.position = "top",
        title.hjust = 0.5,
        label.position = "bottom"
      ),
      linetype = guide_legend(
        order = 2,
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          colour    = c("black", "black"),
          linewidth = c(lw_outline, lw_timber),
          fill      = NA
        )
      )
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_line(color = "grey75", linewidth = 0.20),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.title = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.spacing.y = unit(0.08, "cm"),
      legend.margin = margin(t = 2, r = 2, b = 2, l = 2),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9)
    )
}

## 3. LOAD SHARED BASEMAP DATA ===================================================

### 3.1 Load timber line source layers -------------------------------------------

timber_qc_raw <- safe_read_sf(timber_limit_qc_path, crs_map)
timber_on_raw <- safe_read_sf(timber_limit_on_path, crs_map)

### 3.2 Load local Natural Earth layers ------------------------------------------

ocean <- st_read(ocean_path, quiet = TRUE) |>
  st_transform(crs_map)

lakes <- st_read(lakes_path, quiet = TRUE) |>
  st_transform(crs_map)

### 3.3 Load country and province boundaries -------------------------------------

countries <- ne_countries(scale = "medium", returnclass = "sf") |>
  filter(admin %in% c("Canada", "United States of America")) |>
  st_transform(crs_map)

provinces <- ne_states(country = "canada", returnclass = "sf") |>
  filter(name_en %in% c("Ontario", "Quebec")) |>
  st_transform(crs_map)

qc <- provinces |> dplyr::filter(name_en == "Quebec")
on <- provinces |> dplyr::filter(name_en == "Ontario")

### 3.4 Build QC-ON border line --------------------------------------------------

qc_on_border <- st_intersection(st_boundary(qc), st_boundary(on)) |>
  st_collection_extract("LINESTRING")

### 3.5 Define initial processing extent -----------------------------------------

roi_geom0 <- st_union(provinces) |>
  st_buffer(dist = 200000)

### 3.6 Build timber allocation overlay ------------------------------------------

timber_limit <- build_timber_limit(
  timber_qc_raw = timber_qc_raw,
  timber_on_raw = timber_on_raw,
  qc = qc,
  on = on
)

## 4. RENDER HSI MAP ==============================================================

### 4.1 Load HSI rasters ---------------------------------------------------------

if ("HSI" %in% run_map_types) {
  required_hsi <- c(r_path_1985, r_path_2000, r_path_2020)
  
  missing_hsi <- required_hsi[!file.exists(required_hsi)]
  
  if (length(missing_hsi) > 0) {
    stop(
      "The following HSI raster(s) were not found:\n",
      paste(missing_hsi, collapse = "\n")
    )
  }
  
  r1985 <- rast(r_path_1985)
  r2000 <- rast(r_path_2000)
  r2020 <- rast(r_path_2020)
  
  r_stack_hsi <- c(r1985, r2000, r2020)
  names(r_stack_hsi) <- c("1985", "2000", "2020")
  
  ### 4.1.1 Crop and downsample for plotting -------------------------------------
  
  r_stack_hsi <- crop(r_stack_hsi, vect(roi_geom0))
  r_plot_hsi  <- downsample_for_plot(r_stack_hsi, target_cells = target_cells_per_layer)
  
  ### 4.1.2 Build footprint and map extent ---------------------------------------
  
  footprint_hsi <- build_footprint_objects(r_plot_hsi)
  
  ### 4.1.3 Build HSI map --------------------------------------------------------
  
  p_hsi <- build_three_panel_map(
    r_plot      = r_plot_hsi,
    provinces   = provinces,
    countries   = countries,
    ocean       = ocean,
    lakes       = lakes,
    qc_on_border = qc_on_border,
    r_outline   = footprint_hsi$r_outline,
    timber_limit = timber_limit,
    bb          = footprint_hsi$bb,
    crs_map     = crs_map,
    fill_scale  = build_hsi_fill_scale(hsi_raster_cols)
  )
  
  ### 4.1.4 Export HSI map -------------------------------------------------------
  
  ggsave(
    filename = file.path(output_dir, "hsi_3panel_final.png"),
    plot = p_hsi,
    width = 10,
    height = 7,
    dpi = 600
  )
}

## 5. RENDER dHSI MAP ============================================================

### 5.1 Load dHSI rasters --------------------------------------------------------

if ("dHSI" %in% run_map_types) {
  required_dhsi <- c(r_path_85_20, r_path_85_00, r_path_00_20)
  
  missing_dhsi <- required_dhsi[!file.exists(required_dhsi)]
  
  if (length(missing_dhsi) > 0) {
    stop(
      "The following dHSI raster(s) were not found:\n",
      paste(missing_dhsi, collapse = "\n")
    )
  }
  
  r85_20 <- rast(r_path_85_20)
  r85_00 <- rast(r_path_85_00)
  r00_20 <- rast(r_path_00_20)
  
  r_stack_dhsi <- c(r85_20, r85_00, r00_20)
  names(r_stack_dhsi) <- c("1985-2020", "1985-2000", "2000-2020")
  
  ### 5.1.1 Crop and downsample for plotting -------------------------------------
  
  r_stack_dhsi <- crop(r_stack_dhsi, vect(roi_geom0))
  r_plot_dhsi  <- downsample_for_plot(r_stack_dhsi, target_cells = target_cells_per_layer)
  
  ### 5.1.2 Build footprint and map extent ---------------------------------------
  
  footprint_dhsi <- build_footprint_objects(r_plot_dhsi)
  
  ### 5.1.3 Build dHSI map -------------------------------------------------------
  
  p_dhsi <- build_three_panel_map(
    r_plot       = r_plot_dhsi,
    provinces    = provinces,
    countries    = countries,
    ocean        = ocean,
    lakes        = lakes,
    qc_on_border = qc_on_border,
    r_outline    = footprint_dhsi$r_outline,
    timber_limit = timber_limit,
    bb           = footprint_dhsi$bb,
    crs_map      = crs_map,
    fill_scale   = build_dhsi_fill_scale(dhsi_raster_cols)
  )
  
  ### 5.1.4 Export dHSI map ------------------------------------------------------
  
  ggsave(
    filename = file.path(output_dir, "dhsi_3panel_no_inset.png"),
    plot = p_dhsi,
    width = 10,
    height = 7,
    dpi = 600
  )
}

