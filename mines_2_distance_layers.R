# mines_2_distance_layers.R ----------

# Description:
# This script creates mine-distance rasters, residual suitability rasters,
# and mine-density rasters for the 1985, 2000, and 2020 mine snapshots.
#
# Required inputs:
# - data/mines/processed/1985_mines.csv
# - data/mines/processed/2000_mines.csv
# - data/mines/processed/2020_mines.csv
# - data/roads/paved/paved_roads_zoi_residual_suitability.tif (as template)
#
# Expected outputs:
# - data/mines/processed/distance_layers/dist2mine_masked_1985.tif
# - data/mines/processed/distance_layers/dist2mine_masked_2000.tif
# - data/mines/processed/distance_layers/dist2mine_masked_2020.tif
# - data/mines/processed/RS_layers/RS_masked_1985.tif
# - data/mines/processed/RS_layers/RS_masked_2000.tif
# - data/mines/processed/RS_layers/RS_masked_2020.tif
# - data/mines/processed/density_layers/mine_density_1985.tif
# - data/mines/processed/density_layers/mine_density_2000.tif
# - data/mines/processed/density_layers/mine_density_2020.tif

## 1. SETUP ======================================================================

### 1.1 Define project root ------------------------------------------------------

# This makes the script portable so the packaged project can be stored anywhere.
# The script tries to identify the folder that contains the "scripts" and "data"
# directories, then uses that as the project root for all file paths.

get_project_root <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  
  # Case 1: script run from file path
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(dirname(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))))
  }
  
  # Case 2: script run interactively in RStudio
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    editor_path <- tryCatch(
      rstudioapi::getSourceEditorContext()$path,
      error = function(e) ""
    )
    
    if (nzchar(editor_path)) {
      return(dirname(dirname(normalizePath(editor_path, winslash = "/", mustWork = TRUE))))
    }
  }
  
  # Case 3: fallback if path is available through sys.frames()
  script_path <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  
  if (!is.null(script_path)) {
    return(dirname(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))))
  }
  
  # Case 4: final fallback to current working directory if it looks like project root
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

### 1.2 Define scratch directory -------------------------------------------------

# Set a dedicated temporary folder so both terra and parallel workers
# use the same scratch location inside the packaged project (or an existing TMPDIR
# if one is already defined for the session).

scratch <- Sys.getenv("TMPDIR", unset = file.path(project_root, "temp"))
scratch <- normalizePath(scratch, winslash = "/", mustWork = FALSE)

if (!dir.exists(scratch)) {
  dir.create(scratch, recursive = TRUE)
}

Sys.setenv(TMPDIR = scratch, TMP = scratch, TEMP = scratch)

### 1.3 Load required packages ---------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(future.apply)
})

# Keep terra's own temporary-file location aligned with the scratch directory above
terraOptions(tempdir = scratch)

## 2. DEFINE HELPER FUNCTIONS ===================================================

### 2.1 Create distance and residual suitability layers --------------------------

# This function takes a mine point dataset and creates:
# 1) a raster of distance to the nearest mine, and
# 2) a residual suitability raster derived from that distance
#
# - converts mine coordinates to sf if needed
# - transforms points to the same CRS as the template raster
# - rasterizes mine locations
# - calculates distance to mines
# - crops/resamples/masks to match the template grid
# - applies the Leblond et al. 2014 residual suitability function

process_dist2mines <- function(mine_data,
                               buffered_template,
                               template_layer,
                               output_dist2mine,
                               output_RS) {
  
  # Convert tabular mine coordinates to sf points if needed
  if (!inherits(mine_data, "sf")) {
    mine_data <- st_as_sf(mine_data, coords = c("longitude", "latitude"), crs = 4326)
  }
  
  # Reproject the mine points so they line up with the raster grid
  mine_data <- st_transform(mine_data, terra::crs(template_layer))
  
  # Rasterize mine points using the buffered template.
  # Each mine cell receives a value of 1; all other cells are NA.
  mine_rast <- rasterize(vect(mine_data), buffered_template, field = 1, background = NA)
  
  # Compute distance to the nearest mine for all raster cells
  tmpfile <- tempfile(tmpdir = tempdir(), fileext = ".tif")
  dist_rast <- distance(mine_rast, filename = tmpfile, overwrite = TRUE)
  
  # Bring the raster back to the exact study grid used by the template layer
  dist_rast <- crop(dist_rast, template_layer)
  dist_rast <- resample(dist_rast, template_layer, method = "near")
  dist_rast <- mask(
    dist_rast,
    template_layer,
    filename = output_dist2mine,
    overwrite = TRUE
  )
  
  # Extend + crop is retained here to force an identical grid geometry
  dist_rast <- crop(extend(dist_rast, template_layer), template_layer)
  
  # Clamp distance at 5 km, then apply the residual suitability equation
  # I.e., distances beyond 5 km are truncated before computing suitability.
  dist_clamped <- clamp(dist_rast, 0, 5000, values = TRUE)
  RS <- (-2e-6 * dist_clamped^2) + (0.0279 * dist_clamped) + 5.9313
  
  # Save the residual suitability raster
  writeRaster(RS, output_RS, overwrite = TRUE, datatype = "FLT4S")
  
  invisible(RS)
}

### 2.2 Run focal density on one buffered tile ----------------------------------

# Note: The tile includes a buffer so the focal window is calculated correctly at 
# the tile edges. After the focal sum is calculated, only the tile core is kept.
#
# Output units:
# mines per km²

tile_focal_density <- function(tile_file,
                               win,
                               radius,
                               tile_core_ext,
                               scratch = tempdir()) {
  
  Sys.setenv(TMPDIR = scratch, TMP = scratch, TEMP = scratch)
  
  # Read one raster tile
  r <- rast(tile_file)
  
  # Calculate focal sum using the supplied circular window
  tmpfile <- tempfile(tmpdir = scratch, fileext = ".tif")
  tmp <- focal(
    r,
    w = win,
    fun = "sum",
    na.policy = "omit",
    filename = tmpfile,
    overwrite = TRUE
  )
  
  # Crop back to the tile core and convert counts to density (mines per km²)
  tmp_core <- crop(tmp, tile_core_ext) / (pi * (radius / 1000)^2)
  
  # Save the processed tile
  outfile <- file.path(
    dirname(tile_file),
    sub("\\.tif$", "_dens.tif", basename(tile_file))
  )
  
  writeRaster(tmp_core, outfile, datatype = "FLT4S", overwrite = TRUE)
  
  outfile
}

### 2.3 Build mine density raster in parallel -----------------------------------

# - rasterizes mine points to a count raster
# - splits the raster into overlapping tiles
# - runs a focal density calculation on each tile in parallel
# - mosaics the tiles back together
# - resamples/masks/crops to the final template grid

mine_density_parallel <- function(mine_data,
                                  template_layer,
                                  radius,
                                  tile_size,
                                  ncores = max(1, parallel::detectCores() - 1),
                                  final_file) {
  
  scratch <- tempdir()
  tprefix <- tempfile(tmpdir = scratch, fileext = ".tif")
  
  # Convert mine coordinates to sf points if needed
  if (!inherits(mine_data, "sf")) {
    mine_data <- st_as_sf(mine_data, coords = c("longitude", "latitude"), crs = 4326)
  }
  
  # Reproject mine points to match the raster CRS
  mine_data <- st_transform(mine_data, terra::crs(template_layer))
  
  # Rasterize mine points to a count raster
  count_rast <- rasterize(
    vect(mine_data),
    template_layer,
    field = 1,
    fun = "sum",
    background = 0
  )
  
  # Determine how many cells of buffer are needed around each tile so the focal
  # operation is correct at tile edges
  buf_cells <- ceiling(radius / res(count_rast)[1])
  
  # Create overlapping tiles and store their extents
  tile_files <- makeTiles(
    count_rast,
    y = tile_size,
    buffer = buf_cells,
    filename = tprefix,
    overwrite = TRUE
  )
  
  tile_ext <- getTileExtents(count_rast, y = tile_size, buffer = buf_cells)
  
  # Create a circular focal window and convert it to a binary mask
  win <- focalMat(count_rast, radius, "circle")
  win[win != 0] <- 1
  win[win == 0] <- NA
  
  # Process tiles in parallel
  future::plan(future::multisession, workers = ncores)
  
  core_tiles <- future.apply::future_lapply(
    seq_along(tile_files),
    function(i) {
      tile_focal_density(
        tile_file = tile_files[i],
        win = win,
        radius = radius,
        tile_core_ext = ext(tile_ext[i, 1:4]),
        scratch = scratch
      )
    },
    future.seed = TRUE
  )
  
  # Stitch all processed tiles back together into one raster
  mosaic_rast <- mosaic(
    sprc(lapply(core_tiles, rast)),
    fun = "first",
    filename = tempfile(tmpdir = scratch, fileext = ".tif"),
    datatype = "FLT4S",
    overwrite = TRUE
  )
  
  # Bring the mosaic back to the exact template geometry
  mosaic_rast <- resample(mosaic_rast, template_layer, method = "near")
  mosaic_rast <- mask(
    mosaic_rast,
    template_layer,
    filename = tempfile(tmpdir = scratch, fileext = ".tif"),
    overwrite = TRUE
  )
  
  # Extend + crop retained to force identical final grid
  mosaic_rast <- crop(extend(mosaic_rast, template_layer), template_layer)
  
  mosaic_rast <- mask(
    mosaic_rast,
    template_layer,
    filename = tempfile(tmpdir = scratch, fileext = ".tif"),
    overwrite = TRUE
  )
  
  # Save final density raster
  writeRaster(
    mosaic_rast,
    final_file,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      NAflag = -9999
    )
  )
  
  invisible(mosaic_rast)
}

## 3. DEFINE PATHS ===============================================================

### 3.1 Input file paths --------------------------------------------------------

mines_1985_path <- file.path(project_root, "data", "mines", "processed", "1985_mines.csv")
mines_2000_path <- file.path(project_root, "data", "mines", "processed", "2000_mines.csv")
mines_2020_path <- file.path(project_root, "data", "mines", "processed", "2020_mines.csv")

road_layer_path <- file.path(
  project_root,
  "data", "roads", "paved", "paved_roads_zoi_residual_suitability.tif"
)

### 3.2 Output folder paths -----------------------------------------------------

distance_dir <- file.path(project_root, "data", "mines", "processed", "distance_layers")
rs_dir       <- file.path(project_root, "data", "mines", "processed", "RS_layers")
density_dir  <- file.path(project_root, "data", "mines", "processed", "density_layers")

# Create output folders if needed
dir.create(distance_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rs_dir,       recursive = TRUE, showWarnings = FALSE)
dir.create(density_dir,  recursive = TRUE, showWarnings = FALSE)

## 4. READ INPUT DATA ============================================================

### 4.1 Import mine snapshots and road raster -----------------------------------

# These three mine CSVs were created by mines_1_preprocessing.R
mines1985 <- read.csv(mines_1985_path)
mines2000 <- read.csv(mines_2000_path)
mines2020 <- read.csv(mines_2020_path)

# This road raster is used as the template grid for all outputs
road_layer <- rast(road_layer_path)

## 5. PREPARE TEMPLATE LAYERS ===================================================

### 5.1 Create buffered road template -------------------------------------------

# A buffer is added around the road raster extent before computing distance to
# mines. This helps avoid edge problems during the distance calculation.

extent_sf <- st_as_sf(as.polygons(ext(road_layer)))
buffer_extent <- st_buffer(extent_sf, dist = 5100) |> ext()
roads_buffered <- crop(road_layer, buffer_extent, extend = TRUE)

## 6. CREATE DISTANCE AND RESIDUAL SUITABILITY LAYERS ===========================

### 6.1 Run 1985 distance and residual suitability -------------------------------

process_dist2mines(
  mine_data = mines1985,
  buffered_template = roads_buffered,
  template_layer = road_layer,
  output_dist2mine = file.path(distance_dir, "dist2mine_masked_1985.tif"),
  output_RS = file.path(rs_dir, "RS_masked_1985.tif")
)

### 6.2 Run 2000 distance and residual suitability -------------------------------

process_dist2mines(
  mine_data = mines2000,
  buffered_template = roads_buffered,
  template_layer = road_layer,
  output_dist2mine = file.path(distance_dir, "dist2mine_masked_2000.tif"),
  output_RS = file.path(rs_dir, "RS_masked_2000.tif")
)

### 6.3 Run 2020 distance and residual suitability -------------------------------

process_dist2mines(
  mine_data = mines2020,
  buffered_template = roads_buffered,
  template_layer = road_layer,
  output_dist2mine = file.path(distance_dir, "dist2mine_masked_2020.tif"),
  output_RS = file.path(rs_dir, "RS_masked_2020.tif")
)

## 7. CREATE DENSITY LAYERS =====================================================

### 7.1 Define number of CPU cores to use ---------------------------------------

# Use all but one available core
ncores_to_use <- max(1, parallel::detectCores() - 1)

### 7.2 Run 1985 mine density ---------------------------------------------------

mine_density_parallel(
  mine_data = mines1985,
  template_layer = road_layer,
  radius = 1000,
  tile_size = 2000,
  ncores = ncores_to_use,
  final_file = file.path(density_dir, "mine_density_1985.tif")
)

### 7.3 Run 2000 mine density ---------------------------------------------------

mine_density_parallel(
  mine_data = mines2000,
  template_layer = road_layer,
  radius = 1000,
  tile_size = 2000,
  ncores = ncores_to_use,
  final_file = file.path(density_dir, "mine_density_2000.tif")
)

### 7.4 Run 2020 mine density ---------------------------------------------------

mine_density_parallel(
  mine_data = mines2020,
  template_layer = road_layer,
  radius = 1000,
  tile_size = 2000,
  ncores = ncores_to_use,
  final_file = file.path(density_dir, "mine_density_2020.tif")
)

