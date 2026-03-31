# habitat_3_create_composite.R ----------

# Description:
# Run this script after habitat_2_classification.R.
# It combines the raw binary habitat class tiles into final composite habitat tiles.
#
# Required inputs:
# - data/habitat/outputs/raw/*.tif created by habitat_2_classification.R
#
# Expected outputs:
# - data/habitat/outputs/compositeHabitat/compositeHabitat_1985_tile*.tif
# - data/habitat/outputs/compositeHabitat/compositeHabitat_2000_tile*.tif
# - data/habitat/outputs/compositeHabitat/compositeHabitat_2020_tile*.tif
#
# Composite class codes:
# 1 = natural disturbance
# 2 = harvest 0 to 5 years
# 3 = harvest 6 to 20 years
# 4 = young conifer
# 5 = mature conifer
# 6 = open woodland
# 7 = wetland
# 8 = regenerating stand

# NOTE: These output files are only used in the postprocessing stage (postprocess_4_habitat_proportions.R), 
#       not in the focal sum calculations (the next step, 'habitat_4_focal.R', uses binary tiles from previous step). 

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
  library(stringr)
  library(parallel)
})

### 1.3 User inputs --------------------------------------------------------------

# Choose which snapshot year(s) to run.
# Use any combination of: 1985, 2000, 2020.
run_years <- c(1985, 2000, 2020)

# Input and output folders.
# Keep these as-is unless the project structure was changed.
base_dir      <- "data/habitat/outputs"
raw_dir       <- file.path(base_dir, "raw")
composite_dir <- file.path(base_dir, "compositeHabitat")

# Number of CPU cores to use for parallel tile processing.
ncores <- 4

# Worker raster settings.
# Keep these as-is unless you intentionally want to change runtime behaviour.
worker_memfrac <- 0.25
worker_gdal_cachemax <- "300"
worker_omp_threads <- "1"
worker_num_threads <- 1

### 1.4 Create output folder -----------------------------------------------------

dir.create(composite_dir, recursive = TRUE, showWarnings = FALSE)

## 2. DEFINE HELPER FUNCTIONS ====================================================

### 2.1 Extract tile ID from a file name -----------------------------------------

# Returns values such as "tile1", "tile2", etc.
get_tile_num <- function(x) {
  stringr::str_extract(x, "tile\\d+")
}

### 2.2 Build one composite habitat tile -----------------------------------------

# This combines the raw class tiles for one tile and one year into a single output raster.

makeCompositeHabitat <- function(tileList, year, outputPath) {
  
  # Keep only files for the selected year.
  tileList <- tileList[grepl(as.character(year), tileList)]
  
  if (length(tileList) == 0) {
    stop("Incorrect year: no tiles matched the selected year.")
  }
  
  # Multiplies a binary raster by its composite class code.
  multi <- function(x, by) x * by
  
  # Temporary rasters are reused between steps.
  tempFile1 <- tempfile(fileext = ".tif")
  tempFile2 <- tempfile(fileext = ".tif")
  tempFile3 <- tempfile(fileext = ".tif")
  tempFile4 <- tempfile(fileext = ".tif")
  
  ### Step 1: fire and harvest ---------------------------------------------------
  
  fire <- terra::rast(tileList[grepl("naturalDisturbance", tileList)])
  
  hvY <- terra::app(
    x = terra::rast(tileList[grepl("harvest_0to5_", tileList)]),
    fun = multi,
    by = 2,
    filename = tempFile1
  )
  
  hvO <- terra::app(
    x = terra::rast(tileList[grepl("harvest_6to20_", tileList)]),
    fun = multi,
    by = 3,
    filename = tempFile2
  )
  
  comp <- sum(hvY, hvO, fire, filename = tempFile3, na.rm = TRUE)
  rm(hvY, hvO, fire)
  gc()
  
  ### Step 2: conifers -----------------------------------------------------------
  
  conY <- terra::app(
    x = terra::rast(tileList[grepl("youngConifer", tileList)]),
    fun = multi,
    by = 4,
    filename = tempFile1,
    overwrite = TRUE
  )
  
  conO <- terra::app(
    x = terra::rast(tileList[grepl("matureConifer", tileList)]),
    fun = multi,
    by = 5,
    filename = tempFile2,
    overwrite = TRUE
  )
  
  comp <- sum(conY, conO, comp, na.rm = TRUE, filename = tempFile4)
  rm(conY, conO)
  gc()
  
  ### Step 3: open woodland and wetland ------------------------------------------
  
  woodland <- terra::app(
    x = terra::rast(tileList[grepl("openWoodland", tileList)]),
    fun = multi,
    by = 6,
    filename = tempFile1,
    overwrite = TRUE
  )
  
  wet <- terra::app(
    x = terra::rast(tileList[grepl("wetland", tileList)]),
    fun = multi,
    by = 7,
    filename = tempFile2,
    overwrite = TRUE
  )
  
  comp <- sum(
    comp, woodland, wet,
    na.rm = TRUE,
    filename = tempFile3,
    overwrite = TRUE
  )
  rm(woodland, wet)
  gc()
  
  ### Step 4: regenerating stands and final write --------------------------------
  
  regen <- terra::app(
    x = terra::rast(tileList[grepl("regeneratingStand", tileList)]),
    fun = multi,
    by = 8,
    filename = tempFile1,
    overwrite = TRUE
  )
  
  tileNum <- stringr::str_extract(tileList[1], "tile\\d+")
  outFile <- file.path(outputPath, paste0("compositeHabitat_", year, "_", tileNum, ".tif"))
  
  comp <- sum(comp, regen, na.rm = TRUE, filename = outFile, overwrite = TRUE)
  rm(regen)
  gc()
  
  # Quick check for unexpected overlap between classes.
  mm <- terra::global(comp, "max", na.rm = TRUE)[1, 1]
  
  if (!is.na(mm) && mm > 8) {
    warning("Composite max > 8 in ", basename(outFile), " (possible class overlap).")
  }
  
  outFile
}

## 3. FIND AVAILABLE RAW TILES ===================================================

### 3.1 List all raw habitat class tiles -----------------------------------------

all_raw <- list.files(raw_dir, pattern = "\\.tif$", full.names = TRUE)

### 3.2 Extract available tile IDs -----------------------------------------------

tile_ids <- sort(unique(na.omit(get_tile_num(all_raw))))

## 4. SET UP PARALLEL PROCESSING =================================================

### 4.1 Start cluster ------------------------------------------------------------

cl <- parallel::makeCluster(ncores, outfile = "")
on.exit(parallel::stopCluster(cl), add = TRUE)

### 4.2 Load packages and settings on workers ------------------------------------

parallel::clusterEvalQ(cl, {
  library(terra)
  library(stringr)
  
  terraOptions(
    tempDir = tempdir(),
    memfrac = worker_memfrac,
    numThreads = worker_num_threads
  )
  
  Sys.setenv(
    OMP_NUM_THREADS = worker_omp_threads,
    GDAL_CACHEMAX = worker_gdal_cachemax
  )
})

## 5. BUILD COMPOSITE HABITAT TILES ==============================================

### 5.1 Run each selected year ---------------------------------------------------

for (yr in run_years) {
  message("\n==== YEAR ", yr, " ====")
  
  results <- parallel::parLapplyLB(
    cl,
    tile_ids,
    function(tile_tag, year, raw_dir, composite_dir, .makeCompositeHabitat) {
      
      # Find all raw class files for the current tile.
      tile_files <- list.files(
        raw_dir,
        pattern = paste0("_", tile_tag, "\\.tif$"),
        full.names = TRUE
      )
      
      if (!length(tile_files)) {
        return(sprintf("skip %s (no raw files for tile)", tile_tag))
      }
      
      # Keep processing even if one tile fails.
      tryCatch(
        .makeCompositeHabitat(
          tileList = tile_files,
          year = year,
          outputPath = composite_dir
        ),
        error = function(e) sprintf("error on %s: %s", tile_tag, conditionMessage(e))
      )
    },
    year = yr,
    raw_dir = raw_dir,
    composite_dir = composite_dir,
    .makeCompositeHabitat = makeCompositeHabitat
  )
  
  if (length(results)) {
    print(results)
  }
}

# NEXT STEP -< Continue to habitat_4_focal.R 

