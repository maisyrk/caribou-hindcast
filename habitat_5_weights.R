# habitat_5_weights.R ----------

# Description:
# Run this script after habitat_4_focal_v2.R.
# It combines the focal habitat proportion rasters into one weighted habitat suitability raster per year.
#
# Required inputs:
# - data/habitat/outputs/focalMosaic/prop/focal_<class>_<year>.tif
# - data/habitat/outputs/focalMosaic/denom/denom_year_<year>.tif
#
# Expected outputs:
# - data/habitat/outputs/focalMosaic/weighted/weighted_sum_<year>.tif
#
# These are the final weighted habitat layers that are used in the suitability model ('model_1_calculate_hsm.R')
#
# Notes:
# - Input focal rasters are expected to be scaled from 0 to 1000.
# - The output weighted raster stays on that same 0 to 1000 scale.

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
})

### 1.3 User inputs --------------------------------------------------------------

# Choose which snapshot year(s) to run.
# Use any combination of: 1985, 2000, 2020.
run_years <- c(1985, 2000, 2020)

# Input and output folders.
# Keep these as-is unless the project structure was changed.
base_dir  <- "data/habitat/outputs"
prop_dir  <- file.path(base_dir, "focalMosaic", "prop")
denom_dir <- file.path(base_dir, "focalMosaic", "denom")
out_dir   <- file.path(base_dir, "focalMosaic", "weighted")

# Scale used by the focal outputs.
# Keep this at 1000 unless the focal script output scale was changed.
scale_factor <- 1000L

# Class weights.
# Names must match the <class> portion of: focal_<class>_<year>.tif
weights <- c(
  harvest_0to5       = 0.04,
  harvest_6to20      = 0.04,
  matureConifer      = 0.25,
  naturalDisturbance = 0.06,
  openWoodland       = 0.22,
  regeneratingStand  = 0.06,
  youngConifer       = 0.19,
  wetland            = 0.14
)

# Normalize weights before combining?
# TRUE keeps the weighted sum on the same 0..1000 scale even if weights do not sum to 1.
normalize_weights <- TRUE

# terra memory and write settings.
# Keep these values unless you intentionally want to change runtime behaviour.
terra_memfrac <- 0.6

wopt_fast <- list(
  datatype = "INT2U",
  gdal = c(
    "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512",
    "COMPRESS=ZSTD", "NUM_THREADS=1", "BIGTIFF=YES"
  )
)

### 1.4 Create output folder -----------------------------------------------------

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

### 1.5 Apply terra settings -----------------------------------------------------

terraOptions(memfrac = terra_memfrac, progress = 1)

### 1.6 Normalize weights if requested -------------------------------------------

if (normalize_weights) {
  sw <- sum(weights)
  if (sw <= 0) {
    stop("Sum of weights must be > 0")
  }
  weights <- weights / sw
}

## 2. DEFINE HELPER FUNCTIONS ====================================================

### 2.1 Build expected focal file paths ------------------------------------------

prop_file_for <- function(cls, yr) {
  file.path(prop_dir, sprintf("focal_%s_%d.tif", cls, yr))
}

denom_file_for <- function(yr) {
  file.path(denom_dir, sprintf("denom_year_%d.tif", yr))
}

### 2.2 Create weighted raster for one year --------------------------------------

# This reads in the available focal class rasters for a year,
# applies the weights, and writes the final weighted sum (using denom layer as a mask/template).

make_weighted_for_year <- function(year) {
  message("=== Weighted combine for year ", year, " ===")
  
  # The denominator is used only to restore NA outside valid habitat.
  d_path <- denom_file_for(year)
  if (!file.exists(d_path)) {
    stop("Missing denominator for year ", year, ": ", d_path)
  }
  
  d <- rast(d_path)
  
  # Keep only class rasters that exist for this year.
  classes <- names(weights)
  paths   <- vapply(classes, prop_file_for, character(1L), yr = year)
  exists  <- file.exists(paths)
  
  if (!all(exists)) {
    missing <- classes[!exists]
    warning(
      "Missing class raster(s) for year ", year, ": ",
      paste(missing, collapse = ", ")
    )
  }
  
  classes   <- classes[exists]
  weights_y <- weights[exists]
  paths     <- paths[exists]
  
  if (!length(classes)) {
    stop("No class rasters found for year ", year, " matching provided weights.")
  }
  
  # Start the weighted sum accumulator.
  acc <- NULL
  template <- d
  
  for (i in seq_along(classes)) {
    cls <- classes[i]
    w   <- weights_y[i]
    fp  <- paths[i]
    
    r <- rast(fp)
    
    if (!compareGeom(r, template, stopOnError = FALSE)) {
      r <- resample(r, template, method = "near")
    }
    
    # Treat NA as 0 during the weighted sum.
    term <- classify(r, rbind(c(NA, NA, 0))) * w
    
    if (is.null(acc)) {
      acc <- term
    } else {
      acc <- acc + term
    }
    
    rm(r, term)
    gc()
  }
  
  # Input rasters are already scaled 0..1000, so the weighted sum stays on that scale.
  out_scaled <- round(acc)
  
  # Restore NA outside valid habitat.
  out_scaled <- mask(out_scaled, d, maskvalue = 0)
  
  out_path <- file.path(out_dir, sprintf("weighted_sum_%d.tif", year))
  writeRaster(
    out_scaled, out_path,
    overwrite = TRUE,
    datatype = "INT2U",
    gdal = wopt_fast$gdal
  )
  
  message("Wrote: ", out_path)
  
  # Quick range check in the console.
  stats <- try(global(out_scaled, c("min", "max"), na.rm = TRUE), silent = TRUE)
  if (!inherits(stats, "try-error")) {
    message(sprintf("Range (scaled): min=%s max=%s", stats[1, 1], stats[2, 1]))
  }
  
  invisible(out_path)
}

## 3. RUN WEIGHTED COMBINE =======================================================

### 3.1 Build weighted rasters for selected years --------------------------------

outs <- lapply(run_years, make_weighted_for_year)

message("Done. Outputs in: ", out_dir)

# NEXT STEP -> Continue to model_1_calculate_hsm.R
