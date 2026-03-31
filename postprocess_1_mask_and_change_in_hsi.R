# postprocess_1_mask_and_change_in_hsi.R ----------

# Description:
# Run this script after model_1_calculate_hsm.R.
# It first creates the masked HSI rasters needed for comparison, then calculates
# change-in-HSI rasters between the selected years.
#
# Required inputs:
# - HSI_output/HSI_1985.tif
# - HSI_output/HSI_2000.tif
# - HSI_output/HSI_2020.tif
#
# Expected outputs:
# - HSI_output/HSI_1985_masked.tif
# - HSI_output/HSI_2000_masked.tif
# - HSI_output/HSI_2020_masked.tif
# - HSI_output/dHSI_85-20.tif
# - HSI_output/dHSI_85-00.tif
# - HSI_output/dHSI_00-20.tif
#
# Notes:
# - 1985 is used as the reference mask.
# - Values greater than 0 but smaller than the threshold are set to NA in the
#   1985 raster before it is used as the mask template for later years.

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
})

### 1.3 User inputs --------------------------------------------------------------

# Choose which snapshot year(s) to run.
# Use any combination of: 1985, 2000, 2020.
run_years <- c(1985, 2000, 2020)

# Reference year used to build the valid mask.
# Keep this as 1985 to match the current workflow.
reference_year <- 1985

# Values greater than 0 but smaller than this threshold are set to NA
# in the reference HSI raster before masking the later years.
tiny_value_threshold <- 0.0001

### 1.4 Create output folder -----------------------------------------------------

dir.create("HSI_output", recursive = TRUE, showWarnings = FALSE)

## 2. DEFINE HELPERS =============================================================

### 2.1 Build HSI file paths -----------------------------------------------------

hsi_file_for <- function(year) {
  file.path("HSI_output", paste0("HSI_", year, ".tif"))
}

hsi_masked_file_for <- function(year) {
  file.path("HSI_output", paste0("HSI_", year, "_masked.tif"))
}

### 2.2 Check that required HSI rasters exist ------------------------------------

required_hsi <- vapply(run_years, hsi_file_for, character(1))

missing_hsi <- required_hsi[!file.exists(required_hsi)]

if (length(missing_hsi) > 0) {
  stop(
    "The following required HSI rasters were not found:\n",
    paste(missing_hsi, collapse = "\n")
  )
}

### 2.3 Format year suffix for output names --------------------------------------

year_suffix <- function(year) {
  substr(as.character(year), 3, 4)
}

## 3. CREATE MASKED HSI RASTERS ==================================================

### 3.1 Create the reference masked raster ---------------------------------------

# This applies the threshold only to the reference year.
# Zero values are left unchanged. Only values > 0 and < threshold become NA.

ref_raw_path    <- hsi_file_for(reference_year)
ref_masked_path <- hsi_masked_file_for(reference_year)

hsi_ref <- rast(ref_raw_path)

hsi_ref_masked <- hsi_ref
hsi_ref_masked[hsi_ref_masked > 0 & hsi_ref_masked < tiny_value_threshold] <- NA

writeRaster(
  hsi_ref_masked,
  filename  = ref_masked_path,
  datatype  = "FLT4S",
  overwrite = TRUE,
  wopt      = list(gdal = c("COMPRESS=LZW", "BIGTIFF=YES"))
)

### 3.2 Mask the other selected years to the reference raster --------------------

# Later years are masked to the 1985 masked raster so all comparisons use the
# same valid map extent.

other_years <- setdiff(run_years, reference_year)

for (year_i in other_years) {
  
  hsi_i <- rast(hsi_file_for(year_i))
  
  if (!compareGeom(hsi_ref_masked, hsi_i, stopOnError = FALSE)) {
    stop(
      "Raster geometry does not match between:\n",
      ref_masked_path, "\n",
      hsi_file_for(year_i)
    )
  }
  
  hsi_i_masked <- mask(hsi_i, hsi_ref_masked)
  
  writeRaster(
    hsi_i_masked,
    filename  = hsi_masked_file_for(year_i),
    datatype  = "FLT4S",
    overwrite = TRUE,
    wopt      = list(gdal = c("COMPRESS=LZW", "BIGTIFF=YES"))
  )
}

## 4. CALCULATE CHANGE IN HSI ====================================================

### 4.1 Load masked rasters needed for differencing ------------------------------

masked_hsi <- lapply(run_years, function(y) rast(hsi_masked_file_for(y)))
names(masked_hsi) <- as.character(run_years)

### 4.2 Check geometry across masked rasters -------------------------------------

if (all(c("1985", "2000") %in% names(masked_hsi))) {
  stopifnot(compareGeom(masked_hsi[["1985"]], masked_hsi[["2000"]], stopOnError = FALSE))
}

if (all(c("1985", "2020") %in% names(masked_hsi))) {
  stopifnot(compareGeom(masked_hsi[["1985"]], masked_hsi[["2020"]], stopOnError = FALSE))
}

### 4.3 Write difference rasters -------------------------------------------------

# Each output is written only if both years needed for that comparison are included.

if (all(c("1985", "2020") %in% names(masked_hsi))) {
  writeRaster(
    masked_hsi[["2020"]] - masked_hsi[["1985"]],
    filename  = file.path("HSI_output", "dHSI_85-20.tif"),
    datatype  = "FLT4S",
    overwrite = TRUE
  )
}

if (all(c("1985", "2000") %in% names(masked_hsi))) {
  writeRaster(
    masked_hsi[["2000"]] - masked_hsi[["1985"]],
    filename  = file.path("HSI_output", "dHSI_85-00.tif"),
    datatype  = "FLT4S",
    overwrite = TRUE
  )
}

if (all(c("2000", "2020") %in% names(masked_hsi))) {
  writeRaster(
    masked_hsi[["2020"]] - masked_hsi[["2000"]],
    filename  = file.path("HSI_output", "dHSI_00-20.tif"),
    datatype  = "FLT4S",
    overwrite = TRUE
  )
}

# NEXT STEP -> Continue to postprocessing_1_mask_and_change_in_hsi.R
