# model_1_calculate_hsm.R ----------

# Description:
# Run this script after the habitat, road, and mine steps are complete.
# It calculates the final habitat suitability map for one or more snapshot years.
#
# Required inputs:
# - data/mines/processed/density_layers/mine_density_<year>.tif
# - data/roads/paved/paved_road_density.tif
# - data/roads/unpaved/unpaved_road_density_<year>.tif
# - data/roads/unpaved/unpaved_roads_zoi_residual_suitability_<year>.tif
# - data/roads/paved/paved_roads_zoi_residual_suitability.tif
# - data/mines/processed/RS_layers/RS_masked_<year>.tif
# - data/habitat/outputs/focalMosaic/weighted/weighted_sum_<year>.tif
#
# Expected outputs:
# - data/human_infra_<year>.tif
# - HSI_output/HSI_<year>.tif
#
# Note: The paved road density and paved road residual-suitability layers are not year-specific
#       in the current workflow.

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

# terra memory settings.
# Keep these values unless you intentionally want to change runtime behaviour.
terra_memfrac <- 0.75
terra_progress <- 1

### 1.4 Create output folder -----------------------------------------------------

dir.create("HSI_output", recursive = TRUE, showWarnings = FALSE)

## 2. DEFINE HELPER FUNCTIONS ====================================================

### 2.1 Align one raster layer to another ----------------------------------------

# Checks whether the two raster grids already match.
# If not, it resamples the input layer to the template grid.

align_layers <- function(template_layer,
                         layer_to_align,
                         method  = "near",
                         outfile = NULL,
                         verbose = TRUE) {
  
  stopifnot(
    inherits(template_layer, "SpatRaster"),
    inherits(layer_to_align, "SpatRaster")
  )
  
  if (terra::compareGeom(layer_to_align, template_layer, stopOnError = FALSE)) {
    if (verbose) message("Grids already identical - no action needed.")
    aligned <- layer_to_align
    
  } else {
    if (verbose) message("Geometry differs - realigning ...")
    
    has_align <- "align" %in% names(formals(terra::resample))
    
    if (has_align) {
      aligned <- terra::resample(
        layer_to_align,
        template_layer,
        method = method,
        align  = TRUE,
        snap   = "out"
      )
    } else {
      padded  <- terra::extend(layer_to_align, template_layer)
      padded  <- terra::crop(padded, template_layer)
      aligned <- terra::resample(padded, template_layer, method = method)
    }
    
    if (verbose) {
      message(
        "   -> fixed: ",
        terra::compareGeom(aligned, template_layer, stopOnError = FALSE)
      )
    }
  }
  
  if (!is.null(outfile)) {
    terra::writeRaster(aligned, outfile, overwrite = TRUE)
    if (verbose) message("Written to: ", outfile)
  }
  
  aligned
}

### 2.2 Build file paths for one year --------------------------------------------

# This keeps all year-specific paths in one place.

build_paths_for_year <- function(year) {
  list(
    mines_density   = file.path("data", "mines", "processed", "density_layers", paste0("mine_density_", year, ".tif")),
    paved_density   = file.path("data", "roads", "paved", "paved_road_density.tif"),
    unpaved_density = file.path("data", "roads", "unpaved", paste0("unpaved_road_density_", year, ".tif")),
    
    dist_unpaved = file.path("data", "roads", "unpaved", paste0("unpaved_roads_zoi_residual_suitability_", year, ".tif")),
    dist_paved   = file.path("data", "roads", "paved", "paved_roads_zoi_residual_suitability.tif"),
    dist_mines   = file.path("data", "mines", "processed", "RS_layers", paste0("RS_masked_", year, ".tif")),
    
    habitat = file.path("data", "habitat", "outputs", "focalMosaic", "weighted", paste0("weighted_sum_", year, ".tif")),
    
    human_infra_out = file.path("data", paste0("human_infra_", year, ".tif")),
    hsi_out         = file.path("HSI_output", paste0("HSI_", year, ".tif"))
  )
}

### 2.3 Check that required files exist ------------------------------------------

check_required_files <- function(paths) {
  required_inputs <- c(
    paths$mines_density,
    paths$paved_density,
    paths$unpaved_density,
    paths$dist_unpaved,
    paths$dist_paved,
    paths$dist_mines,
    paths$habitat
  )
  
  missing_files <- required_inputs[!file.exists(required_inputs)]
  
  if (length(missing_files) > 0) {
    stop(
      "The following required input files were not found:\n",
      paste(missing_files, collapse = "\n")
    )
  }
}

## 3. RUN HSM MODEL ==============================================================

### 3.1 Process each selected year -----------------------------------------------

for (year_i in run_years) {
  
  message("\n============================================================")
  message("Running HSM model for year: ", year_i)
  message("============================================================")
  
  paths <- build_paths_for_year(year_i)
  check_required_files(paths)
  
  ## 3.1.1 Create human infrastructure layer =====================================
  
  ### Load infrastructure layers -------------------------------------------------
  
  mines   <- rast(paths$mines_density)
  paved   <- rast(paths$paved_density)
  unpaved <- rast(paths$unpaved_density)
  
  ### Align mine layer to paved road template ------------------------------------
  
  mines <- align_layers(
    template_layer = paved,
    layer_to_align = mines,
    outfile        = NULL
  )
  
  ### Calculate human infrastructure ---------------------------------------------
  
  # Same weighted sum as the uploaded script.
  human_infrastructure <- (0.35 * paved) + (0.30 * unpaved) + (0.35 * mines)
  
  ### Write human infrastructure layer -------------------------------------------
  
  writeRaster(
    human_infrastructure,
    paths$human_infra_out,
    overwrite = TRUE,
    datatype  = "FLT4S"
  )
  
  rm(mines, paved, unpaved)
  gc()
  
  ## 3.1.2 Load final equation inputs ============================================
  
  ### Load distance and infrastructure layers ------------------------------------
  
  dist_unpaved <- rast(paths$dist_unpaved)
  dist_paved   <- rast(paths$dist_paved)
  dist_mines   <- rast(paths$dist_mines)
  
  human_infrastructure <- rast(paths$human_infra_out)
  
  ### Load habitat layer ---------------------------------------------------------
  
  habitat <- rast(paths$habitat)
  
  ## 3.1.3 Align model inputs ====================================================
  
  ### Align mine residual-suitability layer --------------------------------------
  
  dist_mines <- align_layers(
    template_layer = dist_paved,
    layer_to_align = dist_mines,
    outfile        = NULL
  )
  
  ### Align habitat layer --------------------------------------------------------
  
  # Habitat is projected to the distance-layer grid using bilinear resampling.
  habitat_aligned <- project(
    habitat,
    dist_paved,
    method    = "bilinear",
    filename  = file.path(tempdir(), paste0("habitat_aligned_", year_i, ".tif")),
    overwrite = TRUE,
    gdal      = c("COMPRESS=LZW")
  )
  
  ### Align human infrastructure layer -------------------------------------------
  
  human_infrastructure <- align_layers(
    template_layer = dist_paved,
    layer_to_align = human_infrastructure,
    outfile        = NULL
  )
  
  ## 3.1.4 Calculate raw HSM =====================================================
  
  terraOptions(
    memfrac  = terra_memfrac,
    progress = terra_progress
  )
  
  # Same equation as the uploaded script.
  hsm_raw <- (0.55 * habitat_aligned) *
    (dist_unpaved * dist_paved * dist_mines) -
    (0.45 * human_infrastructure)
  
  ## 3.1.5 Standardize to 0-1 ====================================================
  
  mm      <- minmax(hsm_raw)
  min_val <- mm[1]
  max_val <- mm[2]
  
  # Standardize using the raster min and max for this year.
  hsm_final <- (hsm_raw - min_val) / (max_val - min_val)
  
  plot(hsm_final)
  
  ## 3.1.6 Write final HSI layer =================================================
  
  writeRaster(
    hsm_final,
    filename  = paths$hsi_out,
    datatype  = "FLT4S",
    overwrite = TRUE,
    wopt      = list(gdal = c("COMPRESS=LZW", "BIGTIFF=YES"))
  )
  
  ### Clean up before next year --------------------------------------------------
  
  rm(
    dist_unpaved, dist_paved, dist_mines,
    human_infrastructure, habitat, habitat_aligned,
    hsm_raw, hsm_final, mm, min_val, max_val, paths
  )
  gc()
}

# terra::purgeTmp()

# NEXT STEP -> Continue to postprocess_1_mask_and_change_in_hsi.R