# habitat_2_classification.R ----------

# Description:
# Run this script after habitat_1_preprocessing.R.
# It builds the raw binary habitat tiles used later to create the composite habitat maps.
#
# Required inputs:
# - data/habitat/GIS/tiles/*.tif created by habitat_1_preprocessing.R
# - tiled disturbance rasters
# - tiled SCANFI (land cover / canopy / deciduous / age inputs)
#
# Expected outputs:
# - data/habitat/outputs/raw/youngConifer*.tif
# - data/habitat/outputs/raw/matureConifer*.tif
# - data/habitat/outputs/raw/wetland*.tif
# - data/habitat/outputs/raw/naturalDisturbance*.tif
# - data/habitat/outputs/raw/harvest_0to5_*.tif
# - data/habitat/outputs/raw/harvest_6to20_*.tif
# - data/habitat/outputs/raw/openWoodland*.tif
# - data/habitat/outputs/raw/regeneratingStand*.tif

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
  library(data.table)
  library(stringr)
  library(reproducible)
  library(parallel)
})

### 1.3 User inputs --------------------------------------------------------------

# Choose which snapshot year(s) to run: use any combination of: 1985, 2000, 2020.
run_years <- c(1985, 2000, 2020)

# Main input and output folders.
# Keep these as-is unless the project structure was changed.
tile_dir <- "data/habitat/GIS/tiles"
raw_dir  <- "data/habitat/outputs/raw"
comp_dir <- "data/habitat/outputs/compositeHabitat"

# data.table thread setting.
dt_threads <- 4

# terra memory/progress settings.
# Keep these values unless you intentionally want different runtime behaviour.
terra_memfrac <- 0.5
terra_progress <- 1

# GDAL / thread settings used during raster work.
gdal_cachemax <- "1024"
gdal_num_threads <- "ALL_CPUS"
omp_num_threads <- "1"

### 1.4 Create output folders ----------------------------------------------------

checkPath(raw_dir,  create = TRUE)
checkPath(comp_dir, create = TRUE)

### 1.5 Apply global settings ----------------------------------------------------

setDTthreads(dt_threads)

terraOptions(
  memfrac = terra_memfrac,
  progress = terra_progress
)

Sys.setenv(
  GDAL_CACHEMAX    = gdal_cachemax,
  GDAL_NUM_THREADS = gdal_num_threads,
  OMP_NUM_THREADS  = omp_num_threads
)

## 2. DEFINE HELPER FUNCTIONS ====================================================

### 2.1 Name files by tile ID ----------------------------------------------------

# Builds a named vector using tile IDs such as tile1, tile2, etc.
# Stops if a tile is missing from a file name or if duplicate tile files are found.

named_by_tile <- function(files) {
  ids <- stringr::str_extract(files, "tile\\d+")
  
  if (anyNA(ids)) {
    stop(
      "Some files did not contain tile\\d+: ",
      paste(files[is.na(ids)], collapse = "\n")
    )
  }
  
  if (anyDuplicated(ids)) {
    dups <- unique(ids[duplicated(ids)])
    stop(
      "Duplicate tile files detected for: ",
      paste(dups, collapse = ", "),
      "\nFix by deleting extras or tightening list.files() patterns."
    )
  }
  
  setNames(files, ids)
}

### 2.2 Get the age tiles for one year -------------------------------------------

# Matches the tiled age rasters for the selected year.

age_map_by_year <- function(yr) {
  patt <- switch(
    as.character(yr),
    "1985" = "age1985_final.*\\.tif$",
    "2000" = "age2000_final.*\\.tif$",
    "2020" = "age2020_hybrid.*\\.tif$",
    paste0("age", yr, ".*\\.tif$")
  )
  
  named_by_tile(list.files(tile_dir, patt, full.names = TRUE))
}

### 2.3 Combine last disturbance year --------------------------------------------

# Builds a single last-disturbance-year raster for the chosen base year.
# Any disturbance year after the base year is ignored.

combine_last_disturbance_year <- function(baseYear, preFile = NULL, mainFile = NULL) {
  to_rast <- function(x) {
    if (is.null(x)) return(NULL)
    if (inherits(x, "SpatRaster")) return(x)
    if (is.character(x) && length(x) == 1 && nzchar(x) && file.exists(x)) {
      return(terra::rast(x))
    }
    NULL
  }
  
  pre  <- to_rast(preFile)
  main <- to_rast(mainFile)
  
  if (is.null(pre) && is.null(main)) return(NULL)
  if (is.null(pre)) {
    return(terra::ifel(
      main <= baseYear, main, NA,
      filename = tempfile(fileext = ".tif"),
      overwrite = TRUE
    ))
  }
  if (is.null(main)) {
    return(terra::ifel(
      pre <= baseYear, pre, NA,
      filename = tempfile(fileext = ".tif"),
      overwrite = TRUE
    ))
  }
  
  terra::compareGeom(pre, main)
  
  pre_c <- terra::ifel(
    pre <= baseYear, pre, NA,
    filename = tempfile(fileext = ".tif"),
    overwrite = TRUE
  )
  
  main_c <- terra::ifel(
    main <= baseYear, main, NA,
    filename = tempfile(fileext = ".tif"),
    overwrite = TRUE
  )
  
  terra::lapp(
    c(pre_c, main_c),
    fun = function(a, b) pmax(a, b, na.rm = TRUE),
    filename = tempfile(fileext = ".tif"),
    overwrite = TRUE
  )
}

## 3. DEFINE RAW TILE CLASSIFIERS ===============================================

### 3.1 Mature conifers ----------------------------------------------------------

# Writes two masks:
# - young conifer
# - mature conifer

MatureConifers <- function(ageFile, coverFile, decidFile,
                           distYearPreFile = NULL,
                           distYearMainFile = NULL,
                           baseYear) {
  
  id     <- stringr::str_extract(ageFile, "tile\\d+")
  msgHdr <- paste0("MatureConifers ", baseYear, " | ", id, " : ")
  
  age   <- terra::rast(ageFile)
  cover <- terra::rast(coverFile)
  decid <- terra::rast(decidFile)
  
  terra::compareGeom(decid, cover, age)
  message(msgHdr, "geometry check complete")
  
  # Keep conifer-dominated pixels.
  message(msgHdr, "filter: ≥25 % conifer")
  ids <- which(decid[] < 75)
  rm(decid); gc()
  
  # Keep pixels with enough canopy cover.
  message(msgHdr, "filter: ≥25 % canopy cover")
  cov_vals <- cover[][ids]
  ids <- ids[cov_vals >= 25]
  rm(cover, cov_vals); gc()
  
  # Keep pixels with no recent disturbance in the last 20 years.
  d <- combine_last_disturbance_year(baseYear, distYearPreFile, distYearMainFile)
  if (!is.null(d)) {
    message(msgHdr, "filter: ≥20 years since last disturbance")
    lastYr <- d[][ids]
    keep <- is.na(lastYr) | ((baseYear - lastYr) > 20)
    ids <- ids[keep]
    rm(d, lastYr, keep); gc()
  }
  
  # Split the remaining pixels into young vs mature conifer.
  age_vals <- age[][ids]
  young_ids  <- ids[age_vals > 49 & age_vals < 70]
  mature_ids <- ids[age_vals >= 70]
  rm(age_vals, ids); gc()
  
  writeMask <- function(ids, stub) {
    vals <- rep(NA_integer_, terra::ncell(age))
    vals[ids] <- 1L
    out <- terra::setValues(age, vals)
    fn  <- file.path(raw_dir, sprintf("%s%d_%s.tif", stub, baseYear, id))
    terra::writeRaster(out, filename = fn, datatype = "INT1U", overwrite = TRUE)
    invisible(fn)
  }
  
  message(msgHdr, "writing young and mature masks")
  writeMask(young_ids,  "youngConifer")
  writeMask(mature_ids, "matureConifer")
  message(msgHdr, "done")
}

### 3.2 Wetlands -----------------------------------------------------------------

Wetlands <- function(ageFile, landPosFile, lccFile,
                     canopyFile, decidFile,
                     youngFile, matureFile,
                     baseYear) {
  
  id  <- stringr::str_extract(ageFile, "tile\\d+")
  hdr <- paste0("Wetlands ", baseYear, " | ", id, " : ")
  
  age     <- terra::rast(ageFile)
  landPos <- terra::rast(landPosFile)
  lcc     <- terra::rast(lccFile)
  canopy  <- terra::rast(canopyFile)
  decid   <- terra::rast(decidFile)
  youngR  <- terra::rast(youngFile)
  matureR <- terra::rast(matureFile)
  
  terra::compareGeom(age, landPos, lcc, canopy, decid, youngR, matureR, stopOnError = TRUE)
  message(hdr, "geometry check complete")
  
  # Start with wetland land-position pixels.
  sel <- (landPos == 5)
  
  # Remove pixels already assigned as young or mature conifer.
  sel <- sel & is.na(youngR) & is.na(matureR)
  
  # Remove water.
  sel <- sel & (lcc != 8)
  
  # Remove forest-like wetland pixels using the current rule.
  sel <- sel & !(decid <= 75 & age >= 50 & canopy >= 30)
  
  out <- terra::ifel(sel, 1L, NA_integer_)
  
  outName <- file.path(raw_dir, sprintf("wetland%d_%s.tif", baseYear, id))
  terra::writeRaster(out, outName, datatype = "INT1U", overwrite = TRUE)
  
  # Stop if the written tile is empty.
  p <- terra::global(!is.na(out), "mean", na.rm = FALSE)[1, 1]
  if (is.na(p) || p == 0) {
    stop(hdr, "QA FAIL: wetland tile wrote as empty (all NA). File: ", outName)
  }
  
  message(hdr, sprintf("mask written (non-NA fraction = %.4f)", p))
  invisible(outName)
}

### 3.3 Natural disturbance ------------------------------------------------------

NaturalDisturbance <- function(baseYear,
                               dTypePreFile, dYearPreFile,
                               dTypeMainFile, dYearMainFile,
                               lccRefFile, lcc2020File,
                               wetRefFile) {
  
  id <- stringr::str_extract(dTypePreFile, "tile\\d+")
  
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    GDAL_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1"
  )
  
  dTypePre  <- terra::rast(dTypePreFile)
  dYearPre  <- terra::rast(dYearPreFile)
  dTypeMain <- terra::rast(dTypeMainFile)
  dYearMain <- terra::rast(dYearMainFile)
  lccRef    <- terra::rast(lccRefFile)
  wetRef    <- terra::rast(wetRefFile)
  
  terra::compareGeom(dTypePre, dYearPre, dTypeMain, dYearMain, lccRef, wetRef)
  
  delta   <- 20L
  pre_lo  <- baseYear - delta
  pre_hi  <- min(baseYear, 1984L)
  main_lo <- max(1985L, baseYear - delta)
  main_hi <- baseYear
  
  FIRE_CODES_PRE  <- 2L
  FIRE_CODES_MAIN <- 1L
  
  wopt_i1 <- list(
    datatype = "INT1U",
    gdal = c(
      "COMPRESS=ZSTD", "ZSTD_LEVEL=9", "TILED=YES",
      "BLOCKXSIZE=512", "BLOCKYSIZE=512", "BIGTIFF=YES"
    )
  )
  
  and_masks <- function(rs, outfile) {
    terra::lapp(
      rs,
      fun = function(...) {
        v <- Reduce(`*`, list(...))
        v[v > 0L] <- 1L
        as.integer(v)
      },
      filename = outfile,
      overwrite = TRUE,
      wopt = wopt_i1
    )
  }
  
  or_masks <- function(a, b, outfile) {
    terra::lapp(
      c(a, b),
      fun = function(x, y) as.integer(pmax(x, y, na.rm = TRUE)),
      filename = outfile,
      overwrite = TRUE,
      wopt = wopt_i1
    )
  }
  
  pre_type  <- terra::ifel(dTypePre == FIRE_CODES_PRE, 1L, 0L, filename = tempfile(fileext = ".tif"), overwrite = TRUE, wopt = wopt_i1)
  pre_notna <- terra::ifel(!is.na(dYearPre),          1L, 0L, filename = tempfile(fileext = ".tif"), overwrite = TRUE, wopt = wopt_i1)
  pre_ge    <- terra::ifel(dYearPre >= pre_lo,        1L, 0L, filename = tempfile(fileext = ".tif"), overwrite = TRUE, wopt = wopt_i1)
  pre_le    <- terra::ifel(dYearPre <= pre_hi,        1L, 0L, filename = tempfile(fileext = ".tif"), overwrite = TRUE, wopt = wopt_i1)
  
  pre_fire01 <- and_masks(c(pre_type, pre_notna, pre_ge, pre_le), outfile = tempfile(fileext = ".tif"))
  rm(pre_type, pre_notna, pre_ge, pre_le); gc()
  
  main_type  <- terra::ifel(dTypeMain == FIRE_CODES_MAIN, 1L, 0L, filename = tempfile(fileext = ".tif"), overwrite = TRUE, wopt = wopt_i1)
  main_notna <- terra::ifel(!is.na(dYearMain),           1L, 0L, filename = tempfile(fileext = ".tif"), overwrite = TRUE, wopt = wopt_i1)
  main_ge    <- terra::ifel(dYearMain >= main_lo,        1L, 0L, filename = tempfile(fileext = ".tif"), overwrite = TRUE, wopt = wopt_i1)
  main_le    <- terra::ifel(dYearMain <= main_hi,        1L, 0L, filename = tempfile(fileext = ".tif"), overwrite = TRUE, wopt = wopt_i1)
  
  main_fire01 <- and_masks(c(main_type, main_notna, main_ge, main_le), outfile = tempfile(fileext = ".tif"))
  rm(main_type, main_notna, main_ge, main_le); gc()
  
  isFire01 <- or_masks(pre_fire01, main_fire01, outfile = tempfile(fileext = ".tif"))
  rm(pre_fire01, main_fire01); gc()
  
  nonWet01 <- terra::app(
    wetRef,
    fun = function(w) as.integer(is.na(w)),
    filename = tempfile(fileext = ".tif"),
    overwrite = TRUE,
    wopt = wopt_i1
  )
  
  habitat01 <- terra::lapp(
    lccRef,
    fun = function(x) as.integer(!is.na(x) & x != 8L),
    filename = tempfile(fileext = ".tif"),
    overwrite = TRUE,
    wopt = wopt_i1
  )
  
  out_path <- file.path(raw_dir, sprintf("naturalDisturbance%d_%s.tif", baseYear, id))
  sel01 <- and_masks(c(isFire01, nonWet01, habitat01), outfile = tempfile(fileext = ".tif"))
  
  terra::lapp(
    sel01,
    fun = function(z) {
      z[z == 0L] <- NA_integer_
      z
    },
    filename = out_path,
    overwrite = TRUE,
    wopt = wopt_i1
  )
  
  invisible(out_path)
}

### 3.4 Recent harvest -----------------------------------------------------------

RecentHarvest <- function(baseYear,
                          dTypePreFile, dYearPreFile,
                          dTypeMainFile, dYearMainFile,
                          wetRefFile) {
  
  id  <- stringr::str_extract(dTypePreFile, "tile\\d+")
  hdr <- paste0("Harvest ", baseYear, " | ", id, " : ")
  
  dTypePre  <- terra::rast(dTypePreFile)
  dYearPre  <- terra::rast(dYearPreFile)
  dTypeMain <- terra::rast(dTypeMainFile)
  dYearMain <- terra::rast(dYearMainFile)
  wet       <- terra::rast(wetRefFile)
  
  terra::compareGeom(dTypePre, dYearPre, dTypeMain, dYearMain, wet)
  message(hdr, "geometry check passed")
  
  wopt_i2 <- list(
    datatype = "INT2U",
    gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512")
  )
  
  wopt_i1 <- list(
    datatype = "INT1U",
    gdal = c("COMPRESS=LZW", "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512")
  )
  
  HARVEST_CODES_MAIN <- c(2)
  HARVEST_CODES_PRE  <- c(3)
  
  yr_pre <- terra::ifel(
    (dTypePre %in% HARVEST_CODES_PRE) & !is.na(dYearPre) & (dYearPre <= baseYear),
    dYearPre, NA,
    filename = tempfile(fileext = ".tif"),
    overwrite = TRUE,
    wopt = wopt_i2
  )
  
  yr_main <- terra::ifel(
    (dTypeMain %in% HARVEST_CODES_MAIN) & !is.na(dYearMain) & (dYearMain <= baseYear),
    dYearMain, NA,
    filename = tempfile(fileext = ".tif"),
    overwrite = TRUE,
    wopt = wopt_i2
  )
  
  yr_any <- terra::mosaic(
    yr_pre, yr_main,
    fun = "max",
    filename = tempfile(fileext = ".tif"),
    overwrite = TRUE,
    wopt = wopt_i2
  )
  
  yrs <- terra::app(
    yr_any,
    fun = function(x) as.integer(baseYear - x),
    filename = tempfile(fileext = ".tif"),
    overwrite = TRUE,
    wopt = wopt_i2
  )
  
  nonWet <- terra::app(
    wet,
    fun = function(w) as.integer(is.na(w)),
    filename = tempfile(fileext = ".tif"),
    overwrite = TRUE,
    wopt = wopt_i1
  )
  
  mstack <- c(yrs, nonWet)
  
  terra::app(
    mstack,
    fun = function(v) {
      z <- as.integer(!is.na(v[, 1]) & v[, 1] >= 0L & v[, 1] <= 5L & v[, 2] == 1L)
      z[z == 0L] <- NA_integer_
      z
    },
    filename = file.path(raw_dir, sprintf("harvest_0to5_%d_%s.tif", baseYear, id)),
    overwrite = TRUE,
    wopt = wopt_i1
  )
  
  terra::app(
    mstack,
    fun = function(v) {
      z <- as.integer(!is.na(v[, 1]) & v[, 1] > 5L & v[, 1] <= 20L & v[, 2] == 1L)
      z[z == 0L] <- NA_integer_
      z
    },
    filename = file.path(raw_dir, sprintf("harvest_6to20_%d_%s.tif", baseYear, id)),
    overwrite = TRUE,
    wopt = wopt_i1
  )
  
  message(hdr, "masks written")
}

### 3.5 Open woodland ------------------------------------------------------------

OpenWoodlandTile <- function(baseYear,
                             ageFile, canopyFile, decidFile, posFile,
                             hv0_5File, hv6_20File, fireFile) {
  
  id  <- stringr::str_extract(ageFile, "tile\\d+")
  hdr <- paste0("Open woodland ", baseYear, " | ", id, " : ")
  
  age    <- terra::rast(ageFile)
  canopy <- terra::rast(canopyFile)
  decid  <- terra::rast(decidFile)
  pos    <- terra::rast(posFile)
  hv0_5  <- terra::rast(hv0_5File)
  hv6_20 <- terra::rast(hv6_20File)
  fire   <- terra::rast(fireFile)
  
  terra::compareGeom(canopy, decid, age)
  terra::compareGeom(fire, hv0_5, age)
  terra::compareGeom(hv6_20, age)
  message(hdr, "geometry check passed")
  
  sel <- (decid[] < 50) &
    (age[] > 49) &
    (canopy[] < 25) &
    is.na(fire[]) &
    is.na(hv0_5[]) &
    is.na(hv6_20[]) &
    (pos[] != 5)
  
  v <- rep(NA_integer_, terra::ncell(age))
  v[which(sel)] <- 1L
  
  out <- terra::setValues(age, v)
  outFile <- file.path(raw_dir, paste0("openWoodland", baseYear, "_", id, ".tif"))
  terra::writeRaster(out, filename = outFile, datatype = "INT1U", overwrite = TRUE)
  
  message(hdr, "mask written")
}

### 3.6 Regenerating stands ------------------------------------------------------

RegeneratingForest <- function(baseYear,
                               ageFile, decidFile, posFile, canopyFile,
                               distYearPreFile = NULL,
                               distYearMainFile = NULL,
                               hv0_5File = NULL,
                               hv6_20File = NULL,
                               fireFile  = NULL) {
  
  id  <- stringr::str_extract(ageFile, "tile\\d+")
  hdr <- paste0("Regenerating stands ", baseYear, " | ", id, " : ")
  
  age    <- terra::rast(ageFile)
  decid  <- terra::rast(decidFile)
  pos    <- terra::rast(posFile)
  canopy <- terra::rast(canopyFile)
  
  terra::compareGeom(age, decid, pos, canopy)
  message(hdr, "geometry check complete")
  
  # Start with stands older than 20 years.
  sel <- (age[] > 20)
  
  # Remove pixels with a disturbance in the last 20 years.
  d <- combine_last_disturbance_year(baseYear, distYearPreFile, distYearMainFile)
  if (!is.null(d)) {
    sel <- sel & (is.na(d[]) | (baseYear - d[] > 20))
  }
  
  # Apply the remaining stand rules.
  sel <- sel & (pos[] != 5)
  sel <- sel & !(age[] > 50 & canopy[] < 25)
  sel <- sel & ((decid[] > 75) | (decid[] <= 75 & age[] < 50))
  
  # Exclude pixels already flagged as fire or harvest.
  if (!is.null(hv0_5File) && file.exists(hv0_5File) &&
      !is.null(hv6_20File) && file.exists(hv6_20File) &&
      !is.null(fireFile)  && file.exists(fireFile)) {
    
    hv0_5 <- terra::rast(hv0_5File)
    hv6_20 <- terra::rast(hv6_20File)
    fire  <- terra::rast(fireFile)
    
    terra::compareGeom(age, hv0_5, hv6_20, fire)
    sel <- sel & is.na(hv0_5[]) & is.na(hv6_20[]) & is.na(fire[])
  }
  
  v <- rep(NA_integer_, terra::ncell(decid))
  v[which(sel)] <- 1L
  
  out <- terra::setValues(decid, v)
  outFile <- file.path(raw_dir, paste0("regeneratingStand", baseYear, "_", id, ".tif"))
  terra::writeRaster(out, filename = outFile, datatype = "INT1U", overwrite = TRUE)
  
  message(hdr, "mask written")
}

## 4. LOAD SHARED TILE INPUTS ====================================================

### 4.1 Disturbance tile maps ----------------------------------------------------

d_year_pre <- named_by_tile(list.files(
  tile_dir,
  "^canlad_1965_1984_disturbanceYear_tile\\d+\\.tif$",
  full.names = TRUE
))

d_year_main <- named_by_tile(list.files(
  tile_dir,
  "1985_2020_YRt2.*tile\\d+\\.tif$",
  full.names = TRUE
))

d_type_pre <- named_by_tile(list.files(
  tile_dir,
  "^canlad_1965_1984_disturbanceType_tile\\d+\\.tif$",
  full.names = TRUE
))

d_type_main <- named_by_tile(list.files(
  tile_dir,
  "1985_2020_TYPE.*tile\\d+\\.tif$",
  full.names = TRUE
))

## 5. RUN RAW CLASSIFICATION DRIVERS ============================================

### 5.1 Mature conifers ----------------------------------------------------------

for (yr in run_years) {
  
  message("\n==============================")
  message("YEAR ", yr, " | Mature conifers")
  message("==============================")
  
  ageMap  <- age_map_by_year(yr)
  coverYr <- named_by_tile(list.files(tile_dir, paste0("att_closure.*", yr), full.names = TRUE))
  decidYr <- named_by_tile(list.files(tile_dir, paste0("prcB.*", yr), full.names = TRUE))
  
  ids_core <- Reduce(intersect, list(names(ageMap), names(coverYr), names(decidYr)))
  ids_dist <- union(names(d_year_pre), names(d_year_main))
  ids_mc   <- intersect(ids_core, ids_dist)
  
  for (id in ids_mc) {
    MatureConifers(
      ageFile = ageMap[[id]],
      coverFile = coverYr[[id]],
      decidFile = decidYr[[id]],
      distYearPreFile = d_year_pre[[id]],
      distYearMainFile = d_year_main[[id]],
      baseYear = yr
    )
  }
}

### 5.2 Wetlands -----------------------------------------------------------------

for (yr in run_years) {
  
  message("\n==============================")
  message("YEAR ", yr, " | Wetlands")
  message("==============================")
  
  landPosMap <- named_by_tile(list.files(
    tile_dir,
    paste0("^SCANFI_att_land_pos_S_", yr, "_v1_1_tile\\d+\\.tif$"),
    full.names = TRUE
  ))
  
  lccMap <- named_by_tile(list.files(
    tile_dir,
    paste0("^SCANFI_att_nfiLandCover_S_", yr, "_v1_1_tile\\d+\\.tif$"),
    full.names = TRUE
  ))
  
  ageMap  <- age_map_by_year(yr)
  coverYr <- named_by_tile(list.files(tile_dir, paste0("att_closure.*", yr), full.names = TRUE))
  decidYr <- named_by_tile(list.files(tile_dir, paste0("prcB.*", yr), full.names = TRUE))
  
  canopyMap <- coverYr
  yConMap   <- named_by_tile(list.files(raw_dir, paste0("youngConifer", yr), full.names = TRUE))
  mConMap   <- named_by_tile(list.files(raw_dir, paste0("matureConifer", yr), full.names = TRUE))
  
  ids_wet <- Reduce(intersect, list(
    names(ageMap),
    names(landPosMap),
    names(lccMap),
    names(canopyMap),
    names(decidYr),
    names(yConMap),
    names(mConMap)
  ))
  
  for (id in ids_wet) {
    Wetlands(
      ageFile = ageMap[[id]],
      landPosFile = landPosMap[[id]],
      lccFile = lccMap[[id]],
      canopyFile = canopyMap[[id]],
      decidFile = decidYr[[id]],
      youngFile = yConMap[[id]],
      matureFile = mConMap[[id]],
      baseYear = yr
    )
  }
}

### 5.3 Natural disturbance ------------------------------------------------------

for (yr in run_years) {
  
  message("\n==============================")
  message("YEAR ", yr, " | Natural disturbance (fire)")
  message("==============================")
  
  lccMap <- named_by_tile(list.files(
    tile_dir,
    paste0("^SCANFI_att_nfiLandCover_S_", yr, "_v1_1_tile\\d+\\.tif$"),
    full.names = TRUE
  ))
  
  wetMap <- named_by_tile(list.files(
    raw_dir,
    paste0("^wetland", yr, "_tile\\d+\\.tif$"),
    full.names = TRUE
  ))
  
  ids_fire <- Reduce(intersect, list(
    names(d_type_pre), names(d_year_pre),
    names(d_type_main), names(d_year_main),
    names(lccMap),
    names(wetMap)
  ))
  
  for (id in ids_fire) {
    NaturalDisturbance(
      baseYear = yr,
      dTypePreFile = d_type_pre[[id]],
      dYearPreFile = d_year_pre[[id]],
      dTypeMainFile = d_type_main[[id]],
      dYearMainFile = d_year_main[[id]],
      lccRefFile = lccMap[[id]],
      lcc2020File = lccMap[[id]],
      wetRefFile = wetMap[[id]]
    )
    gc()
    try(terra::tmpFiles(current = TRUE, remove = TRUE), silent = TRUE)
  }
}

### 5.4 Harvest ------------------------------------------------------------------

for (yr in run_years) {
  
  message("\n==============================")
  message("YEAR ", yr, " | Harvest")
  message("==============================")
  
  wetMap <- named_by_tile(list.files(
    raw_dir,
    paste0("^wetland", yr, "_tile\\d+\\.tif$"),
    full.names = TRUE
  ))
  
  ids_hv <- Reduce(intersect, list(
    names(d_type_pre), names(d_year_pre),
    names(d_type_main), names(d_year_main),
    names(wetMap)
  ))
  
  for (id in ids_hv) {
    RecentHarvest(
      baseYear = yr,
      dTypePreFile = d_type_pre[[id]],
      dYearPreFile = d_year_pre[[id]],
      dTypeMainFile = d_type_main[[id]],
      dYearMainFile = d_year_main[[id]],
      wetRefFile = wetMap[[id]]
    )
  }
}

### 5.5 Open woodland ------------------------------------------------------------

for (yr in run_years) {
  
  message("\n==============================")
  message("YEAR ", yr, " | Open woodland")
  message("==============================")
  
  posMap_att <- named_by_tile(list.files(
    tile_dir,
    paste0("^SCANFI_att_land_pos_S_", yr, "_v1_1_tile\\d+\\.tif$"),
    full.names = TRUE
  ))
  
  ageMap  <- age_map_by_year(yr)
  coverYr <- named_by_tile(list.files(tile_dir, paste0("att_closure.*", yr), full.names = TRUE))
  decidYr <- named_by_tile(list.files(tile_dir, paste0("prcB.*", yr), full.names = TRUE))
  
  canopyMap <- coverYr
  hv05Map   <- named_by_tile(list.files(raw_dir, paste0("harvest_0to5_", yr), full.names = TRUE))
  hv620Map  <- named_by_tile(list.files(raw_dir, paste0("harvest_6to20_", yr), full.names = TRUE))
  fireMap   <- named_by_tile(list.files(raw_dir, paste0("naturalDisturbance", yr, "_"), full.names = TRUE))
  
  ids_wo <- Reduce(intersect, list(
    names(ageMap),
    names(posMap_att),
    names(canopyMap),
    names(decidYr),
    names(hv05Map),
    names(hv620Map),
    names(fireMap)
  ))
  
  for (id in ids_wo) {
    OpenWoodlandTile(
      baseYear = yr,
      ageFile = ageMap[[id]],
      canopyFile = canopyMap[[id]],
      decidFile = decidYr[[id]],
      posFile = posMap_att[[id]],
      hv0_5File = hv05Map[[id]],
      hv6_20File = hv620Map[[id]],
      fireFile = fireMap[[id]]
    )
  }
}

### 5.6 Regenerating stands ------------------------------------------------------

for (yr in run_years) {
  
  message("\n==============================")
  message("YEAR ", yr, " | Regenerating stands")
  message("==============================")
  
  posMap_att <- named_by_tile(list.files(
    tile_dir,
    paste0("^SCANFI_att_land_pos_S_", yr, "_v1_1_tile\\d+\\.tif$"),
    full.names = TRUE
  ))
  
  ageMap  <- age_map_by_year(yr)
  coverYr <- named_by_tile(list.files(tile_dir, paste0("att_closure.*", yr), full.names = TRUE))
  decidYr <- named_by_tile(list.files(tile_dir, paste0("prcB.*", yr), full.names = TRUE))
  
  canopyMap <- coverYr
  hv05Map   <- named_by_tile(list.files(raw_dir, paste0("harvest_0to5_", yr), full.names = TRUE))
  hv620Map  <- named_by_tile(list.files(raw_dir, paste0("harvest_6to20_", yr), full.names = TRUE))
  fireMap   <- named_by_tile(list.files(raw_dir, paste0("naturalDisturbance", yr, "_"), full.names = TRUE))
  
  ids_reg <- Reduce(intersect, list(
    names(ageMap),
    names(decidYr),
    names(posMap_att),
    names(canopyMap),
    names(d_year_pre),
    names(d_year_main),
    names(hv05Map),
    names(hv620Map),
    names(fireMap)
  ))
  
  for (id in ids_reg) {
    RegeneratingForest(
      baseYear = yr,
      ageFile = ageMap[[id]],
      decidFile = decidYr[[id]],
      posFile = posMap_att[[id]],
      canopyFile = canopyMap[[id]],
      distYearPreFile = d_year_pre[[id]],
      distYearMainFile = d_year_main[[id]],
      hv0_5File = hv05Map[[id]],
      hv6_20File = hv620Map[[id]],
      fireFile = fireMap[[id]]
    )
  }
}

# NEXT STEP -> Continue to habitat_3_create_composite.R 


