# habitat_1_preprocessing.R ----------

# Description:
# Run this script first for the habitat workflow.
# It prepares the habitat input rasters by:
# - checking the required local raster inputs
# - cropping and masking them to the study area
# - checking that processed rasters share the same grid
# - tiling the processed rasters
# - preparing the managed forest and pre-1985 fire inputs used later in habitat classification
#
# Required inputs:
# - data/habitat/raw_layers/ (local habitat raster inputs listed below)
# - data/habitat/GIS/SA.shp
# - data/habitat/GIS/managed_forest/Canada_MFv2017.tif
#   or automatic download fallback if missing
# - data/habitat/GIS/NFDB/NFDB_poly_20210707.shp
#   or automatic download/unzip fallback if missing
#
# Expected outputs:
# - data/habitat/GIS/Eastern_*.tif
# - data/habitat/GIS/tiles/*.tif
# - data/habitat/GIS/Eastern_ManagedForest.tif
# - data/habitat/GIS/NFDB_raster_1965_1985.tif

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
  library(SpaDES.tools)
  library(reproducible)
  library(terra)
  library(magrittr)
  library(stringr)
  library(sf)
  library(httr)
})

### 1.3 User inputs --------------------------------------------------------------

# Choose which snapshot years to prepare.
# Remove years you do not want to run.
run_years <- c(1985, 2000, 2020)

# Folder containing the raw habitat rasters listed below.
# Keep this as-is unless the raw inputs were moved.
local_data_dir <- "data/habitat/raw_layers"

# Tiling settings.
# Keep these values unless you intentionally want a different tile layout.
nx <- 4
ny <- 2
buffer_px <- c(35, 35)

# Disturbance-prep template raster used later in this script.
# Note: This just verifies that all raster grids are aligned, you can use any SCANFI file for this. 
age_template_for_disturbance <- "data/habitat/GIS/Eastern_SCANFI_att_TIDO_age_SW_2020_v1.tif"

### 1.4 Create output folders ----------------------------------------------------

checkPath("data/habitat/outputs/raw", create = TRUE)
checkPath("data/habitat/GIS/tiles", create = TRUE)

### 1.5 Define download helpers --------------------------------------------------

# Checks whether a file exists locally. If not, it tries one or more download URLs.
download_file_if_missing <- function(dest_path, urls, label, manual_url) {
  dir.create(dirname(dest_path), recursive = TRUE, showWarnings = FALSE)
  
  if (file.exists(dest_path) && file.info(dest_path)$size > 0) {
    message(label, " found locally.")
    return(dest_path)
  }
  
  message(label, " not found locally. Attempting download...")
  
  for (u in urls) {
    ok <- FALSE
    
    try({
      resp <- httr::GET(
        u,
        httr::write_disk(dest_path, overwrite = TRUE),
        httr::timeout(300)
      )
      
      ok <- httr::status_code(resp) >= 200 &&
        httr::status_code(resp) < 300 &&
        file.exists(dest_path) &&
        file.info(dest_path)$size > 0
    }, silent = TRUE)
    
    if (ok) {
      message(label, " downloaded successfully.")
      return(dest_path)
    }
    
    if (file.exists(dest_path) && isTRUE(file.info(dest_path)$size == 0)) {
      file.remove(dest_path)
    }
  }
  
  stop(
    label, " could not be downloaded automatically.\n",
    "Download manually from: ", manual_url, "\n",
    "Save it here: ", dest_path
  )
}

# Checks whether a shapefile exists locally. If not, it downloads a zip and extracts it.
download_zip_and_extract_if_missing <- function(zip_url,
                                                zip_dest,
                                                shp_path,
                                                label,
                                                manual_url = zip_url) {
  dir.create(dirname(zip_dest), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(shp_path), recursive = TRUE, showWarnings = FALSE)
  
  if (file.exists(shp_path)) {
    message(label, " found locally.")
    return(shp_path)
  }
  
  if (!file.exists(zip_dest) || file.info(zip_dest)$size == 0) {
    message(label, " not found locally. Attempting download...")
    
    ok <- FALSE
    
    try({
      resp <- httr::GET(
        zip_url,
        httr::write_disk(zip_dest, overwrite = TRUE),
        httr::timeout(300)
      )
      
      ok <- httr::status_code(resp) >= 200 &&
        httr::status_code(resp) < 300 &&
        file.exists(zip_dest) &&
        file.info(zip_dest)$size > 0
    }, silent = TRUE)
    
    if (!ok) {
      if (file.exists(zip_dest) && isTRUE(file.info(zip_dest)$size == 0)) {
        file.remove(zip_dest)
      }
      
      stop(
        label, " could not be downloaded automatically.\n",
        "Download manually from: ", manual_url, "\n",
        "Save the zip here: ", zip_dest, "\n",
        "Then extract it into: ", dirname(shp_path)
      )
    }
  }
  
  unzip(zip_dest, exdir = dirname(shp_path))
  
  if (!file.exists(shp_path)) {
    stop(
      label, " zip was downloaded, but the expected shapefile was not found after extraction.\n",
      "Check the contents of: ", dirname(shp_path)
    )
  }
  
  message(label, " is ready.")
  shp_path
}

## 2. DEFINE INPUT FILES =========================================================

### 2.1 Build the list of local habitat rasters ----------------------------------

age_files <- c(
  "1985" = "stand_age_1985_final.tif",
  "2000" = "stand_age_2000_final.tif",
  "2020" = "stand_age_2020_hybrid.tif"
)

closure_files <- c(
  "1985" = "SCANFI_att_closure_S_1985_v1_1.tif",
  "2000" = "SCANFI_att_closure_S_2000_v1_1.tif",
  "2020" = "SCANFI_att_closure_S_2020_v1_1.tif"
)

broadleaf_files <- c(
  "1985" = "SCANFI_sps_prcB_S_1985_v1_1.tif",
  "2000" = "SCANFI_sps_prcB_S_2000_v1_1.tif",
  "2020" = "SCANFI_sps_prcB_S_2020_v1_1.tif"
)

land_position_files <- c(
  "1985" = "SCANFI_att_land_pos_S_1985_v1_1.tif",
  "2000" = "SCANFI_att_land_pos_S_2000_v1_1.tif",
  "2020" = "SCANFI_att_land_pos_S_2020_v1_1.tif"
)

landcover_files <- c(
  "1985" = "SCANFI_att_nfiLandCover_S_1985_v1_1.tif",
  "2000" = "SCANFI_att_nfiLandCover_S_2000_v1_1.tif",
  "2020" = "SCANFI_att_nfiLandCover_S_2020_v1_1.tif"
)

disturbance_files <- c(
  "CanLaD_Latest_1985_2020_YRt2.tif",
  "CanLaD_Latest_1985_2020_TYPE.tif",
  "canlad_1965_1984_disturbanceType.tif",
  "canlad_1965_1984_disturbanceYear.tif"
)

selected_years_chr <- as.character(run_years)

gis_files <- c(
  unname(age_files[selected_years_chr]),
  unname(closure_files[selected_years_chr]),
  unname(broadleaf_files[selected_years_chr]),
  unname(land_position_files[selected_years_chr]),
  unname(landcover_files[selected_years_chr]),
  disturbance_files
)

canlad_and_scanfi <- file.path(local_data_dir, gis_files)
names(canlad_and_scanfi) <- gis_files

### 2.2 Check that required local rasters exist ---------------------------------

missing_local_files <- canlad_and_scanfi[!file.exists(canlad_and_scanfi)]

if (length(missing_local_files) > 0) {
  stop(
    "The following required input files were not found in ",
    local_data_dir, ":\n",
    paste(basename(missing_local_files), collapse = "\n")
  )
}

## 3. DEFINE STUDY AREA ==========================================================

# ## helper: download once, read with terra::vect
# download_unzip_vect <- function(url, subDir) {
#   dir <- file.path("data/habitat/GIS", subDir)
#   reproducible::checkPath(dir, create = TRUE)
#   
#   zipFile <- file.path(dir, basename(url))
#   if (!file.exists(zipFile)) {
#     message("Downloading ", basename(url))
#     GET(url, write_disk(zipFile, overwrite = TRUE))  
#   }
#   
#   shp <- unzip(zipFile, list = TRUE)$Name
#   shp <- shp[grepl("\\.shp$", shp, ignore.case = TRUE)][1]   
#   
#   if (!file.exists(file.path(dir, shp)))
#     unzip(zipFile, exdir = dir)
#   
#   terra::vect(file.path(dir, shp))
# }
# 
# ## Get provincial boundaries
# Canada <- download_unzip_vect(
#   url    = paste0("https://www12.statcan.gc.ca/census-recensement/2011/",
#                   "geo/bound-limit/files-fichiers/2016/lpr_000b16a_e.zip"),
#   subDir = "Canada"
# )
# 
# QuebecOntario <- Canada[Canada$PRUID %in% c("24", "35"), ]
# QuebecOntario <- terra::aggregate(QuebecOntario)
# 
# ## Get ecozones
# Ecozones <- download_unzip_vect(
#   url    = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip",
#   subDir = "Ecozones"
# )
# Ecozones <- Ecozones[Ecozones$ZONE_NAME %in%
#                        c("Taiga Shield", "Boreal Shield",
#                          "MixedWood Plain", "Hudson Plain")]
# 
# Ecozones <- terra::aggregate(Ecozones)
# sf_proj_network(TRUE)        
# ## Get the intersection
# QuebecOntario <- project(QuebecOntario, crs(Ecozones))
# SA <- terra::intersect(QuebecOntario, Ecozones)

### NOTE: ecozones became corrupted when reprojecting, so I did the intersection in Arc for now
     
### 3.1 Read the prepared study area (fallback option) --------------------------

SA <- terra::vect("data/habitat/GIS/SA.shp")  # This shapefile should be included 
                                              # in the packaged raw data
SA <- terra::aggregate(SA)

## 4. CROP AND MASK INPUT RASTERS ===============================================

### 4.1 Define processed output names --------------------------------------------

processed <- file.path(
  "data/habitat/GIS",
  paste0("Eastern_", basename(canlad_and_scanfi))
)

need_processing <- processed[!file.exists(processed)]

### 4.2 Crop and mask any rasters not yet processed ------------------------------

if (length(need_processing) > 0) {
  message("Cropping and masking ", length(need_processing), " raster(s) ...")
  
  bulk_post_process <- function(infile, out_name, SA) {
    r <- rast(infile)
    d_type <- "INT2S"
    SA <- project(SA, crs(r))
    tmp <- terra::crop(r, SA)
    
    terra::mask(
      tmp, SA,
      filename = out_name,
      datatype = d_type,
      overwrite = TRUE
    )
    
    gc()
  }
  
  mapply(
    bulk_post_process,
    infile = canlad_and_scanfi[basename(processed) %in% basename(need_processing)],
    out_name = need_processing,
    MoreArgs = list(SA = SA)
  )
}

## 5. CHECK RASTER GRID CONSISTENCY =============================================

### 5.1 List processed rasters ---------------------------------------------------

eastern <- list.files(
  "data/habitat/GIS",
  pattern = "^Eastern_.*\\.tif$",
  full.names = TRUE
)

### 5.2 Define raster grid signature helper --------------------------------------

grid_sig <- function(r) {
  c(
    nrow = nrow(r), ncol = ncol(r),
    xmin = xmin(r), xmax = xmax(r),
    ymin = ymin(r), ymax = ymax(r),
    resx = res(r)[1], resy = res(r)[2],
    crs  = crs(r)
  )
}

### 5.3 Find the modal grid ------------------------------------------------------

sigs <- lapply(eastern, function(f) grid_sig(rast(f)))
sig_mat <- do.call(rbind, sigs)

grid_tbl <- as.data.frame(sig_mat)
key_cols <- names(grid_tbl)

template_idx <- which.max(
  vapply(
    split(seq_len(nrow(grid_tbl)), interaction(grid_tbl[key_cols], drop = TRUE)),
    length,
    integer(1)
  )
)

template <- rast(eastern[template_idx])
message("Template grid = ", basename(eastern[template_idx]))

### 5.4 Realign any outlier rasters ----------------------------------------------

needs_fix <- which(
  !duplicated(grid_tbl[key_cols]) &
    seq_len(nrow(grid_tbl)) != template_idx
)

if (length(needs_fix)) {
  message("Fixing ", length(needs_fix), " raster(s) with mismatched grids ...")
  
  for (i in needs_fix) {
    bad_file <- eastern[i]
    r <- rast(bad_file)
    
    fixed <- if (!compareGeom(template, r, stopOnError = FALSE, crs = TRUE)) {
      terra::project(r, template, method = "near", mask = TRUE)
    } else {
      terra::resample(r, template, method = "near")
    }
    
    writeRaster(fixed, bad_file, overwrite = TRUE, datatype = "INT2S")
    message("Success: ", basename(bad_file), " aligned")
  }
}

### 5.5 Re-run the grid check ----------------------------------------------------

cell_counts <- vapply(eastern, function(f) ncell(rast(f)), numeric(1))

if (length(unique(cell_counts)) > 1) {
  stop("GIS error - unequal raster dimensions remain after auto-fix.")
} else {
  message("All cropped rasters now share identical dimensions.")
}


## 6. TILE PROCESSED RASTERS ====================================================

### 6.1 Set terra temp options ---------------------------------------------------

scratch_dir <- file.path(getwd(), "temp")
reproducible::checkPath(scratch_dir, create = TRUE)

terraOptions(
  tempdir = scratch_dir,
  memfrac = 0.30,
  todisk = TRUE,
  progress = 1
)

### 6.2 Define tile settings -----------------------------------------------------

tile_dir <- "data/habitat/GIS/tiles"
checkPath(tile_dir, create = TRUE)

### 6.3 Find rasters that still need tiling --------------------------------------

processed <- list.files(
  "data/habitat/GIS",
  pattern = "^Eastern_.*\\.tif$",
  full.names = TRUE
)

filenames_no_ext <- basename(processed) %>%
  str_remove("\\.tif$") %>%
  str_remove("^Eastern_")

not_tiled <- vapply(
  filenames_no_ext,
  FUN = function(stem) {
    length(list.files(
      path = tile_dir,
      pattern = paste0("^", stem, "_.*\\.tif$")
    ))
  },
  FUN.VALUE = integer(1)
)

missing_tiles <- processed[not_tiled < (nx * ny)]

### 6.4 Tile remaining rasters ---------------------------------------------------

if (length(missing_tiles) == 0L) {
  message("All rasters already tiled (", nx, " x ", ny, "). Step complete.")
} else {
  message("Tiling ", length(missing_tiles), " raster(s) ...")
  
  lapply(missing_tiles, function(to_tile) {
    message("  - ", basename(to_tile))
    r <- rast(to_tile)
    
    splitRaster(
      r,
      nx = nx,
      ny = ny,
      buffer = buffer_px,
      rType = "INT2S",
      path = tile_dir,
      fExt = ".tif"
    )
  })
  
  message("Tiling complete.")
}

## 7. DISTURBANCE MAPPING PREPROCESS ============================================

### 7.1 Reload study area --------------------------------------------------------

SA <- terra::vect("data/habitat/GIS/SA.shp")
SA <- terra::aggregate(SA)

### 7.2 Read the disturbance template raster -------------------------------------

ageRTM <- rast(age_template_for_disturbance)
SA <- terra::project(SA, ageRTM)

### 7.3 Prepare the managed forest raster ----------------------------------------

managed_forest_path <- file.path(
  project_root,
  "data", "habitat", "GIS", "managed_forest", "Canada_MFv2017.tif"
)

managed_forest_manual_url <- paste0(
  "https://drive.google.com/file/d/",
  "1W2EiRtHj_81ZyKk5opqMkRqCA1tRMMvB/view?usp=share_link"
)

managed_forest_download_urls <- c(
  paste0(
    "https://drive.google.com/uc?export=download&id=",
    "1W2EiRtHj_81ZyKk5opqMkRqCA1tRMMvB"
  ),
  managed_forest_manual_url
)

download_file_if_missing(
  dest_path = managed_forest_path,
  urls = managed_forest_download_urls,
  label = "Managed forest raster",
  manual_url = managed_forest_manual_url
)

ManagedForest <- rast(managed_forest_path)

# Before continuing, verify that values are:
# 11 = long-term tenure
# 12 = short-term tenure
# 13 = other
# 20 = Protected areas
# 31 = Federal reserve
# 32 = Indian Reserve
# 33 = Restricted
# 40 = Treaty and Settlement
# 50 = Private forests

ManagedForest <- terra::project(ManagedForest, ageRTM)
ManagedForest <- terra::mask(
  ManagedForest, SA,
  filename = "data/habitat/GIS/Eastern_ManagedForest.tif",
  overwrite = TRUE
)

gc()

### 7.4 Tile the managed forest raster -------------------------------------------

splitRaster(
  ManagedForest,
  nx = nx,
  ny = ny,
  buffer = buffer_px,
  rType = "INT1U",
  path = "data/habitat/GIS/tiles",
  fExt = ".tif"
)

rm(ManagedForest)

### 7.5 Read and filter the NFDB fire polygons -----------------------------------

nfdb_zip_url <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"

nfdb_zip_path <- file.path(
  project_root,
  "data", "habitat", "GIS", "NFDB", "NFDB_poly.zip"
)

nfdb_shp_path <- file.path(
  project_root,
  "data", "habitat", "GIS", "NFDB", "NFDB_poly_20210707.shp"
)

download_zip_and_extract_if_missing(
  zip_url = nfdb_zip_url,
  zip_dest = nfdb_zip_path,
  shp_path = nfdb_shp_path,
  label = "NFDB fire polygons",
  manual_url = nfdb_zip_url
)

NFDB <- terra::vect(nfdb_shp_path)

NFDB <- NFDB[NFDB$YEAR > 1964 & NFDB$YEAR < 1986, ]
NFDB <- terra::project(NFDB, SA)
NFDB <- terra::crop(NFDB, SA)

NFDB <- sf::st_as_sf(NFDB)
NFDB <- st_cast(NFDB, to = "MULTIPOLYGON")
NFDB <- st_transform(NFDB, crs(ageRTM))

### 7.6 Rasterize the NFDB fire layer --------------------------------------------

NFDBras <- rasterize(
  NFDB,
  ageRTM,
  field = "YEAR",
  filename = "data/habitat/GIS/NFDB_raster_1965_1985.tif",
  overwrite = TRUE
)

rm(NFDBras)
gc()

NFDBras <- rast("data/habitat/GIS/NFDB_raster_1965_1985.tif")

### 7.7 Tile the NFDB raster if needed -------------------------------------------

fire_tiles <- list.files(
  "data/habitat/GIS/tiles",
  pattern = "NFDB",
  full.names = TRUE
)

if (length(fire_tiles) < 1) {
  names(NFDBras) <- "Eastern_NFDB"
  
  SpaDES.tools::splitRaster(
    NFDBras,
    nx = nx,
    ny = ny,
    buffer = buffer_px,
    rType = "INT2U",
    path = "data/habitat/GIS/tiles",
    fExt = ".tif"
  )
}

fire_tiles <- list.files(
  "data/habitat/GIS/tiles",
  pattern = "NFDB",
  full.names = TRUE
)

# NEXT STEP -> Continue to habitat_2_classification.R 

