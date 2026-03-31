# habitat_4_focal.R ----------

# Description:
# This script creates focal habitat proportions for each class and year using the raw binary habitat masks.
#
# Required inputs:
# - data/habitat/outputs/raw/*.tif created by habitat_2_classification.R
#
# Expected outputs:
# - data/habitat/outputs/focalMosaic/denom/denom_year_<year>.tif
# - data/habitat/outputs/focalMosaic/num/num_<class>_<year>.tif
# - data/habitat/outputs/focalMosaic/prop/focal_<class>_<year>.tif
#
# Notes:
# - This script uses the raw binary class masks, not the composite habitat tiles.
# - The focal proportion rasters created here are used by habitat_5_weights.R.

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
  library(parallel)
  library(stringr)
  library(tools)
})

### 1.3 User inputs --------------------------------------------------------------

# Choose which snapshot year(s) to run.
# Use any combination of: 1985, 2000, 2020.
run_years <- c(1985, 2000, 2020)

# Main input and output folders.
# Keep these as-is unless the project structure was changed.
base_dir  <- "data/habitat/outputs"
input_dir <- file.path(base_dir, "raw")
out_base  <- file.path(base_dir, "focalMosaic")
denom_dir <- file.path(out_base, "denom")
num_dir   <- file.path(out_base, "num")
prop_dir  <- file.path(out_base, "prop")

# Focal window radius, in metres.
focal_radius_m <- 1000

# Output scaling factor.
# Focal proportions are multiplied by this value before writing.
scale_factor <- 1000L

# Which phases to run.
# Use c("denom","num","combine") for a full run.
# Use a subset, such as c("combine"), to resume from existing intermediate files.
run_phases <- c("denom", "num", "combine")

# Parallel settings.
# ncores_focal = cores used inside each focal run
# nworkers_class = number of classes processed at the same time
ncores_focal   <- 3
nworkers_class <- 1

# Memory and write settings.
# Keep these values unless you intentionally want to change runtime behaviour.
memfrac_each <- 0.6 / (ncores_focal + 1)

wopt_fast <- list(
  datatype = "INT2U",
  gdal = c(
    "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512",
    "COMPRESS=ZSTD", "NUM_THREADS=1", "BIGTIFF=YES"
  )
)

### 1.4 Create output folders ----------------------------------------------------

dir.create(out_base,  showWarnings = FALSE, recursive = TRUE)
dir.create(denom_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(num_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(prop_dir,  showWarnings = FALSE, recursive = TRUE)

### 1.5 Apply terra settings -----------------------------------------------------

terraOptions(memfrac = memfrac_each, progress = 1)

## 2. DISCOVER INPUTS ============================================================

### 2.1 List raw habitat mask files ----------------------------------------------

all_files <- list.files(input_dir, pattern = "\\.tif$", full.names = TRUE)
stopifnot(length(all_files) > 0)

### 2.2 Define filename parsers --------------------------------------------------

extract_year <- function(x) {
  y1 <- str_match(x, "(\\d{4})(?=_tile\\d+)")[, 2]
  y2 <- str_match(x, "(\\d{4})(?!.*\\d)")[, 2]
  as.integer(ifelse(!is.na(y1), y1, y2))
}

extract_tile <- function(x) {
  as.integer(str_match(x, "tile(\\d+)")[, 2])
}

extract_class <- function(x) {
  sub("(_)?\\d{4}_tile\\d+\\.tif$", "", basename(x))
}

### 2.3 Build input index --------------------------------------------------------

meta <- data.frame(
  file  = all_files,
  year  = extract_year(all_files),
  tile  = extract_tile(all_files),
  class = extract_class(all_files),
  stringsAsFactors = FALSE
)

meta <- subset(meta, !is.na(year) & !is.na(tile) & year %in% run_years)
stopifnot(nrow(meta) > 0)

### 2.4 Build the focal kernel ---------------------------------------------------

# This creates a circular focal window using the grid of the first input raster.
# This converts the kernel to raw 1/0 counts, which is important for the
# later numerator and denominator sums.

ref_rast <- rast(meta$file[1])

focal_matrix <- terra::focalMat(
  x = ref_rast,
  d = focal_radius_m,
  type = "circle"
)

focal_matrix <- (focal_matrix > 0) * 1L

message(
  "Inputs: ", nrow(meta), " rasters | years {",
  paste(sort(unique(meta$year)), collapse = ", "),
  "} | tiles {",
  paste(sort(unique(meta$tile)), collapse = ", "),
  "}."
)

## 3. DEFINE HELPER FUNCTIONS ====================================================

### 3.1 Run focal in parallel with halo tiling -----------------------------------

# This splits a raster into overlapping stripes, runs focal() on each stripe,
# then mosaics the results back together.

parallel_focal <- function(x, w,
                           fun = "sum",
                           na.rm = TRUE,
                           expand = FALSE,
                           ncores = 2,
                           filename,
                           wopt = list(datatype = "INT2U"),
                           chunks = max(12L, ncores * 6L)) {
  
  stopifnot(
    inherits(x, "SpatRaster"),
    is.matrix(w),
    chunks >= 1,
    ncores >= 1,
    is.character(filename) && nchar(filename) > 0
  )
  
  # Make sure workers can reopen the raster from disk.
  x_path <- sources(x)
  
  if (length(x_path) == 0L || !file.exists(x_path[1])) {
    x_path <- tempfile(pattern = "pfoc_src_", fileext = ".tif", tmpdir = tempdir())
    writeRaster(
      x, x_path,
      overwrite = TRUE,
      datatype = "FLT4S",
      gdal = c(
        "TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512",
        "COMPRESS=ZSTD", "NUM_THREADS=1", "BIGTIFF=YES"
      )
    )
  } else {
    x_path <- x_path[1]
  }
  
  xr0  <- rast(x_path)
  R    <- nrow(xr0)
  halo <- floor((nrow(w) - 1L) / 2L)
  
  brks   <- unique(round(seq(0, R, length.out = max(1L, chunks) + 1)))
  ranges <- Map(function(a, b) c(a + 1L, b), head(brks, -1), tail(brks, -1))
  ranges <- Filter(function(rr) rr[1] <= rr[2], ranges)
  
  focal_one <- function(rr, x_path, w, fun, na.rm, expand, wopt, halo) {
    xr    <- rast(x_path)
    e     <- ext(xr)
    yres  <- res(xr)[2]
    ymax0 <- ymax(e)
    R     <- nrow(xr)
    
    rmin   <- rr[1]
    rmax   <- rr[2]
    rmin_h <- max(1L, rmin - halo)
    rmax_h <- min(R,  rmax + halo)
    
    ymin_h <- ymax0 - rmax_h * yres
    ymax_h <- ymax0 - (rmin_h - 1L) * yres
    
    xi <- crop(xr, ext(xmin(e), xmax(e), ymin_h, ymax_h), snap = "out")
    
    tmp1 <- file.path(tempdir(), sprintf("pfoc_%d_%d.tif", rmin, rmax))
    yi <- terra::focal(
      x = xi,
      w = w,
      fun = fun,
      na.rm = na.rm,
      expand = expand,
      filename = tmp1,
      overwrite = TRUE,
      wopt = wopt
    )
    
    # Remove the halo before writing the output stripe.
    ymin <- ymax0 - rmax * yres
    ymax <- ymax0 - (rmin - 1L) * yres
    yi2  <- crop(yi, ext(xmin(e), xmax(e), ymin, ymax), snap = "out")
    
    tmp2 <- file.path(tempdir(), sprintf("pfoc_out_%d_%d.tif", rmin, rmax))
    writeRaster(
      yi2, tmp2,
      overwrite = TRUE,
      datatype = wopt$datatype,
      gdal = wopt$gdal
    )
    
    tmp2
  }
  
  if (length(ranges) == 1L || ncores == 1L) {
    parts <- lapply(
      ranges, focal_one,
      x_path, w, fun, na.rm, expand, wopt, halo
    )
  } else {
    cl <- parallel::makeCluster(min(ncores, length(ranges)))
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    
    parallel::clusterExport(
      cl,
      varlist = c("x_path", "w", "fun", "na.rm", "expand", "wopt", "halo", "memfrac_each"),
      envir = environment()
    )
    
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(terra))
      terraOptions(memfrac = memfrac_each, progress = 0)
      NULL
    })
    
    parts <- parallel::parLapplyLB(
      cl, ranges, focal_one,
      x_path, w, fun, na.rm, expand, wopt, halo
    )
  }
  
  rs <- lapply(parts, rast)
  
  if (length(rs) == 1L) {
    writeRaster(
      rs[[1]], filename = filename,
      overwrite = TRUE,
      datatype = wopt$datatype,
      gdal = wopt$gdal
    )
    return(rast(filename))
  } else {
    return(do.call(
      terra::mosaic,
      c(
        rs,
        list(
          fun = "first",
          filename = filename,
          overwrite = TRUE,
          wopt = wopt
        )
      )
    ))
  }
}

### 3.2 Mosaic one class for one year --------------------------------------------

# This mosaics all tiles for one habitat class and one year into one year-wide raster.

mosaic_class_for_year <- function(class_files,
                                  wopt = list(datatype = "INT2U"),
                                  fun = "max",
                                  template = NULL,
                                  tmp_dir = tempdir()) {
  
  stopifnot(length(class_files) >= 1L, all(file.exists(class_files)))
  
  tmp <- tempfile(pattern = "mos_", tmpdir = tmp_dir, fileext = ".tif")
  
  rs <- lapply(class_files, function(f) {
    r <- try(rast(f), silent = TRUE)
    if (inherits(r, "try-error")) {
      stop("Failed to open ", f)
    }
    r
  })
  
  mos <- do.call(
    terra::mosaic,
    c(
      rs,
      list(
        fun = fun,
        filename = tmp,
        overwrite = TRUE,
        wopt = wopt
      )
    )
  )
  
  if (!is.null(template)) {
    mos <- terra::resample(mos, template, method = "near")
  }
  
  mos
}

### 3.3 Optional smoke test ------------------------------------------------------

# Runs only in interactive sessions.
# This is a small check that the focal code still works before a full run.

if (interactive()) {
  message("Running smoke test ...")
  
  k <- nrow(focal_matrix)
  need <- k + 20
  
  nrows_use <- min(need, nrow(ref_rast))
  ncols_use <- min(need, ncol(ref_rast))
  
  rtest <- ref_rast[1:nrows_use, 1:ncols_use, drop = FALSE]
  
  tf <- tempfile(fileext = ".tif")
  dummy <- parallel_focal(rtest, focal_matrix, filename = tf, ncores = 2)
  
  message("Smoke test OK: ", tf)
}

## 4. DEFINE PROCESSING PHASES ===================================================

### 4.1 Compute the denominator raster -------------------------------------------

# This builds a year-wide valid-habitat mask across all classes, then focal-sums it.
# The denominator is the number of non-NA habitat cells in the focal window.

compute_year_denominator <- function(df_year, focal_matrix, out_dir, wopt, ncores = 2) {
  y <- unique(df_year$year)
  stopifnot(length(y) == 1L)
  
  out <- file.path(out_dir, sprintf("denom_year_%d.tif", y))
  if (file.exists(out)) {
    return(out)
  }
  
  class_map     <- split(df_year$file, df_year$class)
  class_mosaics <- lapply(class_map, mosaic_class_for_year, wopt = wopt)
  
  template_year <- class_mosaics[[1]]
  class_mosaics <- lapply(class_mosaics, function(r) {
    resample(r, template_year, method = "near")
  })
  
  cm_stack    <- rast(class_mosaics)
  any_present <- app(cm_stack, function(v) as.integer(any(!is.na(v))))
  valid_mask  <- any_present
  valid_mask[valid_mask == 0] <- NA
  
  parallel_focal(
    x = valid_mask,
    w = focal_matrix,
    fun = "sum",
    na.rm = TRUE,
    expand = FALSE,
    ncores = ncores,
    filename = out,
    wopt = modifyList(wopt, list(datatype = "INT2U")),
    chunks = max(12L, ncores * 6L)
  )
  
  out
}

### 4.2 Compute numerator rasters ------------------------------------------------

# This builds one year-wide focal count raster per habitat class.

compute_year_numerators <- function(df_year, denom_path, focal_matrix, out_dir, wopt,
                                    ncores_focal = 2, nworkers_class = 1) {
  y <- unique(df_year$year)
  stopifnot(length(y) == 1L)
  
  d <- rast(denom_path)
  class_files_map <- split(df_year$file, df_year$class)
  classes <- names(class_files_map)
  
  work_fun <- function(cls) {
    out <- file.path(out_dir, sprintf("num_%s_%d.tif", cls, y))
    if (file.exists(out)) {
      return(out)
    }
    
    class_mos <- mosaic_class_for_year(class_files_map[[cls]], wopt = wopt, template = d)
    
    # Convert the class raster to 1 / 0 / NA before focal-summing.
    x01 <- ifel(class_mos == 1, 1, NA)
    x01[is.na(x01)] <- 0
    x01 <- mask(x01, d)
    
    tmp <- file.path(
      tempdir(),
      paste0("num_", tools::md5sum(paste(class_files_map[[cls]], collapse = "|")), ".tif")
    )
    
    num <- parallel_focal(
      x = x01,
      w = focal_matrix,
      fun = "sum",
      na.rm = TRUE,
      expand = FALSE,
      ncores = ncores_focal,
      filename = tmp,
      wopt = modifyList(wopt, list(datatype = "INT2U"))
    )
    
    if (!compareGeom(num, d, stopOnError = FALSE)) {
      num <- resample(num, d, method = "near")
    }
    
    writeRaster(
      num, out,
      overwrite = TRUE,
      datatype = "INT2U",
      gdal = wopt$gdal
    )
    
    out
  }
  
  if (nworkers_class > 1L) {
    cl <- makeCluster(min(nworkers_class, length(classes)))
    on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
    
    clusterExport(
      cl,
      varlist = c(
        "class_files_map", "y", "focal_matrix", "wopt", "denom_path",
        "mosaic_class_for_year", "parallel_focal", "ncores_focal", "d"
      ),
      envir = environment()
    )
    
    clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(terra))
      terraOptions(memfrac = 0.6, progress = 0)
      NULL
    })
    
    res <- parLapplyLB(cl, classes, work_fun)
  } else {
    res <- lapply(classes, work_fun)
  }
  
  unlist(res, use.names = FALSE)
}

### 4.3 Combine numerators and denominator ---------------------------------------

# This converts each numerator to a focal proportion raster and writes the final outputs.

combine_num_denom <- function(year, denom_dir, num_dir, out_dir, scale_factor, wopt) {
  d_path <- file.path(denom_dir, sprintf("denom_year_%d.tif", year))
  if (!file.exists(d_path)) {
    stop("Missing denominator: ", d_path)
  }
  
  d <- rast(d_path)
  
  num_files <- list.files(
    num_dir,
    pattern = sprintf("^num_.*_%d\\.tif$", year),
    full.names = TRUE
  )
  
  if (!length(num_files)) {
    stop("No numerators found for year ", year, " in ", num_dir)
  }
  
  for (nf in num_files) {
    cls <- sub(sprintf("^num_(.*)_%d\\.tif$", year), "\\1", basename(nf))
    
    out <- file.path(out_dir, sprintf("focal_%s_%d.tif", cls, year))
    if (file.exists(out)) {
      message("Skipping existing: ", out)
      next
    }
    
    n <- rast(nf)
    
    if (!compareGeom(n, d, stopOnError = FALSE)) {
      n <- resample(n, d, method = "near")
    }
    
    ratio <- (n / d) * scale_factor
    ratio <- mask(ratio, d > 0)
    ratio <- round(ratio)
    
    ratio[ratio < 0]     <- 0
    ratio[ratio > 65534] <- 65534
    
    writeRaster(
      ratio, out,
      overwrite = TRUE,
      datatype = "INT2U",
      gdal = wopt$gdal
    )
    
    message("Wrote: ", out)
  }
  
  invisible(TRUE)
}

## 5. RUN FOCAL PROCESSING =======================================================

### 5.1 Process each selected year -----------------------------------------------

for (y in sort(unique(meta$year))) {
  message("=== YEAR ", y, " ===")
  
  df_y <- meta[meta$year == y, , drop = FALSE]
  
  if ("denom" %in% run_phases) {
    message("Phase A: Denominator (parallel focal ×", ncores_focal, ")...")
    denom_path <- compute_year_denominator(
      df_y, focal_matrix, denom_dir, wopt_fast,
      ncores = ncores_focal
    )
    message("  Denominator ready: ", denom_path)
  } else {
    denom_path <- file.path(denom_dir, sprintf("denom_year_%d.tif", y))
    if (!file.exists(denom_path)) {
      stop("Denominator missing; run phase 'denom' first.")
    }
  }
  
  if ("num" %in% run_phases) {
    message(
      "Phase B: Numerators (", nworkers_class,
      " class worker(s); parallel focal ×", ncores_focal, ")..."
    )
    
    nouts <- compute_year_numerators(
      df_y, denom_path, focal_matrix, num_dir, wopt_fast,
      ncores_focal = ncores_focal,
      nworkers_class = nworkers_class
    )
    
    message("  Numerators written: ", length(nouts))
  }
  
  if ("combine" %in% run_phases) {
    message("Phase C: Combine...")
    combine_num_denom(y, denom_dir, num_dir, prop_dir, scale_factor, wopt_fast)
  }
}

message("Done. Outputs:")
message("  Denominator(s): ", denom_dir)
message("  Numerator(s):   ", num_dir)
message("  Proportions:    ", prop_dir, "  (INT2U, ×", scale_factor, ")")

# NEXT STEP -> Continue to habitat_5_weights.R 

