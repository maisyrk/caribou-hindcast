# postprocess_4_extract_herd_values.R ----------

# Description:
# Run this script to extract:
# - habitat proportions within each herd polygon, by year
# - non-habitat model predictors within each herd polygon, by year
#
# Then merge both into one final wide model matrix.
#
# Required inputs:
# - data/habitat/outputs/compositeHabitat/compositeHabitat_<year>_tile1.tif ... tile8.tif
# - data/post_processing/herds/Aire_repartition_populations_caribouForestier.shp
# - data/post_processing/herds/Caribou_range_boundary.shp
# - data/mines/processed/density_layers/mine_density_<year>.tif
# - data/mines/processed/RS_layers/RS_masked_<year>.tif
# - data/roads/unpaved/unpaved_road_density_<year>.tif
# - data/roads/unpaved/unpaved_roads_zoi_residual_suitability_<year>.tif
# - data/roads/paved/paved_road_density.tif
# - data/roads/paved/paved_roads_zoi_residual_suitability.tif
#
# Expected outputs:
# - data/post_processing/habitat_props/habitat_proportions_long.csv
# - data/post_processing/habitat_props/habitat_proportions_matrix_wide.csv
# - data/post_processing/model_matrix/nonhab_means_long.csv
# - data/post_processing/model_matrix/nonhab_means_matrix_wide.csv
# - data/post_processing/model_matrix/final_model_matrix_wide.csv
#
# Notes:
# - Habitat values are extracted from the composite habitat tiles.
# - Non-habitat values are mean raster values by herd polygon.
# - The final merge key is herd_id_safe + "_" + year.

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
  library(dplyr)
  library(stringr)
  library(readr)
  library(rlang)
})

options(stringsAsFactors = FALSE)

### 1.3 User inputs --------------------------------------------------------------

# Choose which snapshot year(s) to run.
# Use any combination of: 1985, 2000, 2020.
run_years <- c(1985, 2000, 2020)

# Choose which non-habitat predictors to extract.
# Valid options:
# "dist_paved", "dist_unpaved", "dist_mines",
# "dens_paved", "dens_unpaved", "dens_mines"
selected_nonhab_predictors <- c("dens_unpaved")

# Keep the full 6-column non-habitat schema?
# TRUE  = unselected predictors are added as NA columns
# FALSE = only selected predictors are kept
keep_unselected_as_na <- TRUE

# Add use_<predictor> columns to the non-habitat outputs?
add_use_flags <- FALSE

# Fill value for unselected predictors when keep_unselected_as_na = TRUE
unselected_fill_value <- NA_real_

# Extraction settins 
# Only counts pixels fully inside the herd boundary = faster processing
habitat_touches <- FALSE
nonhab_touches <- FALSE

# Resume non-habitat extraction from saved year checkpoints if available?
resume_nonhab_checkpoints <- TRUE

# Simplify herd polygons before extraction.
herd_simplify_tolerance_m <- 60

# terra settings
terra_memfrac <- 0.7
terra_progress <- 1
terra_threads <- max(1, parallel::detectCores(logical = FALSE) - 1)


### 1.4 Define folders and terra temp settings -----------------------------------

temp_dir <- file.path(project_root, "temp")
dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

terraOptions(
  tempdir = temp_dir,
  progress = terra_progress,
  memfrac = terra_memfrac,
  threads = terra_threads
)

try(terra::setGDALconfig("GDAL_CACHEMAX", "4096"), silent = TRUE)
try(terra::setGDALconfig("GDAL_NUM_THREADS", "ALL_CPUS"), silent = TRUE)

## 2. DEFINE PATHS AND LOOKUPS ===================================================

### 2.1 Define main folders ------------------------------------------------------

hab_dir <- file.path(project_root, "data", "habitat", "outputs", "compositeHabitat")

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

habitat_out_dir <- file.path(project_root, "data", "post_processing", "habitat_props")
model_out_dir   <- file.path(project_root, "data", "post_processing", "model_matrix")
checkpoint_dir  <- file.path(model_out_dir, "checkpoints_nonhab")

dir.create(habitat_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(model_out_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(checkpoint_dir,  showWarnings = FALSE, recursive = TRUE)

### 2.2 Define log file ----------------------------------------------------------

log_file <- file.path(
  temp_dir,
  paste0("extract_herd_values_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
)

### 2.3 Define habitat legend ----------------------------------------------------

hab_vals <- 1:8

hab_codes <- c(
  "fire", "harv_y", "harv_o", "conf_y",
  "conf_o", "open_w", "wetland", "regen"
)
names(hab_codes) <- as.character(hab_vals)

### 2.4 Define non-habitat raster paths ------------------------------------------

mine_dens_tmpl <- file.path(
  project_root, "data", "mines", "processed", "density_layers",
  "mine_density_%d.tif"
)

mine_dist_tmpl <- file.path(
  project_root, "data", "mines", "processed", "RS_layers",
  "RS_masked_%d.tif"
)

unpaved_dens_tmpl <- file.path(
  project_root, "data", "roads", "unpaved",
  "unpaved_road_density_%d.tif"
)

unpaved_dist_tmpl <- file.path(
  project_root, "data", "roads", "unpaved",
  "unpaved_roads_zoi_residual_suitability_%d.tif"
)

paved_dens_file <- file.path(
  project_root, "data", "roads", "paved",
  "paved_road_density.tif"
)

paved_dist_file <- file.path(
  project_root, "data", "roads", "paved",
  "paved_roads_zoi_residual_suitability.tif"
)

predictor_defs <- list(
  dens_mines   = function(yr) sprintf(mine_dens_tmpl, yr),
  dist_mines   = function(yr) sprintf(mine_dist_tmpl, yr),
  dens_unpaved = function(yr) sprintf(unpaved_dens_tmpl, yr),
  dist_unpaved = function(yr) sprintf(unpaved_dist_tmpl, yr),
  dens_paved   = function(yr) paved_dens_file,
  dist_paved   = function(yr) paved_dist_file
)

predictor_order_all <- c(
  "dist_paved", "dist_unpaved", "dist_mines",
  "dens_paved", "dens_unpaved", "dens_mines"
)

predictor_order_out <- if (keep_unselected_as_na) {
  predictor_order_all
} else {
  selected_nonhab_predictors
}

selection_tag <- selected_nonhab_predictors |>
  sort() |>
  paste(collapse = "__")

## 3. DEFINE TIME LOG HELPERS ====================================================

### 3.1 Write log message --------------------------------------------------------

log_msg <- function(..., level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", ts, "] [", level, "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

### 3.2 Run a named step with timing ---------------------------------------------

run_step <- function(step_name, expr) {
  log_msg("START: ", step_name)
  t0 <- Sys.time()
  
  out <- tryCatch(
    eval.parent(substitute(expr)),
    error = function(e) {
      log_msg("FAILED: ", step_name, " | ", conditionMessage(e), level = "ERROR")
      
      trace_txt <- tryCatch(
        paste(capture.output(rlang::last_trace()), collapse = "\n"),
        error = function(e2) paste(capture.output(traceback()), collapse = "\n")
      )
      
      log_msg("TRACE:\n", trace_txt, level = "ERROR")
      stop(e)
    }
  )
  
  dt <- difftime(Sys.time(), t0, units = "secs")
  log_msg("END: ", step_name, " (", round(as.numeric(dt), 1), " s)")
  out
}

log_msg("Script started.")
log_msg("Log file: ", log_file)
log_msg("run_years = ", paste(run_years, collapse = ", "))
log_msg("selected_nonhab_predictors = ", paste(selected_nonhab_predictors, collapse = ", "))
log_msg("keep_unselected_as_na = ", keep_unselected_as_na)
log_msg("add_use_flags = ", add_use_flags)

## 4. DEFINE HELPER FUNCTIONS ====================================================

### 4.1 Make a safe herd ID ------------------------------------------------------

make_safe_id <- function(x) {
  x %>%
    stringr::str_to_lower() %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("&", "and") %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("^_+|_+$", "")
}

### 4.2 Read and prepare herd polygons -------------------------------------------

read_and_prep_herds <- function(qc_shp, on_shp, target_crs) {
  
  log_msg(" [HERDS] Reading QC and ON herd shapefiles ...")
  QC <- terra::vect(qc_shp)
  ON <- terra::vect(on_shp)
  
  log_msg(" [HERDS] Keeping only herd name fields ...")
  QC <- QC[, "NOM_POP", drop = FALSE]
  ON <- ON[, "RANGE_NAME", drop = FALSE]
  
  QC$herd_name <- QC$NOM_POP
  ON$herd_name <- ON$RANGE_NAME
  QC$region <- "QC"
  ON$region <- "ON"
  
  QC$NOM_POP <- NULL
  ON$RANGE_NAME <- NULL
  
  log_msg(" [HERDS] Validating geometries ...")
  QC <- try(terra::makeValid(QC), silent = TRUE) |> (\(x) if (inherits(x, "try-error")) QC else x)()
  ON <- try(terra::makeValid(ON), silent = TRUE) |> (\(x) if (inherits(x, "try-error")) ON else x)()
  
  log_msg(" [HERDS] Projecting to raster CRS ...")
  QCp <- terra::project(QC, target_crs)
  ONp <- terra::project(ON, target_crs)
  
  QCp$herd_id <- paste0("QC_", QCp$herd_name)
  ONp$herd_id <- paste0("ON_", ONp$herd_name)
  
  QCp$herd_id_safe <- paste0("QC_", make_safe_id(QCp$herd_name))
  ONp$herd_id_safe <- paste0("ON_", make_safe_id(ONp$herd_name))
  
  herds <- rbind(QCp, ONp)
  herds$poly_id <- seq_len(nrow(herds))
  
  log_msg(" [HERDS] Total herds loaded: ", nrow(herds))
  herds
}

### 4.3 Find the 8 composite habitat tiles for one year --------------------------

get_tiles_for_year <- function(year, hab_dir) {
  pat <- paste0("^compositeHabitat_", year, "_tile[1-8]\\.tif$")
  f <- list.files(hab_dir, pattern = pat, full.names = TRUE)
  f <- sort(f)
  
  if (length(f) != 8) {
    stop(
      "Expected 8 tiles for year ", year, " but found ", length(f), ":\n",
      paste(f, collapse = "\n")
    )
  }
  
  f
}

### 4.4 Build one habitat mosaic -------------------------------------------------

# Uses a VRT first. If that fails, writes a real mosaic into temp/.

build_year_mosaic <- function(tile_files, year, temp_dir) {
  vrt_file <- file.path(temp_dir, paste0("compositeHabitat_", year, "_mosaic.vrt"))
  
  log_msg(" [MOSAIC] Trying VRT for year ", year, " ...")
  r <- try(terra::vrt(tile_files, filename = vrt_file, overwrite = TRUE), silent = TRUE)
  
  if (inherits(r, "try-error")) {
    log_msg(" [MOSAIC] VRT failed. Falling back to terra::mosaic().", level = "WARN")
    
    out_tif <- file.path(temp_dir, paste0("compositeHabitat_", year, "_mosaic.tif"))
    s <- terra::sprc(tile_files)
    r <- terra::mosaic(s, filename = out_tif, overwrite = TRUE)
  }
  
  r
}

### 4.5 Extract habitat proportions for one year ---------------------------------

habitat_proportions_for_year <- function(r, herds, year, hab_vals, hab_codes, touches = FALSE) {
  
  log_msg(" [EXTRACT] Counting habitat categories for year ", year, " ...")
  
  counts <- terra::extract(
    r, herds,
    fun = table,
    na.rm = TRUE,
    touches = touches
  )
  
  counts$ID <- as.integer(counts$ID)
  
  nm <- names(counts)
  suffix_num <- stringr::str_extract(nm, "(?<=_)[0-9]+$|(?<=\\.)[0-9]+$|^[0-9]+$")
  
  cat_cols <- nm[!is.na(suffix_num) & suffix_num %in% as.character(hab_vals)]
  
  new_names <- nm
  new_names[match(cat_cols, nm)] <- suffix_num[match(cat_cols, nm)]
  names(counts) <- new_names
  
  for (v in as.character(hab_vals)) {
    if (!v %in% names(counts)) counts[[v]] <- 0
  }
  
  for (v in as.character(hab_vals)) {
    counts[[v]][is.na(counts[[v]])] <- 0
  }
  
  totals <- rowSums(counts[, as.character(hab_vals), drop = FALSE])
  
  props <- counts
  for (v in as.character(hab_vals)) {
    props[[v]] <- ifelse(totals > 0, props[[v]] / totals, NA_real_)
  }
  
  herd_key <- data.frame(
    ID = as.integer(herds$poly_id),
    herd_id = herds$herd_id,
    herd_id_safe = herds$herd_id_safe,
    herd = herds$herd_name,
    region = herds$region,
    stringsAsFactors = FALSE
  )
  
  props <- props %>%
    left_join(herd_key, by = "ID") %>%
    mutate(year = year) %>%
    select(region, herd, herd_id, herd_id_safe, year, all_of(as.character(hab_vals)))
  
  rename_map <- setNames(as.character(hab_vals), hab_codes)
  rename_map <- rename_map[rename_map %in% names(props)]
  
  props <- props %>% rename(!!!rename_map)
  
  props
}

### 4.6 Extract mean values from one non-habitat raster --------------------------

extract_mean_layer <- function(r_path, layer_name, herds, touches = FALSE) {
  
  if (!file.exists(r_path)) {
    stop("Missing raster file: ", r_path)
  }
  
  log_msg("   [RASTER] Loading: ", r_path)
  r <- terra::rast(r_path)
  names(r) <- layer_name
  
  out <- terra::extract(
    r, herds,
    fun = mean,
    na.rm = TRUE,
    touches = touches
  )
  
  out$ID <- as.integer(out$ID)
  
  if (!layer_name %in% names(out)) {
    val_cols <- setdiff(names(out), "ID")
    if (length(val_cols) != 1) {
      stop("Unexpected extract output columns: ", paste(names(out), collapse = ", "))
    }
    names(out)[names(out) == val_cols] <- layer_name
  }
  
  rm(r)
  gc()
  
  out
}

### 4.7 Build a representative non-habitat CRS -----------------------------------

detect_nonhab_crs <- function(years, selected_predictors, predictor_defs) {
  
  for (yr in years) {
    for (nm in selected_predictors) {
      p <- predictor_defs[[nm]](yr)
      if (file.exists(p)) {
        return(terra::crs(terra::rast(p)))
      }
    }
  }
  
  stop("Could not find any selected non-habitat raster to read a CRS from.")
}

## 5. EXTRACT HABITAT PROPORTIONS ================================================

### 5.1 Determine habitat raster CRS ---------------------------------------------

hab_crs <- run_step("STEP 5.1 — Determine habitat raster CRS", {
  sample_tile <- get_tiles_for_year(run_years[1], hab_dir)[1]
  terra::crs(terra::rast(sample_tile))
})

### 5.2 Read and prepare herds for habitat extraction ----------------------------

herds_hab <- run_step("STEP 5.2 — Read and prepare herds for habitat extraction", {
  read_and_prep_herds(qc_shp, on_shp, hab_crs)
})

### 5.3 Simplify habitat herd polygons -------------------------------------------

herds_hab <- run_step("STEP 5.3 — Simplify habitat herd polygons", {
  terra::simplifyGeom(
    herds_hab,
    tolerance = herd_simplify_tolerance_m,
    preserveTopology = TRUE
  )
})

### 5.4 Loop through years and extract habitat proportions -----------------------

all_years_long <- list()

run_step("STEP 5.4 — Loop through years for habitat proportions", {
  
  for (yr in run_years) {
    
    run_step(paste0("YEAR ", yr, " — habitat mosaic + extract"), {
      
      log_msg(" [FILES] Collecting habitat tiles for year ", yr, " ...")
      tiles <- get_tiles_for_year(yr, hab_dir)
      
      r_year <- run_step(paste0(" YEAR ", yr, " — Build habitat mosaic"), {
        build_year_mosaic(tiles, yr, temp_dir)
      })
      
      props_yr <- run_step(paste0(" YEAR ", yr, " — Extract habitat proportions"), {
        habitat_proportions_for_year(
          r = r_year,
          herds = herds_hab,
          year = yr,
          hab_vals = hab_vals,
          hab_codes = hab_codes,
          touches = habitat_touches
        )
      })
      
      all_years_long[[as.character(yr)]] <- props_yr
      
      rm(r_year)
      gc()
    })
  }
})

### 5.5 Build long and wide habitat tables ---------------------------------------

props_long <- run_step("STEP 5.5 — Combine habitat results", {
  bind_rows(all_years_long)
})

props_wide <- run_step("STEP 5.6 — Build habitat wide matrix", {
  props_long %>%
    mutate(herd_year = paste0(herd_id_safe, "_", year)) %>%
    select(herd_year, all_of(hab_codes)) %>%
    arrange(herd_year)
})

## 6. EXTRACT NON-HABITAT MEANS ==================================================

### 6.1 Determine non-habitat raster CRS -----------------------------------------

nonhab_crs <- run_step("STEP 6.1 — Determine non-habitat raster CRS", {
  detect_nonhab_crs(run_years, selected_nonhab_predictors, predictor_defs)
})

### 6.2 Read and prepare herds for non-habitat extraction ------------------------

herds_nonhab <- run_step("STEP 6.2 — Read and prepare herds for non-habitat extraction", {
  read_and_prep_herds(qc_shp, on_shp, nonhab_crs)
})

### 6.3 Simplify non-habitat herd polygons ---------------------------------------

herds_nonhab <- run_step("STEP 6.3 — Simplify non-habitat herd polygons", {
  terra::simplifyGeom(
    herds_nonhab,
    tolerance = herd_simplify_tolerance_m,
    preserveTopology = TRUE
  )
})

### 6.4 Build herd key for non-habitat joins -------------------------------------

herd_key_nonhab <- data.frame(
  ID = as.integer(herds_nonhab$poly_id),
  herd_id_safe = herds_nonhab$herd_id_safe,
  herd_id = herds_nonhab$herd_id,
  herd = herds_nonhab$herd_name,
  region = herds_nonhab$region,
  stringsAsFactors = FALSE
)

### 6.5 Loop through years and extract selected non-habitat predictors -----------

nonhab_long <- run_step("STEP 6.5 — Extract non-habitat predictor means", {
  
  all_years <- list()
  
  for (yr in run_years) {
    
    yr_ckpt <- file.path(
      checkpoint_dir,
      paste0("nonhab_means_", yr, "_", make_safe_id(selection_tag), ".csv")
    )
    
    if (resume_nonhab_checkpoints && file.exists(yr_ckpt)) {
      log_msg(" [RESUME] Found checkpoint for year ", yr, ": ", yr_ckpt)
      yr_df <- readr::read_csv(yr_ckpt, show_col_types = FALSE)
      all_years[[as.character(yr)]] <- yr_df
      next
    }
    
    yr_df <- run_step(paste0("YEAR ", yr, " — non-habitat extraction"), {
      
      selected_paths <- lapply(selected_nonhab_predictors, function(nm) predictor_defs[[nm]](yr))
      names(selected_paths) <- selected_nonhab_predictors
      
      missing <- names(selected_paths)[!file.exists(unlist(selected_paths))]
      
      if (length(missing) > 0) {
        stop(
          "Missing selected raster files for year ", yr, ": ",
          paste(missing, collapse = ", "),
          "\nPaths:\n",
          paste(unlist(selected_paths), collapse = "\n")
        )
      }
      
      layer_dfs <- lapply(names(selected_paths), function(nm) {
        extract_mean_layer(
          r_path = selected_paths[[nm]],
          layer_name = nm,
          herds = herds_nonhab,
          touches = nonhab_touches
        )
      })
      names(layer_dfs) <- names(selected_paths)
      
      yr_out <- Reduce(function(a, b) dplyr::left_join(a, b, by = "ID"), layer_dfs)
      
      yr_out <- yr_out %>%
        left_join(herd_key_nonhab, by = "ID") %>%
        mutate(year = yr)
      
      if (keep_unselected_as_na) {
        for (nm in predictor_order_all) {
          if (!nm %in% names(yr_out)) {
            yr_out[[nm]] <- unselected_fill_value
          }
        }
      }
      
      if (add_use_flags) {
        for (nm in predictor_order_all) {
          flag <- paste0("use_", nm)
          yr_out[[flag]] <- nm %in% selected_nonhab_predictors
        }
      }
      
      keep_cols <- c(
        "region", "herd", "herd_id", "herd_id_safe", "year",
        predictor_order_out
      )
      
      if (add_use_flags) {
        keep_cols <- c(keep_cols, paste0("use_", predictor_order_out))
      }
      
      yr_out %>% select(all_of(keep_cols))
    })
    
    run_step(paste0("YEAR ", yr, " — save non-habitat checkpoint"), {
      readr::write_csv(yr_df, yr_ckpt)
      log_msg(" [CHECKPOINT] Saved: ", yr_ckpt)
    })
    
    all_years[[as.character(yr)]] <- yr_df
  }
  
  bind_rows(all_years)
})

### 6.6 Build wide non-habitat matrix --------------------------------------------

nonhab_wide <- run_step("STEP 6.6 — Build non-habitat wide matrix", {
  
  out <- nonhab_long %>%
    mutate(herd_year = paste0(herd_id_safe, "_", year)) %>%
    select(herd_year, all_of(predictor_order_out)) %>%
    arrange(herd_year)
  
  if (add_use_flags) {
    flag_cols <- paste0("use_", predictor_order_out)
    
    out <- nonhab_long %>%
      mutate(herd_year = paste0(herd_id_safe, "_", year)) %>%
      select(herd_year, all_of(predictor_order_out), all_of(flag_cols)) %>%
      arrange(herd_year)
  }
  
  out
})

## 7. BUILD FINAL MODEL MATRIX ===================================================

### 7.1 Merge habitat and non-habitat values -------------------------------------

final_matrix <- run_step("STEP 7.1 — Merge habitat and non-habitat matrices", {
  props_wide %>%
    left_join(nonhab_wide, by = "herd_year") %>%
    arrange(herd_year)
})

## 8. SAVE OUTPUTS ===============================================================

### 8.1 Save habitat outputs -----------------------------------------------------

run_step("STEP 8.1 — Save habitat outputs", {
  readr::write_csv(
    props_long,
    file.path(habitat_out_dir, "habitat_proportions_long.csv")
  )
  
  readr::write_csv(
    props_wide,
    file.path(habitat_out_dir, "habitat_proportions_matrix_wide.csv")
  )
})

### 8.2 Save non-habitat outputs -------------------------------------------------

run_step("STEP 8.2 — Save non-habitat outputs", {
  readr::write_csv(
    nonhab_long,
    file.path(model_out_dir, "nonhab_means_long.csv")
  )
  
  readr::write_csv(
    nonhab_wide,
    file.path(model_out_dir, "nonhab_means_matrix_wide.csv")
  )
})

### 8.3 Save final model matrix --------------------------------------------------

run_step("STEP 8.3 — Save final model matrix", {
  readr::write_csv(
    final_matrix,
    file.path(model_out_dir, "final_model_matrix_wide.csv")
  )
})

## 9. SANITY CHECKS ==============================================================

### 9.1 Check row counts ---------------------------------------------------------

run_step("STEP 9.1 — Run sanity checks", {
  
  expected_rows <- length(run_years) * nrow(herds_hab)
  
  stopifnot(nrow(props_long) == expected_rows)
  stopifnot(nrow(props_wide) == expected_rows)
  stopifnot(nrow(final_matrix) == expected_rows)
  
  stopifnot(!anyDuplicated(props_wide$herd_year))
  stopifnot(!anyDuplicated(nonhab_wide$herd_year))
  stopifnot(!anyDuplicated(final_matrix$herd_year))
  
  stopifnot(all(hab_codes %in% names(props_wide)))
  stopifnot(all(predictor_order_out %in% names(final_matrix)))
})

log_msg("Checks passed.")
log_msg("Script completed successfully.")

