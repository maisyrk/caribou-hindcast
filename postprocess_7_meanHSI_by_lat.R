# postprocess_7_meanHSI_by_lat.R ----------

# Description:
# This script creates a figure that summarizes mean HSI by latitude for 
# 1985, 2000, and 2020 and overlays a shaded area showing how many raster cells 
# are touched by the cleaned timber allocation line in each latitude bin.
#
# Required inputs:
# - HSI_output/HSI_1985_masked.tif
# - HSI_output/HSI_2000_masked.tif
# - HSI_output/HSI_2020_masked.tif
# - map_data/forest_limit/LIM_FOR_ATT_2018.shp
# - map_data/forest_limit/northern_extent_managed_forest.shp
#
# Expected outputs:
# - HSI_output/mean_HSI_by_latitude_bin_<...>.csv
# - HSI_output/timber_line_contact_by_latitude_bin_<...>.csv
# - HSI_output/timber_limit_clean_combined.gpkg
# - HSI_output/mean_HSI_by_latitude_bin<...>_timber_line_contact_area.png
# - HSI_output/mean_HSI_by_latitude_bin<...>_timber_block.png

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
    "Open the packaged project folder as the working directory, ",
    "or run this script from that project root."
  )
}

project_root <- get_project_root()
setwd(project_root)

### 1.2 Load required packages ---------------------------------------------------

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(ggplot2)
  library(dplyr)
  library(rnaturalearth)
})

### 1.3 User inputs --------------------------------------------------------------

# HSI rasters
hsi_1985_path <- file.path(project_root, "HSI_output", "HSI_1985_masked.tif")
hsi_2000_path <- file.path(project_root, "HSI_output", "HSI_2000_masked.tif")
hsi_2020_path <- file.path(project_root, "HSI_output", "HSI_2020_masked.tif")

# Timber allocation limits
timber_limit_qc_path <- file.path(project_root, "map_data", "forest_limit", "LIM_FOR_ATT_2018.shp")
timber_limit_on_path <- file.path(project_root, "map_data", "forest_limit", "northern_extent_managed_forest.shp")

# Latitude bin width:
# 1.0 = full degree
# 0.1 = every tenth of a degree
lat_bin_width <- 0.1

# Aggregation factor for HSI mean extraction:
# 1  = no aggregation
# 10 = lighter test run
# 20 = lighter test run
agg_factor <- 1

# Number of rows per HSI processing chunk
chunk_nrows <- 64

# Reuse existing HSI summary CSV?
reuse_existing_hsi_csv <- TRUE

# Template raster for timber-line cell contact counting
hsi_template_path <- hsi_2020_path

# Save the cleaned combined timber line?
save_clean_timber_line <- TRUE

# Timber-line 'densification factor', for simplifying the line geometry
densify_factor <- 0.5

### 1.4 Define output files ------------------------------------------------------

lat_bin_label <- gsub("\\.", "p", as.character(lat_bin_width))

out_csv <- file.path(
  project_root,
  "HSI_output",
  paste0(
    "mean_HSI_by_latitude_bin_",
    lat_bin_label,
    "_agg", agg_factor,
    ".csv"
  )
)

out_contact_csv <- file.path(
  project_root,
  "HSI_output",
  paste0(
    "timber_line_contact_by_latitude_bin_",
    lat_bin_label,
    "_template_",
    tools::file_path_sans_ext(basename(hsi_template_path)),
    ".csv"
  )
)

out_clean_line_gpkg <- file.path(
  project_root,
  "HSI_output",
  "timber_limit_clean_combined.gpkg"
)

out_png <- file.path(
  project_root,
  "HSI_output",
  paste0(
    "mean_HSI_by_latitude_bin",
    lat_bin_label,
    "_agg", agg_factor,
    "_timber_line_contact_area.png"
  )
)

out_png_block <- file.path(
  project_root,
  "HSI_output",
  paste0(
    "mean_HSI_by_latitude_bin",
    lat_bin_label,
    "_agg", agg_factor,
    "_timber_block.png"
  )
)

### 1.5 Create output folder and terra options -----------------------------------

dir.create(file.path(project_root, "HSI_output"), recursive = TRUE, showWarnings = FALSE)

terraOptions(
  memfrac  = 0.7,
  progress = 3,
  tempdir  = tempdir()
)

## 2. DEFINE HELPER FUNCTIONS ====================================================

### 2.1 Process one raster by latitude -------------------------------------------

process_hsi_by_latitude <- function(r_path,
                                    year_label,
                                    lat_bin_width,
                                    agg_factor,
                                    chunk_nrows) {
  
  cat("\n----------------------------------------\n")
  cat("Processing:", year_label, "\n")
  cat("Raster:", r_path, "\n")
  cat("----------------------------------------\n")
  
  r <- rast(r_path)
  
  if (nlyr(r) != 1) {
    stop("Raster must contain exactly one layer: ", r_path)
  }
  
  use_weights <- (agg_factor > 1)
  
  # Optional aggregation for lighter runs.
  if (use_weights) {
    cat("Aggregating raster by factor =", agg_factor, "...\n")
    
    r_mean <- aggregate(
      r,
      fact      = agg_factor,
      fun       = mean,
      na.rm     = TRUE,
      filename  = tempfile(pattern = paste0("agg_mean_", year_label, "_"), fileext = ".tif"),
      overwrite = TRUE
    )
    
    valid_mask <- ifel(is.na(r), 0, 1)
    
    r_n <- aggregate(
      valid_mask,
      fact      = agg_factor,
      fun       = sum,
      na.rm     = TRUE,
      filename  = tempfile(pattern = paste0("agg_n_", year_label, "_"), fileext = ".tif"),
      overwrite = TRUE
    )
    
  } else {
    cat("Using native resolution (no aggregation)...\n")
    r_mean <- r
    r_n <- NULL
  }
  
  nr <- nrow(r_mean)
  nc <- ncol(r_mean)
  xs <- xFromCol(r_mean, 1:nc)
  
  row_starts <- seq(1, nr, by = chunk_nrows)
  n_chunks <- length(row_starts)
  
  chunk_results <- vector("list", n_chunks)
  
  readStart(r_mean)
  if (use_weights) readStart(r_n)
  
  on.exit({
    readStop(r_mean)
    if (use_weights) readStop(r_n)
  }, add = TRUE)
  
  for (i in seq_along(row_starts)) {
    
    row_start <- row_starts[i]
    nrows_i   <- min(chunk_nrows, nr - row_start + 1)
    rows_i    <- row_start:(row_start + nrows_i - 1)
    
    cat("  HSI chunk", i, "of", n_chunks, "\n")
    
    vals <- readValues(r_mean, row = row_start, nrows = nrows_i, mat = FALSE)
    
    if (use_weights) {
      wts  <- readValues(r_n, row = row_start, nrows = nrows_i, mat = FALSE)
      keep <- !is.na(vals) & (wts > 0)
    } else {
      keep <- !is.na(vals)
    }
    
    if (!any(keep)) next
    
    valid_idx <- which(keep)
    
    row_in_chunk <- ((valid_idx - 1) %/% nc) + 1
    col_idx      <- ((valid_idx - 1) %%  nc) + 1
    
    ys_valid <- yFromRow(r_mean, rows_i[row_in_chunk])
    
    # Convert y to latitude.
    if (is.lonlat(r_mean)) {
      lat <- ys_valid
    } else {
      xs_valid <- xs[col_idx]
      coords_valid <- cbind(xs_valid, ys_valid)
      coords_ll <- terra::project(coords_valid, from = crs(r_mean), to = "EPSG:4326")
      lat <- coords_ll[, 2]
    }
    
    lat_lower <- floor(lat / lat_bin_width) * lat_bin_width
    
    if (use_weights) {
      d <- data.frame(
        lat_lower = lat_lower,
        sum_hsi   = vals[keep] * wts[keep],
        n_cells   = wts[keep]
      )
    } else {
      d <- data.frame(
        lat_lower = lat_lower,
        sum_hsi   = vals[keep],
        n_cells   = 1L
      )
    }
    
    chunk_results[[i]] <- aggregate(
      cbind(sum_hsi, n_cells) ~ lat_lower,
      data = d,
      FUN = sum
    )
  }
  
  chunk_results <- Filter(Negate(is.null), chunk_results)
  
  if (length(chunk_results) == 0) {
    stop("No valid data found in raster: ", r_path)
  }
  
  combined <- do.call(rbind, chunk_results)
  
  final <- aggregate(
    cbind(sum_hsi, n_cells) ~ lat_lower,
    data = combined,
    FUN  = sum
  )
  
  final$mean_hsi <- final$sum_hsi / final$n_cells
  final$latitude <- final$lat_lower + (lat_bin_width / 2)
  final$year     <- factor(year_label, levels = c("1985", "2000", "2020"))
  
  final <- final[order(final$latitude), c("year", "latitude", "mean_hsi", "n_cells")]
  
  final
}

### 2.2 Convert geometry to simple lines -----------------------------------------

as_lines_sf <- function(x) {
  
  x <- st_make_valid(x)
  g <- st_geometry(x)
  gt <- unique(as.character(st_geometry_type(g, by_geometry = TRUE)))
  
  # If polygons, convert to boundaries.
  if (all(gt %in% c("POLYGON", "MULTIPOLYGON"))) {
    g <- st_boundary(g)
  }
  
  g_line <- suppressWarnings(st_collection_extract(g, "LINESTRING"))
  
  if (length(g_line) == 0) {
    g_line <- suppressWarnings(st_cast(g, "MULTILINESTRING", warn = FALSE))
    g_line <- suppressWarnings(st_cast(g_line, "LINESTRING", warn = FALSE))
  }
  
  g_line <- g_line[!st_is_empty(g_line)]
  
  if (length(g_line) == 0) {
    stop("Could not convert geometry to lines.")
  }
  
  st_sf(geometry = g_line)
}

### 2.3 Get Ontario and Quebec polygons from Natural Earth -----------------------

get_on_qc_from_ne <- function(target_crs) {
  
  provinces <- ne_states(country = "canada", returnclass = "sf") |>
    filter(name_en %in% c("Ontario", "Quebec")) |>
    st_transform(target_crs) |>
    st_make_valid()
  
  qc <- provinces |>
    filter(name_en == "Quebec") |>
    summarise(geometry = st_union(geometry), .groups = "drop") |>
    st_make_valid()
  
  on <- provinces |>
    filter(name_en == "Ontario") |>
    summarise(geometry = st_union(geometry), .groups = "drop") |>
    st_make_valid()
  
  list(qc = qc, on = on)
}

### 2.4 Clip one timber line to one province polygon -----------------------------

clip_line_to_province <- function(line_path, province_sf, target_crs) {
  
  ln <- st_read(line_path, quiet = TRUE)
  ln <- st_make_valid(ln)
  ln <- st_transform(ln, target_crs)
  ln <- st_make_valid(ln)
  ln <- as_lines_sf(ln)
  
  province_sf <- st_transform(province_sf, target_crs)
  province_sf <- st_make_valid(province_sf)
  
  ln <- suppressWarnings(st_crop(ln, st_bbox(province_sf)))
  clipped <- suppressWarnings(st_intersection(ln, province_sf))
  clipped <- clipped[!st_is_empty(clipped), , drop = FALSE]
  
  if (nrow(clipped) == 0) {
    stop("No timber line remained after clipping: ", line_path)
  }
  
  as_lines_sf(clipped)
}

### 2.5 Build one clean combined timber line -------------------------------------

build_clean_timber_line <- function(qc_line_path,
                                    on_line_path,
                                    qc_sf,
                                    on_sf,
                                    target_crs,
                                    bbox_poly = NULL) {
  
  cat("\n----------------------------------------\n")
  cat("Building clean combined timber line\n")
  cat("----------------------------------------\n")
  
  qc_line <- clip_line_to_province(qc_line_path, qc_sf, target_crs)
  on_line <- clip_line_to_province(on_line_path, on_sf, target_crs)
  
  combined_geom <- c(st_geometry(qc_line), st_geometry(on_line))
  clean_geom <- st_union(st_combine(combined_geom))
  clean_geom <- st_make_valid(clean_geom)
  
  clean_line <- st_sf(geometry = clean_geom)
  clean_line <- as_lines_sf(clean_line)
  
  if (!is.null(bbox_poly)) {
    clean_line <- suppressWarnings(st_crop(clean_line, st_bbox(bbox_poly)))
    clean_line <- clean_line[!st_is_empty(clean_line), , drop = FALSE]
  }
  
  clean_line
}

### 2.6 Build the reference raster object ----------------------------------------

build_reference_grid <- function(r_path) {
  cat("\n----------------------------------------\n")
  cat("Using one template raster for timber-line cell contacts\n")
  cat("----------------------------------------\n")
  
  rast(r_path)
}

### 2.7 Summarize timber-line contact counts by latitude -------------------------

summarize_line_contact_by_latitude <- function(ref_raster,
                                               clean_line_sf,
                                               lat_bin_width = 0.1,
                                               densify_factor = 0.5) {
  
  cat("\n----------------------------------------\n")
  cat("Summarizing timber-line contact by latitude\n")
  cat("----------------------------------------\n")
  
  clean_line_sf <- st_transform(clean_line_sf, st_crs(crs(ref_raster)))
  clean_line_sf <- st_make_valid(clean_line_sf)
  
  if (is.lonlat(ref_raster)) {
    stop("This script expects a projected raster CRS, not lon/lat.")
  }
  
  step_len <- min(res(ref_raster), na.rm = TRUE) * densify_factor
  
  if (!is.finite(step_len) || step_len <= 0) {
    stop("Could not determine a valid densification step length.")
  }
  
  dense_line <- st_segmentize(clean_line_sf, dfMaxLength = step_len)
  
  coords <- st_coordinates(dense_line)
  
  if (nrow(coords) == 0) {
    stop("No coordinates found after densifying the timber line.")
  }
  
  xy_line <- coords[, c("X", "Y"), drop = FALSE]
  
  touched_cells <- cellFromXY(ref_raster, xy_line)
  touched_cells <- unique(touched_cells[!is.na(touched_cells)])
  
  if (length(touched_cells) == 0) {
    stop("No raster cells were touched by the timber line.")
  }
  
  touched_vals <- ref_raster[touched_cells]
  touched_cells <- touched_cells[!is.na(touched_vals)]
  
  if (length(touched_cells) == 0) {
    stop("No valid study-area raster cells were touched by the timber line.")
  }
  
  xy_cells <- xyFromCell(ref_raster, touched_cells)
  
  if (is.lonlat(ref_raster)) {
    lat <- xy_cells[, 2]
  } else {
    xy_ll <- terra::project(xy_cells, from = crs(ref_raster), to = "EPSG:4326")
    lat <- xy_ll[, 2]
  }
  
  lat_lower <- floor(lat / lat_bin_width) * lat_bin_width
  
  out <- aggregate(
    x  = list(touched_cells = rep(1L, length(lat_lower))),
    by = list(lat_lower = lat_lower),
    FUN = sum
  )
  
  out$latitude <- out$lat_lower + (lat_bin_width / 2)
  out <- out[order(out$latitude), c("latitude", "touched_cells")]
  
  out
}

## 3. OPTIONAL GEOMETRY CHECKS ===================================================

### 3.1 Compare raster geometry across years -------------------------------------

r1985 <- rast(hsi_1985_path)
r2000 <- rast(hsi_2000_path)
r2020 <- rast(hsi_2020_path)

if (!compareGeom(r1985, r2000, stopOnError = FALSE)) {
  warning("1985 and 2000 rasters do not have identical geometry. Double-check alignment.")
}

if (!compareGeom(r1985, r2020, stopOnError = FALSE)) {
  warning("1985 and 2020 rasters do not have identical geometry. Double-check alignment.")
}

## 4. CREATE OR RELOAD HSI LATITUDE SUMMARY ======================================

### 4.1 Reuse existing summary if requested --------------------------------------

if (reuse_existing_hsi_csv && file.exists(out_csv)) {
  
  cat("\nReading previously saved latitude summary:\n", out_csv, "\n")
  plot_df <- read.csv(out_csv, stringsAsFactors = FALSE)
  plot_df$year <- factor(plot_df$year, levels = c("1985", "2000", "2020"))
  
} else {
  
  ### 4.1.1 Process 1985 ---------------------------------------------------------
  
  res_1985 <- process_hsi_by_latitude(
    r_path        = hsi_1985_path,
    year_label    = "1985",
    lat_bin_width = lat_bin_width,
    agg_factor    = agg_factor,
    chunk_nrows   = chunk_nrows
  )
  
  ### 4.1.2 Process 2000 ---------------------------------------------------------
  
  res_2000 <- process_hsi_by_latitude(
    r_path        = hsi_2000_path,
    year_label    = "2000",
    lat_bin_width = lat_bin_width,
    agg_factor    = agg_factor,
    chunk_nrows   = chunk_nrows
  )
  
  ### 4.1.3 Process 2020 ---------------------------------------------------------
  
  res_2020 <- process_hsi_by_latitude(
    r_path        = hsi_2020_path,
    year_label    = "2020",
    lat_bin_width = lat_bin_width,
    agg_factor    = agg_factor,
    chunk_nrows   = chunk_nrows
  )
  
  ### 4.1.4 Combine and save -----------------------------------------------------
  
  plot_df <- rbind(res_1985, res_2000, res_2020)
  
  write.csv(plot_df, out_csv, row.names = FALSE)
  cat("\nSaved extracted latitude summaries to:\n", out_csv, "\n")
}

## 5. BUILD CLEAN TIMBER LINE ====================================================

### 5.1 Build reference grid and extent polygon ----------------------------------

ref_raster <- build_reference_grid(hsi_template_path)

ref_bbox_poly <- st_as_sfc(
  st_bbox(
    c(
      xmin = xmin(ref_raster),
      ymin = ymin(ref_raster),
      xmax = xmax(ref_raster),
      ymax = ymax(ref_raster)
    ),
    crs = st_crs(crs(ref_raster))
  )
)

### 5.2 Get Ontario and Quebec province boundaries -------------------------------

prov_list <- get_on_qc_from_ne(target_crs = st_crs(crs(ref_raster)))
qc <- prov_list$qc
on <- prov_list$on

### 5.3 Build the cleaned combined line ------------------------------------------

timber_line_clean <- build_clean_timber_line(
  qc_line_path = timber_limit_qc_path,
  on_line_path = timber_limit_on_path,
  qc_sf        = qc,
  on_sf        = on,
  target_crs   = st_crs(crs(ref_raster)),
  bbox_poly    = ref_bbox_poly
)

### 5.4 Optionally save the cleaned line -----------------------------------------

if (save_clean_timber_line) {
  if (file.exists(out_clean_line_gpkg)) file.remove(out_clean_line_gpkg)
  st_write(timber_line_clean, out_clean_line_gpkg, quiet = TRUE)
  cat("\nSaved clean timber line to:\n", out_clean_line_gpkg, "\n")
}

## 6. SUMMARIZE TIMBER-LINE CONTACTS =============================================

### 6.1 Count touched cells by latitude ------------------------------------------

contact_df <- summarize_line_contact_by_latitude(
  ref_raster      = ref_raster,
  clean_line_sf   = timber_line_clean,
  lat_bin_width   = lat_bin_width,
  densify_factor  = densify_factor
)

timber_block <- data.frame(
  xmin = min(contact_df$latitude, na.rm = TRUE) - lat_bin_width / 2,
  xmax = max(contact_df$latitude, na.rm = TRUE) + lat_bin_width / 2,
  label = "Timber allocation limit"
)

### 6.2 Save contact summary -----------------------------------------------------

write.csv(contact_df, out_contact_csv, row.names = FALSE)
cat("\nSaved timber-line contact summary to:\n", out_contact_csv, "\n")

## 7. BUILD MAIN CONTACT-AREA PLOT ===============================================

### 7.1 Prepare y-axis scaling ---------------------------------------------------

p <- ggplot()

y_rng  <- range(plot_df$mean_hsi, na.rm = TRUE)
y_min  <- y_rng[1]
y_max  <- y_rng[2]
y_span <- y_max - y_min

if (y_span == 0) y_span <- 0.05

max_contact <- max(contact_df$touched_cells, na.rm = TRUE)
if (!is.finite(max_contact) || max_contact <= 0) max_contact <- 1

contact_df$contact_y <- y_min + (contact_df$touched_cells / max_contact) * y_span

### 7.2 Build plot ---------------------------------------------------------------

p <- p +
  geom_ribbon(
    data = contact_df,
    aes(x = latitude, ymin = y_min, ymax = contact_y),
    inherit.aes = FALSE,
    fill = "grey65",
    alpha = 0.3,
    color = NA
  ) +
  geom_line(
    data = plot_df,
    aes(x = latitude, y = mean_hsi, color = year, group = year),
    linewidth = 1,
    lineend = "round"
  ) +
  scale_color_manual(
    values = c(
      "1985" = "#1b9e77",
      "2000" = "#d95f02",
      "2020" = "#7570b3"
    ),
    name = "Year"
  ) +
  scale_x_continuous(
    name = "Latitude (°)",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    name = "Mean HSI",
    expand = expansion(mult = c(0.02, 0.04)),
    sec.axis = sec_axis(
      trans = ~ ((. - y_min) / y_span) * max_contact,
      name = "Overlap with Northern Timber Limit (# pixels)"
    )
  ) +
  guides(
    color = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 1,
      byrow = TRUE
    )
  ) +
  labs(
    x = "Latitude (°)",
    y = "Mean HSI"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.y.right = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    legend.position = c(0.98, 0.04),
    legend.justification = c(1, 0),
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = scales::alpha("white", 0.5), color = NA),
    legend.key = element_blank(),
    plot.margin = margin(8, 8, 8, 8)
  )

print(p)

## 8. BUILD ALTERNATIVE BLOCK PLOT ===============================================

### 8.1 Build block-based plot ---------------------------------------------------

p2 <- ggplot() +
  geom_rect(
    data = timber_block,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = label),
    inherit.aes = FALSE,
    alpha = 0.3,
    color = NA
  ) +
  geom_line(
    data = plot_df,
    aes(x = latitude, y = mean_hsi, color = year, group = year),
    linewidth = 1,
    lineend = "round"
  ) +
  scale_color_manual(
    values = c(
      "1985" = "#1b9e77",
      "2000" = "#d95f02",
      "2020" = "#7570b3"
    ),
    name = "Year"
  ) +
  scale_fill_manual(
    values = c("Timber allocation limit" = "grey65"),
    labels = c("Timber allocation limit" = "Timber Allocation Northern Limit"),
    name = NULL
  ) +
  scale_x_continuous(
    name = "Latitude (°)",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    name = "Mean HSI",
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  guides(
    color = guide_legend(
      order = 1,
      title.position = "top",
      title.hjust = 0,
      ncol = 1,
      byrow = TRUE
    ),
    fill = guide_legend(
      order = 2,
      title.position = "top",
      title.hjust = 0,
      ncol = 1,
      byrow = TRUE
    )
  ) +
  labs(
    x = "Latitude (°)",
    y = "Mean HSI"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_blank(),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_text(face = "bold"),
    legend.background = element_rect(fill = scales::alpha("white", 0.5), color = NA),
    legend.key = element_blank(),
    legend.spacing.y = unit(2, "pt"),
    plot.margin = margin(8, 8, 8, 8)
  )

print(p2)

## 9. SAVE FIGURES ===============================================================

### 9.1 Save main plot -----------------------------------------------------------

ggsave(
  filename = out_png,
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)

cat("\nSaved figure to:\n", out_png, "\n")
cat("\nDone.\n")

### 9.2 Save alternative block plot ----------------------------------------------

ggsave(
  filename = out_png_block,
  plot = p2,
  width = 10,
  height = 6,
  dpi = 300
)

cat("\nSaved alternative block figure to:\n", out_png_block, "\n")

