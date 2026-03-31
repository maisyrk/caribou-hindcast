# postprocess_6_PCoA.R ----------

# Description:
# Run this script to create the habitat-only PCoA trajectory figure.
# It facets the plot by province and position relative to the timber line.
#
# Required inputs:
# - data/post_processing/model_matrix/final_model_matrix_wide.csv
# - tables/herd_timber_line_position_by_herd_combined_ON_QC.csv
#
# Expected outputs:
# - plots/pcoa_habitat_only_fixed_labels_Qc_ON_by_timber_position.png
#
# Notes:
# - Ordination is habitat only (Hellinger).
# - Road density is shown only through point size.
# - The herd position table comes from postprocess_5_herd_spatial_classification.R.
# - Actual ordination only starts at step 6. 

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
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(vegan)
  library(ggrepel)
  library(grid)
})

### 1.3 User inputs --------------------------------------------------------------

in_csv <- file.path(
  project_root,
  "data", "post_processing", "model_matrix",
  "final_model_matrix_wide.csv"
)

in_position_csv <- file.path(
  project_root,
  "tables",
  "herd_timber_line_position_by_herd_combined_ON_QC.csv"
)

years_keep <- c(1985, 2000, 2020)
label_year <- 2020

# Final habitat names used in the ordination / envfit
hab_cols <- c(
  "fire", "harv_0-5", "harv_6-20", "young_con",
  "old_con", "open_wood", "wet", "regen"
)

# Numeric habitat columns from older matrix versions
hab_num_map <- c(
  "1" = "fire",
  "2" = "harv_0-5",
  "3" = "harv_6-20",
  "4" = "young_con",
  "5" = "old_con",
  "6" = "open_wood",
  "7" = "wet",
  "8" = "regen"
)

# Habitat columns from the refactored extraction script
hab_refactored_map <- c(
  "fire"    = "fire",
  "harv_y"  = "harv_0-5",
  "harv_o"  = "harv_6-20",
  "conf_y"  = "young_con",
  "conf_o"  = "old_con",
  "open_w"  = "open_wood",
  "wetland" = "wet",
  "regen"   = "regen"
)

# All possible anthropogenic variables in the model matrix
anthro_cols <- c("dist_paved", "dist_unpaved", "dens_paved", "dens_unpaved", "dens_mines")

# Anthropogenic variable shown by point size
road_var <- "dens_unpaved"
road_legend_title <- "Density of unpaved roads"

# Optional envfit filtering
show_only_strong_vectors <- FALSE
r2_threshold <- 0.20

# Invisible obstacle points for ggrepel
n_path_pts  <- 12
n_arrow_pts <- 10

out_plot <- file.path(
  project_root,
  "plots",
  "pcoa_habitat_only_fixed_labels_Qc_ON_by_timber_position.png"
)

# Pretty labels for envfit vectors
envfit_label_map <- c(
  "old_con"   = "Old Con",
  "young_con" = "Young Con",
  "open_wood" = "Open W",
  "wet"       = "Wet",
  "regen"     = "Regen",
  "harv_0-5"  = "Harvest 0-5",
  "harv_6-20" = "Harvest 6-20",
  "fire"      = "Fire"
)

# Herd label -> herd code map
herd_code_map <- c(
  "Assinica" = "ASN",
  "Basse C Te Nord" = "BCN",
  "Caniapiscau" = "CAN",
  "Charlevoix" = "CVX",
  "Detour" = "DET",
  "Joir River" = "JR",
  "Lac Joseph" = "LJ",
  "Manicouagan" = "MAN",
  "Nottaway" = "NOT",
  "Outardes" = "OUT",
  "Pipmuacan" = "PIP",
  "Secteur Baie James" = "SBJ",
  "Secteur Matamec" = "SM",
  "T Miscamie" = "TEM",
  "Val D Or" = "VdO",
  "Berens" = "BER",
  "Brightsand" = "BRS",
  "Churchill" = "CHU",
  "Discontinuous Distribution" = "DD",
  "James Bay" = "JB",
  "Kesagami" = "KES",
  "Kinloch" = "KIN",
  "Lake Superior Coast" = "LSC",
  "Missisa" = "MIS",
  "Nipigon" = "NIP",
  "Ozhiski" = "OZH",
  "Pagwachuan" = "PAG",
  "Spirit" = "SPR",
  "Swan" = "SL",
  "Sydney" = "SYD"
)

# Herd codes to exclude
exclude_herd_codes <- c("JR", "LJ")

# Facet order
position_levels_raw  <- c("North (>75%)", "South (>75%)", "Line straddler")
position_levels_plot <- c("North (>75%)", "South (>75%)", "On the line")

### 1.4 Create output folder -----------------------------------------------------

dir.create(dirname(out_plot), recursive = TRUE, showWarnings = FALSE)

## 2. DEFINE HELPER FUNCTIONS ====================================================

### 2.1 Normalize a text key -----------------------------------------------------

normalize_key <- function(x) {
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT")
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "", x, perl = TRUE)
  x
}

### 2.2 Standardize habitat column names -----------------------------------------

# This allows the script to work with either:
# - older numbered habitat columns (1:8), or
# - the refactored extraction column names.

standardize_habitat_columns <- function(dat_raw, hab_cols, hab_num_map, hab_refactored_map) {
  
  if (all(hab_cols %in% names(dat_raw))) {
    return(dat_raw)
  }
  
  if (all(names(hab_num_map) %in% names(dat_raw))) {
    dat_raw <- dat_raw %>%
      rename(!!!setNames(names(hab_num_map), unname(hab_num_map)))
  }
  
  if (!all(hab_cols %in% names(dat_raw)) &&
      all(names(hab_refactored_map) %in% names(dat_raw))) {
    dat_raw <- dat_raw %>%
      rename(!!!setNames(names(hab_refactored_map), unname(hab_refactored_map)))
  }
  
  missing_hab <- setdiff(hab_cols, names(dat_raw))
  
  if (length(missing_hab) > 0) {
    stop(
      "The following required habitat columns are missing after renaming: ",
      paste(missing_hab, collapse = ", "),
      "\nAvailable columns are:\n",
      paste(names(dat_raw), collapse = ", ")
    )
  }
  
  dat_raw
}

## 3. LOAD AND CLEAN INPUT DATA ==================================================

### 3.1 Read the model matrix ----------------------------------------------------

dat_raw <- readr::read_csv(in_csv, show_col_types = FALSE)

### 3.2 Standardize habitat columns ----------------------------------------------

dat_raw <- standardize_habitat_columns(
  dat_raw = dat_raw,
  hab_cols = hab_cols,
  hab_num_map = hab_num_map,
  hab_refactored_map = hab_refactored_map
)

### 3.3 Drop anthropogenic columns that are all NA -------------------------------

anthro_present <- intersect(anthro_cols, names(dat_raw))

anthro_all_na <- anthro_present[
  vapply(dat_raw[anthro_present], function(x) all(is.na(x)), logical(1))
]

if (length(anthro_all_na) > 0) {
  message(
    "Dropping anthropogenic columns with all NA values: ",
    paste(anthro_all_na, collapse = ", ")
  )
  dat_raw <- dat_raw %>% select(-all_of(anthro_all_na))
}

## 4. PARSE PROVINCE, HERD, AND YEAR =============================================

### 4.1 Split herd_year into components ------------------------------------------

m <- stringr::str_match(dat_raw$herd_year, "^(ON|QC)_(.*)_(\\d{4})$")

dat <- dat_raw %>%
  mutate(
    province_code = m[, 2],
    herd_core     = m[, 3],
    year          = as.integer(m[, 4])
  ) %>%
  filter(!is.na(province_code), !is.na(herd_core), !is.na(year)) %>%
  filter(year %in% years_keep) %>%
  mutate(
    herd_id = paste0(province_code, "_", herd_core),
    herd_label = herd_core %>%
      str_replace_all("_", " ") %>%
      str_to_title(),
    province = factor(
      province_code,
      levels = c("QC", "ON"),
      labels = c("Quebec", "Ontario")
    )
  ) %>%
  arrange(province_code, herd_id, year)

if (anyNA(dat$province_code) | anyNA(dat$herd_core) | anyNA(dat$year)) {
  stop("Could not parse some herd_year values. Expected format: ON_xxx_YYYY or QC_xxx_YYYY.")
}


### 4.2 Check duplicate herd-year rows -------------------------------------------

dups <- dat %>% count(herd_id, year) %>% filter(n > 1)

if (nrow(dups) > 0) {
  print(dups)
  stop("Duplicate rows found for some herd_id-year combinations.")
}

### 4.3 Add herd-code labels ----------------------------------------------------

dat <- dat %>%
  mutate(
    herd_code = recode(herd_label, !!!herd_code_map, .default = herd_label),
    herd_label_key = normalize_key(herd_label)
  )

### 4.4 Optionally remove selected herds -----------------------------------------

if (length(exclude_herd_codes) > 0) {
  dat <- dat %>%
    filter(!herd_code %in% exclude_herd_codes)
}

### 4.5 Print any unmatched herd labels ------------------------------------------

unmatched_herds <- dat %>%
  distinct(herd_label, herd_code) %>%
  filter(herd_label == herd_code) %>%
  arrange(herd_label)

if (nrow(unmatched_herds) > 0) {
  message("The following herd labels were not found in herd_code_map:")
  print(unmatched_herds, n = Inf)
}

## 5. LOAD HERD POSITION TABLE ===================================================

### 5.1 Read the timber-line position table --------------------------------------

pos_tbl_raw <- readr::read_csv(in_position_csv, show_col_types = FALSE)

### 5.2 Standardize join fields --------------------------------------------------

pos_tbl <- pos_tbl_raw %>%
  transmute(
    province_code = province,
    herd_code,
    herd_label_key = normalize_key(herd_label),
    position_class_raw = position_class
  ) %>%
  distinct()

### 5.3 Join position class by province and herd code ----------------------------

dat <- dat %>%
  left_join(
    pos_tbl %>%
      select(province_code, herd_code, position_class_raw) %>%
      distinct(),
    by = c("province_code", "herd_code")
  )

### 5.4 Fallback join by province and normalized herd label ----------------------

missing_pos_idx <- is.na(dat$position_class_raw)

if (any(missing_pos_idx)) {
  
  dat_missing <- dat[missing_pos_idx, ] %>%
    select(province_code, herd_label_key)
  
  fallback_tbl <- pos_tbl %>%
    select(province_code, herd_label_key, position_class_raw) %>%
    distinct()
  
  dat_fallback <- left_join(
    dat_missing,
    fallback_tbl,
    by = c("province_code", "herd_label_key")
  )
  
  dat$position_class_raw[missing_pos_idx] <- dat_fallback$position_class_raw
}

### 5.5 Final check for unmatched herds ------------------------------------------

still_missing <- dat %>%
  filter(is.na(position_class_raw)) %>%
  distinct(province_code, herd_label, herd_code)

if (nrow(still_missing) > 0) {
  print(still_missing, n = Inf)
  stop("Some herds could not be matched to the timber-line position table.")
}

### 5.6 Format final facet variable ----------------------------------------------

dat <- dat %>%
  mutate(
    position_class_raw = factor(position_class_raw, levels = position_levels_raw),
    position_class = recode(
      as.character(position_class_raw),
      "North (>75%)"   = "North (>75%)",
      "South (>75%)"   = "South (>75%)",
      "Line straddler" = "On the line"
    ),
    position_class = factor(position_class, levels = position_levels_plot)
  )

herd_lookup_table <- dat %>%
  distinct(province, herd_label, herd_code, position_class) %>%
  arrange(province, position_class, herd_label)

print(herd_lookup_table, n = Inf)

## 6. BUILD HABITAT-ONLY ORDINATION ==============================================

### 6.1 Prepare habitat matrix ---------------------------------------------------

hab_dat <- dat %>%
  select(all_of(hab_cols))

if (any(is.na(hab_dat))) {
  stop("One or more habitat columns contain NA values.")
}

hab_hel <- vegan::decostand(hab_dat, method = "hellinger")
X_hab <- as.matrix(hab_hel)

### 6.2 Build distance matrix and PCoA -------------------------------------------

d <- dist(X_hab, method = "euclidean")
pcoa <- cmdscale(d, k = 2, eig = TRUE)

eig_pos <- pcoa$eig[pcoa$eig > 0]
var_expl <- round(100 * eig_pos / sum(eig_pos), 1)

scores_df <- as_tibble(pcoa$points, .name_repair = "minimal") %>%
  setNames(c("Axis1", "Axis2")) %>%
  bind_cols(
    dat %>%
      select(
        herd_id, herd_label, herd_code,
        province, province_code,
        position_class,
        year, herd_year,
        all_of(road_var)
      )
  ) %>%
  mutate(road_size = .data[[road_var]])

## 7. FIT HABITAT VECTORS ========================================================

### 7.1 Run envfit ---------------------------------------------------------------

set.seed(1)
fit <- envfit(pcoa$points, X_hab, permutations = 999)

vec <- as.data.frame(scores(fit, display = "vectors"))
vec$var  <- rownames(vec)
vec$r2   <- fit$vectors$r
vec$pval <- fit$vectors$pvals

mul <- vegan::ordiArrowMul(fit, display = "vectors")

vec <- vec %>%
  mutate(
    Axis1 = Dim1 * mul,
    Axis2 = Dim2 * mul
  ) %>%
  arrange(desc(r2))

vec_plot <- vec

if (show_only_strong_vectors) {
  vec_plot <- vec_plot %>% filter(r2 >= r2_threshold)
}

### 7.2 Build adjusted vector labels ---------------------------------------------

vec_lbl <- vec_plot %>%
  mutate(
    display_label = recode(var, !!!envfit_label_map),
    lab_x = Axis1 * 1.08,
    lab_y = Axis2 * 1.08,
    lab_x = case_when(
      var == "harv_0-5"  ~ lab_x - 0.12,
      var == "harv_6-20" ~ lab_x - 0.12,
      var == "young_con" ~ lab_x + 0.01,
      TRUE ~ lab_x
    ),
    lab_y = case_when(
      var == "fire"      ~ lab_y + 0.02,
      var == "young_con" ~ lab_y - 0.03,
      var == "old_con"   ~ lab_y - 0.01,
      TRUE ~ lab_y
    )
  )

## 8. BUILD SHARED AXIS LIMITS ===================================================

### 8.1 Calculate symmetric limits -----------------------------------------------

pad_frac <- 0.001

x_vals <- c(scores_df$Axis1, vec_plot$Axis1, vec_lbl$lab_x)
y_vals <- c(scores_df$Axis2, vec_plot$Axis2, vec_lbl$lab_y)

x_half <- max(abs(x_vals), na.rm = TRUE) * (1 + pad_frac)
y_half <- max(abs(y_vals), na.rm = TRUE) * (1 + pad_frac)

xlim_sym <- c(-x_half, x_half)
ylim_sym <- c(-y_half, y_half)

## 9. BUILD REPEL OBSTACLE POINTS ================================================

### 9.1 Real herd label points ---------------------------------------------------

label_df <- scores_df %>% filter(year == label_year)

label_points <- label_df %>%
  transmute(
    Axis1,
    Axis2,
    province,
    position_class,
    herd_code
  )

### 9.2 Trajectory obstacle points -----------------------------------------------

traj_segments <- scores_df %>%
  arrange(herd_id, year) %>%
  group_by(herd_id, province, position_class) %>%
  mutate(
    x_next = lead(Axis1),
    y_next = lead(Axis2)
  ) %>%
  ungroup() %>%
  filter(!is.na(x_next), !is.na(y_next))

traj_obstacles <- purrr::pmap_dfr(
  list(
    traj_segments$Axis1,
    traj_segments$Axis2,
    traj_segments$x_next,
    traj_segments$y_next,
    traj_segments$province,
    traj_segments$position_class
  ),
  function(x1, y1, x2, y2, province, position_class) {
    tibble(
      Axis1 = seq(x1, x2, length.out = n_path_pts),
      Axis2 = seq(y1, y2, length.out = n_path_pts),
      province = province,
      position_class = position_class,
      herd_code = ""
    )
  }
)

### 9.3 Panel grid for repeated envfit obstacles ---------------------------------

panel_grid <- tidyr::expand_grid(
  province = factor(c("Quebec", "Ontario"), levels = levels(scores_df$province)),
  position_class = factor(position_levels_plot, levels = levels(scores_df$position_class))
)

### 9.4 Envfit arrow obstacle points ---------------------------------------------

arrow_base <- purrr::pmap_dfr(
  list(vec_plot$Axis1, vec_plot$Axis2),
  function(x2, y2) {
    tibble(
      Axis1 = seq(0, x2, length.out = n_arrow_pts),
      Axis2 = seq(0, y2, length.out = n_arrow_pts)
    )
  }
)

arrow_obstacles <- tidyr::crossing(
  panel_grid,
  arrow_base
) %>%
  mutate(herd_code = "")

### 9.5 Envfit label obstacle points ---------------------------------------------

label_base <- vec_lbl %>%
  transmute(
    Axis1 = lab_x,
    Axis2 = lab_y
  )

envfit_label_obstacles <- tidyr::crossing(
  panel_grid,
  label_base
) %>%
  mutate(herd_code = "")

### 9.6 Combine repel points -----------------------------------------------------

repel_df <- bind_rows(
  label_points,
  traj_obstacles,
  arrow_obstacles,
  envfit_label_obstacles
)

## 10. BUILD AND SAVE PLOT =======================================================

### 10.1 Build the PCoA figure --------------------------------------------------

p <- ggplot(scores_df, aes(Axis1, Axis2, colour = factor(year), group = herd_id)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  
  # Habitat vectors
  geom_segment(
    data = vec_plot,
    aes(x = 0, y = 0, xend = Axis1, yend = Axis2),
    inherit.aes = FALSE,
    colour = "grey50",
    arrow = arrow(length = unit(0.15, "inches")),
    linewidth = 0.8
  ) +
  geom_text(
    data = vec_lbl,
    aes(x = lab_x, y = lab_y, label = display_label),
    inherit.aes = FALSE,
    colour = "grey50",
    size = 3
  ) +
  
  # Herd trajectories and points
  geom_path(alpha = 0.45) +
  geom_point(aes(size = road_size), alpha = 0.9) +
  
  # Herd labels
  ggrepel::geom_text_repel(
    data = repel_df,
    aes(x = Axis1, y = Axis2, label = herd_code),
    inherit.aes = FALSE,
    colour = "black",
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.45,
    point.padding = 0.25,
    force = 4,
    force_pull = 0.6,
    min.segment.length = 0.1,
    segment.color = "black",
    segment.size = 0.3,
    segment.alpha = 0.7,
    seed = 123
  ) +
  
  facet_grid(
    rows = vars(province),
    cols = vars(position_class),
    drop = FALSE
  ) +
  coord_equal(xlim = xlim_sym, ylim = ylim_sym) +
  
  scale_colour_manual(
    values = c("1985" = "royalblue", "2000" = "yellowgreen", "2020" = "gold"),
    name = "Year"
  ) +
  scale_size_continuous(
    name = road_legend_title,
    breaks = c(0, 0.25, 0.50, 0.75)
  ) +
  
  labs(
    x = paste0("PCoA 1 (", var_expl[1], "%)"),
    y = paste0("PCoA 2 (", var_expl[2], "%)")
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold")
  )

### 10.2 Print and save ---------------------------------------------------------

print(p)

ggsave(
  filename = out_plot,
  plot = p,
  width = 15,
  height = 10,
  dpi = 300
)

message("Saved: ", out_plot)

