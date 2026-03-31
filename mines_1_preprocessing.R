# mines_1_preprocessing.R ----------

# Description:
# This script reads the raw MinCan mine database, filters it to the study area
# (Quebec and Ontario), identifies which mines were active in each snapshot year,
# and exports one CSV file per year for downstream raster processing.
#
# Required inputs:
# - data/mines/raw/MinCan _Past and Present Productive Mines of Canada, 1950-2022_March2024.xlsx
#   OR, if that file is not present locally, the downloadable MinCan source file from figshare
#
# Expected outputs:
# - data/mines/processed/1985_mines.csv
# - data/mines/processed/1990_mines.csv
# - data/mines/processed/1995_mines.csv
# - data/mines/processed/2000_mines.csv
# - data/mines/processed/2005_mines.csv
# - data/mines/processed/2010_mines.csv
# - data/mines/processed/2015_mines.csv
# - data/mines/processed/2020_mines.csv
# - data/mines/processed/merged_mine_data.csv
# - plots/timeline_of_mines_activity2.jpeg

## 1. SETUP ======================================================================

### 1.1 Define project root ------------------------------------------------------

# This function makes the script portable.
# It tries to identify the project root automatically so the packaged folder
# can be placed anywhere on a computer and still run correctly.
#
# Expected packaged structure:
# caribou_hindcast_packaged/
#   ├─ R/
#   ├─ data/
#   └─ plots/

get_project_root <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  
  # Case 1: script is being run with Rscript or sourced from a file path
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    return(dirname(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))))
  }
  
  # Case 2: script is being run interactively in RStudio
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
    editor_path <- tryCatch(
      rstudioapi::getSourceEditorContext()$path,
      error = function(e) ""
    )
    
    if (nzchar(editor_path)) {
      return(dirname(dirname(normalizePath(editor_path, winslash = "/", mustWork = TRUE))))
    }
  }
  
  # Case 3: fallback if script path is available through sys.frames()
  script_path <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  
  if (!is.null(script_path)) {
    return(dirname(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))))
  }
  
  # Case 4: final fallback to working directory, but only if it looks like
  # the project root (i.e., contains both "scripts" and "data" folders)
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

### 1.2 Load required packages ---------------------------------------------------

suppressPackageStartupMessages({
  library(readxl)
  library(httr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

### 1.3 Define helper function --------------------------------------------------

# This helper function takes the long-format mine table and returns only the mines
# that were active in a given year.
#
# A mine is treated as active if:
# - its opening year is less than or equal to the snapshot year, and
# - its closing year is greater than or equal to the snapshot year

mines_active_in <- function(data, year) {
  data %>%
    filter(open <= year, close >= year) %>%
    distinct(ID, .keep_all = TRUE)  # Removes duplicates that can arise when a mine
                                    # has multiple operating periods in the source data.
}

### 1.4 Define file paths -------------------------------------------------------

raw_mines_path <- file.path(
  project_root,
  "data", "mines", "raw",
  "MinCan _Past and Present Productive Mines of Canada, 1950-2022_March2024.xlsx"
)

processed_dir <- file.path(project_root, "data", "mines", "processed")
plots_dir     <- file.path(project_root, "plots")

# Create output folders if they do not already exist.
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

## 2. READ INPUT DATA ============================================================

### 2.1 Import MinCan dataset ---------------------------------------------------

# Checks whether the Excel file exists locally in the packaged
# project folder. If not, it downloads the source file temporarily.

raw_mines_url <- "https://figshare.com/ndownloader/files/45011833"

if (file.exists(raw_mines_path)) {
  mincan_raw <- read_excel(raw_mines_path, sheet = "Data")
} else {
  temp_file <- tempfile(fileext = ".xlsx")
  GET(raw_mines_url, write_disk(temp_file, overwrite = TRUE))
  mincan_raw <- read_excel(temp_file, sheet = "Data")
}

# Quick column check for the imported data structure
colnames(mincan_raw)

## 3. FILTER AND RESHAPE MINE DATA ==============================================

### 3.1 Keep study area and required columns ------------------------------------

# - Keep only Quebec and Ontario mines because those are the study provinces
# - Convert opening/closing year fields to numeric where possible
# - Keep only mine name, coordinates, and operating-year fields
# - Create a unique ID for each record so mines can be tracked cleanly through
#   later filtering and reshaping

mincan_filtered <- mincan_raw %>%
  filter(province %in% c("Quebec", "Ontario")) %>%   # study area filter
  mutate(across(starts_with("close"), ~ suppressWarnings(as.numeric(.x)))) %>%
  mutate(across(starts_with("open"),  ~ suppressWarnings(as.numeric(.x)))) %>%
  select(
    namemine,
    latitude,
    longitude,
    open1,
    open2,
    open3,
    close1,
    close2,
    close3
  ) %>%
  mutate(ID = row_number())

# Optional check to confirm the filtered table looks reasonable
summary(mincan_filtered)
head(mincan_filtered)

### 3.2 Convert opening/closing periods to long format --------------------------

# In the source file, each mine can have up to three open/close periods
# (open1/close1, open2/close2, open3/close3).
#
# This block converts those paired columns into a long format where each row is
# one operating period, i.e., "Was this mine active in year X?"
#
# - "open" in a closing-year field is treated as Inf, meaning still active
# - missing closing years are also treated as Inf
# - rows with no opening year are removed because they do not define a valid period

mincan_long <- mincan_filtered %>%
  mutate(
    across(
      starts_with("close"),
      ~ case_when(
        .x == "open" ~ Inf,
        is.na(.x)    ~ Inf,
        TRUE         ~ as.numeric(.x)
      )
    )
  ) %>%
  pivot_longer(
    cols = matches("^(open|close)[1-3]$"),
    names_to = c(".value", "period"),
    names_pattern = "(open|close)([1-3])"
  ) %>%
  filter(!is.na(open))

## 4. CREATE SNAPSHOT-SPECIFIC MINE TABLES ======================================

### 4.1 Define snapshot years ---------------------------------------------------

snapshot_years <- c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020)

### 4.2 Extract active mines for each snapshot year -----------------------------

# For each year, keep only the mines that were operating during that year.
# The result is stored as a named list, where each element corresponds to one
# snapshot year.

snapshots <- map(
  snapshot_years,
  ~ mines_active_in(mincan_long, .x)
) %>%
  setNames(paste0("mines_", snapshot_years))
list2env(snapshots, .GlobalEnv)

### 4.3 Export one CSV per snapshot year ----------------------------------------

# Each snapshot is written to disk

for (i in seq_along(snapshots)) {
  year_i <- snapshot_years[i]
  filename <- file.path(processed_dir, paste0(year_i, "_mines.csv"))
  
  write.csv(
    snapshots[[i]],
    file = filename,
    row.names = FALSE
  )
}

## 5. OPTIONAL VISUAL CHECKS ====================================================

### 5.1 Plot number of active mines by year -------------------------------------

mine_counts <- sapply(snapshots, nrow)

counts_df <- data.frame(
  TimePeriod = as.character(snapshot_years),
  MineCount = mine_counts
)

ggplot(counts_df, aes(x = TimePeriod, y = MineCount)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  geom_text(aes(label = MineCount), vjust = -0.5) +
  labs(x = "Year", y = "Number of Mines") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### 5.2 Create mine activity timeline plot --------------------------------------

# This timeline plot shows the operating periods of mines through time and marks
# the snapshot years with dashed vertical lines.

mincan_filtered <- mincan_filtered %>%
  mutate(
    namemine = ifelse(
      namemine == "NA",
      paste0("NA_", row_number()[namemine == "NA"]),
      namemine
    )
  )

timeline_plot <- ggplot() +
  geom_segment(
    data = mincan_long,
    aes(x = open, xend = close, y = namemine, yend = namemine),
    color = "grey",
    size = 2
  ) +
  scale_color_identity() +
  labs(
    title = "Timeline of Mines Activity",
    x = "Year",
    y = "Mine Name"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 6, angle = 0, hjust = 1)) +
  geom_vline(
    xintercept = c(1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020),
    linetype = "dashed",
    color = "black"
  ) +
  scale_x_continuous(limits = c(1950, 2024))

# Save the timeline for reference
ggsave(
  filename = file.path(plots_dir, "timeline_of_mines_activity2.jpeg"),
  plot = timeline_plot,
  width = 15,
  height = 40,
  units = "in",
  dpi = 300
)

## 6. CREATE MERGED TABLE FOR MAP PLOTTING ======================================

### 6.1 Read all processed mine CSVs --------------------------------------------

# Create one combined table containing all yearly mine snapshots

folder_path <- processed_dir
file_list <- list.files(
  path = folder_path,
  pattern = "\\.csv$",
  full.names = TRUE
)

df_list <- list()

for (file in file_list) {
  # The first four characters of each file name are the snapshot year
  snapshot_year <- substr(basename(file), 1, 4)
  
  # Read the file and attach its year as a new column
  df <- read.csv(file)
  df <- df %>% mutate(snapshot_year = snapshot_year)
  
  # Add the data frame to the list
  df_list <- append(df_list, list(df))
}

### 6.2 Merge and save combined mine table --------------------------------------

merged_df <- bind_rows(df_list)

write.csv(
  merged_df,
  file = file.path(processed_dir, "merged_mine_data.csv"),
  row.names = FALSE
)

