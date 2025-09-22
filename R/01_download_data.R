# =============================================================================
# 01_download_data.R
# Author : Romain Loup
# Object : Download data
# Date   : 2025-09-22
# =============================================================================

# ---- Packages ----
# 'pacman' pour load/install all needed packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  here, fs, glue, withr, # chemins & env
  readr, readxl, yaml,   # I/O
  dplyr, tidyr, purrr, stringr, tibble, # wrangling
  ggplot2, ggrepel, scales, patchwork,  # viz
  sf, lwgeom, units,                     # géo
  Matrix, transport,                     # maths
  stats, utils                           # base
)

# Unzip big file for GitHub
zip_files = list.files(path = "data", pattern="\\.zip$", full.names = TRUE)
unzip(zip_files, exdir="data")

# Names of each csv file
files = list.files(path = "data", pattern="\\.csv$", full.names = TRUE)

# List of data frames, indexed by file name
data_list <- lapply(files, read.csv)
names(data_list) <- tools::file_path_sans_ext(basename(files))

names(data_list)
# "detailed_language_2024" "distance_mat_2024"      "IFD_tax_classes_2024"
# "IFD_tax_wealth_2024"    "time_mat_2024"          "vote_info_2024"
# "vote_nb_valid_2024"     "vote_nb_yes_2024"       "vote_theme_names"
# "vote_yes_2024"          "f"

# Shapefiles for mapping
ch_2024 = st_read("data/ch_2024/ch_2024.shp")
lakes = st_read("data/ch_2024/SMV25_lakes.shp")
