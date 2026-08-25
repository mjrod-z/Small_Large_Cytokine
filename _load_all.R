# R/_load_all.R

required_pkgs <- c(
  "here", "dplyr", "tidyr", "readr", "stringr",
  "purrr", "tibble", "ggplot2", "ggnewscale", "ggpattern"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(ggnewscale)
  library(ggpattern)
})

project_root <- here::here()

message("Loading project from: ", project_root)

# 1) Config first
config_path <- here::here("R", "config.R")
if (!file.exists(config_path)) {
  stop("Missing config file: ", config_path)
}
source(config_path, local = FALSE)

# 2) Source function files in explicit order
function_files <- c(
  here::here("R", "functions_data.R"),
  here::here("R", "functions_plots.R"),
  here::here("R", "functions_analysis.R")
)

missing_files <- function_files[!file.exists(function_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required R scripts:\n",
    paste0(" - ", missing_files, collapse = "\n")
  )
}

for (f in function_files) {
  message("Sourcing: ", basename(f))
  source(f, local = FALSE)
}

# 3) Verify critical functions/constants exist
required_objects <- c(
  "filter_genes",
  "normalize_counts",
  "save_plot",
  "save_table",
  "GENE_BACKGROUND_THRESHOLD",
  "PATH_DATA_RAW",
  "PATH_DATA_PROCESSED",
  "PATH_OUTPUT_FIGS",
  "PATH_OUTPUT_TABLES",
  "PATH_OUTPUT_REPORT"
)

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), inherits = TRUE)]
if (length(missing_objects) > 0) {
  stop(
    "Bootstrap completed, but required objects are still missing:\n",
    paste0(" - ", missing_objects, collapse = "\n")
  )
}

message("✓ Project configuration and functions loaded")