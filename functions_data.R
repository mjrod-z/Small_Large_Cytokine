# R/functions_data.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_-]", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

ensure_dir <- function(path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  invisible(path)
}

# Add these helpers to functions_data.R (below standardize_exposure is fine)

parse_race_age <- function(df) {
  # Handles a single combined column like "RACE, AGE" containing "Hispanic, 20"
  race_age_col <- names(df)[tolower(gsub("\\s+", "", names(df))) %in% c("race,age", "raceage")]
  if (length(race_age_col) == 0) return(df)
  
  col <- race_age_col[1]
  parts <- stringr::str_split_fixed(as.character(df[[col]]), ",\\s*", 2)
  
  # Only fill missing RACE/AGE if not already present/populated
  if (!"RACE" %in% names(df)) df$RACE <- NA_character_
  if (!"AGE"  %in% names(df)) df$AGE  <- NA_real_
  
  race_new <- trimws(parts[, 1])
  age_new  <- suppressWarnings(as.numeric(trimws(parts[, 2])))
  
  df$RACE <- dplyr::if_else(is.na(df$RACE) | df$RACE == "", race_new, as.character(df$RACE))
  df$AGE  <- dplyr::if_else(is.na(df$AGE), age_new, as.numeric(df$AGE))
  
  df
}

derive_hormone_from_exposure <- function(exposure_chr, default_hormone = "NONE") {
  x <- as.character(exposure_chr)
  has_e2 <- grepl("1\\s*nM\\s*(B-)?Estradiol|\\+\\s*E2", x, ignore.case = TRUE)
  ifelse(has_e2, "Estradiol", default_hormone)
}

standardize_exposure <- function(x) {
  x <- trimws(as.character(x))
  
  # Remove dose units text and "Smoldering" prefix
  x <- gsub("µg/cm2|ug/cm2", "", x, ignore.case = TRUE)
  x <- gsub("^Smoldering\\s+", "", x, ignore.case = TRUE)
  
  # Remove hormone/paraffin suffixes
  x <- gsub(",?\\s*1\\s*nM\\s*B-Estradiol\\s*$", "", x, ignore.case = TRUE)
  x <- gsub(",?\\s*1\\s*nM\\s*Estradiol\\s*$",   "", x, ignore.case = TRUE)
  x <- gsub(",?\\s*\\+\\s*E2\\s*$",              "", x, ignore.case = TRUE)
  x <- gsub(",?\\s*\\+\\s*Parafin\\s*$",         "", x, ignore.case = TRUE)
  x <- trimws(x)
  
  # Collapse spaces
  x <- gsub("\\s+", " ", x)
  
  # Canonical mappings
  x <- gsub("^Red Oak (5|25)$",    "RedOak_\\1",      x, ignore.case = TRUE)
  x <- gsub("^Eucalyptus (5|25)$", "Eucalyptus_\\1",  x, ignore.case = TRUE)
  x <- gsub("^Pine (5|25)$",       "Pine_\\1",        x, ignore.case = TRUE)
  x <- gsub("^Peat (5|25)$",       "Peat_\\1",        x, ignore.case = TRUE)
  
  x <- gsub("^PBS Control$",       "PBS_Control",      x, ignore.case = TRUE)
  x <- gsub("^PBS$",               "PBS_Control",      x, ignore.case = TRUE)
  x <- gsub("^Untreated Control$", "Untreated_Control",x, ignore.case = TRUE)
  x <- gsub("^Untreated$",         "Untreated_Control",x, ignore.case = TRUE)
  x <- gsub("^Unexposed$",         "Untreated_Control",x, ignore.case = TRUE)
  
  # Final cleanup
  x <- gsub(" ", "_", x)
  x <- gsub("^Red_Oak_(5|25)$", "RedOak_\\1", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

clean_sala_metadata <- function(metadata) {
  md <- metadata
  
  # Normalize sample IDs if present
  if ("SAMPLEID" %in% names(md)) {
    md$SAMPLEID <- standardize_sample_id(md$SAMPLEID)
  }
  
  # Split combined race/age field if present
  md <- parse_race_age(md)
  
  # Derive HORMONE from exposure text if HORMONE missing/blank
  if (!"HORMONE" %in% names(md)) md$HORMONE <- NA_character_
  if ("EXPOSURE" %in% names(md)) {
    inferred_hormone <- derive_hormone_from_exposure(md$EXPOSURE)
    md$HORMONE <- dplyr::if_else(
      is.na(md$HORMONE) | trimws(md$HORMONE) == "",
      inferred_hormone,
      as.character(md$HORMONE)
    )
    md$EXPOSURE <- standardize_exposure(md$EXPOSURE)
  }
  
  # Standardize key factors
  if ("SEX" %in% names(md)) {
    md$SEX <- toupper(trimws(as.character(md$SEX)))
  }
  if ("CELLTYPE" %in% names(md)) {
    md$CELLTYPE <- toupper(trimws(as.character(md$CELLTYPE)))
  }
  if ("HORMONE" %in% names(md)) {
    md$HORMONE <- dplyr::case_when(
      grepl("^E2$|Estradiol", md$HORMONE, ignore.case = TRUE) ~ "Estradiol",
      TRUE ~ "NONE"
    )
  }
  
  md
}

# ============================================================================
# SALA-SPECIFIC DATA FUNCTIONS
# ============================================================================

average_nonzero_by_sample <- function(data) {
  # Average replicate measurements within each sample.
  # All numeric values (including zeros) are included in the mean so that the
  # downstream LOD-imputation step operates on unbiased sample averages.
  # NA values are excluded via na.rm = TRUE.
  
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  id_col <- names(data)[!names(data) %in% numeric_cols][1]
  
  if (is.na(id_col)) {
    stop("average_nonzero_by_sample(): no sample ID column found")
  }
  
  data %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(numeric_cols),
                    ~ mean(., na.rm = TRUE)),
      .groups = "drop"
    )
}

build_sala_full <- function(averaged_data, metadata) {
  # Merge averaged MSD data with metadata and apply factor levels
  # Remove PLATE only if it exists
  
  merged <- averaged_data %>%
    dplyr::left_join(metadata, by = "SAMPLEID") %>%
    dplyr::mutate(
      CELLTYPE    = factor(CELLTYPE,    levels = CELLTYPE_LEVELS),
      HORMONE     = factor(HORMONE,     levels = HORMONE_LEVELS),
      SEX         = factor(SEX,         levels = SEX_LEVELS),
      EXPOSURE    = standardize_exposure(EXPOSURE),
      PATIENTCODE = factor(PATIENTCODE)
    )
  
  # Only remove PLATE if it exists
  if ("PLATE" %in% names(merged)) {
    merged <- merged %>% dplyr::select(-PLATE)
  }
  
  merged
}

apply_factor_spec <- function(data, 
                              celltype_levels = CELLTYPE_LEVELS,
                              hormone_levels = HORMONE_LEVELS,
                              sex_levels = SEX_LEVELS) {
  # NOTE: This is an alias for the factor-coercion already done inside
  # build_sala_full(). Prefer calling build_sala_full() directly.
  data %>%
    dplyr::mutate(
      CELLTYPE = factor(CELLTYPE, levels = celltype_levels),
      HORMONE  = factor(HORMONE,  levels = hormone_levels),
      SEX      = factor(SEX,      levels = sex_levels)
    )
}

# ============================================================================
# STANDARDIZATION FUNCTIONS
# ============================================================================

#' Standardize sample IDs (remove hyphens, fix prefixes)
standardize_sample_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\.fastq.*$", "", x, ignore.case = TRUE)
  x <- gsub("\\.bam$", "", x, ignore.case = TRUE)
  x <- gsub("[[:space:]]+", "", x)
  # preserve SALA_1 exactly; do not drop underscores or leading letters
  x
}
#' Standardize exposure names (remove E2 suffix, clean formatting)
standardize_exposure <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  
  # unify separators first
  x <- gsub("_+", " ", x)
  x <- gsub("\\s+", " ", x)
  
  # remove hormone/paraffin suffixes
  x <- gsub(",?\\s*1\\s*nM\\s*(B-)?Estradiol\\s*$", "", x, ignore.case = TRUE)
  x <- gsub("\\+\\s*E2\\s*$", "", x, ignore.case = TRUE)
  x <- gsub("\\+\\s*Parafin\\s*$", "", x, ignore.case = TRUE)
  
  # remove 'Smoldering' prefix and dose units robustly
  x <- gsub("^Smoldering\\s+", "", x, ignore.case = TRUE)
  x <- gsub("\\s*[µu]g\\s*/\\s*cm2\\s*$", "", x, ignore.case = TRUE)
  
  x <- trimws(gsub("\\s+", " ", x))
  
  # canonical mappings
  x <- gsub("^Red Oak (5|25)$", "RedOak_\\1", x, ignore.case = TRUE)
  x <- gsub("^Eucalyptus (5|25)$", "Eucalyptus_\\1", x, ignore.case = TRUE)
  x <- gsub("^Pine (5|25)$", "Pine_\\1", x, ignore.case = TRUE)
  x <- gsub("^Peat (5|25)$", "Peat_\\1", x, ignore.case = TRUE)
  
  x <- gsub("^PBS Control$", "PBS_Control", x, ignore.case = TRUE)
  x <- gsub("^PBS$", "PBS_Control", x, ignore.case = TRUE)
  
  x <- gsub("^Untreated Control$", "Untreated_Control", x, ignore.case = TRUE)
  x <- gsub("^Untreated$", "Untreated_Control", x, ignore.case = TRUE)
  x <- gsub("^Unexposed$", "Untreated_Control", x, ignore.case = TRUE)
  
  x <- gsub("^10X Lysis Buffer$", "LYSIS_Buffer", x, ignore.case = TRUE)
  
  # fallback underscore format
  x <- gsub(" ", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  
  x
}

# ============================================================================
# 1. NORMALIZE COUNTS
# ============================================================================

normalize_counts <- function(count_matrix) {
  count_matrix <- as.matrix(count_matrix)
  
  if (!is.numeric(count_matrix)) {
    stop("normalize_counts(): count_matrix must be numeric.")
  }
  
  size_factors <- colSums(count_matrix) / mean(colSums(count_matrix))
  normalized <- sweep(count_matrix, 2, size_factors, "/")
  
  return(normalized)
}

# ============================================================================
# 2. FILTER GENES
# ============================================================================

filter_genes <- function(data, background_threshold = GENE_BACKGROUND_THRESHOLD) {
  if (!"Geneid" %in% names(data)) {
    stop("filter_genes(): data must contain a 'Geneid' column.")
  }
  
  counts <- data %>%
    dplyr::select(-Geneid) %>%
    as.matrix()
  
  if (!is.numeric(counts)) {
    stop("filter_genes(): all columns except 'Geneid' must be numeric counts.")
  }
  
  median_counts  <- apply(counts, 1, median, na.rm = TRUE)
  overall_median <- median(median_counts, na.rm = TRUE)
  
  # Guard: if overall_median is 0 (sparse data), fall back to a count > 0 filter
  # to avoid the threshold collapsing to "> 0 in ≥20% of samples" silently.
  if (!is.finite(overall_median) || overall_median == 0) {
    warning("filter_genes(): overall median is 0 or non-finite; ",
            "applying a simple presence filter (count > 0 in ≥20% of samples).")
    keep_genes <- apply(counts > 0, 1, function(x) sum(x, na.rm = TRUE) > ncol(counts) * 0.20)
  } else {
    keep_genes <- apply(
      counts > (overall_median * background_threshold),
      1,
      function(x) sum(x, na.rm = TRUE) > ncol(counts) * 0.20
    )
  }
  
  filtered_data <- data[keep_genes, , drop = FALSE]
  return(filtered_data)
}

# ============================================================================
# 3. PREPARE METADATA
# ============================================================================

prepare_metadata_seq <- function(metadata_full, sample_ids) {
  id_cols <- c("SAMPLEID", "Sample_ID", "SampleID")
  id_col <- id_cols[id_cols %in% names(metadata_full)][1]
  
  if (is.na(id_col)) {
    stop("prepare_metadata_seq(): no sample ID column found.")
  }
  
  metadata_full[[id_col]] <- standardize_sample_id(metadata_full[[id_col]])
  sample_ids <- standardize_sample_id(sample_ids)
  
  metadata_seq <- metadata_full %>%
    dplyr::filter(.data[[id_col]] %in% sample_ids) %>%
    dplyr::mutate(
      dplyr::across(where(is.character), as.character)
    )
  
  if ("EXPOSURE" %in% names(metadata_seq)) {
    metadata_seq <- metadata_seq %>%
      dplyr::mutate(EXPOSURE = standardize_exposure(EXPOSURE))
  }
  
  metadata_seq
}

# ============================================================================
# 4. CLASSIFY DEGS
# ============================================================================

classify_degs <- function(data, fc_cutoff = LOG2FC_CUTOFF, p_cutoff = ADJ_P_CUTOFF) {
  data %>%
    dplyr::mutate(DEG = dplyr::case_when(
      adj.P.Val <= p_cutoff & log2FC >=  fc_cutoff ~ "UP",
      adj.P.Val <= p_cutoff & log2FC <= -fc_cutoff ~ "DOWN",
      TRUE ~ "NO"
    ))
}

# ============================================================================
# 5. READ GSEA HTML
# ============================================================================

read_gsea_html <- function(html_file) {
  if (!requireNamespace("rvest", quietly = TRUE)) {
    stop("read_gsea_html(): package 'rvest' is required.")
  }
  
  page <- rvest::read_html(html_file)
  tables <- rvest::html_table(page)
  
  if (length(tables) > 0) {
    results <- tables[[1]]
  } else {
    results <- data.frame()
  }
  
  return(results)
}

# ============================================================================
# 6. FILTER GSEA RESULTS
# ============================================================================

filter_gsea_results <- function(gsea_results, fdr_cutoff = FDR_THRESHOLD) {
  gsea_results %>%
    dplyr::mutate(Significant = ifelse(FDR < fdr_cutoff, "Yes", "No")) %>%
    dplyr::arrange(FDR)
}

# ============================================================================
# 7. PREPARE DEG EXPORT
# ============================================================================

prepare_deg_export <- function(deg_data) {
  deg_data %>%
    dplyr::select(Gene = gene_name, log2FC, adj.P.Val) %>%
    dplyr::mutate(
      FC = 2^log2FC,
      log10_padj = -log10(adj.P.Val)
    ) %>%
    dplyr::filter(adj.P.Val < ADJ_P_CUTOFF) %>%
    dplyr::arrange(dplyr::desc(abs(log2FC)))
}

cat("✓ Data / RNA-seq functions loaded\n")
