#!/usr/bin/env Rscript

candidate_paths <- c(
  "functions_data.R",
  file.path("..", "functions_data.R")
)
data_path <- candidate_paths[file.exists(candidate_paths)][1]
stopifnot(
  "Cannot locate functions_data.R; run this script from the repository root (or tests/)." =
    !is.na(data_path)
)
helper_env <- new.env(parent = globalenv())
sys.source(data_path, envir = helper_env)

count_matrix <- matrix(
  c(
    10, 0, 0, 0,
    0, 0, 10, 0,
    0, 0, 0, 0
  ),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(c("gene_lae", "gene_na_airway", "gene_zero"), paste0("S", 1:4))
)
sample_metadata <- data.frame(
  SAMPLENAME = paste0("S", 1:4),
  AIRWAY = c("LAE", "LAE", NA, NA),
  HORMONE = c("NONE", "NONE", "E2", "E2"),
  stringsAsFactors = FALSE
)
thresholds <- list(
  strata_cols = c("AIRWAY", "HORMONE"),
  cpm_threshold = 1,
  min_sample_fraction = 0.5,
  min_total_count = 10
)

with_qc <- do.call(
  helper_env$filter_rnaseq_count_matrix,
  c(list(count_matrix, sample_metadata), thresholds)
)
stopifnot(
  "Genes passing CPM in a 50% stratum and total-count threshold are included" =
    identical(
      with_qc$gene_status$included[match(rownames(count_matrix), with_qc$gene_status$GENEID)],
      c(TRUE, TRUE, FALSE)
    ),
  "NA-valued airway strata are represented in filtering QC" =
    any(is.na(with_qc$stratum_qc$AIRWAY)),
  "Zero-library samples produce finite QC fractions" =
    all(is.finite(with_qc$stratum_qc$fraction_above_cpm))
)

without_qc <- do.call(
  helper_env$filter_rnaseq_count_matrix,
  c(list(count_matrix, sample_metadata), thresholds, list(include_qc = FALSE))
)
stopifnot(
  "include_qc = FALSE returns only the inclusion vector" =
    identical(names(without_qc), "keep_genes"),
  "QC-free filtering preserves inclusion results" =
    identical(without_qc$keep_genes, with_qc$keep_genes)
)

cat("All RNA-seq filter sanity checks passed.\n")
