# =============================================================================
# functions_data.R
# Data preparation helpers for SALA and HPIV3 protein workflows
# =============================================================================

pick_existing_column <- function(data, candidates, label) {
  match <- intersect(candidates, names(data))
  if (length(match) == 0) {
    stop(label, " column not found. Expected one of: ",
         paste(candidates, collapse = ", "))
  }
  match[[1]]
}

assert_required_columns <- function(data, required, data_name) {
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0) {
    stop(
      data_name, " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  invisible(TRUE)
}

priority_factor <- function(x, priority_levels = character(), sort_remaining = FALSE) {
  x_chr <- trimws(as.character(x))
  observed <- unique(x_chr[!is.na(x_chr) & nzchar(x_chr)])
  remaining <- setdiff(observed, priority_levels)
  if (sort_remaining) {
    suppressWarnings({
      remaining_num <- as.numeric(remaining)
    })
    numeric_remaining <- remaining[!is.na(remaining_num)]
    other_remaining <- remaining[is.na(remaining_num)]
    remaining <- c(
      numeric_remaining[order(as.numeric(numeric_remaining))],
      sort(other_remaining)
    )
  }
  factor(x_chr, levels = c(priority_levels[priority_levels %in% observed], remaining))
}

apply_factor_spec <- function(data,
                              celltype_levels = CELLTYPE_LEVELS,
                              hormone_levels = HORMONE_LEVELS,
                              sex_levels = SEX_LEVELS) {
  out <- data

  if ("CELLTYPE" %in% names(out)) {
    out$CELLTYPE <- priority_factor(out$CELLTYPE, celltype_levels)
  }
  if ("HORMONE" %in% names(out)) {
    out$HORMONE <- priority_factor(out$HORMONE, hormone_levels)
  }
  if ("SEX" %in% names(out)) {
    out$SEX <- priority_factor(out$SEX, sex_levels)
  }

  out
}

average_nonzero_by_sample <- function(data) {
  stopifnot(is.data.frame(data))

  metadata_candidates <- c(
    "SAMPLEID", "Sample_ID", "SampleID", "SAMPLENAME", "PATIENTCODE",
    "CELLTYPE", "AIRWAY", "EXPOSURE", "INFECTION", "HORMONE", "SEX",
    "GROUP", "CONCENTRATION", "TIMEPOINT", "SMOKER", "PLATE"
  )
  measure_cols <- setdiff(
    names(data)[vapply(data, is.numeric, logical(1))],
    metadata_candidates
  )

  if (length(measure_cols) == 0) {
    return(data)
  }

  group_cols <- setdiff(names(data), measure_cols)

  data %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(measure_cols), ~ dplyr::na_if(.x, 0))) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(measure_cols),
        ~ if (all(is.na(.x))) NA_real_ else mean(.x, na.rm = TRUE)
      ),
      .groups = "drop"
    )
}

build_sala_full <- function(sala_avg, metadata) {
  stopifnot(is.data.frame(sala_avg), is.data.frame(metadata))

  sample_col_data <- pick_existing_column(sala_avg, c("SAMPLEID", "Sample_ID", "SampleID"),
                                          "SALA sample ID")
  sample_col_meta <- pick_existing_column(metadata, c("SAMPLEID", "Sample_ID", "SampleID"),
                                          "SALA metadata sample ID")

  sala_std <- sala_avg %>% dplyr::rename(SAMPLEID = dplyr::all_of(sample_col_data))
  meta_std <- metadata %>% dplyr::rename(SAMPLEID = dplyr::all_of(sample_col_meta))

  joined <- sala_std %>%
    dplyr::left_join(meta_std, by = "SAMPLEID", suffix = c("", ".meta"))

  merge_fields <- c(
    "PATIENTCODE", "CELLTYPE", "AIRWAY", "EXPOSURE", "HORMONE", "SEX",
    "SMOKER", "TIMEPOINT", "GROUP", "CONCENTRATION", "PLATE"
  )

  for (field in merge_fields) {
    meta_field <- paste0(field, ".meta")
    if (meta_field %in% names(joined)) {
      if (field %in% names(joined)) {
        joined[[field]] <- dplyr::coalesce(joined[[field]], joined[[meta_field]])
      } else {
        joined[[field]] <- joined[[meta_field]]
      }
    }
  }

  if (!"CELLTYPE" %in% names(joined) && "AIRWAY" %in% names(joined)) {
    joined$CELLTYPE <- joined$AIRWAY
  }

  joined <- joined %>%
    dplyr::select(-dplyr::matches("\\.meta$")) %>%
    apply_factor_spec()

  unmatched <- setdiff(sala_std$SAMPLEID, meta_std$SAMPLEID)
  if (length(unmatched) > 0) {
    warning(length(unmatched), " SALA sample IDs did not match metadata.")
  }

  joined
}

extract_hpiv3_sample_components <- function(sample_names) {
  sample_names <- trimws(as.character(sample_names))
  match_mat <- stringr::str_match(sample_names, "^(.*)_([^_]+)$")
  bad <- is.na(match_mat[, 1]) | !nzchar(match_mat[, 2]) | !nzchar(match_mat[, 3])

  if (any(bad)) {
    stop(
      "Unable to extract SAMPLEID/TIMEPOINT from SAMPLENAME values: ",
      paste(utils::head(unique(sample_names[bad]), 10), collapse = ", ")
    )
  }

  tibble::tibble(
    SAMPLENAME = sample_names,
    SAMPLEID = match_mat[, 2],
    TIMEPOINT = match_mat[, 3]
  )
}

normalize_hpiv3_infection <- function(x, source_name = "INFECTION") {
  x_chr <- toupper(trimws(as.character(x)))
  x_chr[x_chr %in% c("", "NA")] <- NA_character_
  x_chr[x_chr %in% c("NO", "NONE")] <- "NONE"

  unexpected <- sort(unique(stats::na.omit(x_chr[!x_chr %in% c("NONE", "HPIV3")])))
  if (length(unexpected) > 0) {
    stop(
      "Unexpected ", source_name, " values: ",
      paste(unexpected, collapse = ", "),
      ". Expected only NONE/NO or HPIV3."
    )
  }

  x_chr
}

normalize_hpiv3_hormone <- function(x) {
  x_chr <- toupper(trimws(as.character(x)))
  x_chr[x_chr %in% c("", "NA")] <- NA_character_

  unexpected <- sort(unique(stats::na.omit(x_chr[!x_chr %in% c("NONE", "E2")])))
  if (length(unexpected) > 0) {
    stop(
      "Unexpected HORMONE values: ",
      paste(unexpected, collapse = ", "),
      ". Expected only NONE or E2."
    )
  }

  factor(x_chr, levels = c("NONE", "E2"))
}

coerce_hpiv3_protein_values <- function(data, protein_cols,
                                        llod_table = NULL, ulod_table = NULL) {
  qc_rows <- vector("list", length(protein_cols))
  out <- data
  llod_map <- if (!is.null(llod_table) &&
                  all(c("Analyte", "LLOD") %in% names(llod_table))) {
    stats::setNames(as.numeric(llod_table$LLOD), as.character(llod_table$Analyte))
  } else {
    numeric()
  }
  ulod_map <- if (!is.null(ulod_table) &&
                  all(c("PROTEIN", "ULOD") %in% names(ulod_table))) {
    stats::setNames(as.numeric(ulod_table$ULOD), as.character(ulod_table$PROTEIN))
  } else {
    numeric()
  }

  for (i in seq_along(protein_cols)) {
    col <- protein_cols[[i]]
    raw_chr <- as.character(out[[col]])
    trimmed <- trimws(raw_chr)
    blank_or_na <- is.na(raw_chr) | trimmed == ""
    left_censored <- grepl("(?i)^\\s*<\\s*llod\\s*$", raw_chr, perl = TRUE)
    right_censored <- grepl("(?i)^\\s*>\\s*ulod\\s*$", raw_chr, perl = TRUE)
    censored <- left_censored | right_censored

    cleaned_chr <- raw_chr
    cleaned_chr[blank_or_na | censored] <- NA_character_

    numeric_vals <- suppressWarnings(as.numeric(cleaned_chr))
    other_non_numeric <- !is.na(cleaned_chr) & is.na(numeric_vals)
    quantified <- numeric_vals[is.finite(numeric_vals) & numeric_vals > 0]
    left_limit <- if (col %in% names(llod_map)) unname(llod_map[col]) else NA_real_
    if (length(left_limit) == 0L || !is.finite(left_limit) || left_limit <= 0) {
      left_limit <- if (length(quantified) > 0L) min(quantified) else NA_real_
      left_method <- "minimum quantified value proxy"
    } else {
      left_method <- "analyte LLOD reference"
    }
    right_limit <- if (col %in% names(ulod_map)) unname(ulod_map[col]) else NA_real_
    if (length(right_limit) == 0L || !is.finite(right_limit) || right_limit <= 0) {
      right_limit <- if (length(quantified) > 0L) max(quantified) else NA_real_
      right_method <- "maximum quantified value proxy"
    } else {
      right_method <- "analyte ULOD reference"
    }
    if (is.finite(left_limit)) numeric_vals[left_censored] <- left_limit / sqrt(2)
    if (is.finite(right_limit)) numeric_vals[right_censored] <- right_limit * sqrt(2)

    out[[col]] <- numeric_vals
    qc_rows[[i]] <- tibble::tibble(
      PROTEIN = col,
      n_rows = length(raw_chr),
      n_left_censored = sum(left_censored, na.rm = TRUE),
      n_right_censored = sum(right_censored, na.rm = TRUE),
      n_left_imputed = sum(left_censored & is.finite(left_limit), na.rm = TRUE),
      n_right_imputed = sum(right_censored & is.finite(right_limit), na.rm = TRUE),
      n_censored_unimputed = sum(censored & is.na(numeric_vals), na.rm = TRUE),
      censoring_method = paste(
        if (any(left_censored)) left_method else NA_character_,
        if (any(right_censored)) right_method else NA_character_,
        sep = "; "
      ),
      n_censored_to_na = sum(censored & is.na(numeric_vals), na.rm = TRUE),
      n_other_non_numeric_to_na = sum(other_non_numeric, na.rm = TRUE),
      n_missing_after_conversion = sum(is.na(numeric_vals))
    )
  }

  qc <- dplyr::bind_rows(qc_rows)
  qc$censoring_method <- gsub("^NA; |; NA$", "", qc$censoring_method)
  bad_numeric <- qc %>% dplyr::filter(n_other_non_numeric_to_na > 0)
  if (nrow(bad_numeric) > 0) {
    warning(
      "Non-numeric protein values were coerced to NA in: ",
      paste0(
        bad_numeric$PROTEIN,
        " (", bad_numeric$n_other_non_numeric_to_na, ")",
        collapse = ", "
      )
    )
  }

  list(data = out, qc = qc)
}

build_hpiv3_analysis_data <- function(
    protein_path = here::here(PATH_DATA_RAW, "HPIV3_nomic_subset.csv"),
    metadata_path = here::here(PATH_DATA_RAW, "HPIV3_metadata.csv")) {

  if (!file.exists(protein_path)) {
    stop("Missing HPIV3 protein data file: ", protein_path)
  }
  if (!file.exists(metadata_path)) {
    stop(
      "Missing HPIV3 metadata file: ", metadata_path,
      ". Expected a CSV named 'HPIV3_metadata.csv' in ",
      dirname(metadata_path), "."
    )
  }

  build_hpiv3_rnaseq_data <- function(
      count_path = here::here(PATH_DATA_RAW, "HPIV3_RNA_data.csv"),
      metadata_path = here::here(PATH_DATA_RAW, "HPIV3_metadata.csv")) {
    if (!file.exists(count_path)) {
      stop("Missing HPIV3 RNA-seq count file: ", count_path)
    }
    if (!file.exists(metadata_path)) {
      stop("Missing HPIV3 metadata file: ", metadata_path)
    }

    counts_raw <- readr::read_csv(
      count_path,
      show_col_types = FALSE,
      col_types = readr::cols(.default = readr::col_character())
    )
    metadata_raw <- readr::read_csv(
      metadata_path,
      show_col_types = FALSE,
      col_types = readr::cols(.default = readr::col_character())
    )
    if (ncol(counts_raw) < 2 || names(counts_raw)[[1]] != "GENEID") {
      stop("HPIV3 RNA-seq counts must have GENEID as the first column and at least one sample column.")
    }
    assert_required_columns(
      metadata_raw,
      c("SAMPLEID", "AIRWAY", "PATIENTCODE", "EXPOSURE", "HORMONE", "INFECTION", "SEX"),
      "HPIV3 metadata"
    )

    sample_names <- names(counts_raw)[-1]
    if (anyNA(sample_names) || any(!nzchar(trimws(sample_names))) || anyDuplicated(sample_names)) {
      stop("HPIV3 RNA-seq sample column names must be non-empty and unique.")
    }
    if (anyDuplicated(counts_raw$GENEID) || any(is.na(counts_raw$GENEID)) ||
        any(!nzchar(trimws(counts_raw$GENEID)))) {
      stop("HPIV3 RNA-seq GENEID values must be non-empty and unique.")
    }

    sample_parts <- extract_hpiv3_sample_components(sample_names)
    counts <- counts_raw %>%
      dplyr::mutate(GENEID = trimws(as.character(GENEID))) %>%
      dplyr::mutate(dplyr::across(
        dplyr::all_of(sample_names),
        ~ suppressWarnings(as.numeric(.x))
      ))
    count_values <- unlist(counts[sample_names], use.names = FALSE)
    if (anyNA(count_values) || any(!is.finite(count_values)) ||
        any(count_values < 0) || any(count_values != floor(count_values))) {
      stop("HPIV3 RNA-seq sample columns must contain finite, non-negative integer raw counts.")
    }

    metadata_std <- metadata_raw %>%
      dplyr::mutate(
        SAMPLEID = trimws(as.character(SAMPLEID)),
        AIRWAY = trimws(as.character(AIRWAY)),
        PATIENTCODE = trimws(as.character(PATIENTCODE)),
        EXPOSURE = trimws(as.character(EXPOSURE)),
        INFECTION = normalize_hpiv3_infection(INFECTION, "metadata INFECTION"),
        HORMONE = normalize_hpiv3_hormone(HORMONE),
        SEX = priority_factor(toupper(trimws(as.character(SEX))), c("F", "M"))
      )
    if (anyDuplicated(metadata_std$SAMPLEID)) {
      stop("HPIV3 metadata must contain at most one row per SAMPLEID.")
    }

    unmatched_count_ids <- tibble::tibble(
      SAMPLENAME = sort(setdiff(sample_parts$SAMPLENAME, sample_names))
    )
    unmatched_metadata_ids <- tibble::tibble(
      SAMPLEID = sort(setdiff(unique(metadata_std$SAMPLEID), unique(sample_parts$SAMPLEID)))
    )
    unmatched_sample_ids <- tibble::tibble(
      SAMPLEID = sort(setdiff(unique(sample_parts$SAMPLEID), unique(metadata_std$SAMPLEID)))
    )
    if (nrow(unmatched_sample_ids) > 0 || nrow(unmatched_metadata_ids) > 0) {
      warning(
        "HPIV3 RNA-seq SAMPLEID mismatches detected: ",
        nrow(unmatched_sample_ids), " count-only, ",
        nrow(unmatched_metadata_ids), " metadata-only."
      )
    }

    sample_data <- sample_parts %>%
      dplyr::left_join(metadata_std, by = "SAMPLEID") %>%
      dplyr::mutate(
        TIMEPOINT = priority_factor(TIMEPOINT, c("24", "72"), sort_remaining = TRUE),
        AIRWAY = priority_factor(AIRWAY),
        SEX = priority_factor(SEX, c("F", "M")),
        CELLTYPE = AIRWAY
      ) %>%
      apply_factor_spec()

    count_matrix <- as.matrix(counts[, sample_names, drop = FALSE])
    storage.mode(count_matrix) <- "numeric"
    rownames(count_matrix) <- counts$GENEID
    colnames(count_matrix) <- sample_names

    long_counts <- counts %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(sample_names),
        names_to = "SAMPLENAME",
        values_to = "COUNT"
      ) %>%
      dplyr::left_join(sample_data, by = "SAMPLENAME") %>%
      dplyr::select(
        SAMPLENAME, SAMPLEID, TIMEPOINT, AIRWAY, CELLTYPE, PATIENTCODE,
        EXPOSURE, HORMONE, INFECTION, SEX, GENEID, COUNT,
        dplyr::any_of(c("AGE", "RACE", "SMOKER"))
      )

    list(
      data = long_counts,
      sample_data = sample_data,
      count_matrix = count_matrix,
      gene_cols = rownames(count_matrix),
      qc = list(
        unmatched_count_ids = unmatched_count_ids,
        unmatched_sample_ids = unmatched_sample_ids,
        unmatched_metadata_ids = unmatched_metadata_ids
      )
    )
  }

  filter_hpiv3_rnaseq_genes <- function(data, min_total_count = GENE_BACKGROUND_THRESHOLD) {
    stopifnot(is.data.frame(data), all(c("GENEID", "COUNT") %in% names(data)))
    if (length(min_total_count) != 1 || is.na(min_total_count) || min_total_count < 0) {
      stop("min_total_count must be a single non-negative number.")
    }
    qc <- data %>%
      dplyr::group_by(GENEID) %>%
      dplyr::summarise(
        total_count = sum(COUNT, na.rm = TRUE),
        n_samples = dplyr::n_distinct(SAMPLENAME),
        included_background = total_count >= min_total_count,
        .groups = "drop"
      )
    list(
      data = data %>% dplyr::semi_join(
        qc %>% dplyr::filter(included_background) %>% dplyr::select(GENEID),
        by = "GENEID"
      ),
      qc = qc
    )
  }

  compute_hpiv3_rnaseq_filter_qc <- function(
      data,
      cpm_threshold = 1,
      min_sample_fraction = 0.5,
      strata_cols = c("AIRWAY", "HORMONE", "TIMEPOINT")) {
    stopifnot(is.data.frame(data), all(c("GENEID", "SAMPLENAME", "COUNT") %in% names(data)))
    if (length(cpm_threshold) != 1 || is.na(cpm_threshold) || cpm_threshold < 0 ||
        length(min_sample_fraction) != 1 || is.na(min_sample_fraction) ||
        min_sample_fraction < 0 || min_sample_fraction > 1) {
      stop("CPM threshold must be non-negative and sample fraction must be between 0 and 1.")
    }
    missing_strata <- setdiff(strata_cols, names(data))
    if (length(missing_strata) > 0) {
      stop("RNA-seq filter data is missing stratum columns: ", paste(missing_strata, collapse = ", "))
    }

    library_sizes <- data %>%
      dplyr::distinct(SAMPLENAME, .data[[strata_cols[[1]]]], .data[[if (length(strata_cols) > 1) strata_cols[[2]] else strata_cols[[1]]]])
    sample_libraries <- data %>%
      dplyr::group_by(SAMPLENAME) %>%
      dplyr::summarise(library_size = sum(COUNT, na.rm = TRUE), .groups = "drop")
    cpm_data <- data %>%
      dplyr::left_join(sample_libraries, by = "SAMPLENAME") %>%
      dplyr::mutate(CPM = dplyr::if_else(library_size > 0, COUNT / library_size * 1e6, 0))

    stratum_qc <- cpm_data %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c("GENEID", strata_cols)))) %>%
      dplyr::summarise(
        n_samples = dplyr::n_distinct(SAMPLENAME),
        n_samples_above_threshold = dplyr::n_distinct(SAMPLENAME[CPM >= cpm_threshold]),
        fraction_above_threshold = n_samples_above_threshold / n_samples,
        passes_stratum = fraction_above_threshold >= min_sample_fraction,
        .groups = "drop"
      )
    gene_status <- stratum_qc %>%
      dplyr::group_by(GENEID) %>%
      dplyr::summarise(
        max_fraction_above_threshold = max(fraction_above_threshold, na.rm = TRUE),
        included_completeness = any(passes_stratum),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        cpm_threshold = cpm_threshold,
        min_sample_fraction = min_sample_fraction
      )

    list(qc = stratum_qc, gene_status = gene_status)
  }

  normalize_hpiv3_rnaseq <- function(count_matrix, method = c("TMM", "median_ratio")) {
    method <- match.arg(method)
    count_matrix <- as.matrix(count_matrix)
    if (is.null(rownames(count_matrix)) || is.null(colnames(count_matrix)) ||
        anyNA(count_matrix) || any(count_matrix < 0) ||
        any(count_matrix != floor(count_matrix))) {
      stop("count_matrix must be a non-negative integer matrix with gene and sample names.")
    }
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop("The DESeq2 package is required for HPIV3 RNA-seq VST normalization.")
    }

    if (method == "TMM") {
      if (!requireNamespace("edgeR", quietly = TRUE)) {
        stop("The edgeR package is required for TMM normalization.")
      }
      dge <- edgeR::DGEList(counts = count_matrix)
      dge <- edgeR::calcNormFactors(dge, method = "TMM")
      size_factors <- dge$samples$lib.size * dge$samples$norm.factors
      size_factors <- size_factors / exp(mean(log(size_factors)))
    } else {
      size_factors <- DESeq2::estimateSizeFactorsForMatrix(count_matrix)
    }

    names(size_factors) <- colnames(count_matrix)
    normalized_counts <- sweep(count_matrix, 2, size_factors, "/")
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(count_matrix),
      colData = S4Vectors::DataFrame(row.names = colnames(count_matrix)),
      design = ~ 1
    )
    DESeq2::sizeFactors(dds) <- size_factors
    vst_values <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = TRUE))

    list(
      method = method,
      size_factors = tibble::tibble(SAMPLENAME = names(size_factors), size_factor = as.numeric(size_factors)),
      normalized_counts = normalized_counts,
      vst = vst_values
    )
  }

  protein_raw <- readr::read_csv(
    protein_path,
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )
  metadata_raw <- readr::read_csv(
    metadata_path,
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )

  assert_required_columns(
    protein_raw,
    c("SAMPLENAME", "PLATE", "INFECTION", "AIRWAY"),
    "HPIV3 protein data"
  )
  assert_required_columns(
    metadata_raw,
    c("SAMPLEID", "AIRWAY", "PATIENTCODE", "EXPOSURE", "HORMONE",
      "INFECTION", "AGE", "SEX", "RACE", "SMOKER"),
    "HPIV3 metadata"
  )

  sample_parts <- extract_hpiv3_sample_components(protein_raw$SAMPLENAME)
  protein_data <- dplyr::bind_cols(protein_raw, sample_parts %>% dplyr::select(-SAMPLENAME))

  non_protein_cols <- c("SAMPLENAME", "PLATE", "INFECTION", "AIRWAY", "SAMPLEID", "TIMEPOINT")
  protein_cols <- setdiff(names(protein_data), non_protein_cols)
  if (length(protein_cols) == 0) {
    stop("HPIV3 protein data did not contain any protein concentration columns.")
  }

  numeric_conversion <- coerce_hpiv3_protein_values(
    protein_data,
    protein_cols,
    llod_table = if (exists("cytokine_llod", inherits = TRUE)) cytokine_llod else NULL
  )
  protein_data <- numeric_conversion$data

  protein_data <- protein_data %>%
    dplyr::mutate(
      INFECTION = normalize_hpiv3_infection(INFECTION, "protein-data INFECTION"),
      AIRWAY = trimws(as.character(AIRWAY))
    )

  metadata_std <- metadata_raw %>%
    dplyr::mutate(
      SAMPLEID = trimws(as.character(SAMPLEID)),
      AIRWAY = trimws(as.character(AIRWAY)),
      INFECTION = normalize_hpiv3_infection(INFECTION, "metadata INFECTION"),
      HORMONE = normalize_hpiv3_hormone(HORMONE),
      SEX = priority_factor(toupper(trimws(as.character(SEX))), c("F", "M")),
      EXPOSURE = trimws(as.character(EXPOSURE)),
      PATIENTCODE = trimws(as.character(PATIENTCODE)),
      SMOKER = trimws(as.character(SMOKER)),
      AGE = trimws(as.character(AGE)),
      RACE = trimws(as.character(RACE))
    )

  joined <- protein_data %>%
    dplyr::left_join(metadata_std, by = "SAMPLEID", suffix = c("", ".meta"))

  mismatch_airway <- joined %>%
    dplyr::filter(!is.na(AIRWAY), !is.na(AIRWAY.meta), AIRWAY != AIRWAY.meta) %>%
    dplyr::transmute(SAMPLEID, field = "AIRWAY", protein_value = AIRWAY, metadata_value = AIRWAY.meta)

  mismatch_infection <- joined %>%
    dplyr::filter(!is.na(INFECTION), !is.na(INFECTION.meta), INFECTION != INFECTION.meta) %>%
    dplyr::transmute(SAMPLEID, field = "INFECTION", protein_value = INFECTION, metadata_value = INFECTION.meta)

  mismatch_table <- dplyr::bind_rows(mismatch_airway, mismatch_infection)
  if (nrow(mismatch_table) > 0) {
    warning(nrow(mismatch_table), " HPIV3 metadata mismatches detected for AIRWAY/INFECTION.")
  }

  unmatched_protein_ids <- tibble::tibble(
    SAMPLEID = sort(setdiff(unique(protein_data$SAMPLEID), unique(metadata_std$SAMPLEID)))
  )
  unmatched_metadata_ids <- tibble::tibble(
    SAMPLEID = sort(setdiff(unique(metadata_std$SAMPLEID), unique(protein_data$SAMPLEID)))
  )

  if (nrow(unmatched_protein_ids) > 0 || nrow(unmatched_metadata_ids) > 0) {
    warning(
      "HPIV3 SAMPLEID mismatches detected: ",
      nrow(unmatched_protein_ids), " protein-only, ",
      nrow(unmatched_metadata_ids), " metadata-only."
    )
  }

  joined <- joined %>%
    dplyr::mutate(
      INFECTION = dplyr::coalesce(INFECTION, INFECTION.meta),
      AIRWAY = dplyr::coalesce(AIRWAY, AIRWAY.meta),
      HORMONE = HORMONE,
      TIMEPOINT = priority_factor(TIMEPOINT, c("24", "72"), sort_remaining = TRUE),
      AIRWAY = priority_factor(AIRWAY),
      SEX = priority_factor(SEX, c("F", "M")),
      CELLTYPE = AIRWAY
    ) %>%
    dplyr::select(
      SAMPLENAME, SAMPLEID, TIMEPOINT, PLATE, AIRWAY, CELLTYPE, INFECTION,
      PATIENTCODE, EXPOSURE, HORMONE, AGE, SEX, RACE, SMOKER,
      dplyr::all_of(protein_cols)
    )

  list(
    data = joined,
    protein_cols = protein_cols,
    qc = list(
      numeric_conversion = numeric_conversion$qc,
      unmatched_protein_ids = unmatched_protein_ids,
      unmatched_metadata_ids = unmatched_metadata_ids,
      metadata_mismatches = mismatch_table
    )
  )
}

# Subset HPIV3 data by hormone mode: "both", "E2_only", or "exclude_E2".
filter_hpiv3_hormone_mode <- function(data, mode = c("both", "E2_only", "exclude_E2")) {
  mode <- match.arg(mode)
  out <- switch(
    mode,
    both = data,
    E2_only = dplyr::filter(data, as.character(HORMONE) == "E2"),
    exclude_E2 = dplyr::filter(data, as.character(HORMONE) != "E2")
  )
  if (nrow(out) == 0L) stop("No HPIV3 rows remain after hormone subsetting (mode = ", mode, ")")
  out
}
