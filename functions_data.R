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

coerce_hpiv3_protein_values <- function(data, protein_cols) {
  token_pattern <- "(?i)^\\s*(<\\s*llod|>\\s*ulod)\\s*$"
  qc_rows <- vector("list", length(protein_cols))
  out <- data

  for (i in seq_along(protein_cols)) {
    col <- protein_cols[[i]]
    raw_chr <- as.character(out[[col]])
    trimmed <- trimws(raw_chr)
    blank_or_na <- is.na(raw_chr) | trimmed == ""
    censored <- grepl(token_pattern, raw_chr, perl = TRUE)

    cleaned_chr <- raw_chr
    cleaned_chr[blank_or_na | censored] <- NA_character_

    numeric_vals <- suppressWarnings(as.numeric(cleaned_chr))
    other_non_numeric <- !is.na(cleaned_chr) & is.na(numeric_vals)

    out[[col]] <- numeric_vals
    qc_rows[[i]] <- tibble::tibble(
      PROTEIN = col,
      n_rows = length(raw_chr),
      n_censored_to_na = sum(censored, na.rm = TRUE),
      n_other_non_numeric_to_na = sum(other_non_numeric, na.rm = TRUE),
      n_missing_after_conversion = sum(is.na(numeric_vals))
    )
  }

  qc <- dplyr::bind_rows(qc_rows)
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

  numeric_conversion <- coerce_hpiv3_protein_values(protein_data, protein_cols)
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
