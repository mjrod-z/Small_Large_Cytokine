#!/usr/bin/env Rscript
# Dependency-free structural checks for the HPIV3 ranked-effect workflow.

candidate_roots <- c(".", "..")
root <- candidate_roots[vapply(
  candidate_roots,
  function(path) file.exists(file.path(path, "functions_analysis.R")),
  logical(1)
)][1]
stopifnot(!is.na(root))

analysis_path <- file.path(root, "functions_analysis.R")
data_path <- file.path(root, "functions_data.R")
plots_path <- file.path(root, "functions_plots.R")
report_path <- file.path(root, "02_hpiv3_analysis.Rmd")
rnaseq_report_path <- file.path(root, "03_rnaseq_hpiv3_analysis.Rmd")
loader_path <- file.path(root, "_load_all.R")
analysis_text <- paste(readLines(analysis_path, warn = FALSE), collapse = "\n")
plots_text <- paste(readLines(plots_path, warn = FALSE), collapse = "\n")
report_text <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
rnaseq_report_text <- paste(readLines(rnaseq_report_path, warn = FALSE), collapse = "\n")
loader_text <- paste(readLines(loader_path, warn = FALSE), collapse = "\n")

data_text <- paste(readLines(data_path, warn = FALSE), collapse = "\n")

for (path in c(analysis_path, data_path, plots_path, loader_path)) parse(file = path)
rmd_lines <- readLines(report_path, warn = FALSE)
ranked_chunk_starts <- which(startsWith(rmd_lines, "```{r hpiv3-ranked-"))
stopifnot(length(ranked_chunk_starts) == 4L)
for (ranked_chunk_start in ranked_chunk_starts) {
  ranked_chunk_end <- which(
    seq_along(rmd_lines) > ranked_chunk_start & grepl("^```$", rmd_lines)
  )[1]
  stopifnot(!is.na(ranked_chunk_end))
  parse(text = rmd_lines[(ranked_chunk_start + 1L):(ranked_chunk_end - 1L)])
}
# An explicit "exposure" ranked-bar/volcano chunk must exist alongside the
# hormone and infection chunks, so exposure PNGs are actually generated.
stopifnot(
  "A dedicated hpiv3-ranked-exposure-bars chunk must exist" =
    any(startsWith(rmd_lines, "```{r hpiv3-ranked-exposure-bars"))
)
euler_chunk_start <- match("```{r hpiv3-exposure-euler}", rmd_lines)
euler_chunk_end <- if (is.na(euler_chunk_start)) {
  NA_integer_
} else {
  which(seq_along(rmd_lines) > euler_chunk_start & rmd_lines == "```")[1]
}
stopifnot(!is.na(euler_chunk_start), !is.na(euler_chunk_end))
parse(text = rmd_lines[(euler_chunk_start + 1L):(euler_chunk_end - 1L)])
model_chunk_start <- match("```{r hpiv3-noninfection-models}", rmd_lines)
model_chunk_end <- which(
  seq_along(rmd_lines) > model_chunk_start & rmd_lines == "```"
)[1]
stopifnot(!is.na(model_chunk_end))
parse(text = rmd_lines[(model_chunk_start + 1L):(model_chunk_end - 1L)])

exposure_function_start <- regexpr(
  "fit_hpiv3_exposure_models <- function",
  analysis_text,
  fixed = TRUE
)
exposure_function_end <- regexpr(
  "summarize_hpiv3_strata <- function",
  analysis_text,
  fixed = TRUE
)
exposure_function_text <- substr(
  analysis_text,
  exposure_function_start,
  exposure_function_end - 1L
)

stopifnot(
  "Hormone models must test E2 against NONE" =
    grepl('control_level = "NONE"', analysis_text) &&
      grepl('case_level = "E2"', analysis_text),
  "Hormone models must stratify by airway/timepoint/exposure/infection" =
    grepl("distinct\\(AIRWAY, TIMEPOINT, EXPOSURE, INFECTION\\)", analysis_text),
  "Hormone models must filter and report the infection stratum" =
    grepl("INFECTION == infection_i", analysis_text, fixed = TRUE) &&
      grepl('INFECTION = as.character\\(infection_i\\)', analysis_text) &&
      grepl("group_by\\(AIRWAY, TIMEPOINT, EXPOSURE, SEX, INFECTION\\)", analysis_text),
  "Hormone models must emit pooled and sex-specific groups" =
    grepl('sex_groups <- c\\("All", intersect\\(c\\("F", "M"\\)', analysis_text),
  "Hormone models must retain mixed and fallback linear model paths" =
    grepl('dat, log2_value ~ HORMONE, pairing_factor = "HORMONE"', analysis_text, fixed = TRUE) &&
      grepl("fit_hpiv3_paired_model <- function", analysis_text, fixed = TRUE),
  "Exposure models must create pooled All-sex input rows" =
    grepl('mutate\\(SEX = "All"\\)', analysis_text),
  "Exposure contrasts must compare treatments with PBS only" =
    grepl('method = "trt.vs.ctrl"', exposure_function_text) &&
      grepl('match\\("PBS_Control", exposure_levels\\)', exposure_function_text) &&
      !grepl('method = "pairwise"', exposure_function_text),
  "Exposure models must stratify by airway/hormone/timepoint/sex/infection so exposure rows survive filtering" =
    grepl(
      "distinct\\(AIRWAY, HORMONE, TIMEPOINT, SEX, INFECTION\\)",
      exposure_function_text
    ),
  "All HPIV3 models must share paired-first model selection and avoid unpaired fallback when pairing is available" =
    grepl("fit_hpiv3_paired_model <- function", analysis_text, fixed = TRUE) &&
      grepl("lme4::lmer(random_formula", analysis_text, fixed = TRUE) &&
      grepl("stats::lm(paired_formula", analysis_text, fixed = TRUE) &&
      grepl("paired model failed; unpaired fallback not used", analysis_text, fixed = TRUE) &&
      grepl("paired_with_patientcode", exposure_function_text, fixed = TRUE) &&
      grepl('"lm_donor_fixed"', analysis_text, fixed = TRUE) &&
      grepl('"lm_unpaired"', analysis_text, fixed = TRUE) &&
      grepl("n_paired_donors", exposure_function_text, fixed = TRUE),
  "Exposure models must re-check per-level minimum counts after same-infection PBS donor filtering" =
    grepl("paired_exposure_counts <- table(dat_fit$EXPOSURE)", exposure_function_text, fixed = TRUE) &&
      grepl("after PBS-donor filtering", exposure_function_text, fixed = TRUE),
  "Mixed-model contrasts must use explicit Kenward-Roger df and report contrast df" =
    grepl('HPIV3_LMER_DF_METHOD <- "kenward-roger"', analysis_text, fixed = TRUE) &&
      grepl("arguments$lmer.df <- HPIV3_LMER_DF_METHOD", analysis_text, fixed = TRUE) &&
      grepl("df_method = NA_character_", analysis_text, fixed = TRUE) &&
      grepl("stats_row$df", analysis_text, fixed = TRUE) &&
      grepl('"pbkrtest"', loader_text, fixed = TRUE),
  "A formal exposure-by-infection model and result export must support differential-response claims" =
    grepl("fit_hpiv3_exposure_infection_interactions <- function", analysis_text, fixed = TRUE) &&
      grepl('interaction = c("trt.vs.ctrl", "revpairwise")', analysis_text, fixed = TRUE) &&
      grepl("hpiv3_exposure_infection_interaction_results.csv", report_text, fixed = TRUE) &&
      grepl("interaction_q.value", report_text, fixed = TRUE),
  "Censored HPIV3 measurements must be substituted and auditable rather than all converted to missing" =
    grepl("n_left_imputed", data_text, fixed = TRUE) &&
      grepl("left_limit / sqrt(2)", data_text, fixed = TRUE) &&
      grepl("right_limit * sqrt(2)", data_text, fixed = TRUE) &&
      grepl("hpiv3_censoring_qc.csv", report_text, fixed = TRUE) &&
      grepl("unique membership does not establish a difference", report_text, fixed = TRUE),
  "Exposure baseline audit must report donor pairing and reject unpaired fits when pairing is possible" =
    grepl("n_donors_with_pbs", analysis_text, fixed = TRUE) &&
      grepl("n_exposure_obs_paired", analysis_text, fixed = TRUE) &&
      grepl("donor pairing was possible but an unpaired model was used", analysis_text, fixed = TRUE),
  "Ranked bars must use estimated effects and direction/significance fills" =
    grepl("geom_col\\(width = 0.75\\)", plots_text) &&
      grepl('"Higher, significant"', plots_text) &&
      grepl('"Lower, significant"', plots_text),
  "Ranked plot helper must be present in bootstrap checks" =
    grepl('"plot_hpiv3_ranked_bars"', loader_text),
  "Volcano plot helper must be present in bootstrap checks" =
    grepl('"plot_hpiv3_volcano"', loader_text),
  "Volcano plots must color by the same significance/direction scheme" =
    grepl("plot_hpiv3_volcano <- function", plots_text) &&
      grepl('"Higher, significant" = UP_COLOR_DEFAULT', plots_text, fixed = TRUE) &&
      grepl('"Lower, significant" = DOWN_COLOR_DEFAULT', plots_text, fixed = TRUE) &&
      grepl("ggrepel::geom_text_repel", plots_text),
  "HPIV3 volcano plots must draw a per-panel significance boundary (no log2FC cutoff), UP/DOWN corner counts, and bordered theme" =
    grepl("hpiv3_volcano_thresholds(plot_df, panel_group)", plots_text, fixed = TRUE) &&
      grepl("ggplot2::aes(yintercept = threshold_y)", plots_text, fixed = TRUE) &&
      !grepl("xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF),\n      color = \"grey50\", linetype = \"dashed\", linewidth = 0.5\n    ) +\n    ggplot2::geom_hline(\n      yintercept = -log10(ADJ_P_CUTOFF),\n      color = \"grey50\", linetype = \"dashed\", linewidth = 0.5\n    ) +\n    ggplot2::geom_point(alpha = 0.7, size = 2)", plots_text, fixed = TRUE) &&
      grepl('paste0("UP: ", n_up)', plots_text, fixed = TRUE) &&
      grepl('paste0("DOWN: ", n_down)', plots_text, fixed = TRUE) &&
      grepl("panel.border = ggplot2::element_rect(", plots_text, fixed = TRUE),
  "Plot helpers must define shared title-wrapping and short dimension-tag utilities" =
    grepl("wrap_plot_text <- function", plots_text, fixed = TRUE) &&
      grepl("hpiv3_title_theme <- function", plots_text, fixed = TRUE) &&
      grepl("hpiv3_dimension_tag <- function", plots_text, fixed = TRUE) &&
      grepl("HPIV3_DIMENSION_TAGS", plots_text, fixed = TRUE),
  "Ranked bar and volcano plots must wrap/resize their titles" =
    grepl("plot_title <- wrap_plot_text", plots_text, fixed = TRUE) &&
      grepl("hpiv3_title_theme()", plots_text, fixed = TRUE),
  "Report must create all three requested ranked comparison types" =
    all(vapply(
      c('"exposure"', '"hormone"', '"infection"'),
      grepl,
      logical(1),
      x = report_text,
      fixed = TRUE
    )),
  "Report must generate ranked exposure, hormone, and infection models" =
    all(vapply(
      c("fit_hpiv3_exposure_models", "fit_hpiv3_hormone_models",
        "fit_hpiv3_infection_models", "hpiv3_ranked_bar_plot_index"),
      grepl,
      logical(1),
      x = report_text,
      fixed = TRUE
    )),
  "Report must save a companion volcano plot alongside each ranked bar chart" =
    grepl("plot_hpiv3_volcano(", report_text, fixed = TRUE) &&
      grepl('"volcano"', report_text, fixed = TRUE) &&
      grepl("hpiv3_volc_", report_text, fixed = TRUE),
  "Exposure ranked-bar/volcano chunk must save PNGs for the exposure comparison" =
    grepl("save_hpiv3_ranked_comparisons(", report_text, fixed = TRUE) &&
      grepl('hpiv3_ranked_exposure_results,', report_text, fixed = TRUE) &&
      grepl('hpiv3_ranked_exposure_levels,', report_text, fixed = TRUE) &&
      grepl('"hpiv3_rb_", category_label, "_"', report_text, fixed = TRUE) &&
      grepl('"hpiv3_volc_", category_label, "_"', report_text, fixed = TRUE),
  "Ranked bar/volcano filenames must use the short, shared dimension-tag convention" =
    grepl("stratum_tag <- hpiv3_dimension_tag(", report_text, fixed = TRUE) &&
      !grepl("_hormone-", report_text, fixed = TRUE) &&
      !grepl("_timepoint-", report_text, fixed = TRUE),
  "Ranked plots must filter and encode airway, hormone, timepoint, and infection strata" =
    grepl('c\\("AIRWAY", "HORMONE", "TIMEPOINT", "INFECTION"\\)', report_text) &&
      grepl("chart_filters$SEX <- sex_value", report_text, fixed = TRUE),
  "Euler and heatmap outputs must reuse the shared short dimension-tag naming convention" =
    grepl("euler_stem <- paste0(", report_text, fixed = TRUE) &&
      grepl("hpiv3_dimension_tag(", report_text, fixed = TRUE) &&
      grepl("heatmap_tag <- hpiv3_dimension_tag(", report_text, fixed = TRUE) &&
      grepl('"hpiv3_heatmap_matrix_"', report_text, fixed = TRUE) &&
      grepl('"hpiv3_heatmap_"', report_text, fixed = TRUE),
  "Report must expose exposure/Euler diagnostics and a configurable Euler minimum" =
    grepl("HPIV3_EULER_MIN_SIG", report_text, fixed = TRUE) &&
      grepl("hpiv3_euler_skipped_strata.csv", report_text, fixed = TRUE) &&
      grepl("Exposure summary:", report_text, fixed = TRUE) &&
      !grepl("message(\n     \"Protein Euler plot skipped", report_text, fixed = TRUE) &&
      grepl("prepare_hpiv3_ranked_exposure_results(", report_text, fixed = TRUE),
  "Protein Euler outputs must use infection membership and export CSVs" =
    grepl('group_var = "INFECTION"', report_text, fixed = TRUE) &&
      grepl("plot_unique_protein_euler(", report_text, fixed = TRUE) &&
      grepl('"_membership.csv"', report_text, fixed = TRUE) &&
      grepl("at least two significant proteins per infection stratum", report_text, fixed = TRUE),
  "Session summary must describe the new short-tag ranked-bar/volcano/euler/heatmap filenames" =
    grepl("hpiv3_rb_<comparison>", report_text, fixed = TRUE) &&
      grepl("hpiv3_volc_<comparison>", report_text, fixed = TRUE) &&
      grepl("hpiv3_heatmap_matrix_aw-", report_text, fixed = TRUE) &&
      grepl("hpiv3_heatmap_aw-", report_text, fixed = TRUE),
  "RNA-seq report must reuse the shared short dimension-tag naming convention for Euler outputs" =
    grepl("hpiv3_dimension_tag(", rnaseq_report_text, fixed = TRUE) &&
      !grepl('"_airway-"', rnaseq_report_text, fixed = TRUE) &&
      grepl("hpiv3_euler_rnaseq_<peat|pine>_aw-", rnaseq_report_text, fixed = TRUE),
  "RNA-seq report must still generate EXPOSURE-family volcano/Euler outputs alongside HORMONE/INFECTION" =
    grepl('family == "EXPOSURE"', rnaseq_report_text, fixed = TRUE) &&
      grepl("plot_volcano_deg(", rnaseq_report_text, fixed = TRUE) &&
      grepl('RNA_FAMILIES <- c("EXPOSURE", "HORMONE", "INFECTION"', rnaseq_report_text, fixed = TRUE)
)

# Optional runtime regression checks for censor substitution, post-baseline
# count checks, pairing fallback, and the exposure-by-infection contrast.
runtime_pkgs <- c("dplyr", "tibble", "lme4", "emmeans", "pbkrtest")
if (all(vapply(runtime_pkgs, requireNamespace, logical(1), quietly = TRUE))) {
  suppressPackageStartupMessages(library(dplyr))
  data_env <- new.env(parent = globalenv())
  sys.source(data_path, envir = data_env)
  converted <- data_env$coerce_hpiv3_protein_values(
    data.frame(
      LLOD_REF = c("<LLOD", "2"),
      LLOD_PROXY = c("<LLOD", "4"),
      ULOD_PROXY = c(">ULOD", "8"),
      NO_BOUND = c("<LLOD", NA),
      stringsAsFactors = FALSE
    ),
    c("LLOD_REF", "LLOD_PROXY", "ULOD_PROXY", "NO_BOUND"),
    llod_table = data.frame(Analyte = "LLOD_REF", LLOD = 1),
    ulod_table = NULL
  )
  stopifnot(
    "Reference LLOD must be imputed at LLOD/sqrt(2)" =
      isTRUE(all.equal(converted$data$LLOD_REF[[1]], 1 / sqrt(2))),
    "Missing LLOD must use the minimum quantified value proxy" =
      isTRUE(all.equal(converted$data$LLOD_PROXY[[1]], 4 / sqrt(2))),
    "Missing ULOD must use the maximum quantified value proxy" =
      isTRUE(all.equal(converted$data$ULOD_PROXY[[1]], 8 * sqrt(2))),
    "Censored values without any usable bound must remain missing" =
      is.na(converted$data$NO_BOUND[[1]]) &&
        converted$qc$n_censored_unimputed[[4]] == 1L
  )

  model_env <- new.env(parent = globalenv())
  model_exprs <- parse(file = analysis_path)
  needed_models <- c(
    "HPIV3_LMER_DF_METHOD", "analysis_emmeans", "fit_hpiv3_paired_model",
    "fit_hpiv3_infection_models", "fit_hpiv3_sex_models",
    "fit_hpiv3_hormone_models", "fit_hpiv3_exposure_models",
    "fit_hpiv3_exposure_infection_interactions"
  )
  for (expr in model_exprs) {
    if (
      is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% needed_models
    ) {
      eval(expr, envir = model_env)
    }
  }
  interaction_fixture <- expand.grid(
    PATIENTCODE = paste0("D", seq_len(8)),
    EXPOSURE = c("PBS_Control", "Peat_25"),
    INFECTION = c("NONE", "HPIV3"),
    AIRWAY = "LAE", HORMONE = "NONE", TIMEPOINT = "24", SEX = "M",
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      donor_effect = as.numeric(factor(PATIENTCODE)) / 20,
      log2_value = 5 + donor_effect +
        0.5 * (EXPOSURE == "Peat_25") +
        0.7 * (INFECTION == "HPIV3") +
        0.9 * (EXPOSURE == "Peat_25" & INFECTION == "HPIV3") +
        0.02 * (as.numeric(interaction(PATIENTCODE, EXPOSURE, INFECTION)) %% 3 - 1),
      PROTEIN = 2^log2_value
    ) %>%
    dplyr::bind_rows(
      data.frame(
        PATIENTCODE = "D9", EXPOSURE = "Peat_25", INFECTION = "HPIV3",
        AIRWAY = "LAE", HORMONE = "NONE", TIMEPOINT = "24", SEX = "M",
        PROTEIN = 2^7, stringsAsFactors = FALSE
      )
    )
  interaction_results <- model_env$fit_hpiv3_exposure_infection_interactions(
    interaction_fixture, "PROTEIN",
    protein_status = data.frame(
      PROTEIN = "PROTEIN", included = TRUE, exclusion_reason = NA_character_
    ),
    pseudocount = 0.01, alpha_q = 0.05, min_nonmissing_per_group = 2L,
    pbs_level = "PBS_Control"
  )
  modeled_interaction <- interaction_results %>%
    dplyr::filter(SEX == "M", model_status %in% c("modeled", "modeled_paired_fallback"))
  stopifnot(
    "Interaction model must yield a positive differential exposure response" =
      nrow(modeled_interaction) == 1L && modeled_interaction$estimate[[1]] > 0.5,
    "Mixed-model interaction output must report Kenward-Roger df" =
      modeled_interaction$df_method[[1]] == "kenward-roger" &&
        is.finite(modeled_interaction$df[[1]]),
    "Interaction fit must exclude exposure observations lacking same-infection PBS" =
      modeled_interaction$n_samples[[1]] == 32L
  )
  protein_status <- data.frame(
    PROTEIN = "PROTEIN", included = TRUE, exclusion_reason = NA_character_
  )
  infection_results <- model_env$fit_hpiv3_infection_models(
    interaction_fixture %>% dplyr::filter(EXPOSURE == "PBS_Control"),
    "PROTEIN", protein_status = protein_status,
    pseudocount = 0.01, alpha_q = 0.05
  )
  exposure_results <- model_env$fit_hpiv3_exposure_models(
    interaction_fixture, "PROTEIN", protein_status = protein_status,
    pseudocount = 0.01, alpha_q = 0.05
  )
  hormone_fixture <- expand.grid(
    PATIENTCODE = paste0("D", seq_len(8)),
    HORMONE = c("NONE", "E2"), AIRWAY = "LAE", TIMEPOINT = "24",
    EXPOSURE = "PBS_Control", INFECTION = "NONE", SEX = "M",
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      donor = as.numeric(factor(PATIENTCODE)),
      log2_value = 5 + donor / 20 + 0.4 * (HORMONE == "E2") +
        0.02 * (donor %% 3 - 1) * ifelse(HORMONE == "E2", 1, -1),
      PROTEIN = 2^log2_value
    )
  hormone_results <- model_env$fit_hpiv3_hormone_models(
    hormone_fixture, "PROTEIN", protein_status = protein_status,
    pseudocount = 0.01, alpha_q = 0.05
  )
  sex_fixture <- expand.grid(
    PATIENTCODE = paste0("D", seq_len(8)),
    SEX = c("F", "M"), AIRWAY = "LAE", HORMONE = "NONE",
    TIMEPOINT = "24", EXPOSURE = "PBS_Control", INFECTION = "NONE",
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      donor = as.numeric(factor(PATIENTCODE)),
      log2_value = 5 + donor / 20 + 0.3 * (SEX == "M") +
        0.02 * (donor %% 3 - 1) * ifelse(SEX == "M", 1, -1),
      PROTEIN = 2^log2_value
    )
  sex_results <- model_env$fit_hpiv3_sex_models(
    sex_fixture, "PROTEIN", protein_status = protein_status,
    pseudocount = 0.01, alpha_q = 0.05
  )
  for (model_output in list(
    infection_results, exposure_results, hormone_results, sex_results
  )) {
    mixed_rows <- model_output %>% dplyr::filter(model_type == "lmer")
    stopifnot(
      "Every mixed-model contrast must use Kenward-Roger and report finite df" =
        nrow(mixed_rows) > 0L &&
          all(mixed_rows$df_method == "kenward-roger") &&
          all(is.finite(mixed_rows$df))
    )
  }

  thin_fixture <- data.frame(
    PATIENTCODE = c("D1", "D2", "D3", "D1", "D4"),
    EXPOSURE = c("PBS_Control", "PBS_Control", "PBS_Control", "Peat_25", "Peat_25"),
    INFECTION = "NONE", AIRWAY = "LAE", HORMONE = "NONE", TIMEPOINT = "24",
    SEX = "M", PROTEIN = c(1, 2, 3, 4, 5), stringsAsFactors = FALSE
  )
  thin_result <- model_env$fit_hpiv3_exposure_models(
    thin_fixture, "PROTEIN",
    protein_status = data.frame(
      PROTEIN = "PROTEIN", included = TRUE, exclusion_reason = NA_character_
    ),
    pseudocount = 0.01, alpha_q = 0.05, min_nonmissing_per_group = 2L
  )
  stopifnot(
    "Exposure counts must be re-checked after PBS-donor filtering" =
      all(grepl("after PBS-donor filtering", thin_result$failure_reason))
  )
  unpaired_fixture <- expand.grid(
    PATIENTCODE = paste0("D", seq_len(8)),
    EXPOSURE = c("PBS_Control", "Peat_25"),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(
      (as.numeric(factor(PATIENTCODE)) <= 4 & EXPOSURE == "PBS_Control") |
        (as.numeric(factor(PATIENTCODE)) > 4 & EXPOSURE == "Peat_25")
    ) %>%
    dplyr::mutate(log2_value = seq_len(dplyr::n()))
  unpaired <- model_env$fit_hpiv3_paired_model(
    unpaired_fixture,
    log2_value ~ EXPOSURE, pairing_factor = "EXPOSURE"
  )
  single_paired <- model_env$fit_hpiv3_paired_model(
    interaction_fixture %>% dplyr::filter(PATIENTCODE == "D1") %>%
      dplyr::mutate(PAIR_LEVEL = interaction(EXPOSURE, INFECTION)),
    log2_value ~ EXPOSURE, pairing_factor = "PAIR_LEVEL"
  )
  fallback_fixture <- interaction_fixture %>%
    dplyr::filter(PATIENTCODE != "D9") %>%
    dplyr::mutate(
      log2_value = 5 + 0.5 * (EXPOSURE == "Peat_25") +
        0.7 * (INFECTION == "HPIV3") +
        0.9 * (EXPOSURE == "Peat_25" & INFECTION == "HPIV3") +
        0.05 * (as.numeric(factor(PATIENTCODE)) %% 3 - 1) *
          ifelse(as.numeric(factor(interaction(EXPOSURE, INFECTION))) %% 2 == 0, 1, -1),
      PAIR_LEVEL = interaction(EXPOSURE, INFECTION)
    )
  paired_fallback <- model_env$fit_hpiv3_paired_model(
    fallback_fixture, log2_value ~ EXPOSURE * INFECTION,
    pairing_factor = "PAIR_LEVEL"
  )
  stopifnot(
    "Unpaired lm is reserved for data with no cross-level donor pairing" =
      unpaired$model_type == "lm_unpaired" && !unpaired$paired,
    "One paired donor must not silently trigger an unpaired fallback" =
      is.null(single_paired$fit) && single_paired$paired &&
        single_paired$model_status == "failed",
    "Singular mixed models must use the donor-fixed paired fallback" =
      paired_fallback$model_type == "lm_donor_fixed" && paired_fallback$paired
  )
}

# ---- Behavioural checks (need ggplot2/dplyr/tibble/ggrepel; skipped otherwise) ----
behaviour_pkgs <- c("ggplot2", "dplyr", "tibble", "ggrepel")
if (all(vapply(behaviour_pkgs, requireNamespace, logical(1), quietly = TRUE))) {
  suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })
  ALPHA_Q <- 0.05
  UP_COLOR_DEFAULT <- "#D7191C"
  DOWN_COLOR_DEFAULT <- "#2C7BB6"
  eval(parse(file = plots_path))
  PBS_LEVEL <- "PBS_Control"
  prep_fn_start <- regexpr("prepare_hpiv3_ranked_exposure_results <- function", analysis_text, fixed = TRUE)
  prep_fn_end <- regexpr("summarize_hpiv3_strata <- function", analysis_text, fixed = TRUE)
  eval(parse(text = substr(analysis_text, prep_fn_start, prep_fn_end - 1L)))

  # BH q within one stratum; two panels (infection NONE/HPIV3) with different boundaries.
  make_panel <- function(infection, p, n_sig_expected) {
    q <- p.adjust(p, method = "fdr")
    tibble::tibble(
      PROTEIN = paste0(infection, "_", seq_along(p)),
      INFECTION = infection, SEX = "All", comparison = "EXPOSURE",
      contrast = "Peat_25 - PBS_Control", target_exposure = "Peat_25",
      estimate = seq(-1, 1, length.out = length(p)),
      p.value = p, q.value = q, significant = q < ALPHA_Q
    )
  }
  panel_none <- make_panel("NONE", c(0, 1e-6, 1e-4, 0.03, 0.2, 0.7), 3)
  panel_hpiv3 <- make_panel("HPIV3", c(1e-3, 0.002, 0.04, 0.5, 0.9), 2)
  vol_data <- dplyr::bind_rows(panel_none, panel_hpiv3)

  groups <- vol_data$INFECTION
  thr <- hpiv3_volcano_thresholds(vol_data, groups)
  stopifnot(
    "Threshold table must have one row per panel" = nrow(thr) == 2L,
    "Threshold must be finite even when p = 0 appears in a panel" =
      all(is.finite(thr$threshold_y)),
    "p = 0 must give a finite -log10(p)" = is.finite(hpiv3_neg_log10_p(0)),
    "NA/Inf p-values must stay NA" = all(is.na(hpiv3_neg_log10_p(c(NA, Inf, NaN))))
  )
  for (g in thr$.panel_group) {
    sub <- vol_data[vol_data$INFECTION == g, ]
    y <- hpiv3_neg_log10_p(sub$p.value)
    ty <- thr$threshold_y[thr$.panel_group == g]
    stopifnot(
      "Threshold line must sit at or below every significant point" = all(y[sub$significant] >= ty - 1e-12),
      "Threshold line must sit above every non-significant point" = all(y[!sub$significant] < ty)
    )
  }
  # The raw-p boundary is NOT -log10(ALPHA_Q) in general.
  stopifnot("Boundary must derive from significant proteins, not -log10(ALPHA_Q)" =
    !isTRUE(all.equal(thr$threshold_y[thr$.panel_group == "NONE"], -log10(ALPHA_Q))))

  # Panels without significant hits get no boundary.
  none_sig <- hpiv3_volcano_thresholds(
    dplyr::mutate(panel_hpiv3, significant = FALSE), panel_hpiv3$INFECTION
  )
  stopifnot("No significant proteins -> NA threshold" = is.na(none_sig$threshold_y))

  # The plot must contain a dotted hline layer whose y values equal the thresholds.
  vol_single <- panel_none
  filters <- list(INFECTION = "NONE", SEX = "All", target_exposure = "Peat_25")
  vol_plot <- plot_hpiv3_volcano(vol_single, "exposure", filters = filters)
  hline_layers <- Filter(function(l) inherits(l$geom, "GeomHline"), vol_plot$layers)
  stopifnot("Volcano must have exactly one boundary hline layer" = length(hline_layers) == 1L)
  vline_layers <- Filter(function(l) inherits(l$geom, "GeomVline"), vol_plot$layers)
  stopifnot("Volcano must not draw a log2FC vline" = length(vline_layers) == 0L)
  hline_y <- hline_layers[[1]]$data$threshold_y
  built <- ggplot2::ggplot_build(vol_plot)
  pts <- built$data[[which(vapply(vol_plot$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1)))[1]]]
  stopifnot(
    "Boundary hline y must match hpiv3_volcano_thresholds()" =
      isTRUE(all.equal(hline_y, thr$threshold_y[thr$.panel_group == "NONE"])),
    "Plotted y of p = 0 must be finite" = all(is.finite(pts$y))
  )
  # Facetted volcano also builds.
  invisible(ggplot2::ggplot_build(plot_hpiv3_volcano(
    vol_data, "exposure", filters = list(SEX = "All", target_exposure = "Peat_25")
  )))
  # Ranked bars build for exposure results too.
  invisible(ggplot2::ggplot_build(plot_hpiv3_ranked_bars(vol_single, "exposure", filters = filters)))

  # Exposure contrast parsing: sign, label, PBS-first flip, non-Peat/Pine drop.
  raw_exp <- tibble::tibble(
    comparison = "EXPOSURE",
    contrast = c("Peat_25 - PBS_Control", "PBS_Control - Pine_25", "Eucalyptus_25 - PBS_Control", NA),
    estimate = c(1, 2, 3, NA)
  )
  prep <- prepare_hpiv3_ranked_exposure_results(raw_exp)
  stopifnot(
    "Only Peat/Pine rows survive" = identical(prep$target_exposure, c("Peat_25", "Pine_25")),
    "PBS-first contrast flips the estimate sign" = identical(prep$estimate, c(1, -2)),
    "Contrast labels are normalised to <target> - PBS_Control" =
      identical(prep$contrast, c("Peat_25 - PBS_Control", "Pine_25 - PBS_Control"))
  )
} else {
  cat("NOTE: skipping behavioural volcano checks (missing:",
      paste(behaviour_pkgs[!vapply(behaviour_pkgs, requireNamespace, logical(1), quietly = TRUE)], collapse = ", "), ")\n")
}

cat("PASS: HPIV3 PBS-referenced contrasts, infection-stratified models, ranked plots, naming, and Euler wiring checks passed.\n")
