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
plots_path <- file.path(root, "functions_plots.R")
report_path <- file.path(root, "02_hpiv3_analysis.Rmd")
loader_path <- file.path(root, "_load_all.R")
analysis_text <- paste(readLines(analysis_path, warn = FALSE), collapse = "\n")
plots_text <- paste(readLines(plots_path, warn = FALSE), collapse = "\n")
report_text <- paste(readLines(report_path, warn = FALSE), collapse = "\n")
loader_text <- paste(readLines(loader_path, warn = FALSE), collapse = "\n")

for (path in c(analysis_path, plots_path, loader_path)) parse(file = path)
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
euler_chunk_start <- match("```{r hpiv3-exposure-euler}", rmd_lines)
euler_chunk_end <- if (is.na(euler_chunk_start)) {
  NA_integer_
} else {
  which(seq_along(rmd_lines) > euler_chunk_start & rmd_lines == "```")[1]
}
stopifnot(!is.na(euler_chunk_start), !is.na(euler_chunk_end))
parse(text = rmd_lines[(euler_chunk_start + 1L):(euler_chunk_end - 1L)])

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
    grepl("lme4::lmer\\(log2_value ~ HORMONE", analysis_text) &&
      grepl("stats::lm\\(log2_value ~ HORMONE", analysis_text),
  "Exposure models must create pooled All-sex input rows" =
    grepl('mutate\\(SEX = "All"\\)', analysis_text),
  "Exposure contrasts must compare treatments with PBS only" =
    grepl('method = "trt.vs.ctrl"', exposure_function_text) &&
      grepl('match\\("PBS_Control", exposure_levels\\)', exposure_function_text) &&
      !grepl('method = "pairwise"', exposure_function_text),
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
      grepl("hpiv3_volcano_", report_text, fixed = TRUE),
  "Ranked plots must filter and encode airway, hormone, timepoint, and infection strata" =
    grepl('c\\("AIRWAY", "HORMONE", "TIMEPOINT", "INFECTION"\\)', report_text) &&
      grepl("chart_filters\\$SEX <- sex_value", report_text, fixed = TRUE) &&
      grepl('"_sex-", tolower\\(sex_value\\)', report_text),
  "Protein Euler outputs must use infection membership and export CSVs" =
    grepl('group_var = "INFECTION"', report_text, fixed = TRUE) &&
      grepl("plot_unique_protein_euler(", report_text, fixed = TRUE) &&
      grepl('"_membership.csv"', report_text, fixed = TRUE) &&
      grepl("at least two significant proteins per infection stratum", report_text, fixed = TRUE)
)

cat("PASS: HPIV3 PBS-referenced contrasts, infection-stratified models, ranked plots, and Euler wiring checks passed.\n")
