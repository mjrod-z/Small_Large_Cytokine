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

stopifnot(
  "Hormone models must test E2 against NONE" =
    grepl('control_level = "NONE"', analysis_text) &&
      grepl('case_level = "E2"', analysis_text),
  "Hormone models must pool across hormone and stratify by airway/timepoint/exposure" =
    grepl("distinct\\(AIRWAY, TIMEPOINT, EXPOSURE\\)", analysis_text),
  "Hormone models must emit pooled and sex-specific groups" =
    grepl('sex_groups <- c\\("All", intersect\\(c\\("F", "M"\\)', analysis_text),
  "Hormone models must retain mixed and fallback linear model paths" =
    grepl("lme4::lmer\\(log2_value ~ HORMONE", analysis_text) &&
      grepl("stats::lm\\(log2_value ~ HORMONE", analysis_text),
  "Exposure models must create pooled All-sex input rows" =
    grepl('mutate\\(SEX = "All"\\)', analysis_text),
  "Ranked bars must use estimated effects and direction/significance fills" =
    grepl("geom_col\\(width = 0.75\\)", plots_text) &&
      grepl('"Higher, significant"', plots_text) &&
      grepl('"Lower, significant"', plots_text),
  "Ranked plot helper must be present in bootstrap checks" =
    grepl('"plot_hpiv3_ranked_bars"', loader_text),
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
    ))
)

cat("PASS: HPIV3 hormone, pooled exposure, ranked plot, and report wiring checks passed.\n")
