#!/usr/bin/env Rscript
# Focused checks for RNA-seq sex-effect plots and profile export wiring.

candidate_plot_paths <- c("functions_plots.R", file.path("..", "functions_plots.R"))
plot_path <- candidate_plot_paths[file.exists(candidate_plot_paths)][1]
stopifnot(
  "Cannot locate functions_plots.R; run this script from the repository root (or tests/)." =
    !is.na(plot_path)
)
candidate_rmd_paths <- c(
  "03_rnaseq_hpiv3_analysis.Rmd",
  file.path("..", "03_rnaseq_hpiv3_analysis.Rmd")
)
rmd_path <- candidate_rmd_paths[file.exists(candidate_rmd_paths)][1]
stopifnot(
  "Cannot locate 03_rnaseq_hpiv3_analysis.Rmd; run this script from the repository root (or tests/)." =
    !is.na(rmd_path)
)
candidate_loader_paths <- c("_load_all.R", file.path("..", "_load_all.R"))
loader_path <- candidate_loader_paths[file.exists(candidate_loader_paths)][1]
stopifnot(
  "Cannot locate _load_all.R; run this script from the repository root (or tests/)." =
    !is.na(loader_path)
)

required_packages <- c("dplyr", "tidyr", "ggplot2", "knitr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

plot_env <- new.env(parent = globalenv())
sys.source(plot_path, envir = plot_env)
stopifnot(
  "The sex-comparison helper must be defined" =
    exists("plot_hpiv3_rnaseq_sex_comparison", envir = plot_env)
)

toy_results <- expand.grid(
  GENEID = c("GENE_A", "GENE_B", "GENE_C"),
  sex_stratum = c("M", "F"),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    airway = "LAE",
    family = "INFECTION",
    contrast_label = "HPIV3-NONE within PBS",
    logFC = dplyr::case_when(
      .data$GENEID == "GENE_A" & .data$sex_stratum == "M" ~ 3,
      .data$GENEID == "GENE_A" & .data$sex_stratum == "F" ~ 1,
      .data$GENEID == "GENE_B" & .data$sex_stratum == "M" ~ -2,
      .data$GENEID == "GENE_B" & .data$sex_stratum == "F" ~ -1,
      .data$GENEID == "GENE_C" & .data$sex_stratum == "M" ~ 0.5,
      TRUE ~ 0.25
    ),
    significant = abs(.data$logFC) >= 1
  )

plot <- plot_env$plot_hpiv3_rnaseq_sex_comparison(
  toy_results, airway = "LAE", family = "INFECTION", top_n = 2L
)
stopifnot(
  "The helper must return a ggplot with two sex points per selected gene" =
    inherits(plot, "ggplot") && nrow(plot$data) == 4L,
  "The helper must use paired-effect segments" =
    nrow(plot$layers[[1]]$data) == 2L
)

male_only_plot <- plot_env$plot_hpiv3_rnaseq_sex_comparison(
  dplyr::filter(toy_results, .data$sex_stratum == "M"),
  airway = "LAE",
  family = "INFECTION"
)
stopifnot(
  "Missing female results must produce an informative empty plot" =
    inherits(male_only_plot, "ggplot") &&
      grepl("not both available", male_only_plot$labels$subtitle)
)

one_paired_gene <- toy_results %>%
  dplyr::filter(.data$GENEID == "GENE_A" |
                  (.data$GENEID == "GENE_B" & .data$sex_stratum == "M"))
one_gene_plot <- plot_env$plot_hpiv3_rnaseq_sex_comparison(
  one_paired_gene, airway = "LAE", family = "INFECTION"
)
stopifnot(
  "Fewer than two paired genes must produce an informative empty plot" =
    inherits(one_gene_plot, "ggplot") &&
      grepl("Fewer than two genes", one_gene_plot$labels$subtitle)
)

rmd_text <- paste(readLines(rmd_path, warn = FALSE), collapse = "\n")
loader_text <- paste(readLines(loader_path, warn = FALSE), collapse = "\n")
volcano_position <- regexpr("```\\{r volcano-and-overlap\\}", rmd_text)
sex_section_position <- regexpr("# Sex-stratified effect comparisons", rmd_text)
profile_section_position <- regexpr("# Full gene-by-contrast profile matrix", rmd_text)
stopifnot(
  "New report sections must follow the existing volcano output section" =
    volcano_position > 0 &&
      sex_section_position > volcano_position &&
      profile_section_position > sex_section_position,
  "The report must export a wide logFC profile with interpretable condition columns" =
    grepl("hpiv3_rnaseq_full_gene_contrast_profile\\.csv", rmd_text) &&
      grepl("profile_column = paste\\(", rmd_text, fixed = FALSE) &&
      grepl("tidyr::pivot_wider\\(names_from = profile_column, values_from = logFC\\)", rmd_text),
  "The project loader must require the new plot helper" =
    grepl("plot_hpiv3_rnaseq_sex_comparison", loader_text, fixed = TRUE)
)

report_lines <- readLines(rmd_path, warn = FALSE)
profile_chunk_start <- match("```{r full-gene-contrast-profile}", report_lines)
profile_chunk_end <- if (is.na(profile_chunk_start)) {
  NA_integer_
} else {
  closing_fences <- which(
    seq_along(report_lines) > profile_chunk_start & report_lines == "```"
  )
  if (length(closing_fences) == 0) NA_integer_ else closing_fences[[1]]
}
stopifnot(
  "Cannot locate the full gene-by-contrast profile chunk" =
    !is.na(profile_chunk_start) && !is.na(profile_chunk_end)
)
profile_toy <- data.frame(
  GENEID = c("GENE_A", "GENE_B", "GENE_A", "GENE_C"),
  airway = "LAE",
  family = "INFECTION",
  sex_stratum = c("M", "M", "F", "F"),
  contrast_label = "INFECTION=HPIV3-NONE | within EXPOSURE=PBS",
  logFC = c(2, -1, 1.5, 0.2),
  significant = c(TRUE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)
profile_env <- new.env(parent = globalenv())
profile_env$rnaseq_model_results <- profile_toy
profile_env$RNA_OUTPUT_SUBDIR <- "rna_seq"
profile_env$save_table <- function(data, filename) {
  profile_env$saved_profile <- data
  profile_env$saved_filename <- filename
  invisible(TRUE)
}
eval(
  parse(text = report_lines[seq.int(profile_chunk_start + 1L, profile_chunk_end - 1L)]),
  envir = profile_env
)
expected_m_column <- paste(
  "LAE", "INFECTION", "M", profile_toy$contrast_label[[1]], sep = "__"
)
expected_f_column <- paste(
  "LAE", "INFECTION", "F", profile_toy$contrast_label[[1]], sep = "__"
)
stopifnot(
  "The profile export must contain every observed gene and condition" =
    nrow(profile_env$saved_profile) == 3L &&
      all(c("GENEID", expected_m_column, expected_f_column) %in%
            names(profile_env$saved_profile)),
  "Untested gene-condition pairs must be NA" =
    is.na(profile_env$saved_profile[[expected_m_column]][
      profile_env$saved_profile$GENEID == "GENE_C"
    ]) &&
      is.na(profile_env$saved_profile[[expected_f_column]][
        profile_env$saved_profile$GENEID == "GENE_B"
      ]),
  "The matrix must be written to the required RNA-seq filename" =
    identical(
      profile_env$saved_filename,
      file.path("rna_seq", "hpiv3_rnaseq_full_gene_contrast_profile.csv")
    )
)

cat("PASS: RNA-seq sex-effect plot and profile export checks passed.\n")
