#!/usr/bin/env Rscript
# Lightweight, dependency-free validation fixture for the DREAM/voom
# explicit-effects model builder in `03_rnaseq_hpiv3_analysis.Rmd`.
#
# This script does NOT require edgeR/limma/variancePartition. It validates
# the *design* side of the refactor directly with base R:
#
#   1. It extracts the actual contrast-weight helper functions
#      (`observed_factor_levels`, `factor_contrast_options`,
#      `contrast_label_for`, `build_contrast_weights`) verbatim from the Rmd,
#      so this check stays in sync with the real implementation.
#   2. It reproduces an unbalanced HORMONE x EXPOSURE sample layout like the
#      one described in the issue: EXPOSURE = "PEAT" has zero HORMONE = "E2"
#      samples.
#   3. It confirms that a HORMONE main-effect model scoped only to HORMONE,
#      with EXPOSURE included as a plain additive adjustment term (exactly
#      what `get_fit()` builds for `TESTS = "HORMONE"`), is estimable (full
#      rank) and yields a `HORMONE=E2-NONE` contrast -- without requiring
#      every EXPOSURE x HORMONE cell to be populated.
#   4. For contrast, it shows that the OLD behavior -- crossing every
#      varying factor (HORMONE x EXPOSURE) into one saturated GROUP -- is
#      rank-deficient for this same data, which is exactly the bug this
#      refactor fixes.
#   5. It also confirms an explicitly requested interaction test
#      (`TESTS = "HORMONE:EXPOSURE"`) correctly reports the contrasts that
#      ARE estimable, and fails only the specific contrast that genuinely
#      requires the missing HORMONE=E2 / EXPOSURE=PEAT cell.
#
# Run with: Rscript tests/rnaseq_dream_design_sanity_check.R

candidate_paths <- c(
  "03_rnaseq_hpiv3_analysis.Rmd",
  file.path("..", "03_rnaseq_hpiv3_analysis.Rmd")
)
rmd_path <- candidate_paths[file.exists(candidate_paths)][1]
stopifnot(
  "Cannot locate 03_rnaseq_hpiv3_analysis.Rmd; run this script from the repository root (or tests/)." =
    !is.na(rmd_path)
)
rmd_lines <- readLines(rmd_path, warn = FALSE)
rmd_text <- paste(rmd_lines, collapse = "\n")

selector_start <- match("```{r analysis-selector}", rmd_lines)
selector_end <- if (is.na(selector_start)) {
  NA_integer_
} else {
  closing_fences <- which(seq_along(rmd_lines) > selector_start & rmd_lines == "```")
  if (length(closing_fences) == 0) NA_integer_ else closing_fences[[1]]
}
stopifnot(
  "Cannot locate the analysis-selector chunk in the Rmd" =
    !is.na(selector_start) && !is.na(selector_end)
)
selector_env <- new.env(parent = globalenv())
eval(parse(text = rmd_lines[seq.int(selector_start + 1L, selector_end - 1L)]),
     envir = selector_env)

expected_family_specs <- list(
  EXPOSURE = list(
    infection_levels = "NONE",
    tests = "EXPOSURE",
    report_within = character()
  ),
  EXPOSURE_INFECTED = list(
    infection_levels = "HPIV3",
    tests = "EXPOSURE",
    report_within = character()
  ),
  HORMONE = list(
    infection_levels = "NONE",
    tests = "HORMONE",
    report_within = "EXPOSURE"
  ),
  INFECTION = list(
    infection_levels = c("NONE", "HPIV3"),
    tests = "INFECTION",
    report_within = "EXPOSURE"
  )
)
stopifnot(
  "Selector must define all four expected families" =
    identical(selector_env$RNA_FAMILIES, names(expected_family_specs)),
  "Family specs must preserve existing behavior and restrict EXPOSURE_INFECTED to HPIV3" =
    all(vapply(names(expected_family_specs), function(family) {
      identical(selector_env$family_spec(family), expected_family_specs[[family]])
    }, logical(1))),
  "Exposure contrasts must remain the same three pairwise comparisons" =
    identical(selector_env$EXPOSURE_CONTRASTS,
              c("PEAT-PBS", "PINE-PBS", "PEAT-PINE"))
)

campaign_specs <- selector_env$RNA_CAMPAIGN_SPECS
expected_spec_keys <- as.vector(outer(
  selector_env$RNA_AIRWAY_LEVELS,
  selector_env$RNA_SEX_STRATA,
  paste,
  sep = "_"
))
expected_spec_keys <- unlist(lapply(names(expected_family_specs), function(family) {
  paste(expected_spec_keys, family, sep = "_")
}), use.names = FALSE)
actual_spec_keys <- vapply(campaign_specs, `[[`, character(1), "spec_tag")
stopifnot(
  "Campaign grid must contain one spec per airway, sex stratum, and family" =
    length(campaign_specs) == 24L &&
    setequal(actual_spec_keys, expected_spec_keys),
  "Every EXPOSURE_INFECTED spec must select HPIV3 and use exposure contrasts" =
    all(vapply(Filter(function(spec) spec$FAMILY == "EXPOSURE_INFECTED",
                      campaign_specs), function(spec) {
      identical(spec$SUBSET$INFECTION, "HPIV3") &&
        identical(spec$TESTS, "EXPOSURE") &&
        identical(spec$REPORT_WITHIN, character()) &&
        identical(spec$EXPOSURE_CONTRASTS,
                  c("PEAT-PBS", "PINE-PBS", "PEAT-PINE"))
    }, logical(1)))
)

extract_block <- function(text, start_anchor, end_anchor) {
  start_i <- regexpr(start_anchor, text, fixed = TRUE)
  stopifnot("start anchor not found in Rmd" = start_i > 0)
  end_i <- regexpr(end_anchor, text, fixed = TRUE)
  stopifnot("end anchor not found in Rmd" = end_i > 0)
  substr(text, start_i, end_i - 1L)
}

helper_source <- extract_block(
  rmd_text,
  "# BEGIN contrast-weight helpers",
  "# END contrast-weight helpers"
)

# Selector constants referenced by the extracted helpers.
RNA_HORMONE_LEVELS <- c("NONE", "E2")
RNA_EXPOSURE_LEVELS <- c("PBS", "PEAT", "PINE")
RNA_FACTOR_LEVELS <- list(HORMONE = RNA_HORMONE_LEVELS, EXPOSURE = RNA_EXPOSURE_LEVELS)
EXPOSURE_CONTRASTS <- c("PEAT-PBS", "PINE-PBS", "PEAT-PINE")

eval(parse(text = helper_source), envir = globalenv())
stopifnot(
  is.function(observed_factor_levels),
  is.function(factor_contrast_options),
  is.function(contrast_label_for),
  is.function(build_contrast_weights)
)

# ---------------------------------------------------------------------------
# Unbalanced HORMONE x EXPOSURE sample layout (mirrors the issue's example):
# EXPOSURE = "PEAT" never co-occurs with HORMONE = "E2".
# ---------------------------------------------------------------------------
model_data <- data.frame(
  HORMONE = factor(
    c("NONE", "NONE", "NONE", "E2", "E2", "E2", "NONE", "NONE", "E2", "E2", "NONE", "E2"),
    levels = RNA_HORMONE_LEVELS
  ),
  EXPOSURE = factor(
    c("PBS", "PBS", "PEAT", "PBS", "PBS", "PINE", "PEAT", "PINE", "PINE", "PBS", "PINE", "PBS"),
    levels = RNA_EXPOSURE_LEVELS
  ),
  stringsAsFactors = FALSE
)
cat("Observed HORMONE x EXPOSURE cell counts:\n")
print(table(model_data$HORMONE, model_data$EXPOSURE))
stopifnot(
  "Fixture must leave EXPOSURE=PEAT with zero HORMONE=E2 samples" =
    sum(model_data$HORMONE == "E2" & model_data$EXPOSURE == "PEAT") == 0
)

# ---------------------------------------------------------------------------
# 1) NEW behavior: TESTS = "HORMONE" -> GROUP scoped to HORMONE only, with
#    EXPOSURE as a plain additive adjustment term (as built by get_fit() in
#    the Rmd for a single-factor TESTS entry).
# ---------------------------------------------------------------------------
core_factors <- "HORMONE"
group_raw <- interaction(model_data[[core_factors]], drop = TRUE, sep = "_")
model_data$GROUP <- group_raw
group_cells <- unique(model_data[, c(core_factors, "GROUP"), drop = FALSE])

design <- stats::model.matrix(~ 0 + GROUP + EXPOSURE, data = model_data)
cat("\nScoped HORMONE-only design columns:", paste(colnames(design), collapse = ", "), "\n")
stopifnot(
  "Additive HORMONE + EXPOSURE design must be full rank despite the missing cell" =
    qr(design)$rank == ncol(design)
)

effect_levels <- observed_factor_levels(model_data, "HORMONE")
vectors <- list(factor_contrast_options("HORMONE", effect_levels)[[1]])
weights_result <- build_contrast_weights(
  effect_factors = "HORMONE", vectors = vectors, conditions = character(),
  model_factors = core_factors, group_cells = group_cells
)
stopifnot(
  "HORMONE=E2-NONE main-effect contrast must be estimable under the additive model" =
    !is.null(weights_result$weights)
)
cat("HORMONE=E2-NONE contrast weights (additive model):\n")
print(weights_result$weights)
cat("PASS: HORMONE main effect is estimable with EXPOSURE as an additive adjustment,\n",
    "      even though EXPOSURE=PEAT has no HORMONE=E2 samples.\n", sep = "")

# ---------------------------------------------------------------------------
# 2) OLD behavior for contrast: crossing every varying factor
#    (HORMONE x EXPOSURE) into one saturated GROUP is rank-deficient for this
#    same data -- this is exactly the bug the refactor fixes.
# ---------------------------------------------------------------------------
saturated_raw <- interaction(
  model_data$HORMONE, model_data$EXPOSURE, drop = FALSE, lex.order = TRUE, sep = "_"
)
saturated_design <- stats::model.matrix(~ 0 + saturated_raw)
cat("\nOld saturated-GROUP design rank:", qr(saturated_design)$rank,
    "of", ncol(saturated_design), "columns (drop = FALSE keeps the empty cell)\n")
stopifnot(
  "Saturated cross of HORMONE x EXPOSURE must be rank-deficient for this unbalanced layout" =
    qr(saturated_design)$rank < ncol(saturated_design)
)
cat("CONFIRMED: the saturated all-factors-crossed GROUP design is rank-deficient,\n",
    "           demonstrating why the old model builder could reject this comparison.\n", sep = "")

# ---------------------------------------------------------------------------
# 3) Explicit interaction test: TESTS = "HORMONE:EXPOSURE". GROUP is scoped to
#    just the two interacting factors with only observed combinations
#    (drop = TRUE), so unrelated contrasts still succeed, and only the
#    contrast that truly requires the missing cell fails.
# ---------------------------------------------------------------------------
core_factors2 <- c("HORMONE", "EXPOSURE")
interaction_raw <- interaction(
  model_data[core_factors2], drop = TRUE, lex.order = TRUE, sep = "_"
)
model_data$GROUP2 <- interaction_raw
group_cells2 <- unique(model_data[, c(core_factors2, "GROUP2"), drop = FALSE])
names(group_cells2)[names(group_cells2) == "GROUP2"] <- "GROUP"
design2 <- stats::model.matrix(~ 0 + GROUP2, data = model_data)
stopifnot(
  "Observed-cells-only interaction GROUP (drop = TRUE) must be full rank" =
    qr(design2)$rank == ncol(design2)
)

# A DiD contrast entirely within PBS/PINE (both HORMONE levels present) must
# be estimable:
pbs_pine_weights <- build_contrast_weights(
  effect_factors = "HORMONE",
  vectors = list(stats::setNames(c(1, -1), c("E2", "NONE"))),
  conditions = c(EXPOSURE = "PINE"),
  model_factors = core_factors2, group_cells = group_cells2
)
stopifnot(
  "HORMONE simple effect within EXPOSURE=PINE must be estimable" =
    !is.null(pbs_pine_weights$weights)
)

# A simple effect that truly requires the missing HORMONE=E2/EXPOSURE=PEAT
# cell must fail with a clear reason, not a silent/incorrect result:
peat_weights <- build_contrast_weights(
  effect_factors = "HORMONE",
  vectors = list(stats::setNames(c(1, -1), c("E2", "NONE"))),
  conditions = c(EXPOSURE = "PEAT"),
  model_factors = core_factors2, group_cells = group_cells2
)
stopifnot(
  "HORMONE simple effect within EXPOSURE=PEAT must fail (cell genuinely missing)" =
    is.null(peat_weights$weights)
)
cat("\nPASS: explicit HORMONE:EXPOSURE interaction test still returns estimable\n",
    "      simple effects (e.g. within EXPOSURE=PINE) and correctly/only skips\n",
    "      the one contrast requiring the genuinely missing HORMONE=E2 x\n",
    "      EXPOSURE=PEAT cell (reason: '", peat_weights$reason, "').\n", sep = "")

cat("\nAll RNA-seq DREAM design sanity checks passed.\n")
