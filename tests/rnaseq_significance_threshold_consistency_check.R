#!/usr/bin/env Rscript
# Dependency-free validation fixture for the global RNA-seq significance
# threshold wiring introduced to resolve the `ALPHA_Q` / `ADJ_P_CUTOFF` /
# `FDR_THRESHOLD` inconsistency.
#
# This script does NOT require dplyr/limma/variancePartition. It checks:
#   1. `config.R` defines a single canonical BH-FDR cutoff (`ALPHA_Q`), that
#      `ADJ_P_CUTOFF` is derived from it (not a duplicated literal), that the
#      conflicting unused `FDR_THRESHOLD <- 0.25` has been removed, and that
#      `LOG2FC_CUTOFF` remains the canonical effect-size cutoff.
#   2. The `significant` / `deg_class` rule actually implemented in
#      `03_rnaseq_hpiv3_analysis.Rmd` requires BOTH BH-FDR < cutoff AND
#      |log2FC| >= cutoff, and that `fdr_significant` (FDR-only) never
#      disagrees with `significant` in a way that would let a FDR-only hit
#      with a small effect size count as significant.
#
# Run with: Rscript tests/rnaseq_significance_threshold_consistency_check.R

candidate_config_paths <- c("config.R", file.path("..", "config.R"))
config_path <- candidate_config_paths[file.exists(candidate_config_paths)][1]
stopifnot(
  "Cannot locate config.R; run this script from the repository root (or tests/)." =
    !is.na(config_path)
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
candidate_plot_paths <- c("functions_plots.R", file.path("..", "functions_plots.R"))
plot_path <- candidate_plot_paths[file.exists(candidate_plot_paths)][1]
stopifnot(
    "Cannot locate functions_plots.R; run this script from the repository root (or tests/)." =
      !is.na(plot_path)
)

# ── 1. config.R: single canonical BH-FDR cutoff ──────────────────────────────
config_env <- new.env()
sys.source(config_path, envir = config_env)

stopifnot(
  "ALPHA_Q must be defined in config.R" = exists("ALPHA_Q", envir = config_env),
  "ALPHA_Q must equal 0.05" = identical(config_env$ALPHA_Q, 0.05),
  "ADJ_P_CUTOFF must be defined in config.R" =
    exists("ADJ_P_CUTOFF", envir = config_env),
  "ADJ_P_CUTOFF must equal ALPHA_Q (not a duplicated literal)" =
    identical(config_env$ADJ_P_CUTOFF, config_env$ALPHA_Q),
  "LOG2FC_CUTOFF must be defined in config.R" =
    exists("LOG2FC_CUTOFF", envir = config_env),
  "LOG2FC_CUTOFF must equal 1.0" = identical(config_env$LOG2FC_CUTOFF, 1.0),
  "FDR_THRESHOLD must be removed from config.R" =
    !exists("FDR_THRESHOLD", envir = config_env)
)

# The source text itself must show ADJ_P_CUTOFF referencing ALPHA_Q, not a
# re-typed `0.05` literal.
config_text <- paste(readLines(config_path, warn = FALSE), collapse = "\n")
stopifnot(
  "config.R must set ADJ_P_CUTOFF <- ALPHA_Q" =
    grepl("ADJ_P_CUTOFF\\s*<-\\s*ALPHA_Q", config_text),
  "config.R must not define FDR_THRESHOLD" =
    !grepl("FDR_THRESHOLD", config_text)
)

cat("PASS: config.R exposes a single canonical ALPHA_Q/ADJ_P_CUTOFF pair and",
    "LOG2FC_CUTOFF, with no conflicting FDR_THRESHOLD.\n")

# ── 2. Rmd: significant/deg_class rule combines FDR AND effect size ─────────
rmd_text <- paste(readLines(rmd_path, warn = FALSE), collapse = "\n")
plot_text <- paste(readLines(plot_path, warn = FALSE), collapse = "\n")

stopifnot(
  "Rmd must derive RNA_ADJ_P_CUTOFF from the global ADJ_P_CUTOFF" =
    grepl("RNA_ADJ_P_CUTOFF\\s*<-\\s*ADJ_P_CUTOFF", rmd_text),
  "Rmd must derive RNA_LOG2FC_CUTOFF from the global LOG2FC_CUTOFF" =
    grepl("RNA_LOG2FC_CUTOFF\\s*<-\\s*LOG2FC_CUTOFF", rmd_text),
  "Rmd must define an FDR-only fdr_significant field" =
    grepl("fdr_significant\\s*=", rmd_text),
  "Rmd `significant` must be derived from fdr_significant" =
    grepl("significant\\s*=\\s*fdr_significant", rmd_text),
  "Rmd `significant` must also require the log2FC cutoff (not FDR alone)" =
    grepl("significant\\s*=\\s*fdr_significant[^\\n]*RNA_LOG2FC_CUTOFF", rmd_text),
  "Shared volcano helper must draw the horizontal BH-FDR cutoff line" =
    grepl("yintercept\\s*=\\s*-log10\\(ADJ_P_CUTOFF\\)", plot_text),
  "Shared volcano helper must use fixed facet y scales" =
    grepl("facet_wrap\\([^\\n]*scales\\s*=\\s*\"fixed\"", plot_text),
  "Shared volcano y-axis range must be configured globally" =
    grepl("VOLCANO_Y_LIMITS\\s*<-\\s*c\\(0,\\s*50\\)", config_text)
)

subtitle_start <- regexpr("subtitle\\s*=\\s*sprintf\\(", rmd_text)
stopifnot("Volcano subtitle must be built with sprintf()" = subtitle_start > 0)
subtitle_block <- substr(rmd_text, subtitle_start, subtitle_start + 400L)
subtitle_end <- regexpr("\\)\\s*,", subtitle_block)
stopifnot("Could not locate end of volcano subtitle sprintf() call" = subtitle_end > 0)
subtitle_block <- substr(subtitle_block, 1, subtitle_end)

stopifnot(
  "Volcano subtitle must mention FDR" = grepl("FDR", subtitle_block),
  "Volcano subtitle must reference RNA_ADJ_P_CUTOFF" =
    grepl("RNA_ADJ_P_CUTOFF", subtitle_block),
  "Volcano subtitle must reference RNA_LOG2FC_CUTOFF" =
    grepl("RNA_LOG2FC_CUTOFF", subtitle_block)
)

cat("PASS: Rmd wires RNA_ADJ_P_CUTOFF/RNA_LOG2FC_CUTOFF into `fdr_significant`",
    "and the combined `significant` rule; the shared volcano helper applies the",
    "global FDR threshold and fixed y-axis range.\n")

# ── 3. Reproduce the mutate() rule on synthetic data with base R ────────────
# Mirrors (without requiring dplyr) the logic at the `dream-models` chunk:
#   fdr_significant = !is.na(q.value) & q.value < RNA_ADJ_P_CUTOFF
#   significant      = fdr_significant & abs(logFC) >= RNA_LOG2FC_CUTOFF
#   deg_class        = Positive / Negative / Not significant
rna_adj_p_cutoff <- config_env$ADJ_P_CUTOFF
rna_log2fc_cutoff <- config_env$LOG2FC_CUTOFF

toy <- data.frame(
  GENEID = c("A", "B", "C", "D", "E"),
  q.value = c(0.01, 0.01, 0.01, NA, 0.04),
  logFC = c(2.0, 0.3, -3.0, 1.5, -1.0)
)

fdr_significant <- !is.na(toy$q.value) & toy$q.value < rna_adj_p_cutoff
significant <- fdr_significant & abs(toy$logFC) >= rna_log2fc_cutoff
deg_class <- ifelse(
  significant & toy$logFC >= rna_log2fc_cutoff, "Positive",
  ifelse(significant & toy$logFC <= -rna_log2fc_cutoff, "Negative", "Not significant")
)

# Gene A: FDR-significant AND |logFC| >= 1 -> significant, Positive.
stopifnot(identical(fdr_significant[1], TRUE), identical(significant[1], TRUE))
stopifnot(identical(deg_class[1], "Positive"))

# Gene B: FDR-significant but |logFC| < 1 -> fdr_significant TRUE,
# significant FALSE (this is exactly the disagreement the fix removes).
stopifnot(identical(fdr_significant[2], TRUE), identical(significant[2], FALSE))
stopifnot(identical(deg_class[2], "Not significant"))

# Gene C: FDR-significant and large negative logFC -> Negative.
stopifnot(identical(significant[3], TRUE), identical(deg_class[3], "Negative"))

# Gene D: missing q-value must never be treated as significant.
stopifnot(identical(fdr_significant[4], FALSE), identical(significant[4], FALSE))

# Gene E: FDR-significant with logFC exactly at the cutoff -> significant.
stopifnot(identical(fdr_significant[5], TRUE), identical(significant[5], TRUE))
stopifnot(identical(deg_class[5], "Negative"))

# `significant` must never be TRUE unless `deg_class` is a real direction,
# and `deg_class == "Not significant"` must never coincide with
# `significant == TRUE`.
stopifnot(all(significant == (deg_class != "Not significant")))
# `significant` must always imply `fdr_significant` (effect-size filtering
# only removes genes, never adds them).
stopifnot(all(!significant | fdr_significant))

cat("PASS: `fdr_significant` (FDR-only) and `significant` (FDR + effect size)",
    "agree with `deg_class` on synthetic data, including the",
    "FDR-significant-but-small-effect disagreement case this fix resolves.\n")

cat("All RNA-seq significance threshold consistency checks passed.\n")
