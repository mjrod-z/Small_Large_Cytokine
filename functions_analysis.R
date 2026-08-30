# =============================================================================
# functions_analysis.R
# Statistical models: LMER, ART, screening, significance tables
# =============================================================================

suppressPackageStartupMessages({
  library(lme4)
  library(emmeans)
  library(ARTool)
})

# ── Hormone label recoding ────────────────────────────────────────────────────

recode_hormone_labels <- function(df, exp_short, e2_label, hormone_levels) {
  df %>%
    dplyr::mutate(
      HORMONE = dplyr::case_match(
        as.character(HORMONE),
        "NONE" ~ exp_short,
        "Estradiol" ~ e2_label,
        .default = as.character(HORMONE)
      ),
      HORMONE = factor(HORMONE, levels = hormone_levels)
    )
}

# ── log2FC long format ────────────────────────────────────────────────────────

make_log2fc_long <- function(df, cytokine_cols,
                             pbs_level         = PBS_LEVEL,
                             exclude_exposures = EXCLUDE_EXPOSURES,
                             pseudocount       = PSEUDOCOUNT) {
  d <- df %>%
    dplyr::filter(!EXPOSURE %in% exclude_exposures) %>%
    dplyr::mutate(
      SEX         = as.character(SEX),
      CELLTYPE    = factor(CELLTYPE),
      HORMONE     = factor(HORMONE),
      EXPOSURE    = factor(EXPOSURE),
      PATIENTCODE = factor(PATIENTCODE)
    )
  
  long <- d %>%
    tidyr::pivot_longer(cols = dplyr::all_of(cytokine_cols),
                        names_to = "CYTOKINE", values_to = "VALUE") %>%
    dplyr::filter(!is.na(VALUE)) %>%
    dplyr::mutate(log2_val = log2(as.numeric(VALUE) + pseudocount))
  
  pbs <- long %>%
    dplyr::filter(EXPOSURE == pbs_level) %>%
    dplyr::select(PATIENTCODE, SEX, CELLTYPE, HORMONE, CYTOKINE,
                  log2_pbs = log2_val)
  
  long %>%
    dplyr::filter(EXPOSURE != pbs_level) %>%
    dplyr::left_join(pbs,
                     by = c("PATIENTCODE","SEX","CELLTYPE","HORMONE","CYTOKINE")) %>%
    dplyr::mutate(log2FC = log2_val - log2_pbs)
}

# ── LMER screening: one exposure vs PBS, per CELLTYPE × HORMONE × SEX ────────
# NOTE: screen_one_exposure_lmer_log2() is SUPERSEDED and no longer called in
# 01_sala_analysis.Rmd.  The single source of statistical truth is now the
# exposure_lmer_pairwise() pipeline run inside run_lmer_chunk(), whose output
# is returned as lmer_plot_data and used for both significance tables AND plots.
# This function is kept here for reference only.

screen_one_exposure_lmer_log2 <- function(df, cytokine_cols, target_exposure,
                                          ctrl_level        = "PBS_Control",
                                          exclude_exposures = character(0),
                                          pseudocount       = 1e-6,
                                          group             = c("All","F","M"),
                                          alpha             = 0.1,
                                          emmeans_weights   = c("equal","proportional")) {
  group           <- match.arg(group)
  emmeans_weights <- match.arg(emmeans_weights)
  
  d0 <- df %>%
    dplyr::filter(!EXPOSURE %in% exclude_exposures,
                  EXPOSURE %in% c(ctrl_level, target_exposure)) %>%
    dplyr::mutate(
      EXPOSURE    = factor(EXPOSURE),
      PATIENTCODE = factor(PATIENTCODE),
      CELLTYPE    = factor(CELLTYPE),
      HORMONE     = factor(HORMONE)
    )
  
  if (group %in% c("F","M"))
    d0 <- d0 %>% dplyr::filter(as.character(SEX) == group)
  
  if (ctrl_level %in% levels(d0$EXPOSURE))
    d0$EXPOSURE <- relevel(d0$EXPOSURE, ref = ctrl_level)
  
  # Empty result template, used whenever nothing survives screening
  # (e.g. too few samples in this SEX/CELLTYPE/HORMONE stratum). Returning
  # this instead of an empty bind_rows() output avoids a downstream
  # `dplyr::group_by()` error on a 0-row/0-column tibble.
  empty_result <- tibble::tibble(
    SEX      = character(),
    CELLTYPE = character(),
    HORMONE  = character(),
    EXPOSURE = character(),
    CYTOKINE = character(),
    estimate = numeric(),
    SE       = numeric(),
    p.value  = numeric(),
    q        = numeric(),
    sig      = logical()
  )
  
  get_contrast <- function(fit) {
    emm  <- emmeans::emmeans(fit, ~ EXPOSURE, weights = emmeans_weights)
    levs <- levels(emmeans::summary(emm)$EXPOSURE)
    v    <- if (identical(levs, c(ctrl_level, target_exposure))) c(-1, 1) else c(1, -1)
    contrast_list        <- list(v)
    names(contrast_list) <- paste0(target_exposure, " - ", ctrl_level)
    emmeans::contrast(emm, method = contrast_list, adjust = "none")
  }
  
  combos        <- d0 %>% dplyr::distinct(CELLTYPE, HORMONE)
  cytokine_cols <- intersect(cytokine_cols, names(d0))
  out           <- list()
  
  for (i in seq_len(nrow(combos))) {
    ct   <- combos$CELLTYPE[i]
    ho   <- combos$HORMONE[i]
    dsub <- d0 %>% dplyr::filter(CELLTYPE == ct, HORMONE == ho)
    if (!all(c(ctrl_level, target_exposure) %in% unique(dsub$EXPOSURE))) next
    
    for (cyt in cytokine_cols) {
      dat <- dsub %>% dplyr::filter(!is.na(.data[[cyt]]))
      if (nrow(dat) < 3) next
      if (!all(c(ctrl_level, target_exposure) %in% unique(dat$EXPOSURE))) next
      
      dat  <- dat %>% dplyr::mutate(resp = log2(.data[[cyt]] + pseudocount))
      form <- resp ~ EXPOSURE + (1 | PATIENTCODE)
      
      fit <- try(lme4::lmer(form, data = dat), silent = TRUE)
      if (inherits(fit, "try-error")) next
      # Check for convergence / singular-fit warnings stored in the fit object
      if (length(lme4::isSingular(fit)) > 0 && lme4::isSingular(fit)) {
        warning("Singular fit for ", cyt, " in CELLTYPE=", ct, " HORMONE=", ho,
                "; estimates may be unreliable.")
      }
      
      con <- try(get_contrast(fit), silent = TRUE)
      if (inherits(con, "try-error")) next
      
      s <- as.data.frame(summary(con))
      out[[paste(group, ct, ho, cyt, sep = "|")]] <- tibble::tibble(
        SEX      = ifelse(group == "All", "All", group),
        CELLTYPE = as.character(ct),
        HORMONE  = as.character(ho),
        EXPOSURE = target_exposure,
        CYTOKINE = cyt,
        estimate = s$estimate[1],
        SE       = s$SE[1],
        p.value  = s$p.value[1]
      )
    }
  }
  
  result <- dplyr::bind_rows(out)
  
  if (nrow(result) == 0) {
    warning("screen_one_exposure_lmer_log2(): no cytokines survived screening for ",
            "group='", group, "', target_exposure='", target_exposure,
            "'. Returning an empty result.")
    return(empty_result)
  }
  
  result %>%
    dplyr::group_by(SEX, CELLTYPE, HORMONE, EXPOSURE) %>%
    dplyr::mutate(q   = p.adjust(p.value, method = "fdr"),
                  sig = q < alpha) %>%
    dplyr::ungroup()
}

# ── LMER pairwise: exposure vs PBS ────────────────────────────────────────────

exposure_lmer_pairwise <- function(data, group = "All", adjust_method = "fdr",
                                   response_columns = NULL,
                                   ctrl_level = "PBS_Control") {
  if (group != "All" && "SEX" %in% names(data))
    data <- data %>% dplyr::filter(SEX == group)
  
  stopifnot("EXPOSURE" %in% names(data), "PATIENTCODE" %in% names(data))
  
  data <- data %>%
    dplyr::mutate(
      EXPOSURE    = factor(EXPOSURE),
      PATIENTCODE = factor(PATIENTCODE)
    )
  
  if (ctrl_level %in% levels(data$EXPOSURE))
    data$EXPOSURE <- relevel(data$EXPOSURE, ref = ctrl_level)
  
  if (is.null(response_columns))
    response_columns <- names(data)[sapply(data, is.numeric)]
  
  results_list <- list()
  
  for (resp in response_columns) {
    if (!resp %in% names(data)) next
    df <- data %>% dplyr::filter(!is.na(.data[[resp]]))
    if (length(unique(df$EXPOSURE)) < 2) next
    
    model <- try(
      lme4::lmer(as.formula(paste(resp, "~ EXPOSURE + (1|PATIENTCODE)")),
                 data = df),
      silent = TRUE)
    if (inherits(model, "try-error")) { warning("Model failed for ", resp); next }
    
    emm      <- emmeans::emmeans(model, ~ EXPOSURE, weights = "equal")
    ctrl_idx <- which(levels(data$EXPOSURE) == ctrl_level)
    pairwise <- emmeans::contrast(emm, method = "trt.vs.ctrl",
                                  ref = ctrl_idx, adjust = adjust_method)
    pairwise_df          <- as.data.frame(summary(pairwise))
    pairwise_df$response <- resp
    results_list[[resp]] <- pairwise_df
  }
  
  dplyr::bind_rows(results_list)
}

# ── LMER pairwise: sex × exposure interaction ─────────────────────────────────

interaction_lmer_pairwise <- function(data, group = "All", adjust_method = "fdr",
                                      response_columns = NULL,
                                      ctrl_level = "PBS_Control") {
  if (group != "All" && "SEX" %in% names(data))
    data <- data %>% dplyr::filter(SEX == group)
  
  stopifnot("EXPOSURE" %in% names(data), "PATIENTCODE" %in% names(data),
            "SEX"      %in% names(data))
  
  data <- data %>%
    dplyr::mutate(
      EXPOSURE    = factor(EXPOSURE),
      PATIENTCODE = factor(PATIENTCODE),
      SEX         = factor(SEX)
    )
  
  if (ctrl_level %in% levels(data$EXPOSURE))
    data$EXPOSURE <- relevel(data$EXPOSURE, ref = ctrl_level)
  
  if (is.null(response_columns))
    response_columns <- names(data)[sapply(data, is.numeric)]
  
  results_list <- list()
  
  for (resp in response_columns) {
    if (!resp %in% names(data)) next
    df <- data %>% dplyr::filter(!is.na(.data[[resp]]))
    if (length(unique(df$EXPOSURE)) < 2 || length(unique(df$SEX)) < 2) next
    
    model <- try(
      lme4::lmer(
        as.formula(paste(resp, "~ EXPOSURE * SEX + (1|PATIENTCODE)")),
        data = df),
      silent = TRUE)
    if (inherits(model, "try-error")) {
      warning("Interaction model failed for ", resp); next
    }
    
    emm_exp  <- emmeans::emmeans(model, ~ EXPOSURE, weights = "equal")
    ctrl_idx <- which(levels(data$EXPOSURE) == ctrl_level)
    pw_exp   <- emmeans::contrast(emm_exp, "trt.vs.ctrl",
                                  ref = ctrl_idx, adjust = adjust_method)
    exp_df   <- as.data.frame(summary(pw_exp))
    exp_df$type     <- "Exposure_vs_Control"
    exp_df$response <- resp
    
    emm_int  <- emmeans::emmeans(model, ~ SEX | EXPOSURE)
    pw_int   <- emmeans::contrast(emm_int, "pairwise",
                                  simple = "SEX", combine = TRUE)
    int_df   <- as.data.frame(summary(pw_int))
    int_df$type     <- "Sex_within_Exposure"
    int_df$response <- resp
    
    results_list[[resp]] <- dplyr::bind_rows(exp_df, int_df)
  }
  
  dplyr::bind_rows(results_list)
}

# ── ART pairwise ──────────────────────────────────────────────────────────────

exposure_art_pairwise <- function(data, group = "All", adjust_method = "fdr",
                                  response_columns = NULL,
                                  ctrl_level = "PBS_Control") {
  filtered_data <- switch(group,
                          "M"   = data %>% dplyr::filter(EXPOSURE != "Untreated_Control", SEX == "M"),
                          "F"   = data %>% dplyr::filter(EXPOSURE != "Untreated_Control", SEX == "F"),
                          "All" = data %>% dplyr::filter(EXPOSURE != "Untreated_Control"),
                          stop("Invalid group. Choose 'All', 'M', or 'F'.")
  ) %>%
    dplyr::mutate(
      SEX         = factor(SEX),
      EXPOSURE    = factor(EXPOSURE),
      PATIENTCODE = factor(PATIENTCODE)
    )
  
  if (is.null(response_columns))
    response_columns <- names(filtered_data)[sapply(filtered_data, is.numeric)]
  
  results_list <- list()
  
  for (response in response_columns) {
    if (!response %in% names(filtered_data)) {
      warning("Column ", response, " not found"); next
    }
    
    formula <- if (group == "All") {
      as.formula(paste(response, "~ EXPOSURE * SEX + (1|PATIENTCODE)"))
    } else {
      as.formula(paste(response, "~ EXPOSURE + (1|PATIENTCODE)"))
    }
    
    m.art <- try(ARTool::art(formula, data = filtered_data), silent = TRUE)
    if (inherits(m.art, "try-error")) {
      warning("ART failed for ", response); next
    }
    
    anova_res   <- anova(m.art)
    exposure_p  <- {
      r <- anova_res[grepl("^EXPOSURE$", anova_res[[1]], ignore.case = TRUE), ]
      if (nrow(r) > 0) r[["Pr(>F)"]][1] else NA
    }
    interaction_p <- if (group == "All") {
      r <- anova_res[grepl("EXPOSURE:SEX", anova_res[[1]], ignore.case = TRUE), ]
      if (nrow(r) > 0) r[["Pr(>F)"]][1] else NA
    } else NA
    
    # art.con() is the correct way to get pairwise contrasts from an ART model
    # (artlm.con() + emmeans() produces unreliable d.f. for interaction contrasts).
    pairwise_con <- try(
      ARTool::art.con(m.art, "EXPOSURE", adjust = adjust_method),
      silent = TRUE
    )
    if (inherits(pairwise_con, "try-error")) {
      warning("art.con() failed for ", response); next
    }
    ctrl_idx <- which(levels(filtered_data$EXPOSURE) == ctrl_level)
    pairwise <- emmeans::contrast(pairwise_con, "trt.vs.ctrl",
                                  ref = ctrl_idx, adjust = adjust_method)
    
    res               <- as.data.frame(pairwise)
    res$response      <- response
    res$exposure_p    <- exposure_p
    res$interaction_p <- interaction_p
    results_list[[response]] <- res
  }
  
  dplyr::bind_rows(results_list)
}

# ── Convenience wrapper: run all 4 CELLTYPE × HORMONE strata ─────────────────

run_lmer_chunk <- function(label, celltype_filter, hormone_filter,
                           sala_full     = NULL,
                           llod_table    = cytokine_llod,
                           zero_co       = ZERO_CUTOFF,
                           alpha_q       = ALPHA_Q,
                           trend_a       = TREND_ALPHA,
                           all_cytokines = NULL) {
  
  cat("\n── LMER:", label, "──\n")
  
  if (is.null(sala_full))
    stop("sala_full must be provided as a data frame containing the full SALA dataset")
  
  d_sub <- sala_full %>%
    dplyr::filter(CELLTYPE == celltype_filter, HORMONE == hormone_filter)
  
  imp_res      <- impute_lod_sqrt2(d_sub, cytokine_llod = llod_table,
                                   zero_cutoff = zero_co)
  d_imp        <- imp_res$data
  valid_cyts   <- imp_res$valid_cytokines
  
  d_filt <- d_imp %>%
    dplyr::filter(EXPOSURE != "Untreated_Control") %>%
    dplyr::mutate(EXPOSURE    = factor(EXPOSURE),
                  SEX         = factor(SEX),
                  PATIENTCODE = factor(PATIENTCODE))
  
  msd_sum <- summarize_to_wide(d_filt, measure_vars = valid_cyts)
  pbs_ctl <- msd_sum %>%
    dplyr::select(Measurement, `PBS_Control`) %>%
    dplyr::arrange(Measurement)
  
  lmer_All <- exposure_lmer_pairwise(d_filt, "All", "fdr", valid_cyts)
  lmer_F   <- exposure_lmer_pairwise(d_filt, "F",   "fdr", valid_cyts)
  lmer_M   <- exposure_lmer_pairwise(d_filt, "M",   "fdr", valid_cyts)
  lmer_int <- interaction_lmer_pairwise(d_filt, "All", "fdr", valid_cyts)
  
  # Convert LMER results to plot-compatible format (used by cytokine dotplots
  # and bar plots instead of the retired screen_one_exposure_lmer_log2()).
  lmer_plot_data <- lmer_results_to_plot_format(
    lmer_All, lmer_F, lmer_M,
    celltype = celltype_filter,
    hormone  = hormone_filter,
    alpha    = alpha_q
  )
  
  out_csv <- paste0("MSD_SALA_", label, "_lmer.csv")
  sig_tbl <- build_lmer_sig_table(
    lmer_All      = lmer_All,
    lmer_F        = lmer_F,
    lmer_M        = lmer_M,
    msd_summary   = msd_sum,
    pbs_control   = pbs_ctl,
    cytokine_llod = llod_table,
    out_filename  = out_csv,
    alpha         = alpha_q,
    trend_alpha   = trend_a,
    all_cytokines = all_cytokines
  )
  
  cat("  Saved:", out_csv, "\n")
  # Return significance table, interaction results, and plot-ready LMER data
  invisible(list(sig_table = sig_tbl, lmer_interaction = lmer_int,
                 lmer_plot_data = lmer_plot_data))
}

build_lmer_sig_table <- function(lmer_All, lmer_F, lmer_M,
                                 msd_summary, pbs_control,
                                 cytokine_llod,
                                 out_filename,
                                 alpha         = 0.1,
                                 trend_alpha   = 0.2,
                                 epsilon       = 1e-5,
                                 all_cytokines = NULL) {
  
  # Step 6: format with group suffixes
  fmt <- function(df, grp) {
    col <- paste0("response_", grp)
    df %>%
      dplyr::mutate(
        !!col    := paste0(response, "_", grp),
        EXPOSURE  = stringr::str_replace(contrast, " - PBS_Control", "")
      )
  }
  combined <- dplyr::bind_rows(fmt(lmer_All, "All"),
                               fmt(lmer_F,   "F"),
                               fmt(lmer_M,   "M"))
  
  # Step 7: apply FDR correction per group, then assign direction-aware stars
  # Stars are based on the FDR-adjusted q-value, not the raw p.value.
    stars_df <- combined %>%
    dplyr::group_by(response) %>%
    dplyr::mutate(q = p.adjust(p.value, method = "fdr")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(stars = dplyr::case_when(
      q < 0.001 ~ "***",
      q < 0.01  ~ "**",
      q < 0.05  ~ "*",
      TRUE      ~ ""
    ))
  
  # Steps 8-9: pivot stars long
  stars_long <- stars_df %>%
    tidyr::pivot_longer(
      cols           = dplyr::any_of(c("response_All","response_F","response_M")),
      names_to       = "Group",
      values_to      = "Measurement",
      values_drop_na = TRUE
    ) %>%
    dplyr::select(Measurement, EXPOSURE, stars) %>%
    tidyr::pivot_wider(names_from = EXPOSURE, values_from = stars) %>%
    tidyr::pivot_longer(cols = -Measurement,
                        names_to = "Exposure", values_to = "Stars") %>%
    dplyr::mutate(Exposure = stringr::str_trim(
      gsub(" - PBS_Control$", "", Exposure)))
  
  # Step 10: summary without PBS
  wide_no_pbs <- msd_summary %>%
    dplyr::select(-`PBS_Control`) %>%
    dplyr::arrange(Measurement) %>%
    dplyr::mutate(Analyte_Base = gsub("_.*", "", Measurement)) %>%
    dplyr::left_join(cytokine_llod, by = c("Analyte_Base" = "Analyte"))
  
  if (!"LLOD" %in% colnames(wide_no_pbs))
    stop("LLOD column not found. Check cytokine_llod data.")
  
  # Step 11: pivot long + numeric values
  long_vals <- wide_no_pbs %>%
    tidyr::pivot_longer(
      cols      = -c(Measurement, Analyte_Base, LLOD),
      names_to  = "Exposure",
      values_to = "Value_full"
    ) %>%
    dplyr::mutate(
      Value_num = as.numeric(sub("^(\\d*\\.?\\d+).*", "\\1", Value_full))
    )
  
  # Step 12: merge stars + LLOD check
  long_vals <- long_vals %>%
    dplyr::left_join(stars_long, by = c("Measurement","Exposure")) %>%
    dplyr::mutate(Value_final = ifelse(
      !is.na(Stars) & !is.na(Value_num) & Value_num >= LLOD,
      paste0(Value_full, " ", Stars),
      Value_full
    ))
  
  # Step 13: pivot wide
  wide_final <- long_vals %>%
    dplyr::select(Measurement, Exposure, Value_final, LLOD, Analyte_Base) %>%
    tidyr::pivot_wider(names_from = Exposure, values_from = Value_final)
  
  # Step 14: add PBS, apply epsilon, export
  pbs_numeric <- as.numeric(sub(" \u00b1.*", "", pbs_control$`PBS_Control`))
  
  sig_table <- wide_final %>%
    dplyr::left_join(pbs_control, by = "Measurement") %>%
    dplyr::mutate(dplyr::across(
      dplyr::where(is.character) &
        !dplyr::all_of(c("Measurement","PBS_Control")),
      ~ ifelse(pbs_numeric < epsilon, gsub("[*#~]+", "", .), .)
    )) %>%
    dplyr::select(Measurement, `PBS_Control`, dplyr::everything(),
                  -Analyte_Base, -LLOD)

  # Append placeholder rows for raw-panel cytokines not reached by the model
  # (absent from cytokine_llod or dropped by zero-cutoff filter). These rows
  # will have NA for all exposure columns and PBS_Control so they are clearly
  # identifiable as untested in the exported CSV.
  if (!is.null(all_cytokines)) {
    tested_cyts <- gsub("_(All|F|M)$", "", sig_table$Measurement) %>% unique()
    missing_cyts <- setdiff(all_cytokines, tested_cyts)
    if (length(missing_cyts) > 0) {
      # Build one "_All" stub row per missing cytokine with NA exposure values
      stub_rows <- tibble::tibble(
        Measurement = paste0(missing_cyts, "_All"),
        PBS_Control = NA_character_
      )
      sig_table <- dplyr::bind_rows(sig_table, stub_rows)
    }
  }

  save_table(sig_table, out_filename)
  invisible(sig_table)
}

# ── Convert exposure_lmer_pairwise() output to plot-compatible format ─────────
# Used internally by run_lmer_chunk() to produce lmer_plot_data.
# FDR correction matches build_lmer_sig_table() exactly:
#   group_by(CYTOKINE) %>% p.adjust(p.value, method = "fdr")
# so q-values are IDENTICAL between significance tables and plots.

lmer_results_to_plot_format <- function(lmer_All, lmer_F, lmer_M,
                                        celltype, hormone,
                                        alpha = ALPHA_Q) {
  parse_rows <- function(df, sex_label) {
    df %>%
      dplyr::mutate(
        SEX      = sex_label,
        CELLTYPE = celltype,
        HORMONE  = hormone,
        EXPOSURE = stringr::str_replace(contrast, " - PBS_Control$", ""),
        CYTOKINE = response
      ) %>%
      dplyr::select(SEX, CELLTYPE, HORMONE, EXPOSURE, CYTOKINE,
                    estimate, SE, p.value)
  }
  
  combined <- dplyr::bind_rows(
    parse_rows(lmer_All, "All"),
    parse_rows(lmer_F,   "F"),
    parse_rows(lmer_M,   "M")
  )
  
  # FDR correction matching build_lmer_sig_table(): group by CYTOKINE
  # across all sex groups and exposures together.
  combined %>%
    dplyr::group_by(CYTOKINE) %>%
    dplyr::mutate(q = p.adjust(p.value, method = "fdr")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sig = q < alpha)
}

# ── LLOD imputation + summary helpers ─────────────────────────────────────────

impute_lod_sqrt2 <- function(input_data, cols = NULL, cytokine_llod, zero_cutoff = ZERO_CUTOFF) {
  if (is.null(cols)) {
    cols <- intersect(names(input_data), cytokine_llod$Analyte)
  }
  
  llod_map <- stats::setNames(cytokine_llod$LLOD, cytokine_llod$Analyte)
  valid_cytokines <- cols[vapply(cols, function(col) {
    vals <- suppressWarnings(as.numeric(input_data[[col]]))
    if (all(is.na(vals))) return(FALSE)
    
    if ("SEX" %in% names(input_data)) {
      max_zero <- input_data %>%
        dplyr::mutate(.val = vals) %>%
        dplyr::group_by(SEX) %>%
        dplyr::summarise(p_zero = mean(is.na(.val) | .val <= 0), .groups = "drop") %>%
        dplyr::summarise(max_p = max(p_zero, na.rm = TRUE)) %>%
        dplyr::pull(max_p)
      return(is.finite(max_zero) && max_zero <= zero_cutoff)
    }
    
    mean(is.na(vals) | vals <= 0) <= zero_cutoff
  }, logical(1))]
  
  skipped_cytokines <- setdiff(cols, valid_cytokines)
  
  out <- input_data
  for (col in valid_cytokines) {
    llod_val <- llod_map[[col]]
    if (!is.null(llod_val) && is.finite(llod_val)) {
      vals <- suppressWarnings(as.numeric(out[[col]]))
      vals[!is.na(vals) & vals < llod_val] <- llod_val / sqrt(2)
      out[[col]] <- vals
    }
  }
  
  list(data = out, valid_cytokines = valid_cytokines, skipped_cytokines = skipped_cytokines)
}

summarize_to_wide <- function(data, measure_vars) {
  long <- data %>%
    dplyr::select(dplyr::any_of(c("EXPOSURE", "SEX", measure_vars))) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(measure_vars),
                        names_to = "CYTOKINE", values_to = "VALUE")
  
  by_sex <- long %>%
    dplyr::group_by(CYTOKINE, SEX, EXPOSURE) %>%
    dplyr::summarise(
      mu = mean(VALUE, na.rm = TRUE),
      sd = stats::sd(VALUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Measurement = paste0(CYTOKINE, "_", SEX))
  
  all_rows <- long %>%
    dplyr::group_by(CYTOKINE, EXPOSURE) %>%
    dplyr::summarise(
      mu = mean(VALUE, na.rm = TRUE),
      sd = stats::sd(VALUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Measurement = paste0(CYTOKINE, "_All"))
  
  dplyr::bind_rows(by_sex, all_rows) %>%
    dplyr::mutate(
      mu = ifelse(is.nan(mu), NA_real_, mu),
      sd = ifelse(is.na(sd) | is.nan(sd), 0, sd),
      Value = sprintf("%.3f \u00b1 %.3f", mu, sd)
    ) %>%
    dplyr::select(Measurement, EXPOSURE, Value) %>%
    tidyr::pivot_wider(names_from = EXPOSURE, values_from = Value) %>%
    dplyr::arrange(Measurement)
}
