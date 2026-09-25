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
    levs <- levels(emm)[["EXPOSURE"]]
    ctrl_idx <- match(ctrl_level, levs)
    target_idx <- match(target_exposure, levs)
    if (is.na(ctrl_idx) || is.na(target_idx)) {
      stop(
        "Expected contrast levels not found in emmeans results for EXPOSURE. ",
        "Requested control='", ctrl_level,
        "', target='", target_exposure,
        "'. Available levels: ",
        paste(levs, collapse = ", ")
      )
    }
    v <- rep(0, length(levs))
    v[ctrl_idx] <- -1
    v[target_idx] <- 1
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


compute_hpiv3_missingness_qc <- function(data, protein_cols,
                                         threshold = ZERO_CUTOFF,
                                         strata_cols = c("AIRWAY", "HORMONE", "TIMEPOINT"),
                                         sex_col = "SEX") {
  stopifnot(is.data.frame(data))
  protein_cols <- intersect(protein_cols, names(data))

  build_scope <- function(group_cols, scope_name) {
    long <- data %>%
      dplyr::select(dplyr::any_of(c(group_cols, "INFECTION", protein_cols))) %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(protein_cols),
        names_to = "PROTEIN",
        values_to = "VALUE"
      ) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c("PROTEIN", group_cols)))) %>%
      dplyr::summarise(
        n_samples = dplyr::n(),
        n_observed = sum(!is.na(VALUE)),
        n_missing = sum(is.na(VALUE)),
        pct_missing = n_missing / n_samples,
        n_none = sum(INFECTION == "NONE" & !is.na(VALUE), na.rm = TRUE),
        n_hpiv3 = sum(INFECTION == "HPIV3" & !is.na(VALUE), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(scope = scope_name)

    if (length(group_cols) == 0) {
      long <- long %>% dplyr::mutate(scope_label = "overall")
    } else {
      long <- long %>%
        tidyr::unite("scope_label", dplyr::all_of(group_cols), sep = " | ", remove = FALSE)
    }

    long
  }

  overall_qc <- build_scope(character(0), "overall")
  stratum_qc <- build_scope(strata_cols, "analysis_stratum")
  sex_qc <- if (sex_col %in% names(data)) {
    build_scope(c(strata_cols, sex_col), "analysis_stratum_by_sex")
  } else {
    tibble::tibble()
  }

  protein_status <- stratum_qc %>%
    dplyr::group_by(PROTEIN) %>%
    dplyr::summarise(
      max_pct_missing = max(pct_missing, na.rm = TRUE),
      any_exceeds_threshold = any(pct_missing > threshold),
      all_missing = all(n_observed == 0),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      included = !all_missing & !any_exceeds_threshold,
      exclusion_reason = dplyr::case_when(
        all_missing ~ "all values missing",
        any_exceeds_threshold ~ paste0(
          "missingness exceeds ", formatC(threshold * 100, format = "f", digits = 0),
          "% in at least one AIRWAY/HORMONE/TIMEPOINT stratum"
        ),
        TRUE ~ NA_character_
      )
    )

  qc <- dplyr::bind_rows(overall_qc, stratum_qc, sex_qc) %>%
    dplyr::left_join(
      protein_status %>% dplyr::select(PROTEIN, included, exclusion_reason, max_pct_missing),
      by = "PROTEIN"
    ) %>%
    dplyr::mutate(missingness_threshold = threshold)

  list(qc = qc, protein_status = protein_status)
}

fit_hpiv3_infection_models <- function(data, protein_cols,
                                       protein_status = NULL,
                                       pseudocount = PSEUDOCOUNT,
                                       alpha_q = ALPHA_Q,
                                       min_nonmissing_per_group = 2L,
                                       min_total_nonmissing = 3L,
                                       control_level = "NONE",
                                       case_level = "HPIV3") {
  stopifnot(is.data.frame(data))
  protein_cols <- intersect(protein_cols, names(data))
  strata <- data %>% dplyr::distinct(AIRWAY, HORMONE, TIMEPOINT, EXPOSURE)
  sex_groups <- c("All", intersect(c("F", "M"), unique(as.character(data$SEX))))

  if (is.null(protein_status)) {
    protein_status <- tibble::tibble(
      PROTEIN = protein_cols,
      included = TRUE,
      exclusion_reason = NA_character_
    )
  }
  protein_status <- protein_status %>%
    dplyr::select(PROTEIN, included, exclusion_reason)

  get_contrast <- function(fit) {
    emm <- emmeans::emmeans(fit, ~ INFECTION, weights = "equal")
    infection_levels <- levels(emm)[["INFECTION"]]
    ctrl_idx <- match(control_level, infection_levels)
    if (is.na(ctrl_idx)) {
      stop(
        "Control level '", control_level,
        "' not found in emmeans results for INFECTION. Available levels: ",
        paste(infection_levels, collapse = ", ")
      )
    }
    emmeans::contrast(emm, method = "trt.vs.ctrl", ref = ctrl_idx, adjust = "none")
  }

  results <- list()
  idx <- 1L

  for (row_idx in seq_len(nrow(strata))) {
    airway_i <- strata$AIRWAY[[row_idx]]
    hormone_i <- strata$HORMONE[[row_idx]]
    timepoint_i <- strata$TIMEPOINT[[row_idx]]
    exposure_i <- strata$EXPOSURE[[row_idx]]

    stratum_data <- data %>%
      dplyr::filter(
        AIRWAY == airway_i,
        HORMONE == hormone_i,
        TIMEPOINT == timepoint_i,
        EXPOSURE == exposure_i
      )

    for (sex_group in sex_groups) {
      sex_data <- if (sex_group == "All") {
        stratum_data
      } else {
        stratum_data %>% dplyr::filter(as.character(SEX) == sex_group)
      }

      for (protein in protein_cols) {
        protein_rule <- protein_status %>% dplyr::filter(PROTEIN == protein)
        included <- if (nrow(protein_rule) == 1) isTRUE(protein_rule$included[[1]]) else TRUE
        exclusion_reason <- if (nrow(protein_rule) == 1) protein_rule$exclusion_reason[[1]] else NA_character_

        dat <- sex_data %>%
          dplyr::transmute(
            PATIENTCODE = PATIENTCODE,
            INFECTION = factor(INFECTION, levels = c(control_level, case_level)),
            VALUE = .data[[protein]]
          ) %>%
          dplyr::filter(!is.na(PATIENTCODE), !is.na(INFECTION), !is.na(VALUE))

        n_obs <- nrow(dat)
        n_none <- sum(dat$INFECTION == control_level)
        n_hpiv3 <- sum(dat$INFECTION == case_level)
        n_donors <- dplyr::n_distinct(dat$PATIENTCODE)
        repeated_donor <- anyDuplicated(as.character(dat$PATIENTCODE)) > 0

        result_row <- tibble::tibble(
          AIRWAY = as.character(airway_i),
          HORMONE = as.character(hormone_i),
          TIMEPOINT = as.character(timepoint_i),
          EXPOSURE = as.character(exposure_i),
          SEX = sex_group,
          PROTEIN = protein,
          contrast = paste(case_level, "-", control_level),
          n_samples = n_obs,
          n_none = n_none,
          n_hpiv3 = n_hpiv3,
          n_donors = n_donors,
          used_random_intercept = FALSE,
          model_type = NA_character_,
          model_status = "skipped",
          failure_reason = NA_character_,
          estimate = NA_real_,
          SE = NA_real_,
          p.value = NA_real_,
          singular_fit = NA
        )

        if (!included) {
          result_row$failure_reason <- exclusion_reason
          results[[idx]] <- result_row
          idx <- idx + 1L
          next
        }
        if (n_obs < min_total_nonmissing) {
          result_row$failure_reason <- paste0("fewer than ", min_total_nonmissing, " non-missing observations")
          results[[idx]] <- result_row
          idx <- idx + 1L
          next
        }
        if (n_none < min_nonmissing_per_group || n_hpiv3 < min_nonmissing_per_group) {
          result_row$failure_reason <- paste0(
            "insufficient per-infection observations (NONE=", n_none,
            ", HPIV3=", n_hpiv3, ")"
          )
          results[[idx]] <- result_row
          idx <- idx + 1L
          next
        }
        if (length(unique(stats::na.omit(dat$INFECTION))) < 2) {
          result_row$failure_reason <- "both infection levels not present"
          results[[idx]] <- result_row
          idx <- idx + 1L
          next
        }

        dat <- dat %>% dplyr::mutate(log2_value = log2(as.numeric(VALUE) + pseudocount))

        fit <- NULL
        used_lmer <- FALSE
        singular_fit <- FALSE
        lmer_error <- NA_character_
        lm_error <- NA_character_

        if (n_donors >= 2 && repeated_donor) {
          fit <- tryCatch(
            lme4::lmer(log2_value ~ INFECTION + (1 | PATIENTCODE), data = dat),
            error = function(e) {
              lmer_error <<- conditionMessage(e)
              NULL
            }
          )
          if (!is.null(fit)) {
            singular_fit <- isTRUE(tryCatch(lme4::isSingular(fit), error = function(...) FALSE))
            if (!singular_fit) {
              used_lmer <- TRUE
            }
          }
        }

        if (!used_lmer) {
          fit <- tryCatch(
            stats::lm(log2_value ~ INFECTION, data = dat),
            error = function(e) {
              lm_error <<- conditionMessage(e)
              NULL
            }
          )
          if (is.null(fit)) {
            base_reason <- if (singular_fit) {
              "singular mixed model and fallback linear model failed"
            } else {
              "model fitting failed"
            }
            detail_reason <- lm_error
            if (is.na(detail_reason) || !nzchar(detail_reason)) {
              detail_reason <- lmer_error
            }
            result_row$failure_reason <- if (is.na(detail_reason) || !nzchar(detail_reason)) {
              base_reason
            } else {
              paste0(base_reason, ": ", detail_reason)
            }
            result_row$singular_fit <- singular_fit
            results[[idx]] <- result_row
            idx <- idx + 1L
            next
          }
        }

        contrast_error <- NA_character_
        contrast <- tryCatch(
          get_contrast(fit),
          error = function(e) {
            contrast_error <<- conditionMessage(e)
            NULL
          }
        )
        if (is.null(contrast)) {
          result_row$failure_reason <- if (is.na(contrast_error) || !nzchar(contrast_error)) {
            "emmeans contrast failed"
          } else {
            paste0("emmeans contrast failed: ", contrast_error)
          }
          result_row$singular_fit <- singular_fit
          results[[idx]] <- result_row
          idx <- idx + 1L
          next
        }

        stats_row <- as.data.frame(summary(contrast))
        result_row$used_random_intercept <- used_lmer
        result_row$model_type <- if (used_lmer) "lmer" else "lm"
        result_row$model_status <- if (used_lmer) "modeled" else "modeled_fallback"
        result_row$failure_reason <- if (!used_lmer && singular_fit) {
          "singular mixed model; used linear-model fallback"
        } else if (!used_lmer && !(n_donors >= 2 && repeated_donor)) {
          "random intercept unsupported; used linear-model fallback"
        } else {
          NA_character_
        }
        result_row$estimate <- stats_row$estimate[[1]]
        result_row$SE <- stats_row$SE[[1]]
        result_row$p.value <- stats_row$p.value[[1]]
        result_row$singular_fit <- singular_fit

        results[[idx]] <- result_row
        idx <- idx + 1L
      }
    }
  }

  out <- dplyr::bind_rows(results) %>%
    dplyr::group_by(AIRWAY, HORMONE, TIMEPOINT, EXPOSURE, SEX) %>%
    dplyr::mutate(
      q.value = {
        q_vals <- rep(NA_real_, dplyr::n())
        keep <- which(!is.na(p.value))
        if (length(keep) > 0) {
          q_vals[keep] <- p.adjust(p.value[keep], method = "fdr")
        }
        q_vals
      },
      significant = !is.na(q.value) & q.value < alpha_q
    ) %>%
    dplyr::ungroup()

  out
}

fit_hpiv3_sex_models <- function(data, protein_cols,
                                 protein_status = NULL,
                                 pseudocount = PSEUDOCOUNT,
                                 alpha_q = ALPHA_Q,
                                 min_nonmissing_per_group = 2L,
                                 min_total_nonmissing = 3L) {
  stopifnot(is.data.frame(data))
  protein_cols <- intersect(protein_cols, names(data))
  strata <- data %>% dplyr::distinct(AIRWAY, HORMONE, TIMEPOINT, EXPOSURE, INFECTION)

  if (is.null(protein_status)) {
    protein_status <- tibble::tibble(
      PROTEIN = protein_cols,
      included = TRUE,
      exclusion_reason = NA_character_
    )
  }
  protein_status <- protein_status %>%
    dplyr::select(PROTEIN, included, exclusion_reason)

  results <- list()
  idx <- 1L

  for (row_idx in seq_len(nrow(strata))) {
    airway_i <- strata$AIRWAY[[row_idx]]
    hormone_i <- strata$HORMONE[[row_idx]]
    timepoint_i <- strata$TIMEPOINT[[row_idx]]
    exposure_i <- strata$EXPOSURE[[row_idx]]
    infection_i <- strata$INFECTION[[row_idx]]

    stratum_data <- data %>%
      dplyr::filter(
        AIRWAY == airway_i,
        HORMONE == hormone_i,
        TIMEPOINT == timepoint_i,
        EXPOSURE == exposure_i,
        INFECTION == infection_i
      )

    for (protein in protein_cols) {
      protein_rule <- protein_status %>% dplyr::filter(PROTEIN == protein)
      included <- if (nrow(protein_rule) == 1) isTRUE(protein_rule$included[[1]]) else TRUE
      exclusion_reason <- if (nrow(protein_rule) == 1) protein_rule$exclusion_reason[[1]] else NA_character_

      dat_raw <- stratum_data %>%
        dplyr::transmute(
          PATIENTCODE = PATIENTCODE,
          SEX = as.character(SEX),
          VALUE = .data[[protein]]
        ) %>%
        dplyr::filter(!is.na(PATIENTCODE), !is.na(SEX), !is.na(VALUE))

      sex_order <- c(intersect(c("F", "M"), unique(dat_raw$SEX)), sort(setdiff(unique(dat_raw$SEX), c("F", "M"))))
      dat <- dat_raw %>%
        dplyr::mutate(SEX = factor(SEX, levels = sex_order))

      sex_counts <- table(dat$SEX)
      n_obs <- nrow(dat)
      n_donors <- dplyr::n_distinct(dat$PATIENTCODE)
      repeated_donor <- anyDuplicated(as.character(dat$PATIENTCODE)) > 0
      counts_label <- if (length(sex_counts) == 0) {
        NA_character_
      } else {
        paste0(names(sex_counts), "=", as.integer(sex_counts), collapse = "; ")
      }

      result_row <- tibble::tibble(
        AIRWAY = as.character(airway_i),
        HORMONE = as.character(hormone_i),
        TIMEPOINT = as.character(timepoint_i),
        EXPOSURE = as.character(exposure_i),
        INFECTION = as.character(infection_i),
        PROTEIN = protein,
        comparison = "SEX",
        contrast = NA_character_,
        n_samples = n_obs,
        n_groups = length(sex_counts),
        group_counts = counts_label,
        n_donors = n_donors,
        used_random_intercept = FALSE,
        model_type = NA_character_,
        model_status = "skipped",
        failure_reason = NA_character_,
        estimate = NA_real_,
        SE = NA_real_,
        p.value = NA_real_,
        singular_fit = NA
      )

      if (!included) {
        result_row$failure_reason <- exclusion_reason
        results[[idx]] <- result_row
        idx <- idx + 1L
        next
      }
      if (n_obs < min_total_nonmissing) {
        result_row$failure_reason <- paste0("fewer than ", min_total_nonmissing, " non-missing observations")
        results[[idx]] <- result_row
        idx <- idx + 1L
        next
      }

      valid_levels <- names(sex_counts)[as.integer(sex_counts) >= min_nonmissing_per_group]
      dat_fit <- dat %>%
        dplyr::filter(as.character(SEX) %in% valid_levels) %>%
        droplevels()

      if (length(unique(as.character(dat_fit$SEX))) < 2) {
        result_row$failure_reason <- paste0(
          "insufficient per-sex observations (required >= ", min_nonmissing_per_group,
          " per level; observed: ", counts_label, ")"
        )
        results[[idx]] <- result_row
        idx <- idx + 1L
        next
      }

      dat_fit <- dat_fit %>% dplyr::mutate(log2_value = log2(as.numeric(VALUE) + pseudocount))

      fit <- NULL
      used_lmer <- FALSE
      singular_fit <- FALSE
      lmer_error <- NA_character_
      lm_error <- NA_character_

      if (n_donors >= 2 && repeated_donor) {
        fit <- tryCatch(
          lme4::lmer(log2_value ~ SEX + (1 | PATIENTCODE), data = dat_fit),
          error = function(e) {
            lmer_error <<- conditionMessage(e)
            NULL
          }
        )
        if (!is.null(fit)) {
          singular_fit <- isTRUE(tryCatch(lme4::isSingular(fit), error = function(...) FALSE))
          if (!singular_fit) {
            used_lmer <- TRUE
          }
        }
      }

      if (!used_lmer) {
        fit <- tryCatch(
          stats::lm(log2_value ~ SEX, data = dat_fit),
          error = function(e) {
            lm_error <<- conditionMessage(e)
            NULL
          }
        )
        if (is.null(fit)) {
          base_reason <- if (singular_fit) {
            "singular mixed model and fallback linear model failed"
          } else {
            "model fitting failed"
          }
          detail_reason <- lm_error
          if (is.na(detail_reason) || !nzchar(detail_reason)) {
            detail_reason <- lmer_error
          }
          result_row$failure_reason <- if (is.na(detail_reason) || !nzchar(detail_reason)) {
            base_reason
          } else {
            paste0(base_reason, ": ", detail_reason)
          }
          result_row$singular_fit <- singular_fit
          results[[idx]] <- result_row
          idx <- idx + 1L
          next
        }
      }

      contrast_error <- NA_character_
      contrast <- tryCatch(
        {
          emm <- emmeans::emmeans(fit, ~ SEX, weights = "equal")
          emmeans::contrast(emm, method = "pairwise", adjust = "none")
        },
        error = function(e) {
          contrast_error <<- conditionMessage(e)
          NULL
        }
      )

      if (is.null(contrast)) {
        result_row$failure_reason <- if (is.na(contrast_error) || !nzchar(contrast_error)) {
          "emmeans contrast failed"
        } else {
          paste0("emmeans contrast failed: ", contrast_error)
        }
        result_row$singular_fit <- singular_fit
        results[[idx]] <- result_row
        idx <- idx + 1L
        next
      }

      stats_rows <- as.data.frame(summary(contrast))
      if (nrow(stats_rows) == 0) {
        result_row$failure_reason <- "no pairwise sex contrasts available"
        result_row$singular_fit <- singular_fit
        results[[idx]] <- result_row
        idx <- idx + 1L
        next
      }

      for (k in seq_len(nrow(stats_rows))) {
        out_row <- result_row
        out_row$contrast <- as.character(stats_rows$contrast[[k]])
        out_row$used_random_intercept <- used_lmer
        out_row$model_type <- if (used_lmer) "lmer" else "lm"
        out_row$model_status <- if (used_lmer) "modeled" else "modeled_fallback"
        out_row$failure_reason <- if (!used_lmer && singular_fit) {
          "singular mixed model; used linear-model fallback"
        } else if (!used_lmer && !(n_donors >= 2 && repeated_donor)) {
          "random intercept unsupported; used linear-model fallback"
        } else {
          NA_character_
        }
        out_row$estimate <- stats_rows$estimate[[k]]
        out_row$SE <- stats_rows$SE[[k]]
        out_row$p.value <- stats_rows$p.value[[k]]
        out_row$singular_fit <- singular_fit
        results[[idx]] <- out_row
        idx <- idx + 1L
      }
    }
  }

  dplyr::bind_rows(results) %>%
    # FDR family: all proteins/sex-contrasts within each biological stratum.
    dplyr::group_by(AIRWAY, HORMONE, TIMEPOINT, EXPOSURE, INFECTION) %>%
    dplyr::mutate(
      q.value = {
        q_vals <- rep(NA_real_, dplyr::n())
        keep <- which(!is.na(p.value))
        if (length(keep) > 0) {
          q_vals[keep] <- p.adjust(p.value[keep], method = "fdr")
        }
        q_vals
      },
      significant = !is.na(q.value) & q.value < alpha_q
    ) %>%
    dplyr::ungroup()
}

fit_hpiv3_exposure_models <- function(data, protein_cols,
                                      protein_status = NULL,
                                      pseudocount = PSEUDOCOUNT,
                                      alpha_q = ALPHA_Q,
                                      min_nonmissing_per_group = 2L,
                                      min_total_nonmissing = 3L) {
  stopifnot(is.data.frame(data))
  protein_cols <- intersect(protein_cols, names(data))
  strata <- data %>% dplyr::distinct(AIRWAY, HORMONE, TIMEPOINT, SEX, INFECTION)

  if (is.null(protein_status)) {
    protein_status <- tibble::tibble(
      PROTEIN = protein_cols,
      included = TRUE,
      exclusion_reason = NA_character_
    )
  }
  protein_status <- protein_status %>%
    dplyr::select(PROTEIN, included, exclusion_reason)

  results <- list()
  idx <- 1L

  for (row_idx in seq_len(nrow(strata))) {
    airway_i <- strata$AIRWAY[[row_idx]]
    hormone_i <- strata$HORMONE[[row_idx]]
    timepoint_i <- strata$TIMEPOINT[[row_idx]]
    sex_i <- strata$SEX[[row_idx]]
    infection_i <- strata$INFECTION[[row_idx]]

    stratum_data <- data %>%
      dplyr::filter(
        AIRWAY == airway_i,
        HORMONE == hormone_i,
        TIMEPOINT == timepoint_i,
        SEX == sex_i,
        INFECTION == infection_i
      )

    for (protein in protein_cols) {
      protein_rule <- protein_status %>% dplyr::filter(PROTEIN == protein)
      included <- if (nrow(protein_rule) == 1) isTRUE(protein_rule$included[[1]]) else TRUE
      exclusion_reason <- if (nrow(protein_rule) == 1) protein_rule$exclusion_reason[[1]] else NA_character_

      dat_raw <- stratum_data %>%
        dplyr::transmute(
          PATIENTCODE = PATIENTCODE,
          EXPOSURE = as.character(EXPOSURE),
          VALUE = .data[[protein]]
        ) %>%
        dplyr::filter(!is.na(PATIENTCODE), !is.na(EXPOSURE), !is.na(VALUE))

      exp_order <- c(
        intersect("PBS_Control", unique(dat_raw$EXPOSURE)),
        sort(setdiff(unique(dat_raw$EXPOSURE), "PBS_Control"))
      )
      dat <- dat_raw %>%
        dplyr::mutate(EXPOSURE = factor(EXPOSURE, levels = exp_order))

      exposure_counts <- table(dat$EXPOSURE)
      n_obs <- nrow(dat)
      n_donors <- dplyr::n_distinct(dat$PATIENTCODE)
      repeated_donor <- anyDuplicated(as.character(dat$PATIENTCODE)) > 0
      counts_label <- if (length(exposure_counts) == 0) {
        NA_character_
      } else {
        paste0(names(exposure_counts), "=", as.integer(exposure_counts), collapse = "; ")
      }

      result_row <- tibble::tibble(
        AIRWAY = as.character(airway_i),
        HORMONE = as.character(hormone_i),
        TIMEPOINT = as.character(timepoint_i),
        SEX = as.character(sex_i),
        INFECTION = as.character(infection_i),
        PROTEIN = protein,
        comparison = "EXPOSURE",
        contrast = NA_character_,
        n_samples = n_obs,
        n_groups = length(exposure_counts),
        group_counts = counts_label,
        n_donors = n_donors,
        used_random_intercept = FALSE,
        model_type = NA_character_,
        model_status = "skipped",
        failure_reason = NA_character_,
        estimate = NA_real_,
        SE = NA_real_,
        p.value = NA_real_,
        singular_fit = NA
      )

      if (!included) {
        result_row$failure_reason <- exclusion_reason
        results[[idx]] <- result_row
        idx <- idx + 1L
        next
      }
      if (n_obs < min_total_nonmissing) {
        result_row$failure_reason <- paste0("fewer than ", min_total_nonmissing, " non-missing observations")
        results[[idx]] <- result_row
        idx <- idx + 1L
        next
      }

      valid_levels <- names(exposure_counts)[as.integer(exposure_counts) >= min_nonmissing_per_group]
      dat_fit <- dat %>%
        dplyr::filter(as.character(EXPOSURE) %in% valid_levels) %>%
        droplevels()

      if (length(unique(as.character(dat_fit$EXPOSURE))) < 2) {
        result_row$failure_reason <- paste0(
          "insufficient per-exposure observations (required >= ", min_nonmissing_per_group,
          " per level; observed: ", counts_label, ")"
        )
        results[[idx]] <- result_row
        idx <- idx + 1L
        next
      }

      dat_fit <- dat_fit %>% dplyr::mutate(log2_value = log2(as.numeric(VALUE) + pseudocount))

      fit <- NULL
      used_lmer <- FALSE
      singular_fit <- FALSE
      lmer_error <- NA_character_
      lm_error <- NA_character_

      if (n_donors >= 2 && repeated_donor) {
        fit <- tryCatch(
          lme4::lmer(log2_value ~ EXPOSURE + (1 | PATIENTCODE), data = dat_fit),
          error = function(e) {
            lmer_error <<- conditionMessage(e)
            NULL
          }
        )
        if (!is.null(fit)) {
          singular_fit <- isTRUE(tryCatch(lme4::isSingular(fit), error = function(...) FALSE))
          if (!singular_fit) {
            used_lmer <- TRUE
          }
        }
      }

      if (!used_lmer) {
        fit <- tryCatch(
          stats::lm(log2_value ~ EXPOSURE, data = dat_fit),
          error = function(e) {
            lm_error <<- conditionMessage(e)
            NULL
          }
        )
        if (is.null(fit)) {
          base_reason <- if (singular_fit) {
            "singular mixed model and fallback linear model failed"
          } else {
            "model fitting failed"
          }
          detail_reason <- lm_error
          if (is.na(detail_reason) || !nzchar(detail_reason)) {
            detail_reason <- lmer_error
          }
          result_row$failure_reason <- if (is.na(detail_reason) || !nzchar(detail_reason)) {
            base_reason
          } else {
            paste0(base_reason, ": ", detail_reason)
          }
          result_row$singular_fit <- singular_fit
          results[[idx]] <- result_row
          idx <- idx + 1L
          next
        }
      }

      contrast_error <- NA_character_
      contrast <- tryCatch(
        {
          emm <- emmeans::emmeans(fit, ~ EXPOSURE, weights = "equal")
          emmeans::contrast(emm, method = "pairwise", adjust = "none")
        },
        error = function(e) {
          contrast_error <<- conditionMessage(e)
          NULL
        }
      )

      if (is.null(contrast)) {
        result_row$failure_reason <- if (is.na(contrast_error) || !nzchar(contrast_error)) {
          "emmeans contrast failed"
        } else {
          paste0("emmeans contrast failed: ", contrast_error)
        }
        result_row$singular_fit <- singular_fit
        results[[idx]] <- result_row
        idx <- idx + 1L
        next
      }

      stats_rows <- as.data.frame(summary(contrast))
      if (nrow(stats_rows) == 0) {
        result_row$failure_reason <- "no pairwise exposure contrasts available"
        result_row$singular_fit <- singular_fit
        results[[idx]] <- result_row
        idx <- idx + 1L
        next
      }

      for (k in seq_len(nrow(stats_rows))) {
        out_row <- result_row
        out_row$contrast <- as.character(stats_rows$contrast[[k]])
        out_row$used_random_intercept <- used_lmer
        out_row$model_type <- if (used_lmer) "lmer" else "lm"
        out_row$model_status <- if (used_lmer) "modeled" else "modeled_fallback"
        out_row$failure_reason <- if (!used_lmer && singular_fit) {
          "singular mixed model; used linear-model fallback"
        } else if (!used_lmer && !(n_donors >= 2 && repeated_donor)) {
          "random intercept unsupported; used linear-model fallback"
        } else {
          NA_character_
        }
        out_row$estimate <- stats_rows$estimate[[k]]
        out_row$SE <- stats_rows$SE[[k]]
        out_row$p.value <- stats_rows$p.value[[k]]
        out_row$singular_fit <- singular_fit
        results[[idx]] <- out_row
        idx <- idx + 1L
      }
    }
  }

  dplyr::bind_rows(results) %>%
    # FDR family: all proteins/exposure-contrasts within each biological stratum.
    dplyr::group_by(AIRWAY, HORMONE, TIMEPOINT, SEX, INFECTION) %>%
    dplyr::mutate(
      q.value = {
        q_vals <- rep(NA_real_, dplyr::n())
        keep <- which(!is.na(p.value))
        if (length(keep) > 0) {
          q_vals[keep] <- p.adjust(p.value[keep], method = "fdr")
        }
        q_vals
      },
      significant = !is.na(q.value) & q.value < alpha_q
    ) %>%
    dplyr::ungroup()
}

summarize_hpiv3_strata <- function(data, protein_cols) {
  stopifnot(is.data.frame(data))
  protein_cols <- intersect(protein_cols, names(data))

  long <- data %>%
    dplyr::select(AIRWAY, HORMONE, TIMEPOINT, EXPOSURE, SEX, INFECTION, dplyr::all_of(protein_cols)) %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(protein_cols),
      names_to = "PROTEIN",
      values_to = "VALUE"
    )

  pooled <- long %>%
    dplyr::mutate(
      SEX = "All",
      EXPOSURE = as.character(EXPOSURE)
    )

  dplyr::bind_rows(long, pooled) %>%
    dplyr::group_by(AIRWAY, HORMONE, TIMEPOINT, EXPOSURE, SEX, INFECTION, PROTEIN) %>%
    dplyr::summarise(
      n = sum(!is.na(VALUE)),
      mean_pg_ml = mean(VALUE, na.rm = TRUE),
      sd_pg_ml = stats::sd(VALUE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      mean_pg_ml = ifelse(is.nan(mean_pg_ml), NA_real_, mean_pg_ml),
      sd_pg_ml = ifelse(is.na(sd_pg_ml) | is.nan(sd_pg_ml), 0, sd_pg_ml),
      summary_pg_ml = ifelse(
        is.na(mean_pg_ml),
        NA_character_,
        sprintf("%.3f ± %.3f", mean_pg_ml, sd_pg_ml)
      )
    )
}

summarize_unique_significant_proteins <- function(
    model_results,
    group_var,
    group_levels = NULL,
    filters = list(),
    alpha_q = ALPHA_Q,
    membership_source = c("column", "contrast"),
    contrast_col = "contrast",
    estimate_col = "estimate"
) {
  membership_source <- match.arg(membership_source)
  stopifnot(is.data.frame(model_results), group_var %in% names(model_results))
  stopifnot("PROTEIN" %in% names(model_results))
  stopifnot(is.list(filters))

  d <- model_results
  for (col in names(filters)) {
    stopifnot(col %in% names(d))
    d <- d[as.character(d[[col]]) %in% as.character(filters[[col]]), , drop = FALSE]
  }

  empty_membership <- tibble::tibble(
    PROTEIN = character(),
    sig_sets = list(),
    n_sig_groups = integer(),
    membership_class = character()
  )
  empty_membership_long <- tibble::tibble(
    PROTEIN = character(),
    sig_group = character(),
    n_sig_groups = integer(),
    membership_class = character()
  )
  empty_summary <- tibble::tibble(
    membership_class = character(),
    n_proteins = integer()
  )

  if (nrow(d) == 0) {
    if (is.null(group_levels)) {
      group_levels <- character()
    } else {
      group_levels <- as.character(group_levels)
    }
    return(list(
      membership = empty_membership,
      membership_long = empty_membership_long,
      summary = empty_summary,
      group_var = group_var,
      group_levels = group_levels,
      filters = filters,
      alpha_q = alpha_q
    ))
  }

  if (!("significant" %in% names(d))) {
    if ("q.value" %in% names(d)) {
      d$significant <- !is.na(d$q.value) & d$q.value < alpha_q
    } else {
      d$significant <- FALSE
    }
  }

  if (membership_source == "column") {
    if (is.null(group_levels)) {
      group_levels <- sort(unique(as.character(d[[group_var]])))
    } else {
      group_levels <- as.character(group_levels)
      d <- d[as.character(d[[group_var]]) %in% group_levels, , drop = FALSE]
    }

    d <- d %>%
      dplyr::mutate(
        significant = dplyr::coalesce(significant, FALSE),
        !!group_var := as.character(.data[[group_var]])
      )

    membership <- d %>%
      dplyr::group_by(PROTEIN) %>%
      dplyr::summarise(
        sig_sets = list(sort(unique(.data[[group_var]][significant]))),
        .groups = "drop"
      )
  } else {
    stopifnot(contrast_col %in% names(d), estimate_col %in% names(d))

    parsed_levels <- d %>%
      dplyr::transmute(
        .contrast = as.character(.data[[contrast_col]]),
        .split = stringr::str_split(.contrast, " - ", n = 2)
      ) %>%
      dplyr::mutate(
        group_first = purrr::map_chr(.split, ~ if (length(.x) == 2) .x[[1]] else NA_character_),
        group_second = purrr::map_chr(.split, ~ if (length(.x) == 2) .x[[2]] else NA_character_)
      ) %>%
      dplyr::select(group_first, group_second)

    if (is.null(group_levels)) {
      group_levels <- sort(unique(c(parsed_levels$group_first, parsed_levels$group_second)))
      group_levels <- group_levels[!is.na(group_levels)]
    } else {
      group_levels <- as.character(group_levels)
    }

    d <- d %>%
      dplyr::mutate(
        significant = dplyr::coalesce(significant, FALSE),
        .contrast = as.character(.data[[contrast_col]]),
        .split = stringr::str_split(.contrast, " - ", n = 2),
        group_first = purrr::map_chr(.split, ~ if (length(.x) == 2) .x[[1]] else NA_character_),
        group_second = purrr::map_chr(.split, ~ if (length(.x) == 2) .x[[2]] else NA_character_),
        sig_group = dplyr::case_when(
          !significant ~ NA_character_,
          is.na(.data[[estimate_col]]) ~ NA_character_,
          .data[[estimate_col]] > 0 ~ group_first,
          .data[[estimate_col]] < 0 ~ group_second,
          TRUE ~ NA_character_
        )
      )

    if (length(group_levels) > 0) {
      d <- d %>%
        dplyr::mutate(sig_group = ifelse(sig_group %in% group_levels, sig_group, NA_character_))
    }

    membership <- d %>%
      dplyr::group_by(PROTEIN) %>%
      dplyr::summarise(
        sig_sets = list(sort(unique(stats::na.omit(sig_group)))),
        .groups = "drop"
      )
  }

  membership <- membership %>%
    dplyr::mutate(
      n_sig_groups = vapply(sig_sets, length, integer(1)),
      unique_group = purrr::map_chr(
        sig_sets,
        function(x) if (length(x) > 0) x[[1]] else NA_character_
      ),
      membership_class = dplyr::case_when(
        n_sig_groups == 0 ~ "Not significant",
        n_sig_groups == 1 ~ paste0("Unique: ", unique_group),
        TRUE ~ paste0("Shared (", n_sig_groups, "-way)")
      )
    ) %>%
    dplyr::select(-unique_group)

  membership_long <- membership %>%
    dplyr::transmute(
      PROTEIN = PROTEIN,
      sig_group = purrr::map(
        sig_sets,
        function(x) if (length(x) == 0) NA_character_ else x
      ),
      n_sig_groups = n_sig_groups,
      membership_class = membership_class
    ) %>%
    tidyr::unnest_longer(sig_group)

  summary_counts <- membership %>%
    dplyr::count(membership_class, name = "n_proteins") %>%
    dplyr::arrange(dplyr::desc(n_proteins), membership_class)

  list(
    membership = membership,
    membership_long = membership_long,
    summary = summary_counts,
    group_var = group_var,
    group_levels = group_levels,
    filters = filters,
    alpha_q = alpha_q
  )
}
