# =============================================================================
# functions_plots.R
# All ggplot2 helpers: dotplots, bar plots, histograms, boxplots,
# gene expression bars, GSEA barplots
# =============================================================================

UP_COLOR_DEFAULT   <- "#D7191C"
DOWN_COLOR_DEFAULT <- "#2C7BB6"
ESTRADIOL_FILL_DEFAULT <- "#E8F5E9"

# ── Label helpers ───────────────────────────────────────────────────────────────

short_exposure_label <- function(exposure) {
  exposure %>%
    gsub("_", " ", ., fixed = TRUE) %>%
    sub(" \\d+$", "", .)
}

sig_label_from_q <- function(q, alpha_q = ALPHA_Q) {
  dplyr::case_when(
    is.na(q)   ~ "",
    q < 0.001  ~ "***",
    q < 0.01   ~ "**",
    q < alpha_q ~ "*",
    TRUE       ~ ""
  )
}

# ── Size rescaling helper ─────────────────────────────────────────────────────

add_rescaled_size <- function(df, value_col,
                              clamp_quantile = 0.95,
                              transform = c("sqrt","log1p","none"),
                              out_col = "size_val") {
  transform <- match.arg(transform)
  mag  <- abs(df[[value_col]])
  cap  <- stats::quantile(mag, probs = clamp_quantile, na.rm = TRUE, names = FALSE)
  if (!is.finite(cap) || cap <= 0) cap <- max(mag, na.rm = TRUE)
  mag2 <- pmin(mag, cap)
  df[[out_col]] <- switch(transform,
                          sqrt  = sqrt(mag2),
                          log1p = log1p(mag2),
                          none  = mag2)
  df
}

# ── Cytokine dotplot (pooled or sex-stratified) ───────────────────────────────
# Call once per exposure inside a loop; returns a combined cowplot grid.
# Requires: plot_df_pooled, plot_df_sex, cyt_order, cyt_levels,
#           border_color, panel_fill_none, e2_label, hormone_levels
#           (all set in the calling chunk)

make_dotplot <- function(df, x_var, title_str,
                         cyt_levels, border_color,
                         up_color   = UP_COLOR_DEFAULT,
                         down_color = DOWN_COLOR_DEFAULT,
                         show_y = TRUE) {
  
  fill_scale  <- ggplot2::scale_fill_manual(
    name   = "Direction\n(vs matched PBS)",
    values = c("Up" = up_color, "Down" = down_color, "Zero" = "grey70"),
    na.value = "grey70")
  color_scale <- ggplot2::scale_color_manual(
    name   = "Direction\n(vs matched PBS)",
    values = c("Up" = up_color, "Down" = down_color, "Zero" = "grey70"),
    na.value = "grey70")
  size_scale  <- ggplot2::scale_size_continuous(
    name  = "Effect size\n(sqrt|mean log2FC|)",
    range = c(1, 8), limits = c(0, NA))
  
  ref_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = NA, colour = border_color,
                                               linewidth = 1.2),
      panel.ontop      = TRUE,
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = border_color,
                                               colour = border_color, linewidth = 1.2),
      strip.text       = ggplot2::element_text(color = "white", face = "bold", size = 11),
      axis.title.x     = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(face = "bold", size = 11),
      axis.text.y      = ggplot2::element_text(face = "bold", size = 10),
      axis.title.y     = ggplot2::element_text(size = 11),
      legend.position  = "right",
      plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA),
      plot.title       = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5)
    )
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[x_var]], y = CYTOKINE)) +
    ggplot2::geom_tile(ggplot2::aes(fill = panel_fill), width = Inf, height = 1,
                       alpha = 1, show.legend = FALSE) +
    ggplot2::scale_fill_identity() +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_hline(yintercept = seq_along(cyt_levels),
                        color = "white", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = seq_along(unique(df[[x_var]])),
                        color = "white", linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(fill = direction, color = direction,
                                     size = size_val),
                        shape = 21, stroke = 1.5, alpha = 0.9) +
    fill_scale + color_scale + size_scale +
    ggplot2::guides(
      color = "none",
      fill  = ggplot2::guide_legend(title = "Direction\n(vs matched PBS)"),
      size  = ggplot2::guide_legend(title = "Effect size\n(sqrt|mean log2FC|)")) +
    ggplot2::facet_grid(rows = ggplot2::vars(CELLTYPE),
                        cols = ggplot2::vars(HORMONE), drop = FALSE) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = c(0.8, 0.8)),
                              drop = TRUE) +
    ggplot2::scale_y_discrete(drop = FALSE) +
    ref_theme +
    ggplot2::labs(title = title_str, y = "Cytokine")
  
  if (!show_y)
    p <- p + ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank())
  p
}

make_cytokine_barplot <- function(summary_raw, d_raw, cyt,
                                  target_exposure, exp_short,
                                  border_color) {
  
  y_max     <- max(summary_raw$mu + summary_raw$sd, na.rm = TRUE)
  y_bracket <- y_max * 1.10
  bg_df     <- summary_raw %>% dplyr::distinct(CELLTYPE, HORMONE, panel_fill)
  
  bracket_df <- summary_raw %>%
    dplyr::filter(EXPOSURE == target_exposure, sig_label != "") %>%
    dplyr::mutate(
      x_left   = as.numeric(x_label) - 1,
      x_right  = as.numeric(x_label),
      x_center = as.numeric(x_label) - 0.5,
      y_top    = y_bracket
    )
  
  bar_theme <- ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = NA, colour = border_color,
                                               linewidth = 1),
      panel.ontop      = TRUE,
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = border_color,
                                               colour = border_color, linewidth = 1),
      strip.text       = ggplot2::element_text(color = "white", face = "bold", size = 11),
      axis.title.x     = ggplot2::element_blank(),
      axis.text.x      = ggplot2::element_text(face = "bold", size = 9,
                                               angle = 30, hjust = 1),
      axis.text.y      = ggplot2::element_text(face = "bold", size = 10),
      axis.title.y     = ggplot2::element_text(size = 11),
      legend.position  = "right",
      plot.background  = ggplot2::element_rect(fill = "transparent", colour = NA),
      plot.title       = ggplot2::element_text(face = "bold", size = 13, hjust = 0.5)
    )
  
  p <- ggplot2::ggplot(summary_raw, ggplot2::aes(x = x_label, y = mu)) +
    ggplot2::geom_rect(
      data = bg_df, ggplot2::aes(fill = panel_fill),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      inherit.aes = FALSE, show.legend = FALSE) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_hline(yintercept = pretty(c(0, y_max)),
                        color = "white", linewidth = 0.5) +
    ggnewscale::new_scale_fill() +
    
    ggplot2::geom_col(
      data = summary_raw %>% dplyr::filter(SEX == "M"),
      ggplot2::aes(x = x_label, y = mu, fill = bar_fill, color = bar_color),
      width = 0.9, linewidth = 0.7, alpha = 1) +
    
    ggplot2::geom_col(
      data = summary_raw %>% dplyr::filter(SEX == "F"),
      ggplot2::aes(x = x_label, y = mu, fill = bar_fill, color = bar_color),
      width = 0.9, linewidth = 0.7, alpha = 0) +
    
    ggplot2::scale_fill_identity(guide = "none") +
    ggplot2::scale_color_identity(guide = "none") +
    
    ggplot2::geom_errorbar(
      ggplot2::aes(x = x_label, ymin = mu - sd, ymax = mu + sd),
      width = 0.2, linewidth = 0.6, color = "grey20") +
    
    { if (nrow(bracket_df) > 0) list(
      ggplot2::geom_segment(
        data = bracket_df,
        ggplot2::aes(x = x_left, xend = x_right,
                     y = y_top, yend = y_top, color = sig_color),
        linewidth = 0.7, inherit.aes = FALSE, show.legend = FALSE),
      ggplot2::geom_segment(
        data = bracket_df,
        ggplot2::aes(x = x_left, xend = x_left,
                     y = y_top, yend = y_top - y_max * 0.04,
                     color = sig_color),
        linewidth = 0.7, inherit.aes = FALSE, show.legend = FALSE),
      ggplot2::geom_segment(
        data = bracket_df,
        ggplot2::aes(x = x_right, xend = x_right,
                     y = y_top, yend = y_top - y_max * 0.04,
                     color = sig_color),
        linewidth = 0.7, inherit.aes = FALSE, show.legend = FALSE),
      ggplot2::geom_text(
        data = bracket_df,
        ggplot2::aes(x = x_center,
                     y = y_top + y_max * 0.03,
                     label = sig_label, color = sig_color),
        size = 6, fontface = "bold", hjust = 0.5, vjust = 0,
        inherit.aes = FALSE, show.legend = FALSE)
    ) else list() } +
    
    ggplot2::geom_jitter(
      data = d_raw,
      ggplot2::aes(x = x_label, y = VALUE,
                   fill = point_fill, color = point_color),
      shape = 21, width = 0.12, size = 2.5,
      stroke = 1.0, alpha = 0.95, inherit.aes = FALSE) +
    ggplot2::scale_fill_identity(guide  = "none") +
    ggplot2::scale_color_identity(guide = "none") +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.6)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.18))) +
    ggplot2::facet_grid(rows = ggplot2::vars(CELLTYPE),
                        cols = ggplot2::vars(HORMONE)) +
    bar_theme +
    ggplot2::labs(
      title = paste0(cyt, " — ", exp_short,
                     " vs Control | M = filled, F = transparent"),
      y = paste0(cyt, " (pg/mL)"))
  
  p
}

make_cytokine_barplot_single_combo <- function(
    df, cytokine, target_exposure, sex, airway,
    pbs_level      = PBS_LEVEL,
    up_color       = UP_COLOR_DEFAULT,
    down_color     = DOWN_COLOR_DEFAULT,
    estradiol_fill = ESTRADIOL_FILL_DEFAULT,
    sig_data       = NULL
) {
  stopifnot(cytokine %in% names(df))

  exp_short       <- short_exposure_label(target_exposure)
  e2_label        <- paste(exp_short, "+ E2")
  hormone_levels  <- c(exp_short, e2_label)
  border_color    <- EXPOSURE_COLORS_DEEP[[target_exposure]]
  panel_fill_none <- EXPOSURE_COLORS_LIGHT[[target_exposure]]

  d_raw <- dplyr::bind_rows(
    df %>%
      dplyr::filter(SEX == sex, CELLTYPE == airway,
                    EXPOSURE == target_exposure,
                    HORMONE %in% c("NONE", "Estradiol")),
    df %>%
      dplyr::filter(SEX == sex, CELLTYPE == airway,
                    EXPOSURE == pbs_level,
                    HORMONE %in% c("NONE", "Estradiol"))
  ) %>%
    dplyr::mutate(
      HORMONE  = dplyr::recode(as.character(HORMONE),
                               "NONE" = exp_short, "Estradiol" = e2_label),
      HORMONE  = factor(HORMONE, levels = hormone_levels),
      SEX      = factor(as.character(SEX),      levels = c("M", "F")),
      CELLTYPE = factor(as.character(CELLTYPE),  levels = c("LAE", "SAE")),
      VALUE    = as.numeric(.data[[cytokine]]),
      EXPOSURE = factor(as.character(EXPOSURE),  levels = c(pbs_level, target_exposure))
    )

  if (nrow(d_raw) == 0) return(NULL)

  male_pt_fill <- grDevices::adjustcolor(border_color, alpha.f = 1,
                                         red.f = 1.3, green.f = 1.3, blue.f = 1.3)

  d_raw <- d_raw %>%
    dplyr::mutate(
      x_label = factor(
        dplyr::case_when(
          EXPOSURE == pbs_level       ~ "Control",
          EXPOSURE == target_exposure ~ exp_short,
          TRUE ~ as.character(EXPOSURE)
        ),
        levels = c("Control", exp_short)
      ),
      point_fill = dplyr::case_when(
        SEX == "M" & EXPOSURE == pbs_level       ~ "grey80",
        SEX == "M" & EXPOSURE == target_exposure ~ male_pt_fill,
        TRUE                                     ~ "white"
      ),
      point_color = dplyr::case_when(
        EXPOSURE == pbs_level       ~ "grey40",
        EXPOSURE == target_exposure ~ border_color,
        TRUE                        ~ "grey40"
      )
    )

  summary_raw <- d_raw %>%
    dplyr::group_by(SEX, EXPOSURE, HORMONE, CELLTYPE) %>%
    dplyr::summarise(
      mu = mean(VALUE, na.rm = TRUE),
      sd = stats::sd(VALUE, na.rm = TRUE),
      n  = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      panel_fill = ifelse(as.character(HORMONE) == e2_label,
                          estradiol_fill, panel_fill_none),
      bar_fill   = border_color,
      bar_color  = border_color,
      direction  = dplyr::case_when(
        mu > 0 ~ "Up", mu < 0 ~ "Down", TRUE ~ "Zero"
      ),
      sig_color  = dplyr::case_when(
        direction == "Up"   ~ up_color,
        direction == "Down" ~ down_color,
        TRUE                ~ "grey50"
      ),
      x_label = factor(
        dplyr::case_when(
          EXPOSURE == pbs_level       ~ "Control",
          EXPOSURE == target_exposure ~ exp_short,
          TRUE ~ as.character(EXPOSURE)
        ),
        levels = c("Control", exp_short)
      ),
      q = NA_real_,
      sig = FALSE,
      sig_label = ""
    )

  if (!is.null(sig_data)) {
    required_sig_cols <- c("CELLTYPE", "HORMONE", "SEX", "EXPOSURE", "CYTOKINE", "q", "sig")
    missing_sig_cols <- setdiff(required_sig_cols, names(sig_data))
    if (length(missing_sig_cols) > 0) {
      stop("sig_data is missing required columns: ", paste(missing_sig_cols, collapse = ", "))
    }

    sig_data <- sig_data %>%
      dplyr::mutate(
        SEX = factor(as.character(SEX), levels = c("M", "F")),
        CELLTYPE = factor(as.character(CELLTYPE), levels = c("LAE", "SAE")),
        HORMONE = factor(
          dplyr::recode(as.character(HORMONE),
                        "NONE" = exp_short, "Estradiol" = e2_label),
          levels = hormone_levels
        ),
        EXPOSURE = factor(as.character(EXPOSURE), levels = c(pbs_level, target_exposure)),
        q = as.numeric(q),
        sig = dplyr::coalesce(sig, q < ALPHA_Q, FALSE)
      ) %>%
      dplyr::filter(
        CYTOKINE == cytokine,
        EXPOSURE == target_exposure,
        SEX == sex,
        CELLTYPE == airway
      ) %>%
      dplyr::group_by(CELLTYPE, HORMONE, SEX, EXPOSURE) %>%
      dplyr::summarise(
        q = suppressWarnings(min(q, na.rm = TRUE)),
        sig = any(sig %in% TRUE, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        q = dplyr::if_else(is.infinite(q), NA_real_, q),
        sig_label = sig_label_from_q(q, alpha_q = ALPHA_Q)
      )

    summary_raw <- summary_raw %>%
      dplyr::left_join(
        sig_data %>%
          dplyr::select(CELLTYPE, HORMONE, SEX, EXPOSURE, q, sig, sig_label) %>%
          dplyr::rename(q_sig = q, sig_sig = sig, sig_label_sig = sig_label),
        by = c("CELLTYPE", "HORMONE", "SEX", "EXPOSURE")
      ) %>%
      dplyr::mutate(
        sig_label = dplyr::coalesce(sig_label_sig, sig_label),
        sig       = dplyr::coalesce(sig_sig, sig, FALSE),
        q         = dplyr::coalesce(q_sig, q)
      ) %>%
      dplyr::select(-q_sig, -sig_sig, -sig_label_sig)
  }

  make_cytokine_barplot(
    summary_raw     = summary_raw,
    d_raw           = d_raw,
    cyt             = cytokine,
    target_exposure = target_exposure,
    exp_short       = exp_short,
    border_color    = border_color
  ) +
    ggplot2::labs(
      title = paste0(cytokine, " — ", target_exposure,
                     " | SEX=", sex, " | AIRWAY=", airway)
    )
}

plot_all_cytokine_single_combo_bars <- function(
    df,
    cytokines,
    exposures  = ALL_HIGH_DOSE_EXPOSURES,
    sexes      = c("M", "F"),
    airways    = c("LAE", "SAE"),
    sig_data   = NULL,
    save       = TRUE,
    out_subdir = "single_combo_barplots",
    width      = 8,
    height     = 5,
    dpi        = 150
) {
  out <- list()

  if (!is.null(sig_data)) {
    required_sig_cols <- c("CELLTYPE", "HORMONE", "SEX", "EXPOSURE", "CYTOKINE", "q", "sig")
    missing_sig_cols <- setdiff(required_sig_cols, names(sig_data))
    if (length(missing_sig_cols) > 0) {
      stop("sig_data is missing required columns: ", paste(missing_sig_cols, collapse = ", "))
    }
  }

  if (save) {
    ensure_dir(here::here(PATH_OUTPUT_FIGS, out_subdir))
  }

  for (cyt in cytokines) {
    for (exp in exposures) {
      for (sx in sexes) {
        for (aw in airways) {
          sig_subset <- NULL
          if (!is.null(sig_data)) {
            sig_subset <- sig_data %>%
              dplyr::filter(CYTOKINE == cyt, EXPOSURE == exp, SEX == sx, CELLTYPE == aw)
          }

          p <- make_cytokine_barplot_single_combo(
            df              = df,
            cytokine        = cyt,
            target_exposure = exp,
            sex             = sx,
            airway          = aw,
            sig_data        = sig_subset
          )

          key      <- paste(cyt, exp, sx, aw, sep = "__")
          out[[key]] <- p

          if (!is.null(p)) {
            print(p)
            if (save) {
              fn <- paste0(
                "barplot_", safe_name(cyt), "_",
                safe_name(exp), "_SEX_", sx, "_AIRWAY_", aw, ".png"
              )
              save_plot(
                filename = file.path(out_subdir, fn),
                plot     = p,
                width    = width,
                height   = height,
                dpi      = dpi,
                bg       = "white"
              )
            }
          }
        }
      }
    }
  }

  invisible(out)
}

# ── Save helpers ──────────────────────────────────────────────────────────────

save_plot <- function(filename, plot = ggplot2::last_plot(),
                      width = 10, height = 8, dpi = 300,
                      path = PATH_OUTPUT_FIGURES, bg = "white") {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)

  ggplot2::ggsave(
    filename = file.path(path, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    bg = bg
  )

  cat("✓ Saved plot:", file.path(path, filename), "\n")
  invisible(TRUE)
}

save_table <- function(data, filename, path = PATH_OUTPUT_TABLES) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)

  readr::write_csv(data, file = file.path(path, filename))

  cat("✓ Saved table:", file.path(path, filename), "\n")
  invisible(TRUE)
}

# ── Per-cytokine histogram ────────────────────────────────────────────────────

make_cytokine_histogram <- function(d_raw, cyt, target_exposure) {
  ggplot2::ggplot(
    d_raw %>% dplyr::filter(EXPOSURE == target_exposure),
    ggplot2::aes(x = VALUE, fill = SEX)
  ) +
    ggplot2::geom_histogram(position = "identity", alpha = 0.65,
                            bins = 20, color = "grey40") +
    ggplot2::facet_grid(rows = ggplot2::vars(CELLTYPE),
                        cols = ggplot2::vars(HORMONE)) +
    ggplot2::scale_fill_manual(values = c("M" = "lemonchiffon", "F" = "white"),
                               name = "Sex") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = paste("Raw Values —", cyt, "|", target_exposure),
      x = paste(cyt, "(pg/mL)"), y = "Count")
}

# ── Boxplot helpers (legacy plotting functions) ───────────────────────────────

sex_exposure_plot <- function(data, y_variable) {
  p <- data %>%
    dplyr::filter(SEX %in% c("M","F")) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = EXPOSURE, y = .data[[y_variable]],
      color = EXPOSURE, fill = factor(SEX), pattern = CONCENTRATION)) +
    ggplot2::geom_boxplot(linewidth = 1) +
    ggpattern::geom_boxplot_pattern(
      position           = ggplot2::position_dodge(preserve = "single"),
      pattern_fill       = "black", pattern_angle = 45,
      pattern_density    = 0.2, pattern_spacing = 0.025,
      pattern_key_scale_factor = 0.6) +
    ggplot2::guides(
      fill  = ggplot2::guide_legend(override.aes = list(pattern = "none")),
      color = ggplot2::guide_legend(override.aes = list(pattern = "none"))) +
    ggpattern::scale_pattern_manual(
      values = c(HIGH = "stripe", LOW = "none", NONE = "none")) +
    ggplot2::geom_jitter(
      position = ggplot2::position_jitterdodge(), size = 2, alpha = 1) +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(y = y_variable, title = paste(y_variable, "and Sex")) +
    ggplot2::scale_fill_manual(values = c("M" = "lemonchiffon", "F" = "thistle1")) +
    ggplot2::scale_color_manual(values = color_exposure) +
    ggplot2::theme(
      axis.title   = ggplot2::element_text(size = 20, face = "bold"),
      axis.text.x  = ggplot2::element_blank(),
      axis.text.y  = ggplot2::element_text(face = "bold", color = "black", size = 20),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line    = ggplot2::element_line(colour = "black", linewidth = 1),
      panel.background = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::scale_x_discrete(limits = c(
      "PBS_Control","Untreated_Control","Peat_5","Peat_25","Pine_5","Pine_25",
      "Eucalyptus_5","Eucalyptus_25","RedOak_5","RedOak_25"))
  p
}

exposure_plot <- function(data, title, y_variable) {
  p <- data %>%
    ggplot2::ggplot(ggplot2::aes(
      x = EXPOSURE, y = .data[[y_variable]],
      color = EXPOSURE, pattern = CONCENTRATION)) +
    ggplot2::geom_boxplot(linewidth = 1) +
    ggpattern::geom_boxplot_pattern(
      position           = ggplot2::position_dodge(preserve = "single"),
      pattern_fill       = "black", pattern_angle = 45,
      pattern_density    = 0.2, pattern_spacing = 0.025,
      pattern_key_scale_factor = 0.6) +
    ggplot2::guides(
      fill  = ggplot2::guide_legend(override.aes = list(pattern = "none")),
      color = ggplot2::guide_legend(override.aes = list(pattern = "none"))) +
    ggpattern::scale_pattern_manual(
      values = c(HIGH = "stripe", LOW = "none", NONE = "none")) +
    ggplot2::geom_jitter(
      position = ggplot2::position_jitterdodge(), size = 2, alpha = 1) +
    ggplot2::scale_fill_manual(values = c("M" = "lemonchiffon", "F" = "thistle1")) +
    ggplot2::scale_color_manual(values = color_exposure) +
    ggplot2::theme(
      legend.position  = "right",
      axis.title       = ggplot2::element_text(size = 20, face = "bold"),
      axis.text.x      = ggplot2::element_blank(),
      axis.text.y      = ggplot2::element_text(face = "bold", color = "black", size = 20),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(colour = "black", linewidth = 1),
      panel.background = ggplot2::element_blank(),
      axis.ticks.x     = ggplot2::element_blank()) +
    ggplot2::labs(y = y_variable, title = title) +
    ggplot2::scale_x_discrete(limits = c(
      "Untreated_Control","Peat_5","Peat_25","Pine_5","Pine_25",
      "Eucalyptus_5","Eucalyptus_25","RedOak_5","RedOak_25")) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray1")
  p
}

# ── Gene expression bar plots (RNA-seq) ──────────────────────────────────────

plot_gene_expression_bars <- function(vst_data, metadata, genes) {
  long_dat <- vst_data %>%
    dplyr::filter(Geneid %in% genes) %>%
    tidyr::pivot_longer(cols = -Geneid,
                        names_to = "SAMPLEID", values_to = "Expression") %>%
    dplyr::left_join(metadata, by = "SAMPLEID") %>%
    dplyr::mutate(EXPOSURE = factor(EXPOSURE))
  
  ggplot2::ggplot(long_dat,
                  ggplot2::aes(x = EXPOSURE, y = Expression, fill = SEX)) +
    ggplot2::stat_summary(fun = mean, geom = "bar",
                          position = ggplot2::position_dodge(0.8), width = 0.7) +
    ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar",
                          position = ggplot2::position_dodge(0.8), width = 0.25) +
    ggplot2::scale_fill_manual(values = c("M" = "steelblue", "F" = "salmon")) +
    ggplot2::facet_wrap(~ Geneid, scales = "free_y") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text       = ggplot2::element_text(face = "bold"),
      panel.border     = ggplot2::element_rect(color = "black", fill = NA,
                                               linewidth = 0.5)) +
    ggplot2::labs(y = "VST Expression", x = NULL)
}

# ── GSEA barplot helper (RNA-seq) ─────────────────────────────────────────────

plot_gsea_barplot <- function(gsea_data, n_top = 12, facet_by = "sample_name",
                              color_low = "darkgreen", color_high = "gray",
                              fdr_limit = 0.35) {
  plot_dat <- gsea_data %>%
    dplyr::mutate(
      clean_path = gsub("^HALLMARK_", "", Pathway),
      clean_path = gsub("_", " ", clean_path)
    ) %>%
    dplyr::group_by(.data[[facet_by]]) %>%
    dplyr::arrange(FDR_qval) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::ungroup()
  
  ggplot2::ggplot(plot_dat,
                  ggplot2::aes(x = reorder(clean_path, NES), y = NES, fill = FDR_qval)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_gradient(low = color_low, high = color_high,
                                 limits = c(0, fdr_limit)) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)),
                        scales = "free_y") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(face = "bold"),
      axis.text.y      = ggplot2::element_text(face = "bold"),
      strip.text       = ggplot2::element_text(size = 12, face = "bold"),
      panel.border     = ggplot2::element_rect(color = "black", fill = NA,
                                               linewidth = 1),
      panel.background = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()) +
    ggplot2::labs(x = "Pathway",
                  y = "Normalized Enrichment Score",
                  fill = "FDR q-val")
}
