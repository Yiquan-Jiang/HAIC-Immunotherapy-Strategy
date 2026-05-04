#!/usr/bin/env Rscript
# =============================================================================
# make_publication_figures.R
#
# Publication-ready figures for RCS non-linear interaction analysis.
#
# Comparisons (two PSM pairs):
#   Pair 02: HAIC_then_I  vs HAIC_alone
#   Pair 06: HAIC_then_I+T vs HAIC_alone
#
# Biomarkers (5): AFP, PIVKA, SII, PLR, NLR
# Timepoints (3): Baseline, Pre-IT, Pre-IT Change Rate
#
# Outputs (in ./ relative to this script):
#   Fig_Main_RCS_5indicators.{pdf,png,tiff}             5 x 6 matrix (main text)
#   Fig_Supp_RCS_Baseline.{pdf,png,tiff}                5 x 2, baseline
#   Fig_Supp_RCS_PreIT.{pdf,png,tiff}                   5 x 2, pre-IT
#   Fig_Supp_RCS_PreIT_ChangeRate.{pdf,png,tiff}        5 x 2, pre-IT change
#
# Survival endpoint: Landmark OS (42-day landmark), which is the IT-appropriate
# primary analysis. Total-OS versions are also saved under total_os/.
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(rms)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(grid)
})

# -- Paths --------------------------------------------------------------------
args_all <- commandArgs(trailingOnly = FALSE)
fa <- args_all[grepl("^--file=", args_all)]
SCRIPT_DIR <- if (length(fa)) {
  dirname(normalizePath(sub("^--file=", "", fa[1]), winslash = "/"))
} else {
  getwd()
}
DATA_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/")
OUT_ROOT <- file.path(SCRIPT_DIR, "psm")
OUT_LM   <- file.path(OUT_ROOT, "landmark")
OUT_TOT  <- file.path(OUT_ROOT, "total_os")
dir.create(OUT_LM,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TOT, showWarnings = FALSE, recursive = TRUE)

# -- Analysis knobs -----------------------------------------------------------
RMS_RCS_NK <- suppressWarnings(as.integer(Sys.getenv("RMS_RCS_NK", "3")))
if (is.na(RMS_RCS_NK) || RMS_RCS_NK < 3L) RMS_RCS_NK <- 3L
N_BOOT <- suppressWarnings(as.integer(Sys.getenv("RMS_RCS_N_BOOT", "300")))
if (is.na(N_BOOT) || N_BOOT < 2L) N_BOOT <- 300L
MIN_N  <- 40L
LANDMARK_MONTHS <- 42 / 30.44

cat("SCRIPT_DIR:", SCRIPT_DIR,
    "\nDATA_DIR:", DATA_DIR,
    "\nRMS_RCS_NK:", RMS_RCS_NK, "| N_BOOT:", N_BOOT, "\n\n")

# -- Pairs of interest --------------------------------------------------------
PAIRS <- list(
  list(id = "02", file = "cohort_02_HAIC_alone_vs_HAIC_then_I.csv",
       group1 = "HAIC_alone", group2 = "HAIC_then_I",
       label  = expression(bold("HAIC" %->% "I  vs  HAIC alone"))),
  list(id = "06", file = "cohort_06_HAIC_alone_vs_HAIC_then_IT.csv",
       group1 = "HAIC_alone", group2 = "HAIC_then_I+T",
       label  = expression(bold("HAIC" %->% "I+T  vs  HAIC alone")))
)

# -- Matrix configuration: 5 biomarkers x 3 timepoints -----------------------
INDICATORS <- list(
  list(row = "AFP",   cols = list(
    list(col = "afp",               type = "static",  log = TRUE),
    list(col = "afp_pre_it",        type = "static",  log = TRUE),
    list(col = "afp_change_pre_it", type = "dynamic", log = FALSE))),
  list(row = "PIVKA", cols = list(
    list(col = "pivka",               type = "static",  log = TRUE),
    list(col = "pivka_pre_it",        type = "static",  log = TRUE),
    list(col = "pivka_change_pre_it", type = "dynamic", log = FALSE))),
  list(row = "SII",   cols = list(
    list(col = "sii_bl",              type = "static",  log = TRUE),
    list(col = "sii_pre_it",          type = "static",  log = TRUE),
    list(col = "sii_change_pre_it",   type = "dynamic", log = FALSE))),
  list(row = "PLR",   cols = list(
    list(col = "plr_bl",              type = "static",  log = FALSE),
    list(col = "plr_pre_it",          type = "static",  log = FALSE),
    list(col = "plr_change_pre_it",   type = "dynamic", log = FALSE))),
  list(row = "NLR",   cols = list(
    list(col = "nlr_bl",              type = "static",  log = FALSE),
    list(col = "nlr_pre_it",          type = "static",  log = FALSE),
    list(col = "nlr_change_pre_it",   type = "dynamic", log = FALSE)))
)

TIMEPOINT_HEADERS <- c("Baseline", "Pre-IT", "Pre-IT Change Rate (%)")

# -- Helpers: anova p, curve prediction ---------------------------------------
parse_anova_p <- function(cell) {
  x <- as.character(cell)[1]
  if (is.na(x) || !nzchar(x)) return(NA_real_)
  if (grepl("^<", x)) return(suppressWarnings(as.numeric(sub("^<\\.?", "", x))))
  suppressWarnings(as.numeric(x))
}

extract_rms_anova_p <- function(fit) {
  a  <- anova(fit); rn <- rownames(a)
  pcol <- if ("P" %in% colnames(a)) "P" else ncol(a)
  int_p    <- NA_real_; nonlin_p <- NA_real_
  idx_int <- grep("^trt \\* [^\\(]+  \\(Factor\\+Higher Order Factors\\)$", rn)
  if (length(idx_int)) int_p <- parse_anova_p(a[idx_int[1], pcol])
  idx_nl <- grep("Nonlinear Interaction : f(A,B) vs. AB", rn, fixed = TRUE)
  if (length(idx_nl)) nonlin_p <- parse_anova_p(a[idx_nl[1], pcol])
  list(int_p = int_p, nonlin_p = nonlin_p)
}

fit_boot_cph_rms <- function(df_b, surv_time_col, nk) {
  dd_b <- suppressWarnings(datadist(df_b[, c("trt", "rcsx"), drop = FALSE]))
  options(datadist = dd_b)
  fml <- as.formula(paste0(
    "Surv(", surv_time_col, ", death_status) ~ trt * rcs(rcsx, ", nk, ")"))
  cph(fml, data = df_b, x = TRUE, y = TRUE, robust = FALSE)
}

predict_hr_curve <- function(fit, x_grid) {
  nd1 <- data.frame(trt = 1, rcsx = x_grid)
  nd0 <- data.frame(trt = 0, rcsx = x_grid)
  exp(as.numeric(predict(fit, nd1, type = "lp")) -
      as.numeric(predict(fit, nd0, type = "lp")))
}

find_crossings <- function(df_plot, hr_ref) {
  cps <- numeric(0)
  for (i in seq_len(nrow(df_plot) - 1)) {
    y1 <- df_plot$hr[i] - hr_ref; y2 <- df_plot$hr[i + 1] - hr_ref
    if (y1 * y2 < 0) {
      frac <- (hr_ref - df_plot$hr[i]) / (df_plot$hr[i + 1] - df_plot$hr[i])
      cps  <- c(cps, df_plot$x[i] + frac * (df_plot$x[i + 1] - df_plot$x[i]))
    }
  }
  cps
}

# -- Derive missing columns ----------------------------------------------------
derive_cols <- function(df) {
  if (!"sii_bl" %in% names(df) && all(c("plt", "nlr_bl") %in% names(df)))
    df$sii_bl <- df$plt * df$nlr_bl
  if (!"nlr_pre_it" %in% names(df) && all(c("neut_pre_it", "lymph_pre_it") %in% names(df)))
    df$nlr_pre_it <- df$neut_pre_it / ifelse(df$lymph_pre_it == 0, NA, df$lymph_pre_it)
  if (!"plr_pre_it" %in% names(df) && all(c("plt_pre_it", "lymph_pre_it") %in% names(df)))
    df$plr_pre_it <- df$plt_pre_it / ifelse(df$lymph_pre_it == 0, NA, df$lymph_pre_it)
  if (!"sii_pre_it" %in% names(df) && all(c("plt_pre_it", "neut_pre_it", "lymph_pre_it") %in% names(df)))
    df$sii_pre_it <- df$plt_pre_it * df$neut_pre_it / ifelse(df$lymph_pre_it == 0, NA, df$lymph_pre_it)
  map <- list(
    list(pre = "afp_pre_it",   bl = "afp",     out = "afp_change_pre_it"),
    list(pre = "pivka_pre_it", bl = "pivka",   out = "pivka_change_pre_it"),
    list(pre = "sii_pre_it",   bl = "sii_bl",  out = "sii_change_pre_it"),
    list(pre = "plr_pre_it",   bl = "plr_bl",  out = "plr_change_pre_it"),
    list(pre = "nlr_pre_it",   bl = "nlr_bl",  out = "nlr_change_pre_it"))
  for (m in map) {
    if (!(m$out %in% names(df)) && all(c(m$pre, m$bl) %in% names(df))) {
      bl_v <- as.numeric(df[[m$bl]]); pv <- as.numeric(df[[m$pre]])
      df[[m$out]] <- ifelse(!is.na(bl_v) & !is.na(pv) & bl_v != 0,
                            (pv - bl_v) / abs(bl_v) * 100, NA)
    }
  }
  if (!"os_lm" %in% names(df))
    df$os_lm <- df$os_months - LANDMARK_MONTHS
  df
}

# -- Single panel (publication styling) ----------------------------------------
make_panel <- function(df, var_col, var_x_label, surv_time_col, plot_type,
                       log_transform = FALSE, show_y_axis = TRUE,
                       show_x_title = TRUE) {
  df <- df[!is.na(df[[var_col]]), ]
  n  <- nrow(df)
  if (n < MIN_N) {
    return(ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0("Insufficient data\n(n=", n, ")"),
               size = 3, color = "grey50") +
      theme_void())
  }

  if (plot_type == "static" && log_transform) {
    df$rcsx <- log1p(pmax(df[[var_col]], 0))
  } else {
    df$rcsx <- df[[var_col]]
  }
  x_vals <- df$rcsx
  nk <- RMS_RCS_NK
  dd <- suppressWarnings(datadist(df[, c("trt", "rcsx"), drop = FALSE]))
  options(datadist = dd)

  out <- tryCatch({
    fml <- as.formula(paste0("Surv(", surv_time_col,
                             ", death_status) ~ trt * rcs(rcsx, ", nk, ")"))
    fit <- cph(fml, data = df, x = TRUE, y = TRUE, robust = FALSE)
    ap  <- extract_rms_anova_p(fit)

    xlim   <- quantile(df$rcsx, c(0.05, 0.95))
    x_grid <- seq(xlim[1], xlim[2], length.out = 200)
    hr_c   <- predict_hr_curve(fit, x_grid)

    set.seed(42)
    boot <- matrix(NA_real_, N_BOOT, length(x_grid))
    for (b in seq_len(N_BOOT)) {
      idx  <- sample.int(nrow(df), replace = TRUE)
      tryCatch({
        fit_b <- fit_boot_cph_rms(df[idx, ], surv_time_col, nk)
        boot[b, ] <- predict_hr_curve(fit_b, x_grid)
      }, error = function(e) NULL)
    }
    hr_lo <- apply(boot, 2, quantile, 0.025, na.rm = TRUE)
    hr_hi <- apply(boot, 2, quantile, 0.975, na.rm = TRUE)
    eps <- 0.05
    hr_lo[is.na(hr_lo)] <- eps
    hr_hi[is.na(hr_hi)] <- pmax(10, hr_c, na.rm = TRUE)
    hr_lo <- pmax(hr_lo, eps)
    hr_hi <- pmax(hr_hi, hr_lo * 1.001)
    hr_c  <- pmax(hr_c,  eps)

    HR_REFS    <- c(1.0, 0.85, 0.7)
    REF_COLORS <- c("#333333", "#0072B2", "#DC0000")
    REF_LT     <- c("dashed", "dotdash", "dotdash")
    REF_LW     <- c(0.7, 0.55, 0.55)

    plot_df    <- data.frame(x = x_grid, hr = hr_c, hr_lo = hr_lo, hr_hi = hr_hi)
    cross_list <- lapply(HR_REFS, function(r) find_crossings(plot_df, r))

    ok_rug <- is.finite(x_vals) & x_vals >= xlim[1] & x_vals <= xlim[2]
    rug_x  <- x_vals[ok_rug]

    p <- ggplot(plot_df, aes(x = x)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.25, ymax = 1,
               fill = "#E8F5E9", alpha = 0.30) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 2.5,
               fill = "#FFF3E0", alpha = 0.30) +
      geom_ribbon(aes(ymin = hr_lo, ymax = hr_hi),
                  fill = "#3C5488", alpha = 0.20) +
      geom_line(aes(y = hr), color = "#3C5488", linewidth = 0.9)

    for (ri in seq_along(HR_REFS)) {
      p <- p + geom_hline(yintercept = HR_REFS[ri], linetype = REF_LT[ri],
                          color = REF_COLORS[ri], linewidth = REF_LW[ri],
                          alpha = 0.9)
    }

    cross_y_labels <- c(1.45, 0.57, 0.41)
    for (ri in seq_along(HR_REFS)) {
      cps <- cross_list[[ri]]
      if (length(cps) > 0) {
        for (cp in cps) {
          if (plot_type == "static" && log_transform) {
            lab <- sprintf("%.3g", expm1(cp))
          } else if (plot_type == "dynamic") {
            lab <- sprintf("%.0f%%", cp)
          } else {
            lab <- sprintf("%.3g", cp)
          }
          p <- p +
            geom_vline(xintercept = cp, linetype = "dotted",
                       color = REF_COLORS[ri], linewidth = 0.45, alpha = 0.8) +
            annotate("point", x = cp, y = HR_REFS[ri],
                     color = REF_COLORS[ri], size = 1.4, shape = 18) +
            annotate("label", x = cp, y = cross_y_labels[ri], label = lab,
                     color = REF_COLORS[ri], size = 2.5, fontface = "bold",
                     fill = "white", label.size = 0.15,
                     label.padding = unit(0.10, "lines"))
        }
      }
    }

    if (plot_type == "dynamic") {
      p <- p + geom_vline(xintercept = 0, linetype = "dotted",
                          color = "#888888", linewidth = 0.45)
    }

    p_sub <- sprintf("P[int]=%s | P[nl]=%s | n=%d",
       ifelse(is.na(ap$int_p),    "NA", sprintf("%.3f", ap$int_p)),
       ifelse(is.na(ap$nonlin_p), "NA", sprintf("%.3f", ap$nonlin_p)), n)

    p <- p +
      geom_rug(data = data.frame(xr = rug_x), aes(x = xr),
               sides = "b", alpha = 0.3, color = "#333333",
               length = unit(1.8, "mm"), linewidth = 0.3,
               inherit.aes = FALSE) +
      labs(title = NULL, subtitle = p_sub,
           x = if (show_x_title) var_x_label else NULL,
           y = NULL) +
      scale_y_log10(breaks = c(0.3, 0.5, 0.7, 1.0, 1.5, 2.0),
                    labels = c("0.3", "0.5", "0.7", "1.0", "1.5", "2.0")) +
      coord_cartesian(ylim = c(0.25, 2.5), xlim = xlim) +
      theme_bw(base_size = 9, base_family = "sans") +
      theme(
        plot.subtitle = element_text(size = 7, hjust = 0.5, color = "grey30",
                                     margin = margin(t = 1, b = 2)),
        axis.text     = element_text(size = 7, color = "black"),
        axis.title.x  = element_text(size = 8, margin = margin(t = 2)),
        axis.title.y  = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
        panel.border  = element_rect(color = "black", linewidth = 0.4),
        plot.margin   = margin(3, 4, 2, 3)
      )

    if (!show_y_axis) {
      p <- p + theme(axis.text.y = element_blank(),
                     axis.ticks.y = element_blank())
    }

    if (plot_type == "static" && log_transform) {
      raw_range <- expm1(xlim)
      cand <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000,
                10000, 20000, 50000, 100000, 200000, 500000)
      rt <- cand[cand >= raw_range[1] & cand <= raw_range[2]]
      if (length(rt) < 3) {
        rt <- pretty(raw_range, n = 4); rt <- rt[rt > 0]
      }
      if (length(rt) > 5)
        rt <- rt[seq(1, length(rt), length.out = 5)]
      fmt <- function(v) ifelse(v >= 1000, sprintf("%gK", v/1000),
                          ifelse(v >= 1, sprintf("%g", v), sprintf("%.1f", v)))
      p <- p + scale_x_continuous(breaks = log1p(rt), labels = fmt(rt),
                                   limits = xlim,
                                   expand = expansion(mult = c(0.02, 0.02)))
    } else {
      p <- p + scale_x_continuous(limits = xlim,
                                   expand = expansion(mult = c(0.03, 0.03)))
    }
    p
  }, error = function(e) {
    cat("    ERROR [", var_col, "]:", conditionMessage(e), "\n")
    ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Model failed",
                        size = 3, color = "red") + theme_void()
  })

  out
}

# -- Text grobs ---------------------------------------------------------------
txt_grob <- function(label, size = 10, face = "plain", col = "black",
                     rot = 0, fill = NA, pad = 0.2) {
  g <- textGrob(label, rot = rot,
                gp = gpar(fontsize = size, fontface = face, col = col,
                          fontfamily = "sans"))
  if (!is.na(fill)) {
    bg <- rectGrob(gp = gpar(fill = fill, col = NA))
    g  <- grobTree(bg, g)
  }
  g
}

panel_label_grob <- function(label) {
  grobTree(
    rectGrob(gp = gpar(fill = "white", col = NA)),
    textGrob(label, x = 0.02, y = 0.5, just = c("left", "center"),
             gp = gpar(fontsize = 13, fontface = "bold", col = "black",
                       fontfamily = "sans"))
  )
}

# -- Build 5 x 3 block for one pair --------------------------------------------
build_block <- function(df, surv_time_col, show_x_titles = TRUE) {
  plots <- vector("list", length(INDICATORS) * 3)
  for (ri in seq_along(INDICATORS)) {
    ind <- INDICATORS[[ri]]
    for (ci in seq_along(ind$cols)) {
      cfg <- ind$cols[[ci]]
      idx <- (ri - 1) * 3 + ci
      xlab <- if (cfg$type == "dynamic") {
        "Change (%)"
      } else if (cfg$log) {
        paste0(ind$row, " (log-scale)")
      } else {
        ind$row
      }
      cat(sprintf("    [%s] %s\n", ind$row, cfg$col))
      if (cfg$col %in% names(df)) {
        plots[[idx]] <- make_panel(
          df, var_col = cfg$col, var_x_label = xlab,
          surv_time_col = surv_time_col,
          plot_type = cfg$type, log_transform = cfg$log,
          show_y_axis = (ci == 1), show_x_title = show_x_titles)
      } else {
        plots[[idx]] <- ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "N/A",
                   size = 3, color = "grey60") + theme_void()
      }
    }
  }
  plots
}

# ===========================================================================
# Figure builders
# ===========================================================================

# ---- Main figure: 5 x 6 (two pairs side by side) --------------------------
build_main_figure <- function(data_list, surv_time_col) {
  blocks <- list(
    build_block(data_list[[1]], surv_time_col, show_x_titles = TRUE),
    build_block(data_list[[2]], surv_time_col, show_x_titles = TRUE)
  )

  # 30 panel plots
  plots <- c()
  for (ri in seq_along(INDICATORS)) {
    for (pi in seq_along(PAIRS)) {
      for (ci in 1:3) {
        idx <- (ri - 1) * 3 + ci
        plots <- c(plots, list(blocks[[pi]][[idx]]))
      }
    }
  }

  # Convert to grobs
  plot_grobs <- lapply(plots, ggplotGrob)

  # Headers
  pair_headers <- lapply(PAIRS, function(p)
    txt_grob(sprintf("%s  vs  %s",
                     gsub("_", " ", p$group2),
                     gsub("_", " ", p$group1)),
             size = 11, face = "bold", col = "white", fill = "#2E4E7E"))

  sub_headers <- lapply(TIMEPOINT_HEADERS, function(h)
    txt_grob(h, size = 10, face = "bold", col = "#1A1A1A"))

  row_labels <- lapply(INDICATORS, function(ind)
    txt_grob(ind$row, size = 12, face = "bold", rot = 90, col = "#1A1A1A"))

  panel_label_A <- panel_label_grob("A")
  panel_label_B <- panel_label_grob("B")

  # Layout:
  # Row 1 : [      ] [ A header pair1                        ] [ B header pair2                        ]
  # Row 2 : [      ] [BL] [PreIT] [PreIT change] [BL] [PreIT] [PreIT change]
  # Row 3 : [AFP  ] [panel] [panel] [panel] [panel] [panel] [panel]
  # ...
  # Row 7 : [NLR  ] [panel] [panel] [panel] [panel] [panel] [panel]
  #
  # Column widths: 0.6cm row-label; then six equal panel columns.
  n_panels_total <- length(plot_grobs)              # 30
  # grob layout indices
  # slots: 1..30 panels, 31..32 pair headers, 33..35 sub-pair1 subheaders,
  # 36..38 sub-pair2 subheaders, 39..43 row labels, 44 panel label A, 45 panel label B
  all_grobs <- c(plot_grobs,
                 pair_headers,
                 sub_headers,      # pair 1
                 sub_headers,      # pair 2
                 row_labels,
                 list(panel_label_A, panel_label_B))

  n_rows <- length(INDICATORS)      # 5
  n_cols <- 6
  # layout matrix: (2 header rows + 5 indicator rows) x (1 label col + 6 panel cols) = 7 x 7
  lm <- matrix(NA_integer_, nrow = 2 + n_rows, ncol = 1 + n_cols)

  # Row 1 (pair headers): cols 2-4 = pair1 header (idx 31), cols 5-7 = pair2 header (idx 32)
  lm[1, 2:4] <- 31L
  lm[1, 5:7] <- 32L
  # Row 2 (sub headers): 3 pair1 subheaders (33,34,35) + 3 pair2 subheaders (36,37,38)
  lm[2, 2:4] <- c(33L, 34L, 35L)
  lm[2, 5:7] <- c(36L, 37L, 38L)
  # Rows 3-7: row label + 6 panels
  for (ri in seq_len(n_rows)) {
    lm[ri + 2, 1] <- 38L + ri    # indices 39..43
    start <- (ri - 1) * n_cols + 1L
    lm[ri + 2, 2:7] <- start:(start + 5L)
  }
  # Panel labels overlaid via annotation grob - we'll place them as additional row
  # Actually, panel-label grobs (A,B) can replace nothing slot-wise; put them above headers.
  # Simpler: prepend a narrow row above row 1 carrying A / B labels at the left of each block
  lm <- rbind(matrix(NA_integer_, nrow = 1, ncol = ncol(lm)), lm)
  lm[1, 2] <- 44L   # "A"
  lm[1, 5] <- 45L   # "B"

  widths  <- unit.c(unit(0.55, "cm"), unit(rep(1, n_cols), "null"))
  heights <- unit.c(
    unit(0.55, "cm"),   # A / B panel labels
    unit(0.70, "cm"),   # pair headers
    unit(0.55, "cm"),   # sub-headers
    unit(rep(1, n_rows), "null")
  )

  arrangeGrob(
    grobs = all_grobs, layout_matrix = lm,
    widths = widths, heights = heights,
    bottom = textGrob(
      sprintf(paste0("Hazard ratio (treatment vs. HAIC alone) from PSM-matched ",
                     "Cox model with restricted cubic splines (knots=%d), ",
                     "bootstrap 95%% CI (B=%d). Landmark: %d days."),
              RMS_RCS_NK, N_BOOT,
              ifelse(surv_time_col == "os_lm", 42L, 0L)),
      gp = gpar(fontsize = 8, col = "grey30", fontfamily = "sans"))
  )
}

# ---- Supplementary figure: 5 x 2 (one timepoint x two pairs) --------------
build_supp_figure <- function(data_list, surv_time_col, tp_index, tp_title) {
  # tp_index: 1 = baseline, 2 = pre-IT, 3 = pre-IT change
  plots <- vector("list", length(INDICATORS) * 2)
  for (ri in seq_along(INDICATORS)) {
    ind <- INDICATORS[[ri]]
    cfg <- ind$cols[[tp_index]]
    xlab <- if (cfg$type == "dynamic") {
      "Change (%)"
    } else if (cfg$log) {
      paste0(ind$row, " (log-scale)")
    } else {
      ind$row
    }
    for (pi in seq_along(PAIRS)) {
      idx <- (ri - 1) * 2 + pi
      df  <- data_list[[pi]]
      if (cfg$col %in% names(df)) {
        plots[[idx]] <- make_panel(
          df, var_col = cfg$col, var_x_label = xlab,
          surv_time_col = surv_time_col,
          plot_type = cfg$type, log_transform = cfg$log,
          show_y_axis = (pi == 1), show_x_title = TRUE)
      } else {
        plots[[idx]] <- ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "N/A") + theme_void()
      }
    }
  }

  plot_grobs <- lapply(plots, ggplotGrob)

  pair_headers <- lapply(PAIRS, function(p)
    txt_grob(sprintf("%s  vs  %s",
                     gsub("_", " ", p$group2),
                     gsub("_", " ", p$group1)),
             size = 11, face = "bold", col = "white", fill = "#2E4E7E"))

  row_labels <- lapply(INDICATORS, function(ind)
    txt_grob(ind$row, size = 12, face = "bold", rot = 90))

  all_grobs <- c(plot_grobs, pair_headers, row_labels)

  n_rows <- length(INDICATORS)
  n_cols <- 2
  # (1 header row + 5 indicator rows) x (1 label col + 2 panel cols)
  lm <- matrix(NA_integer_, nrow = 1 + n_rows, ncol = 1 + n_cols)
  lm[1, 2] <- length(plot_grobs) + 1L                # pair 1
  lm[1, 3] <- length(plot_grobs) + 2L                # pair 2
  for (ri in seq_len(n_rows)) {
    lm[ri + 1, 1] <- length(plot_grobs) + 2L + ri
    start <- (ri - 1) * n_cols + 1L
    lm[ri + 1, 2:3] <- start:(start + 1L)
  }

  widths  <- unit.c(unit(0.55, "cm"), unit(rep(1, n_cols), "null"))
  heights <- unit.c(unit(0.70, "cm"), unit(rep(1, n_rows), "null"))

  arrangeGrob(
    grobs = all_grobs, layout_matrix = lm,
    widths = widths, heights = heights,
    top = textGrob(tp_title, gp = gpar(fontsize = 13, fontface = "bold",
                                       fontfamily = "sans")),
    bottom = textGrob(
      sprintf(paste0("Hazard ratio (treatment vs. HAIC alone) from PSM-matched ",
                     "Cox model with restricted cubic splines (knots=%d), ",
                     "bootstrap 95%% CI (B=%d)."),
              RMS_RCS_NK, N_BOOT),
      gp = gpar(fontsize = 8, col = "grey30", fontfamily = "sans"))
  )
}

# -- Save helper (PDF + PNG + TIFF) -------------------------------------------
save_fig <- function(grob, out_dir, base_name, width, height) {
  ggsave(file.path(out_dir, paste0(base_name, ".pdf")),
         grob, width = width, height = height, device = "pdf",
         useDingbats = FALSE)
  ggsave(file.path(out_dir, paste0(base_name, ".png")),
         grob, width = width, height = height, dpi = 600, device = "png")
  tryCatch({
    ggsave(file.path(out_dir, paste0(base_name, ".tiff")),
           grob, width = width, height = height, dpi = 600,
           device = "tiff", compression = "lzw")
  }, error = function(e) cat("  (TIFF save skipped:", conditionMessage(e), ")\n"))
  cat("  saved:", base_name, ".pdf/.png/.tiff\n")
}

# ===========================================================================
# MAIN
# ===========================================================================

# 1. Load data for both pairs
data_list <- lapply(PAIRS, function(p) {
  csv_path <- file.path(DATA_DIR, p$file)
  cat("\nLoading:", p$file, "\n")
  if (!file.exists(csv_path)) {
    stop("Missing: ", csv_path, "\nRun build_cohort_psm.py first")
  }
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  df <- derive_cols(df)
  cat("  Rows:", nrow(df),
      "| trt counts:", paste(table(df$trt), collapse = "/"),
      "\n")
  df
})

# 2. Both timescales
for (ts in list(
  list(col = "os_lm",     out = OUT_LM,  label = "Landmark OS (42 days)"),
  list(col = "os_months", out = OUT_TOT, label = "Total OS")
)) {
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat("Timescale:", ts$label, "\n")
  cat(strrep("=", 72), "\n")

  # Filter rows with positive follow-up
  data_list_f <- lapply(data_list, function(d) {
    if (ts$col == "os_lm") d[!is.na(d$os_lm) & d$os_lm > 0, ]
    else                   d[!is.na(d$os_months) & d$os_months > 0, ]
  })

  # ------- Main figure (5 x 6) -------
  cat("\n[Main figure: 5 indicators x 2 cohorts x 3 timepoints]\n")
  main_fig <- build_main_figure(data_list_f, ts$col)
  save_fig(main_fig, ts$out, "Fig_Main_RCS_5indicators",
           width = 17.5, height = 10.0)

  # ------- Supplementary figures (5 x 2, one per timepoint) -------
  supp_specs <- list(
    list(idx = 1, name = "Fig_Supp_RCS_Baseline",
         title = "Baseline biomarker x treatment interaction (RCS)"),
    list(idx = 2, name = "Fig_Supp_RCS_PreIT",
         title = "Pre-IT biomarker x treatment interaction (RCS)"),
    list(idx = 3, name = "Fig_Supp_RCS_PreIT_ChangeRate",
         title = "Pre-IT change rate x treatment interaction (RCS)")
  )
  for (ss in supp_specs) {
    cat("\n[Supp figure:", ss$title, "]\n")
    sfig <- build_supp_figure(data_list_f, ts$col, ss$idx, ss$title)
    save_fig(sfig, ts$out, ss$name, width = 8.2, height = 13.5)
  }
}

cat("\n=== DONE. Outputs at:", SCRIPT_DIR, "===\n")
