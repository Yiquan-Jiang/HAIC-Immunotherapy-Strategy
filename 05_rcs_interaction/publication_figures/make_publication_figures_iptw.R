#!/usr/bin/env Rscript
# =============================================================================
# make_publication_figures_iptw.R
#
# Publication-ready figures for RCS non-linear interaction analysis
# using ROUTE B (composite cohort + IPTW-weighted Cox).
#
# Comparisons (two composite cohorts):
#   THEN_I  : HAIC_then_I   vs HAIC_alone  (trt_haic_then_i)
#   THEN_IT : HAIC_then_I+T vs HAIC_alone  (trt_haic_then_it)
#
# Biomarkers (5): AFP, PIVKA, SII, PLR, NLR
# Timepoints (3): Baseline, Pre-IT, Pre-IT Change Rate
#
# Outputs (in ./iptw/{landmark,total_os}/):
#   Fig_Main_RCS_5indicators.{pdf,png,tiff}              (wider, landscape)
#   Fig_Supp_RCS_Baseline.{pdf,png,tiff}
#   Fig_Supp_RCS_PreIT.{pdf,png,tiff}
#   Fig_Supp_RCS_PreIT_ChangeRate.{pdf,png,tiff}
#
# Primary endpoint: 42-day landmark OS (IT-appropriate).
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(rms)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(glmnet)
  library(grid)
})

# -- Paths --------------------------------------------------------------------
args_all <- commandArgs(trailingOnly = FALSE)
fa <- args_all[grepl("^--file=", args_all)]
SCRIPT_DIR <- if (length(fa)) {
  dirname(normalizePath(sub("^--file=", "", fa[1]), winslash = "/"))
} else getwd()

COMPOSITE_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", "afp_pivka_composite"),
                               winslash = "/")
OUT_ROOT <- file.path(SCRIPT_DIR, "iptw")
OUT_LM   <- file.path(OUT_ROOT, "landmark")
OUT_TOT  <- file.path(OUT_ROOT, "total_os")
dir.create(OUT_LM,  showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_TOT, showWarnings = FALSE, recursive = TRUE)

# -- Analysis knobs -----------------------------------------------------------
RMS_RCS_NK <- suppressWarnings(as.integer(Sys.getenv("RMS_RCS_NK", "3")))
if (is.na(RMS_RCS_NK) || RMS_RCS_NK < 3L) RMS_RCS_NK <- 3L
N_BOOT <- suppressWarnings(as.integer(Sys.getenv("RMS_RCS_N_BOOT", "300")))
if (is.na(N_BOOT) || N_BOOT < 2L) N_BOOT <- 300L
MIN_N  <- 50L
LANDMARK_MONTHS <- 42 / 30.44

cat("SCRIPT_DIR:", SCRIPT_DIR,
    "\nCOMPOSITE_DIR:", COMPOSITE_DIR,
    "\nRMS_RCS_NK:", RMS_RCS_NK, "| N_BOOT:", N_BOOT, "\n\n")

# -- Cohorts of interest (Route B: composite + IPTW) --------------------------
COHORTS <- list(
  list(id = "THEN_I",  csv = "composite_THEN_I_cohort.csv",
       trt_col = "trt_haic_then_i",
       trt_label = "HAIC_then_I",   ctrl_label = "HAIC_alone"),
  list(id = "THEN_IT", csv = "composite_THEN_IT_cohort.csv",
       trt_col = "trt_haic_then_it",
       trt_label = "HAIC_then_I+T", ctrl_label = "HAIC_alone")
)

# -- 5 biomarkers x 3 timepoints ---------------------------------------------
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

# -- Propensity score variables (ridge logistic) ------------------------------
PS_VARS_STATIC <- c("albi_bl", "inr", "plt",
                    "tumor_max_diameter_cm", "tumor_count_enc",
                    "pvtt_grade", "hvtt_binary", "ivc_ra_binary", "ascites_score_enc",
                    "log_afp_bl", "log_pivka_bl",
                    "metastasis_binary", "lymph_node_binary",
                    "neut_bl", "lymph_bl", "mono_bl")
PS_VARS_DYN_EXTRA <- c("afp_pre_it", "pivka_pre_it",
                       "albi_pre_it", "neut_pre_it", "ast_pre_it", "alt_pre_it",
                       "afp_change_pre_it", "pivka_change_pre_it")

# -- IPTW: ridge logistic PS, stabilized weights, clipped to [0.05, 0.95] -----
compute_iptw <- function(df_in, trt_col, ps_vars_base) {
  ps_use <- ps_vars_base[ps_vars_base %in% names(df_in)]
  ps_use <- ps_use[sapply(ps_use,
                          function(v) sd(df_in[[v]], na.rm = TRUE) > 0)]
  df_cc <- df_in[complete.cases(df_in[, ps_use, drop = FALSE]), ]
  if (nrow(df_cc) < 20) return(df_cc[0, , drop = FALSE])
  X_sc <- scale(as.matrix(df_cc[, ps_use, drop = FALSE]))
  y    <- df_cc[[trt_col]]
  cv_fit  <- cv.glmnet(X_sc, y, family = "binomial", alpha = 0, nfolds = 5)
  ps_prob <- as.numeric(predict(cv_fit, newx = X_sc,
                                s = "lambda.min", type = "response"))
  ps_prob <- pmin(pmax(ps_prob, 0.05), 0.95)
  p_treat <- mean(y)
  df_cc$sw <- ifelse(df_cc[[trt_col]] == 1,
                     p_treat / ps_prob, (1 - p_treat) / (1 - ps_prob))
  df_cc
}

# -- Helpers: anova p, curve prediction ---------------------------------------
parse_anova_p <- function(cell) {
  x <- as.character(cell)[1]
  if (is.na(x) || !nzchar(x)) return(NA_real_)
  if (grepl("^<", x)) return(suppressWarnings(as.numeric(sub("^<\\.?", "", x))))
  suppressWarnings(as.numeric(x))
}

extract_rms_anova_p <- function(fit, trt_col) {
  a  <- anova(fit); rn <- rownames(a)
  pcol <- if ("P" %in% colnames(a)) "P" else ncol(a)
  int_p    <- NA_real_; nonlin_p <- NA_real_
  idx_int <- grep(paste0("^", trt_col, " \\* [^\\(]+  \\(Factor\\+Higher Order Factors\\)$"), rn)
  if (length(idx_int)) int_p <- parse_anova_p(a[idx_int[1], pcol])
  idx_nl <- grep("Nonlinear Interaction : f(A,B) vs. AB", rn, fixed = TRUE)
  if (length(idx_nl)) nonlin_p <- parse_anova_p(a[idx_nl[1], pcol])
  list(int_p = int_p, nonlin_p = nonlin_p)
}

fit_boot_cph_iptw <- function(df_b, surv_time_col, nk, trt_col) {
  dd_b <- suppressWarnings(datadist(df_b[, c(trt_col, "rcsx"), drop = FALSE]))
  options(datadist = dd_b)
  fml <- as.formula(paste0(
    "Surv(", surv_time_col, ", death_status) ~ ", trt_col,
    " * rcs(rcsx, ", nk, ")"))
  cph(fml, data = df_b, weights = sw, x = TRUE, y = TRUE, robust = FALSE)
}

predict_hr_curve <- function(fit, x_grid, trt_col) {
  nd1 <- setNames(data.frame(1L, x_grid), c(trt_col, "rcsx"))
  nd0 <- setNames(data.frame(0L, x_grid), c(trt_col, "rcsx"))
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

# -- Derive missing columns ---------------------------------------------------
derive_cols <- function(df) {
  if (!"sii_bl" %in% names(df) && all(c("plt", "nlr_bl") %in% names(df)))
    df$sii_bl <- df$plt * df$nlr_bl
  if (!"nlr_pre_it" %in% names(df) && all(c("neut_pre_it", "lymph_pre_it") %in% names(df)))
    df$nlr_pre_it <- df$neut_pre_it / ifelse(df$lymph_pre_it == 0, NA, df$lymph_pre_it)
  if (!"plr_pre_it" %in% names(df) && all(c("plt_pre_it", "lymph_pre_it") %in% names(df)))
    df$plr_pre_it <- df$plt_pre_it / ifelse(df$lymph_pre_it == 0, NA, df$lymph_pre_it)
  if (!"sii_pre_it" %in% names(df) && all(c("plt_pre_it", "neut_pre_it", "lymph_pre_it") %in% names(df)))
    df$sii_pre_it <- df$plt_pre_it * df$neut_pre_it / ifelse(df$lymph_pre_it == 0, NA, df$lymph_pre_it)
  if (!"albi_pre_it" %in% names(df) && all(c("tbil_pre_it", "alb_pre_it") %in% names(df)))
    df$albi_pre_it <- ifelse(df$tbil_pre_it > 0 & !is.na(df$tbil_pre_it) & !is.na(df$alb_pre_it),
                             0.66 * log10(df$tbil_pre_it) - 0.085 * df$alb_pre_it, NA)

  change_map <- list(
    list(pre = "afp_pre_it",   bl = "afp",     out = "afp_change_pre_it"),
    list(pre = "pivka_pre_it", bl = "pivka",   out = "pivka_change_pre_it"),
    list(pre = "sii_pre_it",   bl = "sii_bl",  out = "sii_change_pre_it"),
    list(pre = "plr_pre_it",   bl = "plr_bl",  out = "plr_change_pre_it"),
    list(pre = "nlr_pre_it",   bl = "nlr_bl",  out = "nlr_change_pre_it"))
  for (m in change_map) {
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

# -- Single panel (IPTW, publication styling) ---------------------------------
make_panel <- function(df, var_col, var_x_label, surv_time_col, plot_type,
                       log_transform = FALSE, trt_col,
                       ps_vars_base, show_y_axis = TRUE, show_x_title = TRUE) {
  df_w <- compute_iptw(df, trt_col, ps_vars_base)
  df_w <- df_w[!is.na(df_w[[var_col]]), ]
  n <- nrow(df_w)
  if (n < MIN_N) {
    return(ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0("Insufficient data\n(n=", n, ")"),
               size = 3, color = "grey50") + theme_void())
  }

  if (plot_type == "static" && log_transform) {
    df_w$rcsx <- log1p(pmax(df_w[[var_col]], 0))
  } else {
    df_w$rcsx <- df_w[[var_col]]
  }
  x_vals <- df_w$rcsx
  nk <- RMS_RCS_NK
  dd <- suppressWarnings(datadist(df_w[, c(trt_col, "rcsx"), drop = FALSE]))
  options(datadist = dd)

  out <- tryCatch({
    fml <- as.formula(paste0("Surv(", surv_time_col,
                             ", death_status) ~ ", trt_col,
                             " * rcs(rcsx, ", nk, ")"))
    fit <- cph(fml, data = df_w, weights = sw,
               x = TRUE, y = TRUE, robust = FALSE)
    ap  <- extract_rms_anova_p(fit, trt_col)

    xlim   <- quantile(df_w$rcsx, c(0.05, 0.95))
    x_grid <- seq(xlim[1], xlim[2], length.out = 200)
    hr_c   <- predict_hr_curve(fit, x_grid, trt_col)

    set.seed(42)
    boot <- matrix(NA_real_, N_BOOT, length(x_grid))
    for (b in seq_len(N_BOOT)) {
      idx  <- sample.int(nrow(df_w), replace = TRUE)
      tryCatch({
        fit_b <- fit_boot_cph_iptw(df_w[idx, ], surv_time_col, nk, trt_col)
        boot[b, ] <- predict_hr_curve(fit_b, x_grid, trt_col)
      }, error = function(e) NULL)
    }
    hr_lo <- apply(boot, 2, quantile, 0.025, na.rm = TRUE)
    hr_hi <- apply(boot, 2, quantile, 0.975, na.rm = TRUE)
    eps <- 0.05
    hr_lo[is.na(hr_lo)] <- eps
    hr_hi[is.na(hr_hi)] <- pmax(10, hr_c, na.rm = TRUE)
    hr_lo <- pmax(hr_lo, eps); hr_hi <- pmax(hr_hi, hr_lo * 1.001)
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
           x = if (show_x_title) var_x_label else NULL, y = NULL) +
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
      if (length(rt) < 3) { rt <- pretty(raw_range, n = 4); rt <- rt[rt > 0] }
      if (length(rt) > 5) rt <- rt[seq(1, length(rt), length.out = 5)]
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
                     rot = 0, fill = NA) {
  g <- textGrob(label, rot = rot,
                gp = gpar(fontsize = size, fontface = face, col = col,
                          fontfamily = "sans"))
  if (!is.na(fill)) {
    g <- grobTree(rectGrob(gp = gpar(fill = fill, col = NA)), g)
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

# -- Build a 5 x 3 block for one cohort ---------------------------------------
build_block <- function(df, surv_time_col, trt_col, show_x_titles = TRUE) {
  plots <- vector("list", length(INDICATORS) * 3)
  for (ri in seq_along(INDICATORS)) {
    ind <- INDICATORS[[ri]]
    for (ci in seq_along(ind$cols)) {
      cfg <- ind$cols[[ci]]
      idx <- (ri - 1) * 3 + ci
      ps_base <- if (cfg$type == "dynamic")
                   c(PS_VARS_STATIC, PS_VARS_DYN_EXTRA)
                 else PS_VARS_STATIC
      xlab <- if (cfg$type == "dynamic") "Change (%)" else
              if (cfg$log) paste0(ind$row, " (log-scale)") else ind$row
      cat(sprintf("    [%s] %s\n", ind$row, cfg$col))
      if (cfg$col %in% names(df)) {
        plots[[idx]] <- make_panel(
          df, var_col = cfg$col, var_x_label = xlab,
          surv_time_col = surv_time_col, plot_type = cfg$type,
          log_transform = cfg$log, trt_col = trt_col,
          ps_vars_base = ps_base,
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

# -- Main figure: 5 x 6 (two cohorts side-by-side), landscape -----------------
build_main_figure <- function(data_list, surv_time_col) {
  blocks <- lapply(seq_along(COHORTS), function(i)
    build_block(data_list[[i]], surv_time_col, COHORTS[[i]]$trt_col,
                show_x_titles = TRUE))

  # Order: AFP row -> cohort1 c1-c3, cohort2 c1-c3 ; PIVKA row ... etc.
  plots <- list()
  for (ri in seq_along(INDICATORS)) {
    for (pi in seq_along(COHORTS)) {
      for (ci in 1:3) {
        plots <- c(plots, list(blocks[[pi]][[(ri - 1) * 3 + ci]]))
      }
    }
  }
  plot_grobs <- lapply(plots, ggplotGrob)

  cohort_headers <- lapply(COHORTS, function(c)
    txt_grob(sprintf("%s  vs  %s",
                     gsub("_", " ", c$trt_label),
                     gsub("_", " ", c$ctrl_label)),
             size = 11, face = "bold", col = "white", fill = "#2E4E7E"))
  sub_headers <- lapply(TIMEPOINT_HEADERS, function(h)
    txt_grob(h, size = 10, face = "bold", col = "#1A1A1A"))
  row_labels <- lapply(INDICATORS, function(ind)
    txt_grob(ind$row, size = 12, face = "bold", rot = 90, col = "#1A1A1A"))
  panel_A <- panel_label_grob("A"); panel_B <- panel_label_grob("B")

  n_rows <- length(INDICATORS); n_cols <- 6
  all_grobs <- c(plot_grobs, cohort_headers,
                 sub_headers, sub_headers,
                 row_labels, list(panel_A, panel_B))

  lm <- matrix(NA_integer_, nrow = 2 + n_rows, ncol = 1 + n_cols)
  lm[1, 2:4] <- 31L; lm[1, 5:7] <- 32L
  lm[2, 2:4] <- c(33L, 34L, 35L); lm[2, 5:7] <- c(36L, 37L, 38L)
  for (ri in seq_len(n_rows)) {
    lm[ri + 2, 1] <- 38L + ri
    s <- (ri - 1) * n_cols + 1L
    lm[ri + 2, 2:7] <- s:(s + 5L)
  }
  lm <- rbind(matrix(NA_integer_, nrow = 1, ncol = ncol(lm)), lm)
  lm[1, 2] <- 44L; lm[1, 5] <- 45L

  widths  <- unit.c(unit(0.55, "cm"), unit(rep(1, n_cols), "null"))
  heights <- unit.c(
    unit(0.55, "cm"), unit(0.70, "cm"), unit(0.55, "cm"),
    unit(rep(1, n_rows), "null"))

  arrangeGrob(
    grobs = all_grobs, layout_matrix = lm,
    widths = widths, heights = heights,
    bottom = textGrob(
      sprintf(paste0("Hazard ratio (treatment vs. HAIC alone) from ",
                     "IPTW-weighted Cox model with restricted cubic splines ",
                     "(knots=%d); propensity score from ridge logistic ",
                     "regression, stabilized weights clipped to [0.05, 0.95]; ",
                     "bootstrap 95%% CI (B=%d)."),
              RMS_RCS_NK, N_BOOT),
      gp = gpar(fontsize = 8, col = "grey30", fontfamily = "sans")))
}

# -- Supplementary figure: 5 x 2 (one timepoint x two cohorts) ---------------
build_supp_figure <- function(data_list, surv_time_col, tp_index, tp_title) {
  plots <- vector("list", length(INDICATORS) * 2)
  for (ri in seq_along(INDICATORS)) {
    ind <- INDICATORS[[ri]]
    cfg <- ind$cols[[tp_index]]
    ps_base <- if (cfg$type == "dynamic")
                 c(PS_VARS_STATIC, PS_VARS_DYN_EXTRA)
               else PS_VARS_STATIC
    xlab <- if (cfg$type == "dynamic") "Change (%)" else
            if (cfg$log) paste0(ind$row, " (log-scale)") else ind$row
    for (pi in seq_along(COHORTS)) {
      idx <- (ri - 1) * 2 + pi
      df  <- data_list[[pi]]
      trt_col <- COHORTS[[pi]]$trt_col
      cat(sprintf("    [%s | %s] %s\n", ind$row, COHORTS[[pi]]$id, cfg$col))
      if (cfg$col %in% names(df)) {
        plots[[idx]] <- make_panel(
          df, var_col = cfg$col, var_x_label = xlab,
          surv_time_col = surv_time_col, plot_type = cfg$type,
          log_transform = cfg$log, trt_col = trt_col,
          ps_vars_base = ps_base,
          show_y_axis = (pi == 1), show_x_title = TRUE)
      } else {
        plots[[idx]] <- ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "N/A") + theme_void()
      }
    }
  }
  plot_grobs <- lapply(plots, ggplotGrob)

  cohort_headers <- lapply(COHORTS, function(c)
    txt_grob(sprintf("%s  vs  %s",
                     gsub("_", " ", c$trt_label),
                     gsub("_", " ", c$ctrl_label)),
             size = 11, face = "bold", col = "white", fill = "#2E4E7E"))
  row_labels <- lapply(INDICATORS, function(ind)
    txt_grob(ind$row, size = 12, face = "bold", rot = 90))

  all_grobs <- c(plot_grobs, cohort_headers, row_labels)
  n_rows <- length(INDICATORS); n_cols <- 2
  lm <- matrix(NA_integer_, nrow = 1 + n_rows, ncol = 1 + n_cols)
  lm[1, 2] <- length(plot_grobs) + 1L
  lm[1, 3] <- length(plot_grobs) + 2L
  for (ri in seq_len(n_rows)) {
    lm[ri + 1, 1] <- length(plot_grobs) + 2L + ri
    s <- (ri - 1) * n_cols + 1L
    lm[ri + 1, 2:3] <- s:(s + 1L)
  }
  widths  <- unit.c(unit(0.55, "cm"), unit(rep(1, n_cols), "null"))
  heights <- unit.c(unit(0.70, "cm"), unit(rep(1, n_rows), "null"))

  arrangeGrob(
    grobs = all_grobs, layout_matrix = lm,
    widths = widths, heights = heights,
    top = textGrob(tp_title, gp = gpar(fontsize = 13, fontface = "bold",
                                       fontfamily = "sans")),
    bottom = textGrob(
      sprintf(paste0("Hazard ratio from IPTW-weighted Cox model with ",
                     "restricted cubic splines (knots=%d); bootstrap 95%% CI ",
                     "(B=%d)."), RMS_RCS_NK, N_BOOT),
      gp = gpar(fontsize = 8, col = "grey30", fontfamily = "sans")))
}

# -- Save helper --------------------------------------------------------------
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
  }, error = function(e) cat("  (TIFF skipped:", conditionMessage(e), ")\n"))
  cat("  saved:", base_name, ".pdf/.png/.tiff\n")
}

# ===========================================================================
# MAIN
# ===========================================================================
data_list <- lapply(COHORTS, function(c) {
  csv_path <- file.path(COMPOSITE_DIR, c$csv)
  cat("\nLoading:", c$csv, "\n")
  if (!file.exists(csv_path))
    stop("Missing: ", csv_path,
         "\nRun afp_pivka_composite/00_build_composite_cohorts.py first")
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  df <- derive_cols(df)
  cat("  Rows:", nrow(df),
      "| ", c$trt_col, " counts:",
      paste(table(df[[c$trt_col]]), collapse = "/"), "\n")
  df
})

for (ts in list(
  list(col = "os_lm",     out = OUT_LM,  label = "Landmark OS (42 days)"),
  list(col = "os_months", out = OUT_TOT, label = "Total OS")
)) {
  cat("\n", strrep("=", 72), "\n", sep = "")
  cat("Timescale:", ts$label, "\n")
  cat(strrep("=", 72), "\n")

  data_list_f <- lapply(data_list, function(d) {
    if (ts$col == "os_lm") d[!is.na(d$os_lm) & d$os_lm > 0, ]
    else                   d[!is.na(d$os_months) & d$os_months > 0, ]
  })

  # ------ Main (wider, landscape) ------
  cat("\n[Main figure: 5 indicators x 2 cohorts x 3 timepoints]\n")
  main_fig <- build_main_figure(data_list_f, ts$col)
  save_fig(main_fig, ts$out, "Fig_Main_RCS_5indicators",
           width = 17.5, height = 10.0)

  # ------ Supplementary (3 figures) ------
  for (ss in list(
    list(idx = 1, name = "Fig_Supp_RCS_Baseline",
         title = "Baseline biomarker x treatment interaction (RCS, IPTW)"),
    list(idx = 2, name = "Fig_Supp_RCS_PreIT",
         title = "Pre-IT biomarker x treatment interaction (RCS, IPTW)"),
    list(idx = 3, name = "Fig_Supp_RCS_PreIT_ChangeRate",
         title = "Pre-IT change rate x treatment interaction (RCS, IPTW)")
  )) {
    cat("\n[Supp figure:", ss$title, "]\n")
    sfig <- build_supp_figure(data_list_f, ts$col, ss$idx, ss$title)
    save_fig(sfig, ts$out, ss$name, width = 8.2, height = 13.5)
  }
}

cat("\n=== DONE. IPTW outputs at:", OUT_ROOT, "===\n")
