#!/usr/bin/env Rscript
# =============================================================================
# 02_rcs_matrix_panel.R
#
# 8×5 Matrix Panel (IPTW weighted):
#   Rows: AFP, PIVKA, PIV, SII, NLR, PLR, MONOCYTE, ALBI
#   Cols: Baseline, Pre-HAIC-3, Pre-IT, Pre-HAIC-3 Change Rate, Pre-IT Change Rate
#
# 对每个 cohort × 每种时间尺度(landmark / total_os) 各输出一张大图
#
# 用法: Rscript 02_rcs_matrix_panel.R ALL
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

SCRIPT_DIR <- Sys.getenv("AFP_PIVKA_SCRIPT_DIR", unset = "")
if (!nzchar(SCRIPT_DIR)) {
  tryCatch({
    args_all <- commandArgs(trailingOnly = FALSE)
    fa <- args_all[grepl("^--file=", args_all)]
    SCRIPT_DIR <- if (length(fa)) {
      dirname(normalizePath(sub("^--file=", "", fa[1]), winslash = "/"))
    } else {
      getwd()
    }
  }, error = function(e) {
    SCRIPT_DIR <<- getwd()
  })
}

cohort_key <- Sys.getenv("AFP_PIVKA_COHORT", unset = "")
if (!nzchar(cohort_key)) {
  tryCatch({
    ck <- commandArgs(trailingOnly = TRUE)
    cohort_key <- if (length(ck) > 0) toupper(ck[1]) else "ALL"
  }, error = function(e) {
    cohort_key <<- "ALL"
  })
} else {
  cohort_key <- toupper(cohort_key)
}

COHORT_REGISTRY <- list(
  THEN_IT = list(
    csv = "composite_THEN_IT_cohort.csv",
    trt_col = "trt_haic_then_it",
    trt_label = "HAIC_then_I+T",
    ctrl_label = "HAIC_alone"
  ),
  THEN_I = list(
    csv = "composite_THEN_I_cohort.csv",
    trt_col = "trt_haic_then_i",
    trt_label = "HAIC_then_I",
    ctrl_label = "HAIC_alone"
  ),
  THEN_T = list(
    csv = "composite_THEN_T_cohort.csv",
    trt_col = "trt_haic_then_t",
    trt_label = "HAIC_then_T",
    ctrl_label = "HAIC_alone"
  ),
  T_CONC = list(
    csv = "composite_T_CONC_cohort.csv",
    trt_col = "trt_haic_t_conc",
    trt_label = "HAIC+T_concurrent",
    ctrl_label = "HAIC_alone"
  ),
  I_CONC = list(
    csv = "composite_I_CONC_cohort.csv",
    trt_col = "trt_haic_i_conc",
    trt_label = "HAIC+I_concurrent",
    ctrl_label = "HAIC_alone"
  ),
  IT_CONC = list(
    csv = "composite_IT_CONC_cohort.csv",
    trt_col = "trt_haic_it_conc",
    trt_label = "HAIC+I+T_concurrent",
    ctrl_label = "HAIC_alone"
  )
)

nk_env <- Sys.getenv("RMS_RCS_NK", "3")
RMS_RCS_NK <- suppressWarnings(as.integer(nk_env))
if (length(RMS_RCS_NK) != 1L || is.na(RMS_RCS_NK) || RMS_RCS_NK < 3L) RMS_RCS_NK <- 3L
nb_env <- Sys.getenv("RMS_RCS_N_BOOT", "")
N_BOOT <- if (nzchar(nb_env)) as.integer(nb_env) else 200L
if (length(N_BOOT) != 1L || is.na(N_BOOT) || N_BOOT < 2L) N_BOOT <- 200L
MIN_N <- 50L

FIG_WIDTH  <- 14.7
FIG_HEIGHT <- 17.7

cat("SCRIPT_DIR:", SCRIPT_DIR, "\nRMS_RCS_NK:", RMS_RCS_NK,
    "| N_BOOT:", N_BOOT, "| cohort_key:", cohort_key, "\n\n")

# ── 7×5 Matrix 指标配置 ──────────────────────────────────────────────────────
MATRIX_CONFIG <- list(
  list(
    row_label = "AFP",
    cols = list(
      list(col = "afp",               label = "AFP\nBaseline",                    type = "static", log_transform = TRUE),
      list(col = "afp_pre3",          label = "AFP\nPre-HAIC-3",                 type = "static", log_transform = TRUE),
      list(col = "afp_pre_it",        label = "AFP\nPre-IT",                     type = "static", log_transform = TRUE),
      list(col = "afp_change_pre3",   label = "AFP\nPre-HAIC-3 Change Rate",    type = "dynamic"),
      list(col = "afp_change_pre_it", label = "AFP\nPre-IT Change Rate",         type = "dynamic")
    )
  ),
  list(
    row_label = "PIVKA",
    cols = list(
      list(col = "pivka",               label = "PIVKA\nBaseline",                  type = "static", log_transform = TRUE),
      list(col = "pivka_pre3",          label = "PIVKA\nPre-HAIC-3",               type = "static", log_transform = TRUE),
      list(col = "pivka_pre_it",        label = "PIVKA\nPre-IT",                   type = "static", log_transform = TRUE),
      list(col = "pivka_change_pre3",   label = "PIVKA\nPre-HAIC-3 Change Rate",  type = "dynamic"),
      list(col = "pivka_change_pre_it", label = "PIVKA\nPre-IT Change Rate",       type = "dynamic")
    )
  ),
  list(
    row_label = "PIV",
    cols = list(
      list(col = "piv_bl",              label = "PIV\nBaseline",                    type = "static", log_transform = TRUE),
      list(col = "piv_pre3",            label = "PIV\nPre-HAIC-3",                 type = "static", log_transform = TRUE),
      list(col = "piv_pre_it",          label = "PIV\nPre-IT",                     type = "static", log_transform = TRUE),
      list(col = "piv_change_pre3",     label = "PIV\nPre-HAIC-3 Change Rate",    type = "dynamic"),
      list(col = "piv_change_pre_it",   label = "PIV\nPre-IT Change Rate",         type = "dynamic")
    )
  ),
  list(
    row_label = "SII",
    cols = list(
      list(col = "sii_bl",              label = "SII\nBaseline",                    type = "static", log_transform = TRUE),
      list(col = "sii_pre3",            label = "SII\nPre-HAIC-3",                 type = "static", log_transform = TRUE),
      list(col = "sii_pre_it",          label = "SII\nPre-IT",                     type = "static", log_transform = TRUE),
      list(col = "sii_change_pre3",     label = "SII\nPre-HAIC-3 Change Rate",    type = "dynamic"),
      list(col = "sii_change_pre_it",   label = "SII\nPre-IT Change Rate",         type = "dynamic")
    )
  ),
  list(
    row_label = "NLR",
    cols = list(
      list(col = "nlr_bl",              label = "NLR\nBaseline",                    type = "static", log_transform = FALSE),
      list(col = "nlr_pre3",            label = "NLR\nPre-HAIC-3",                 type = "static", log_transform = FALSE),
      list(col = "nlr_pre_it",          label = "NLR\nPre-IT",                     type = "static", log_transform = FALSE),
      list(col = "nlr_change_pre3",     label = "NLR\nPre-HAIC-3 Change Rate",    type = "dynamic"),
      list(col = "nlr_change_pre_it",   label = "NLR\nPre-IT Change Rate",         type = "dynamic")
    )
  ),
  list(
    row_label = "PLR",
    cols = list(
      list(col = "plr_bl",              label = "PLR\nBaseline",                    type = "static", log_transform = FALSE),
      list(col = "plr_pre3",            label = "PLR\nPre-HAIC-3",                 type = "static", log_transform = FALSE),
      list(col = "plr_pre_it",          label = "PLR\nPre-IT",                     type = "static", log_transform = FALSE),
      list(col = "plr_change_pre3",     label = "PLR\nPre-HAIC-3 Change Rate",    type = "dynamic"),
      list(col = "plr_change_pre_it",   label = "PLR\nPre-IT Change Rate",         type = "dynamic")
    )
  ),
  list(
    row_label = "MONOCYTE",
    cols = list(
      list(col = "mono_bl",             label = "Monocyte\nBaseline",               type = "static", log_transform = FALSE),
      list(col = "mono_pre3",           label = "Monocyte\nPre-HAIC-3",            type = "static", log_transform = FALSE),
      list(col = "mono_pre_it",         label = "Monocyte\nPre-IT",                type = "static", log_transform = FALSE),
      list(col = "mono_change_pre3",    label = "Monocyte\nPre-HAIC-3 Change Rate", type = "dynamic"),
      list(col = "mono_change_pre_it",  label = "Monocyte\nPre-IT Change Rate",     type = "dynamic")
    )
  ),
  list(
    row_label = "ALBI",
    cols = list(
      list(col = "albi_bl",             label = "ALBI\nBaseline",                   type = "static", log_transform = FALSE),
      list(col = "albi_pre3",           label = "ALBI\nPre-HAIC-3",                type = "static", log_transform = FALSE),
      list(col = "albi_pre_it",         label = "ALBI\nPre-IT",                    type = "static", log_transform = FALSE),
      list(col = "albi_change_pre3",    label = "ALBI\nPre-HAIC-3 Change Rate",   type = "dynamic"),
      list(col = "albi_change_pre_it",  label = "ALBI\nPre-IT Change Rate",        type = "dynamic")
    )
  )
)

COL_HEADERS <- c("Baseline", "Pre-HAIC-3", "Pre-IT", "Pre-HAIC-3\nChange Rate", "Pre-IT\nChange Rate")

# ── IPTW ─────────────────────────────────────────────────────────────────────
compute_iptw <- function(df_input, var_col, trt_col, ps_vars_base) {
  ps_vars_use <- ps_vars_base[ps_vars_base %in% names(df_input)]
  ps_vars_use <- ps_vars_use[sapply(ps_vars_use,
    function(v) sd(df_input[[v]], na.rm = TRUE) > 0)]

  df_cc <- df_input[complete.cases(df_input[, ps_vars_use, drop = FALSE]), ]
  if (nrow(df_cc) < 20) return(df_cc[0, , drop = FALSE])
  X_sc <- scale(as.matrix(df_cc[, ps_vars_use, drop = FALSE]))
  y    <- df_cc[[trt_col]]
  cv_fit  <- cv.glmnet(X_sc, y, family = "binomial", alpha = 0, nfolds = 5)
  ps_prob <- as.numeric(predict(cv_fit, newx = X_sc, s = "lambda.min", type = "response"))
  ps_prob <- pmin(pmax(ps_prob, 0.05), 0.95)
  p_treat <- mean(y)
  df_cc$sw <- ifelse(df_cc[[trt_col]] == 1, p_treat / ps_prob, (1 - p_treat) / (1 - ps_prob))
  df_cc
}

# ── anova P 值提取 ──────────────────────────────────────────────────────────
parse_anova_p <- function(cell) {
  x <- as.character(cell)[1]
  if (is.na(x) || !nzchar(x)) return(NA_real_)
  if (grepl("^<", x)) {
    z <- sub("^<\\.?", "", x)
    return(suppressWarnings(as.numeric(z)))
  }
  suppressWarnings(as.numeric(x))
}

extract_rms_anova_p <- function(fit, trt_col) {
  a <- anova(fit)
  rn <- rownames(a)
  pcol <- if ("P" %in% colnames(a)) "P" else ncol(a)
  int_p <- NA_real_
  nonlin_p <- NA_real_
  idx_int <- grep(paste0("^", trt_col, " \\* [^\\(]+  \\(Factor\\+Higher Order Factors\\)$"), rn)
  if (length(idx_int)) int_p <- parse_anova_p(a[idx_int[1], pcol])
  idx_nl <- grep("Nonlinear Interaction : f(A,B) vs. AB", rn, fixed = TRUE)
  if (length(idx_nl)) nonlin_p <- parse_anova_p(a[idx_nl[1], pcol])
  list(int_p = int_p, nonlin_p = nonlin_p)
}

# ── 核心模型函数 ────────────────────────────────────────────────────────────
fit_boot_cph_rms <- function(df_b, surv_time_col, nk, trt_col) {
  dd_b <- suppressWarnings(datadist(df_b[, c(trt_col, "rcsx"), drop = FALSE]))
  options(datadist = dd_b)
  fml <- as.formula(paste0(
    "Surv(", surv_time_col, ", death_status) ~ ", trt_col, " * rcs(rcsx, ", nk, ")"
  ))
  cph(fml, data = df_b, weights = sw, x = TRUE, y = TRUE, robust = FALSE)
}

predict_hr_curve <- function(fit, x_grid, trt_col) {
  nd1 <- setNames(data.frame(1L, x_grid), c(trt_col, "rcsx"))
  nd0 <- setNames(data.frame(0L, x_grid), c(trt_col, "rcsx"))
  lp1 <- as.numeric(predict(fit, nd1, type = "lp"))
  lp0 <- as.numeric(predict(fit, nd0, type = "lp"))
  exp(lp1 - lp0)
}

find_crossings <- function(df_plot, hr_ref) {
  cps <- numeric(0)
  for (i in seq_len(nrow(df_plot) - 1)) {
    y1 <- df_plot$hr[i]     - hr_ref
    y2 <- df_plot$hr[i + 1] - hr_ref
    if (y1 * y2 < 0) {
      frac <- (hr_ref - df_plot$hr[i]) / (df_plot$hr[i + 1] - df_plot$hr[i])
      cps  <- c(cps, df_plot$x[i] + frac * (df_plot$x[i + 1] - df_plot$x[i]))
    }
  }
  cps
}

# ── 单个子图绘制（IPTW 版，精简适配 matrix panel）──────────────────────────
make_panel_plot <- function(df_sub, var_col, var_label, surv_time_col,
                            plot_type, log_transform = FALSE,
                            trt_col, trt_label, ctrl_label, ps_vars_base) {
  # IPTW 加权
  df_sub <- compute_iptw(df_sub, var_col, trt_col, ps_vars_base)
  df_sub <- df_sub[!is.na(df_sub[[var_col]]), ]
  n <- nrow(df_sub)
  if (n < MIN_N) {
    return(ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = paste0(var_label, "\nn=", n, " < ", MIN_N),
               size = 3.5, color = "grey50") +
      theme_void() +
      theme(plot.margin = margin(2, 2, 2, 2)))
  }

  if (plot_type == "static" && log_transform) {
    df_sub$rcsx <- log1p(pmax(df_sub[[var_col]], 0))
    inv_transform <- function(x) expm1(x)
  } else {
    df_sub$rcsx <- df_sub[[var_col]]
    inv_transform <- identity
  }
  x_vals <- df_sub$rcsx

  nk <- RMS_RCS_NK
  dd <- suppressWarnings(datadist(df_sub[, c(trt_col, "rcsx"), drop = FALSE]))
  options(datadist = dd)

  result <- tryCatch({
    fml <- as.formula(paste0(
      "Surv(", surv_time_col, ", death_status) ~ ", trt_col, " * rcs(rcsx, ", nk, ")"
    ))
    fit <- cph(fml, data = df_sub, weights = sw, x = TRUE, y = TRUE, robust = FALSE)
    ap  <- extract_rms_anova_p(fit, trt_col)

    xlim     <- quantile(df_sub$rcsx, c(0.05, 0.95))
    x_grid   <- seq(xlim[1], xlim[2], length.out = 200)
    hr_curve <- predict_hr_curve(fit, x_grid, trt_col)

    set.seed(42)
    hr_boot <- matrix(NA_real_, N_BOOT, length(x_grid))
    for (b in seq_len(N_BOOT)) {
      idx  <- sample.int(nrow(df_sub), replace = TRUE)
      df_b <- df_sub[idx, ]
      tryCatch({
        fit_b <- fit_boot_cph_rms(df_b, surv_time_col, nk, trt_col)
        hr_boot[b, ] <- predict_hr_curve(fit_b, x_grid, trt_col)
      }, error = function(e) NULL)
    }
    hr_lo <- apply(hr_boot, 2, quantile, 0.025, na.rm = TRUE)
    hr_hi <- apply(hr_boot, 2, quantile, 0.975, na.rm = TRUE)
    eps_hr <- 0.05
    hr_lo[is.na(hr_lo)] <- eps_hr
    hr_hi[is.na(hr_hi)] <- pmax(10, hr_curve, na.rm = TRUE)
    hr_lo    <- pmax(hr_lo, eps_hr)
    hr_hi    <- pmax(hr_hi, hr_lo * 1.001)
    hr_curve <- pmax(hr_curve, eps_hr)

    HR_REFS    <- c(1.0, 0.9, 0.8, 0.7)
    REF_COLORS <- c("#333333", "#0072B2", "#DC0000", "#7E6148")
    REF_LT     <- c("dashed", "dotdash", "dotdash", "dotdash")
    REF_LW     <- c(0.8, 0.6, 0.6, 0.6)

    plot_df    <- data.frame(x = x_grid, hr = hr_curve, hr_lo = hr_lo, hr_hi = hr_hi)
    cross_list <- lapply(HR_REFS, function(r) find_crossings(plot_df, r))

    fmt_cp <- function(cp) {
      val <- inv_transform(cp)
      if (abs(val) >= 1000) sprintf("%.0f", val)
      else if (abs(val) >= 10) sprintf("%.1f", val)
      else sprintf("%.2f", val)
    }

    ok_rug <- is.finite(x_vals) & x_vals >= xlim[1] & x_vals <= xlim[2]
    rug_x  <- x_vals[ok_rug]
    if (length(rug_x) < 1L && any(is.finite(x_vals))) rug_x <- x_vals[is.finite(x_vals)]

    p <- ggplot(plot_df, aes(x = x)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.25, ymax = 1,
               fill = "#E8F5E9", alpha = 0.35) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = 2.5,
               fill = "#FFF3E0", alpha = 0.35) +
      geom_ribbon(aes(ymin = hr_lo, ymax = hr_hi), fill = "#3C5488", alpha = 0.18) +
      geom_line(aes(y = hr), color = "#3C5488", linewidth = 1.2)

    for (ri in seq_along(HR_REFS)) {
      p <- p + geom_hline(yintercept = HR_REFS[ri], linetype = REF_LT[ri],
                          color = REF_COLORS[ri], linewidth = REF_LW[ri], alpha = 0.9)
    }

    cross_y_labels <- c(1.38, 0.60, 0.50, 0.38)
    for (ri in seq_along(HR_REFS)) {
      cps <- cross_list[[ri]]
      if (length(cps) > 0) {
        for (cp in cps) {
          if (plot_type == "dynamic") {
            lab <- sprintf("%.1f%%", cp)
          } else {
            lab <- fmt_cp(cp)
          }
          p <- p +
            geom_vline(xintercept = cp, linetype = "dotted",
                       color = REF_COLORS[ri], linewidth = 0.5, alpha = 0.8) +
            annotate("point", x = cp, y = HR_REFS[ri],
                     color = REF_COLORS[ri], size = 1.5, shape = 18) +
            annotate("label", x = cp, y = cross_y_labels[ri], label = lab,
                     color = REF_COLORS[ri], size = 2.2, fontface = "bold",
                     fill = "white", label.size = 0.2,
                     label.padding = unit(0.1, "lines"))
        }
      }
    }

    if (plot_type == "dynamic") {
      p <- p + geom_vline(xintercept = 0, linetype = "dotted", color = "#888888", linewidth = 0.5)
    }

    p <- p +
      geom_rug(data = data.frame(xr = rug_x),
               aes(x = xr), sides = "b", alpha = 0.3, color = "#333333",
               length = unit(1.5, "mm"), linewidth = 0.3, inherit.aes = FALSE) +
      labs(
        title    = var_label,
        subtitle = sprintf("Int.p=%s | NL.p=%s | n=%d",
          ifelse(is.na(ap$int_p),    "NA", sprintf("%.3f", ap$int_p)),
          ifelse(is.na(ap$nonlin_p), "NA", sprintf("%.3f", ap$nonlin_p)), n),
        x = NULL, y = NULL
      ) +
      scale_y_log10(breaks = c(0.3, 0.5, 0.7, 0.9, 1.0, 1.5, 2.0),
                    labels = c("0.3", "0.5", "0.7", "0.9", "1.0", "1.5", "2.0")) +
      coord_cartesian(ylim = c(0.25, 2.5), xlim = xlim) +
      theme_bw(base_size = 9) +
      theme(
        plot.title    = element_text(size = 8, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 6.5, hjust = 0.5, color = "grey40"),
        axis.text     = element_text(size = 6.5),
        axis.title    = element_blank(),
        plot.margin   = margin(3, 4, 3, 4)
      )

    if (plot_type == "static" && log_transform) {
      raw_range <- expm1(xlim)
      candidate_ticks <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000,
                           10000, 20000, 50000, 100000, 200000, 500000)
      raw_ticks <- candidate_ticks[candidate_ticks >= raw_range[1] & candidate_ticks <= raw_range[2]]
      if (length(raw_ticks) < 3) {
        raw_ticks <- pretty(raw_range, n = 4)
        raw_ticks <- raw_ticks[raw_ticks > 0]
      }
      if (length(raw_ticks) > 6)
        raw_ticks <- raw_ticks[seq(1, length(raw_ticks), length.out = min(6, length(raw_ticks)))]
      log_ticks <- log1p(raw_ticks)
      fmt_label <- function(v) ifelse(v >= 1000, sprintf("%gK", v / 1000),
                                      ifelse(v >= 1, sprintf("%g", v), sprintf("%.1f", v)))
      p <- p + scale_x_continuous(breaks = log_ticks, labels = fmt_label(raw_ticks),
                                   limits = xlim,
                                   expand = ggplot2::expansion(mult = c(0.02, 0.02)))
    } else {
      p <- p + scale_x_continuous(limits = xlim,
                                   expand = ggplot2::expansion(mult = c(0.03, 0.03)),
                                   labels = function(x) {
                                     vals <- inv_transform(x)
                                     ifelse(abs(vals) >= 1000, sprintf("%.0f", vals),
                                       ifelse(abs(vals) >= 10, sprintf("%.0f", vals),
                                         sprintf("%.1f", vals)))
                                   })
    }

    p
  }, error = function(e) {
    cat("  ERROR [", var_col, "]:", conditionMessage(e), "\n")
    ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste0(var_label, "\nERROR"),
               size = 3, color = "red") +
      theme_void() +
      theme(plot.margin = margin(2, 2, 2, 2))
  })

  result
}

# ── 组装 7×5 Matrix Panel ────────────────────────────────────────────────────
build_matrix_panel <- function(df, surv_time_col, trt_col, trt_label, ctrl_label,
                               ps_vars_static, ps_vars_dynamic, title_main) {
  n_rows <- length(MATRIX_CONFIG)
  n_cols <- 5
  plot_list <- vector("list", n_rows * n_cols)

  for (ri in seq_along(MATRIX_CONFIG)) {
    row_cfg <- MATRIX_CONFIG[[ri]]
    for (ci in seq_along(row_cfg$cols)) {
      col_cfg <- row_cfg$cols[[ci]]
      idx <- (ri - 1) * n_cols + ci

      ps_use <- if (col_cfg$type == "dynamic") ps_vars_dynamic else ps_vars_static

      if (col_cfg$col %in% names(df)) {
        cat(sprintf("  [%d/%d] %s ...\n", idx, n_rows * n_cols, col_cfg$col))
        plot_list[[idx]] <- make_panel_plot(
          df_sub        = df,
          var_col       = col_cfg$col,
          var_label     = col_cfg$label,
          surv_time_col = surv_time_col,
          plot_type     = col_cfg$type,
          log_transform = if (!is.null(col_cfg$log_transform)) col_cfg$log_transform else FALSE,
          trt_col       = trt_col,
          trt_label     = trt_label,
          ctrl_label    = ctrl_label,
          ps_vars_base  = ps_use
        )
      } else {
        plot_list[[idx]] <- ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = paste0(col_cfg$label, "\nN/A"),
                   size = 3, color = "grey60") +
          theme_void() +
          theme(plot.margin = margin(2, 2, 2, 2))
      }
    }
  }

  # 列标题
  col_header_grobs <- lapply(COL_HEADERS, function(h) {
    textGrob(h, gp = gpar(fontsize = 10, fontface = "bold", col = "#333333"))
  })

  # 行标题
  row_label_grobs <- lapply(MATRIX_CONFIG, function(rc) {
    textGrob(rc$row_label, rot = 90,
             gp = gpar(fontsize = 11, fontface = "bold", col = "#333333"))
  })

  # layout matrix: (n_rows+1) × (n_cols+1)
  layout_mat <- matrix(NA, nrow = n_rows + 1, ncol = n_cols + 1)

  col_start <- length(plot_list) + 1
  for (ci in seq_len(n_cols)) {
    layout_mat[1, ci + 1] <- col_start + ci - 1
  }
  row_start <- col_start + n_cols
  for (ri in seq_len(n_rows)) {
    layout_mat[ri + 1, 1] <- row_start + ri - 1
  }
  for (ri in seq_len(n_rows)) {
    for (ci in seq_len(n_cols)) {
      layout_mat[ri + 1, ci + 1] <- (ri - 1) * n_cols + ci
    }
  }

  all_grobs <- c(plot_list,
                 lapply(col_header_grobs, function(g) g),
                 lapply(row_label_grobs, function(g) g))

  all_grobs <- lapply(all_grobs, function(g) {
    if (inherits(g, "gg")) ggplotGrob(g) else g
  })

  combined <- arrangeGrob(
    grobs  = all_grobs,
    layout_matrix = layout_mat,
    widths  = unit(c(0.6, rep(1, n_cols)), c("cm", rep("null", n_cols))),
    heights = unit(c(1.0, rep(1, n_rows)), c("cm", rep("null", n_rows))),
    top     = textGrob(title_main, gp = gpar(fontsize = 14, fontface = "bold")),
    bottom  = textGrob(
      sprintf("HR (%s vs %s) | IPTW-weighted RCS(nk=%d) | Bootstrap CI (n=%d)",
              trt_label, ctrl_label, RMS_RCS_NK, N_BOOT),
      gp = gpar(fontsize = 8, col = "grey50"))
  )

  combined
}

# ── 运行一个队列 ────────────────────────────────────────────────────────────
run_cohort <- function(key) {
  reg <- COHORT_REGISTRY[[key]]
  data_csv <- file.path(SCRIPT_DIR, reg$csv)
  if (!file.exists(data_csv)) {
    cat("Data file not found:", data_csv, "- SKIP\n")
    return(invisible(NULL))
  }
  df0 <- read.csv(data_csv, stringsAsFactors = FALSE)
  cat("\n", strrep("#", 70), "\n", sep = "")
  cat("COHORT:", key, "| File:", reg$csv, "| nrow =", nrow(df0),
      "| trt:", reg$trt_col, "\n")
  cat(strrep("#", 70), "\n")

  if (!reg$trt_col %in% names(df0)) {
    cat("  ERROR: trt_col", reg$trt_col, "not in data\n")
    return(invisible(NULL))
  }

  if (!"os_lm" %in% names(df0)) df0$os_lm <- df0$os_months - 42 / 30.44

  # ── 衍生列 ──────────────────────────────────────────────────────────────
  if (!"sii_bl" %in% names(df0) && all(c("plt", "nlr_bl") %in% names(df0)))
    df0$sii_bl <- df0$plt * df0$nlr_bl
  if (!"piv_bl" %in% names(df0) && all(c("mono_bl", "plt", "nlr_bl") %in% names(df0)))
    df0$piv_bl <- df0$mono_bl * df0$plt * df0$nlr_bl
  if (!"sii_pre3" %in% names(df0) && all(c("plt_pre3", "nlr_pre3") %in% names(df0)))
    df0$sii_pre3 <- df0$plt_pre3 * df0$nlr_pre3
  if (!"piv_pre3" %in% names(df0) && all(c("mono_pre3", "plt_pre3", "nlr_pre3") %in% names(df0)))
    df0$piv_pre3 <- df0$mono_pre3 * df0$plt_pre3 * df0$nlr_pre3

  # Pre-IT 衍生列
  if (!"nlr_pre_it" %in% names(df0) && all(c("neut_pre_it", "lymph_pre_it") %in% names(df0)))
    df0$nlr_pre_it <- df0$neut_pre_it / ifelse(df0$lymph_pre_it == 0, NA, df0$lymph_pre_it)
  if (!"plr_pre_it" %in% names(df0) && all(c("plt_pre_it", "lymph_pre_it") %in% names(df0)))
    df0$plr_pre_it <- df0$plt_pre_it / ifelse(df0$lymph_pre_it == 0, NA, df0$lymph_pre_it)
  if (!"sii_pre_it" %in% names(df0) && all(c("plt_pre_it", "neut_pre_it", "lymph_pre_it") %in% names(df0)))
    df0$sii_pre_it <- df0$plt_pre_it * df0$neut_pre_it / ifelse(df0$lymph_pre_it == 0, NA, df0$lymph_pre_it)
  if (!"piv_pre_it" %in% names(df0) && all(c("mono_pre_it", "plt_pre_it", "neut_pre_it", "lymph_pre_it") %in% names(df0)))
    df0$piv_pre_it <- df0$mono_pre_it * df0$plt_pre_it * df0$neut_pre_it / ifelse(df0$lymph_pre_it == 0, NA, df0$lymph_pre_it)
  if (!"albi_pre_it" %in% names(df0) && all(c("tbil_pre_it", "alb_pre_it") %in% names(df0)))
    df0$albi_pre_it <- ifelse(df0$tbil_pre_it > 0 & !is.na(df0$tbil_pre_it) & !is.na(df0$alb_pre_it),
                              0.66 * log10(df0$tbil_pre_it) - 0.085 * df0$alb_pre_it, NA)

  # 变化率
  change_maps <- list(
    list(pre = "afp_pre3",     bl = "afp",     out = "afp_change_pre3"),
    list(pre = "pivka_pre3",   bl = "pivka",   out = "pivka_change_pre3"),
    list(pre = "piv_pre3",     bl = "piv_bl",  out = "piv_change_pre3"),
    list(pre = "sii_pre3",     bl = "sii_bl",  out = "sii_change_pre3"),
    list(pre = "nlr_pre3",     bl = "nlr_bl",  out = "nlr_change_pre3"),
    list(pre = "plr_pre3",     bl = "plr_bl",  out = "plr_change_pre3"),
    list(pre = "mono_pre3",    bl = "mono_bl", out = "mono_change_pre3"),
    list(pre = "albi_pre3",    bl = "albi_bl", out = "albi_change_pre3"),
    list(pre = "afp_pre_it",   bl = "afp",     out = "afp_change_pre_it"),
    list(pre = "pivka_pre_it", bl = "pivka",   out = "pivka_change_pre_it"),
    list(pre = "piv_pre_it",   bl = "piv_bl",  out = "piv_change_pre_it"),
    list(pre = "sii_pre_it",   bl = "sii_bl",  out = "sii_change_pre_it"),
    list(pre = "nlr_pre_it",   bl = "nlr_bl",  out = "nlr_change_pre_it"),
    list(pre = "plr_pre_it",   bl = "plr_bl",  out = "plr_change_pre_it"),
    list(pre = "mono_pre_it",  bl = "mono_bl", out = "mono_change_pre_it"),
    list(pre = "albi_pre_it",  bl = "albi_bl", out = "albi_change_pre_it")
  )
  for (m in change_maps) {
    if (!(m$out %in% names(df0)) && all(c(m$pre, m$bl) %in% names(df0))) {
      bl_vals  <- as.numeric(df0[[m$bl]])
      pre_vals <- as.numeric(df0[[m$pre]])
      df0[[m$out]] <- ifelse(!is.na(bl_vals) & !is.na(pre_vals) & bl_vals != 0,
                             (pre_vals - bl_vals) / abs(bl_vals) * 100, NA)
    }
  }

  df_lm    <- df0[!is.na(df0$os_lm) & df0$os_lm > 0, ]
  df_total <- df0[!is.na(df0$os_months) & df0$os_months > 0, ]

  # PS 变量
  ps_vars_full <- c("albi_bl", "inr", "plt",
                    "tumor_max_diameter_cm", "tumor_count_enc",
                    "pvtt_grade", "hvtt_binary", "ivc_ra_binary", "ascites_score_enc",
                    "log_afp_bl", "log_pivka_bl",
                    "metastasis_binary", "lymph_node_binary",
                    "neut_bl", "lymph_bl", "mono_bl")
  ps_vars_dynamic_extra <- c(
    "afp_pre_it", "pivka_pre_it",
    "albi_pre_it", "neut_pre_it", "ast_pre_it", "alt_pre_it",
    "afp_change_pre_it", "pivka_change_pre_it"
  )

  PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", "..", ".."), winslash = "/")
  base_out <- file.path(PROJECT_ROOT, "output", "step1_rcs_interaction", "afp_pivka_composite", key)

  for (arm_info in list(
    list(df = df_lm,    subdir = "landmark",  surv = "os_lm",      tag = "Landmark OS"),
    list(df = df_total, subdir = "total_os",  surv = "os_months",  tag = "Total OS")
  )) {
    df_arm <- arm_info$df
    out_dir <- file.path(base_out, arm_info$subdir)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    ps_static  <- ps_vars_full[ps_vars_full %in% names(df_arm)]
    ps_static  <- ps_static[sapply(ps_static, function(v) sd(df_arm[[v]], na.rm = TRUE) > 0)]
    ps_dynamic <- c(ps_vars_full, ps_vars_dynamic_extra)
    ps_dynamic <- ps_dynamic[ps_dynamic %in% names(df_arm)]
    ps_dynamic <- ps_dynamic[sapply(ps_dynamic, function(v) sd(df_arm[[v]], na.rm = TRUE) > 0)]

    cat(sprintf("\n  Building matrix panel [%s | %s] ...\n", key, arm_info$tag))
    panel <- build_matrix_panel(
      df            = df_arm,
      surv_time_col = arm_info$surv,
      trt_col       = reg$trt_col,
      trt_label     = reg$trt_label,
      ctrl_label    = reg$ctrl_label,
      ps_vars_static  = ps_static,
      ps_vars_dynamic = ps_dynamic,
      title_main    = sprintf("%s vs %s — RCS Matrix Panel (%s, IPTW)",
                              reg$trt_label, reg$ctrl_label, arm_info$tag)
    )
    ggsave(file.path(out_dir, "rcs_matrix_panel_8x5.pdf"),
           panel, width = FIG_WIDTH, height = FIG_HEIGHT, device = "pdf")
    ggsave(file.path(out_dir, "rcs_matrix_panel_8x5.png"),
           panel, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = 300)
    cat(sprintf("  Saved: %s/rcs_matrix_panel_8x5.pdf/.png\n", arm_info$subdir))
  }
}

# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════
if (cohort_key == "ALL") {
  for (k in names(COHORT_REGISTRY)) run_cohort(k)
} else if (cohort_key %in% names(COHORT_REGISTRY)) {
  run_cohort(cohort_key)
} else {
  stop("Unknown cohort key: ", cohort_key,
       "\nUse one of: ", paste(names(COHORT_REGISTRY), collapse = ", "), ", ALL",
       call. = FALSE)
}

cat("\n=== Matrix panel generation completed ===\n")
