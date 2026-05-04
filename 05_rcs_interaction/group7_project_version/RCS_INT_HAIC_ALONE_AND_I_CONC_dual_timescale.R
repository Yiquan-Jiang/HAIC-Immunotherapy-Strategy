#!/usr/bin/env Rscript
# =============================================================================
# RCS_INT_HAIC_ALONE_AND_I_CONC — rms::rcs() + Overlap Weighting rms::cph()，双时间尺度
#
# 队列: TIDY main_group 为 HAIC_alone vs HAIC+I_concurrent；暴露项
#   Surv(...) ~ trt_haic_i_conc * rcs(rcsx, nk)
#
# 依赖: survival, rms, Hmisc, ggplot2, dplyr, glmnet, gridExtra, grid
#
# 数据: RCS_INT_HAIC_ALONE_AND_I_CONC_cohort.csv（先运行 RCS_INT_HAIC_ALONE_AND_I_CONC_build_cohort.py）
# 输出: output/RCS_INT_HAIC_ALONE_AND_I_CONC/{landmark,total_os}/
# 环境变量: RMS_RCS_DATA_CSV, RMS_RCS_OUT_DIR, RMS_RCS_N_BOOT, RMS_RCS_NK（可选）
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

args_all <- commandArgs(trailingOnly = FALSE)
fa <- args_all[grepl("^--file=", args_all)]
SCRIPT_DIR <- if (length(fa)) {
  dirname(normalizePath(sub("^--file=", "", fa[1]), winslash = "/"))
} else {
  getwd()
}

DATA_CSV <- Sys.getenv(
  "RMS_RCS_DATA_CSV",
  unset = normalizePath(file.path(SCRIPT_DIR, "data", "RCS_INT_HAIC_ALONE_AND_I_CONC_cohort.csv"), mustWork = FALSE)
)
BASE_OUT <- Sys.getenv(
  "RMS_RCS_OUT_DIR",
  unset = file.path(SCRIPT_DIR, "output", "RCS_INT_HAIC_ALONE_AND_I_CONC")
)
dir.create(BASE_OUT, showWarnings = FALSE, recursive = TRUE)
BASE_OUT <- normalizePath(BASE_OUT, winslash = "/", mustWork = TRUE)

# rcs(x, nk)：rms 要求 nk>=3；默认固定 nk=3（可用 RMS_RCS_NK 覆盖）
nk_env <- Sys.getenv("RMS_RCS_NK", "3")
RMS_RCS_NK <- suppressWarnings(as.integer(nk_env))
if (length(RMS_RCS_NK) != 1L || is.na(RMS_RCS_NK) || RMS_RCS_NK < 3L) RMS_RCS_NK <- 3L
nb_env <- Sys.getenv("RMS_RCS_N_BOOT", "")
N_BOOT <- if (nzchar(nb_env)) as.integer(nb_env) else 200L
if (length(N_BOOT) != 1L || is.na(N_BOOT) || N_BOOT < 2L) N_BOOT <- 200L
MIN_N <- 50L

# 治疗变量列名与标签（可通过环境变量覆盖，默认 HAIC_alone vs HAIC+I_concurrent）
TRT_COL   <- Sys.getenv("RMS_RCS_TRT_COL",   unset = "trt_haic_i_conc")
TRT_LABEL <- Sys.getenv("RMS_RCS_TRT_LABEL", unset = "HAIC+I_concurrent")
CTRL_LABEL<- Sys.getenv("RMS_RCS_CTRL_LABEL",unset = "HAIC_alone")

cat("TRT_COL:", TRT_COL, "| TRT_LABEL:", TRT_LABEL, "| CTRL_LABEL:", CTRL_LABEL, "\n")

# ── Overlap Weighting（glmnet ridge PS）────────────────────────────────────────
impute_ps_vars <- function(df, ps_vars_full) {
  for (v in ps_vars_full) {
    if (v %in% names(df) && any(is.na(df[[v]]))) {
      med_val <- median(df[[v]], na.rm = TRUE)
      df[[v]][is.na(df[[v]])] <- med_val
    }
  }
  df
}

compute_overlap_weight_static <- function(df_input, var_col, ps_vars_base) {
  # 静态变量分析：使用完整ps_vars_base，不排除任何变量
  ps_vars_use <- ps_vars_base[ps_vars_base %in% names(df_input)]
  ps_vars_use <- ps_vars_use[sapply(ps_vars_use,
    function(v) sd(df_input[[v]], na.rm = TRUE) > 0)]

  df_cc <- df_input[complete.cases(df_input[, ps_vars_use, drop = FALSE]), ]
  if (nrow(df_cc) < 20) return(df_cc[0, , drop = FALSE])
  X_sc <- scale(as.matrix(df_cc[, ps_vars_use, drop = FALSE]))
  y    <- df_cc[[TRT_COL]]
  cv_fit  <- cv.glmnet(X_sc, y, family = "binomial", alpha = 0, nfolds = 5)
  ps_prob <- as.numeric(predict(cv_fit, newx = X_sc, s = "lambda.min", type = "response"))
  ps_prob <- pmin(pmax(ps_prob, 0.05), 0.95)
  
  # Overlap weighting: w = (1-ps)^trt * ps^(1-trt)
  # 标准化权重使其均值为1
  df_cc$sw_raw <- (1 - ps_prob)^y * ps_prob^(1 - y)
  df_cc$sw <- df_cc$sw_raw / mean(df_cc$sw_raw)
  df_cc
}

compute_overlap_weight_dynamic <- function(df_input, var_col, ps_vars_base) {
  # 动态变量分析：使用完整ps_vars_base + 额外纳入变化率和pre-HAIC-3 ALBI
  extra_vars <- c("pivka_change_pre3", "afp_change_pre3", "albi_pre3")
  ps_vars_extended <- c(ps_vars_base, extra_vars)
  ps_vars_use <- ps_vars_extended[ps_vars_extended %in% names(df_input)]
  ps_vars_use <- ps_vars_use[sapply(ps_vars_use,
    function(v) sd(df_input[[v]], na.rm = TRUE) > 0)]

  df_cc <- df_input[complete.cases(df_input[, ps_vars_use, drop = FALSE]), ]
  if (nrow(df_cc) < 20) return(df_cc[0, , drop = FALSE])
  X_sc <- scale(as.matrix(df_cc[, ps_vars_use, drop = FALSE]))
  y    <- df_cc[[TRT_COL]]
  cv_fit  <- cv.glmnet(X_sc, y, family = "binomial", alpha = 0, nfolds = 5)
  ps_prob <- as.numeric(predict(cv_fit, newx = X_sc, s = "lambda.min", type = "response"))
  ps_prob <- pmin(pmax(ps_prob, 0.05), 0.95)
  
  # Overlap weighting: w = (1-ps)^trt * ps^(1-trt)
  df_cc$sw_raw <- (1 - ps_prob)^y * ps_prob^(1 - y)
  df_cc$sw <- df_cc$sw_raw / mean(df_cc$sw_raw)
  df_cc
}

# ── 从 rms::anova(fit) 提取治疗×暴露交互与非线性交互 P 值 ────────────────────
# 注意：anova 表子行名常带前导空格（缩进层级），勿用 ^ 锚定行首匹配 Nonlinear Interaction
parse_anova_p <- function(cell) {
  x <- as.character(cell)[1]
  if (is.na(x) || !nzchar(x)) return(NA_real_)
  if (grepl("^<", x)) {
    z <- sub("^<\\.?", "", x)
    return(suppressWarnings(as.numeric(z)))
  }
  suppressWarnings(as.numeric(x))
}

extract_rms_anova_p <- function(fit) {
  a <- anova(fit)
  rn <- rownames(a)
  pcol <- if ("P" %in% colnames(a)) "P" else ncol(a)
  int_p <- NA_real_
  nonlin_p <- NA_real_
  idx_int <- grep(paste0("^", TRT_COL, " \\* [^\\(]+  \\(Factor\\+Higher Order Factors\\)$"), rn)
  if (length(idx_int)) int_p <- parse_anova_p(a[idx_int[1], pcol])
  idx_nl <- grep("Nonlinear Interaction : f(A,B) vs. AB", rn, fixed = TRUE)
  if (length(idx_nl)) nonlin_p <- parse_anova_p(a[idx_nl[1], pcol])
  list(int_p = int_p, nonlin_p = nonlin_p)
}

choose_nk <- function(n) {
  RMS_RCS_NK
}

fit_boot_cph_rms <- function(df_b, surv_time_col, nk) {
  dd_b <- suppressWarnings(datadist(df_b[, c(TRT_COL, "rcsx"), drop = FALSE]))
  options(datadist = dd_b)
  fml <- as.formula(paste0(
    "Surv(", surv_time_col, ", death_status) ~ ", TRT_COL, " * rcs(rcsx, ", nk, ")"
  ))
  cph(fml, data = df_b, weights = sw, x = TRUE, y = TRUE, robust = FALSE)
}

predict_hr_curve <- function(fit, x_grid) {
  nd1 <- setNames(data.frame(1L, x_grid), c(TRT_COL, "rcsx"))
  nd0 <- setNames(data.frame(0L, x_grid), c(TRT_COL, "rcsx"))
  lp1 <- as.numeric(predict(fit, nd1, type = "lp"))
  lp0 <- as.numeric(predict(fit, nd0, type = "lp"))
  exp(lp1 - lp0)
}

plot_rms_rcs_static <- function(df_sub, var_col, var_label, surv_time_col, caption_core,
                                title_tag, log_transform = FALSE, ps_vars_base = NULL) {
  df_sub <- compute_overlap_weight_static(df_sub, var_col, ps_vars_base)
  df_sub <- df_sub[!is.na(df_sub[[var_col]]), ]
  n <- nrow(df_sub)
  cat(sprintf("\n--- [%s] %s: n=%d ---\n", title_tag, var_label, n))
  if (n < MIN_N) { cat("  SKIP: n < MIN_N\n"); return(NULL) }

  if (log_transform) {
    df_sub$rcsx <- log1p(pmax(df_sub[[var_col]], 0))
    x_vals_raw <- df_sub[[var_col]]
    x_vals <- df_sub$rcsx
    x_label_plot <- sprintf("%s (log scale)", var_label)
  } else {
    df_sub$rcsx <- df_sub[[var_col]]
    x_vals_raw <- df_sub[[var_col]]
    x_vals <- x_vals_raw
    x_label_plot <- var_label
  }

  nk <- choose_nk(n)
  dd <- suppressWarnings(datadist(df_sub[, c(TRT_COL, "rcsx"), drop = FALSE]))
  options(datadist = dd)

  tryCatch({
    fml <- as.formula(paste0(
      "Surv(", surv_time_col, ", death_status) ~ ", TRT_COL, " * rcs(rcsx, ", nk, ")"
    ))
    fit <- cph(fml, data = df_sub, weights = sw, x = TRUE, y = TRUE, robust = FALSE)
    ap <- extract_rms_anova_p(fit)
    int_p <- ap$int_p
    nonlin_p <- ap$nonlin_p

    xlim <- quantile(df_sub$rcsx, c(0.05, 0.95))
    x_grid <- seq(xlim[1], xlim[2], length.out = 200)
    hr_curve <- predict_hr_curve(fit, x_grid)

    set.seed(42)
    hr_boot <- matrix(NA_real_, N_BOOT, length(x_grid))
    for (b in seq_len(N_BOOT)) {
      idx  <- sample.int(nrow(df_sub), replace = TRUE)
      df_b <- df_sub[idx, ]
      nk_b <- choose_nk(nrow(df_b))
      tryCatch({
        fit_b <- fit_boot_cph_rms(df_b, surv_time_col, nk_b)
        hr_boot[b, ] <- predict_hr_curve(fit_b, x_grid)
      }, error = function(e) NULL)
    }
    hr_lo <- apply(hr_boot, 2, quantile, 0.025, na.rm = TRUE)
    hr_hi <- apply(hr_boot, 2, quantile, 0.975, na.rm = TRUE)
    eps_hr <- 0.05
    hr_lo[is.na(hr_lo)] <- eps_hr
    hr_hi[is.na(hr_hi)] <- pmax(10, hr_curve, na.rm = TRUE)
    hr_lo <- pmax(hr_lo, eps_hr)
    hr_hi <- pmax(hr_hi, hr_lo * 1.001)
    hr_curve <- pmax(hr_curve, eps_hr)

    HR_REF_10 <- 1
    HR_REF_09 <- 0.9
    HR_REF_08 <- 0.8
    HR_REF_07 <- 0.7
    plot_df <- data.frame(x = x_grid, hr = hr_curve, hr_lo = hr_lo, hr_hi = hr_hi)
    find_crossings <- function(df_plot, hr_ref) {
      cps <- numeric(0)
      for (i in seq_len(nrow(df_plot) - 1)) {
        y1 <- df_plot$hr[i] - hr_ref
        y2 <- df_plot$hr[i + 1] - hr_ref
        if (y1 * y2 < 0) {
          frac <- (hr_ref - df_plot$hr[i]) / (df_plot$hr[i + 1] - df_plot$hr[i])
          cps <- c(cps, df_plot$x[i] + frac * (df_plot$x[i + 1] - df_plot$x[i]))
        }
      }
      cps
    }
    cross_pts_10 <- find_crossings(plot_df, HR_REF_10)
    cross_pts_09 <- find_crossings(plot_df, HR_REF_09)
    cross_pts_08 <- find_crossings(plot_df, HR_REF_08)
    cross_pts_07 <- find_crossings(plot_df, HR_REF_07)
    y_lo_vis <- 0.25
    y_hi_vis <- 2.5

    fmt_x_cross <- function(cp) {
      if (log_transform) sprintf("%.3g", expm1(cp)) else sprintf("%.3g", cp)
    }

    # rug：与 x 轴数据同尺度；避免 NA/Inf；log 图后续 scale_x 须带 limits 否则 rug 易被裁没
    ok_rug <- is.finite(x_vals) & is.finite(xlim[1]) & is.finite(xlim[2]) &
      x_vals >= xlim[1] & x_vals <= xlim[2]
    rug_x <- x_vals[ok_rug]
    if (length(rug_x) < 1L && any(is.finite(x_vals))) {
      rug_x <- x_vals[is.finite(x_vals)]
    }

    p <- ggplot(plot_df, aes(x = x)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = y_lo_vis, ymax = 1,
               fill = "#E8F5E9", alpha = 0.35) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = y_hi_vis,
               fill = "#FFF3E0", alpha = 0.35) +
      geom_ribbon(aes(ymin = hr_lo, ymax = hr_hi), fill = "#3C5488", alpha = 0.18) +
      geom_line(aes(y = hr), color = "#3C5488", linewidth = 1.8) +
      geom_hline(yintercept = HR_REF_10, linetype = "dashed", color = "#333333", linewidth = 1.0) +
      geom_hline(yintercept = HR_REF_09, linetype = "dotdash", color = "#0072B2",
                 linewidth = 0.8, alpha = 0.8) +
      geom_hline(yintercept = HR_REF_08, linetype = "dotdash", color = "#DC0000",
                 linewidth = 0.8, alpha = 0.8) +
      geom_hline(yintercept = HR_REF_07, linetype = "dotdash", color = "#7E6148",
                 linewidth = 0.8, alpha = 0.8) +
      annotate("text", x = xlim[2] - diff(xlim) * 0.02, y = HR_REF_10,
               label = "HR=1.0", color = "#333333",
               size = 3.5, hjust = 1, vjust = -0.5, fontface = "bold") +
      annotate("text", x = xlim[2] - diff(xlim) * 0.02, y = HR_REF_09,
               label = sprintf("HR=%.1f", HR_REF_09), color = "#0072B2",
               size = 3.5, hjust = 1, vjust = 1.5, fontface = "bold") +
      annotate("text", x = xlim[2] - diff(xlim) * 0.02, y = HR_REF_08,
               label = sprintf("HR=%.1f", HR_REF_08), color = "#DC0000",
               size = 3.5, hjust = 1, vjust = -0.5, fontface = "bold") +
      annotate("text", x = xlim[2] - diff(xlim) * 0.02, y = HR_REF_07,
               label = sprintf("HR=%.1f", HR_REF_07), color = "#7E6148",
               size = 3.5, hjust = 1, vjust = 1.5, fontface = "bold")

    if (length(cross_pts_10) > 0) {
      for (cp in cross_pts_10) {
        lab <- fmt_x_cross(cp)
        p <- p +
          geom_vline(xintercept = cp, linetype = "dotted", color = "#333333",
                     linewidth = 0.7, alpha = 0.8) +
          annotate("point", x = cp, y = HR_REF_10, color = "#333333",
                   size = 3, shape = 18) +
          annotate("label", x = cp, y = 1.38, label = lab, color = "#333333",
                    size = 3.5, fontface = "bold", fill = "white",
                    label.size = 0.3, label.padding = unit(0.15, "lines"))
      }
    }
    if (length(cross_pts_09) > 0) {
      for (cp in cross_pts_09) {
        lab <- fmt_x_cross(cp)
        p <- p +
          geom_vline(xintercept = cp, linetype = "dotted", color = "#0072B2",
                     linewidth = 0.7, alpha = 0.8) +
          annotate("point", x = cp, y = HR_REF_09, color = "#0072B2",
                   size = 3, shape = 18) +
          annotate("label", x = cp, y = 0.66, label = lab, color = "#0072B2",
                    size = 3.5, fontface = "bold", fill = "white",
                    label.size = 0.3, label.padding = unit(0.15, "lines"))
      }
    }
    if (length(cross_pts_08) > 0) {
      for (cp in cross_pts_08) {
        lab <- fmt_x_cross(cp)
        p <- p +
          geom_vline(xintercept = cp, linetype = "dotted", color = "#DC0000",
                     linewidth = 0.7, alpha = 0.8) +
          annotate("point", x = cp, y = HR_REF_08, color = "#DC0000",
                   size = 3, shape = 18) +
          annotate("label", x = cp, y = 0.52, label = lab, color = "#DC0000",
                    size = 3.5, fontface = "bold", fill = "white",
                    label.size = 0.3, label.padding = unit(0.15, "lines"))
      }
    }
    if (length(cross_pts_07) > 0) {
      for (cp in cross_pts_07) {
        lab <- fmt_x_cross(cp)
        p <- p +
          geom_vline(xintercept = cp, linetype = "dotted", color = "#7E6148",
                     linewidth = 0.7, alpha = 0.8) +
          annotate("point", x = cp, y = HR_REF_07, color = "#7E6148",
                   size = 3, shape = 18) +
          annotate("label", x = cp, y = 0.40, label = lab, color = "#7E6148",
                    size = 3.3, fontface = "bold", fill = "white",
                    label.size = 0.3, label.padding = unit(0.15, "lines"))
      }
    }

    p <- p +
      geom_rug(data = data.frame(xr = rug_x),
               aes(x = xr), sides = "b", alpha = 0.42, color = "#333333",
               length = unit(2.5, "mm"), linewidth = 0.4, inherit.aes = FALSE) +
      annotate("text", x = xlim[1] + diff(xlim) * 0.03, y = 0.50,
               label = paste0("Favors\n", TRT_LABEL), color = "#2E7D32",
               size = 3.8, hjust = 0, fontface = "italic") +
      annotate("text", x = xlim[1] + diff(xlim) * 0.03, y = 1.65,
               label = paste0("Favors\n", CTRL_LABEL), color = "#E65100",
               size = 3.8, hjust = 0, fontface = "italic") +
      labs(
        title = sprintf("%s %s", title_tag, var_label),
        subtitle = sprintf(
          "rms::rcs(nk=%d) Overlap-cph | Int. p=%s | Nonlin.int. p=%s | n=%d",
          nk,
          ifelse(is.na(int_p), "NA", sprintf("%.3f", int_p)),
          ifelse(is.na(nonlin_p), "NA", sprintf("%.3f", nonlin_p)), n),
        x = x_label_plot,
        y = sprintf("HR (%s vs %s)", TRT_LABEL, CTRL_LABEL),
        caption = paste0(caption_core,
          " | Harrell RCS | P5\u2013P95 | 95% bootstrap CI | Gray: HR=1 | Blue: 0.9 | Red: 0.8 | Brown: 0.7")
      ) +
      scale_y_log10(breaks = c(0.3, 0.5, 0.7, 0.9, 1.0, 1.5, 2.0),
                    labels = c("0.3", "0.5", "0.7", "0.9", "1.0", "1.5", "2.0")) +
      coord_cartesian(ylim = c(0.25, 2.5), xlim = xlim) +
      theme_bw(base_size = 13)

    if (log_transform) {
      raw_range <- expm1(xlim)
      candidate_ticks <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000,
                           10000, 20000, 50000, 100000, 200000, 500000)
      raw_ticks <- candidate_ticks[candidate_ticks >= raw_range[1] & candidate_ticks <= raw_range[2]]
      if (length(raw_ticks) < 4) {
        raw_ticks <- pretty(raw_range, n = 6)
        raw_ticks <- raw_ticks[raw_ticks > 0]
      }
      if (length(raw_ticks) > 8) {
        raw_ticks <- raw_ticks[seq(1, length(raw_ticks), length.out = min(8, length(raw_ticks)))]
      }
      log_ticks <- log1p(raw_ticks)
      fmt_label <- function(v) {
        ifelse(v >= 1000, sprintf("%gK", v / 1000),
               ifelse(v >= 1, sprintf("%g", v), sprintf("%.1f", v)))
      }
      p <- p + scale_x_continuous(
        breaks = log_ticks, labels = fmt_label(raw_ticks),
        name = x_label_plot,
        limits = xlim,
        expand = ggplot2::expansion(mult = c(0.02, 0.02))
      )
    } else {
      # 非 log：须显式 limits+expand，否则 rug 常与默认 x scale 冲突而被裁没（如 mono、NLR）
      p <- p + scale_x_continuous(
        limits = xlim,
        expand = ggplot2::expansion(mult = c(0.04, 0.04))
      )
    }

    cross_pts10_raw <- if (log_transform) expm1(cross_pts_10) else cross_pts_10
    cross_pts09_raw <- if (log_transform) expm1(cross_pts_09) else cross_pts_09
    cross_pts08_raw <- if (log_transform) expm1(cross_pts_08) else cross_pts_08
    cross_pts07_raw <- if (log_transform) expm1(cross_pts_07) else cross_pts_07
    hr_first <- hr_curve[1]
    hr_last  <- hr_curve[length(hr_curve)]
    xlim_raw <- if (log_transform) expm1(xlim) else xlim

    return(list(plot = p, int_p = int_p, nonlin_p = nonlin_p, n = n, nk = nk,
                cross_pts10_raw = cross_pts10_raw, cross_pts09_raw = cross_pts09_raw,
                cross_pts08_raw = cross_pts08_raw, cross_pts07_raw = cross_pts07_raw,
                hr_first = hr_first, hr_last = hr_last, xlim_raw = xlim_raw))
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    return(NULL)
  })
}

plot_rms_rcs_dynamic <- function(df_sub, var_col, var_label, surv_time_col, caption_core,
                                 title_tag, ps_vars_base = NULL) {
  df_sub <- compute_overlap_weight_dynamic(df_sub, var_col, ps_vars_base)
  df_sub <- df_sub[!is.na(df_sub[[var_col]]), ]
  n <- nrow(df_sub)
  cat(sprintf("\n--- [%s] %s: n=%d ---\n", title_tag, var_label, n))
  if (n < MIN_N) { cat("  SKIP: n < MIN_N\n"); return(NULL) }

  df_sub$rcsx <- df_sub[[var_col]]
  x_vals <- df_sub$rcsx
  x_label_plot <- sprintf("%s (%%)", var_label)
  nk <- choose_nk(n)
  dd <- suppressWarnings(datadist(df_sub[, c(TRT_COL, "rcsx"), drop = FALSE]))
  options(datadist = dd)
  xlim <- quantile(x_vals, c(0.05, 0.95))

  tryCatch({
    fml <- as.formula(paste0(
      "Surv(", surv_time_col, ", death_status) ~ ", TRT_COL, " * rcs(rcsx, ", nk, ")"
    ))
    fit <- cph(fml, data = df_sub, weights = sw, x = TRUE, y = TRUE, robust = FALSE)
    ap <- extract_rms_anova_p(fit)
    int_p <- ap$int_p
    nonlin_p <- ap$nonlin_p

    x_grid <- seq(xlim[1], xlim[2], length.out = 200)
    hr_curve <- predict_hr_curve(fit, x_grid)

    set.seed(42)
    hr_boot <- matrix(NA_real_, N_BOOT, length(x_grid))
    for (b in seq_len(N_BOOT)) {
      idx  <- sample.int(nrow(df_sub), replace = TRUE)
      df_b <- df_sub[idx, ]
      nk_b <- choose_nk(nrow(df_b))
      tryCatch({
        fit_b <- fit_boot_cph_rms(df_b, surv_time_col, nk_b)
        hr_boot[b, ] <- predict_hr_curve(fit_b, x_grid)
      }, error = function(e) NULL)
    }
    hr_lo <- apply(hr_boot, 2, quantile, 0.025, na.rm = TRUE)
    hr_hi <- apply(hr_boot, 2, quantile, 0.975, na.rm = TRUE)
    eps_hr <- 0.05
    hr_lo[is.na(hr_lo)] <- eps_hr
    hr_hi[is.na(hr_hi)] <- pmax(10, hr_curve, na.rm = TRUE)
    hr_lo <- pmax(hr_lo, eps_hr)
    hr_hi <- pmax(hr_hi, hr_lo * 1.001)
    hr_curve <- pmax(hr_curve, eps_hr)

    HR_REF_10 <- 1
    HR_REF_09 <- 0.9
    HR_REF_08 <- 0.8
    HR_REF_07 <- 0.7
    plot_df <- data.frame(x = x_grid, hr = hr_curve, hr_lo = hr_lo, hr_hi = hr_hi)
    find_crossings <- function(df_plot, hr_ref) {
      cps <- numeric(0)
      for (i in seq_len(nrow(df_plot) - 1)) {
        y1 <- df_plot$hr[i] - hr_ref
        y2 <- df_plot$hr[i + 1] - hr_ref
        if (y1 * y2 < 0) {
          frac <- (hr_ref - df_plot$hr[i]) / (df_plot$hr[i + 1] - df_plot$hr[i])
          cps <- c(cps, df_plot$x[i] + frac * (df_plot$x[i + 1] - df_plot$x[i]))
        }
      }
      cps
    }
    cross_pts_10 <- find_crossings(plot_df, HR_REF_10)
    cross_pts_09 <- find_crossings(plot_df, HR_REF_09)
    cross_pts_08 <- find_crossings(plot_df, HR_REF_08)
    cross_pts_07 <- find_crossings(plot_df, HR_REF_07)
    y_lo_vis <- 0.25
    y_hi_vis <- 2.5

    cp_disp <- function(cp) {
      ifelse(abs(cp) >= 100, sprintf("%.0f%%", cp),
        ifelse(abs(cp) >= 10, sprintf("%.1f%%", cp), sprintf("%.1f%%", cp)))
    }

    ok_rug_d <- is.finite(x_vals) & is.finite(xlim[1]) & is.finite(xlim[2]) &
      x_vals >= xlim[1] & x_vals <= xlim[2]
    rug_x_d <- x_vals[ok_rug_d]
    if (length(rug_x_d) < 1L && any(is.finite(x_vals))) {
      rug_x_d <- x_vals[is.finite(x_vals)]
    }

    p <- ggplot(plot_df, aes(x = x)) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = y_lo_vis, ymax = 1,
               fill = "#E8F5E9", alpha = 0.35) +
      annotate("rect", xmin = -Inf, xmax = Inf, ymin = 1, ymax = y_hi_vis,
               fill = "#FFF3E0", alpha = 0.35) +
      geom_ribbon(aes(ymin = hr_lo, ymax = hr_hi), fill = "#3C5488", alpha = 0.18) +
      geom_line(aes(y = hr), color = "#3C5488", linewidth = 1.8) +
      geom_hline(yintercept = HR_REF_10, linetype = "dashed", color = "#333333", linewidth = 1.0) +
      geom_hline(yintercept = HR_REF_09, linetype = "dotdash", color = "#0072B2",
                 linewidth = 0.8, alpha = 0.8) +
      geom_hline(yintercept = HR_REF_08, linetype = "dotdash", color = "#DC0000",
                 linewidth = 0.8, alpha = 0.8) +
      geom_hline(yintercept = HR_REF_07, linetype = "dotdash", color = "#7E6148",
                 linewidth = 0.8, alpha = 0.8) +
      annotate("text", x = xlim[2] - diff(xlim) * 0.02, y = HR_REF_10,
               label = "HR=1.0", color = "#333333",
               size = 3.5, hjust = 1, vjust = -0.5, fontface = "bold") +
      annotate("text", x = xlim[2] - diff(xlim) * 0.02, y = HR_REF_09,
               label = sprintf("HR=%.1f", HR_REF_09), color = "#0072B2",
               size = 3.5, hjust = 1, vjust = 1.5, fontface = "bold") +
      annotate("text", x = xlim[2] - diff(xlim) * 0.02, y = HR_REF_08,
               label = sprintf("HR=%.1f", HR_REF_08), color = "#DC0000",
               size = 3.5, hjust = 1, vjust = -0.5, fontface = "bold") +
      annotate("text", x = xlim[2] - diff(xlim) * 0.02, y = HR_REF_07,
               label = sprintf("HR=%.1f", HR_REF_07), color = "#7E6148",
               size = 3.5, hjust = 1, vjust = 1.5, fontface = "bold")

    if (length(cross_pts_10) > 0) {
      for (cp in cross_pts_10) {
        p <- p +
          geom_vline(xintercept = cp, linetype = "dotted", color = "#333333",
                     linewidth = 0.7, alpha = 0.8) +
          annotate("point", x = cp, y = HR_REF_10, color = "#333333",
                   size = 3, shape = 18) +
          annotate("label", x = cp, y = 1.38, label = cp_disp(cp), color = "#333333",
                    size = 3.5, fontface = "bold", fill = "white",
                    label.size = 0.3, label.padding = unit(0.15, "lines"))
      }
    }
    if (length(cross_pts_09) > 0) {
      for (cp in cross_pts_09) {
        p <- p +
          geom_vline(xintercept = cp, linetype = "dotted", color = "#0072B2",
                     linewidth = 0.7, alpha = 0.8) +
          annotate("point", x = cp, y = HR_REF_09, color = "#0072B2",
                   size = 3, shape = 18) +
          annotate("label", x = cp, y = 0.66, label = cp_disp(cp), color = "#0072B2",
                    size = 3.5, fontface = "bold", fill = "white",
                    label.size = 0.3, label.padding = unit(0.15, "lines"))
      }
    }
    if (length(cross_pts_08) > 0) {
      for (cp in cross_pts_08) {
        p <- p +
          geom_vline(xintercept = cp, linetype = "dotted", color = "#DC0000",
                     linewidth = 0.7, alpha = 0.8) +
          annotate("point", x = cp, y = HR_REF_08, color = "#DC0000",
                   size = 3, shape = 18) +
          annotate("label", x = cp, y = 0.52, label = cp_disp(cp), color = "#DC0000",
                    size = 3.5, fontface = "bold", fill = "white",
                    label.size = 0.3, label.padding = unit(0.15, "lines"))
      }
    }
    if (length(cross_pts_07) > 0) {
      for (cp in cross_pts_07) {
        p <- p +
          geom_vline(xintercept = cp, linetype = "dotted", color = "#7E6148",
                     linewidth = 0.7, alpha = 0.8) +
          annotate("point", x = cp, y = HR_REF_07, color = "#7E6148",
                   size = 3, shape = 18) +
          annotate("label", x = cp, y = 0.40, label = cp_disp(cp), color = "#7E6148",
                    size = 3.3, fontface = "bold", fill = "white",
                    label.size = 0.3, label.padding = unit(0.15, "lines"))
      }
    }

    p <- p +
      geom_vline(xintercept = 0, linetype = "dotted", color = "#888888", linewidth = 0.6) +
      geom_rug(data = data.frame(xr = rug_x_d),
               aes(x = xr), sides = "b", alpha = 0.42, color = "#333333",
               length = unit(2.5, "mm"), linewidth = 0.4, inherit.aes = FALSE) +
      annotate("text", x = xlim[1] + diff(xlim) * 0.03, y = 0.50,
               label = paste0("Favors\n", TRT_LABEL), color = "#2E7D32",
               size = 3.8, hjust = 0, fontface = "italic") +
      annotate("text", x = xlim[1] + diff(xlim) * 0.03, y = 1.65,
               label = paste0("Favors\n", CTRL_LABEL), color = "#E65100",
               size = 3.8, hjust = 0, fontface = "italic") +
      annotate("text", x = 0, y = 0.28, label = "No change",
               color = "#888888", size = 3.5, hjust = 0.5, fontface = "italic") +
      labs(
        title = sprintf("%s %s", title_tag, var_label),
        subtitle = sprintf(
          "rms::rcs(nk=%d) Overlap-cph | Int. p=%s | Nonlin.int. p=%s | n=%d",
          nk,
          ifelse(is.na(int_p), "NA", sprintf("%.3f", int_p)),
          ifelse(is.na(nonlin_p), "NA", sprintf("%.3f", nonlin_p)), n),
        x = x_label_plot, y = sprintf("HR (%s vs %s)", TRT_LABEL, CTRL_LABEL),
        caption = paste0(caption_core,
          " | Harrell RCS | P5\u2013P95 | 95% bootstrap CI | Gray: HR=1 | Blue: 0.9 | Red: 0.8 | Brown: 0.7")
      ) +
      scale_y_log10(breaks = c(0.3, 0.5, 0.7, 0.9, 1.0, 1.5, 2.0),
                    labels = c("0.3", "0.5", "0.7", "0.9", "1.0", "1.5", "2.0")) +
      coord_cartesian(ylim = c(0.25, 2.5), xlim = xlim) +
      theme_bw(base_size = 13)

    p <- p + scale_x_continuous(
      limits = xlim,
      expand = ggplot2::expansion(mult = c(0.04, 0.04))
    )

    list(plot = p, int_p = int_p, nonlin_p = nonlin_p, n = n, nk = nk,
         cross_pts10_raw = cross_pts_10, cross_pts09_raw = cross_pts_09,
         cross_pts08_raw = cross_pts_08, cross_pts07_raw = cross_pts_07,
         hr_first = hr_curve[1], hr_last = hr_curve[length(hr_curve)], xlim_raw = xlim)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    return(NULL)
  })
}

vars_config_static <- list(
  list(col = "afp", label = "AFP (baseline, ng/mL)", log_transform = TRUE,  group = "baseline"),
  list(col = "pivka", label = "PIVKA-II (baseline, mAU/mL)", log_transform = TRUE,  group = "baseline"),
  list(col = "nlr_bl", label = "NLR (baseline)", log_transform = FALSE, group = "baseline"),
  list(col = "albi_bl", label = "ALBI score (baseline)", log_transform = FALSE, group = "baseline"),
  list(col = "alb_bl", label = "Albumin (baseline, g/L)", log_transform = FALSE, group = "baseline"),
  list(col = "tbil_bl", label = "Total bilirubin (baseline, \u03bcmol/L)", log_transform = TRUE,  group = "baseline"),
  list(col = "plr_bl", label = "PLR (baseline)", log_transform = FALSE, group = "baseline"),
  list(col = "sii_bl", label = "SII (PLT\u00d7NLR, baseline)", log_transform = TRUE,  group = "baseline"),
  list(col = "piv_bl", label = "PIV (Mono\u00d7PLT\u00d7NLR, baseline)", log_transform = TRUE,  group = "baseline"),
  list(col = "neut_bl", label = "Neutrophil (baseline, 10^9/L)", log_transform = FALSE, group = "baseline"),
  list(col = "lymph_bl", label = "Lymphocyte (baseline, 10^9/L)", log_transform = FALSE, group = "baseline"),
  list(col = "mono_bl", label = "Monocyte (baseline, 10^9/L)", log_transform = FALSE, group = "baseline"),
  list(col = "plt", label = "Platelet (baseline, 10^9/L)", log_transform = TRUE,  group = "baseline"),
  list(col = "tumor_max_diameter_cm", label = "Tumor diameter (baseline, cm)", log_transform = FALSE, group = "baseline"),
  list(col = "afp_pre3", label = "AFP (pre-HAIC-3, ng/mL)", log_transform = TRUE,  group = "pre3_static"),
  list(col = "pivka_pre3", label = "PIVKA-II (pre-HAIC-3, mAU/mL)", log_transform = TRUE,  group = "pre3_static"),
  list(col = "nlr_pre3", label = "NLR (pre-HAIC-3)", log_transform = FALSE, group = "pre3_static"),
  list(col = "albi_pre3", label = "ALBI score (pre-HAIC-3)", log_transform = FALSE, group = "pre3_static"),
  list(col = "alb_pre3", label = "Albumin (pre-HAIC-3, g/L)", log_transform = FALSE, group = "pre3_static"),
  list(col = "tbil_pre3", label = "Total bilirubin (pre-HAIC-3, \u03bcmol/L)", log_transform = TRUE,  group = "pre3_static"),
  list(col = "plr_pre3", label = "PLR (pre-HAIC-3)", log_transform = FALSE, group = "pre3_static"),
  list(col = "mono_pre3", label = "Monocyte (pre-HAIC-3, 10^9/L)", log_transform = FALSE, group = "pre3_static"),
  list(col = "plt_pre3", label = "Platelet (pre-HAIC-3, 10^9/L)", log_transform = TRUE,  group = "pre3_static"),
  list(col = "sii_pre3", label = "SII (PLT\u00d7NLR, pre-HAIC-3)", log_transform = TRUE,  group = "pre3_static"),
  list(col = "piv_pre3", label = "PIV (Mono\u00d7PLT\u00d7NLR, pre-HAIC-3)", log_transform = TRUE,  group = "pre3_static"),
  list(col = "neut_pre3", label = "Neutrophil (pre-HAIC-3, 10^9/L)", log_transform = FALSE, group = "pre3_static"),
  list(col = "lymph_pre3", label = "Lymphocyte (pre-HAIC-3, 10^9/L)", log_transform = FALSE, group = "pre3_static")
)

vars_config_dynamic <- list(
  list(col = "nlr_change_pre3", label = "NLR change rate"),
  list(col = "sii_change_pre3", label = "SII change rate"),
  list(col = "piv_change_pre3", label = "PIV change rate"),
  list(col = "afp_change_pre3", label = "AFP change rate"),
  list(col = "pivka_change_pre3", label = "PIVKA-II change rate"),
  list(col = "neut_change_pre3", label = "Neutrophil change rate"),
  list(col = "plr_change_pre3", label = "PLR change rate"),
  list(col = "albi_change_pre3", label = "ALBI change rate"),
  list(col = "alb_change_pre3", label = "Albumin change rate"),
  list(col = "tbil_change_pre3", label = "TBIL change rate"),
  list(col = "lymph_change_pre3", label = "Lymphocyte change rate"),
  list(col = "mono_change_pre3", label = "Monocyte change rate"),
  list(col = "plt_change_pre3", label = "Platelet change rate")
)

run_one_arm <- function(df_arm, arm_folder, surv_time_col, caption_core, title_tag) {
  OUT_DIR <- file.path(BASE_OUT, arm_folder)
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("ARM:", arm_folder, "| Surv(", surv_time_col, ", death_status) | n =", nrow(df_arm), "\n")
  cat(strrep("=", 70), "\n")

  ps_vars_full <- c("albi_bl", "alb_bl", "tbil_bl", "inr", "plt",
                    "tumor_max_diameter_cm", "tumor_count_enc",
                    "pvtt_grade", "ascites_score_enc",
                    "log_afp_bl", "log_pivka_bl",
                    "metastasis_binary", "lymph_node_binary",
                    "neut_bl", "lymph_bl", "mono_bl")
  df <- impute_ps_vars(df_arm, ps_vars_full)
  ps_vars_available <- ps_vars_full[ps_vars_full %in% names(df)]
  ps_vars_available <- ps_vars_available[sapply(ps_vars_available,
    function(v) sd(df[[v]], na.rm = TRUE) > 0)]

  results_static <- list()
  plots_static   <- list()
  for (cfg in vars_config_static) {
    if (!(cfg$col %in% names(df))) next
    res <- plot_rms_rcs_static(df, cfg$col, cfg$label, surv_time_col, caption_core, title_tag,
                               log_transform = cfg$log_transform, ps_vars_base = ps_vars_available)
    if (!is.null(res)) {
      results_static[[cfg$col]] <- c(res, list(label = cfg$label, group = cfg$group))
      plots_static[[cfg$col]]   <- res$plot
    }
  }

  bl_keys <- names(plots_static)[sapply(names(plots_static),
    function(k) results_static[[k]]$group == "baseline")]
  p3_keys <- names(plots_static)[sapply(names(plots_static),
    function(k) results_static[[k]]$group == "pre3_static")]
  if (length(bl_keys) >= 2) {
    comb <- gridExtra::grid.arrange(
      grobs = plots_static[bl_keys], ncol = 2,
      top = grid::textGrob(
        paste(title_tag, sprintf("Harrell RCS — Baseline (%s vs %s, rms)", CTRL_LABEL, TRT_LABEL)),
        gp = grid::gpar(fontsize = 14, fontface = "bold")))
    h_bl <- ceiling(length(bl_keys) / 2) * 4.2
    ggsave(file.path(OUT_DIR, "rms_rcs_static_baseline_combined.pdf"), comb, width = 11.5, height = h_bl, device = "pdf")
    ggsave(file.path(OUT_DIR, "rms_rcs_static_baseline_combined.png"), comb, width = 11.5, height = h_bl, dpi = 300)
  }
  if (length(p3_keys) >= 2) {
    comb <- gridExtra::grid.arrange(
      grobs = plots_static[p3_keys], ncol = 2,
      top = grid::textGrob(
        paste(title_tag, sprintf("Harrell RCS — Pre-HAIC-3 (%s vs %s, rms)", CTRL_LABEL, TRT_LABEL)),
        gp = grid::gpar(fontsize = 14, fontface = "bold")))
    h_p3 <- ceiling(length(p3_keys) / 2) * 4.2
    ggsave(file.path(OUT_DIR, "rms_rcs_static_pre3_combined.pdf"), comb, width = 11.5, height = h_p3, device = "pdf")
    ggsave(file.path(OUT_DIR, "rms_rcs_static_pre3_combined.png"), comb, width = 11.5, height = h_p3, dpi = 300)
  }

  if (length(names(results_static)) > 0) {
    summary_table <- do.call(rbind, lapply(names(results_static), function(k) {
      r <- results_static[[k]]
      safe <- function(x) ifelse(is.null(x) || length(x) == 0 || is.na(x[1]), NA_real_, x[1])
      data.frame(variable = k, label = r$label, group = r$group, n = r$n, nk = r$nk,
                 interaction_p = round(safe(r$int_p), 4), nonlinear_interaction_p = round(safe(r$nonlin_p), 4),
                 stringsAsFactors = FALSE)
    }))
    write.csv(summary_table, file.path(OUT_DIR, "rms_rcs_static_anova_summary.csv"), row.names = FALSE)
    cross_table <- do.call(rbind, lapply(names(results_static), function(k) {
      r <- results_static[[k]]
      cp10 <- r$cross_pts10_raw; cp09 <- r$cross_pts09_raw
      cp08 <- r$cross_pts08_raw; cp07 <- r$cross_pts07_raw
      data.frame(variable = k, label = r$label, group = r$group,
                 n_crosspoints_hr10 = length(cp10),
                 cross1_hr10 = ifelse(length(cp10) >= 1, round(cp10[1], 4), NA_real_),
                 cross2_hr10 = ifelse(length(cp10) >= 2, round(cp10[2], 4), NA_real_),
                 n_crosspoints_hr09 = length(cp09),
                 cross1_hr09 = ifelse(length(cp09) >= 1, round(cp09[1], 4), NA_real_),
                 cross2_hr09 = ifelse(length(cp09) >= 2, round(cp09[2], 4), NA_real_),
                 n_crosspoints_hr08 = length(cp08),
                 cross1_hr08 = ifelse(length(cp08) >= 1, round(cp08[1], 4), NA_real_),
                 cross2_hr08 = ifelse(length(cp08) >= 2, round(cp08[2], 4), NA_real_),
                 n_crosspoints_hr07 = length(cp07),
                 cross1_hr07 = ifelse(length(cp07) >= 1, round(cp07[1], 4), NA_real_),
                 cross2_hr07 = ifelse(length(cp07) >= 2, round(cp07[2], 4), NA_real_),
                 hr_at_xmin = round(r$hr_first, 4), hr_at_xmax = round(r$hr_last, 4),
                 xmin_raw = round(r$xlim_raw[1], 4), xmax_raw = round(r$xlim_raw[2], 4),
                 stringsAsFactors = FALSE)
    }))
    write.csv(cross_table, file.path(OUT_DIR, "rms_rcs_static_hr07_hr08_crosspoints.csv"), row.names = FALSE)
  }

  results_dyn <- list()
  plots_dyn   <- list()
  for (cfg in vars_config_dynamic) {
    if (!(cfg$col %in% names(df))) next
    res <- plot_rms_rcs_dynamic(df, cfg$col, cfg$label, surv_time_col, caption_core, title_tag,
                                 ps_vars_base = ps_vars_available)
    if (!is.null(res)) {
      results_dyn[[cfg$col]] <- c(res, list(label = cfg$label))
      plots_dyn[[cfg$col]]   <- res$plot
    }
  }
  if (length(plots_dyn) >= 2) {
    dyn_plot_order <- vapply(vars_config_dynamic, function(z) z$col, character(1))
    dyn_plot_order <- dyn_plot_order[dyn_plot_order %in% names(plots_dyn)]
    comb <- gridExtra::grid.arrange(
      grobs = plots_dyn[dyn_plot_order], ncol = 2,
      top = grid::textGrob(
        paste(title_tag, sprintf("Harrell RCS — Dynamic change %% (%s vs %s, rms)", CTRL_LABEL, TRT_LABEL)),
        gp = grid::gpar(fontsize = 14, fontface = "bold")))
    h <- ceiling(length(plots_dyn) / 2) * 4.2
    ggsave(file.path(OUT_DIR, "rms_rcs_dynamic_combined.pdf"), comb, width = 11.5, height = h, device = "pdf")
    ggsave(file.path(OUT_DIR, "rms_rcs_dynamic_combined.png"), comb, width = 11.5, height = h, dpi = 300)
  }
  if (length(names(results_dyn)) > 0) {
    summary_dyn <- do.call(rbind, lapply(names(results_dyn), function(k) {
      r <- results_dyn[[k]]
      safe <- function(x) ifelse(is.null(x) || length(x) == 0 || is.na(x[1]), NA_real_, x[1])
      data.frame(variable = k, label = r$label, n = r$n, nk = r$nk,
                 interaction_p = round(safe(r$int_p), 4), nonlinear_interaction_p = round(safe(r$nonlin_p), 4),
                 stringsAsFactors = FALSE)
    }))
    write.csv(summary_dyn, file.path(OUT_DIR, "rms_rcs_dynamic_anova_summary.csv"), row.names = FALSE)
    cross_dyn <- do.call(rbind, lapply(names(results_dyn), function(k) {
      r <- results_dyn[[k]]
      cp10 <- r$cross_pts10_raw; cp09 <- r$cross_pts09_raw
      cp08 <- r$cross_pts08_raw; cp07 <- r$cross_pts07_raw
      data.frame(variable = k, label = r$label, group = "dynamic",
                 n_crosspoints_hr10 = length(cp10),
                 cross1_hr10 = ifelse(length(cp10) >= 1, round(cp10[1], 4), NA_real_),
                 cross2_hr10 = ifelse(length(cp10) >= 2, round(cp10[2], 4), NA_real_),
                 n_crosspoints_hr09 = length(cp09),
                 cross1_hr09 = ifelse(length(cp09) >= 1, round(cp09[1], 4), NA_real_),
                 cross2_hr09 = ifelse(length(cp09) >= 2, round(cp09[2], 4), NA_real_),
                 n_crosspoints_hr08 = length(cp08),
                 cross1_hr08 = ifelse(length(cp08) >= 1, round(cp08[1], 4), NA_real_),
                 cross2_hr08 = ifelse(length(cp08) >= 2, round(cp08[2], 4), NA_real_),
                 n_crosspoints_hr07 = length(cp07),
                 cross1_hr07 = ifelse(length(cp07) >= 1, round(cp07[1], 4), NA_real_),
                 cross2_hr07 = ifelse(length(cp07) >= 2, round(cp07[2], 4), NA_real_),
                 hr_at_xmin = round(r$hr_first, 4), hr_at_xmax = round(r$hr_last, 4),
                 xmin_raw = round(r$xlim_raw[1], 4), xmax_raw = round(r$xlim_raw[2], 4),
                 stringsAsFactors = FALSE)
    }))
    write.csv(cross_dyn, file.path(OUT_DIR, "rms_rcs_dynamic_hr07_hr08_crosspoints.csv"), row.names = FALSE)
  }
  cat("Saved outputs under:", OUT_DIR, "\n")
}

# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════
if (!file.exists(DATA_CSV)) {
  stop("Data file not found: ", DATA_CSV,
       "\n先在本目录运行: python3 RCS_INT_HAIC_ALONE_AND_I_CONC_build_cohort.py\n或设置 RMS_RCS_DATA_CSV",
       call. = FALSE)
}

df0 <- read.csv(DATA_CSV, stringsAsFactors = FALSE)
cat("Loaded:", DATA_CSV, "| nrow =", nrow(df0), "\n")
if (!TRT_COL %in% names(df0)) {
  stop(paste0("Column '", TRT_COL, "' missing in data. Check RMS_RCS_TRT_COL or rebuild cohort CSV."), call. = FALSE)
}

if (!"sii_bl" %in% names(df0) && all(c("plt", "nlr_bl") %in% names(df0)))
  df0$sii_bl <- df0$plt * df0$nlr_bl
if (!"piv_bl" %in% names(df0) && all(c("mono_bl", "plt", "nlr_bl") %in% names(df0)))
  df0$piv_bl <- df0$mono_bl * df0$plt * df0$nlr_bl
if (!"sii_pre3" %in% names(df0) && all(c("plt_pre3", "nlr_pre3") %in% names(df0)))
  df0$sii_pre3 <- df0$plt_pre3 * df0$nlr_pre3
if (!"piv_pre3" %in% names(df0) && all(c("mono_pre3", "plt_pre3", "nlr_pre3") %in% names(df0)))
  df0$piv_pre3 <- df0$mono_pre3 * df0$plt_pre3 * df0$nlr_pre3

if (!"os_lm" %in% names(df0)) {
  df0$os_lm <- df0$os_months - 42 / 30.44
}

df_landmark <- df0[!is.na(df0$os_lm) & df0$os_lm > 0, ]
df_total    <- df0[!is.na(df0$os_months) & df0$os_months > 0, ]

cat("Landmark subset (os_lm > 0):", nrow(df_landmark), "rows\n")
cat("Total OS subset (os_months > 0):", nrow(df_total), "rows\n")

run_one_arm(
  df_landmark, "landmark", "os_lm",
  caption_core = "42-day landmark residual OS (primary time axis)",
  title_tag = "[Landmark]"
)
run_one_arm(
  df_total, "total_os", "os_months",
  caption_core = "Total OS from baseline (sensitivity)",
  title_tag = "[Total OS]"
)

cat("\n=== All arms completed. Output root:", BASE_OUT, "===\n")
