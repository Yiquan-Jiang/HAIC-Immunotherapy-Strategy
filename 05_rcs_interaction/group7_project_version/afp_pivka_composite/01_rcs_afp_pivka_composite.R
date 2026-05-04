#!/usr/bin/env Rscript
# =============================================================================
# AFP-PIVKA 组合指标 RCS 非线性交互分析
#
# 参数化脚本：通过命令行参数或环境变量选择队列
#   Rscript 01_rcs_afp_pivka_composite.R THEN_IT
#   Rscript 01_rcs_afp_pivka_composite.R THEN_I
#   Rscript 01_rcs_afp_pivka_composite.R THEN_T
#   Rscript 01_rcs_afp_pivka_composite.R T_CONC
#   Rscript 01_rcs_afp_pivka_composite.R ALL    # 全部运行
#
# 模型: Surv(time, death_status) ~ trt * rcs(composite_var, nk)
# 方法: IPTW 加权 Cox (rms::cph) + Bootstrap 95% CI
# 时间尺度: landmark (42-day) + total_os
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
  I_T_CONC = list(
    csv = "composite_I_T_CONC_cohort.csv",
    trt_col = "trt_haic_i_t_conc",
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

cat("SCRIPT_DIR:", SCRIPT_DIR, "\nRMS_RCS_NK:", RMS_RCS_NK,
    "| N_BOOT:", N_BOOT, "| cohort_key:", cohort_key, "\n\n")

# ── IPTW ─────────────────────────────────────────────────────────────────────
impute_ps_vars <- function(df, ps_vars_full) {
  for (v in ps_vars_full) {
    if (v %in% names(df) && any(is.na(df[[v]]))) {
      med_val <- median(df[[v]], na.rm = TRUE)
      df[[v]][is.na(df[[v]])] <- med_val
    }
  }
  df
}

compute_iptw <- function(df_input, var_col, trt_col, ps_vars_base) {
  # 排除映射：
  #   - 静态组合指标：排除原始基线 log 值 + 排除所有动态变量（变化率）
  #   - 动态组合指标：仅排除原始基线 log 值，保留所有变化率作为平衡变量
  #   - 单独变化率：排除自身，保留其他变化率
  exclude_map <- list(
    # 静态组合指标：排除原始基线 log 值 + 排除所有动态变量
    log_afp_pivka_ratio_bl  = c("log_afp_bl", "log_pivka_bl", "afp_change_pre3", "pivka_change_pre3", "albi_change_pre3"),
    log_afp_pivka_sum_bl    = c("log_afp_bl", "log_pivka_bl", "afp_change_pre3", "pivka_change_pre3", "albi_change_pre3"),
    afp_pivka_pc1_bl        = c("log_afp_bl", "log_pivka_bl", "afp_change_pre3", "pivka_change_pre3", "albi_change_pre3"),
    afp_pivka_pc2_bl        = c("log_afp_bl", "log_pivka_bl", "afp_change_pre3", "pivka_change_pre3", "albi_change_pre3"),
    log_afp_pivka_ratio_pre3 = c("log_afp_bl", "log_pivka_bl", "afp_change_pre3", "pivka_change_pre3", "albi_change_pre3"),
    log_afp_pivka_sum_pre3   = c("log_afp_bl", "log_pivka_bl", "afp_change_pre3", "pivka_change_pre3", "albi_change_pre3"),
    afp_pivka_pc1_pre3       = c("log_afp_bl", "log_pivka_bl", "afp_change_pre3", "pivka_change_pre3", "albi_change_pre3"),
    afp_pivka_pc2_pre3       = c("log_afp_bl", "log_pivka_bl", "afp_change_pre3", "pivka_change_pre3", "albi_change_pre3"),
    # 动态组合指标：仅排除原始基线 log 值，保留所有变化率
    afp_pivka_change_diff   = c("log_afp_bl", "log_pivka_bl"),
    afp_pivka_change_sum    = c("log_afp_bl", "log_pivka_bl"),
    afp_pivka_pc1_dyn       = c("log_afp_bl", "log_pivka_bl"),
    afp_pivka_pc2_dyn       = c("log_afp_bl", "log_pivka_bl"),
    # 单独变化率：仅排除自身（保留其他变化率作为平衡变量）
    afp_change_pre3         = c("afp_change_pre3"),
    pivka_change_pre3       = c("pivka_change_pre3"),
    albi_change_pre3        = c("albi_change_pre3")
  )
  to_exclude <- if (var_col %in% names(exclude_map)) exclude_map[[var_col]] else c()
  ps_vars_use <- setdiff(ps_vars_base, to_exclude)
  ps_vars_use <- ps_vars_use[ps_vars_use %in% names(df_input)]
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

# ── 核心绘图函数 ────────────────────────────────────────────────────────────
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

plot_rcs_composite <- function(df_sub, var_col, var_label, surv_time_col, caption_core,
                                title_tag, trt_col, trt_label, ctrl_label,
                                ps_vars_base, is_dynamic = FALSE) {
  df_sub <- compute_iptw(df_sub, var_col, trt_col, ps_vars_base)
  df_sub <- df_sub[!is.na(df_sub[[var_col]]), ]
  n <- nrow(df_sub)
  cat(sprintf("\n--- [%s] %s: n=%d ---\n", title_tag, var_label, n))
  if (n < MIN_N) { cat("  SKIP: n < MIN_N\n"); return(NULL) }

  df_sub$rcsx <- df_sub[[var_col]]
  x_vals <- df_sub$rcsx
  nk <- RMS_RCS_NK
  dd <- suppressWarnings(datadist(df_sub[, c(trt_col, "rcsx"), drop = FALSE]))
  options(datadist = dd)

  tryCatch({
    fml <- as.formula(paste0(
      "Surv(", surv_time_col, ", death_status) ~ ", trt_col, " * rcs(rcsx, ", nk, ")"
    ))
    fit <- cph(fml, data = df_sub, weights = sw, x = TRUE, y = TRUE, robust = FALSE)
    ap <- extract_rms_anova_p(fit, trt_col)
    int_p <- ap$int_p
    nonlin_p <- ap$nonlin_p

    xlim <- quantile(df_sub$rcsx, c(0.05, 0.95))
    x_grid <- seq(xlim[1], xlim[2], length.out = 200)
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

    ok_rug <- is.finite(x_vals) & x_vals >= xlim[1] & x_vals <= xlim[2]
    rug_x <- x_vals[ok_rug]
    if (length(rug_x) < 1L && any(is.finite(x_vals))) rug_x <- x_vals[is.finite(x_vals)]

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

    fmt_cp <- function(cp) sprintf("%.2f", cp)

    if (length(cross_pts_10) > 0) {
      for (cp in cross_pts_10) {
        p <- p +
          geom_vline(xintercept = cp, linetype = "dotted", color = "#333333",
                     linewidth = 0.7, alpha = 0.8) +
          annotate("point", x = cp, y = HR_REF_10, color = "#333333",
                   size = 3, shape = 18) +
          annotate("label", x = cp, y = 1.38, label = fmt_cp(cp), color = "#333333",
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
          annotate("label", x = cp, y = 0.60, label = fmt_cp(cp), color = "#0072B2",
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
          annotate("label", x = cp, y = 0.50, label = fmt_cp(cp), color = "#DC0000",
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
          annotate("label", x = cp, y = 0.38, label = fmt_cp(cp), color = "#7E6148",
                    size = 3.3, fontface = "bold", fill = "white",
                    label.size = 0.3, label.padding = unit(0.15, "lines"))
      }
    }

    if (is_dynamic) {
      p <- p + geom_vline(xintercept = 0, linetype = "dotted", color = "#888888", linewidth = 0.6)
    }

    p <- p +
      geom_rug(data = data.frame(xr = rug_x),
               aes(x = xr), sides = "b", alpha = 0.42, color = "#333333",
               length = unit(2.5, "mm"), linewidth = 0.4, inherit.aes = FALSE) +
      annotate("text", x = xlim[1] + diff(xlim) * 0.03, y = 0.50,
               label = paste0("Favors\n", trt_label), color = "#2E7D32",
               size = 3.8, hjust = 0, fontface = "italic") +
      annotate("text", x = xlim[1] + diff(xlim) * 0.03, y = 1.65,
               label = paste0("Favors\n", ctrl_label), color = "#E65100",
               size = 3.8, hjust = 0, fontface = "italic") +
      labs(
        title = sprintf("%s %s", title_tag, var_label),
        subtitle = sprintf(
          "rms::rcs(nk=%d) IPTW-cph | Int. p=%s | Nonlin.int. p=%s | n=%d",
          nk,
          ifelse(is.na(int_p), "NA", sprintf("%.3f", int_p)),
          ifelse(is.na(nonlin_p), "NA", sprintf("%.3f", nonlin_p)), n),
        x = var_label,
        y = sprintf("HR (%s vs %s)", trt_label, ctrl_label),
        caption = paste0(caption_core,
          " | Harrell RCS | P5\u2013P95 | 95% bootstrap CI | Gray: HR=1 | Blue: 0.9 | Red: 0.8 | Brown: 0.7")
      ) +
      scale_y_log10(breaks = c(0.3, 0.5, 0.7, 0.9, 1.0, 1.5, 2.0),
                    labels = c("0.3", "0.5", "0.7", "0.9", "1.0", "1.5", "2.0")) +
      coord_cartesian(ylim = c(0.25, 2.5), xlim = xlim) +
      scale_x_continuous(limits = xlim,
                         expand = ggplot2::expansion(mult = c(0.04, 0.04))) +
      theme_bw(base_size = 13)

    list(plot = p, int_p = int_p, nonlin_p = nonlin_p, n = n, nk = nk,
         cross_pts10 = cross_pts_10, cross_pts09 = cross_pts_09,
         cross_pts08 = cross_pts_08, cross_pts07 = cross_pts_07,
         hr_first = hr_curve[1], hr_last = hr_curve[length(hr_curve)],
         xlim = xlim)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# ── 变量配置 ─────────────────────────────────────────────────────────────────
vars_static <- list(
  list(col = "log_afp_pivka_ratio_bl",
       label = "log(AFP/PIVKA) ratio (baseline)", group = "baseline"),
  list(col = "log_afp_pivka_sum_bl",
       label = "log(AFP)+log(PIVKA) burden (baseline)", group = "baseline"),
  list(col = "afp_pivka_pc1_bl",
       label = "AFP-PIVKA PC1 shared burden (baseline)", group = "baseline"),
  list(col = "afp_pivka_pc2_bl",
       label = "AFP-PIVKA PC2 discordance (baseline)", group = "baseline"),

  list(col = "log_afp_pivka_ratio_pre3",
       label = "log(AFP/PIVKA) ratio (pre-HAIC-3)", group = "pre3"),
  list(col = "log_afp_pivka_sum_pre3",
       label = "log(AFP)+log(PIVKA) burden (pre-HAIC-3)", group = "pre3"),
  list(col = "afp_pivka_pc1_pre3",
       label = "AFP-PIVKA PC1 shared burden (pre-HAIC-3)", group = "pre3"),
  list(col = "afp_pivka_pc2_pre3",
       label = "AFP-PIVKA PC2 discordance (pre-HAIC-3)", group = "pre3")
)

vars_dynamic <- list(
  # --- 独立变化率（IPTW 时交叉平衡对方变化率）---
  list(col = "afp_change_pre3",
       label = "AFP change rate (%)"),
  list(col = "pivka_change_pre3",
       label = "PIVKA change rate (%)"),
  # --- 组合指标 ---
  list(col = "afp_pivka_change_diff",
       label = "AFP-PIVKA change difference (%)"),
  list(col = "afp_pivka_change_sum",
       label = "AFP+PIVKA change sum (%)"),
  list(col = "afp_pivka_pc1_dyn",
       label = "AFP-PIVKA PC1 dynamic (shared response)"),
  list(col = "afp_pivka_pc2_dyn",
       label = "AFP-PIVKA PC2 dynamic (discordance)")
)

# ── 运行一个队列×时间尺度 ──────────────────────────────────────────────────
run_one_arm <- function(df_arm, out_dir, surv_time_col, caption_core, title_tag,
                        trt_col, trt_label, ctrl_label) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("ARM:", out_dir, "| Surv(", surv_time_col, ") | n =", nrow(df_arm),
      "| trt_col:", trt_col, "\n")
  cat(strrep("=", 70), "\n")

  ps_vars_full <- c("albi_bl", "alb_bl", "tbil_bl", "inr", "plt", "creatinine",
                    "tumor_max_diameter_cm", "tumor_count_enc",
                    "pvtt_grade", "hvtt_binary", "ivc_ra_binary", "ascites_score_enc",
                    "log_afp_bl", "log_pivka_bl",
                    "metastasis_binary", "lymph_node_binary",
                    "neut_bl", "lymph_bl", "mono_bl", "egv_binary",
                    # 新增：pre-HAIC-3 相对 baseline 的变化率
                    "afp_change_pre3", "pivka_change_pre3", "albi_change_pre3")
  df <- impute_ps_vars(df_arm, ps_vars_full)
  ps_avail <- ps_vars_full[ps_vars_full %in% names(df)]
  ps_avail <- ps_avail[sapply(ps_avail, function(v) sd(df[[v]], na.rm = TRUE) > 0)]

  # --- Static ---
  results_static <- list()
  plots_static   <- list()
  for (cfg in vars_static) {
    if (!(cfg$col %in% names(df))) next
    res <- plot_rcs_composite(df, cfg$col, cfg$label, surv_time_col, caption_core,
                               title_tag, trt_col, trt_label, ctrl_label,
                               ps_avail, is_dynamic = FALSE)
    if (!is.null(res)) {
      results_static[[cfg$col]] <- c(res, list(label = cfg$label, group = cfg$group))
      plots_static[[cfg$col]]   <- res$plot
    }
  }

  bl_keys <- names(plots_static)[sapply(names(plots_static),
    function(k) results_static[[k]]$group == "baseline")]
  p3_keys <- names(plots_static)[sapply(names(plots_static),
    function(k) results_static[[k]]$group == "pre3")]

  if (length(bl_keys) >= 2) {
    comb <- gridExtra::grid.arrange(
      grobs = plots_static[bl_keys], ncol = 2,
      top = grid::textGrob(
        paste(title_tag, sprintf("AFP-PIVKA Composite RCS — Baseline (%s vs %s)", ctrl_label, trt_label)),
        gp = grid::gpar(fontsize = 13, fontface = "bold")))
    h_bl <- ceiling(length(bl_keys) / 2) * 4.2
    ggsave(file.path(out_dir, "afp_pivka_static_baseline_combined.pdf"),
           comb, width = 11.5, height = h_bl, device = "pdf")
    ggsave(file.path(out_dir, "afp_pivka_static_baseline_combined.png"),
           comb, width = 11.5, height = h_bl, dpi = 300)
  }
  if (length(p3_keys) >= 2) {
    comb <- gridExtra::grid.arrange(
      grobs = plots_static[p3_keys], ncol = 2,
      top = grid::textGrob(
        paste(title_tag, sprintf("AFP-PIVKA Composite RCS — Pre-HAIC-3 (%s vs %s)", ctrl_label, trt_label)),
        gp = grid::gpar(fontsize = 13, fontface = "bold")))
    h_p3 <- ceiling(length(p3_keys) / 2) * 4.2
    ggsave(file.path(out_dir, "afp_pivka_static_pre3_combined.pdf"),
           comb, width = 11.5, height = h_p3, device = "pdf")
    ggsave(file.path(out_dir, "afp_pivka_static_pre3_combined.png"),
           comb, width = 11.5, height = h_p3, dpi = 300)
  }

  if (length(results_static) > 0) {
    smry <- do.call(rbind, lapply(names(results_static), function(k) {
      r <- results_static[[k]]
      sf <- function(x) ifelse(is.null(x) || length(x) == 0 || is.na(x[1]), NA_real_, x[1])
      data.frame(variable = k, label = r$label, group = r$group, n = r$n, nk = r$nk,
                 interaction_p = round(sf(r$int_p), 4),
                 nonlinear_interaction_p = round(sf(r$nonlin_p), 4),
                 hr_at_xmin = round(r$hr_first, 4), hr_at_xmax = round(r$hr_last, 4),
                 n_cross10 = length(r$cross_pts10), n_cross09 = length(r$cross_pts09),
                 n_cross08 = length(r$cross_pts08), n_cross07 = length(r$cross_pts07),
                 cross1_10 = ifelse(length(r$cross_pts10) >= 1, round(r$cross_pts10[1], 4), NA_real_),
                 cross1_09 = ifelse(length(r$cross_pts09) >= 1, round(r$cross_pts09[1], 4), NA_real_),
                 cross1_08 = ifelse(length(r$cross_pts08) >= 1, round(r$cross_pts08[1], 4), NA_real_),
                 cross1_07 = ifelse(length(r$cross_pts07) >= 1, round(r$cross_pts07[1], 4), NA_real_),
                 stringsAsFactors = FALSE)
    }))
    write.csv(smry, file.path(out_dir, "afp_pivka_static_anova_summary.csv"), row.names = FALSE)
  }

  # --- Dynamic ---
  results_dyn <- list()
  plots_dyn   <- list()
  for (cfg in vars_dynamic) {
    if (!(cfg$col %in% names(df))) next
    res <- plot_rcs_composite(df, cfg$col, cfg$label, surv_time_col, caption_core,
                               title_tag, trt_col, trt_label, ctrl_label,
                               ps_avail, is_dynamic = TRUE)
    if (!is.null(res)) {
      results_dyn[[cfg$col]] <- c(res, list(label = cfg$label))
      plots_dyn[[cfg$col]]   <- res$plot
    }
  }
  if (length(plots_dyn) >= 2) {
    dyn_order <- vapply(vars_dynamic, function(z) z$col, character(1))
    dyn_order <- dyn_order[dyn_order %in% names(plots_dyn)]
    comb <- gridExtra::grid.arrange(
      grobs = plots_dyn[dyn_order], ncol = 2,
      top = grid::textGrob(
        paste(title_tag, sprintf("AFP-PIVKA Composite RCS — Dynamic (%s vs %s)", ctrl_label, trt_label)),
        gp = grid::gpar(fontsize = 13, fontface = "bold")))
    h <- ceiling(length(plots_dyn) / 2) * 4.2
    ggsave(file.path(out_dir, "afp_pivka_dynamic_combined.pdf"),
           comb, width = 11.5, height = h, device = "pdf")
    ggsave(file.path(out_dir, "afp_pivka_dynamic_combined.png"),
           comb, width = 11.5, height = h, dpi = 300)
  }
  if (length(results_dyn) > 0) {
    smry_d <- do.call(rbind, lapply(names(results_dyn), function(k) {
      r <- results_dyn[[k]]
      sf <- function(x) ifelse(is.null(x) || length(x) == 0 || is.na(x[1]), NA_real_, x[1])
      data.frame(variable = k, label = r$label, group = "dynamic", n = r$n, nk = r$nk,
                 interaction_p = round(sf(r$int_p), 4),
                 nonlinear_interaction_p = round(sf(r$nonlin_p), 4),
                 hr_at_xmin = round(r$hr_first, 4), hr_at_xmax = round(r$hr_last, 4),
                 n_cross10 = length(r$cross_pts10), n_cross09 = length(r$cross_pts09),
                 n_cross08 = length(r$cross_pts08), n_cross07 = length(r$cross_pts07),
                 cross1_10 = ifelse(length(r$cross_pts10) >= 1, round(r$cross_pts10[1], 4), NA_real_),
                 cross1_09 = ifelse(length(r$cross_pts09) >= 1, round(r$cross_pts09[1], 4), NA_real_),
                 cross1_08 = ifelse(length(r$cross_pts08) >= 1, round(r$cross_pts08[1], 4), NA_real_),
                 cross1_07 = ifelse(length(r$cross_pts07) >= 1, round(r$cross_pts07[1], 4), NA_real_),
                 stringsAsFactors = FALSE)
    }))
    write.csv(smry_d, file.path(out_dir, "afp_pivka_dynamic_anova_summary.csv"), row.names = FALSE)
  }
  cat("Saved outputs under:", out_dir, "\n")
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

  df_lm    <- df0[!is.na(df0$os_lm) & df0$os_lm > 0, ]
  df_total <- df0[!is.na(df0$os_months) & df0$os_months > 0, ]

  base_out <- file.path(SCRIPT_DIR, "output", key)

  run_one_arm(df_lm,    file.path(base_out, "landmark"), "os_lm",
              "42-day landmark residual OS", "[Landmark]",
              reg$trt_col, reg$trt_label, reg$ctrl_label)
  run_one_arm(df_total, file.path(base_out, "total_os"), "os_months",
              "Total OS from baseline", "[Total OS]",
              reg$trt_col, reg$trt_label, reg$ctrl_label)
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

cat("\n=== AFP-PIVKA Composite RCS analysis completed ===\n")
