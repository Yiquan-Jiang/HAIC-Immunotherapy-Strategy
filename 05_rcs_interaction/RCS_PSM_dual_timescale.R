#!/usr/bin/env Rscript
# =============================================================================
# RCS_PSM_dual_timescale.R
#
# 基于 PSM 匹配队列的 RCS 非线性交互分析（无 IPTW/OW 权重）
#
# 队列: 6个 HAIC_alone vs 其余组的 PSM 匹配对子
# 模型: Surv(...) ~ trt * rcs(rcsx, nk)，无权重
#
# 依赖: survival, rms, Hmisc, ggplot2, dplyr, gridExtra, grid
#
# 数据: cohort_0X_*.csv（先运行 build_cohort_psm.py）
# 输出: <project_root>/output/step1_rcs_interaction/psm/pair_XX/{landmark,total_os}/
# 环境变量: RMS_RCS_N_BOOT, RMS_RCS_NK（可选）
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(rms)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(grid)
})

args_all <- commandArgs(trailingOnly = FALSE)
fa <- args_all[grepl("^--file=", args_all)]
SCRIPT_DIR <- if (length(fa)) {
  dirname(normalizePath(sub("^--file=", "", fa[1]), winslash = "/"))
} else {
  getwd()
}

PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), winslash = "/")
BASE_OUT <- file.path(PROJECT_ROOT, "output", "step1_rcs_interaction", "psm")
dir.create(BASE_OUT, showWarnings = FALSE, recursive = TRUE)
BASE_OUT <- normalizePath(BASE_OUT, winslash = "/", mustWork = TRUE)

nk_env <- Sys.getenv("RMS_RCS_NK", "3")
RMS_RCS_NK <- suppressWarnings(as.integer(nk_env))
if (length(RMS_RCS_NK) != 1L || is.na(RMS_RCS_NK) || RMS_RCS_NK < 3L) RMS_RCS_NK <- 3L

nb_env <- Sys.getenv("RMS_RCS_N_BOOT", "")
N_BOOT <- if (nzchar(nb_env)) as.integer(nb_env) else 200L
if (length(N_BOOT) != 1L || is.na(N_BOOT) || N_BOOT < 2L) N_BOOT <- 200L

MIN_N <- 40L

cat("SCRIPT_DIR:", SCRIPT_DIR, "\nBASE_OUT:", BASE_OUT,
    "\nRMS_RCS_NK:", RMS_RCS_NK, "| N_BOOT:", N_BOOT, "\n\n")

# ── 6个对子配置 ───────────────────────────────────────────────────────────────
PAIRS <- list(
  list(id = "01", file = "cohort_01_HAIC_alone_vs_HAIC_I_conc.csv",
       group1 = "HAIC_alone", group2 = "HAIC+I_concurrent",
       folder = "pair_01_HAIC_alone_vs_HAIC_I_conc"),
  list(id = "02", file = "cohort_02_HAIC_alone_vs_HAIC_then_I.csv",
       group1 = "HAIC_alone", group2 = "HAIC_then_I",
       folder = "pair_02_HAIC_alone_vs_HAIC_then_I"),
  list(id = "03", file = "cohort_03_HAIC_alone_vs_HAIC_T_conc.csv",
       group1 = "HAIC_alone", group2 = "HAIC+T_concurrent",
       folder = "pair_03_HAIC_alone_vs_HAIC_T_conc"),
  list(id = "04", file = "cohort_04_HAIC_alone_vs_HAIC_then_T.csv",
       group1 = "HAIC_alone", group2 = "HAIC_then_T",
       folder = "pair_04_HAIC_alone_vs_HAIC_then_T"),
  list(id = "05", file = "cohort_05_HAIC_alone_vs_HAIC_IT_conc.csv",
       group1 = "HAIC_alone", group2 = "HAIC+I+T_concurrent",
       folder = "pair_05_HAIC_alone_vs_HAIC_IT_conc"),
  list(id = "06", file = "cohort_06_HAIC_alone_vs_HAIC_then_IT.csv",
       group1 = "HAIC_alone", group2 = "HAIC_then_I+T",
       folder = "pair_06_HAIC_alone_vs_HAIC_then_IT")
)

# ── anova P 值解析 ────────────────────────────────────────────────────────────
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
  a  <- anova(fit)
  rn <- rownames(a)
  pcol <- if ("P" %in% colnames(a)) "P" else ncol(a)
  int_p    <- NA_real_
  nonlin_p <- NA_real_
  idx_int <- grep("^trt \\* [^\\(]+  \\(Factor\\+Higher Order Factors\\)$", rn)
  if (length(idx_int)) int_p <- parse_anova_p(a[idx_int[1], pcol])
  idx_nl <- grep("Nonlinear Interaction : f(A,B) vs. AB", rn, fixed = TRUE)
  if (length(idx_nl)) nonlin_p <- parse_anova_p(a[idx_nl[1], pcol])
  list(int_p = int_p, nonlin_p = nonlin_p)
}

# ── 模型拟合辅助 ──────────────────────────────────────────────────────────────
fit_boot_cph_rms <- function(df_b, surv_time_col, nk) {
  dd_b <- suppressWarnings(datadist(df_b[, c("trt", "rcsx"), drop = FALSE]))
  options(datadist = dd_b)
  fml <- as.formula(paste0(
    "Surv(", surv_time_col, ", death_status) ~ trt * rcs(rcsx, ", nk, ")"
  ))
  cph(fml, data = df_b, x = TRUE, y = TRUE, robust = FALSE)
}

predict_hr_curve <- function(fit, x_grid) {
  nd1 <- data.frame(trt = 1, rcsx = x_grid)
  nd0 <- data.frame(trt = 0, rcsx = x_grid)
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

# ── 静态变量绘图函数 ──────────────────────────────────────────────────────────
plot_rms_rcs_static <- function(df_sub, var_col, var_label, surv_time_col,
                                caption_core, title_tag, log_transform = FALSE,
                                group1, group2) {
  df_sub <- df_sub[!is.na(df_sub[[var_col]]), ]
  n <- nrow(df_sub)
  cat(sprintf("\n--- [%s] %s: n=%d ---\n", title_tag, var_label, n))
  if (n < MIN_N) { cat("  SKIP: n < MIN_N\n"); return(NULL) }

  if (log_transform) {
    df_sub$rcsx   <- log1p(pmax(df_sub[[var_col]], 0))
    x_vals        <- df_sub$rcsx
    x_label_plot  <- sprintf("%s (log scale)", var_label)
  } else {
    df_sub$rcsx   <- df_sub[[var_col]]
    x_vals        <- df_sub$rcsx
    x_label_plot  <- var_label
  }

  nk <- RMS_RCS_NK
  dd <- suppressWarnings(datadist(df_sub[, c("trt", "rcsx"), drop = FALSE]))
  options(datadist = dd)

  tryCatch({
    fml <- as.formula(paste0(
      "Surv(", surv_time_col, ", death_status) ~ trt * rcs(rcsx, ", nk, ")"
    ))
    fit <- cph(fml, data = df_sub, x = TRUE, y = TRUE, robust = FALSE)
    ap  <- extract_rms_anova_p(fit)

    xlim   <- quantile(df_sub$rcsx, c(0.05, 0.95))
    x_grid <- seq(xlim[1], xlim[2], length.out = 200)
    hr_curve <- predict_hr_curve(fit, x_grid)

    set.seed(42)
    hr_boot <- matrix(NA_real_, N_BOOT, length(x_grid))
    for (b in seq_len(N_BOOT)) {
      idx  <- sample.int(nrow(df_sub), replace = TRUE)
      df_b <- df_sub[idx, ]
      tryCatch({
        fit_b <- fit_boot_cph_rms(df_b, surv_time_col, nk)
        hr_boot[b, ] <- predict_hr_curve(fit_b, x_grid)
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

    HR_REFS    <- c(1.0, 0.85, 0.7)
    REF_COLORS <- c("#333333", "#0072B2", "#DC0000")
    REF_YTXT   <- c(-0.5, -0.5, 1.5)
    REF_LT     <- c("dashed", "dotdash", "dotdash")
    REF_LW     <- c(1.0, 0.8, 0.8)

    plot_df <- data.frame(x = x_grid, hr = hr_curve, hr_lo = hr_lo, hr_hi = hr_hi)
    cross_list <- lapply(HR_REFS, function(r) find_crossings(plot_df, r))

    fmt_x_cross <- function(cp) {
      if (log_transform) sprintf("%.3g", expm1(cp)) else sprintf("%.3g", cp)
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
      geom_line(aes(y = hr), color = "#3C5488", linewidth = 1.8)

    for (ri in seq_along(HR_REFS)) {
      hr_label <- if (HR_REFS[ri] == 0.85) "HR=0.85" else sprintf("HR=%.1f", HR_REFS[ri])
      p <- p + geom_hline(yintercept = HR_REFS[ri], linetype = REF_LT[ri],
                          color = REF_COLORS[ri], linewidth = REF_LW[ri], alpha = 0.9) +
        annotate("text", x = xlim[2] - diff(xlim) * 0.02, y = HR_REFS[ri],
                 label = hr_label, color = REF_COLORS[ri],
                 size = 3.5, hjust = 1, vjust = REF_YTXT[ri], fontface = "bold")
    }

    cross_y_labels <- c(1.38, 0.60, 0.45)
    for (ri in seq_along(HR_REFS)) {
      cps <- cross_list[[ri]]
      if (length(cps) > 0) {
        for (cp in cps) {
          lab <- fmt_x_cross(cp)
          p <- p +
            geom_vline(xintercept = cp, linetype = "dotted",
                       color = REF_COLORS[ri], linewidth = 0.7, alpha = 0.8) +
            annotate("point", x = cp, y = HR_REFS[ri],
                     color = REF_COLORS[ri], size = 3, shape = 18) +
            annotate("label", x = cp, y = cross_y_labels[ri], label = lab,
                     color = REF_COLORS[ri], size = 3.5, fontface = "bold",
                     fill = "white", label.size = 0.3,
                     label.padding = unit(0.15, "lines"))
        }
      }
    }

    short_title <- sprintf("%s vs %s  %s", group2, group1, var_label)
    p <- p +
      geom_rug(data = data.frame(xr = rug_x),
               aes(x = xr), sides = "b", alpha = 0.42, color = "#333333",
               length = unit(2.5, "mm"), linewidth = 0.4, inherit.aes = FALSE) +
      annotate("text", x = xlim[1] + diff(xlim) * 0.03, y = 0.38,
               label = sprintf("Favors\n%s", group2),
               color = "#2E7D32", size = 3.8, hjust = 0, fontface = "italic") +
      annotate("text", x = xlim[1] + diff(xlim) * 0.03, y = 1.90,
               label = sprintf("Favors\n%s", group1),
               color = "#E65100", size = 3.8, hjust = 0, fontface = "italic") +
      labs(
        title    = short_title,
        subtitle = sprintf(
          "Int. p=%s | Nonlin.int. p=%s | n=%d",
          ifelse(is.na(ap$int_p),    "NA", sprintf("%.3f", ap$int_p)),
          ifelse(is.na(ap$nonlin_p), "NA", sprintf("%.3f", ap$nonlin_p)), n),
        x       = x_label_plot,
        y       = sprintf("HR (%s vs %s)", group2, group1)
      ) +
      scale_y_log10(breaks = c(0.3, 0.5, 0.7, 0.85, 1.0, 1.5, 2.0),
                    labels = c("0.3", "0.5", "0.7", "0.85", "1.0", "1.5", "2.0")) +
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
      if (length(raw_ticks) > 8)
        raw_ticks <- raw_ticks[seq(1, length(raw_ticks), length.out = min(8, length(raw_ticks)))]
      log_ticks <- log1p(raw_ticks)
      fmt_label <- function(v) ifelse(v >= 1000, sprintf("%gK", v / 1000),
                                      ifelse(v >= 1, sprintf("%g", v), sprintf("%.1f", v)))
      p <- p + scale_x_continuous(breaks = log_ticks, labels = fmt_label(raw_ticks),
                                   name = x_label_plot, limits = xlim,
                                   expand = ggplot2::expansion(mult = c(0.02, 0.02)))
    } else {
      p <- p + scale_x_continuous(limits = xlim,
                                   expand = ggplot2::expansion(mult = c(0.04, 0.04)))
    }

    cross_pts_raw <- lapply(seq_along(HR_REFS), function(ri) {
      cps <- cross_list[[ri]]
      if (log_transform) expm1(cps) else cps
    })

    list(plot = p, int_p = ap$int_p, nonlin_p = ap$nonlin_p, n = n, nk = nk,
         cross_pts10_raw = cross_pts_raw[[1]],
         cross_pts085_raw = cross_pts_raw[[2]],
         cross_pts07_raw = cross_pts_raw[[3]],
         hr_first = hr_curve[1], hr_last = hr_curve[length(hr_curve)],
         xlim_raw = if (log_transform) expm1(xlim) else xlim)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# ── 动态变量绘图函数 ──────────────────────────────────────────────────────────
plot_rms_rcs_dynamic <- function(df_sub, var_col, var_label, surv_time_col,
                                  caption_core, title_tag, group1, group2) {
  df_sub <- df_sub[!is.na(df_sub[[var_col]]), ]
  n <- nrow(df_sub)
  cat(sprintf("\n--- [%s] %s: n=%d ---\n", title_tag, var_label, n))
  if (n < MIN_N) { cat("  SKIP: n < MIN_N\n"); return(NULL) }

  df_sub$rcsx  <- df_sub[[var_col]]
  x_vals       <- df_sub$rcsx
  x_label_plot <- sprintf("%s (%%)", var_label)
  nk   <- RMS_RCS_NK
  dd   <- suppressWarnings(datadist(df_sub[, c("trt", "rcsx"), drop = FALSE]))
  options(datadist = dd)
  xlim <- quantile(x_vals, c(0.05, 0.95))

  tryCatch({
    fml <- as.formula(paste0(
      "Surv(", surv_time_col, ", death_status) ~ trt * rcs(rcsx, ", nk, ")"
    ))
    fit <- cph(fml, data = df_sub, x = TRUE, y = TRUE, robust = FALSE)
    ap  <- extract_rms_anova_p(fit)

    x_grid   <- seq(xlim[1], xlim[2], length.out = 200)
    hr_curve <- predict_hr_curve(fit, x_grid)

    set.seed(42)
    hr_boot <- matrix(NA_real_, N_BOOT, length(x_grid))
    for (b in seq_len(N_BOOT)) {
      idx  <- sample.int(nrow(df_sub), replace = TRUE)
      df_b <- df_sub[idx, ]
      tryCatch({
        fit_b <- fit_boot_cph_rms(df_b, surv_time_col, nk)
        hr_boot[b, ] <- predict_hr_curve(fit_b, x_grid)
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

    HR_REFS    <- c(1.0, 0.85, 0.7)
    REF_COLORS <- c("#333333", "#0072B2", "#DC0000")
    REF_YTXT   <- c(-0.5, -0.5, 1.5)
    REF_LT     <- c("dashed", "dotdash", "dotdash")
    REF_LW     <- c(1.0, 0.8, 0.8)

    plot_df    <- data.frame(x = x_grid, hr = hr_curve, hr_lo = hr_lo, hr_hi = hr_hi)
    cross_list <- lapply(HR_REFS, function(r) find_crossings(plot_df, r))

    cp_disp <- function(cp) {
      ifelse(abs(cp) >= 10, sprintf("%.1f%%", cp), sprintf("%.1f%%", cp))
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
      geom_line(aes(y = hr), color = "#3C5488", linewidth = 1.8)

    for (ri in seq_along(HR_REFS)) {
      hr_label <- if (HR_REFS[ri] == 0.85) "HR=0.85" else sprintf("HR=%.1f", HR_REFS[ri])
      p <- p + geom_hline(yintercept = HR_REFS[ri], linetype = REF_LT[ri],
                          color = REF_COLORS[ri], linewidth = REF_LW[ri], alpha = 0.9) +
        annotate("text", x = xlim[2] - diff(xlim) * 0.02, y = HR_REFS[ri],
                 label = hr_label, color = REF_COLORS[ri],
                 size = 3.5, hjust = 1, vjust = REF_YTXT[ri], fontface = "bold")
    }

    cross_y_labels <- c(1.38, 0.60, 0.45)
    for (ri in seq_along(HR_REFS)) {
      cps <- cross_list[[ri]]
      if (length(cps) > 0) {
        for (cp in cps) {
          p <- p +
            geom_vline(xintercept = cp, linetype = "dotted",
                       color = REF_COLORS[ri], linewidth = 0.7, alpha = 0.8) +
            annotate("point", x = cp, y = HR_REFS[ri],
                     color = REF_COLORS[ri], size = 3, shape = 18) +
            annotate("label", x = cp, y = cross_y_labels[ri], label = cp_disp(cp),
                     color = REF_COLORS[ri], size = 3.5, fontface = "bold",
                     fill = "white", label.size = 0.3,
                     label.padding = unit(0.15, "lines"))
        }
      }
    }

    short_title <- sprintf("%s vs %s  %s", group2, group1, var_label)
    p <- p +
      geom_vline(xintercept = 0, linetype = "dotted", color = "#888888", linewidth = 0.6) +
      geom_rug(data = data.frame(xr = rug_x),
               aes(x = xr), sides = "b", alpha = 0.42, color = "#333333",
               length = unit(2.5, "mm"), linewidth = 0.4, inherit.aes = FALSE) +
      annotate("text", x = xlim[1] + diff(xlim) * 0.03, y = 0.38,
               label = sprintf("Favors\n%s", group2),
               color = "#2E7D32", size = 3.8, hjust = 0, fontface = "italic") +
      annotate("text", x = xlim[1] + diff(xlim) * 0.03, y = 1.90,
               label = sprintf("Favors\n%s", group1),
               color = "#E65100", size = 3.8, hjust = 0, fontface = "italic") +
      annotate("text", x = 0, y = 0.28, label = "No change",
               color = "#888888", size = 3.5, hjust = 0.5, fontface = "italic") +
      labs(
        title    = short_title,
        subtitle = sprintf(
          "Int. p=%s | Nonlin.int. p=%s | n=%d",
          ifelse(is.na(ap$int_p),    "NA", sprintf("%.3f", ap$int_p)),
          ifelse(is.na(ap$nonlin_p), "NA", sprintf("%.3f", ap$nonlin_p)), n),
        x       = x_label_plot,
        y       = sprintf("HR (%s vs %s)", group2, group1)
      ) +
      scale_y_log10(breaks = c(0.3, 0.5, 0.7, 0.85, 1.0, 1.5, 2.0),
                    labels = c("0.3", "0.5", "0.7", "0.85", "1.0", "1.5", "2.0")) +
      coord_cartesian(ylim = c(0.25, 2.5), xlim = xlim) +
      theme_bw(base_size = 13) +
      scale_x_continuous(limits = xlim,
                          expand = ggplot2::expansion(mult = c(0.04, 0.04)))

    list(plot = p, int_p = ap$int_p, nonlin_p = ap$nonlin_p, n = n, nk = nk,
         cross_pts10_raw = cross_list[[1]],
         cross_pts085_raw = cross_list[[2]],
         cross_pts07_raw = cross_list[[3]],
         hr_first = hr_curve[1], hr_last = hr_curve[length(hr_curve)], xlim_raw = xlim)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# ── 变量配置 ──────────────────────────────────────────────────────────────────
vars_config_static <- list(
  # ── 肿瘤标志物 ──
  list(col = "afp",       label = "AFP (baseline, ng/mL)",              log_transform = TRUE,  group = "baseline"),
  list(col = "pivka",     label = "PIVKA-II (baseline, mAU/mL)",        log_transform = TRUE,  group = "baseline"),
  # ── 肿瘤特征 ──
  list(col = "tumor_max_diameter_cm", label = "Tumor diameter (baseline, cm)", log_transform = FALSE, group = "baseline"),
  # ── 肝功能 ──
  list(col = "albi_bl",   label = "ALBI score (baseline)",              log_transform = FALSE, group = "baseline"),
  list(col = "alb_bl",    label = "Albumin (baseline, g/L)",            log_transform = FALSE, group = "baseline"),
  list(col = "tbil_bl",   label = "Total bilirubin (baseline, \u03bcmol/L)", log_transform = TRUE,  group = "baseline"),
  list(col = "alt_bl",    label = "ALT (baseline, U/L)",                log_transform = TRUE,  group = "baseline"),
  list(col = "ast_bl",    label = "AST (baseline, U/L)",                log_transform = TRUE,  group = "baseline"),
  # ── 炎症指标 ──
  list(col = "nlr_bl",    label = "NLR (baseline)",                     log_transform = FALSE, group = "baseline"),
  list(col = "plr_bl",    label = "PLR (baseline)",                     log_transform = FALSE, group = "baseline"),
  list(col = "sii_bl",    label = "SII (PLT\u00d7NLR, baseline)",       log_transform = TRUE,  group = "baseline"),
  list(col = "piv_bl",    label = "PIV (Mono\u00d7PLT\u00d7NLR, baseline)", log_transform = TRUE,  group = "baseline"),
  # ── 血细胞计数 ──
  list(col = "neut_bl",   label = "Neutrophil (baseline, 10\u2079/L)",  log_transform = FALSE, group = "baseline"),
  list(col = "lymph_bl",  label = "Lymphocyte (baseline, 10\u2079/L)",  log_transform = FALSE, group = "baseline"),
  list(col = "mono_bl",   label = "Monocyte (baseline, 10\u2079/L)",    log_transform = FALSE, group = "baseline"),
  list(col = "plt",       label = "Platelet (baseline, 10\u2079/L)",    log_transform = TRUE,  group = "baseline"),
  # ── Pre-HAIC-3 肿瘤标志物 ──
  list(col = "afp_pre3",  label = "AFP (pre-HAIC-3, ng/mL)",            log_transform = TRUE,  group = "pre3_static"),
  list(col = "pivka_pre3",label = "PIVKA-II (pre-HAIC-3, mAU/mL)",     log_transform = TRUE,  group = "pre3_static"),
  # ── Pre-HAIC-3 肝功能 ──
  list(col = "albi_pre3", label = "ALBI score (pre-HAIC-3)",            log_transform = FALSE, group = "pre3_static"),
  list(col = "alb_pre3",  label = "Albumin (pre-HAIC-3, g/L)",          log_transform = FALSE, group = "pre3_static"),
  list(col = "tbil_pre3", label = "Total bilirubin (pre-HAIC-3, \u03bcmol/L)", log_transform = TRUE, group = "pre3_static"),
  list(col = "alt_pre3",  label = "ALT (pre-HAIC-3, U/L)",              log_transform = TRUE,  group = "pre3_static"),
  list(col = "ast_pre3",  label = "AST (pre-HAIC-3, U/L)",              log_transform = TRUE,  group = "pre3_static"),
  # ── Pre-HAIC-3 炎症指标 ──
  list(col = "nlr_pre3",  label = "NLR (pre-HAIC-3)",                   log_transform = FALSE, group = "pre3_static"),
  list(col = "plr_pre3",  label = "PLR (pre-HAIC-3)",                   log_transform = FALSE, group = "pre3_static"),
  list(col = "sii_pre3",  label = "SII (PLT\u00d7NLR, pre-HAIC-3)",    log_transform = TRUE,  group = "pre3_static"),
  list(col = "piv_pre3",  label = "PIV (Mono\u00d7PLT\u00d7NLR, pre-HAIC-3)", log_transform = TRUE, group = "pre3_static"),
  # ── Pre-HAIC-3 血细胞计数 ──
  list(col = "neut_pre3", label = "Neutrophil (pre-HAIC-3, 10\u2079/L)",log_transform = FALSE, group = "pre3_static"),
  list(col = "lymph_pre3",label = "Lymphocyte (pre-HAIC-3, 10\u2079/L)",log_transform = FALSE, group = "pre3_static"),
  list(col = "mono_pre3", label = "Monocyte (pre-HAIC-3, 10\u2079/L)",  log_transform = FALSE, group = "pre3_static"),
  list(col = "plt_pre3",  label = "Platelet (pre-HAIC-3, 10\u2079/L)",  log_transform = TRUE,  group = "pre3_static"),

  # ── Pre-IT 肿瘤标志物 ──
  list(col = "afp_pre_it",  label = "AFP (pre-IT, ng/mL)",              log_transform = TRUE,  group = "pre_it_static"),
  list(col = "pivka_pre_it",label = "PIVKA-II (pre-IT, mAU/mL)",        log_transform = TRUE,  group = "pre_it_static"),
  # ── Pre-IT 肝功能 ──
  list(col = "albi_pre_it", label = "ALBI score (pre-IT)",              log_transform = FALSE, group = "pre_it_static"),
  list(col = "alb_pre_it",  label = "Albumin (pre-IT, g/L)",            log_transform = FALSE, group = "pre_it_static"),
  list(col = "tbil_pre_it", label = "Total bilirubin (pre-IT, \u03bcmol/L)", log_transform = TRUE,  group = "pre_it_static"),
  list(col = "alt_pre_it",  label = "ALT (pre-IT, U/L)",                log_transform = TRUE,  group = "pre_it_static"),
  list(col = "ast_pre_it",  label = "AST (pre-IT, U/L)",                log_transform = TRUE,  group = "pre_it_static"),
  # ── Pre-IT 炎症指标 ──
  list(col = "nlr_pre_it",  label = "NLR (pre-IT)",                     log_transform = FALSE, group = "pre_it_static"),
  list(col = "plr_pre_it",  label = "PLR (pre-IT)",                     log_transform = FALSE, group = "pre_it_static"),
  list(col = "sii_pre_it",  label = "SII (PLT\u00d7NLR, pre-IT)",      log_transform = TRUE,  group = "pre_it_static"),
  list(col = "piv_pre_it",  label = "PIV (Mono\u00d7PLT\u00d7NLR, pre-IT)", log_transform = TRUE,  group = "pre_it_static"),
  # ── Pre-IT 血细胞计数 ──
  list(col = "neut_pre_it",  label = "Neutrophil (pre-IT, 10\u2079/L)", log_transform = FALSE, group = "pre_it_static"),
  list(col = "lymph_pre_it", label = "Lymphocyte (pre-IT, 10\u2079/L)", log_transform = FALSE, group = "pre_it_static"),
  list(col = "mono_pre_it",  label = "Monocyte (pre-IT, 10\u2079/L)",   log_transform = FALSE, group = "pre_it_static"),
  list(col = "plt_pre_it",   label = "Platelet (pre-IT, 10\u2079/L)",   log_transform = TRUE,  group = "pre_it_static")
)

vars_config_dynamic <- list(
  # ── 肿瘤标志物变化率 ──
  list(col = "afp_change_pre3",   label = "AFP change rate"),
  list(col = "pivka_change_pre3", label = "PIVKA-II change rate"),
  # ── 肝功能变化率 ──
  list(col = "albi_change_pre3",  label = "ALBI change rate"),
  list(col = "alb_change_pre3",   label = "Albumin change rate"),
  list(col = "tbil_change_pre3",  label = "TBIL change rate"),
  list(col = "alt_change_pre3",   label = "ALT change rate"),
  list(col = "ast_change_pre3",   label = "AST change rate"),
  # ── 炎症指标变化率 ──
  list(col = "nlr_change_pre3",   label = "NLR change rate"),
  list(col = "plr_change_pre3",   label = "PLR change rate"),
  list(col = "sii_change_pre3",   label = "SII change rate"),
  list(col = "piv_change_pre3",   label = "PIV change rate"),
  # ── 血细胞计数变化率 ──
  list(col = "neut_change_pre3",  label = "Neutrophil change rate"),
  list(col = "lymph_change_pre3", label = "Lymphocyte change rate"),
  list(col = "mono_change_pre3",  label = "Monocyte change rate"),
  list(col = "plt_change_pre3",   label = "Platelet change rate"),
  # ── 肿瘤标志物变化率 (pre-IT) ──
  list(col = "afp_change_pre_it",  label = "AFP change rate (pre-IT)"),
  list(col = "pivka_change_pre_it",label = "PIVKA change rate (pre-IT)"),
  # ── 肝功能变化率 (pre-IT) ──
  list(col = "albi_change_pre_it", label = "ALBI change rate (pre-IT)"),
  list(col = "alb_change_pre_it",  label = "Albumin change rate (pre-IT)"),
  list(col = "tbil_change_pre_it", label = "TBIL change rate (pre-IT)"),
  list(col = "alt_change_pre_it",  label = "ALT change rate (pre-IT)"),
  list(col = "ast_change_pre_it",  label = "AST change rate (pre-IT)"),
  # ── 炎症指标变化率 (pre-IT) ──
  list(col = "nlr_change_pre_it",  label = "NLR change rate (pre-IT)"),
  list(col = "plr_change_pre_it",  label = "PLR change rate (pre-IT)"),
  list(col = "sii_change_pre_it",  label = "SII change rate (pre-IT)"),
  list(col = "piv_change_pre_it",  label = "PIV change rate (pre-IT)"),
  # ── 血细胞计数变化率 (pre-IT) ──
  list(col = "neut_change_pre_it", label = "Neutrophil change rate (pre-IT)"),
  list(col = "lymph_change_pre_it",label = "Lymphocyte change rate (pre-IT)"),
  list(col = "mono_change_pre_it", label = "Monocyte change rate (pre-IT)"),
  list(col = "plt_change_pre_it",  label = "Platelet change rate (pre-IT)")
)

# ── 保存汇总表辅助 ────────────────────────────────────────────────────────────
save_summary_csv <- function(results, out_dir, prefix) {
  if (length(results) == 0) return(invisible(NULL))
  safe <- function(x) ifelse(is.null(x) || length(x) == 0 || is.na(x[1]), NA_real_, x[1])

  anova_tbl <- do.call(rbind, lapply(names(results), function(k) {
    r <- results[[k]]
    data.frame(variable = k, label = r$label,
               group = if (!is.null(r$group)) r$group else "dynamic",
               n = r$n, nk = r$nk,
               interaction_p = round(safe(r$int_p), 4),
               nonlinear_interaction_p = round(safe(r$nonlin_p), 4),
               stringsAsFactors = FALSE)
  }))
  write.csv(anova_tbl, file.path(out_dir, paste0(prefix, "_anova_summary.csv")), row.names = FALSE)

  cross_tbl <- do.call(rbind, lapply(names(results), function(k) {
    r <- results[[k]]
    cp10  <- r$cross_pts10_raw
    cp085 <- r$cross_pts085_raw
    cp07  <- r$cross_pts07_raw
    data.frame(variable = k, label = r$label,
               group = if (!is.null(r$group)) r$group else "dynamic",
               n_crosspoints_hr10 = length(cp10),
               cross1_hr10 = ifelse(length(cp10) >= 1, round(cp10[1], 4), NA_real_),
               cross2_hr10 = ifelse(length(cp10) >= 2, round(cp10[2], 4), NA_real_),
               n_crosspoints_hr085 = length(cp085),
               cross1_hr085 = ifelse(length(cp085) >= 1, round(cp085[1], 4), NA_real_),
               cross2_hr085 = ifelse(length(cp085) >= 2, round(cp085[2], 4), NA_real_),
               n_crosspoints_hr07 = length(cp07),
               cross1_hr07 = ifelse(length(cp07) >= 1, round(cp07[1], 4), NA_real_),
               cross2_hr07 = ifelse(length(cp07) >= 2, round(cp07[2], 4), NA_real_),
               hr_at_xmin = round(r$hr_first, 4), hr_at_xmax = round(r$hr_last, 4),
               xmin_raw = round(r$xlim_raw[1], 4), xmax_raw = round(r$xlim_raw[2], 4),
               stringsAsFactors = FALSE)
  }))
  write.csv(cross_tbl, file.path(out_dir, paste0(prefix, "_crosspoints.csv")), row.names = FALSE)
}

# ── 单臂分析函数 ──────────────────────────────────────────────────────────────
run_one_arm <- function(df_arm, arm_folder, surv_time_col, caption_core,
                        title_tag, group1, group2) {
  OUT_DIR <- file.path(arm_folder)
  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("ARM:", arm_folder, "| Surv(", surv_time_col, ", death_status) | n =", nrow(df_arm), "\n")
  cat("    group1 =", group1, "| group2 =", group2, "\n")
  cat(strrep("=", 70), "\n")

  # ── 静态变量 ──────────────────────────────────────────────────────────────
  results_static <- list()
  plots_static   <- list()
  for (cfg in vars_config_static) {
    if (!(cfg$col %in% names(df_arm))) next
    res <- plot_rms_rcs_static(df_arm, cfg$col, cfg$label, surv_time_col,
                               caption_core, title_tag,
                               log_transform = cfg$log_transform,
                               group1 = group1, group2 = group2)
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
        sprintf("%s vs %s  Harrell RCS — Baseline (PSM)", group2, group1),
        gp = grid::gpar(fontsize = 14, fontface = "bold")))
    h_bl <- ceiling(length(bl_keys) / 2) * 4.2
    ggsave(file.path(OUT_DIR, "rcs_psm_static_baseline_combined.pdf"),
           comb, width = 11.5, height = h_bl, device = "pdf")
    ggsave(file.path(OUT_DIR, "rcs_psm_static_baseline_combined.png"),
           comb, width = 11.5, height = h_bl, dpi = 300)
    cat("  Saved: rcs_psm_static_baseline_combined\n")
  }
  if (length(p3_keys) >= 2) {
    comb <- gridExtra::grid.arrange(
      grobs = plots_static[p3_keys], ncol = 2,
      top = grid::textGrob(
        sprintf("%s vs %s  Harrell RCS — Pre-HAIC-3 (PSM)", group2, group1),
        gp = grid::gpar(fontsize = 14, fontface = "bold")))
    h_p3 <- ceiling(length(p3_keys) / 2) * 4.2
    ggsave(file.path(OUT_DIR, "rcs_psm_static_pre3_combined.pdf"),
           comb, width = 11.5, height = h_p3, device = "pdf")
    ggsave(file.path(OUT_DIR, "rcs_psm_static_pre3_combined.png"),
           comb, width = 11.5, height = h_p3, dpi = 300)
    cat("  Saved: rcs_psm_static_pre3_combined\n")
  }
  pit_keys <- names(plots_static)[sapply(names(plots_static),
    function(k) results_static[[k]]$group == "pre_it_static")]
  if (length(pit_keys) >= 2) {
    comb <- gridExtra::grid.arrange(
      grobs = plots_static[pit_keys], ncol = 2,
      top = grid::textGrob(
        sprintf("%s vs %s  Harrell RCS — Pre-IT (PSM)", group2, group1),
        gp = grid::gpar(fontsize = 14, fontface = "bold")))
    h_pit <- ceiling(length(pit_keys) / 2) * 4.2
    ggsave(file.path(OUT_DIR, "rcs_psm_static_pre_it_combined.pdf"),
           comb, width = 11.5, height = h_pit, device = "pdf")
    ggsave(file.path(OUT_DIR, "rcs_psm_static_pre_it_combined.png"),
           comb, width = 11.5, height = h_pit, dpi = 300)
    cat("  Saved: rcs_psm_static_pre_it_combined\n")
  }
  save_summary_csv(results_static, OUT_DIR, "rcs_psm_static")

  # ── 动态变量 ──────────────────────────────────────────────────────────────
  results_dyn <- list()
  plots_dyn   <- list()
  for (cfg in vars_config_dynamic) {
    if (!(cfg$col %in% names(df_arm))) next
    res <- plot_rms_rcs_dynamic(df_arm, cfg$col, cfg$label, surv_time_col,
                                 caption_core, title_tag,
                                 group1 = group1, group2 = group2)
    if (!is.null(res)) {
      results_dyn[[cfg$col]] <- c(res, list(label = cfg$label))
      plots_dyn[[cfg$col]]   <- res$plot
    }
  }
  # Split dynamic into pre3 and pre_it groups
  dyn_pre3_keys <- names(plots_dyn)[grepl("_pre3$", names(plots_dyn))]
  dyn_pre_it_keys <- names(plots_dyn)[grepl("_pre_it$", names(plots_dyn))]

  if (length(dyn_pre3_keys) >= 2) {
    dyn_order <- vapply(vars_config_dynamic, function(z) z$col, character(1))
    dyn_order <- dyn_order[dyn_order %in% dyn_pre3_keys]
    comb <- gridExtra::grid.arrange(
      grobs = plots_dyn[dyn_order], ncol = 2,
      top = grid::textGrob(
        sprintf("%s vs %s  Harrell RCS — Dynamic change %% pre3 (PSM)", group2, group1),
        gp = grid::gpar(fontsize = 14, fontface = "bold")))
    h <- ceiling(length(dyn_order) / 2) * 4.2
    ggsave(file.path(OUT_DIR, "rcs_psm_dynamic_combined.pdf"),
           comb, width = 11.5, height = h, device = "pdf")
    ggsave(file.path(OUT_DIR, "rcs_psm_dynamic_combined.png"),
           comb, width = 11.5, height = h, dpi = 300)
    cat("  Saved: rcs_psm_dynamic_combined\n")
  }
  if (length(dyn_pre_it_keys) >= 2) {
    dyn_order_it <- vapply(vars_config_dynamic, function(z) z$col, character(1))
    dyn_order_it <- dyn_order_it[dyn_order_it %in% dyn_pre_it_keys]
    comb <- gridExtra::grid.arrange(
      grobs = plots_dyn[dyn_order_it], ncol = 2,
      top = grid::textGrob(
        sprintf("%s vs %s  Harrell RCS — Dynamic change %% pre-IT (PSM)", group2, group1),
        gp = grid::gpar(fontsize = 14, fontface = "bold")))
    h_it <- ceiling(length(dyn_order_it) / 2) * 4.2
    ggsave(file.path(OUT_DIR, "rcs_psm_dynamic_pre_it_combined.pdf"),
           comb, width = 11.5, height = h_it, device = "pdf")
    ggsave(file.path(OUT_DIR, "rcs_psm_dynamic_pre_it_combined.png"),
           comb, width = 11.5, height = h_it, dpi = 300)
    cat("  Saved: rcs_psm_dynamic_pre_it_combined\n")
  }
  save_summary_csv(results_dyn, OUT_DIR, "rcs_psm_dynamic")

  cat("Outputs saved under:", OUT_DIR, "\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# Main：循环6个对子
# ═══════════════════════════════════════════════════════════════════════════════
LANDMARK_MONTHS <- 42 / 30.44

for (pair in PAIRS) {
  csv_path <- file.path(SCRIPT_DIR, pair$file)
  if (!file.exists(csv_path)) {
    cat("\n[SKIP] 文件不存在:", csv_path,
        "\n先运行: python3 build_cohort_psm.py\n")
    next
  }

  df0 <- read.csv(csv_path, stringsAsFactors = FALSE)
  cat("\n\n", strrep("#", 70), "\n", sep = "")
  cat("PAIR", pair$id, ":", pair$group1, "vs", pair$group2,
      "| n =", nrow(df0), "\n")
  cat(strrep("#", 70), "\n")

  if (!"trt" %in% names(df0)) {
    cat("  [SKIP] 缺少 trt 列\n"); next
  }

  # 衍生列（若 build_cohort_psm.py 未生成）
  if (!"sii_bl" %in% names(df0) && all(c("plt", "nlr_bl") %in% names(df0)))
    df0$sii_bl <- df0$plt * df0$nlr_bl
  if (!"piv_bl" %in% names(df0) && all(c("mono_bl", "plt", "nlr_bl") %in% names(df0)))
    df0$piv_bl <- df0$mono_bl * df0$plt * df0$nlr_bl
  if (!"sii_pre3" %in% names(df0) && all(c("plt_pre3", "nlr_pre3") %in% names(df0)))
    df0$sii_pre3 <- df0$plt_pre3 * df0$nlr_pre3
  if (!"piv_pre3" %in% names(df0) && all(c("mono_pre3", "plt_pre3", "nlr_pre3") %in% names(df0)))
    df0$piv_pre3 <- df0$mono_pre3 * df0$plt_pre3 * df0$nlr_pre3

  # ── Pre-IT 衍生列 ──
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

  # ── Pre-IT 变化率（相对于 baseline） ──
  pre_it_change_map <- list(
    list(pre = "afp_pre_it",   bl = "afp",     out = "afp_change_pre_it"),
    list(pre = "pivka_pre_it", bl = "pivka",   out = "pivka_change_pre_it"),
    list(pre = "alb_pre_it",   bl = "alb_bl",  out = "alb_change_pre_it"),
    list(pre = "tbil_pre_it",  bl = "tbil_bl", out = "tbil_change_pre_it"),
    list(pre = "alt_pre_it",   bl = "alt_bl",  out = "alt_change_pre_it"),
    list(pre = "ast_pre_it",   bl = "ast_bl",  out = "ast_change_pre_it"),
    list(pre = "albi_pre_it",  bl = "albi_bl", out = "albi_change_pre_it"),
    list(pre = "nlr_pre_it",   bl = "nlr_bl",  out = "nlr_change_pre_it"),
    list(pre = "plr_pre_it",   bl = "plr_bl",  out = "plr_change_pre_it"),
    list(pre = "sii_pre_it",   bl = "sii_bl",  out = "sii_change_pre_it"),
    list(pre = "piv_pre_it",   bl = "piv_bl",  out = "piv_change_pre_it"),
    list(pre = "neut_pre_it",  bl = "neut_bl", out = "neut_change_pre_it"),
    list(pre = "lymph_pre_it", bl = "lymph_bl",out = "lymph_change_pre_it"),
    list(pre = "mono_pre_it",  bl = "mono_bl", out = "mono_change_pre_it"),
    list(pre = "plt_pre_it",   bl = "plt",     out = "plt_change_pre_it")
  )
  for (m in pre_it_change_map) {
    if (!(m$out %in% names(df0)) && all(c(m$pre, m$bl) %in% names(df0))) {
      bl_vals <- as.numeric(df0[[m$bl]])
      pre_vals <- as.numeric(df0[[m$pre]])
      df0[[m$out]] <- ifelse(!is.na(bl_vals) & !is.na(pre_vals) & bl_vals != 0,
                             (pre_vals - bl_vals) / abs(bl_vals) * 100, NA)
    }
  }

  if (!"os_lm" %in% names(df0))
    df0$os_lm <- df0$os_months - LANDMARK_MONTHS

  df_landmark <- df0[!is.na(df0$os_lm)     & df0$os_lm     > 0, ]
  df_total    <- df0[!is.na(df0$os_months)  & df0$os_months > 0, ]
  cat("  Landmark subset (os_lm > 0):", nrow(df_landmark), "\n")
  cat("  Total OS subset:", nrow(df_total), "\n")

  pair_out <- file.path(BASE_OUT, pair$folder)

  run_one_arm(
    df_landmark,
    arm_folder   = file.path(pair_out, "landmark"),
    surv_time_col = "os_lm",
    caption_core  = sprintf("42-day landmark residual OS | %s vs %s", pair$group2, pair$group1),
    title_tag     = sprintf("%s vs %s", pair$group2, pair$group1),
    group1 = pair$group1, group2 = pair$group2
  )
  run_one_arm(
    df_total,
    arm_folder   = file.path(pair_out, "total_os"),
    surv_time_col = "os_months",
    caption_core  = sprintf("Total OS from baseline | %s vs %s", pair$group2, pair$group1),
    title_tag     = sprintf("%s vs %s", pair$group2, pair$group1),
    group1 = pair$group1, group2 = pair$group2
  )
}

cat("\n=== All pairs completed. Output root:", BASE_OUT, "===\n")
