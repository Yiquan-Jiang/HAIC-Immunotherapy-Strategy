#!/usr/bin/env Rscript
# =============================================================================
# rcs_all_pairs_dual_timescale.R
#
# HAIC_alone 与其余 6 个治疗组的两两配对 RCS 非线性生存分析。
# 对每个配对分别运行：
#   Surv(...) ~ trt_compare * rcs(rcsx, nk)
# 使用 glmnet ridge 倾向评分 + 稳定 IPTW 加权 rms::cph()，bootstrap 95% CI。
#
# 6 个配对（均以 HAIC_alone 为对照，trt_compare=0）:
#   1. HAIC_alone vs HAIC_then_I
#   2. HAIC_alone vs HAIC_then_I+T
#   3. HAIC_alone vs HAIC_then_T
#   4. HAIC_alone vs HAIC+I_concurrent
#   5. HAIC_alone vs HAIC+I+T_concurrent
#   6. HAIC_alone vs HAIC+T_concurrent
#
# 依赖: survival, rms, Hmisc, ggplot2, dplyr, glmnet, gridExtra, grid
#
# 数据: 先运行 build_all_pairs_cohorts.py 生成各配对 CSV
# 输出: <project_root>/output/step1_rcs_interaction/iptw/<pair_label>/{landmark,total_os}/
#
# 环境变量（可选覆盖）:
#   RMS_RCS_NK      — rcs 节点数，默认 3
#   RMS_RCS_N_BOOT  — bootstrap 次数，默认 200
#   RCS_OUT_ROOT    — 输出根目录，默认 <项目根>/output/step1_rcs_interaction/iptw
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

# ── 路径与参数 ────────────────────────────────────────────────────────────────
args_all <- commandArgs(trailingOnly = FALSE)
fa <- args_all[grepl("^--file=", args_all)]
SCRIPT_DIR <- if (length(fa)) {
  dirname(normalizePath(sub("^--file=", "", fa[1]), winslash = "/"))
} else {
  getwd()
}

nk_env <- Sys.getenv("RMS_RCS_NK", "3")
RMS_RCS_NK <- suppressWarnings(as.integer(nk_env))
if (length(RMS_RCS_NK) != 1L || is.na(RMS_RCS_NK) || RMS_RCS_NK < 3L) RMS_RCS_NK <- 3L

nb_env <- Sys.getenv("RMS_RCS_N_BOOT", "")
N_BOOT <- if (nzchar(nb_env)) as.integer(nb_env) else 200L
if (length(N_BOOT) != 1L || is.na(N_BOOT) || N_BOOT < 2L) N_BOOT <- 200L

PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), winslash = "/")
OUT_ROOT <- Sys.getenv(
  "RCS_OUT_ROOT",
  unset = file.path(PROJECT_ROOT, "output", "step1_rcs_interaction", "iptw")
)
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)
OUT_ROOT <- normalizePath(OUT_ROOT, winslash = "/", mustWork = TRUE)

MIN_N <- 50L

cat("SCRIPT_DIR:", SCRIPT_DIR, "\nOUT_ROOT:", OUT_ROOT,
    "\nRMS_RCS_NK:", RMS_RCS_NK, "| N_BOOT:", N_BOOT, "\n\n")

# ── 6 个配对配置 ──────────────────────────────────────────────────────────────
# list(csv, compare_label, pair_label)
#   csv           — build_all_pairs_cohorts.py 生成的 CSV 文件名
#   compare_label — 对照组的显示名称（用于图标题、标注）
#   pair_label    — 输出子目录名
PAIR_CONFIGS <- list(
  list(csv = "cohort_HAIC_alone_vs_HAIC_then_I.csv",
       compare_label = "HAIC_then_I",
       pair_label    = "HAIC_alone_vs_HAIC_then_I"),
  list(csv = "cohort_HAIC_alone_vs_HAIC_then_IT.csv",
       compare_label = "HAIC_then_I+T",
       pair_label    = "HAIC_alone_vs_HAIC_then_IT"),
  list(csv = "cohort_HAIC_alone_vs_HAIC_then_T.csv",
       compare_label = "HAIC_then_T",
       pair_label    = "HAIC_alone_vs_HAIC_then_T"),
  list(csv = "cohort_HAIC_alone_vs_HAIC_I_conc.csv",
       compare_label = "HAIC+I_concurrent",
       pair_label    = "HAIC_alone_vs_HAIC_I_conc"),
  list(csv = "cohort_HAIC_alone_vs_HAIC_IT_conc.csv",
       compare_label = "HAIC+I+T_concurrent",
       pair_label    = "HAIC_alone_vs_HAIC_IT_conc"),
  list(csv = "cohort_HAIC_alone_vs_HAIC_T_conc.csv",
       compare_label = "HAIC+T_concurrent",
       pair_label    = "HAIC_alone_vs_HAIC_T_conc")
)

# ── IPTW 函数 ─────────────────────────────────────────────────────────────────
# 静态变量：使用全量 PS 变量，不排除任何共线项
compute_iptw_static <- function(df_input, var_col, ps_vars_base) {
  ps_vars_use <- ps_vars_base[ps_vars_base %in% names(df_input)]
  ps_vars_use <- ps_vars_use[sapply(ps_vars_use,
    function(v) sd(df_input[[v]], na.rm = TRUE) > 0)]

  df_cc <- df_input[complete.cases(df_input[, ps_vars_use, drop = FALSE]), ]
  if (nrow(df_cc) < 20) return(df_cc[0, , drop = FALSE])
  X_sc   <- scale(as.matrix(df_cc[, ps_vars_use, drop = FALSE]))
  y      <- df_cc$trt_compare
  cv_fit <- cv.glmnet(X_sc, y, family = "binomial", alpha = 0, nfolds = 5)
  ps_prob <- as.numeric(predict(cv_fit, newx = X_sc, s = "lambda.min", type = "response"))
  ps_prob <- pmin(pmax(ps_prob, 0.05), 0.95)
  p_treat <- mean(y)
  df_cc$sw <- ifelse(df_cc$trt_compare == 1,
                     p_treat / ps_prob,
                     (1 - p_treat) / (1 - ps_prob))
  df_cc
}

# 动态变量（变化率）：使用全量 PS 变量，不排除任何共线项
compute_iptw_dynamic <- function(df_input, var_col, ps_vars_base) {
  ps_vars_use <- ps_vars_base[ps_vars_base %in% names(df_input)]
  ps_vars_use <- ps_vars_use[sapply(ps_vars_use,
    function(v) sd(df_input[[v]], na.rm = TRUE) > 0)]
  df_cc <- df_input[complete.cases(df_input[, ps_vars_use, drop = FALSE]), ]
  if (nrow(df_cc) < 20) return(df_cc[0, , drop = FALSE])
  X_sc   <- scale(as.matrix(df_cc[, ps_vars_use, drop = FALSE]))
  y      <- df_cc$trt_compare
  cv_fit <- cv.glmnet(X_sc, y, family = "binomial", alpha = 0, nfolds = 5)
  ps_prob <- as.numeric(predict(cv_fit, newx = X_sc, s = "lambda.min", type = "response"))
  ps_prob <- pmin(pmax(ps_prob, 0.05), 0.95)
  p_treat <- mean(y)
  df_cc$sw <- ifelse(df_cc$trt_compare == 1,
                     p_treat / ps_prob,
                     (1 - p_treat) / (1 - ps_prob))
  df_cc
}

# ── ANOVA P 值解析 ─────────────────────────────────────────────────────────────
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
  a   <- anova(fit)
  rn  <- rownames(a)
  pcol <- if ("P" %in% colnames(a)) "P" else ncol(a)
  int_p    <- NA_real_
  nonlin_p <- NA_real_
  idx_int <- grep(
    "^trt_compare \\* [^\\(]+  \\(Factor\\+Higher Order Factors\\)$", rn
  )
  if (length(idx_int)) int_p <- parse_anova_p(a[idx_int[1], pcol])
  idx_nl <- grep("Nonlinear Interaction : f(A,B) vs. AB", rn, fixed = TRUE)
  if (length(idx_nl)) nonlin_p <- parse_anova_p(a[idx_nl[1], pcol])
  list(int_p = int_p, nonlin_p = nonlin_p)
}

# ── RCS 拟合与 bootstrap HR 曲线 ──────────────────────────────────────────────
fit_boot_cph_rms <- function(df_b, surv_time_col, nk) {
  dd_b <- suppressWarnings(
    datadist(df_b[, c("trt_compare", "rcsx"), drop = FALSE])
  )
  options(datadist = dd_b)
  fml <- as.formula(paste0(
    "Surv(", surv_time_col, ", death_status) ~ trt_compare * rcs(rcsx, ", nk, ")"
  ))
  cph(fml, data = df_b, weights = sw, x = TRUE, y = TRUE, robust = FALSE)
}

predict_hr_curve <- function(fit, x_grid) {
  nd1 <- data.frame(trt_compare = 1, rcsx = x_grid)
  nd0 <- data.frame(trt_compare = 0, rcsx = x_grid)
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
  list(col = "plt_pre3",  label = "Platelet (pre-HAIC-3, 10\u2079/L)",  log_transform = TRUE,  group = "pre3_static")
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
  list(col = "plt_change_pre3",   label = "Platelet change rate")
)

# ── 核心绘图函数（静态变量）──────────────────────────────────────────────────
plot_rcs_static <- function(df_sub, var_col, var_label, surv_time_col,
                             caption_core, title_tag, compare_label,
                             log_transform = FALSE, ps_vars_base = NULL) {
  df_sub <- compute_iptw_static(df_sub, var_col, ps_vars_base)
  df_sub <- df_sub[!is.na(df_sub[[var_col]]), ]
  n <- nrow(df_sub)
  cat(sprintf("\n--- [%s] %s: n=%d ---\n", title_tag, var_label, n))
  if (n < MIN_N) { cat("  SKIP: n < MIN_N\n"); return(NULL) }

  if (log_transform) {
    df_sub$rcsx <- log1p(pmax(df_sub[[var_col]], 0))
    x_vals      <- df_sub$rcsx
    x_label_plot <- sprintf("%s (log scale)", var_label)
  } else {
    df_sub$rcsx  <- df_sub[[var_col]]
    x_vals       <- df_sub$rcsx
    x_label_plot <- var_label
  }

  nk <- RMS_RCS_NK
  dd <- suppressWarnings(
    datadist(df_sub[, c("trt_compare", "rcsx"), drop = FALSE])
  )
  options(datadist = dd)

  tryCatch({
    fml <- as.formula(paste0(
      "Surv(", surv_time_col, ", death_status) ~ trt_compare * rcs(rcsx, ", nk, ")"
    ))
    fit <- cph(fml, data = df_sub, weights = sw, x = TRUE, y = TRUE, robust = FALSE)
    ap  <- extract_rms_anova_p(fit)

    xlim    <- quantile(df_sub$rcsx, c(0.05, 0.95))
    x_grid  <- seq(xlim[1], xlim[2], length.out = 200)
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

    short_title <- sprintf("%s vs HAIC_alone  %s", compare_label, var_label)
    p <- p +
      geom_rug(data = data.frame(xr = rug_x),
               aes(x = xr), sides = "b", alpha = 0.42, color = "#333333",
               length = unit(2.5, "mm"), linewidth = 0.4, inherit.aes = FALSE) +
      annotate("text", x = xlim[1] + diff(xlim) * 0.03, y = 0.38,
               label = sprintf("Favors\n%s", compare_label),
               color = "#2E7D32", size = 3.8, hjust = 0, fontface = "italic") +
      annotate("text", x = xlim[1] + diff(xlim) * 0.03, y = 1.90,
               label = "Favors\nHAIC_alone",
               color = "#E65100", size = 3.8, hjust = 0, fontface = "italic") +
      labs(
        title    = short_title,
        subtitle = sprintf(
          "Int. p=%s | Nonlin.int. p=%s | n=%d",
          ifelse(is.na(ap$int_p),    "NA", sprintf("%.3f", ap$int_p)),
          ifelse(is.na(ap$nonlin_p), "NA", sprintf("%.3f", ap$nonlin_p)), n),
        x       = x_label_plot,
        y       = sprintf("HR (%s vs HAIC_alone)", compare_label)
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
      fmt_label <- function(v) ifelse(v >= 1000, sprintf("%gK", v/1000),
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
         cross_pts10_raw  = cross_pts_raw[[1]],
         cross_pts085_raw = cross_pts_raw[[2]],
         cross_pts07_raw  = cross_pts_raw[[3]],
         hr_first = hr_curve[1], hr_last = hr_curve[length(hr_curve)],
         xlim_raw = if (log_transform) expm1(xlim) else xlim)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  })
}

# ── 核心绘图函数（动态变量）──────────────────────────────────────────────────
plot_rcs_dynamic <- function(df_sub, var_col, var_label, surv_time_col,
                              caption_core, title_tag, compare_label,
                              ps_vars_base = NULL) {
  df_sub <- compute_iptw_dynamic(df_sub, var_col, ps_vars_base)
  df_sub <- df_sub[!is.na(df_sub[[var_col]]), ]
  n <- nrow(df_sub)
  cat(sprintf("\n--- [%s] %s: n=%d ---\n", title_tag, var_label, n))
  if (n < MIN_N) { cat("  SKIP: n < MIN_N\n"); return(NULL) }

  df_sub$rcsx  <- df_sub[[var_col]]
  x_vals       <- df_sub$rcsx
  x_label_plot <- sprintf("%s (%%)", var_label)
  nk  <- RMS_RCS_NK
  dd  <- suppressWarnings(
    datadist(df_sub[, c("trt_compare", "rcsx"), drop = FALSE])
  )
  options(datadist = dd)
  xlim <- quantile(x_vals, c(0.05, 0.95))

  tryCatch({
    fml <- as.formula(paste0(
      "Surv(", surv_time_col, ", death_status) ~ trt_compare * rcs(rcsx, ", nk, ")"
    ))
    fit <- cph(fml, data = df_sub, weights = sw, x = TRUE, y = TRUE, robust = FALSE)
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

    short_title <- sprintf("%s vs HAIC_alone  %s", compare_label, var_label)
    p <- p +
      geom_vline(xintercept = 0, linetype = "dotted", color = "#888888", linewidth = 0.6) +
      geom_rug(data = data.frame(xr = rug_x),
               aes(x = xr), sides = "b", alpha = 0.42, color = "#333333",
               length = unit(2.5, "mm"), linewidth = 0.4, inherit.aes = FALSE) +
      annotate("text", x = xlim[1] + diff(xlim) * 0.03, y = 0.38,
               label = sprintf("Favors\n%s", compare_label),
               color = "#2E7D32", size = 3.8, hjust = 0, fontface = "italic") +
      annotate("text", x = xlim[1] + diff(xlim) * 0.03, y = 1.90,
               label = "Favors\nHAIC_alone",
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
        y       = sprintf("HR (%s vs HAIC_alone)", compare_label)
      ) +
      scale_y_log10(breaks = c(0.3, 0.5, 0.7, 0.85, 1.0, 1.5, 2.0),
                    labels = c("0.3", "0.5", "0.7", "0.85", "1.0", "1.5", "2.0")) +
      coord_cartesian(ylim = c(0.25, 2.5), xlim = xlim) +
      theme_bw(base_size = 13) +
      scale_x_continuous(limits = xlim,
                          expand = ggplot2::expansion(mult = c(0.04, 0.04)))

    list(plot = p, int_p = ap$int_p, nonlin_p = ap$nonlin_p, n = n, nk = nk,
         cross_pts10_raw  = cross_list[[1]],
         cross_pts085_raw = cross_list[[2]],
         cross_pts07_raw  = cross_list[[3]],
         hr_first = hr_curve[1], hr_last = hr_curve[length(hr_curve)],
         xlim_raw = xlim)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    NULL
  })
}

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

# ── 单时间轴分析（landmark 或 total_os）──────────────────────────────────────
run_one_timescale <- function(df_arm, out_dir, surv_time_col,
                               caption_core, title_tag, compare_label) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("ARM:", basename(out_dir), "| Surv(", surv_time_col, ", death_status) | n =",
      nrow(df_arm), "\n")
  cat(strrep("=", 70), "\n")

  # 静态变量 IPTW：仅使用基线变量
  ps_vars_static_full <- c(
    "albi_bl", "alb_bl", "tbil_bl", "alt_bl", "ast_bl",
    "inr", "plt", "creatinine",
    "tumor_max_diameter_cm", "tumor_count_enc",
    "pvtt_grade", "hvtt_binary", "ivc_ra_binary", "ascites_score_enc",
    "log_afp_bl", "log_pivka_bl",
    "metastasis_binary", "lymph_node_binary",
    "neut_bl", "lymph_bl", "mono_bl", "egv_binary"
  )
  # 动态变量 IPTW：基线变量 + pre-HAIC-3 变化率（ALBI、AFP、PIVKA）
  ps_vars_dynamic_full <- c(
    ps_vars_static_full,
    "albi_change_pre3", "afp_change_pre3", "pivka_change_pre3",
    "alt_change_pre3", "ast_change_pre3"
  )

  df <- df_arm

  make_avail <- function(vars) {
    v <- vars[vars %in% names(df)]
    v[sapply(v, function(x) sd(df[[x]], na.rm = TRUE) > 0)]
  }
  ps_vars_static_avail  <- make_avail(ps_vars_static_full)
  ps_vars_dynamic_avail <- make_avail(ps_vars_dynamic_full)

  # ── 静态变量 ──
  results_static <- list()
  plots_static   <- list()
  for (cfg in vars_config_static) {
    if (!(cfg$col %in% names(df))) next
    res <- plot_rcs_static(df, cfg$col, cfg$label, surv_time_col,
                            caption_core, title_tag, compare_label,
                            log_transform = cfg$log_transform,
                            ps_vars_base  = ps_vars_static_avail)
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
        sprintf("%s vs HAIC_alone  Harrell RCS — Baseline (IPTW)", compare_label),
        gp = grid::gpar(fontsize = 14, fontface = "bold")))
    h_bl <- ceiling(length(bl_keys) / 2) * 4.2
    ggsave(file.path(out_dir, "rms_rcs_static_baseline_combined.pdf"),
           comb, width = 11.5, height = h_bl, device = "pdf")
    ggsave(file.path(out_dir, "rms_rcs_static_baseline_combined.png"),
           comb, width = 11.5, height = h_bl, dpi = 300)
    cat("  Saved: rms_rcs_static_baseline_combined\n")
  }
  if (length(p3_keys) >= 2) {
    comb <- gridExtra::grid.arrange(
      grobs = plots_static[p3_keys], ncol = 2,
      top = grid::textGrob(
        sprintf("%s vs HAIC_alone  Harrell RCS — Pre-HAIC-3 (IPTW)", compare_label),
        gp = grid::gpar(fontsize = 14, fontface = "bold")))
    h_p3 <- ceiling(length(p3_keys) / 2) * 4.2
    ggsave(file.path(out_dir, "rms_rcs_static_pre3_combined.pdf"),
           comb, width = 11.5, height = h_p3, device = "pdf")
    ggsave(file.path(out_dir, "rms_rcs_static_pre3_combined.png"),
           comb, width = 11.5, height = h_p3, dpi = 300)
    cat("  Saved: rms_rcs_static_pre3_combined\n")
  }

  save_summary_csv(results_static, out_dir, "rms_rcs_static")

  # ── 动态变量 ──
  results_dyn <- list()
  plots_dyn   <- list()
  for (cfg in vars_config_dynamic) {
    if (!(cfg$col %in% names(df))) next
    res <- plot_rcs_dynamic(df, cfg$col, cfg$label, surv_time_col,
                             caption_core, title_tag, compare_label,
                             ps_vars_base = ps_vars_dynamic_avail)
    if (!is.null(res)) {
      results_dyn[[cfg$col]] <- c(res, list(label = cfg$label))
      plots_dyn[[cfg$col]]   <- res$plot
    }
  }
  if (length(plots_dyn) >= 2) {
    dyn_order <- vapply(vars_config_dynamic, function(z) z$col, character(1))
    dyn_order <- dyn_order[dyn_order %in% names(plots_dyn)]
    comb <- gridExtra::grid.arrange(
      grobs = plots_dyn[dyn_order], ncol = 2,
      top = grid::textGrob(
        sprintf("%s vs HAIC_alone  Harrell RCS — Dynamic change %% (IPTW)",
                compare_label),
        gp = grid::gpar(fontsize = 14, fontface = "bold")))
    h <- ceiling(length(plots_dyn) / 2) * 4.2
    ggsave(file.path(out_dir, "rms_rcs_dynamic_combined.pdf"),
           comb, width = 11.5, height = h, device = "pdf")
    ggsave(file.path(out_dir, "rms_rcs_dynamic_combined.png"),
           comb, width = 11.5, height = h, dpi = 300)
    cat("  Saved: rms_rcs_dynamic_combined\n")
  }
  save_summary_csv(results_dyn, out_dir, "rms_rcs_dynamic")

  cat("Saved outputs under:", out_dir, "\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 主循环：遍历 6 个配对
# ═══════════════════════════════════════════════════════════════════════════════
for (pair in PAIR_CONFIGS) {
  csv_path <- file.path(SCRIPT_DIR, pair$csv)
  if (!file.exists(csv_path)) {
    cat("\n[SKIP] CSV not found:", csv_path,
        "\n先运行: python3 build_all_pairs_cohorts.py\n")
    next
  }

  cat("\n", strrep("#", 70), "\n", sep = "")
  cat("# 配对:", pair$pair_label, "\n")
  cat(strrep("#", 70), "\n")

  df0 <- read.csv(csv_path, stringsAsFactors = FALSE)
  cat("Loaded:", csv_path, "| nrow =", nrow(df0), "\n")

  if (!"trt_compare" %in% names(df0)) {
    cat("[SKIP] Column trt_compare missing in", pair$csv, "\n")
    next
  }

  # 补充衍生列（防御性）
  if (!"sii_bl" %in% names(df0) && all(c("plt","nlr_bl") %in% names(df0)))
    df0$sii_bl <- df0$plt * df0$nlr_bl
  if (!"piv_bl" %in% names(df0) && all(c("mono_bl","plt","nlr_bl") %in% names(df0)))
    df0$piv_bl <- df0$mono_bl * df0$plt * df0$nlr_bl
  if (!"sii_pre3" %in% names(df0) && all(c("plt_pre3","nlr_pre3") %in% names(df0)))
    df0$sii_pre3 <- df0$plt_pre3 * df0$nlr_pre3
  if (!"piv_pre3" %in% names(df0) && all(c("mono_pre3","plt_pre3","nlr_pre3") %in% names(df0)))
    df0$piv_pre3 <- df0$mono_pre3 * df0$plt_pre3 * df0$nlr_pre3
  if (!"os_lm" %in% names(df0))
    df0$os_lm <- df0$os_months - 42 / 30.44

  df_landmark <- df0[!is.na(df0$os_lm)     & df0$os_lm     > 0, ]
  df_total    <- df0[!is.na(df0$os_months)  & df0$os_months > 0, ]
  cat("Landmark subset (os_lm > 0):", nrow(df_landmark), "rows\n")
  cat("Total OS subset (os_months > 0):", nrow(df_total), "rows\n")

  pair_out <- file.path(OUT_ROOT, pair$pair_label)

  run_one_timescale(
    df_landmark,
    out_dir       = file.path(pair_out, "landmark"),
    surv_time_col = "os_lm",
    caption_core  = "42-day landmark residual OS (primary time axis)",
    title_tag     = "[Landmark]",
    compare_label = pair$compare_label
  )
  run_one_timescale(
    df_total,
    out_dir       = file.path(pair_out, "total_os"),
    surv_time_col = "os_months",
    caption_core  = "Total OS from baseline (sensitivity)",
    title_tag     = "[Total OS]",
    compare_label = pair$compare_label
  )
}

cat("\n", strrep("=", 70), "\n", sep = "")
cat("All pairs completed. Output root:", OUT_ROOT, "\n")
cat(strrep("=", 70), "\n")
