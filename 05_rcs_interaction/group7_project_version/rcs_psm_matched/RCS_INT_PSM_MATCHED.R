#!/usr/bin/env Rscript
# =============================================================================
# RCS_INT_PSM_MATCHED - 基于PSM匹配队列的RCS非线性交互分析
#
# 与原脚本 RCS_INT_ALL_GROUPS_unified.R 的区别:
#   1. 人群: PSM 1:1 匹配后的队列 (不再使用全队列)
#   2. 加权: 去除 Overlap Weighting / IPTW (PSM已平衡混杂)
#   3. 依赖: 去除 glmnet 包
#
# 支持所有6个治疗组对比（均以 HAIC_alone 为对照组）:
#   1. I_CONC   - HAIC_alone vs HAIC+I_concurrent
#   2. T_CONC   - HAIC_alone vs HAIC+T_concurrent
#   3. I_T_CONC - HAIC_alone vs HAIC+I+T_concurrent
#   4. THEN_I   - HAIC_alone vs HAIC_then_I
#   5. THEN_T   - HAIC_alone vs HAIC_then_T
#   6. THEN_IT  - HAIC_alone vs HAIC_then_I+T
#
# 用法: Rscript RCS_INT_PSM_MATCHED.R [group_name]
#   group_name: I_CONC, T_CONC, I_T_CONC, THEN_I, THEN_T, THEN_IT
#   不指定则运行所有组
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(rms)
  library(ggplot2)
  library(dplyr)
})

# ── 配置所有治疗组对比 ──────────────────────────────────────────────────────
# matched_ids_file: 相对于 PSM_MATCHED_IDS_DIR 的路径
# data_file: 相对于脚本所在目录下 data/ 的路径
GROUPS_CONFIG <- list(
  I_CONC = list(
    trt_label = "HAIC+I_concurrent",
    ctrl_label = "HAIC_alone",
    data_file = "RCS_INT_HAIC_ALONE_AND_I_CONC_cohort.csv",
    matched_ids_file = "matched_ids_01_HAIC_alone_vs_HAIC+I_concurrent.csv",
    output_dir = "I_CONC"
  ),
  T_CONC = list(
    trt_label = "HAIC+T_concurrent",
    ctrl_label = "HAIC_alone",
    data_file = "RCS_INT_HAIC_ALONE_AND_T_CONC_cohort.csv",
    matched_ids_file = "matched_ids_03_HAIC_alone_vs_HAIC+T_concurrent.csv",
    output_dir = "T_CONC"
  ),
  I_T_CONC = list(
    trt_label = "HAIC+I+T_concurrent",
    ctrl_label = "HAIC_alone",
    data_file = "RCS_INT_HAIC_ALONE_AND_I_T_CONC_cohort.csv",
    matched_ids_file = "matched_ids_05_HAIC_alone_vs_HAIC+I+T_concurrent.csv",
    output_dir = "I_T_CONC"
  ),
  THEN_I = list(
    trt_label = "HAIC_then_I",
    ctrl_label = "HAIC_alone",
    data_file = "RCS_INT_HAIC_ALONE_AND_THEN_I_cohort.csv",
    matched_ids_file = "matched_ids_02_HAIC_alone_vs_HAIC_then_I.csv",
    output_dir = "THEN_I"
  ),
  THEN_T = list(
    trt_label = "HAIC_then_T",
    ctrl_label = "HAIC_alone",
    data_file = "RCS_INT_HAIC_ALONE_AND_THEN_T_cohort.csv",
    matched_ids_file = "matched_ids_04_HAIC_alone_vs_HAIC_then_T.csv",
    output_dir = "THEN_T"
  ),
  THEN_IT = list(
    trt_label = "HAIC_then_I+T",
    ctrl_label = "HAIC_alone",
    data_file = "RCS_INT_HAIC_ALONE_AND_THEN_IT_cohort.csv",
    matched_ids_file = "matched_ids_06_HAIC_alone_vs_HAIC_then_I+T.csv",
    output_dir = "THEN_IT"
  )
)

# ── 路径配置 ────────────────────────────────────────────────────────────────
args_all <- commandArgs(trailingOnly = FALSE)
fa <- args_all[grepl("^--file=", args_all)]
SCRIPT_DIR <- if (length(fa)) {
  dirname(normalizePath(sub("^--file=", "", fa[1]), winslash = "/"))
} else {
  getwd()
}

DATA_DIR <- file.path(SCRIPT_DIR, "data")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "output")

# PSM matched_ids 所在目录 (项目级)
PROJECT_ROOT <- file.path(SCRIPT_DIR, "..", "..", "..")
PSM_MATCHED_IDS_DIR <- file.path(PROJECT_ROOT, "results", "psm_balance_tables_complete")

# ── 全局参数 ────────────────────────────────────────────────────────────────
RMS_RCS_NK <- 3L       # RCS 节点数
N_BOOT     <- 200L     # Bootstrap 重抽样次数
MIN_N      <- 50L      # 最小样本量（低于此跳过）
LANDMARK_DAYS <- 42L   # Landmark 分析截止天数

# ── RCS 分析工具函数 ────────────────────────────────────────────────────────

#' 从 rms::anova() 提取交互 p 值
#' 返回: list(int_p, nonlin_p)
extract_rms_anova_p <- function(fit) {
  a <- anova(fit)
  rn <- rownames(a)
  pcol <- if ("P" %in% colnames(a)) "P" else ncol(a)
  
  int_p <- NA_real_
  nonlin_p <- NA_real_
  
  # 总交互 p 值: "treatment * biomarker  (Factor+Higher Order Factors)"
  idx_int <- grep("treatment \\* [^\\(]+  \\(Factor\\+Higher Order Factors\\)$", rn)
  if (length(idx_int)) int_p <- parse_anova_p(a[idx_int[1], pcol])
  
  # 非线性交互 p 值: "Nonlinear Interaction : f(A,B) vs. AB"
  idx_nl <- grep("Nonlinear Interaction : f\\(A,B\\) vs\\. AB", rn)
  if (length(idx_nl)) nonlin_p <- parse_anova_p(a[idx_nl[1], pcol])
  
  list(int_p = int_p, nonlin_p = nonlin_p)
}

#' 解析 ANOVA p 值 (处理 "<0.001" 等格式)
parse_anova_p <- function(cell) {
  x <- as.character(cell)[1]
  if (is.na(x) || !nzchar(x)) return(NA_real_)
  if (grepl("^<", x)) {
    z <- sub("^<\\.?", "", x)
    return(suppressWarnings(as.numeric(z)))
  }
  suppressWarnings(as.numeric(x))
}

#' 预测 HR 曲线: trt=1 vs trt=0
predict_hr_curve <- function(fit, x_grid) {
  nd1 <- setNames(data.frame(1L, x_grid), c("treatment", "rcsx"))
  nd0 <- setNames(data.frame(0L, x_grid), c("treatment", "rcsx"))
  lp1 <- as.numeric(predict(fit, nd1, type = "lp"))
  lp0 <- as.numeric(predict(fit, nd0, type = "lp"))
  exp(lp1 - lp0)
}

# ── 核心: PSM匹配后队列的 RCS 分析 + 绘图 ──────────────────────────────────

#' 在 PSM 匹配队列上拟合 RCS Cox 模型并绘图
#'
#' @param df_sub    PSM匹配后的数据框
#' @param var_col   生物标志物列名
#' @param var_label 显示标签
#' @param surv_time_col 生存时间列名
#' @param trt_label 治疗组标签
#' @param ctrl_label 对照组标签
#' @param log_transform 是否 log1p 转换
#' @param is_dynamic 是否为动态变化率变量
#' @param output_dir 输出目录路径
#'
#' @return list(plot, int_p, nonlin_p, n) 或 NULL (样本量不足)
plot_rms_rcs_psm <- function(df_sub, var_col, var_label, surv_time_col,
                             trt_label, ctrl_label,
                             log_transform = FALSE, is_dynamic = FALSE,
                             output_dir = NULL) {
  
  # 去除缺失值
  df_sub <- df_sub[!is.na(df_sub[[var_col]]) & !is.na(df_sub$treatment) &
                   !is.na(df_sub[[surv_time_col]]) & !is.na(df_sub$death_status), ]
  n <- nrow(df_sub)
  
  if (n < MIN_N) return(NULL)
  
  # 变量转换
  if (log_transform) {
    df_sub$rcsx <- log1p(pmax(df_sub[[var_col]], 0))
    x_label_plot <- sprintf("%s (log scale)", var_label)
  } else {
    df_sub$rcsx <- df_sub[[var_col]]
    x_label_plot <- if (is_dynamic) sprintf("%s (%%)", var_label) else var_label
  }
  
  # datadist 设置
  dd <- suppressWarnings(datadist(df_sub[, c("treatment", "rcsx"), drop = FALSE]))
  options(datadist = dd)
  
  # 拟合 RCS Cox 模型 (无权重)
  fml <- as.formula(paste0(
    "Surv(", surv_time_col, ", death_status) ~ treatment * rcs(rcsx, ", RMS_RCS_NK, ")"
  ))
  
  fit <- cph(fml, data = df_sub, x = TRUE, y = TRUE, robust = FALSE)
  ap <- extract_rms_anova_p(fit)
  
  # 预测 HR 曲线 (在 5th-95th 百分位范围内)
  xlim <- quantile(df_sub$rcsx, c(0.05, 0.95))
  x_grid <- seq(xlim[1], xlim[2], length.out = 200)
  hr_curve <- predict_hr_curve(fit, x_grid)
  
  # Bootstrap 95% CI
  set.seed(42)
  hr_boot <- matrix(NA_real_, N_BOOT, length(x_grid))
  for (b in seq_len(N_BOOT)) {
    idx <- sample.int(n, replace = TRUE)
    df_b <- df_sub[idx, ]
    tryCatch({
      dd_b <- suppressWarnings(datadist(df_b[, c("treatment", "rcsx"), drop = FALSE]))
      options(datadist = dd_b)
      fit_b <- cph(fml, data = df_b, x = TRUE, y = TRUE, robust = FALSE)
      hr_boot[b, ] <- predict_hr_curve(fit_b, x_grid)
    }, error = function(e) NULL)
  }
  
  hr_lo <- apply(hr_boot, 2, quantile, 0.025, na.rm = TRUE)
  hr_hi <- apply(hr_boot, 2, quantile, 0.975, na.rm = TRUE)
  
  # 绘图
  plot_df <- data.frame(x = x_grid, hr = hr_curve, hr_lo = hr_lo, hr_hi = hr_hi)
  
  p <- ggplot(plot_df, aes(x = x)) +
    geom_ribbon(aes(ymin = hr_lo, ymax = hr_hi), fill = "#3C5488", alpha = 0.18) +
    geom_line(aes(y = hr), color = "#3C5488", linewidth = 1.8) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "#333333", linewidth = 1.0) +
    labs(
      title = var_label,
      subtitle = sprintf("PSM Matched | Int. p=%s | Nonlin.int. p=%s | n=%d",
                         ifelse(is.na(ap$int_p), "NA", format_p(ap$int_p)),
                         ifelse(is.na(ap$nonlin_p), "NA", format_p(ap$nonlin_p)), n),
      x = x_label_plot,
      y = sprintf("HR (%s vs %s)", trt_label, ctrl_label)
    ) +
    scale_y_log10(breaks = c(0.3, 0.5, 0.7, 0.9, 1.0, 1.5, 2.0)) +
    coord_cartesian(ylim = c(0.25, 2.5), xlim = xlim) +
    theme_bw(base_size = 13)
  
  # 保存图片
  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # 生成文件名: 去除特殊字符
    safe_name <- gsub("[^a-zA-Z0-9_-]", "_", var_col)
    prefix <- if (grepl("^landmark", surv_time_col)) "landmark" else "total_os"
    
    ggsave(file.path(output_dir, sprintf("%s_%s.png", prefix, safe_name)),
           p, width = 8, height = 6, dpi = 300)
    ggsave(file.path(output_dir, sprintf("%s_%s.pdf", prefix, safe_name)),
           p, width = 8, height = 6)
  }
  
  list(plot = p, int_p = ap$int_p, nonlin_p = ap$nonlin_p, n = n)
}

#' 格式化 p 值
format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# ── 主分析函数 ──────────────────────────────────────────────────────────────

run_group_analysis <- function(group_name, config) {
  
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("分析组:", group_name, "|", config$trt_label, "vs", config$ctrl_label, "\n")
  cat(strrep("=", 80), "\n")
  
  # ── 1. 加载原始 cohort 数据 ──
  cohort_path <- file.path(DATA_DIR, config$data_file)
  if (!file.exists(cohort_path)) {
    cat("错误: cohort 数据不存在:", cohort_path, "\n")
    return(NULL)
  }
  df_cohort <- read.csv(cohort_path, stringsAsFactors = FALSE)
  cat("Cohort 数据:", config$data_file, "| n =", nrow(df_cohort), "\n")
  
  # ── 2. 加载 PSM matched_ids ──
  matched_path <- file.path(PSM_MATCHED_IDS_DIR, config$matched_ids_file)
  if (!file.exists(matched_path)) {
    cat("错误: matched_ids 不存在:", matched_path, "\n")
    return(NULL)
  }
  df_matched <- read.csv(matched_path, stringsAsFactors = FALSE)
  cat("PSM matched_ids:", config$matched_ids_file,
      "| n =", nrow(df_matched), "\n")
  
  # ── 3. 合并: 通过 patient_id 筛选匹配队列 ──
  # matched_ids 中 treatment 列: 0=对照组(HAIC_alone), 1=实验组
  df_psm <- merge(df_cohort, df_matched[, c("patient_id", "treatment")],
                  by = "patient_id", all.x = FALSE, all.y = FALSE)
  cat("PSM 匹配后队列: n =", nrow(df_psm),
      "| 对照组 n =", sum(df_psm$treatment == 0),
      "| 实验组 n =", sum(df_psm$treatment == 1), "\n")
  
  if (nrow(df_psm) < MIN_N) {
    cat("警告: 匹配后样本量不足 (n=", nrow(df_psm), "), 跳过\n", sep = "")
    return(NULL)
  }
  
  # ── 4. 准备生存时间 ──
  if (!"os_lm" %in% names(df_psm)) {
    df_psm$os_lm <- df_psm$os_months - LANDMARK_DAYS / 30.44
  }
  df_landmark <- df_psm[!is.na(df_psm$os_lm) & df_psm$os_lm > 0, ]
  df_total    <- df_psm[!is.na(df_psm$os_months) & df_psm$os_months > 0, ]
  cat("Landmark 队列: n =", nrow(df_landmark), "\n")
  cat("Total OS 队列: n =", nrow(df_total), "\n")
  
  # ── 5. 创建输出目录 ──
  out_dir <- file.path(OUTPUT_DIR, config$output_dir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ── 6. 定义分析变量 ──
  # 静态基线变量 (14 baseline + 12 pre-HAIC-3 = 26)
  vars_static <- list(
    # --- Baseline ---
    list(col = "afp",                  label = "AFP (baseline, ng/mL)",                  log = TRUE,  group = "baseline"),
    list(col = "pivka",                label = "PIVKA-II (baseline, mAU/mL)",            log = TRUE,  group = "baseline"),
    list(col = "nlr_bl",               label = "NLR (baseline)",                         log = FALSE, group = "baseline"),
    list(col = "albi_bl",              label = "ALBI score (baseline)",                  log = FALSE, group = "baseline"),
    list(col = "alb_bl",               label = "Albumin (baseline, g/L)",                log = FALSE, group = "baseline"),
    list(col = "tbil_bl",              label = "Total bilirubin (baseline, \u00b5mol/L)", log = TRUE,  group = "baseline"),
    list(col = "plr_bl",               label = "PLR (baseline)",                         log = FALSE, group = "baseline"),
    list(col = "sii_bl",               label = "SII (PLT\u00d7NLR, baseline)",           log = TRUE,  group = "baseline"),
    list(col = "piv_bl",               label = "PIV (Mono\u00d7PLT\u00d7NLR, baseline)", log = TRUE,  group = "baseline"),
    list(col = "neut_bl",              label = "Neutrophil (baseline, 10^9/L)",          log = FALSE, group = "baseline"),
    list(col = "lymph_bl",             label = "Lymphocyte (baseline, 10^9/L)",          log = FALSE, group = "baseline"),
    list(col = "mono_bl",              label = "Monocyte (baseline, 10^9/L)",            log = FALSE, group = "baseline"),
    list(col = "plt",                  label = "Platelet (baseline, 10^9/L)",            log = TRUE,  group = "baseline"),
    list(col = "tumor_max_diameter_cm", label = "Tumor diameter (baseline, cm)",          log = FALSE, group = "baseline"),
    # --- Pre-HAIC-3 ---
    list(col = "afp_pre3",             label = "AFP (pre-HAIC-3, ng/mL)",               log = TRUE,  group = "pre3_static"),
    list(col = "pivka_pre3",           label = "PIVKA-II (pre-HAIC-3, mAU/mL)",         log = TRUE,  group = "pre3_static"),
    list(col = "nlr_pre3",             label = "NLR (pre-HAIC-3)",                       log = FALSE, group = "pre3_static"),
    list(col = "albi_pre3",            label = "ALBI score (pre-HAIC-3)",                log = FALSE, group = "pre3_static"),
    list(col = "alb_pre3",             label = "Albumin (pre-HAIC-3, g/L)",              log = FALSE, group = "pre3_static"),
    list(col = "tbil_pre3",            label = "Total bilirubin (pre-HAIC-3, \u00b5mol/L)", log = TRUE, group = "pre3_static"),
    list(col = "plr_pre3",             label = "PLR (pre-HAIC-3)",                       log = FALSE, group = "pre3_static"),
    list(col = "mono_pre3",            label = "Monocyte (pre-HAIC-3, 10^9/L)",          log = FALSE, group = "pre3_static"),
    list(col = "plt_pre3",             label = "Platelet (pre-HAIC-3, 10^9/L)",          log = TRUE,  group = "pre3_static"),
    list(col = "sii_pre3",             label = "SII (PLT\u00d7NLR, pre-HAIC-3)",         log = TRUE,  group = "pre3_static"),
    list(col = "piv_pre3",             label = "PIV (Mono\u00d7PLT\u00d7NLR, pre-HAIC-3)", log = TRUE, group = "pre3_static"),
    list(col = "neut_pre3",            label = "Neutrophil (pre-HAIC-3, 10^9/L)",        log = FALSE, group = "pre3_static"),
    list(col = "lymph_pre3",           label = "Lymphocyte (pre-HAIC-3, 10^9/L)",        log = FALSE, group = "pre3_static")
  )

  # 动态变化率变量 (13)
  vars_dynamic <- list(
    list(col = "nlr_change_pre3",    label = "NLR change rate"),
    list(col = "sii_change_pre3",    label = "SII change rate"),
    list(col = "piv_change_pre3",    label = "PIV change rate"),
    list(col = "afp_change_pre3",    label = "AFP change rate"),
    list(col = "pivka_change_pre3",  label = "PIVKA-II change rate"),
    list(col = "neut_change_pre3",   label = "Neutrophil change rate"),
    list(col = "plr_change_pre3",    label = "PLR change rate"),
    list(col = "albi_change_pre3",   label = "ALBI change rate"),
    list(col = "alb_change_pre3",    label = "Albumin change rate"),
    list(col = "tbil_change_pre3",   label = "TBIL change rate"),
    list(col = "lymph_change_pre3",  label = "Lymphocyte change rate"),
    list(col = "mono_change_pre3",   label = "Monocyte change rate"),
    list(col = "plt_change_pre3",    label = "Platelet change rate")
  )
  
  # ── 7. 分析函数 ──
  analyze_subset <- function(df_subset, subset_name, surv_time_col) {
    
    cat("\n---", subset_name, "(n =", nrow(df_subset), ") ---\n")
    
    results <- list()
    
    # 静态变量
    for (v in vars_static) {
      if (!v$col %in% names(df_subset)) next
      cat("  [", v$group, "]", v$label, "...")
      res <- plot_rms_rcs_psm(
        df_subset, v$col, v$label, surv_time_col,
        config$trt_label, config$ctrl_label,
        log_transform = v$log, is_dynamic = FALSE,
        output_dir = out_dir
      )
      if (!is.null(res)) {
        results[[paste0("static_", v$col)]] <- c(res, list(label = v$label, group = v$group))
        cat(" n =", res$n,
            "int.p =", format_p(res$int_p),
            "nonlin.p =", format_p(res$nonlin_p), "\n")
      } else {
        cat(" SKIP (n < ", MIN_N, " or NA)\n")
      }
    }
    
    # 动态变量
    for (v in vars_dynamic) {
      if (!v$col %in% names(df_subset)) next
      cat("  [dynamic]", v$label, "...")
      res <- plot_rms_rcs_psm(
        df_subset, v$col, v$label, surv_time_col,
        config$trt_label, config$ctrl_label,
        log_transform = FALSE, is_dynamic = TRUE,
        output_dir = out_dir
      )
      if (!is.null(res)) {
        results[[paste0("dynamic_", v$col)]] <- c(res, list(label = v$label, group = "dynamic"))
        cat(" n =", res$n,
            "int.p =", format_p(res$int_p),
            "nonlin.p =", format_p(res$nonlin_p), "\n")
      } else {
        cat(" SKIP (n < ", MIN_N, " or NA)\n")
      }
    }
    
    # 保存 ANOVA 摘要 CSV（含 label 和 group 列）
    if (length(results) > 0) {
      summary_df <- do.call(rbind, lapply(names(results), function(k) {
        r <- results[[k]]
        data.frame(
          variable = k,
          label = r$label,
          group = r$group,
          n = r$n,
          interaction_p = round(r$int_p, 4),
          nonlinear_interaction_p = round(r$nonlin_p, 4),
          stringsAsFactors = FALSE
        )
      }))
      
      csv_path <- file.path(out_dir, paste0(subset_name, "_anova_summary.csv"))
      write.csv(summary_df, csv_path, row.names = FALSE)
      cat("  摘要已保存:", csv_path, "\n")
    }
    
    results
  }
  
  # ── 8. 运行两种生存框架的分析 ──
  landmark_results <- analyze_subset(df_landmark, "landmark", "os_lm")
  total_results    <- analyze_subset(df_total,    "total_os", "os_months")
  
  cat("\n输出目录:", out_dir, "\n")
  
  list(landmark = landmark_results, total_os = total_results)
}

# ── 主程序 ──────────────────────────────────────────────────────────────────

main <- function() {
  
  # 验证关键路径
  if (!dir.exists(DATA_DIR)) {
    stop("数据目录不存在: ", DATA_DIR)
  }
  if (!dir.exists(PSM_MATCHED_IDS_DIR)) {
    stop("PSM matched_ids 目录不存在: ", PSM_MATCHED_IDS_DIR)
  }
  
  cat("=== RCS非线性交互分析 (PSM匹配队列) ===\n")
  cat("数据目录:", DATA_DIR, "\n")
  cat("PSM IDs目录:", PSM_MATCHED_IDS_DIR, "\n")
  cat("输出目录:", OUTPUT_DIR, "\n")
  cat("RCS节点数:", RMS_RCS_NK, "| Bootstrap次数:", N_BOOT, "| 最小n:", MIN_N, "\n")
  
  # 解析命令行参数
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) > 0) {
    # 运行指定的组
    group_name <- toupper(args[1])
    if (!group_name %in% names(GROUPS_CONFIG)) {
      cat("错误: 未知组名:", group_name, "\n")
      cat("可用组:", paste(names(GROUPS_CONFIG), collapse = ", "), "\n")
      quit(status = 1)
    }
    
    run_group_analysis(group_name, GROUPS_CONFIG[[group_name]])
    
  } else {
    # 运行所有 6 组
    cat("\n将运行所有6个治疗组对比分析...\n")
    
    for (group_name in names(GROUPS_CONFIG)) {
      tryCatch({
        run_group_analysis(group_name, GROUPS_CONFIG[[group_name]])
      }, error = function(e) {
        cat("\n错误 [", group_name, "]:", conditionMessage(e), "\n")
      })
    }
    
    cat("\n========== 所有分析完成! ==========\n")
    cat("结果保存在:", OUTPUT_DIR, "\n")
  }
}

# 执行
main()
