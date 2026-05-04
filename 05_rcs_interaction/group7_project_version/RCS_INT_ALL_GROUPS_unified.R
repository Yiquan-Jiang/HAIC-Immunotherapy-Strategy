#!/usr/bin/env Rscript
# =============================================================================
# RCS_INT 统一脚本 - 所有治疗组对比的RCS非线性交互分析
# 
# 使用Overlap Weighting方法
# 支持所有6个治疗组对比：
#   1. HAIC_alone vs HAIC+I_concurrent
#   2. HAIC_alone vs HAIC+T_concurrent  
#   3. HAIC_alone vs HAIC+I+T_concurrent
#   4. HAIC_alone vs HAIC_then_I
#   5. HAIC_alone vs HAIC_then_T
#   6. HAIC_alone vs HAIC_then_I+T
#
# 用法: Rscript RCS_INT_ALL_GROUPS_unified.R [group_name]
#   group_name: I_CONC, T_CONC, I_T_CONC, THEN_I, THEN_T, THEN_IT
#   不指定则运行所有组
# =============================================================================

suppressPackageStartupMessages({
  library(survival)
  library(rms)
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(grid)
  library(glmnet)
})

# ── 配置所有治疗组对比 ──────────────────────────────────────────────────────
GROUPS_CONFIG <- list(
  I_CONC = list(
    trt_col = "trt_haic_i_conc",
    trt_label = "HAIC+I_concurrent",
    ctrl_label = "HAIC_alone",
    data_file = "RCS_INT_HAIC_ALONE_AND_I_CONC_cohort.csv",
    output_dir = "RCS_INT_HAIC_ALONE_AND_I_CONC"
  ),
  T_CONC = list(
    trt_col = "trt_haic_t_conc",
    trt_label = "HAIC+T_concurrent",
    ctrl_label = "HAIC_alone",
    data_file = "RCS_INT_HAIC_ALONE_AND_T_CONC_cohort.csv",
    output_dir = "RCS_INT_HAIC_ALONE_AND_T_CONC"
  ),
  I_T_CONC = list(
    trt_col = "trt_haic_i_t_conc",
    trt_label = "HAIC+I+T_concurrent",
    ctrl_label = "HAIC_alone",
    data_file = "RCS_INT_HAIC_ALONE_AND_I_T_CONC_cohort.csv",
    output_dir = "RCS_INT_HAIC_ALONE_AND_I_T_CONC"
  ),
  THEN_I = list(
    trt_col = "trt_haic_then_i",
    trt_label = "HAIC_then_I",
    ctrl_label = "HAIC_alone",
    data_file = "RCS_INT_HAIC_ALONE_AND_THEN_I_cohort.csv",
    output_dir = "RCS_INT_HAIC_ALONE_AND_THEN_I"
  ),
  THEN_T = list(
    trt_col = "trt_haic_then_t",
    trt_label = "HAIC_then_T",
    ctrl_label = "HAIC_alone",
    data_file = "RCS_INT_HAIC_ALONE_AND_THEN_T_cohort.csv",
    output_dir = "RCS_INT_HAIC_ALONE_AND_THEN_T"
  ),
  THEN_IT = list(
    trt_col = "trt_haic_then_it",
    trt_label = "HAIC_then_I+T",
    ctrl_label = "HAIC_alone",
    data_file = "RCS_INT_HAIC_ALONE_AND_THEN_IT_cohort.csv",
    output_dir = "RCS_INT_HAIC_ALONE_AND_THEN_IT"
  )
)

# ── 获取脚本目录 ────────────────────────────────────────────────────────────
args_all <- commandArgs(trailingOnly = FALSE)
fa <- args_all[grepl("^--file=", args_all)]
SCRIPT_DIR <- if (length(fa)) {
  dirname(normalizePath(sub("^--file=", "", fa[1]), winslash = "/"))
} else {
  getwd()
}

# ── 全局参数 ────────────────────────────────────────────────────────────────
RMS_RCS_NK <- 3L
N_BOOT <- 200L
MIN_N <- 50L

# ── Overlap Weighting 函数 ───────────────────────────────────────────────────
impute_ps_vars <- function(df, ps_vars_full) {
  for (v in ps_vars_full) {
    if (v %in% names(df) && any(is.na(df[[v]]))) {
      med_val <- median(df[[v]], na.rm = TRUE)
      df[[v]][is.na(df[[v]])] <- med_val
    }
  }
  df
}

compute_overlap_weight_static <- function(df_input, var_col, ps_vars_base, trt_col) {
  ps_vars_use <- ps_vars_base[ps_vars_base %in% names(df_input)]
  ps_vars_use <- ps_vars_use[sapply(ps_vars_use, function(v) sd(df_input[[v]], na.rm = TRUE) > 0)]
  
  df_cc <- df_input[complete.cases(df_input[, ps_vars_use, drop = FALSE]), ]
  if (nrow(df_cc) < 20) return(df_cc[0, , drop = FALSE])
  
  X_sc <- scale(as.matrix(df_cc[, ps_vars_use, drop = FALSE]))
  y <- df_cc[[trt_col]]
  cv_fit <- cv.glmnet(X_sc, y, family = "binomial", alpha = 0, nfolds = 5)
  ps_prob <- as.numeric(predict(cv_fit, newx = X_sc, s = "lambda.min", type = "response"))
  ps_prob <- pmin(pmax(ps_prob, 0.05), 0.95)
  
  df_cc$sw_raw <- (1 - ps_prob)^y * ps_prob^(1 - y)
  df_cc$sw <- df_cc$sw_raw / mean(df_cc$sw_raw)
  df_cc
}

compute_overlap_weight_dynamic <- function(df_input, var_col, ps_vars_base, trt_col) {
  extra_vars <- c("pivka_change_pre3", "afp_change_pre3", "albi_pre3")
  ps_vars_extended <- c(ps_vars_base, extra_vars)
  ps_vars_use <- ps_vars_extended[ps_vars_extended %in% names(df_input)]
  ps_vars_use <- ps_vars_use[sapply(ps_vars_use, function(v) sd(df_input[[v]], na.rm = TRUE) > 0)]
  
  df_cc <- df_input[complete.cases(df_input[, ps_vars_use, drop = FALSE]), ]
  if (nrow(df_cc) < 20) return(df_cc[0, , drop = FALSE])
  
  X_sc <- scale(as.matrix(df_cc[, ps_vars_use, drop = FALSE]))
  y <- df_cc[[trt_col]]
  cv_fit <- cv.glmnet(X_sc, y, family = "binomial", alpha = 0, nfolds = 5)
  ps_prob <- as.numeric(predict(cv_fit, newx = X_sc, s = "lambda.min", type = "response"))
  ps_prob <- pmin(pmax(ps_prob, 0.05), 0.95)
  
  df_cc$sw_raw <- (1 - ps_prob)^y * ps_prob^(1 - y)
  df_cc$sw <- df_cc$sw_raw / mean(df_cc$sw_raw)
  df_cc
}

# ── RCS 分析函数 ────────────────────────────────────────────────────────────
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

predict_hr_curve <- function(fit, x_grid, trt_col) {
  nd1 <- setNames(data.frame(1L, x_grid), c(trt_col, "rcsx"))
  nd0 <- setNames(data.frame(0L, x_grid), c(trt_col, "rcsx"))
  lp1 <- as.numeric(predict(fit, nd1, type = "lp"))
  lp0 <- as.numeric(predict(fit, nd0, type = "lp"))
  exp(lp1 - lp0)
}

# ── 绘图函数（简化版，核心逻辑与原脚本相同）────────────────────────────────
plot_rms_rcs <- function(df_sub, var_col, var_label, surv_time_col, 
                         trt_col, trt_label, ctrl_label,
                         log_transform = FALSE, ps_vars_base = NULL,
                         is_dynamic = FALSE) {
  
  # 计算权重
  if (is_dynamic) {
    df_sub <- compute_overlap_weight_dynamic(df_sub, var_col, ps_vars_base, trt_col)
  } else {
    df_sub <- compute_overlap_weight_static(df_sub, var_col, ps_vars_base, trt_col)
  }
  
  df_sub <- df_sub[!is.na(df_sub[[var_col]]), ]
  n <- nrow(df_sub)
  
  if (n < MIN_N) return(NULL)
  
  # 数据准备
  if (log_transform) {
    df_sub$rcsx <- log1p(pmax(df_sub[[var_col]], 0))
    x_label_plot <- sprintf("%s (log scale)", var_label)
  } else {
    df_sub$rcsx <- df_sub[[var_col]]
    x_label_plot <- if (is_dynamic) sprintf("%s (%%)", var_label) else var_label
  }
  
  # 拟合模型
  dd <- suppressWarnings(datadist(df_sub[, c(trt_col, "rcsx"), drop = FALSE]))
  options(datadist = dd)
  
  fml <- as.formula(paste0(
    "Surv(", surv_time_col, ", death_status) ~ ", trt_col, " * rcs(rcsx, ", RMS_RCS_NK, ")"
  ))
  
  fit <- cph(fml, data = df_sub, weights = sw, x = TRUE, y = TRUE, robust = FALSE)
  ap <- extract_rms_anova_p(fit, trt_col)
  
  # 预测HR曲线
  xlim <- quantile(df_sub$rcsx, c(0.05, 0.95))
  x_grid <- seq(xlim[1], xlim[2], length.out = 200)
  hr_curve <- predict_hr_curve(fit, x_grid, trt_col)
  
  # Bootstrap CI
  set.seed(42)
  hr_boot <- matrix(NA_real_, N_BOOT, length(x_grid))
  for (b in seq_len(N_BOOT)) {
    idx <- sample.int(nrow(df_sub), replace = TRUE)
    df_b <- df_sub[idx, ]
    tryCatch({
      dd_b <- suppressWarnings(datadist(df_b[, c(trt_col, "rcsx"), drop = FALSE]))
      options(datadist = dd_b)
      fit_b <- cph(fml, data = df_b, weights = sw, x = TRUE, y = TRUE, robust = FALSE)
      hr_boot[b, ] <- predict_hr_curve(fit_b, x_grid, trt_col)
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
      subtitle = sprintf("Overlap Weighting | Int. p=%s | Nonlin.int. p=%s | n=%d",
                         ifelse(is.na(ap$int_p), "NA", sprintf("%.3f", ap$int_p)),
                         ifelse(is.na(ap$nonlin_p), "NA", sprintf("%.3f", ap$nonlin_p)), n),
      x = x_label_plot,
      y = sprintf("HR (%s vs %s)", trt_label, ctrl_label)
    ) +
    scale_y_log10(breaks = c(0.3, 0.5, 0.7, 0.9, 1.0, 1.5, 2.0)) +
    coord_cartesian(ylim = c(0.25, 2.5), xlim = xlim) +
    theme_bw(base_size = 13)
  
  list(plot = p, int_p = ap$int_p, nonlin_p = ap$nonlin_p, n = n)
}

# ── 主分析函数 ──────────────────────────────────────────────────────────────
run_group_analysis <- function(group_name, config) {
  
  cat("\n", strrep("=", 80), "\n", sep = "")
  cat("分析组:", group_name, "|", config$trt_label, "vs", config$ctrl_label, "\n")
  cat(strrep("=", 80), "\n")
  
  # 加载数据
  data_path <- file.path(SCRIPT_DIR, "data", config$data_file)
  if (!file.exists(data_path)) {
    cat("错误: 数据文件不存在:", data_path, "\n")
    return(NULL)
  }
  
  df0 <- read.csv(data_path, stringsAsFactors = FALSE)
  cat("加载数据:", config$data_file, "| n =", nrow(df0), "\n")
  
  # 准备变量
  ps_vars_full <- c("albi_bl", "alb_bl", "tbil_bl", "inr", "plt",
                    "tumor_max_diameter_cm", "tumor_count_enc",
                    "pvtt_grade", "ascites_score_enc",
                    "log_afp_bl", "log_pivka_bl",
                    "metastasis_binary", "lymph_node_binary",
                    "neut_bl", "lymph_bl", "mono_bl")
  
  df <- impute_ps_vars(df0, ps_vars_full)
  ps_vars_available <- ps_vars_full[ps_vars_full %in% names(df)]
  ps_vars_available <- ps_vars_available[sapply(ps_vars_available, 
                                                function(v) sd(df[[v]], na.rm = TRUE) > 0)]
  
  # 创建输出目录
  out_dir <- file.path(SCRIPT_DIR, "output", config$output_dir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 准备数据子集
  if (!"os_lm" %in% names(df)) {
    df$os_lm <- df$os_months - 42 / 30.44
  }
  
  df_landmark <- df[!is.na(df$os_lm) & df$os_lm > 0, ]
  df_total <- df[!is.na(df$os_months) & df$os_months > 0, ]
  
  # 定义变量配置
  vars_static <- list(
    list(col = "afp", label = "AFP (baseline, ng/mL)", log = TRUE),
    list(col = "pivka", label = "PIVKA-II (baseline, mAU/mL)", log = TRUE),
    list(col = "nlr_bl", label = "NLR (baseline)", log = FALSE),
    list(col = "albi_bl", label = "ALBI score (baseline)", log = FALSE),
    list(col = "plr_bl", label = "PLR (baseline)", log = FALSE),
    list(col = "sii_bl", label = "SII (baseline)", log = TRUE),
    list(col = "tumor_max_diameter_cm", label = "Tumor diameter (cm)", log = FALSE)
  )
  
  vars_dynamic <- list(
    list(col = "nlr_change_pre3", label = "NLR change rate"),
    list(col = "afp_change_pre3", label = "AFP change rate"),
    list(col = "pivka_change_pre3", label = "PIVKA-II change rate"),
    list(col = "plr_change_pre3", label = "PLR change rate"),
    list(col = "albi_change_pre3", label = "ALBI change rate")
  )
  
  # 分析函数
  analyze_subset <- function(df_subset, subset_name, surv_time_col) {
    
    cat("\n---", subset_name, "---\n")
    
    results <- list()
    
    # 静态变量
    for (v in vars_static) {
      if (!v$col %in% names(df_subset)) next
      cat("分析:", v$label, "...")
      res <- plot_rms_rcs(df_subset, v$col, v$label, surv_time_col,
                          config$trt_col, config$trt_label, config$ctrl_label,
                          v$log, ps_vars_available, is_dynamic = FALSE)
      if (!is.null(res)) {
        results[[paste0("static_", v$col)]] <- res
        cat(" n =", res$n, "int.p =", sprintf("%.4f", res$int_p), "\n")
      } else {
        cat(" SKIP\n")
      }
    }
    
    # 动态变量
    for (v in vars_dynamic) {
      if (!v$col %in% names(df_subset)) next
      cat("分析:", v$label, "...")
      res <- plot_rms_rcs(df_subset, v$col, v$label, surv_time_col,
                          config$trt_col, config$trt_label, config$ctrl_label,
                          FALSE, ps_vars_available, is_dynamic = TRUE)
      if (!is.null(res)) {
        results[[paste0("dynamic_", v$col)]] <- res
        cat(" n =", res$n, "int.p =", sprintf("%.4f", res$int_p), "\n")
      } else {
        cat(" SKIP\n")
      }
    }
    
    # 保存结果摘要
    if (length(results) > 0) {
      summary_df <- do.call(rbind, lapply(names(results), function(k) {
        r <- results[[k]]
        data.frame(
          variable = k,
          n = r$n,
          interaction_p = round(r$int_p, 4),
          nonlinear_interaction_p = round(r$nonlin_p, 4),
          stringsAsFactors = FALSE
        )
      }))
      
      write.csv(summary_df, 
                file.path(out_dir, paste0(subset_name, "_anova_summary.csv")),
                row.names = FALSE)
    }
    
    results
  }
  
  # 运行分析
  landmark_results <- analyze_subset(df_landmark, "landmark", "os_lm")
  total_results <- analyze_subset(df_total, "total_os", "os_months")
  
  cat("\n输出目录:", out_dir, "\n")
  
  list(landmark = landmark_results, total_os = total_results)
}

# ── 主程序 ──────────────────────────────────────────────────────────────────
main <- function() {
  
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
    # 运行所有组
    cat("\n将运行所有6个治疗组对比分析...\n")
    
    for (group_name in names(GROUPS_CONFIG)) {
      tryCatch({
        run_group_analysis(group_name, GROUPS_CONFIG[[group_name]])
      }, error = function(e) {
        cat("\n错误 [", group_name, "]:", conditionMessage(e), "\n")
      })
    }
    
    cat("\n所有分析完成！\n")
  }
}

# 执行
main()