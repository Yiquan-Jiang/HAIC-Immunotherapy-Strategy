#!/usr/bin/env Rscript
# =============================================================================
# RCS_INT_HAIC_ALONE_AND_THEN_T — rms::rcs() + IPTW 加权 rms::cph()，双时间尺度
#
# 队列: TIDY main_group 为 HAIC_alone vs HAIC_then_T；暴露项
#   Surv(...) ~ trt_haic_then_t * rcs(rcsx, nk)
#
# 依赖: survival, rms, Hmisc, ggplot2, dplyr, glmnet, gridExtra, grid
#
# 数据: RCS_INT_HAIC_ALONE_AND_THEN_T_cohort.csv（先运行 RCS_INT_HAIC_ALONE_AND_THEN_T_build_cohort.py）
# 输出: output/RCS_INT_HAIC_ALONE_AND_THEN_T/{landmark,total_os}/
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
  unset = normalizePath(file.path(SCRIPT_DIR, "data", "RCS_INT_HAIC_ALONE_AND_THEN_T_cohort.csv"), mustWork = FALSE)
)
BASE_OUT <- Sys.getenv(
  "RMS_RCS_OUT_DIR",
  unset = file.path(SCRIPT_DIR, "output", "RCS_INT_HAIC_ALONE_AND_THEN_T")
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

# 治疗变量列名与标签（可通过环境变量覆盖，默认 HAIC_alone vs HAIC_then_T）
TRT_COL   <- Sys.getenv("RMS_RCS_TRT_COL",   unset = "trt_haic_then_t")
TRT_LABEL <- Sys.getenv("RMS_RCS_TRT_LABEL", unset = "HAIC_then_T")
CTRL_LABEL<- Sys.getenv("RMS_RCS_CTRL_LABEL",unset = "HAIC_alone")

cat("TRT_COL:", TRT_COL, "| TRT_LABEL:", TRT_LABEL, "| CTRL_LABEL:", CTRL_LABEL, "\n")

# ── IPTW（glmnet ridge PS + 稳定权重）────────────────────────────────────────
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

