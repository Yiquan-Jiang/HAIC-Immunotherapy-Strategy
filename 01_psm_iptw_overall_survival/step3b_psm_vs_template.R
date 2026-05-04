#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
#
# CBPS-IPTW 多组加权 — WeightIt (Covariate Balancing Propensity Score)
# 估计量: ATE，同时平衡全部 7 个治疗组的协变量分布
# 参照组(HR): HAIC_alone
#
# 输出:
#   results/psm_vs_template/iptw_weights.csv
#   results/psm_vs_template/iptw_balance_summary.csv
#   results/psm_vs_template/survival_gps_final.csv
#   results/psm_vs_template/iptw_km_data.csv         (加权 KM 数据供 Python 绘图)
#   results/psm_vs_template/iptw_global_test.csv      (全模型 Wald P + PH 检验)

library(tidyverse)
library(WeightIt)
library(cobalt)
library(survival)
library(survey)

BASE_DIR   <- "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
DATA_DIR   <- file.path(BASE_DIR, "data")
OUTPUT_DIR <- file.path(BASE_DIR, "results", "psm_vs_template")
LOG_DIR    <- file.path(BASE_DIR, "logs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR,    showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(LOG_DIR, "psm_vs_template.log")
sink(log_file, split = TRUE)

cat("============================================================\n")
cat("CBPS-IPTW 多组加权 — WeightIt\n")
cat("估计量: ATE | 参照组(HR): HAIC_alone\n")
cat("============================================================\n")

REF_GROUP_HR <- "HAIC_alone"

GROUP_ORDER <- c(
  "HAIC_alone", "HAIC+I_concurrent", "HAIC_then_I",
  "HAIC+T_concurrent", "HAIC_then_T",
  "HAIC+I+T_concurrent", "HAIC_then_I+T"
)

# ════════════════════════════════════════════════════════════════════
# 1. 读取数据
# ════════════════════════════════════════════════════════════════════
cat("\n1. 读取数据...\n")

analysis_data <- read_csv(
  file.path(DATA_DIR, "analysis_ready.csv"),
  show_col_types = FALSE
) %>%
  filter(os_months >= 0) %>%
  mutate(
    group = factor(main_group, levels = GROUP_ORDER),
    death_status = case_when(
      death_status %in% c("Yes", "1", "TRUE", "yes") ~ 1L,
      death_status %in% c("No",  "0", "FALSE","no")  ~ 0L,
      TRUE ~ as.integer(as.numeric(death_status))
    )
  ) %>%
  filter(!is.na(group))

cat(sprintf("   有效分析患者: %d\n", nrow(analysis_data)))
cat("   分组分布:\n")
print(table(analysis_data$group))

# ════════════════════════════════════════════════════════════════════
# 2. 构建匹配变量
#    [C5] 删除 pvtt_present（与 pvtt_grade_cat 完全共线）
#    [C4] 有序分类变量转 factor
#    [m2] 删除 pivka_cat（与 pivka_log 冗余）
# ════════════════════════════════════════════════════════════════════
cat("\n2. 构建匹配变量...\n")

analysis_data <- analysis_data %>%
  mutate(
    afp_cat = factor(case_when(
      afp < 20   ~ 0L,
      afp < 400  ~ 1L,
      TRUE       ~ 2L
    )),
    pivka_log = log10(pmax(pivka, 0.01) + 1),
    tbil_log  = log10(pmax(tbil,  0.01) + 1),
    tumor_gt10cm   = if_else(tumor_max_diameter_cm > 10, 1L, 0L),
    tumor_multiple = if_else(tumor_count_category == ">3", 1L, 0L),
    pvtt_grade_cat = factor(case_when(
      pvtt_classification == "Absent" ~ 0L,
      pvtt_classification == "Vp1/2"  ~ 1L,
      pvtt_classification == "Vp3/4"  ~ 2L,
      TRUE                            ~ 0L
    )),
    hvtt_present     = if_else(hvtt == "Yes", 1L, 0L),
    ivc_ra_present   = if_else(ivc_or_ra_thrombus == "Yes", 1L, 0L),
    distant_meta_bin = if_else(distant_metastasis == "Yes", 1L, 0L),
    lymph_meta_bin   = if_else(lymph_node_metastasis == "Yes", 1L, 0L),
    ascites_bin      = if_else(ascites != "Absent", 1L, 0L),
    varices_bin      = if_else(varices == "Yes", 1L, 0L),
    albi_grade_fac   = factor(as.integer(albi_grade))
  )

iptw_formula <- group ~ afp_cat + pivka_log +
  tumor_gt10cm + tumor_multiple +
  pvtt_grade_cat + hvtt_present +
  ivc_ra_present + distant_meta_bin + lymph_meta_bin +
  ascites_bin + varices_bin +
  albi_grade_fac + tbil_log + alb + plt + age +
  tumor_max_diameter_cm + nlr

cov_vars <- c("afp_cat", "pivka_log",
              "tumor_gt10cm", "tumor_multiple",
              "pvtt_grade_cat", "hvtt_present",
              "ivc_ra_present", "distant_meta_bin", "lymph_meta_bin",
              "ascites_bin", "varices_bin",
              "albi_grade_fac", "tbil_log", "alb", "plt", "age",
              "tumor_max_diameter_cm", "nlr")

miss_count <- sum(!complete.cases(analysis_data[, cov_vars]))
if (miss_count > 0) {
  cat(sprintf("   ⚠ %d 行有缺失值，将被剔除\n", miss_count))
  analysis_data <- analysis_data %>% filter(complete.cases(across(all_of(cov_vars))))
  cat(sprintf("   剔除后剩余: %d\n", nrow(analysis_data)))
}

cat("   IPTW 公式: ", deparse(iptw_formula, width.cutoff = 200), "\n")
cat("   变量数: ", length(cov_vars), "\n")
cat("   [C5] 已删除 pvtt_present（与 pvtt_grade_cat 共线）\n")
cat("   [C4] afp_cat, pvtt_grade_cat, albi_grade_fac 已转为 factor\n")
cat("   [m2] 已删除 pivka_cat（与 pivka_log 冗余）\n")

# ════════════════════════════════════════════════════════════════════
# 3. CBPS-IPTW 权重估计 (ATE)
# ════════════════════════════════════════════════════════════════════
cat("\n3. CBPS-IPTW 权重估计 (ATE)...\n")

W <- weightit(
  iptw_formula,
  data     = analysis_data,
  method   = "cbps",
  estimand = "ATE"
)

cat("\n   权重摘要:\n")
print(summary(W))

# ════════════════════════════════════════════════════════════════════
# 4. [C1] 极端权重截断 (1st/99th percentile)
# ════════════════════════════════════════════════════════════════════
cat("\n4. 极端权重截断 (1st/99th percentile)...\n")

raw_weights <- W$weights
q01 <- quantile(raw_weights, 0.01)
q99 <- quantile(raw_weights, 0.99)
cat(sprintf("   截断前: min=%.2f, max=%.2f, mean=%.2f\n",
            min(raw_weights), max(raw_weights), mean(raw_weights)))
cat(sprintf("   截断阈值: P1=%.2f, P99=%.2f\n", q01, q99))

truncated_weights <- pmin(pmax(raw_weights, q01), q99)
n_truncated <- sum(raw_weights != truncated_weights)
cat(sprintf("   截断后: min=%.2f, max=%.2f, mean=%.2f\n",
            min(truncated_weights), max(truncated_weights), mean(truncated_weights)))
cat(sprintf("   被截断的观测数: %d (%.1f%%)\n", n_truncated, 100 * n_truncated / length(raw_weights)))

analysis_data$iptw_weight <- truncated_weights

# ════════════════════════════════════════════════════════════════════
# 5. 平衡诊断（截断后）
# ════════════════════════════════════════════════════════════════════
cat("\n5. 平衡诊断（截断后权重）...\n")

W_trunc <- W
W_trunc$weights <- truncated_weights

bal <- bal.tab(W_trunc, stats = c("m", "v"), thresholds = c(m = 0.1))
cat("\n")
print(bal)

cat("\n   各组有效样本量 (ESS):\n")
for (g in GROUP_ORDER) {
  sub <- analysis_data %>% filter(group == g)
  w   <- sub$iptw_weight
  ess <- sum(w)^2 / sum(w^2)
  cat(sprintf("     %s: N=%d, ESS=%.0f (%.0f%%)\n", g, nrow(sub), ess, 100 * ess / nrow(sub)))
}

# ════════════════════════════════════════════════════════════════════
# 6. 导出权重
# ════════════════════════════════════════════════════════════════════
cat("\n6. 导出 IPTW 权重...\n")

weight_export <- analysis_data %>%
  select(patient_id, main_group, iptw_weight, os_months, death_status)

write_csv(weight_export, file.path(OUTPUT_DIR, "iptw_weights.csv"))
cat(sprintf("   已保存: iptw_weights.csv (%d 例)\n", nrow(weight_export)))

# ════════════════════════════════════════════════════════════════════
# 7. [C2] 全模型加权 Cox 回归 + [M4] PH 假设检验
# ════════════════════════════════════════════════════════════════════
cat("\n7. 全模型加权 Cox 回归（各组 vs HAIC_alone）...\n")

analysis_data$group <- relevel(analysis_data$group, ref = REF_GROUP_HR)

cox_full <- coxph(
  Surv(os_months, death_status) ~ group,
  data    = analysis_data,
  weights = iptw_weight,
  robust  = TRUE
)

cat("\n   全模型 Cox 回归结果:\n")
print(summary(cox_full))

# 全局 Wald 检验
wald_stat <- cox_full$wald.test
wald_df   <- length(coef(cox_full))
wald_p    <- 1 - pchisq(wald_stat, df = wald_df)
cat(sprintf("\n   全局 Wald 检验: chi2=%.2f, df=%d, P=%s\n",
            wald_stat, wald_df,
            ifelse(wald_p < 0.0001, "< 0.0001", sprintf("%.4f", wald_p))))

# [M4] PH 假设检验
cat("\n   [M4] 比例风险假设检验 (Schoenfeld):\n")
ph_test <- cox.zph(cox_full)
print(ph_test)

ph_global_p <- ph_test$table["GLOBAL", "p"]
if (ph_global_p < 0.05) {
  cat(sprintf("\n   ⚠ PH 假设全局检验 P=%.4f < 0.05，比例风险假设可能不成立\n", ph_global_p))
  cat("   建议: 考虑 RMST 分析或分时段 HR\n")
} else {
  cat(sprintf("\n   ✓ PH 假设全局检验 P=%.4f，比例风险假设成立\n", ph_global_p))
}

# 提取各组 HR
others <- GROUP_ORDER[GROUP_ORDER != REF_GROUP_HR]
cox_coefs <- summary(cox_full)$coefficients
cox_ci    <- exp(confint(cox_full))

survival_results <- data.frame()

for (g in others) {
  coef_name <- paste0("group", g)
  if (!(coef_name %in% rownames(cox_coefs))) next

  hr    <- exp(cox_coefs[coef_name, "coef"])
  se    <- cox_coefs[coef_name, "robust se"]
  z_val <- cox_coefs[coef_name, "coef"] / se
  p_raw <- 2 * pnorm(-abs(z_val))
  ci_lo <- cox_ci[coef_name, 1]
  ci_hi <- cox_ci[coef_name, 2]
  n_grp <- sum(analysis_data$group == g)
  n_ref <- sum(analysis_data$group == REF_GROUP_HR)

  survival_results <- rbind(survival_results, data.frame(
    Group    = g,
    N        = n_grp,
    N_ref    = n_ref,
    HR       = round(hr, 3),
    CI_lower = round(ci_lo, 3),
    CI_upper = round(ci_hi, 3),
    P_value  = p_raw,
    stringsAsFactors = FALSE
  ))

  cat(sprintf("   %s (n=%d) vs %s (n=%d): HR=%.3f (%.3f-%.3f) P=%s\n",
              g, n_grp, REF_GROUP_HR, n_ref,
              hr, ci_lo, ci_hi,
              ifelse(p_raw < 0.0001, "< 0.0001", sprintf("%.4f", p_raw))))
}

# [M3] 多重比较校正
cat("\n   [M3] 多重比较校正:\n")
survival_results$P_holm       <- p.adjust(survival_results$P_value, method = "holm")
survival_results$P_fdr        <- p.adjust(survival_results$P_value, method = "fdr")
survival_results$P_bonferroni <- p.adjust(survival_results$P_value, method = "bonferroni")

for (i in 1:nrow(survival_results)) {
  r <- survival_results[i, ]
  cat(sprintf("     %s: P_raw=%s | P_holm=%s | P_fdr=%s | P_bonf=%s\n",
              r$Group,
              ifelse(r$P_value < 0.0001, "<0.0001", sprintf("%.4f", r$P_value)),
              ifelse(r$P_holm < 0.0001, "<0.0001", sprintf("%.4f", r$P_holm)),
              ifelse(r$P_fdr < 0.0001, "<0.0001", sprintf("%.4f", r$P_fdr)),
              ifelse(r$P_bonferroni < 0.0001, "<0.0001", sprintf("%.4f", r$P_bonferroni))))
}

# [m4] P 值格式修正：不再 round 到 4 位导致 0
survival_results$P_value_fmt <- ifelse(
  survival_results$P_value < 0.0001,
  "< 0.0001",
  sprintf("%.4f", survival_results$P_value)
)

write_csv(survival_results, file.path(OUTPUT_DIR, "survival_gps_final.csv"))
cat(sprintf("\n   已保存: survival_gps_final.csv\n"))

# 保存全局检验结果
global_test <- data.frame(
  Test = c("Wald_global", "PH_global"),
  Statistic = c(wald_stat, ph_test$table["GLOBAL", "chisq"]),
  DF = c(wald_df, ph_test$table["GLOBAL", "df"]),
  P_value = c(wald_p, ph_global_p),
  stringsAsFactors = FALSE
)
write_csv(global_test, file.path(OUTPUT_DIR, "iptw_global_test.csv"))
cat("   已保存: iptw_global_test.csv\n")

# ════════════════════════════════════════════════════════════════════
# 8. [M5] 加权 KM 数据导出 (survey::svykm)
# ════════════════════════════════════════════════════════════════════
cat("\n8. 导出加权 KM 数据 (survey::svykm)...\n")

svy_design <- svydesign(
  ids     = ~1,
  weights = ~iptw_weight,
  data    = analysis_data
)

km_all <- list()
for (g in GROUP_ORDER) {
  sub_design <- subset(svy_design, group == g)
  km_fit <- svykm(Surv(os_months, death_status) ~ 1, design = sub_design, se = TRUE)

  km_df <- data.frame(
    group = g,
    time  = km_fit$time,
    surv  = km_fit$surv,
    se    = km_fit$varlog,
    stringsAsFactors = FALSE
  )
  km_df$surv_lower <- pmax(0, km_df$surv * exp(-1.96 * sqrt(km_df$se)))
  km_df$surv_upper <- pmin(1, km_df$surv * exp( 1.96 * sqrt(km_df$se)))

  km_all[[g]] <- km_df
}

km_export <- do.call(rbind, km_all)
write_csv(km_export, file.path(OUTPUT_DIR, "iptw_km_data.csv"))
cat(sprintf("   已保存: iptw_km_data.csv (%d 行)\n", nrow(km_export)))

cat("\n============================================================\n")
cat("CBPS-IPTW 多组加权完成！\n")
cat(sprintf("  方法: CBPS (ATE) + 权重截断 (P1/P99)\n"))
cat(sprintf("  全模型 Cox: group ~ HAIC_alone (ref), robust SE\n"))
cat(sprintf("  总样本: %d\n", nrow(analysis_data)))
cat(sprintf("结果: %s\n", OUTPUT_DIR))
cat("============================================================\n")

sink()
