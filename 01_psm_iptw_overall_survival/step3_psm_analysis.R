#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
#
# PSM 分析 — update_group_7（7组，21组两两对比）
# 数据源: update_group_7/data/analysis_ready.csv
# 输出:
#   results/psm_balance_tables_complete/matched_ids_XX_*.csv
#   results/psm_balance_tables_complete/survival_analysis_final.csv
#   figures/psm_final/survival_curves_final.pdf

library(tidyverse)
library(MatchIt)
library(survival)
library(survminer)

BASE_DIR   <- "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
EIGHT_GROUP <- Sys.getenv("EIGHT_GROUP", "0") == "1"
SFX <- if (EIGHT_GROUP) "_8group" else ""
DATA_CSV <- if (EIGHT_GROUP) "analysis_ready_8group.csv" else "analysis_ready.csv"
DATA_DIR   <- file.path(BASE_DIR, "data")
OUTPUT_DIR <- file.path(BASE_DIR, "results", paste0("psm_balance_tables_complete", SFX))
FIGURE_DIR <- file.path(BASE_DIR, "figures", paste0("psm_final", SFX))
LOG_DIR    <- file.path(BASE_DIR, "logs")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGURE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR,    showWarnings = FALSE, recursive = TRUE)

# 日志
log_file <- file.path(LOG_DIR, "psm_analysis.log")
sink(log_file, split = TRUE)

cat("============================================================\n")
cat("PSM 分析 — update_group_7（7组，21组两两对比）\n")
cat("============================================================\n")

# ── 配色（Okabe-Ito 扩展，7色）──────────────────────────────────
group_colors <- c(
  "HAIC_alone"            = "#0072B2",
  "HAIC+I_concurrent"     = "#E69F00",
  "HAIC_then_I"           = "#009E73",
  "HAIC+T_concurrent"     = "#F0E442",
  "HAIC_then_T"           = "#CC79A7",
  "HAIC+I+T_concurrent"   = "#D55E00",
  "HAIC_then_I+T"         = "#56B4E9",
  "Systemic_I+T"          = "#117733"
)
GROUP_ORDER <- names(group_colors)
if (!EIGHT_GROUP) GROUP_ORDER <- setdiff(GROUP_ORDER, "Systemic_I+T")

# ════════════════════════════════════════════════════════════════════
# 1. 读取数据
# ════════════════════════════════════════════════════════════════════
cat("\n1. 读取数据...\n")

analysis_data <- read_csv(
  file.path(DATA_DIR, DATA_CSV),
  show_col_types = FALSE
) %>%
  filter(os_months >= 0) %>%
  mutate(
    group = factor(main_group, levels = GROUP_ORDER),
    sex_male = if_else(sex == "Male", 1L, 0L, missing = 1L),
    death_status = case_when(
      death_status %in% c("Yes", "1", "TRUE", "yes") ~ 1L,
      death_status %in% c("No",  "0", "FALSE","no")  ~ 0L,
      TRUE ~ as.integer(as.numeric(death_status))
    )
  ) %>%
  filter(!is.na(group))

cat(sprintf("   - 有效分析患者: %d\n", nrow(analysis_data)))
cat("   - 分组分布:\n")
print(table(analysis_data$group))

# ════════════════════════════════════════════════════════════════════
# 2. 构建 PSM 变量
# ════════════════════════════════════════════════════════════════════
cat("\n2. 构建 PSM 匹配变量...\n")

analysis_data <- analysis_data %>%
  mutate(
    afp_cat = case_when(
      afp < 20   ~ 0L,
      afp < 400  ~ 1L,
      TRUE       ~ 2L
    ),
    pivka_cat = case_when(
      pivka < 40   ~ 0L,
      pivka < 400  ~ 1L,
      TRUE         ~ 2L
    ),
    tumor_gt10cm = if_else(tumor_max_diameter_cm > 10, 1L, 0L),
    tumor_multiple = if_else(tumor_count_category == ">3", 1L, 0L),
    pvtt_grade_cat = case_when(
      pvtt_classification == "Absent" ~ 0L,
      pvtt_classification == "Vp1/2"  ~ 1L,
      pvtt_classification == "Vp3/4"  ~ 2L,
      TRUE                            ~ 0L
    ),
    pvtt_present        = if_else(pvtt_classification != "Absent", 1L, 0L),
    hvtt_present        = if_else(hvtt == "Yes", 1L, 0L),
    ivc_ra_present      = if_else(ivc_or_ra_thrombus == "Yes", 1L, 0L),
    distant_meta_bin    = if_else(distant_metastasis == "Yes", 1L, 0L),
    lymph_meta_bin      = if_else(lymph_node_metastasis == "Yes", 1L, 0L),
    ascites_bin         = if_else(ascites != "Absent", 1L, 0L),
    varices_bin         = if_else(varices == "Yes", 1L, 0L),
    albi_grade_num      = as.integer(albi_grade),
    afp_log   = log10(pmax(afp,   0.01) + 1),
    pivka_log = log10(pmax(pivka, 0.01) + 1),
    tbil_log  = log10(pmax(tbil,  0.01) + 1),
    age_std        = as.numeric(scale(age)),
    pivka_std      = as.numeric(scale(pivka_log)),
    tumor_size_std = as.numeric(scale(tumor_max_diameter_cm)),
    tbil_std       = as.numeric(scale(tbil_log)),
    alb_std        = as.numeric(scale(alb)),
    plt_std        = as.numeric(scale(plt)),
    nlr_std        = as.numeric(scale(nlr))
  )

# ════════════════════════════════════════════════════════════════════
# 3. PSM 匹配（21组两两对比）
# ════════════════════════════════════════════════════════════════════
cat("\n3. PSM 匹配（21组两两对比）...\n")

PSM_FORMULA <- treatment ~
  afp_cat + pivka_cat + pivka_std +
  tumor_gt10cm + tumor_multiple +
  pvtt_grade_cat + pvtt_present + hvtt_present +
  ivc_ra_present + distant_meta_bin + lymph_meta_bin +
  ascites_bin + varices_bin +
  albi_grade_num + tbil_std + alb_std + plt_std +
  age_std + tumor_size_std + nlr_std

# 生成所有21组两两对比
all_groups <- GROUP_ORDER
n_groups   <- length(all_groups)
comparisons <- list()
comp_idx <- 1
for (i in 1:(n_groups - 1)) {
  for (j in (i+1):n_groups) {
    g1 <- all_groups[i]
    g2 <- all_groups[j]
    key <- paste0(g1, "_vs_", g2)
    comparisons[[comp_idx]] <- list(
      id     = comp_idx,
      group1 = g1,
      group2 = g2,
      key    = key
    )
    comp_idx <- comp_idx + 1
  }
}

cat(sprintf("   共 %d 组两两对比\n", length(comparisons)))

survival_results  <- data.frame()
matched_data_list <- list()

for (comp in comparisons) {
  i    <- comp$id
  g1   <- comp$group1
  g2   <- comp$group2
  ckey <- comp$key

  cat(sprintf("\n   === 比较 %02d: %s vs %s ===\n", i, g1, g2))

  comp_data <- analysis_data %>%
    filter(group %in% c(g1, g2)) %>%
    mutate(
      treatment   = if_else(group == g2, 1L, 0L),
      group_label = factor(group, levels = c(g1, g2))
    )

  n_treat <- sum(comp_data$treatment == 1)
  n_ctrl  <- sum(comp_data$treatment == 0)
  cat(sprintf("   PSM 前: %s=%d, %s=%d\n", g1, n_ctrl, g2, n_treat))

  # 自适应卡钳值：小样本组用更宽的卡钳
  min_n <- min(n_treat, n_ctrl)
  caliper_val <- if (min_n < 150) 0.25 else if (min_n < 300) 0.15 else 0.10

  # 固定每个对比的随机种子，确保 matchit 对并列距离的 tie-breaking 完全可重复
  set.seed(1000L + i)

  match_result <- tryCatch(
    matchit(
      PSM_FORMULA,
      data        = comp_data,
      method      = "nearest",
      distance    = "glm",
      link        = "logit",
      caliper     = caliper_val,
      std.caliper = TRUE,
      ratio       = 1,
      replace     = FALSE
    ),
    error = function(e) {
      cat(sprintf("   ⚠ matchit 失败: %s\n", conditionMessage(e)))
      NULL
    }
  )

  if (is.null(match_result)) next

  matched_data <- match.data(match_result)
  n1 <- sum(matched_data$treatment == 0)
  n2 <- sum(matched_data$treatment == 1)
  cat(sprintf("   PSM 后: %s=%d, %s=%d (caliper=%.2f)\n", g1, n1, g2, n2, caliper_val))

  matched_data_list[[ckey]] <- matched_data

  # 导出匹配后患者 ID（供 Python KM 脚本使用）
  id_export <- matched_data %>%
    select(patient_id, group_label, treatment, subclass) %>%
    mutate(comparison_key = ckey,
           group1 = g1,
           group2 = g2)
  write_csv(id_export,
            file.path(OUTPUT_DIR, sprintf("matched_ids_%02d_%s.csv", i, ckey)))

  # 生存分析
  # HR 方向约定：HR = h(Group1) / h(Group2)，即 "A_vs_B" 表示 A 相对于 B
  # 故 Cox 模型使用 g2 作为参考组，coef 即 g1 相对 g2 的 log-HR
  surv_fit     <- survfit(Surv(os_months, death_status) ~ group_label, data = matched_data)
  logrank_test <- survdiff(Surv(os_months, death_status) ~ group_label, data = matched_data)
  logrank_p    <- 1 - pchisq(logrank_test$chisq, df = 1)
  median_os    <- summary(surv_fit)$table[, "median"]
  cox_model    <- coxph(Surv(os_months, death_status) ~ relevel(group_label, ref = g2),
                        data = matched_data)
  hr           <- exp(coef(cox_model))
  hr_ci        <- exp(confint(cox_model))

  survival_results <- rbind(survival_results, data.frame(
    Comparison   = sprintf("%02d_%s", i, ckey),
    Group1       = g1,
    Group2       = g2,
    N1_before    = n_ctrl,
    N2_before    = n_treat,
    N1_after     = n1,
    N2_after     = n2,
    Caliper      = caliper_val,
    Median_OS_1  = round(median_os[1], 1),
    Median_OS_2  = round(median_os[2], 1),
    HR           = round(hr, 3),
    CI_lower     = round(hr_ci[1], 3),
    CI_upper     = round(hr_ci[2], 3),
    P_value      = round(logrank_p, 4)
  ))

  cat(sprintf("   中位OS: %.1f vs %.1f mo | HR=%.2f (%.2f-%.2f) | P=%.4f\n",
              median_os[1], median_os[2], hr, hr_ci[1], hr_ci[2], logrank_p))
}

# ════════════════════════════════════════════════════════════════════
# 4. 保存生存分析结果
# ════════════════════════════════════════════════════════════════════
cat("\n4. 保存结果...\n")
write_csv(survival_results, file.path(OUTPUT_DIR, "survival_analysis_final.csv"))
cat(sprintf("   已保存: %s\n", file.path(OUTPUT_DIR, "survival_analysis_final.csv")))
cat(sprintf("   共 %d 组对比完成\n", nrow(survival_results)))

# ════════════════════════════════════════════════════════════════════
# 5. 生成 KM 生存曲线（PDF）
# ════════════════════════════════════════════════════════════════════
cat("\n5. 生成 KM 曲线...\n")

theme_pub_km <- theme_classic(base_size = 9, base_family = "Helvetica") +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    panel.border      = element_blank(),
    axis.line         = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.4, colour = "black"),
    axis.text         = element_text(size = 7,  colour = "#333333"),
    axis.title        = element_text(size = 9,  colour = "#333333"),
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.text       = element_text(size = 6),
    legend.title      = element_text(size = 7),
    legend.position   = "right",
    plot.title        = element_text(size = 10, face = "bold", hjust = 0, colour = "#333333"),
    strip.background  = element_blank(),
    strip.text        = element_text(size = 8,  face = "bold", colour = "#333333")
  )

pdf_file <- file.path(FIGURE_DIR, "survival_curves_final.pdf")
pdf(pdf_file, width = 7.5, height = 5.5)

# PSM 前整体曲线（7组）
surv_fit_all <- survfit(Surv(os_months, death_status) ~ group, data = analysis_data)
p_overall <- ggsurvplot(
  surv_fit_all, data = analysis_data,
  pval = TRUE, pval.method = TRUE, pval.size = 3.5,
  risk.table = TRUE, conf.int = TRUE, conf.int.alpha = 0.08,
  palette          = unname(group_colors),
  legend.labs      = levels(analysis_data$group),
  legend.title     = "Treatment Group",
  xlab             = "Time (months)",
  ylab             = "Overall Survival Probability",
  title            = sprintf("Overall Survival — Before PSM (N=%d)", nrow(analysis_data)),
  risk.table.height = 0.32,
  risk.table.fontsize = 2.5,
  risk.table.y.text   = FALSE,
  ggtheme          = theme_pub_km,
  break.time.by    = 12,
  xlim             = c(0, 60),
  surv.median.line = "hv",
  size             = 0.8
)
print(p_overall)

# PSM 后各对比
for (comp in comparisons) {
  i    <- comp$id
  g1   <- comp$group1
  g2   <- comp$group2
  ckey <- comp$key

  if (is.null(matched_data_list[[ckey]])) next

  matched_data <- matched_data_list[[ckey]]
  n1 <- sum(matched_data$treatment == 0)
  n2 <- sum(matched_data$treatment == 1)

  surv_fit <- survfit(Surv(os_months, death_status) ~ group_label, data = matched_data)

  p <- ggsurvplot(
    surv_fit, data = matched_data,
    pval = TRUE, conf.int = TRUE, conf.int.alpha = 0.08,
    risk.table = TRUE,
    palette          = c(group_colors[g1], group_colors[g2]),
    legend.labs      = c(g1, g2),
    legend.title     = "Treatment Group",
    xlab             = "Time (months)",
    ylab             = "Overall Survival Probability",
    title            = sprintf("Comp %02d: %s vs %s\n(After PSM, n=%d vs %d)",
                               i, g1, g2, n1, n2),
    risk.table.height   = 0.25,
    risk.table.fontsize = 3.0,
    risk.table.y.text   = FALSE,
    ggtheme          = theme_pub_km,
    break.time.by    = 12,
    xlim             = c(0, 60),
    surv.median.line = "hv",
    linetype         = c("solid", "dashed"),
    size             = 1.0
  )
  print(p)
}

dev.off()
cat(sprintf("   已保存: %s\n", pdf_file))

cat("\n============================================================\n")
cat("PSM 分析完成！\n")
cat(sprintf("结果: %s\n", OUTPUT_DIR))
cat(sprintf("图表: %s\n", FIGURE_DIR))
cat("============================================================\n")

sink()
