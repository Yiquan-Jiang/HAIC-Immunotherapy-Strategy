#!/usr/bin/env Rscript
# ============================================================
# 高危亚组分析 — Overlap Weighting + 加权 Cox
# 对比1: HAIC+I_concurrent vs HAIC_then_I
# 对比2: HAIC+I+T_concurrent vs HAIC_then_I+T
# ============================================================
#
# 方法: Overlap Weighting (Li, Morgan & Zaslavsky, JASA 2018)
#   - PS 模型: logistic regression (~12 核心变量)
#   - 权重: w_treat = 1 - PS, w_control = PS
#   - 效应估计: 加权 Cox (robust sandwich SE)
#   - 平衡性: 加权后 SMD < 0.10
#
# 输出 (每组对比加后缀 _I / _IT):
#   results/subgroup_analysis/{tag}/ow_subgroup_results_{tag}.csv
#   results/subgroup_analysis/{tag}/ow_smd_balance_{tag}.csv
#   results/subgroup_analysis/{tag}/ow_weighted_ids_*.csv
#   figures/subgroup/{tag}/ow_love_plot_*.pdf/png

library(tidyverse)
library(WeightIt)
library(cobalt)
library(survival)
library(survey)

BASE_DIR   <- "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
DATA_DIR   <- file.path(BASE_DIR, "data")
RES_DIR    <- file.path(BASE_DIR, "results", "subgroup_analysis")
FIG_DIR    <- file.path(BASE_DIR, "figures", "subgroup")
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("  高危亚组 Overlap Weighting 分析\n")
cat("  对比1: HAIC+I_concurrent vs HAIC_then_I\n")
cat("  对比2: HAIC+I+T_concurrent vs HAIC_then_I+T\n")
cat("============================================================\n")

# ═══════���════════════════════════════════════════════════════════
# 0. ��取全部数据 & 定义对比组
# ══��═════════════���═══════════════════════════════════════════════
cat("\n0. 读取数据...\n")

df_all_raw <- read_csv(file.path(DATA_DIR, "analysis_ready.csv"), show_col_types = FALSE) %>%
  filter(os_months >= 0)

COMPARISONS <- list(
  list(tag = "I",
       treat_group   = "HAIC+I_concurrent",
       control_group = "HAIC_then_I",
       title_short   = "HAIC+Immuno Concurrent vs Delayed"),
  list(tag = "IT",
       treat_group   = "HAIC+I+T_concurrent",
       control_group = "HAIC_then_I+T",
       title_short   = "HAIC+Immuno+TKI Concurrent vs Delayed")
)

prepare_two_group <- function(df_raw, treat_group, control_group) {
  df_raw %>%
    filter(main_group %in% c(treat_group, control_group)) %>%
    mutate(
      treat = if_else(main_group == treat_group, 1L, 0L),
      event = case_when(
        death_status %in% c("Yes", "1", "TRUE") ~ 1L,
        TRUE ~ 0L
      ),
      sex_male         = if_else(sex == "Male", 1L, 0L),
      pvtt_vp34        = if_else(pvtt_classification == "Vp3/4", 1L, 0L),
      pvtt_vp12        = if_else(pvtt_classification == "Vp1/2", 1L, 0L),
      pvtt_present     = if_else(pvtt_classification != "Absent", 1L, 0L),
      hvtt_yes         = if_else(hvtt == "Yes", 1L, 0L),
      ivc_ra_yes       = if_else(ivc_or_ra_thrombus == "Yes", 1L, 0L),
      dist_meta_yes    = if_else(distant_metastasis == "Yes", 1L, 0L),
      lymph_meta_yes   = if_else(lymph_node_metastasis == "Yes", 1L, 0L),
      ascites_yes      = if_else(ascites != "Absent", 1L, 0L),
      varices_yes      = if_else(varices == "Yes", 1L, 0L),
      tumor_gt10       = if_else(tumor_max_diameter_cm > 10, 1L, 0L),
      tumor_multi      = if_else(tumor_count_category == ">3", 1L, 0L),
      afp_high_bin     = if_else(afp_high == "Yes", 1L, 0L),
      pivka_high_bin   = if_else(pivka_high == "Yes", 1L, 0L),
      albi_grade_num   = as.integer(albi_grade),
      high_risk_composite = as.integer(
        ivc_ra_yes | tumor_multi | pvtt_vp34 | dist_meta_yes | tumor_gt10
      )
    )
}

# ════════════════════════════════════════════════════════════════
# 遍历每组对比
# ══���════════════════════════��════════════════════════════════════
for (COMP in COMPARISONS) {

CMP_TAG       <- COMP$tag
TREAT_GROUP   <- COMP$treat_group
CONTROL_GROUP <- COMP$control_group
TITLE_SHORT   <- COMP$title_short

CMP_RES_DIR <- file.path(RES_DIR, CMP_TAG)
CMP_FIG_DIR <- file.path(FIG_DIR, CMP_TAG)
dir.create(CMP_RES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CMP_FIG_DIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("\n######################################################\n"))
cat(sprintf("#  对比组: %s vs %s  (tag=%s)\n", TREAT_GROUP, CONTROL_GROUP, CMP_TAG))
cat(sprintf("######################################################\n"))

df <- prepare_two_group(df_all_raw, TREAT_GROUP, CONTROL_GROUP)

cat(sprintf("   总样本: %d (Concurrent=%d, Delayed=%d)\n",
            nrow(df), sum(df$treat == 1), sum(df$treat == 0)))
cat(sprintf("   事件数: %d (%.1f%%)\n", sum(df$event), mean(df$event) * 100))

# ════════════════════════════════════════════════════════════════
# 2. 定义亚组和 PS 变量
# ════════════════════════════════════════════════════════════════

CORE_PS_VARS <- c(
  "age", "sex_male", "albi_score",
  "tumor_max_diameter_cm", "tumor_multi",
  "pvtt_vp34", "hvtt_yes",
  "dist_meta_yes", "lymph_meta_yes",
  "ascites_yes",
  "afp_high_bin", "pivka_high_bin",
  "nlr", "plt", "alb"
)

SUBGROUPS <- list(
  list(name = "Composite high-risk",
       filter_expr = quote(high_risk_composite == 1),
       exclude_ps = character(0)),
  list(name = "Tumor count >3",
       filter_expr = quote(tumor_multi == 1),
       exclude_ps = c("tumor_multi")),
  list(name = "Tumor diameter >10 cm",
       filter_expr = quote(tumor_gt10 == 1),
       exclude_ps = c("tumor_gt10", "tumor_max_diameter_cm")),
  list(name = "PVTT Vp3/4",
       filter_expr = quote(pvtt_vp34 == 1),
       exclude_ps = c("pvtt_vp34")),
  list(name = "Extrahepatic metastasis",
       filter_expr = quote(dist_meta_yes == 1),
       exclude_ps = c("dist_meta_yes"))
)

# 互补亚组
SUBGROUPS_COMP <- list(
  list(name = "Non-Tumor count >3",
       filter_expr = quote(tumor_multi == 0),
       exclude_ps = c("tumor_multi")),
  list(name = "Non-Tumor diameter >10 cm",
       filter_expr = quote(tumor_gt10 == 0),
       exclude_ps = c("tumor_gt10", "tumor_max_diameter_cm")),
  list(name = "Non-PVTT Vp3/4",
       filter_expr = quote(pvtt_vp34 == 0),
       exclude_ps = c("pvtt_vp34")),
  list(name = "Non-Extrahepatic metastasis",
       filter_expr = quote(dist_meta_yes == 0),
       exclude_ps = c("dist_meta_yes"))
)

ALL_SUBGROUPS <- c(SUBGROUPS, SUBGROUPS_COMP)

# ════════════════════════════════════════════════════════════════
# 3. 交互作用检验 (全人群)
# ════════════════════════════════════════════════════════════════
cat("\n2. 交互作用检验 (全人群 n=%d)...\n", nrow(df))

interaction_vars <- list(
  list(name = "Tumor count >3",       var = "tumor_multi"),
  list(name = "Tumor diameter >10 cm", var = "tumor_gt10"),
  list(name = "PVTT Vp3/4",           var = "pvtt_vp34"),
  list(name = "Extrahepatic metastasis", var = "dist_meta_yes")
)

inter_results <- data.frame()
for (iv in interaction_vars) {
  fml <- as.formula(paste0(
    "Surv(os_months, event) ~ treat * ", iv$var,
    " + age + albi_score + afp_high_bin + nlr"
  ))
  fit <- coxph(fml, data = df)
  inter_term <- paste0("treat:", iv$var)
  coefs <- summary(fit)$coefficients
  if (inter_term %in% rownames(coefs)) {
    p_inter <- coefs[inter_term, "Pr(>|z|)"]
  } else {
    p_inter <- NA
  }
  inter_results <- rbind(inter_results, data.frame(
    Subgroup = iv$name,
    P_interaction = round(p_inter, 4)
  ))
  cat(sprintf("   %s: P_interaction = %.4f\n", iv$name, p_inter))
}

# Composite
fit_comp <- coxph(Surv(os_months, event) ~ treat * high_risk_composite +
                    age + albi_score + afp_high_bin + nlr, data = df)
coefs_comp <- summary(fit_comp)$coefficients
p_comp <- coefs_comp["treat:high_risk_composite", "Pr(>|z|)"]
inter_results <- rbind(inter_results, data.frame(
  Subgroup = "Composite high-risk", P_interaction = round(p_comp, 4)
))
cat(sprintf("   Composite high-risk: P_interaction = %.4f\n", p_comp))

write_csv(inter_results, file.path(CMP_RES_DIR, paste0("ow_interaction_tests_", CMP_TAG, ".csv")))

# ════════════════════════════════════════════════════════════════
# 4. OW 分析函数
# ════════════════════════════════════════════════════════════════

analyze_ow <- function(data, sg_name, exclude_ps, save_love = TRUE) {

  n_treat <- sum(data$treat == 1)
  n_ctrl  <- sum(data$treat == 0)
  n_events <- sum(data$event)

  cat(sprintf("\n   --- %s (n=%d+%d, events=%d) ---\n",
              sg_name, n_treat, n_ctrl, n_events))

  if (n_treat < 10 || n_ctrl < 10 || n_events < 15) {
    cat("   ⚠ 样本量不足，跳过\n")
    return(NULL)
  }

  # PS 变量
  ps_vars <- setdiff(CORE_PS_VARS, exclude_ps)

  # 检查常量列
  ps_vars <- ps_vars[sapply(ps_vars, function(v) sd(data[[v]], na.rm = TRUE) > 0)]

  ps_formula <- as.formula(paste("treat ~", paste(ps_vars, collapse = " + ")))

  # Overlap Weighting
  W <- tryCatch(
    weightit(ps_formula, data = data, method = "glm", estimand = "ATO"),
    error = function(e) {
      cat(sprintf("   ⚠ WeightIt 失败: %s\n", conditionMessage(e)))
      NULL
    }
  )
  if (is.null(W)) return(NULL)

  data$ow_weights <- W$weights

  # 平衡性检查
  bal <- bal.tab(W, stats = "mean.diffs", thresholds = c(m = 0.1), un = TRUE)
  bal_df <- bal$Balance
  bal_df$Variable <- rownames(bal_df)

  n_unbal_before <- sum(abs(bal_df$Diff.Un) > 0.1, na.rm = TRUE)
  n_unbal_after  <- sum(abs(bal_df$Diff.Adj) > 0.1, na.rm = TRUE)
  max_smd_before <- max(abs(bal_df$Diff.Un), na.rm = TRUE)
  max_smd_after  <- max(abs(bal_df$Diff.Adj), na.rm = TRUE)

  cat(sprintf("   OW 前: %d/%d vars SMD>0.1, max=%.3f\n",
              n_unbal_before, nrow(bal_df), max_smd_before))
  cat(sprintf("   OW 后: %d/%d vars SMD>0.1, max=%.3f\n",
              n_unbal_after, nrow(bal_df), max_smd_after))

  # ESS
  ess_treat <- sum(data$ow_weights[data$treat == 1])^2 /
    sum(data$ow_weights[data$treat == 1]^2)
  ess_ctrl  <- sum(data$ow_weights[data$treat == 0])^2 /
    sum(data$ow_weights[data$treat == 0]^2)
  cat(sprintf("   ESS: treat=%.0f, ctrl=%.0f\n", ess_treat, ess_ctrl))

  # Love Plot
  if (save_love) {
    lp <- love.plot(W,
                    stats = "mean.diffs",
                    threshold = 0.1,
                    abs = TRUE,
                    var.order = "unadjusted",
                    colors = c("#999999", "#0072B2"),
                    shapes = c("circle", "triangle"),
                    size = 3,
                    title = sg_name,
                    sample.names = c("Before OW", "After OW"),
                    limits = c(0, 0.6)) +
      ggplot2::theme_bw(base_size = 9, base_family = "Helvetica") +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        legend.position    = "bottom",
        legend.title       = ggplot2::element_blank(),
        plot.title         = ggplot2::element_text(size = 10, face = "bold"),
        axis.text          = ggplot2::element_text(size = 8),
        axis.title         = ggplot2::element_text(size = 9)
      ) +
      ggplot2::labs(x = "Absolute Standardized Mean Difference")

    safe_name <- gsub("[+>≥ ]", "_", gsub("[/]", "_", sg_name))
    lp_path <- file.path(CMP_FIG_DIR, paste0("ow_love_", safe_name))
    ggplot2::ggsave(paste0(lp_path, ".png"), lp, width = 5.5, height = 5, dpi = 300)
    tryCatch(
      ggplot2::ggsave(paste0(lp_path, ".pdf"), lp, width = 5.5, height = 5,
                       device = cairo_pdf),
      error = function(e) cat("   PDF保存失败，已保存PNG\n")
    )
    cat(sprintf("   Love Plot: %s\n", basename(lp_path)))
  }

  # 加权 Cox (robust SE)
  des <- svydesign(ids = ~1, weights = ~ow_weights, data = data)
  cox_fit <- svycoxph(Surv(os_months, event) ~ treat, design = des)
  cox_sum <- summary(cox_fit)

  hr     <- exp(coef(cox_fit)["treat"])
  ci     <- exp(confint(cox_fit)["treat", ])
  ci_lo  <- ci[1]
  ci_hi  <- ci[2]
  p_val  <- cox_sum$coefficients["treat", "Pr(>|z|)"]

  # 90% CI (等价性)
  z90    <- qnorm(0.95)
  se_log <- cox_sum$coefficients["treat", "robust se"]
  ci90_lo <- exp(coef(cox_fit)["treat"] - z90 * se_log)
  ci90_hi <- exp(coef(cox_fit)["treat"] + z90 * se_log)

  equiv <- (ci90_lo >= 0.60) & (ci90_hi <= 1.67)

  # Log-rank (unweighted, for reference)
  lr <- survdiff(Surv(os_months, event) ~ treat, data = data)
  lr_p <- 1 - pchisq(lr$chisq, df = 1)

  # 未调整 HR
  cox_unadj <- coxph(Surv(os_months, event) ~ treat, data = data)
  hr_unadj  <- exp(coef(cox_unadj)["treat"])
  ci_unadj  <- exp(confint(cox_unadj)["treat", ])
  p_unadj   <- summary(cox_unadj)$coefficients["treat", "Pr(>|z|)"]

  # 中位 OS
  sf <- survfit(Surv(os_months, event) ~ treat, data = data)
  med_os <- summary(sf)$table[, "median"]

  cat(sprintf("   未调整: HR=%.2f (%.2f-%.2f) P=%.4f\n",
              hr_unadj, ci_unadj[1], ci_unadj[2], p_unadj))
  cat(sprintf("   OW加权: HR=%.2f (%.2f-%.2f) P=%.4f\n",
              hr, ci_lo, ci_hi, p_val))
  cat(sprintf("   90%%CI: (%.2f-%.2f) → 等价性: %s\n",
              ci90_lo, ci90_hi, ifelse(equiv, "✓", "✗")))

  # 导出加权患者数据 (供 Python KM)
  export <- data %>%
    select(patient_id, treat, os_months, event, ow_weights) %>%
    mutate(subgroup = sg_name)
  safe_fn <- gsub("[+>≥ /]", "_", sg_name)
  write_csv(export, file.path(CMP_RES_DIR, paste0("ow_weighted_ids_", safe_fn, ".csv")))

  # SMD 保存
  smd_export <- bal_df %>%
    mutate(Subgroup = sg_name) %>%
    select(Subgroup, Variable, SMD_before = Diff.Un, SMD_after = Diff.Adj)

  list(
    result = data.frame(
      Subgroup         = sg_name,
      N_concurrent     = n_treat,
      N_then_I         = n_ctrl,
      Events_total     = n_events,
      Median_OS_conc   = round(med_os[2], 1),
      Median_OS_then   = round(med_os[1], 1),
      HR_unadjusted    = round(hr_unadj, 3),
      CI95_lo_unadj    = round(ci_unadj[1], 3),
      CI95_hi_unadj    = round(ci_unadj[2], 3),
      P_unadjusted     = round(p_unadj, 4),
      HR_OW            = round(hr, 3),
      CI95_lo_OW       = round(ci_lo, 3),
      CI95_hi_OW       = round(ci_hi, 3),
      P_OW             = round(p_val, 4),
      CI90_lo_OW       = round(ci90_lo, 3),
      CI90_hi_OW       = round(ci90_hi, 3),
      Equivalence      = equiv,
      N_PS_vars        = length(ps_vars),
      N_unbal_before   = n_unbal_before,
      N_unbal_after    = n_unbal_after,
      Max_SMD_before   = round(max_smd_before, 3),
      Max_SMD_after    = round(max_smd_after, 3),
      ESS_treat        = round(ess_treat, 0),
      ESS_ctrl         = round(ess_ctrl, 0),
      Logrank_P        = round(lr_p, 4),
      stringsAsFactors = FALSE
    ),
    smd = smd_export
  )
}

# ════════════════════════════════════════════════════════════════
# 5. 执行分析
# ════════════════════════════════════════════════════════════════
cat("\n3. 执行 Overlap Weighting 分析...\n")

all_results <- data.frame()
all_smd     <- data.frame()

for (sg in ALL_SUBGROUPS) {
  sg_data <- df %>% filter(!!sg$filter_expr)
  res <- analyze_ow(sg_data, sg$name, sg$exclude_ps)
  if (!is.null(res)) {
    all_results <- rbind(all_results, res$result)
    all_smd     <- rbind(all_smd, res$smd)
  }
}

# 保存
write_csv(all_results, file.path(CMP_RES_DIR, paste0("ow_subgroup_results_", CMP_TAG, ".csv")))
write_csv(all_smd, file.path(CMP_RES_DIR, paste0("ow_smd_balance_", CMP_TAG, ".csv")))

# ════════════════════════════════════════════════════════════════
# 6. 汇总
# ════════════════════════════════════════════════════════════════
cat("\n============================================================\n")
cat("  分析完成！结果汇总\n")
cat("============================================================\n\n")

cat(sprintf("%-35s %22s %8s %18s %6s %5s %5s\n",
            "Subgroup", "OW HR (95%CI)", "P", "90%CI", "Equiv",
            "SMD>0.1", "ESS"))
cat(paste(rep("-", 105), collapse = ""), "\n")

for (i in seq_len(nrow(all_results))) {
  r <- all_results[i, ]
  p_str <- ifelse(r$P_OW < 0.001, "<0.001", sprintf("%.3f", r$P_OW))
  eq_str <- ifelse(r$Equivalence, "✓", "✗")
  cat(sprintf("%-35s %5.2f (%4.2f-%4.2f) %8s (%4.2f-%4.2f) %5s %3d→%d %3.0f+%.0f\n",
              r$Subgroup,
              r$HR_OW, r$CI95_lo_OW, r$CI95_hi_OW,
              p_str,
              r$CI90_lo_OW, r$CI90_hi_OW,
              eq_str,
              r$N_unbal_before, r$N_unbal_after,
              r$ESS_treat, r$ESS_ctrl))
}

cat(sprintf("\n交互作用检验:\n"))
for (i in seq_len(nrow(inter_results))) {
  cat(sprintf("  %-35s P_interaction = %.4f\n",
              inter_results$Subgroup[i], inter_results$P_interaction[i]))
}

cat(sprintf("\n输出文件:\n"))
cat(sprintf("  %s\n", file.path(CMP_RES_DIR, paste0("ow_subgroup_results_", CMP_TAG, ".csv"))))
cat(sprintf("  %s\n", file.path(CMP_RES_DIR, paste0("ow_smd_balance_", CMP_TAG, ".csv"))))
cat(sprintf("  %s\n", file.path(CMP_RES_DIR, paste0("ow_interaction_tests_", CMP_TAG, ".csv"))))
cat(sprintf("  %s (Love Plots)\n", CMP_FIG_DIR))

}  # end for (COMP in COMPARISONS)
