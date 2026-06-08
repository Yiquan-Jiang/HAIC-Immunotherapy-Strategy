#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
#
# 生成论文所需表格 + PSM Love Plot — update_group_7
# Table S0（全队列描述性总表）+ Table 1（7组对比）+ 21组 PSM 平衡表 + 21组 Love Plot
#
# 输出：
#   results/tables/table_s0_overall_descriptive.docx  ← 新增：全队列描述性总表（竖向）
#   results/tables/table1_overall_baseline.docx        ← 7组基线对比（横向）
#   results/tables/tableXX_compXX_*_psm_balance.docx  ← 21张平衡表（竖向）
#   figures/loveplots/XX_loveplot_compXX_*.pdf/png     ← 21组 Love Plot
# 变更：剔除 APTT / Fibrinogen / HBsAg / PTA；所有表格含缩写备注脚注

library(tidyverse)
library(MatchIt)
library(cobalt)
library(gtsummary)
library(flextable)
library(officer)
library(survival)

BASE_DIR   <- "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
EIGHT_GROUP <- Sys.getenv("EIGHT_GROUP", "0") == "1"
SFX <- if (EIGHT_GROUP) "_8group" else ""
DATA_CSV <- if (EIGHT_GROUP) "analysis_ready_8group.csv" else "analysis_ready.csv"
DATA_DIR   <- file.path(BASE_DIR, "data")
TABLE_DIR  <- file.path(BASE_DIR, "results", paste0("tables", SFX))
LOVE_DIR   <- file.path(BASE_DIR, "figures", paste0("loveplots", SFX))
LOG_DIR    <- file.path(BASE_DIR, "logs")
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(LOVE_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(LOG_DIR,   showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(LOG_DIR, "tables_loveplots.log")
sink(log_file, split = TRUE)

cat("============================================================\n")
cat("生成论文表格 + Love Plot — update_group_7\n")
cat("============================================================\n")

group_colors <- c(
  "HAIC_alone"            = "#0072B2",
  "HAIC+I_concurrent"     = "#E69F00",
  "HAIC_then_I"           = "#009E73",
  "HAIC+T_concurrent"     = "#F0E442",
  "HAIC_then_T"           = "#CC79A7",
  "HAIC+I+T_concurrent"   = "#D55E00",
  "HAIC_then_I+T"         = "#56B4E9",
  "Systemic_I+T"          = "#009E73"
)
GROUP_ORDER <- c(
  "HAIC_alone", "HAIC+I_concurrent", "HAIC_then_I",
  "HAIC+T_concurrent", "HAIC_then_T",
  "HAIC+I+T_concurrent", "HAIC_then_I+T"
)
if (EIGHT_GROUP) GROUP_ORDER <- c(GROUP_ORDER, "Systemic_I+T")

# ════════════════════════════════════════════════════════════════════
# 1. 读取并整合数据
# ════════════════════════════════════════════════════════════════════
cat("\n1. 读取数据...\n")

analysis_data <- read_csv(
  file.path(DATA_DIR, DATA_CSV), show_col_types = FALSE
) %>%
  filter(os_months >= 0) %>%
  mutate(
    group = factor(main_group, levels = GROUP_ORDER),
    death_status = case_when(
      death_status %in% c("Yes","1","TRUE") ~ 1L,
      TRUE ~ 0L
    ),
    sex_male            = if_else(sex == "Male", 1L, 0L),
    pvtt_grade_cat      = case_when(
      pvtt_classification == "Absent" ~ 0L,
      pvtt_classification == "Vp1/2"  ~ 1L,
      pvtt_classification == "Vp3/4"  ~ 2L,
      TRUE ~ 0L
    ),
    pvtt_present        = if_else(pvtt_classification != "Absent", 1L, 0L),
    hvtt_present        = if_else(hvtt == "Yes", 1L, 0L),
    ivc_ra_present      = if_else(ivc_or_ra_thrombus == "Yes", 1L, 0L),
    distant_meta_bin    = if_else(distant_metastasis == "Yes", 1L, 0L),
    lymph_meta_bin      = if_else(lymph_node_metastasis == "Yes", 1L, 0L),
    ascites_bin         = if_else(ascites != "Absent", 1L, 0L),
    varices_bin         = if_else(varices == "Yes", 1L, 0L),
    tumor_gt10cm        = if_else(tumor_max_diameter_cm > 10, 1L, 0L),
    tumor_multiple      = if_else(tumor_count_category == ">3", 1L, 0L),
    albi_grade_num      = as.integer(albi_grade),
    child_pugh_grade_num = case_when(
      child_pugh_grade == "A" ~ 1L, child_pugh_grade == "B" ~ 2L, TRUE ~ 1L
    ),
    afp_log   = log10(pmax(afp,   0.01) + 1),
    pivka_log = log10(pmax(pivka, 0.01) + 1),
    tbil_log  = log10(pmax(tbil,  0.01) + 1),
    age_std        = as.numeric(scale(age)),
    afp_std        = as.numeric(scale(afp_log)),
    pivka_std      = as.numeric(scale(pivka_log)),
    tumor_size_std = as.numeric(scale(tumor_max_diameter_cm)),
    tbil_std       = as.numeric(scale(tbil_log)),
    alb_std        = as.numeric(scale(alb)),
    plt_std        = as.numeric(scale(plt)),
    nlr_std        = as.numeric(scale(nlr)),
    albi_std       = as.numeric(scale(albi_score)),
    afp_cat = case_when(afp < 20 ~ 0L, afp < 400 ~ 1L, TRUE ~ 2L),
    pivka_cat = case_when(pivka < 40 ~ 0L, pivka < 400 ~ 1L, TRUE ~ 2L)
  ) %>%
  filter(!is.na(group))

cat(sprintf("   有效患者: %d\n", nrow(analysis_data)))

# ════════════════════════════════════════════════════════════════════
# 2. 定义变量列表
# ════════════════════════════════════════════════════════════════════
# 剔除 APTT / Fibrinogen / HBsAg / PTA（用户要求）
VARS_DEMO    <- c("age", "sex", "etiology")          # 去掉 hbsag
VARS_TUMOR   <- c("tumor_max_diameter_cm", "tumor_size_category",
                   "tumor_count_category", "pvtt_classification",
                   "hvtt", "ivc_or_ra_thrombus",
                   "distant_metastasis", "lymph_node_metastasis")
VARS_LIVER   <- c("alt", "ast", "tbil", "dbil", "alb",
                   "pt", "inr")                       # 去掉 pta / aptt / fbg
VARS_HEMA    <- c("plt", "hb", "wbc", "creatinine")
VARS_MARKERS <- c("afp", "afp_high", "pivka", "pivka_high")
VARS_SCORES  <- c("albi_score", "albi_grade",
                   "child_pugh_score", "child_pugh_grade",
                   "bclc_stage", "ascites", "varices")
VARS_INFLAM  <- c("neut", "lymph", "mono", "nlr", "plr")

ALL_VARS <- c(VARS_DEMO, VARS_TUMOR, VARS_LIVER,
              VARS_HEMA, VARS_MARKERS, VARS_SCORES, VARS_INFLAM)

CONT_VARS <- c("age", "tumor_max_diameter_cm",
                "alt", "ast", "tbil", "dbil", "alb", "pt", "inr",
                "plt", "hb", "wbc", "creatinine",
                "afp", "pivka", "albi_score", "child_pugh_score",
                "neut", "lymph", "mono", "nlr", "plr")

CAT_VARS <- c("sex", "etiology",
               "tumor_size_category", "tumor_count_category",
               "pvtt_classification", "hvtt", "ivc_or_ra_thrombus",
               "distant_metastasis", "lymph_node_metastasis",
               "afp_high", "pivka_high",
               "albi_grade", "child_pugh_grade", "bclc_stage",
               "ascites", "varices")

VAR_LABELS <- list(
  age = "Age (years)", sex = "Sex", etiology = "Etiology",
  tumor_max_diameter_cm = "Max Tumor Diameter (cm)",
  tumor_size_category = "Tumor Size (≤10 / >10 cm)",
  tumor_count_category = "Tumor Number",
  pvtt_classification = "PVTT", hvtt = "HVTT",
  ivc_or_ra_thrombus = "IVC/RA Thrombus",
  distant_metastasis = "Distant Metastasis",
  lymph_node_metastasis = "Lymph Node Metastasis",
  alt = "ALT (U/L)", ast = "AST (U/L)",
  tbil = "Total Bilirubin (umol/L)", dbil = "Direct Bilirubin (umol/L)",
  alb = "Albumin (g/L)", pt = "Prothrombin Time (s)", inr = "INR",
  plt = "Platelet (x10^9/L)", hb = "Hemoglobin (g/L)",
  wbc = "WBC (x10^9/L)", creatinine = "Creatinine (umol/L)",
  afp = "AFP (ng/mL)", afp_high = "AFP > 400 ng/mL",
  pivka = "PIVKA-II (mAU/mL)", pivka_high = "PIVKA-II > 8000 mAU/mL",
  albi_score = "ALBI Score", albi_grade = "ALBI Grade",
  child_pugh_score = "Child-Pugh Score", child_pugh_grade = "Child-Pugh Grade",
  bclc_stage = "BCLC Stage", ascites = "Ascites",
  varices = "Esophagogastric Varices",
  neut = "Neutrophil (x10^9/L)", lymph = "Lymphocyte (x10^9/L)",
  mono = "Monocyte (x10^9/L)", nlr = "NLR", plr = "PLR"
)

VAR_GROUPING <- list(
  "Demographics & Etiology"        = VARS_DEMO,
  "Tumor Characteristics"          = VARS_TUMOR,
  "Liver & Synthetic Function"     = VARS_LIVER,
  "Hematology & Renal Function"    = VARS_HEMA,
  "Tumor Markers"                  = VARS_MARKERS,
  "Liver Reserve Scores & Staging" = VARS_SCORES,
  "Inflammatory Indices"           = VARS_INFLAM
)

# ── 缩写备注脚注（所有表格通用）────────────────────────────────────
ABBREV_FOOTNOTE <- paste(
  "Abbreviations: ALT, alanine aminotransferase; AST, aspartate aminotransferase;",
  "TBIL, total bilirubin; DBIL, direct bilirubin; ALB, albumin; PT, prothrombin time;",
  "INR, international normalized ratio; PLT, platelet count; HB, hemoglobin;",
  "WBC, white blood cell count; AFP, alpha-fetoprotein;",
  "PIVKA-II, protein induced by vitamin K absence or antagonist-II;",
  "ALBI, albumin-bilirubin; BCLC, Barcelona Clinic Liver Cancer;",
  "PVTT, portal vein tumor thrombus; HVTT, hepatic vein tumor thrombus;",
  "IVC/RA, inferior vena cava/right atrium; NLR, neutrophil-to-lymphocyte ratio;",
  "PLR, platelet-to-lymphocyte ratio; PSM, propensity score matching;",
  "HAIC, hepatic arterial infusion chemotherapy;",
  "I, immunotherapy; T, targeted therapy.",
  "Continuous variables are presented as median (IQR); categorical variables as n (%).",
  sep = " "
)

PSM_FORMULA <- treatment ~
  afp_cat + pivka_cat + pivka_std +
  tumor_gt10cm + tumor_multiple +
  pvtt_grade_cat + pvtt_present + hvtt_present +
  ivc_ra_present + distant_meta_bin + lymph_meta_bin +
  ascites_bin + varices_bin +
  albi_grade_num + tbil_std + alb_std + plt_std +
  age_std + tumor_size_std + nlr_std

LOVE_LABELS <- c(
  afp_cat = "AFP category", pivka_cat = "PIVKA-II category",
  pivka_std = "PIVKA-II (std)", tumor_gt10cm = "Tumor > 10 cm",
  tumor_multiple = "Tumor number > 3", pvtt_grade_cat = "PVTT grade",
  pvtt_present = "PVTT present", hvtt_present = "HVTT present",
  ivc_ra_present = "IVC/RA thrombus", distant_meta_bin = "Distant metastasis",
  lymph_meta_bin = "Lymph node metastasis", ascites_bin = "Ascites",
  varices_bin = "Esophagogastric varices", albi_grade_num = "ALBI grade",
  tbil_std = "Total bilirubin (std)", alb_std = "Albumin (std)",
  plt_std = "Platelet (std)", age_std = "Age (std)",
  tumor_size_std = "Tumor diameter (std)", nlr_std = "NLR (std)"
)

# ── 二分类变量（只显示阳性/Yes那一行，减少约1/3行数）────────────────
BINARY_YES_NO  <- c("hvtt", "ivc_or_ra_thrombus",
                    "distant_metastasis", "lymph_node_metastasis",
                    "afp_high", "pivka_high", "varices")
BINARY_POS_NEG <- character(0)   # hbsag 已从变量列表移除

# ── 紧凑型 flextable 样式（统一调用）────────────────────────────────
compact_flextable <- function(ft, font_size = 7.5, footnote = NULL) {
  ft <- ft %>%
    fontsize(size = font_size, part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    padding(padding.top = 1, padding.bottom = 1, part = "body") %>%
    padding(padding.top = 3, padding.bottom = 3, part = "header") %>%
    line_spacing(space = 1, part = "all") %>%
    set_table_properties(layout = "autofit", width = 1) %>%
    theme_booktabs()
  if (!is.null(footnote)) {
    ft <- ft %>%
      add_footer_lines(footnote) %>%
      fontsize(size = 6.5, part = "footer") %>%
      font(fontname = "Arial", part = "footer") %>%
      color(color = "#555555", part = "footer") %>%
      padding(padding.top = 4, padding.bottom = 2, part = "footer")
  }
  ft
}

# ── 页面属性 ─────────────────────────────────────────────────────────
# Table 1（7组，列多）用横向
landscape_section <- prop_section(
  page_size = page_size(orient = "landscape",
                        width  = 11.69, height = 8.27),
  page_margins = page_mar(top = 0.5, bottom = 0.5,
                          left = 0.5, right = 0.5)
)
# 描述性总表 + 21张平衡表用竖向（A4）
portrait_section <- prop_section(
  page_size = page_size(orient = "portrait",
                        width  = 8.27, height = 11.69),
  page_margins = page_mar(top = 0.8, bottom = 0.8,
                          left = 0.8, right = 0.8)
)

# ── 主汇总表生成函数 ─────────────────────────────────────────────────
make_tbl_summary <- function(data, by_var, cont_vars, cat_vars,
                              all_vars, var_labels, var_grouping,
                              add_overall_col = FALSE,
                              add_p_col = TRUE,
                              p_test_cont = "wilcox.test",
                              p_test_cat  = "chisq.test") {
  tbl <- data %>%
    select(all_of(c(by_var, all_vars))) %>%
    tbl_summary(
      by       = all_of(by_var),
      label    = var_labels,
      type     = list(
        all_of(cont_vars) ~ "continuous",
        all_of(cat_vars)  ~ "categorical"
      ),
      statistic = list(
        all_continuous()  ~ "{median} ({p25}, {p75})",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits = list(
        all_continuous()  ~ 1,
        all_categorical() ~ c(0, 1)
      ),
      missing = "no"
    )
  # add_overall / add_p 必须在 remove_row_type 之前调用
  if (add_overall_col) tbl <- tbl %>% add_overall(last = FALSE)
  if (add_p_col) {
    tbl <- tbl %>% add_p(
      test = list(
        all_continuous()  ~ p_test_cont,
        all_categorical() ~ p_test_cat
      ),
      pvalue_fun = function(x) style_pvalue(x, digits = 3)
    )
  }
  # 二分类变量只显示阳性那一行，减少行数（在 add_overall/add_p 之后调用）
  yn_vars  <- intersect(BINARY_YES_NO,  all_vars)
  pn_vars  <- intersect(BINARY_POS_NEG, all_vars)
  if (length(yn_vars) > 0) {
    tbl <- tbl %>%
      remove_row_type(variables = all_of(yn_vars), type = "level",
                      level_value = "No")
  }
  if (length(pn_vars) > 0) {
    tbl <- tbl %>%
      remove_row_type(variables = all_of(pn_vars), type = "level",
                      level_value = "Negative")
  }
  for (grp_name in names(var_grouping)) {
    grp_vars <- intersect(var_grouping[[grp_name]], all_vars)
    if (length(grp_vars) == 0) next
    tbl <- tbl %>%
      add_variable_group_header(header = grp_name, variables = all_of(grp_vars))
  }
  tbl <- tbl %>% bold_labels() %>% modify_header(label = "**Variable**")
  tbl
}

# ════════════════════════════════════════════════════════════════════
# 2. Table S0 — 全队列描述性总表（不分组，含治疗分组变量）
# ════════════════════════════════════════════════════════════════════
cat("\n2. 生成全队列描述性总表（Table S0）...\n")

# 治疗分组标签（英文展示名）
GROUP_LABELS <- c(
  "HAIC_alone"          = "HAIC alone",
  "HAIC+I_concurrent"   = "HAIC + Immunotherapy (concurrent)",
  "HAIC_then_I"         = "HAIC then Immunotherapy (sequential)",
  "HAIC+T_concurrent"   = "HAIC + Targeted therapy (concurrent)",
  "HAIC_then_T"         = "HAIC then Targeted therapy (sequential)",
  "HAIC+I+T_concurrent" = "HAIC + Immunotherapy + Targeted therapy (concurrent)",
  "HAIC_then_I+T"       = "HAIC then Immunotherapy + Targeted therapy (sequential)",
  "Systemic_I+T"        = "Systemic I+T"
)

desc_data <- analysis_data %>%
  mutate(
    afp   = round(afp,   1),
    pivka = round(pivka, 1),
    treatment_group = factor(
      GROUP_LABELS[as.character(main_group)],
      levels = unname(GROUP_LABELS)
    )
  )

VARS_TREAT   <- "treatment_group"
ALL_VARS_S0  <- c(VARS_TREAT, ALL_VARS)
CONT_VARS_S0 <- CONT_VARS
CAT_VARS_S0  <- c("treatment_group", CAT_VARS)
VAR_LABELS_S0 <- c(
  list(treatment_group = "Treatment Group"),
  VAR_LABELS
)
VAR_GROUPING_S0 <- c(
  list("Treatment Pattern" = VARS_TREAT),
  VAR_GROUPING
)

tbl_s0 <- desc_data %>%
  select(all_of(ALL_VARS_S0)) %>%
  tbl_summary(
    label    = VAR_LABELS_S0,
    type     = list(
      all_of(CONT_VARS_S0) ~ "continuous",
      all_of(CAT_VARS_S0)  ~ "categorical"
    ),
    statistic = list(
      all_continuous()  ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous()  ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    missing = "no"
  )

# 二分类变量只显示 Yes 行
yn_s0 <- intersect(BINARY_YES_NO, ALL_VARS_S0)
if (length(yn_s0) > 0) {
  tbl_s0 <- tbl_s0 %>%
    remove_row_type(variables = all_of(yn_s0), type = "level", level_value = "No")
}

# 添加分组标题
for (grp_name in names(VAR_GROUPING_S0)) {
  grp_vars <- intersect(VAR_GROUPING_S0[[grp_name]], ALL_VARS_S0)
  if (length(grp_vars) == 0) next
  tbl_s0 <- tbl_s0 %>%
    add_variable_group_header(header = grp_name, variables = all_of(grp_vars))
}

tbl_s0 <- tbl_s0 %>%
  bold_labels() %>%
  modify_header(label = "**Variable**") %>%
  modify_caption("**Table S0. Baseline Characteristics of the Overall Cohort (N = 3,885)**")

tbl_s0_ft <- as_flex_table(tbl_s0) %>%
  compact_flextable(font_size = 7.5, footnote = ABBREV_FOOTNOTE)

doc_s0 <- read_docx() %>%
  body_add_par("Table S0. Baseline Characteristics of the Overall Cohort", style = "heading 1") %>%
  body_add_flextable(tbl_s0_ft, split = TRUE) %>%
  body_end_block_section(block_section(portrait_section))
print(doc_s0, target = file.path(TABLE_DIR, "table_s0_overall_descriptive.docx"))
cat("   已保存: table_s0_overall_descriptive.docx\n")

# ════════════════════════════════════════════════════════════════════
# 3. Table 1 — 整体基线特征（7组对比）
# ════════════════════════════════════════════════════════════════════
cat("\n3. 生成 Table 1（7组整体基线特征）...\n")

tbl1 <- make_tbl_summary(
  data = analysis_data %>% mutate(afp = round(afp, 1), pivka = round(pivka, 1)),
  by_var = "group", cont_vars = CONT_VARS, cat_vars = CAT_VARS,
  all_vars = ALL_VARS, var_labels = VAR_LABELS, var_grouping = VAR_GROUPING,
  add_overall_col = TRUE, add_p_col = TRUE,
  p_test_cont = "kruskal.test", p_test_cat = "chisq.test"
) %>%
  modify_caption("**Table 1. Baseline Characteristics of All Patients Before PSM (7 Groups)**")

tbl1_ft <- as_flex_table(tbl1) %>%
  compact_flextable(font_size = 7.5, footnote = ABBREV_FOOTNOTE)

doc1 <- read_docx() %>%
  body_add_par("Table 1. Baseline Characteristics Before PSM (7 Groups)", style = "heading 1") %>%
  body_add_flextable(tbl1_ft, split = TRUE) %>%
  body_end_block_section(block_section(landscape_section))
print(doc1, target = file.path(TABLE_DIR, "table1_overall_baseline.docx"))
cat("   已保存: table1_overall_baseline.docx\n")

# ════════════════════════════════════════════════════════════════════
# 4. 21组对比：PSM + Love Plot + 平衡表
# ════════════════════════════════════════════════════════════════════
cat("\n3. 生成 21 组 Love Plot + 平衡表...\n")

all_groups <- GROUP_ORDER
n_groups   <- length(all_groups)
comparisons <- list()
comp_idx <- 1
for (i in 1:(n_groups - 1)) {
  for (j in (i+1):n_groups) {
    comparisons[[comp_idx]] <- list(
      id = comp_idx, group1 = all_groups[i], group2 = all_groups[j],
      key = paste0(all_groups[i], "_vs_", all_groups[j])
    )
    comp_idx <- comp_idx + 1
  }
}

for (comp in comparisons) {
  cid  <- comp$id
  g1   <- comp$group1
  g2   <- comp$group2
  ckey <- comp$key

  cat(sprintf("\n   === Comp %02d: %s vs %s ===\n", cid, g1, g2))

  comp_data <- analysis_data %>%
    filter(group %in% c(g1, g2)) %>%
    mutate(
      treatment   = if_else(group == g2, 1L, 0L),
      group_label = factor(group, levels = c(g1, g2))
    )

  n_treat <- sum(comp_data$treatment == 1)
  n_ctrl  <- sum(comp_data$treatment == 0)
  min_n   <- min(n_treat, n_ctrl)
  cal     <- if (min_n < 150) 0.25 else if (min_n < 300) 0.15 else 0.10

  match_result <- tryCatch(
    matchit(PSM_FORMULA, data = comp_data,
            method = "nearest", distance = "glm", link = "logit",
            caliper = cal, std.caliper = TRUE, ratio = 1, replace = FALSE),
    error = function(e) { cat("   ⚠ matchit 失败:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(match_result)) next

  matched_data <- match.data(match_result)
  n1 <- sum(matched_data$treatment == 0)
  n2 <- sum(matched_data$treatment == 1)
  cat(sprintf("   PSM 后: %s=%d, %s=%d (caliper=%.2f)\n", g1, n1, g2, n2, cal))

  # Love Plot
  love_plot_obj <- love.plot(
    match_result,
    stats        = "mean.diffs",
    threshold    = 0.1,
    abs          = TRUE,
    var.order    = "unadjusted",
    var.names    = LOVE_LABELS,
    colors       = c("#999999", "#0072B2"),
    shapes       = c("circle", "triangle"),
    size         = 3,
    title        = sprintf("Comp %02d: %s vs %s", cid, g1, g2),
    sample.names = c("Before PSM", "After PSM"),
    limits       = c(0, 0.6)
  ) +
    ggplot2::theme_bw(base_size = 9, base_family = "Helvetica") +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position    = "bottom",
      legend.title       = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(size = 10, face = "bold", hjust = 0),
      axis.text          = ggplot2::element_text(size = 8),
      axis.title         = ggplot2::element_text(size = 9)
    ) +
    ggplot2::labs(x = "Absolute Standardized Mean Difference")

  lp_base <- file.path(LOVE_DIR,
    sprintf("%02d_loveplot_comp%02d_%s_vs_%s",
            cid, cid, gsub("[+]","_", g1), gsub("[+]","_", g2)))
  ggplot2::ggsave(paste0(lp_base, ".pdf"), love_plot_obj,
                  width = 5.5, height = 5.5, device = cairo_pdf)
  ggplot2::ggsave(paste0(lp_base, ".png"), love_plot_obj,
                  width = 5.5, height = 5.5, dpi = 300, type = "cairo")
  cat(sprintf("   Love Plot 已保存: %s\n", basename(lp_base)))

  # PSM 前后平衡表
  comp_data_tbl    <- comp_data    %>% mutate(treatment = factor(treatment, levels=c(0L,1L), labels=c(g1,g2)))
  matched_data_tbl <- matched_data %>% mutate(treatment = factor(treatment, levels=c(0L,1L), labels=c(g1,g2)))

  tbl_before <- make_tbl_summary(
    data = comp_data_tbl, by_var = "treatment",
    cont_vars = CONT_VARS, cat_vars = CAT_VARS,
    all_vars = ALL_VARS, var_labels = VAR_LABELS, var_grouping = VAR_GROUPING,
    add_overall_col = FALSE, add_p_col = TRUE
  ) %>%
    modify_spanning_header(all_stat_cols() ~ "**Before PSM**")

  tbl_after <- make_tbl_summary(
    data = matched_data_tbl, by_var = "treatment",
    cont_vars = CONT_VARS, cat_vars = CAT_VARS,
    all_vars = ALL_VARS, var_labels = VAR_LABELS, var_grouping = VAR_GROUPING,
    add_overall_col = FALSE, add_p_col = TRUE
  ) %>%
    bold_labels() %>%
    modify_header(label = "**Variable**") %>%
    modify_spanning_header(all_stat_cols() ~ "**After PSM**")

  tbl_merged <- tbl_merge(
    tbls = list(tbl_before, tbl_after),
    tab_spanner = c(
      sprintf("**Before PSM** (n=%d+%d)", n_ctrl, n_treat),
      sprintf("**After PSM** (n=%d+%d)", n1, n2)
    )
  ) %>%
    modify_caption(sprintf("**Table %d. Baseline Balance: %s vs %s**", cid+1, g1, g2))

  tbl_ft <- as_flex_table(tbl_merged) %>%
    compact_flextable(font_size = 7.5, footnote = ABBREV_FOOTNOTE)

  doc <- read_docx() %>%
    body_add_par(sprintf("Table %d. Baseline Balance: %s vs %s", cid+1, g1, g2),
                 style = "heading 1") %>%
    body_add_flextable(tbl_ft, split = TRUE) %>%
    body_end_block_section(block_section(portrait_section))

  tbl_fname <- file.path(TABLE_DIR,
    sprintf("table%02d_comp%02d_%s_vs_%s_psm_balance.docx",
            cid+1, cid, gsub("[+]","_", g1), gsub("[+]","_", g2)))
  print(doc, target = tbl_fname)
  cat(sprintf("   平衡表已保存: %s\n", basename(tbl_fname)))
}

cat("\n============================================================\n")
cat("全部完成！\n")
cat(sprintf("表格目录: %s\n", TABLE_DIR))
cat(sprintf("  - table_s0_overall_descriptive.docx  (全队列描述性总表)\n"))
cat(sprintf("  - table1_overall_baseline.docx        (7组基线对比，横向)\n"))
cat(sprintf("  - table02~22 各组 PSM 平衡表 (竖向，含脚注)\n"))
cat(sprintf("Love Plot 目录: %s\n", LOVE_DIR))
cat("============================================================\n")

sink()
