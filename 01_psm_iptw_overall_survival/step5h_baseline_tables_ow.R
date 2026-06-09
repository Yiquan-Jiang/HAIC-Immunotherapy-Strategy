#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
#
# step5h — Before/after overlap-weighting baseline (Table-1) balance tables
# =========================================================================
# Companion to step5e/step5f/step5g. For EACH of the 7 HAIC strategies vs Systemic I+T,
# a baseline-characteristics table in the SAME house style as the step6 PSM balance tables
# (gtsummary + flextable, grouped sections, median [IQR] / n (%), per-row SMD), with:
#   - "Before weighting"  = unweighted distributions + SMD
#   - "After OW (ATO)"     = overlap-weighted distributions (survey design) + SMD
# Output (results/ow_vs_systemic_it_8group/tables/):
#   ow_tableXX_<group>_vs_Systemic_I+T_balance.docx
#
# Notes: dbil/creatinine/mono are not collected in the Systemic cohort and are omitted.
# Estimand = ATO (overlap population). Reads ow_weights.csv + ow_forest_data.csv from step5e.

suppressMessages({
  library(dplyr); library(readr); library(survey)
  library(gtsummary); library(flextable); library(officer)
})

BASE_DIR <- "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
DATA_DIR <- file.path(BASE_DIR, "data")
OW_DIR   <- file.path(BASE_DIR, "results", "ow_vs_systemic_it_8group")
TAB_DIR  <- file.path(OW_DIR, "tables")
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

# ── variable sets (step6 ALL_VARS minus dbil/creatinine/mono) ───────────────
VARS_DEMO    <- c("age", "sex", "etiology")
VARS_TUMOR   <- c("tumor_max_diameter_cm", "tumor_size_category", "tumor_count_category",
                  "pvtt_classification", "hvtt", "ivc_or_ra_thrombus",
                  "distant_metastasis", "lymph_node_metastasis")
VARS_LIVER   <- c("alt", "ast", "tbil", "alb", "pt", "inr")
VARS_HEMA    <- c("plt", "hb", "wbc")
VARS_MARKERS <- c("afp", "afp_high", "pivka", "pivka_high")
VARS_SCORES  <- c("albi_score", "albi_grade", "child_pugh_score", "child_pugh_grade",
                  "bclc_stage", "ascites", "varices")
VARS_INFLAM  <- c("neut", "lymph", "nlr", "plr")
ALL_VARS  <- c(VARS_DEMO, VARS_TUMOR, VARS_LIVER, VARS_HEMA, VARS_MARKERS,
               VARS_SCORES, VARS_INFLAM)
CONT_VARS <- c("age", "tumor_max_diameter_cm", "alt", "ast", "tbil", "alb", "pt", "inr",
               "plt", "hb", "wbc", "afp", "pivka", "albi_score", "child_pugh_score",
               "neut", "lymph", "nlr", "plr")
CAT_VARS  <- c("sex", "etiology", "tumor_size_category", "tumor_count_category",
               "pvtt_classification", "hvtt", "ivc_or_ra_thrombus", "distant_metastasis",
               "lymph_node_metastasis", "afp_high", "pivka_high", "albi_grade",
               "child_pugh_grade", "bclc_stage", "ascites", "varices")

VAR_LABELS <- list(
  age = "Age (years)", sex = "Sex", etiology = "Etiology",
  tumor_max_diameter_cm = "Max Tumor Diameter (cm)",
  tumor_size_category = "Tumor Size (≤10 / >10 cm)", tumor_count_category = "Tumor Number",
  pvtt_classification = "PVTT", hvtt = "HVTT", ivc_or_ra_thrombus = "IVC/RA Thrombus",
  distant_metastasis = "Distant Metastasis", lymph_node_metastasis = "Lymph Node Metastasis",
  alt = "ALT (U/L)", ast = "AST (U/L)", tbil = "Total Bilirubin (umol/L)",
  alb = "Albumin (g/L)", pt = "Prothrombin Time (s)", inr = "INR",
  plt = "Platelet (x10^9/L)", hb = "Hemoglobin (g/L)", wbc = "WBC (x10^9/L)",
  afp = "AFP (ng/mL)", afp_high = "AFP > 400 ng/mL", pivka = "PIVKA-II (mAU/mL)",
  pivka_high = "PIVKA-II > 8000 mAU/mL", albi_score = "ALBI Score", albi_grade = "ALBI Grade",
  child_pugh_score = "Child-Pugh Score", child_pugh_grade = "Child-Pugh Grade",
  bclc_stage = "BCLC Stage", ascites = "Ascites", varices = "Esophagogastric Varices",
  neut = "Neutrophil (x10^9/L)", lymph = "Lymphocyte (x10^9/L)", nlr = "NLR", plr = "PLR"
)
VAR_GROUPING <- list(
  "Demographics & Etiology" = VARS_DEMO, "Tumor Characteristics" = VARS_TUMOR,
  "Liver & Synthetic Function" = VARS_LIVER, "Hematology" = VARS_HEMA,
  "Tumor Markers" = VARS_MARKERS, "Liver Reserve Scores & Staging" = VARS_SCORES,
  "Inflammatory Indices" = VARS_INFLAM
)
BINARY_YES_NO <- c("hvtt", "ivc_or_ra_thrombus", "distant_metastasis",
                   "lymph_node_metastasis", "afp_high", "pivka_high", "varices")
ABBREV_FOOTNOTE <- paste(
  "Values are median (IQR) for continuous and n (%) for categorical variables.",
  "SMD, absolute standardized mean difference (|SMD| < 0.1 indicates good balance);",
  "OW, overlap weighting; ATO, average treatment effect in the overlap population;",
  "ESS, effective sample size. dbil/creatinine/monocyte not collected in the systemic",
  "cohort and omitted. Abbreviations as in the main baseline table.")

compact_flextable <- function(ft, font_size = 7.5, footnote = NULL) {
  ft <- ft %>% fontsize(size = font_size, part = "all") %>%
    font(fontname = "Arial", part = "all") %>%
    padding(padding.top = 1, padding.bottom = 1, part = "body") %>%
    padding(padding.top = 3, padding.bottom = 3, part = "header") %>%
    line_spacing(space = 1, part = "all") %>%
    set_table_properties(layout = "autofit", width = 1) %>% theme_booktabs()
  if (!is.null(footnote)) {
    ft <- ft %>% add_footer_lines(footnote) %>%
      fontsize(size = 6.5, part = "footer") %>% font(fontname = "Arial", part = "footer") %>%
      color(color = "#555555", part = "footer")
  }
  ft
}
portrait_section <- prop_section(
  page_size = page_size(orient = "portrait", width = 8.27, height = 11.69),
  page_margins = page_mar(top = 0.8, bottom = 0.8, left = 0.8, right = 0.8))

# tbl builder: x is a data.frame (unweighted) or survey.design (weighted)
make_summary <- function(x, weighted) {
  base <- if (weighted) {
    tbl_svysummary(x, by = arm, include = all_of(ALL_VARS), label = VAR_LABELS,
                   type = list(all_of(CONT_VARS) ~ "continuous", all_of(CAT_VARS) ~ "categorical"),
                   statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                                    all_categorical() ~ "{n} ({p}%)"),
                   digits = list(all_continuous() ~ 1, all_categorical() ~ c(0, 1)),
                   missing = "no")
  } else {
    tbl_summary(x, by = arm, include = all_of(ALL_VARS), label = VAR_LABELS,
                type = list(all_of(CONT_VARS) ~ "continuous", all_of(CAT_VARS) ~ "categorical"),
                statistic = list(all_continuous() ~ "{median} ({p25}, {p75})",
                                 all_categorical() ~ "{n} ({p}%)"),
                digits = list(all_continuous() ~ 1, all_categorical() ~ c(0, 1)),
                missing = "no")
  }
  tbl <- base %>%
    add_difference(everything() ~ "smd",
                   estimate_fun = everything() ~ function(x) style_number(abs(x), digits = 3))
  tbl <- tbl %>% modify_column_hide(c(conf.low, conf.high))
  if (weighted) tbl <- tbl %>% modify_header(all_stat_cols() ~ "**{level}**")
  yn <- intersect(BINARY_YES_NO, ALL_VARS)
  if (length(yn) > 0)
    tbl <- tbl %>% remove_row_type(variables = all_of(yn), type = "level", level_value = "No")
  for (grp in names(VAR_GROUPING)) {
    gv <- intersect(VAR_GROUPING[[grp]], ALL_VARS)
    if (length(gv) > 0)
      tbl <- tbl %>% add_variable_group_header(header = grp, variables = all_of(gv))
  }
  tbl %>% bold_labels() %>% modify_header(label = "**Variable**", estimate = "**SMD**")
}

# ── data ────────────────────────────────────────────────────────────────────
d <- read_csv(file.path(DATA_DIR, "analysis_ready_8group.csv"), show_col_types = FALSE) %>%
  mutate(afp = round(afp, 1), pivka = round(pivka, 1))
wts <- read_csv(file.path(OW_DIR, "ow_weights.csv"), show_col_types = FALSE)
fd  <- read_csv(file.path(OW_DIR, "ow_forest_data.csv"), show_col_types = FALSE)
GLAB <- c("HAIC_alone" = "HAIC alone", "HAIC+I_concurrent" = "HAIC + I (concurrent)",
          "HAIC_then_I" = "HAIC → I (deferred)", "HAIC+T_concurrent" = "HAIC + T (concurrent)",
          "HAIC_then_T" = "HAIC → T (deferred)", "HAIC+I+T_concurrent" = "HAIC + I + T (concurrent)",
          "HAIC_then_I+T" = "HAIC → I + T (deferred)")
HAIC_ORDER <- names(GLAB)

for (i in seq_along(HAIC_ORDER)) {
  g <- HAIC_ORDER[i]
  glab <- GLAB[[g]]
  sub <- wts %>% filter(group == g) %>% select(patient_id, ow) %>%
    inner_join(d, by = "patient_id") %>%
    mutate(arm = factor(if_else(main_group == "Systemic_I+T", "Systemic I+T", glab),
                        levels = c(glab, "Systemic I+T")))
  n1 <- sum(sub$arm == glab); n2 <- sum(sub$arm == "Systemic I+T")
  r  <- fd[fd$group == g, ]

  tbl_before <- make_summary(sub, weighted = FALSE) %>%
    modify_spanning_header(all_stat_cols() ~ "**Before weighting**")
  des <- svydesign(~1, data = sub, weights = ~ow)
  tbl_after <- make_summary(des, weighted = TRUE) %>%
    modify_spanning_header(all_stat_cols() ~ "**After OW (ATO)**")

  tbl_merged <- tbl_merge(
    tbls = list(tbl_before, tbl_after),
    tab_spanner = c(sprintf("**Before weighting** (n=%d+%d)", n1, n2),
                    sprintf("**After OW · ATO** (ESS=%d+%d)", r$ess_haic, r$ess_sys))) %>%
    modify_caption(sprintf("**Baseline balance: %s vs Systemic I+T (overlap weighting)**", glab))

  ft <- as_flex_table(tbl_merged) %>% compact_flextable(font_size = 7.5, footnote = ABBREV_FOOTNOTE)
  fname <- sprintf("ow_table%02d_%s_vs_Systemic_I+T_balance.docx", i, gsub("[+]", "_", g))
  read_docx() %>%
    body_add_par(sprintf("Baseline balance: %s vs Systemic I+T", glab), style = "heading 1") %>%
    body_add_flextable(ft, split = TRUE) %>%
    body_end_block_section(block_section(portrait_section)) %>%
    print(target = file.path(TAB_DIR, fname))
  cat(sprintf("  saved: %s  (before %d+%d -> after ESS %d+%d, max|SMD| %.3f)\n",
              fname, n1, n2, r$ess_haic, r$ess_sys, r$max_smd_adj))
}
cat("\nDone ->", TAB_DIR, "\n")
