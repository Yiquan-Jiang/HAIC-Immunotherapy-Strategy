#!/usr/bin/env Rscript
# ============================================================
# OW 加权后基线平衡表
# 对比1: HAIC+I_concurrent vs HAIC_then_I
# 对比2: HAIC+I+T_concurrent vs HAIC_then_I+T
# 输出：每个亚组的加权均值/比例 + SMD（加权前后对比）
# ============================================================

library(tidyverse)
library(WeightIt)
library(cobalt)

BASE_DIR <- "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
DATA_DIR <- file.path(BASE_DIR, "data")
RES_DIR  <- file.path(BASE_DIR, "results", "subgroup_analysis")
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

cat("============================================================\n")
cat("  OW 加权后基线平衡表生成\n")
cat("============================================================\n\n")

# ── 读取数据 ──────────────────────────────────────────────────
df_all_raw <- read_csv(file.path(DATA_DIR, "analysis_ready.csv"), show_col_types = FALSE) %>%
  filter(os_months >= 0)

COMPARISONS <- list(
  list(tag = "I",
       treat_group   = "HAIC+I_concurrent",
       control_group = "HAIC_then_I"),
  list(tag = "IT",
       treat_group   = "HAIC+I+T_concurrent",
       control_group = "HAIC_then_I+T")
)

prepare_two_group_bt <- function(df_raw, treat_group, control_group) {
  df_raw %>%
    filter(main_group %in% c(treat_group, control_group)) %>%
    mutate(
      treat          = if_else(main_group == treat_group, 1L, 0L),
    event          = case_when(death_status %in% c("Yes","1","TRUE") ~ 1L, TRUE ~ 0L),
    sex_male       = if_else(sex == "Male", 1L, 0L),
    pvtt_vp34      = if_else(pvtt_classification == "Vp3/4", 1L, 0L),
    hvtt_yes       = if_else(hvtt == "Yes", 1L, 0L),
    ivc_ra_yes     = if_else(ivc_or_ra_thrombus == "Yes", 1L, 0L),
    dist_meta_yes  = if_else(distant_metastasis == "Yes", 1L, 0L),
    lymph_meta_yes = if_else(lymph_node_metastasis == "Yes", 1L, 0L),
    ascites_yes    = if_else(ascites != "Absent", 1L, 0L),
    varices_yes    = if_else(varices == "Yes", 1L, 0L),
    tumor_gt10     = if_else(tumor_max_diameter_cm > 10, 1L, 0L),
    tumor_multi    = if_else(tumor_count_category == ">3", 1L, 0L),
    afp_high_bin   = if_else(afp_high == "Yes", 1L, 0L),
    pivka_high_bin = if_else(pivka_high == "Yes", 1L, 0L),
    high_risk_composite = as.integer(ivc_ra_yes | tumor_multi | pvtt_vp34 |
                                     dist_meta_yes | tumor_gt10)
  )
}

# ════════════════════════════════════════════════════════════════
# 遍历每组对比
# ════════════════════════════════════════════════════════════════
for (COMP in COMPARISONS) {

CMP_TAG       <- COMP$tag
TREAT_GROUP   <- COMP$treat_group
CONTROL_GROUP <- COMP$control_group

CMP_RES_DIR <- file.path(RES_DIR, CMP_TAG)
dir.create(CMP_RES_DIR, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("\n######################################################\n"))
cat(sprintf("#  对比组: %s vs %s  (tag=%s)\n", TREAT_GROUP, CONTROL_GROUP, CMP_TAG))
cat(sprintf("######################################################\n"))

df <- prepare_two_group_bt(df_all_raw, TREAT_GROUP, CONTROL_GROUP)

# ── 核心 PS 变量 ───────────────────────────────────────────────
CORE_PS_VARS <- c(
  "age", "sex_male", "albi_score",
  "tumor_max_diameter_cm", "tumor_multi",
  "pvtt_vp34", "hvtt_yes",
  "dist_meta_yes", "lymph_meta_yes",
  "ascites_yes", "afp_high_bin", "pivka_high_bin",
  "nlr", "plt", "alb"
)

# ── 展示变量（临床基线表）────────────────────────────────────
# 格式：list(var_name, display_label, type)  type = "cont" | "bin"
DISPLAY_VARS <- list(
  list("age",               "Age (years)",              "cont"),
  list("sex_male",          "Male sex",                 "bin"),
  list("albi_score",        "ALBI score",               "cont"),
  list("child_pugh_score",  "Child-Pugh score",         "cont"),
  list("alt",               "ALT (U/L)",                "cont"),
  list("ast",               "AST (U/L)",                "cont"),
  list("tbil",              "TBIL (μmol/L)",            "cont"),
  list("alb",               "Albumin (g/L)",            "cont"),
  list("inr",               "INR",                      "cont"),
  list("plt",               "Platelets (x10^9/L)",      "cont"),
  list("creatinine",        "Creatinine (μmol/L)",      "cont"),
  list("afp_high_bin",      "AFP >400 ng/mL",           "bin"),
  list("pivka_high_bin",    "PIVKA-II >8000 mAU/mL",   "bin"),
  list("tumor_max_diameter_cm", "Max tumor diameter (cm)", "cont"),
  list("tumor_gt10",        "Tumor diameter >10 cm",    "bin"),
  list("tumor_multi",       "Tumor count >3",           "bin"),
  list("pvtt_vp34",         "PVTT Vp3/4",               "bin"),
  list("hvtt_yes",          "HVTT",                     "bin"),
  list("dist_meta_yes",     "Extrahepatic metastasis",  "bin"),
  list("lymph_meta_yes",    "Lymph node metastasis",    "bin"),
  list("ascites_yes",       "Ascites",                  "bin"),
  list("varices_yes",       "Esophagogastric varices",  "bin"),
  list("nlr",               "NLR",                      "cont"),
  list("plr",               "PLR",                      "cont")
)

# ── 亚组定义 ──────────────────────────────────────────────────
SUBGROUPS <- list(
  list(name = "Composite high-risk",
       filter_expr = quote(high_risk_composite == 1),
       exclude_ps  = character(0)),
  list(name = "Tumor count >3",
       filter_expr = quote(tumor_multi == 1),
       exclude_ps  = c("tumor_multi")),
  list(name = "Tumor diameter >10 cm",
       filter_expr = quote(tumor_gt10 == 1),
       exclude_ps  = c("tumor_gt10", "tumor_max_diameter_cm")),
  list(name = "PVTT Vp3/4",
       filter_expr = quote(pvtt_vp34 == 1),
       exclude_ps  = c("pvtt_vp34")),
  list(name = "Extrahepatic metastasis",
       filter_expr = quote(dist_meta_yes == 1),
       exclude_ps  = c("dist_meta_yes"))
)

# ── 辅助函数 ──────────────────────────────────────────────────

# 加权均值（连续）或加权比例（二分类）
weighted_stat <- function(x, w, type) {
  if (all(is.na(x))) return(c(val = NA_real_, sd = NA_real_))
  idx <- !is.na(x)
  x_c <- x[idx]; w_c <- w[idx]
  w_c <- w_c / sum(w_c)
  if (type == "cont") {
    mu  <- sum(w_c * x_c)
    vr  <- sum(w_c * (x_c - mu)^2) / (1 - sum(w_c^2))
    return(c(val = mu, sd = sqrt(max(vr, 0))))
  } else {
    p <- sum(w_c * x_c)
    return(c(val = p * 100, sd = NA_real_))
  }
}

# 计算 SMD（加权后）
calc_smd <- function(x, treat, w, type) {
  x1 <- x[treat == 1]; w1 <- w[treat == 1]
  x0 <- x[treat == 0]; w0 <- w[treat == 0]
  if (type == "cont") {
    s1 <- weighted_stat(x1, w1, "cont")
    s0 <- weighted_stat(x0, w0, "cont")
    pool_sd <- sqrt((s1["sd"]^2 + s0["sd"]^2) / 2)
    if (is.na(pool_sd) || pool_sd == 0) return(NA_real_)
    return((s1["val"] - s0["val"]) / pool_sd)
  } else {
    p1 <- weighted_stat(x1, w1, "bin")["val"] / 100
    p0 <- weighted_stat(x0, w0, "bin")["val"] / 100
    denom <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
    if (is.na(denom) || denom == 0) return(NA_real_)
    return((p1 - p0) / denom)
  }
}

# ── 主循环 ────────────────────────────────────────────────────
all_tables <- list()

for (sg in SUBGROUPS) {
  cat(sprintf("\n--- %s ---\n", sg$name))

  sg_data <- df %>% filter(!!sg$filter_expr)
  n_conc  <- sum(sg_data$treat == 1)
  n_then  <- sum(sg_data$treat == 0)
  cat(sprintf("  n = %d+%d\n", n_conc, n_then))

  ps_vars <- setdiff(CORE_PS_VARS, sg$exclude_ps)
  ps_vars <- ps_vars[sapply(ps_vars, function(v) sd(sg_data[[v]], na.rm = TRUE) > 0)]
  ps_formula <- as.formula(paste("treat ~", paste(ps_vars, collapse = " + ")))

  W <- tryCatch(
    weightit(ps_formula, data = sg_data, method = "glm", estimand = "ATO"),
    error = function(e) { cat("  WeightIt failed\n"); NULL }
  )
  if (is.null(W)) next

  sg_data$ow_w <- W$weights

  rows <- list()
  for (dv in DISPLAY_VARS) {
    vname <- dv[[1]]; vlabel <- dv[[2]]; vtype <- dv[[3]]

    if (!vname %in% names(sg_data)) next

    x  <- sg_data[[vname]]
    tr <- sg_data$treat
    w  <- sg_data$ow_w

    # 未加权
    s1_un <- weighted_stat(x[tr == 1], rep(1, sum(tr == 1)), vtype)
    s0_un <- weighted_stat(x[tr == 0], rep(1, sum(tr == 0)), vtype)
    smd_un <- calc_smd(x, tr, rep(1, nrow(sg_data)), vtype)

    # 加权后
    s1_ow <- weighted_stat(x[tr == 1], w[tr == 1], vtype)
    s0_ow <- weighted_stat(x[tr == 0], w[tr == 0], vtype)
    smd_ow <- calc_smd(x, tr, w, vtype)

    if (vtype == "cont") {
      conc_un  <- sprintf("%.1f (%.1f)", s1_un["val"], s1_un["sd"])
      then_un  <- sprintf("%.1f (%.1f)", s0_un["val"], s0_un["sd"])
      conc_ow  <- sprintf("%.1f (%.1f)", s1_ow["val"], s1_ow["sd"])
      then_ow  <- sprintf("%.1f (%.1f)", s0_ow["val"], s0_ow["sd"])
    } else {
      conc_un  <- sprintf("%.1f%%", s1_un["val"])
      then_un  <- sprintf("%.1f%%", s0_un["val"])
      conc_ow  <- sprintf("%.1f%%", s1_ow["val"])
      then_ow  <- sprintf("%.1f%%", s0_ow["val"])
    }

    rows[[length(rows) + 1]] <- data.frame(
      Subgroup        = sg$name,
      Variable        = vlabel,
      Type            = vtype,
      N_concurrent    = n_conc,
      N_then_I        = n_then,
      Concurrent_unw  = conc_un,
      Then_I_unw      = then_un,
      SMD_unweighted  = round(abs(smd_un), 3),
      Concurrent_OW   = conc_ow,
      Then_I_OW       = then_ow,
      SMD_OW          = round(abs(smd_ow), 3),
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) > 0) {
    sg_table <- bind_rows(rows)
    all_tables[[sg$name]] <- sg_table
    cat(sprintf("  变量数: %d, 加权后 SMD>0.1: %d\n",
                nrow(sg_table),
                sum(sg_table$SMD_OW > 0.1, na.rm = TRUE)))
  }
}

# ── 合并并保存 ────────────────────────────────────────────────
full_table <- bind_rows(all_tables)
out_path   <- file.path(CMP_RES_DIR, paste0("ow_balance_table_full_", CMP_TAG, ".csv"))
write_csv(full_table, out_path)
cat(sprintf("\n全部完成，共 %d 行 → %s\n", nrow(full_table), out_path))

# ── 宽表（每个亚组一列）────────────────────────────────────────
# 用于 Python 渲染多列对比表
wide_list <- list()
for (sg_name in names(all_tables)) {
  t <- all_tables[[sg_name]]
  safe <- gsub("[ >≥/]", "_", sg_name)
  col_conc <- paste0(safe, "_Concurrent_OW")
  col_then <- paste0(safe, "_ThenI_OW")
  col_smd  <- paste0(safe, "_SMD_OW")
  sub_wide <- t %>%
    select(Variable, !!col_conc := Concurrent_OW,
           !!col_then := Then_I_OW,
           !!col_smd  := SMD_OW)
  wide_list[[sg_name]] <- sub_wide
}

wide_table <- reduce(wide_list, full_join, by = "Variable")
write_csv(wide_table, file.path(CMP_RES_DIR, paste0("ow_balance_table_wide_", CMP_TAG, ".csv")))
cat(sprintf("宽表已保存 → %s\n",
            file.path(CMP_RES_DIR, paste0("ow_balance_table_wide_", CMP_TAG, ".csv"))))

}  # end for (COMP in COMPARISONS)
