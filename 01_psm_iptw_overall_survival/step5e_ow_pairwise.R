#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
#
# step5e — Overlap Weighting (ATO) of each group vs a reference arm
# =====================================================================
# Focused pairwise overlap-weighting (ATO): for EACH non-reference group we run a binary
# OW contrast against the reference arm REF. Default REF = HAIC_alone, so the 7 rows are
# the 6 other HAIC strategies + Systemic_I+T, each weighted against HAIC alone.
# (Set OW_REF=Systemic_I+T to reproduce the original vs-Systemic track, which fixes the
#  no-HAIC arm's inadequate overlap under the joint 8-group IPTW — ESS 182/570.)
#   - exact mean balance on the PS-model covariates (overlap/equipoise population)
#   - varices DROPPED (uniformly 'No' in the Systemic source -> artifactual imbalance;
#     kept dropped for ALL contrasts so one PS model is shared across the forest)
#   - weighted Cox HR (row group vs REF; HR<1 = row group better than REF)
#   - per-contrast love plot + balance summary + ESS
# Estimand: ATO (average treatment effect in the overlap population).
#
# Outputs (results/ow_vs_<refkey>_8group/, e.g. ow_vs_haic_alone_8group/):
#   ow_forest_data.csv        7 rows: HR/CI/p, ESS, max|SMD| pre/post
#                             (n_haic/ess_haic = row-group n/ESS; n_sys/ess_sys = REF n/ESS)
#   ow_balance_long.csv       per-covariate SMD (unadj/adj) for every contrast
#   ow_km_data.csv            ATO-weighted KM curve data per contrast (for plotting)
# Love plots (figures/ow_vs_<refkey>_8group/): love_<group>.png

suppressMessages({
  library(dplyr); library(readr); library(WeightIt); library(cobalt)
  library(survival); library(ggplot2)
})

BASE_DIR <- "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30/HAIC_NO_TACE_4_TIDY/update_group_7"
DATA_DIR <- file.path(BASE_DIR, "data")

# ── Reference arm (the comparator each row group is weighted against) ───────
# Default HAIC_alone; set OW_REF=Systemic_I+T to reproduce the original vs-Systemic track.
REF <- Sys.getenv("OW_REF", "HAIC_alone")

ALL_GROUPS <- c("HAIC_alone", "HAIC+I_concurrent", "HAIC_then_I", "HAIC+T_concurrent",
                "HAIC_then_T", "HAIC+I+T_concurrent", "HAIC_then_I+T", "Systemic_I+T")
ROW_GROUPS <- setdiff(ALL_GROUPS, REF)   # every other group, contrasted against REF

group_label <- c(
  "HAIC_alone" = "HAIC alone", "HAIC+I_concurrent" = "HAIC + I (concurrent)",
  "HAIC_then_I" = "HAIC -> I (deferred)", "HAIC+T_concurrent" = "HAIC + T (concurrent)",
  "HAIC_then_T" = "HAIC -> T (deferred)", "HAIC+I+T_concurrent" = "HAIC + I + T (concurrent)",
  "HAIC_then_I+T" = "HAIC -> I + T (deferred)", "Systemic_I+T" = "Systemic I+T")
REF_LABEL <- group_label[[REF]]

ref_key <- function(x) {
  if (x == "Systemic_I+T") "systemic_it"
  else if (x == "HAIC_alone") "haic_alone"
  else gsub("[^A-Za-z0-9]+", "_", tolower(x))
}
OUT_DIR  <- file.path(BASE_DIR, "results", paste0("ow_vs_", ref_key(REF), "_8group"))
FIG_DIR  <- file.path(BASE_DIR, "figures", paste0("ow_vs_", ref_key(REF), "_8group"))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Reference arm: %s  ->  %s\n", REF, OUT_DIR))

raw <- read_csv(file.path(DATA_DIR, "analysis_ready_8group.csv"), show_col_types = FALSE)

prep <- raw %>%
  mutate(
    death = if_else(death_status %in% c("Yes","1","TRUE","yes"), 1L, 0L),
    afp_cat = factor(case_when(afp < 20 ~ 0L, afp < 400 ~ 1L, TRUE ~ 2L)),
    pivka_log = log10(pmax(pivka, 0.01) + 1),
    tbil_log  = log10(pmax(tbil,  0.01) + 1),
    tumor_gt10cm   = if_else(tumor_max_diameter_cm > 10, 1L, 0L),
    tumor_multiple = if_else(tumor_count_category == ">3", 1L, 0L),
    pvtt_grade_cat = factor(case_when(
      pvtt_classification == "Absent" ~ 0L, pvtt_classification == "Vp1/2" ~ 1L,
      pvtt_classification == "Vp3/4"  ~ 2L, TRUE ~ 0L)),
    hvtt_present     = if_else(hvtt == "Yes", 1L, 0L),
    ivc_ra_present   = if_else(ivc_or_ra_thrombus == "Yes", 1L, 0L),
    distant_meta_bin = if_else(distant_metastasis == "Yes", 1L, 0L),
    lymph_meta_bin   = if_else(lymph_node_metastasis == "Yes", 1L, 0L),
    ascites_bin      = if_else(ascites != "Absent", 1L, 0L),
    albi_grade_fac   = factor(as.integer(albi_grade))
  ) %>%
  filter(!is.na(os_months), os_months > 0, !is.na(death))

# PS-model covariates (step3b set MINUS varices_bin)
COVS <- c("afp_cat","pivka_log","tumor_gt10cm","tumor_multiple","pvtt_grade_cat",
          "hvtt_present","ivc_ra_present","distant_meta_bin","lymph_meta_bin",
          "ascites_bin","albi_grade_fac","tbil_log","alb","plt","age",
          "tumor_max_diameter_cm","nlr")
ow_formula <- as.formula(paste("treat ~", paste(COVS, collapse = " + ")))

var_labels <- c(afp_cat = "AFP category", pivka_log = "log PIVKA-II",
  tumor_gt10cm = "Tumor >10 cm", tumor_multiple = "Tumor >3 nodules",
  pvtt_grade_cat = "PVTT grade", hvtt_present = "HVTT", ivc_ra_present = "IVC/RA thrombus",
  distant_meta_bin = "Distant metastasis", lymph_meta_bin = "Lymph-node metastasis",
  ascites_bin = "Ascites", albi_grade_fac = "ALBI grade", tbil_log = "log TBIL",
  alb = "Albumin", plt = "Platelets", age = "Age",
  tumor_max_diameter_cm = "Max tumor diameter", nlr = "NLR")

ess <- function(w) sum(w)^2 / sum(w^2)
forest_rows <- list(); bal_rows <- list(); km_rows <- list(); wt_rows <- list()

for (g in ROW_GROUPS) {
  d <- prep %>% filter(main_group %in% c(REF, g)) %>%
    mutate(treat = if_else(main_group == g, 1L, 0L))   # 1 = row group, 0 = reference (REF)
  d <- d %>% filter(if_all(all_of(COVS), ~ !is.na(.)))

  W <- weightit(ow_formula, data = d, method = "glm", estimand = "ATO")
  d$ow <- W$weights
  wt_rows[[g]] <- d %>% transmute(group = g, patient_id, main_group, treat,
                                  os_months, death, ow)

  bt <- bal.tab(W, un = TRUE, stats = "mean.diffs", binary = "std", continuous = "std")$Balance
  bt$covariate <- rownames(bt)
  bal_rows[[g]] <- data.frame(group = g, covariate = bt$covariate,
                              smd_unadj = bt$Diff.Un, smd_adj = bt$Diff.Adj)

  # love plot
  lp <- love.plot(W, stats = "mean.diffs", binary = "std", abs = TRUE,
                  thresholds = c(m = 0.1), var.names = var_labels,
                  title = sprintf("%s vs %s (ATO)", group_label[[g]], REF_LABEL),
                  sample.names = c("Unweighted", "Overlap-weighted")) +
        theme(legend.position = "bottom")
  ggsave(file.path(FIG_DIR, sprintf("love_%s.png", gsub("[^A-Za-z0-9]+","_", g))),
         lp, width = 7, height = 5.2, dpi = 300)

  # weighted Cox: HR of row group vs REF (reference)
  cox <- coxph(Surv(os_months, death) ~ treat, data = d, weights = ow, robust = TRUE)
  sm <- summary(cox)
  hr <- sm$conf.int[1, "exp(coef)"]; lo <- sm$conf.int[1, "lower .95"]
  hi <- sm$conf.int[1, "upper .95"]; pv <- sm$coefficients[1, "Pr(>|z|)"]

  # NOTE column names kept for backward compat: n_haic/ess_haic = ROW-group, n_sys/ess_sys = REF
  forest_rows[[g]] <- data.frame(
    group = g, n_haic = sum(d$treat == 1), n_sys = sum(d$treat == 0),
    ess_haic = round(ess(d$ow[d$treat == 1])), ess_sys = round(ess(d$ow[d$treat == 0])),
    max_smd_unadj = round(max(abs(bt$Diff.Un), na.rm = TRUE), 3),
    max_smd_adj   = round(max(abs(bt$Diff.Adj), na.rm = TRUE), 3),
    HR = hr, CI_lower = lo, CI_upper = hi, p = pv)

  # ATO-weighted KM data (both arms)
  for (tr in c(0, 1)) {
    sub <- d %>% filter(treat == tr)
    fit <- survfit(Surv(os_months, death) ~ 1, data = sub, weights = ow)
    km_rows[[paste(g, tr)]] <- data.frame(
      group = g, arm = if_else(tr == 1L, g, REF),
      time = fit$time, surv = fit$surv, n_risk = fit$n.risk)
  }
  cat(sprintf("%-22s vs %-12s | ESS_ref=%d/%d  max|SMD| %.3f->%.3f  HR=%.2f (%.2f-%.2f) p=%.3g\n",
              g, REF, forest_rows[[g]]$ess_sys, forest_rows[[g]]$n_sys,
              forest_rows[[g]]$max_smd_unadj, forest_rows[[g]]$max_smd_adj, hr, lo, hi, pv))
}

forest_df <- bind_rows(forest_rows)
write_csv(forest_df, file.path(OUT_DIR, "ow_forest_data.csv"))
write_csv(bind_rows(bal_rows), file.path(OUT_DIR, "ow_balance_long.csv"))
write_csv(bind_rows(km_rows),  file.path(OUT_DIR, "ow_km_data.csv"))
write_csv(bind_rows(wt_rows),  file.path(OUT_DIR, "ow_weights.csv"))

cat(sprintf("\nAll contrasts balanced: max post-OW |SMD| across 7 contrasts = %.3f\n",
            max(forest_df$max_smd_adj)))
cat("Saved ->", OUT_DIR, "\n")
