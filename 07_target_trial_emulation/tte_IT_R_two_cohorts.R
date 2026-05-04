#!/usr/bin/env Rscript
# =============================================================================
# Target Trial Emulation (TTE): IT_RULES applied to TWO cohorts
# =============================================================================
# Version    : IT_RULES_two_cohorts_v1
# Purpose    : Applies the Clone-Censor-Weight IT-Rules pipeline from
#              tte_IT_R_core_cohort_3matched.R to TWO cohorts, each with its
#              own output folder:
#                (A) cohort_3matched        — matched_06 + all HAIC+I+T_concurrent
#                                             Rules: PLR > 102.4, AFP drop < 32.5%
#                                             Add-on: immune + target (I+T)
#                (B) cohort_7group_psm02    — matched_02 + all HAIC+I_concurrent
#                                             Rules: PLR > 98.7,  AFP drop < 40%
#                                             Add-on: immune only (I)
# Run        : Rscript tte_IT_R_two_cohorts.R <data_dir>
# Framework  : Clone-Censor-Weight (CCW) with stabilized IPCW
# Display    : Internal arm name 'dynamic' is displayed as "Adaptive On Demand"
#              in user-facing strings (CSV variable tags, cat() messages).
#              Filenames (km_dynamic.csv etc.) are preserved for downstream code.
# =============================================================================

Sys.setenv(LANG = "en_US.UTF-8")
Sys.setlocale("LC_ALL", "en_US.UTF-8")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(survival)
  library(survey)
  library(cobalt)
  library(boot)
})

args_all <- commandArgs(trailingOnly = TRUE)
data_dir <- args_all[1]
if (is.na(data_dir) || data_dir == "") {
  script_path <- tryCatch(
    normalizePath(sys.frames()[[1]]$ofile),
    error = function(e) NULL
  )
  if (!is.null(script_path)) {
    data_dir <- dirname(dirname(script_path))
  } else {
    data_dir <- normalizePath(file.path(getwd(), ".."))
  }
}
data_dir <- normalizePath(data_dir, winslash = "/", mustWork = FALSE)

args_raw <- commandArgs(trailingOnly = FALSE)
fa_raw   <- args_raw[grepl("^--file=", args_raw)]
project_root <- if (length(fa_raw)) {
  normalizePath(file.path(dirname(sub("^--file=", "", fa_raw[1])), "..", ".."), winslash = "/")
} else {
  normalizePath(file.path(data_dir, ".."), winslash = "/")
}
base_out_dir <- file.path(project_root, "output", "step3_tte", "IT_RULES_R_two_cohorts")
dir.create(base_out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Shared parameters ────────────────────────────────────────────────────────
EARLY_GRACE_DAYS       <- 14
DYNAMIC_GRACE_DAYS     <- 90
AFP_POST_TRIGGER_ABS   <- 20
PIVKA_POST_TRIGGER_ABS <- 40
AFP_POST_TRIGGER_RISE  <- 1.3
RMST_TAUS              <- c(12, 18, 24, 36)
N_BOOT                 <- 500
RANDOM_SEED            <- 42
IPCW_TRUNCATION        <- 0.99
EPSILON                <- 0.01
set.seed(RANDOM_SEED)

# ── Display label mapping (internal arm key -> user-facing label) ────────────
DISPLAY_DYN   <- "Adaptive On Demand"
DISPLAY_EARLY <- "Early Combination"

# ── Cohort configurations ────────────────────────────────────────────────────
CONFIGS <- list(
  cohort_3matched = list(
    label            = "HAIC then I+T on-demand (IT_RULES_v2: +SII +PIVKA +LN)",
    out_subdir       = "cohort_3matched",
    matched_ids_csv  = "matched_ids_06_HAIC_alone_vs_HAIC_then_I+T.csv",
    concurrent_group = "HAIC+I+T_concurrent",
    treatment_mode   = "I_and_T",
    afp_trigger_pct  = -32.5,
    plr_trigger      = 102.4,
    sii_trigger      = 390.9,
    nlr_trigger      = NA_real_,
    pivka_trigger_pct = -45.6,
    use_lymph_node   = TRUE
  ),
  cohort_7group_psm02 = list(
    label            = "HAIC then I on-demand (IT_RULES_v2: PLR>98.7, AFP<40%, +NLR +PIVKA +LN)",
    out_subdir       = "cohort_7group_psm02",
    matched_ids_csv  = "matched_ids_02_HAIC_alone_vs_HAIC_then_I.csv",
    concurrent_group = "HAIC+I_concurrent",
    treatment_mode   = "I_only",
    afp_trigger_pct  = -40.0,
    plr_trigger      = 98.7,
    sii_trigger      = NA_real_,
    nlr_trigger      = 2.68,
    pivka_trigger_pct = -51.2,
    use_lymph_node   = TRUE
  )
)

# =============================================================================
# Shared helper functions
# =============================================================================

normalize_pid <- function(x) trimws(as.character(x))

#' Two-layer IT-Rules trigger (v2: AFP, PLR, SII, NLR, PIVKA, meta, lymph node).
#' Layer 1 (pre_haic, cycle >= 3): trigger if any of the activated rules fire:
#'   - afp_chg_pct > afp_trigger_pct                               (always)
#'   - plr_t > plr_trigger                                         (always)
#'   - sii_t > sii_trigger                                         (if !NA)
#'   - nlr   > nlr_trigger                                         (if !NA)
#'   - pivka_chg_pct > pivka_trigger_pct                           (if !NA)
#'   - baseline distant metastasis == 1, fires at cycle 3
#'   - baseline lymph node metastasis == 1, fires at cycle 3       (if map given)
#' Layer 2 (post_haic, untriggered): trigger if afp > AFP_POST_TRIGGER_ABS
#'                                   OR pivka > PIVKA_POST_TRIGGER_ABS.
compute_tte_it_triggers <- function(lon_df, post_df = NULL, bl_afp = NULL,
                                     afp_trigger_pct = -32.5,
                                     plr_trigger = 102.4,
                                     sii_trigger = NA_real_,
                                     nlr_trigger = NA_real_,
                                     pivka_trigger_pct = NA_real_,
                                     metastasis_map = NULL,
                                     lymph_node_map = NULL) {
  req <- c("patient_id", "haic_cycle", "days_from_start", "afp_chg_pct", "plr_t")
  miss <- setdiff(req, names(lon_df))
  if (length(miss)) stop("compute_tte_it_triggers: lon_df missing columns: ",
                         paste(miss, collapse = ", "))
  if (!is.na(sii_trigger)   && !"sii_t"          %in% names(lon_df))
    stop("compute_tte_it_triggers: sii_trigger set but lon_df missing 'sii_t'")
  if (!is.na(nlr_trigger)   && !"nlr"            %in% names(lon_df))
    stop("compute_tte_it_triggers: nlr_trigger set but lon_df missing 'nlr'")
  if (!is.na(pivka_trigger_pct) && !"pivka_chg_pct" %in% names(lon_df))
    stop("compute_tte_it_triggers: pivka_trigger_pct set but lon_df missing 'pivka_chg_pct'")

  has_post <- !is.null(post_df) && nrow(post_df) > 0
  if (has_post) {
    post_req <- c("patient_id", "days_from_start", "afp", "afp_chg_pct")
    post_miss <- setdiff(post_req, names(post_df))
    if (length(post_miss)) stop("compute_tte_it_triggers: post_df missing columns: ",
                                paste(post_miss, collapse = ", "))
  }

  all_pids <- unique(lon_df$patient_id)
  if (has_post)       all_pids <- unique(c(all_pids, post_df$patient_id))
  if (!is.null(bl_afp)) all_pids <- unique(c(all_pids, names(bl_afp)))

  rows <- vector("list", length(all_pids))
  for (i in seq_along(all_pids)) {
    pid <- all_pids[i]
    sub <- lon_df[lon_df$patient_id == pid, , drop = FALSE]
    sub <- sub[order(sub$haic_cycle), , drop = FALSE]
    trig_day  <- NA_integer_
    trig_rule <- NA_character_

    has_meta <- !is.null(metastasis_map) &&
                pid %in% names(metastasis_map) &&
                !is.na(metastasis_map[[pid]]) &&
                metastasis_map[[pid]] == 1L
    has_ln   <- !is.null(lymph_node_map) &&
                pid %in% names(lymph_node_map) &&
                !is.na(lymph_node_map[[pid]]) &&
                lymph_node_map[[pid]] == 1L

    sub_eval <- sub[sub$haic_cycle >= 3L, , drop = FALSE]
    for (j in seq_len(nrow(sub_eval))) {
      ac  <- sub_eval$afp_chg_pct[j]
      plr <- sub_eval$plr_t[j]
      sii <- if ("sii_t"         %in% names(sub_eval)) sub_eval$sii_t[j]         else NA_real_
      nlv <- if ("nlr"           %in% names(sub_eval)) sub_eval$nlr[j]           else NA_real_
      piv <- if ("pivka_chg_pct" %in% names(sub_eval)) sub_eval$pivka_chg_pct[j] else NA_real_
      rule1 <- !is.na(ac)  && ac  > afp_trigger_pct
      rule2 <- !is.na(plr) && plr > plr_trigger
      rule3 <- has_meta && sub_eval$haic_cycle[j] == 3L
      rule4 <- !is.na(sii_trigger)        && !is.na(sii) && sii > sii_trigger
      rule5 <- !is.na(nlr_trigger)        && !is.na(nlv) && nlv > nlr_trigger
      rule6 <- !is.na(pivka_trigger_pct)  && !is.na(piv) && piv > pivka_trigger_pct
      rule7 <- has_ln && sub_eval$haic_cycle[j] == 3L
      if (rule1 || rule2 || rule3 || rule4 || rule5 || rule6 || rule7) {
        trig_day <- as.integer(sub_eval$days_from_start[j])
        parts <- c(if (rule1) "afp",   if (rule2) "plr",
                   if (rule3) "meta",  if (rule4) "sii",
                   if (rule5) "nlr",   if (rule6) "pivka",
                   if (rule7) "ln")
        trig_rule <- paste0(
          "layer1_cycle", sub_eval$haic_cycle[j], "_",
          paste(parts, collapse = "+")
        )
        break
      }
    }

    if (is.na(trig_day) && has_post) {
      post_sub <- post_df[post_df$patient_id == pid, , drop = FALSE]
      post_sub <- post_sub[order(post_sub$days_from_start), , drop = FALSE]
      has_pivka_col <- "pivka" %in% names(post_sub)
      if (nrow(post_sub) >= 1L) {
        for (jp in seq_len(nrow(post_sub))) {
          cur_afp   <- post_sub$afp[jp]
          cur_pivka <- if (has_pivka_col) post_sub$pivka[jp] else NA_real_
          fire_afp   <- !is.na(cur_afp)   && cur_afp   > AFP_POST_TRIGGER_ABS
          fire_pivka <- !is.na(cur_pivka) && cur_pivka > PIVKA_POST_TRIGGER_ABS
          if (fire_afp || fire_pivka) {
            trig_day  <- as.integer(post_sub$days_from_start[jp])
            parts2 <- c(if (fire_afp) "afp", if (fire_pivka) "pivka")
            trig_rule <- paste0("layer2_post_haic_", paste(parts2, collapse = "+"))
            break
          }
        }
      }
    }

    rows[[i]] <- data.frame(
      patient_id = pid,
      trigger_day = trig_day,
      trigger_stage = if (is.na(trig_rule)) "never" else trig_rule,
      stringsAsFactors = FALSE
    )
  }
  bind_rows(rows)
}

fit_stabilized_ipcw <- function(df_arm, arm_name, covar_cols) {
  df <- df_arm
  df$uncensored <- 1L - df$censored
  p_uncens <- mean(df$uncensored)
  if (p_uncens == 1 || p_uncens == 0) {
    cat(sprintf("  %s: no artificial censoring -> IPCW=1.0\n", arm_name))
    return(rep(1.0, nrow(df)))
  }
  avail <- covar_cols[covar_cols %in% names(df)]
  for (v in avail) df[[v]][is.na(df[[v]])] <- median(df[[v]], na.rm = TRUE)
  avail <- avail[sapply(avail, function(v) {
    x <- df[[v]]; !all(is.na(x)) && length(unique(x[!is.na(x)])) > 1
  })]
  time_cuts <- unique(quantile(df$os_m, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
  use_time <- FALSE
  if (length(time_cuts) >= 2) {
    df$time_cat <- tryCatch({
      tc <- cut(df$os_m, breaks = c(-Inf, time_cuts, Inf), include.lowest = TRUE)
      as.numeric(tc)
    }, error = function(e) NULL)
    use_time <- !is.null(df$time_cat) && length(unique(df$time_cat)) > 1
  }
  if (use_time) {
    fml_denom <- as.formula(paste("uncensored ~ time_cat +", paste(avail, collapse = " + ")))
    fml_numer <- uncensored ~ time_cat
  } else {
    fml_denom <- as.formula(paste("uncensored ~", paste(avail, collapse = " + ")))
    fml_numer <- uncensored ~ 1
  }
  m_denom <- tryCatch(
    glm(fml_denom, data = df, family = binomial(link = "logit")),
    error = function(e) tryCatch({
      fml_simple <- as.formula(paste("uncensored ~", paste(avail, collapse = " + ")))
      glm(fml_simple, data = df, family = binomial(link = "logit"))
    }, error = function(e2) {
      cat(sprintf("  %s: IPCW denom model failed, using intercept-only\n", arm_name))
      glm(uncensored ~ 1, data = df, family = binomial(link = "logit"))
    })
  )
  p_denom <- predict(m_denom, type = "response")
  p_denom <- pmin(pmax(p_denom, 0.05), 0.95)
  m_numer <- tryCatch(
    glm(fml_numer, data = df, family = binomial(link = "logit")),
    error = function(e) glm(uncensored ~ 1, data = df, family = binomial(link = "logit"))
  )
  p_numer <- predict(m_numer, type = "response")
  p_numer <- pmin(pmax(p_numer, 0.05), 0.95)
  sw <- p_numer / p_denom
  p_trunc <- quantile(sw, IPCW_TRUNCATION)
  sw <- pmin(sw, p_trunc)
  cat(sprintf("  %s: n_covars=%d (+time), uncensored_rate=%.3f, IPCW [%.3f, %.3f], mean=%.3f (trunc %.0fth)\n",
              arm_name, length(avail), p_uncens, min(sw), max(sw), mean(sw), IPCW_TRUNCATION * 100))
  return(sw)
}

weighted_rmst <- function(time, event, weights, tau) {
  df_tmp <- data.frame(T = time, E = event, W = weights) %>%
    filter(T > 0) %>% arrange(T)
  uniq_t <- sort(unique(df_tmp$T))
  surv <- 1.0
  at_risk_w <- sum(df_tmp$W)
  t_out <- 0.0
  s_out <- 1.0
  for (ti in uniq_t) {
    mask_t <- df_tmp$T == ti
    d_w <- sum(df_tmp$W[mask_t & df_tmp$E == 1])
    n_w <- sum(df_tmp$W[mask_t])
    if (at_risk_w > 0) surv <- surv * (1 - d_w / at_risk_w)
    at_risk_w <- at_risk_w - n_w
    t_out <- c(t_out, ti)
    s_out <- c(s_out, surv)
  }
  mask <- t_out <= tau
  tr <- t_out[mask]; sr <- s_out[mask]
  tr <- c(tr, tau)
  n_pts <- length(tr)
  area <- 0.0
  for (j in seq_len(n_pts - 1)) area <- area + sr[j] * (tr[j + 1] - tr[j])
  list(rmst = area, t_grid = t_out, s_grid = s_out)
}

compute_evalue <- function(hr_est) {
  if (is.na(hr_est)) return(NA)
  h <- if (hr_est >= 1) hr_est else 1 / hr_est
  h + sqrt(h * (h - 1))
}

# =============================================================================
# Main per-cohort pipeline
# =============================================================================

run_cohort <- function(cfg, baseline_all, longitudinal_all,
                       psm_dir, project_root_dir, ar_df) {

  out_dir <- file.path(base_out_dir, cfg$out_subdir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  cat("\n\n")
  cat("######################################################################\n")
  cat(sprintf("## Cohort: %s\n", cfg$out_subdir))
  cat(sprintf("## %s\n", cfg$label))
  cat(sprintf("## treatment_mode=%s, AFP drop<%.1f%% (afp_chg_pct>%.1f), PLR>%.2f\n",
              cfg$treatment_mode, -cfg$afp_trigger_pct, cfg$afp_trigger_pct, cfg$plr_trigger))
  if (!is.na(cfg$sii_trigger))
    cat(sprintf("##   SII (PLT*NLR) > %.2f\n", cfg$sii_trigger))
  if (!is.na(cfg$nlr_trigger))
    cat(sprintf("##   NLR > %.2f\n", cfg$nlr_trigger))
  if (!is.na(cfg$pivka_trigger_pct))
    cat(sprintf("##   PIVKA drop<%.1f%% (pivka_chg_pct>%.1f)\n",
                -cfg$pivka_trigger_pct, cfg$pivka_trigger_pct))
  if (isTRUE(cfg$use_lymph_node))
    cat("##   Lymph-node metastasis (fires at cycle 3 if baseline LN=1)\n")
  cat(sprintf("## Output: %s\n", out_dir))
  cat("######################################################################\n")

  # ── Phase 0b: cohort selection ────────────────────────────────────────────
  cat("\n--- Phase 0b: cohort selection ---\n")

  m_path <- file.path(psm_dir, cfg$matched_ids_csv)
  if (!file.exists(m_path)) stop("Matched IDs file not found: ", m_path)
  m_df <- read.csv(m_path, stringsAsFactors = FALSE)
  m_ids <- unique(normalize_pid(m_df$patient_id))

  # Concurrent patients: prefer baseline.main_group; fall back to analysis_ready
  if ("main_group" %in% names(baseline_all)) {
    conc_ids <- unique(normalize_pid(
      baseline_all$patient_id[baseline_all$main_group == cfg$concurrent_group]
    ))
  } else {
    conc_ids <- character(0)
  }
  if (length(conc_ids) == 0 && !is.null(ar_df)) {
    conc_ids <- unique(normalize_pid(
      ar_df$patient_id[ar_df$main_group == cfg$concurrent_group]
    ))
  }

  cohort_ids <- unique(c(m_ids, conc_ids))
  cat(sprintf("  matched cohort (%s): %d\n", cfg$matched_ids_csv, length(m_ids)))
  cat(sprintf("  %-25s : %d\n", cfg$concurrent_group, length(conc_ids)))
  cat(sprintf("  overlap: %d   union: %d\n",
              length(intersect(m_ids, conc_ids)), length(cohort_ids)))

  out_ids_csv <- file.path(out_dir, sprintf("tte_cohort_%s_IT_ids.csv", cfg$out_subdir))
  write.csv(data.frame(patient_id = sort(cohort_ids)), out_ids_csv, row.names = FALSE)
  cat(sprintf("  Wrote %s\n", out_ids_csv))

  baseline     <- baseline_all     %>% filter(normalize_pid(patient_id) %in% cohort_ids)
  longitudinal <- longitudinal_all %>% filter(normalize_pid(patient_id) %in% cohort_ids)
  cat(sprintf("  baseline rows after filter: %d (unique pids: %d)\n",
              nrow(baseline), n_distinct(baseline$patient_id)))

  # ── Phase 1: Baseline feature engineering ─────────────────────────────────
  cat("\n--- Phase 1: Baseline feature engineering ---\n")
  bl <- baseline %>%
    mutate(
      os_months = os_days / 30.44,
      death_status = as.integer(death_status == 1 | tolower(death_status) == "yes"),
      immune_flag = as.integer(has_immunotherapy == 1 | tolower(has_immunotherapy) == "yes"),
      target_flag = as.integer(has_target_therapy == 1 | tolower(has_target_therapy) == "yes"),
      sex_binary = as.numeric(sex == "\u7537"),
      log_afp_bl = log1p(pmax(afp, 0)),
      log_pivka_bl = log1p(pmax(pivka, 0)),
      afp_high = as.numeric(afp > 400),
      pivka_high = as.numeric(pivka > 8000),
      pvtt_grade = case_when(
        pvtt_classification == "\u65e0" ~ 0L,
        pvtt_classification == "vp1/2" ~ 1L,
        pvtt_classification == "vp3/4" ~ 2L,
        TRUE ~ 0L
      ),
      pvtt_advanced = as.integer(pvtt_grade >= 2),
      hvtt_binary = as.numeric(hvtt == "\u6709"),
      ivc_ra_binary = as.numeric(ivc_or_ra_thrombus == "\u6709"),
      metastasis_binary = as.numeric(
        tolower(trimws(as.character(distant_metastasis))) %in% c("yes", "\u662f", "1")
      ),
      lymph_node_meta_binary = as.numeric(
        tolower(trimws(as.character(lymph_node_metastasis))) %in% c("yes", "\u662f", "1")
      ),
      tumor_count_enc = case_when(
        tumor_count_category == "\u5355\u4e2a" ~ 0,
        tumor_count_category == "2-3\u4e2a" ~ 1,
        tumor_count_category == "\u591a\u53d1" ~ 2,
        TRUE ~ 1
      ),
      tumor_large = as.numeric(tumor_max_diameter_cm > 10),
      albi_bl = coalesce(albi_score_calculated, albi_score),
      albi_grade_enc = pmin(pmax(coalesce(albi_grade_calculated, albi_grade), 1), 3),
      ascites_score_enc = case_when(
        ascites == "\u65e0" ~ 0,
        ascites == "\u5c11\u91cf" ~ 1,
        ascites == "\u4e2d-\u5927\u91cf" ~ 2,
        TRUE ~ 0
      ),
      nlr_bl = nlr,
      plr_bl = plt / ifelse(lymph == 0, NA_real_, lymph),
      neut_bl = neut,
      lymph_bl = lymph
    ) %>%
    filter(os_months > 0, !is.na(death_status))

  # ── treatment_added & actual_treatment_day depend on config ───────────────
  bl <- bl %>% mutate(
    mono_bl = suppressWarnings(as.numeric(mono)),
    piv_bl = plt * mono_bl * neut / pmax(lymph, EPSILON),
    days_to_immune = suppressWarnings(as.numeric(days_haic_to_immune)),
    days_to_target = suppressWarnings(as.numeric(days_haic_to_target))
  )
  if (cfg$treatment_mode == "I_and_T") {
    bl$treatment_added <- as.integer(bl$immune_flag == 1 & bl$target_flag == 1)
    bl$actual_treatment_day <- ifelse(
      bl$treatment_added == 1,
      pmax(coalesce(bl$days_to_immune, 9999), coalesce(bl$days_to_target, 9999)),
      9999)
  } else if (cfg$treatment_mode == "I_only") {
    bl$treatment_added <- as.integer(bl$immune_flag == 1)
    bl$actual_treatment_day <- ifelse(
      bl$treatment_added == 1,
      coalesce(bl$days_to_immune, 9999),
      9999)
  } else {
    stop("Unknown treatment_mode: ", cfg$treatment_mode)
  }

  n_eligible <- nrow(bl)
  cat(sprintf("  Eligible patients: %d\n", n_eligible))
  cat(sprintf("  treatment_added=1 (%s): %d, HAIC-only/partial: %d\n",
              cfg$treatment_mode, sum(bl$treatment_added),
              sum(bl$treatment_added == 0)))

  # ── Phase 2: Propensity score (IPTW sensitivity) ──────────────────────────
  cat("\n--- Phase 2: Propensity score (IPTW sensitivity) ---\n")
  ps_vars <- c("age", "sex_binary", "alt", "ast", "tbil", "alb", "inr", "plt",
               "creatinine", "tumor_max_diameter_cm", "tumor_large",
               "tumor_count_enc", "pvtt_grade", "hvtt_binary", "ivc_ra_binary",
               "metastasis_binary", "ascites_score_enc",
               "log_afp_bl", "afp_high", "log_pivka_bl", "pivka_high",
               "albi_bl", "albi_grade_enc", "neut_bl", "lymph_bl", "nlr_bl", "plr_bl",
               "mono_bl", "piv_bl")
  ps_vars <- ps_vars[ps_vars %in% names(bl)]
  for (v in ps_vars) bl[[v]][is.na(bl[[v]])] <- median(bl[[v]], na.rm = TRUE)

  ps_formula <- as.formula(paste("treatment_added ~", paste(ps_vars, collapse = " + ")))
  ps_model <- glm(ps_formula, data = bl, family = binomial(link = "logit"))
  bl$ps <- predict(ps_model, type = "response")
  bl$ps_clip <- pmin(pmax(bl$ps, 0.05), 0.95)
  p_t <- mean(bl$treatment_added)
  bl$sw_iptw <- ifelse(bl$treatment_added == 1, p_t / bl$ps_clip, (1 - p_t) / (1 - bl$ps_clip))
  cat(sprintf("  PS range: [%.3f, %.3f], mean SW=%.3f\n",
              min(bl$ps), max(bl$ps), mean(bl$sw_iptw)))

  # ── Phase 3: Longitudinal prep ────────────────────────────────────────────
  cat("\n--- Phase 3: Longitudinal data preparation ---\n")
  lon <- longitudinal %>%
    filter(str_detect(timepoint_type, "^pre_haic_\\d+$")) %>%
    mutate(
      haic_cycle = as.integer(str_extract(timepoint_type, "\\d+")),
      haic_date = as.Date(haic_date)
    )

  bl_merge <- bl %>% select(
    patient_id, first_haic_date,
    days_to_immune, days_to_target, treatment_added, actual_treatment_day,
    os_days, os_months, death_status, sw_iptw,
    albi_bl, nlr_bl, neut_bl, lymph_bl
  ) %>% mutate(first_haic_date = as.Date(first_haic_date))

  bl_biomarkers <- bl %>% select(patient_id) %>%
    mutate(
      afp_bl0 = bl$afp,
      nlr_bl0 = bl$nlr,
      lymph_bl0 = bl$lymph,
      pivka_bl0 = suppressWarnings(as.numeric(bl$pivka))
    )
  pre1_tmp <- lon %>% filter(haic_cycle == 1) %>%
    mutate(plr_bl0 = plt / ifelse(lymph == 0, NA_real_, lymph)) %>%
    select(patient_id, plr_bl0)
  bl_biomarkers <- bl_biomarkers %>% left_join(pre1_tmp, by = "patient_id")

  lon <- lon %>%
    inner_join(bl_merge, by = "patient_id", suffix = c("", "_bl")) %>%
    left_join(bl_biomarkers, by = "patient_id")
  for (nm in c("albi_score_calculated", "albi", "albi_score")) {
    if (!nm %in% names(lon)) lon[[nm]] <- NA_real_
  }
  lon <- lon %>% mutate(
    albi_t = suppressWarnings(as.numeric(dplyr::coalesce(
      .data[["albi_score_calculated"]], .data[["albi"]], .data[["albi_score"]]))),
    mono_t = suppressWarnings(as.numeric(mono)),
    days_from_start = as.integer(difftime(haic_date, first_haic_date, units = "days")),
    days_from_start = ifelse(is.na(days_from_start), 0L, days_from_start),
    plr_t = plt / ifelse(lymph == 0, NA_real_, lymph),
    piv_t = plt * mono_t * neut / pmax(lymph, EPSILON),
    sii_t = plt * suppressWarnings(as.numeric(nlr)),
    pivka_num = suppressWarnings(as.numeric(pivka)),
    afp_chg_pct = (afp - afp_bl0) / pmax(afp_bl0, EPSILON) * 100,
    nlr_chg_pct = (nlr - nlr_bl0) / pmax(nlr_bl0, EPSILON) * 100,
    plr_chg_pct = (plr_t - plr_bl0) / pmax(plr_bl0, EPSILON) * 100,
    lymph_chg_pct = (lymph - lymph_bl0) / pmax(lymph_bl0, EPSILON) * 100,
    pivka_chg_pct = (pivka_num - pivka_bl0) / pmax(pivka_bl0, EPSILON) * 100
  )
  cat(sprintf("  Longitudinal rows: %d, unique patients: %d\n",
              nrow(lon), n_distinct(lon$patient_id)))

  cat("  Loading post_haic data for Layer 2 AFP/PIVKA monitoring...\n")
  post_haic_df <- longitudinal %>%
    filter(timepoint_type == "post_haic") %>%
    inner_join(bl_merge %>% select(patient_id, first_haic_date),
               by = "patient_id") %>%
    left_join(bl_biomarkers %>% select(patient_id, afp_bl0), by = "patient_id") %>%
    mutate(
      haic_date = as.Date(haic_date),
      days_from_start = as.integer(difftime(haic_date, first_haic_date, units = "days")),
      days_from_start = ifelse(is.na(days_from_start), 0L, days_from_start),
      afp = suppressWarnings(as.numeric(afp)),
      pivka = suppressWarnings(as.numeric(pivka)),
      afp_chg_pct = (afp - afp_bl0) / pmax(afp_bl0, EPSILON) * 100
    ) %>% arrange(patient_id, days_from_start)
  cat(sprintf("  post_haic rows: %d, unique patients: %d\n",
              nrow(post_haic_df), n_distinct(post_haic_df$patient_id)))

  bl_afp_vec <- setNames(bl$afp, bl$patient_id)
  meta_vec   <- setNames(as.integer(bl$metastasis_binary == 1), bl$patient_id)
  ln_vec     <- if (isTRUE(cfg$use_lymph_node) &&
                    "lymph_node_meta_binary" %in% names(bl)) {
    setNames(as.integer(bl$lymph_node_meta_binary == 1), bl$patient_id)
  } else NULL
  cat(sprintf("  Baseline distant metastasis: %d / %d patients\n",
              sum(meta_vec == 1, na.rm = TRUE), length(meta_vec)))
  if (!is.null(ln_vec))
    cat(sprintf("  Baseline lymph-node metastasis: %d / %d patients\n",
                sum(ln_vec == 1, na.rm = TRUE), length(ln_vec)))

  # ── Phase 4: Trigger classification ───────────────────────────────────────
  cat("\n--- Phase 4: Classify under both strategies (IT_RULES_v2: AFP/PLR/SII/NLR/PIVKA + meta + LN) ---\n")
  trig_tbl <- compute_tte_it_triggers(
    lon, post_haic_df, bl_afp_vec,
    afp_trigger_pct   = cfg$afp_trigger_pct,
    plr_trigger       = cfg$plr_trigger,
    sii_trigger       = cfg$sii_trigger,
    nlr_trigger       = cfg$nlr_trigger,
    pivka_trigger_pct = cfg$pivka_trigger_pct,
    metastasis_map    = meta_vec,
    lymph_node_map    = ln_vec
  )

  pat_df <- bl %>%
    select(patient_id, os_months, death_status, sw_iptw, actual_treatment_day) %>%
    left_join(trig_tbl, by = "patient_id") %>%
    mutate(
      triggered = as.integer(!is.na(trigger_day)),
      trigger_stage = coalesce(trigger_stage, "never"),
      trigger_day_eff = ifelse(is.na(trigger_day), 9999L, trigger_day),
      eff_trigger_day = trigger_day_eff,
      eff_triggered = triggered
    )
  stage_counts <- table(pat_df$trigger_stage)
  cat("  Trigger stage breakdown:\n")
  for (s in sort(names(stage_counts)))
    cat(sprintf("    %-20s: %d\n", s, stage_counts[[s]]))
  cat(sprintf("  Total strategy-triggered: %d\n", sum(pat_df$eff_triggered == 1L)))
  cat(sprintf("  Never triggered:          %d\n", sum(pat_df$eff_triggered == 0L)))
  cat(sprintf("  Early add-on (<=14d):     %d\n", sum(pat_df$actual_treatment_day <= EARLY_GRACE_DAYS)))

  # ── Phase 5: CCW dataset ──────────────────────────────────────────────────
  cat("\n--- Phase 5: Clone-Censor dataset ---\n")

  early <- pat_df %>% mutate(
    arm = "early_combo",
    non_compliant = actual_treatment_day > EARLY_GRACE_DAYS,
    censored = as.integer(non_compliant),
    os_m = os_months,
    event = death_status
  )
  art_censor_m <- EARLY_GRACE_DAYS / 30.44
  censor_needed <- early$non_compliant & early$os_months > art_censor_m
  early$os_m[censor_needed] <- art_censor_m
  early$event[censor_needed] <- 0

  dyn <- pat_df %>% mutate(
    arm = "dynamic",
    censored = 0L,
    os_m = os_months,
    event = death_status
  )
  dc_trig   <- dyn$eff_triggered == 1 & dyn$actual_treatment_day > (dyn$eff_trigger_day + DYNAMIC_GRACE_DAYS)
  dc_untrig <- dyn$eff_triggered == 0 & dyn$actual_treatment_day < 9999
  dc_early  <- dyn$eff_triggered == 1 & dyn$actual_treatment_day < dyn$eff_trigger_day
  dc_mask   <- dc_trig | dc_untrig | dc_early

  cd_trig   <- (dyn$eff_trigger_day + DYNAMIC_GRACE_DAYS) / 30.44
  cd_untrig <- dyn$actual_treatment_day / 30.44
  cd_early  <- dyn$actual_treatment_day / 30.44
  cd_dyn    <- ifelse(dc_trig,  cd_trig,
               ifelse(dc_early, cd_early,
               ifelse(dc_untrig, cd_untrig, 9999)))
  cn_dyn <- dc_mask & dyn$os_months > cd_dyn
  dyn$os_m[cn_dyn]  <- cd_dyn[cn_dyn]
  dyn$event[cn_dyn] <- 0
  dyn$censored[dc_mask] <- 1L

  cat(sprintf("  Case 1 (triggered, added too late >grace):    %d censored\n", sum(dc_trig)))
  cat(sprintf("  Case 2 (untriggered, but received add-on):    %d censored\n", sum(dc_untrig)))
  cat(sprintf("  Case 3 (triggered, added EARLY before trig):  %d censored\n", sum(dc_early)))

  clone_df <- bind_rows(
    dyn   %>% select(patient_id, arm, os_m, event, sw_iptw, censored),
    early %>% select(patient_id, arm, os_m, event, sw_iptw, censored)
  ) %>% filter(os_m > 0)

  for (a in c("dynamic", "early_combo")) {
    sub <- clone_df %>% filter(arm == a)
    disp <- if (a == "dynamic") DISPLAY_DYN else DISPLAY_EARLY
    cat(sprintf("  %s [%s]: n=%d, events=%d, art_censored=%d, median_os=%.1fm\n",
                a, disp, nrow(sub), sum(sub$event),
                sum(sub$censored), median(sub$os_m)))
  }

  # ── Phase 6: IPCW weights ─────────────────────────────────────────────────
  cat("\n--- Phase 6: Stabilized IPCW weights ---\n")
  clone_df <- clone_df %>%
    left_join(bl %>% select(patient_id, all_of(ps_vars)), by = "patient_id")

  pre3_data <- lon %>%
    filter(haic_cycle == 3) %>%
    distinct(patient_id, .keep_all = TRUE) %>%
    transmute(
      patient_id,
      nlr_pre3 = nlr, neut_pre3 = neut, lymph_pre3 = lymph, afp_pre3 = afp,
      mono_pre3 = mono_t, piv_pre3 = piv_t, albi_pre3 = albi_t,
      afp_chg_pre3 = afp_chg_pct, nlr_chg_pre3 = nlr_chg_pct,
      plr_pre3 = plr_t, plr_chg_pre3 = plr_chg_pct, lymph_chg_pre3 = lymph_chg_pct
    )
  clone_df <- clone_df %>% left_join(pre3_data, by = "patient_id")

  mask_e <- clone_df$arm == "early_combo"
  clone_df$sw_ipcw <- NA_real_
  clone_df$sw_ipcw[mask_e] <- fit_stabilized_ipcw(
    clone_df[mask_e, ], "early_combo", ps_vars)

  mask_d <- clone_df$arm == "dynamic"
  dyn_ipcw_vars <- ps_vars
  clone_df$sw_ipcw[mask_d] <- fit_stabilized_ipcw(
    clone_df[mask_d, ], "dynamic", dyn_ipcw_vars)

  clone_df$sw <- clone_df$sw_ipcw
  for (a in c("dynamic", "early_combo")) {
    mask <- clone_df$arm == a
    p_trunc <- quantile(clone_df$sw[mask], IPCW_TRUNCATION)
    clone_df$sw[mask] <- pmin(clone_df$sw[mask], p_trunc)
  }
  clone_df$sw_iptw_ipcw <- clone_df$sw_iptw * clone_df$sw_ipcw
  for (a in c("dynamic", "early_combo")) {
    mask <- clone_df$arm == a
    p_trunc <- quantile(clone_df$sw_iptw_ipcw[mask], IPCW_TRUNCATION)
    clone_df$sw_iptw_ipcw[mask] <- pmin(clone_df$sw_iptw_ipcw[mask], p_trunc)
  }
  for (a in c("dynamic", "early_combo")) {
    sub <- clone_df[clone_df$arm == a, ]
    disp <- if (a == "dynamic") DISPLAY_DYN else DISPLAY_EARLY
    cat(sprintf("  %s [%s]: IPCW-only mean=%.3f [%.3f, %.3f]; IPTW*IPCW mean=%.3f\n",
                a, disp, mean(sub$sw), min(sub$sw), max(sub$sw), mean(sub$sw_iptw_ipcw)))
  }

  # ── Phase 7: Cox + RMST ───────────────────────────────────────────────────
  cat("\n--- Phase 7: Weighted survival analysis ---\n")
  clone_df$A <- as.integer(clone_df$arm == "dynamic")

  cat("  Fitting weighted Cox model (survival::coxph, robust=TRUE)...\n")
  cox_fit <- coxph(Surv(os_m, event) ~ A, data = clone_df,
                   weights = sw, robust = TRUE)
  cox_summ <- summary(cox_fit)
  hr <- exp(coef(cox_fit)["A"])
  ci <- exp(confint(cox_fit)["A", ])
  pv <- cox_summ$coefficients["A", "Pr(>|z|)"]
  cat(sprintf("  HR (%s vs %s) = %.3f (%.3f-%.3f), P=%.4f\n",
              DISPLAY_DYN, DISPLAY_EARLY, hr, ci[1], ci[2], pv))

  cat("  Testing PH assumption (cox.zph)...\n")
  ph_test <- tryCatch(cox.zph(cox_fit), error = function(e) NULL)
  ph_p <- NA_real_
  if (!is.null(ph_test)) {
    ph_p <- ph_test$table["A", "p"]
    cat(sprintf("  PH test (Schoenfeld): P = %.4f %s\n",
                ph_p,
                ifelse(ph_p < 0.05, "[WARNING: PH assumption may be violated]", "[OK]")))
    ph_test_result <- data.frame(
      variable = sprintf("A (%s vs %s)", DISPLAY_DYN, DISPLAY_EARLY),
      chisq = ph_test$table["A", "chisq"],
      df = ph_test$table["A", "df"],
      p_value = ph_p,
      interpretation = ifelse(ph_p < 0.05, "PH assumption may be violated", "PH assumption holds")
    )
    write.csv(ph_test_result, file.path(out_dir, "R_cox_zph_test.csv"), row.names = FALSE)
  } else {
    cat("  PH test could not be computed\n")
  }

  cat("  Computing weighted RMST...\n")
  dyn_sub   <- clone_df %>% filter(arm == "dynamic")
  early_sub <- clone_df %>% filter(arm == "early_combo")

  rmst_results <- list()
  for (tau in RMST_TAUS) {
    r_dyn   <- weighted_rmst(dyn_sub$os_m,   dyn_sub$event,   dyn_sub$sw,   tau)
    r_early <- weighted_rmst(early_sub$os_m, early_sub$event, early_sub$sw, tau)
    delta   <- r_dyn$rmst - r_early$rmst

    unique_pids <- unique(clone_df$patient_id)
    n_pids <- length(unique_pids)
    boot_diffs <- numeric(N_BOOT)
    cat(sprintf("  Bootstrap (tau=%dm, N=%d, with IPCW re-estimation)...\n", tau, N_BOOT))

    for (b in seq_len(N_BOOT)) {
      boot_pids <- sample(unique_pids, n_pids, replace = TRUE)
      pid_counts <- as.data.frame(table(pid = boot_pids), stringsAsFactors = FALSE)
      dup_rows <- list()
      for (rr in seq_len(nrow(pid_counts))) {
        pid_val <- pid_counts$pid[rr]
        cnt <- pid_counts$Freq[rr]
        idx <- which(clone_df$patient_id == pid_val)
        dup_rows[[rr]] <- rep(idx, cnt)
      }
      boot_idx <- unlist(dup_rows)
      boot_df <- clone_df[boot_idx, ]
      boot_df <- boot_df[!is.na(boot_df$os_m) & boot_df$os_m > 0, ]

      for (a_boot in c("dynamic", "early_combo")) {
        mask_boot <- boot_df$arm == a_boot
        boot_arm  <- boot_df[mask_boot, ]
        boot_uncens <- 1L - boot_arm$censored
        p_uncens_boot <- mean(boot_uncens)
        if (p_uncens_boot == 1 || p_uncens_boot == 0 || nrow(boot_arm) < 20) {
          boot_df$sw[mask_boot] <- 1.0; next
        }
        boot_covars <- if (a_boot == "dynamic") dyn_ipcw_vars else ps_vars
        boot_avail  <- boot_covars[boot_covars %in% names(boot_arm)]
        boot_arm$uncensored <- boot_uncens
        for (vv in boot_avail)
          boot_arm[[vv]][is.na(boot_arm[[vv]])] <- median(boot_arm[[vv]], na.rm = TRUE)
        fml_boot <- tryCatch(
          as.formula(paste("uncensored ~", paste(boot_avail, collapse = " + "))),
          error = function(e) NULL)
        if (is.null(fml_boot)) { boot_df$sw[mask_boot] <- 1.0; next }
        m_boot <- tryCatch(
          glm(fml_boot, data = boot_arm, family = binomial(link = "logit")),
          error = function(e) NULL, warning = function(w) NULL)
        if (is.null(m_boot)) { boot_df$sw[mask_boot] <- 1.0; next }
        p_d_boot <- predict(m_boot, type = "response")
        p_d_boot <- pmin(pmax(p_d_boot, 0.05), 0.95)
        sw_boot <- p_uncens_boot / p_d_boot
        pt_boot <- quantile(sw_boot, IPCW_TRUNCATION)
        sw_boot <- pmin(sw_boot, pt_boot)
        boot_df$sw[mask_boot] <- sw_boot
      }

      bd <- boot_df[boot_df$arm == "dynamic", ]
      be <- boot_df[boot_df$arm == "early_combo", ]
      if (nrow(bd) < 10 || nrow(be) < 10) { boot_diffs[b] <- NA; next }
      rd <- tryCatch(weighted_rmst(bd$os_m, bd$event, bd$sw, tau)$rmst, error = function(e) NA)
      re <- tryCatch(weighted_rmst(be$os_m, be$event, be$sw, tau)$rmst, error = function(e) NA)
      boot_diffs[b] <- rd - re
    }
    boot_diffs <- boot_diffs[!is.na(boot_diffs)]
    ci_lo <- quantile(boot_diffs, 0.025)
    ci_hi <- quantile(boot_diffs, 0.975)
    p_boot <- 2 * min(mean(boot_diffs <= 0), mean(boot_diffs >= 0))

    rmst_results[[as.character(tau)]] <- list(
      rmst_dyn = r_dyn$rmst, rmst_early = r_early$rmst,
      delta = delta, ci_lo = ci_lo, ci_hi = ci_hi, p = p_boot
    )
    cat(sprintf("  tau=%dm: %s=%.2fm, %s=%.2fm, delta=%+.2fm (95%%CI %+.2f,%+.2f), P=%.4f\n",
                tau, DISPLAY_DYN, r_dyn$rmst,
                DISPLAY_EARLY, r_early$rmst, delta, ci_lo, ci_hi, p_boot))
  }

  km_dyn   <- weighted_rmst(dyn_sub$os_m,   dyn_sub$event,   dyn_sub$sw,   max(RMST_TAUS) + 6)
  km_early <- weighted_rmst(early_sub$os_m, early_sub$event, early_sub$sw, max(RMST_TAUS) + 6)
  write.csv(data.frame(time = km_dyn$t_grid, surv = km_dyn$s_grid),
            file.path(out_dir, "km_dynamic.csv"), row.names = FALSE)
  write.csv(data.frame(time = km_early$t_grid, surv = km_early$s_grid),
            file.path(out_dir, "km_early_combo.csv"), row.names = FALSE)

  # ── Phase 8: Diagnostics ──────────────────────────────────────────────────
  cat("\n--- Phase 8: Diagnostics ---\n")
  treated <- bl$treatment_added == 1
  smd_before <- numeric(length(ps_vars))
  smd_after  <- numeric(length(ps_vars))
  for (i in seq_along(ps_vars)) {
    v <- ps_vars[i]; x <- bl[[v]]; tr <- as.integer(treated)
    smd_before[i] <- tryCatch({
      col_w_smd(x, treat = tr, weights = NULL, std = TRUE, abs = FALSE)
    }, error = function(e) {
      pool_sd <- sqrt((var(x[treated]) + var(x[!treated])) / 2)
      if (pool_sd < 1e-10) 0 else (mean(x[treated]) - mean(x[!treated])) / pool_sd
    })
    smd_after[i] <- tryCatch({
      col_w_smd(x, treat = tr, weights = bl$sw_iptw, std = TRUE, abs = FALSE)
    }, error = function(e) {
      w_t <- bl$sw_iptw[treated]; w_c <- bl$sw_iptw[!treated]
      x_t <- x[treated]; x_c <- x[!treated]
      wm_t <- weighted.mean(x_t, w_t); wm_c <- weighted.mean(x_c, w_c)
      wv_t <- sum(w_t * (x_t - wm_t)^2) / sum(w_t)
      wv_c <- sum(w_c * (x_c - wm_c)^2) / sum(w_c)
      pool_sd_w <- sqrt((wv_t + wv_c) / 2)
      if (pool_sd_w < 1e-10) 0 else (wm_t - wm_c) / pool_sd_w
    })
  }
  smd_table <- data.frame(
    variable = ps_vars,
    SMD_unadjusted = round(smd_before, 4),
    SMD_after_IPTW = round(smd_after, 4)
  )
  write.csv(smd_table, file.path(out_dir, "R_eTable5_SMD.csv"), row.names = FALSE)

  ess_dyn   <- sum(dyn_sub$sw)^2   / sum(dyn_sub$sw^2)
  ess_early <- sum(early_sub$sw)^2 / sum(early_sub$sw^2)
  cat(sprintf("  ESS: %s=%.0f, %s=%.0f\n",
              DISPLAY_DYN, ess_dyn, DISPLAY_EARLY, ess_early))
  cat(sprintf("  |SMD|>0.1 before: %d, after IPTW: %d\n",
              sum(abs(smd_before) > 0.1), sum(abs(smd_after) > 0.1)))

  # ── Phase 8b: IPCW stability & ESS diagnostics ────────────────────────────
  cat("\n--- Phase 8b: IPCW stability & ESS diagnostics ---\n")
  ipcw_summ <- clone_df %>%
    group_by(arm) %>%
    summarise(
      n          = n(),
      mean_w     = mean(sw, na.rm = TRUE),
      sd_w       = sd(sw,   na.rm = TRUE),
      cv         = sd_w / mean_w,
      min_w      = min(sw,  na.rm = TRUE),
      median_w   = median(sw, na.rm = TRUE),
      p95        = as.numeric(quantile(sw, 0.95, na.rm = TRUE)),
      p99        = as.numeric(quantile(sw, 0.99, na.rm = TRUE)),
      max_w      = max(sw,  na.rm = TRUE),
      pct_gt_10  = mean(sw > 10,  na.rm = TRUE) * 100,
      pct_gt_20  = mean(sw > 20,  na.rm = TRUE) * 100,
      ESS        = sum(sw, na.rm = TRUE)^2 / sum(sw^2, na.rm = TRUE),
      ESS_ratio  = ESS / n,
      uncens_rate = mean(censored == 0),
      .groups    = "drop"
    )
  ess_status <- function(r) ifelse(r > 0.80, "good",
                              ifelse(r > 0.50, "acceptable",
                                ifelse(r > 0.30, "WARNING", "SEVERE")))
  ipcw_summ$ESS_status <- ess_status(ipcw_summ$ESS_ratio)
  write.csv(ipcw_summ, file.path(out_dir, "R_IPCW_diagnostics.csv"), row.names = FALSE)

  wt_quants <- clone_df %>%
    group_by(arm) %>%
    summarise(
      p01 = as.numeric(quantile(sw, 0.01, na.rm = TRUE)),
      p05 = as.numeric(quantile(sw, 0.05, na.rm = TRUE)),
      p25 = as.numeric(quantile(sw, 0.25, na.rm = TRUE)),
      p50 = as.numeric(quantile(sw, 0.50, na.rm = TRUE)),
      p75 = as.numeric(quantile(sw, 0.75, na.rm = TRUE)),
      p95 = as.numeric(quantile(sw, 0.95, na.rm = TRUE)),
      p99 = as.numeric(quantile(sw, 0.99, na.rm = TRUE)),
      p995= as.numeric(quantile(sw, 0.995, na.rm = TRUE)),
      .groups = "drop"
    )
  write.csv(wt_quants, file.path(out_dir, "R_IPCW_weight_quantiles.csv"), row.names = FALSE)

  time_grid <- seq(0, max(clone_df$os_m, na.rm = TRUE), length.out = 40)
  ess_time_rows <- list()
  for (a in unique(clone_df$arm)) {
    for (t in time_grid) {
      sub_w <- clone_df$sw[clone_df$arm == a & clone_df$os_m >= t]
      sub_w <- sub_w[!is.na(sub_w)]
      ess_t <- if (length(sub_w) > 0) sum(sub_w)^2 / sum(sub_w^2) else NA_real_
      n_t   <- length(sub_w)
      ess_time_rows[[length(ess_time_rows) + 1]] <- data.frame(
        arm = a, t = t, n_at_risk = n_t, ESS_at_t = ess_t,
        ESS_ratio_at_t = ifelse(n_t > 0, ess_t / n_t, NA_real_)
      )
    }
  }
  ess_time_df <- bind_rows(ess_time_rows)
  write.csv(ess_time_df, file.path(out_dir, "R_IPCW_time_varying_ESS.csv"), row.names = FALSE)

  cens_curve_rows <- list()
  for (a in unique(clone_df$arm)) {
    sub_a <- clone_df[clone_df$arm == a, ]; n_a <- nrow(sub_a)
    for (t in time_grid) {
      cens_by_t <- sum(sub_a$censored == 1 & sub_a$os_m <= t)
      cens_curve_rows[[length(cens_curve_rows) + 1]] <- data.frame(
        arm = a, t = t, n_censored_by_t = cens_by_t,
        cum_cens_prop = cens_by_t / n_a
      )
    }
  }
  write.csv(bind_rows(cens_curve_rows), file.path(out_dir, "R_IPCW_censoring_curve.csv"), row.names = FALSE)

  red_rows <- list()
  for (i in seq_len(nrow(ipcw_summ))) {
    r <- ipcw_summ[i, ]; flags <- character(0)
    if (abs(r$mean_w - 1.0) > 0.10) flags <- c(flags, sprintf("mean(sw)=%.2f deviates from 1.0", r$mean_w))
    if (r$max_w > 50)               flags <- c(flags, sprintf("max=%.1f (>50)", r$max_w))
    if (r$p99 > 20)                 flags <- c(flags, sprintf("p99=%.1f (>20)", r$p99))
    if (r$cv > 1.0)                 flags <- c(flags, sprintf("CV=%.2f (>1.0)", r$cv))
    if (r$pct_gt_10 > 1)            flags <- c(flags, sprintf("%.2f%% weights >10", r$pct_gt_10))
    if (r$ESS_ratio < 0.30)         flags <- c(flags, sprintf("ESS/N=%.2f SEVERE", r$ESS_ratio))
    else if (r$ESS_ratio < 0.50)    flags <- c(flags, sprintf("ESS/N=%.2f WARNING", r$ESS_ratio))
    red_rows[[i]] <- data.frame(
      arm = r$arm, n_flags = length(flags),
      flags = if (length(flags)) paste(flags, collapse = " | ") else "None",
      stringsAsFactors = FALSE
    )
  }
  write.csv(bind_rows(red_rows), file.path(out_dir, "R_IPCW_red_flags.csv"), row.names = FALSE)

  # ── Phase 9: Sensitivity analyses ─────────────────────────────────────────
  cat("\n--- Phase 9: Sensitivity analyses ---\n")

  run_sens_cox <- function(clone_data, trunc_pct, label) {
    cd <- clone_data
    cd$sw_s <- cd$sw_ipcw
    for (a in c("dynamic", "early_combo")) {
      mask <- cd$arm == a
      p_t <- quantile(cd$sw_s[mask], trunc_pct)
      cd$sw_s[mask] <- pmin(cd$sw_s[mask], p_t)
    }
    cd$A <- as.integer(cd$arm == "dynamic")
    cd <- cd[cd$os_m > 0, ]
    fit <- tryCatch(
      coxph(Surv(os_m, event) ~ A, data = cd, weights = sw_s, robust = TRUE),
      error = function(e) NULL)
    if (is.null(fit)) return(data.frame(analysis=label, HR=NA, CI_lo=NA, CI_hi=NA, P=NA))
    s <- summary(fit)
    data.frame(
      analysis = label,
      HR = round(exp(coef(fit)["A"]), 4),
      CI_lo = round(exp(confint(fit)["A", 1]), 4),
      CI_hi = round(exp(confint(fit)["A", 2]), 4),
      P = round(s$coefficients["A", "Pr(>|z|)"], 4)
    )
  }

  run_grace_sens <- function(bl_data, lon_data, post_data, bl_afp_v,
                              ps_vars_list, grace_d, label,
                              meta_v = NULL, ln_v = NULL) {
    keep_cols <- intersect(
      c("patient_id", "haic_cycle", "days_from_start", "afp_chg_pct",
        "plr_t", "sii_t", "nlr", "pivka_chg_pct"),
      names(lon_data)
    )
    lon_slim <- lon_data[, keep_cols, drop = FALSE]
    post_slim <- if (!is.null(post_data) && nrow(post_data) > 0) {
      post_keep <- intersect(
        c("patient_id", "days_from_start", "afp", "afp_chg_pct", "pivka"),
        names(post_data)
      )
      post_data[, post_keep, drop = FALSE]
    } else NULL
    trig_t <- compute_tte_it_triggers(
      lon_slim, post_slim, bl_afp_v,
      afp_trigger_pct   = cfg$afp_trigger_pct,
      plr_trigger       = cfg$plr_trigger,
      sii_trigger       = cfg$sii_trigger,
      nlr_trigger       = cfg$nlr_trigger,
      pivka_trigger_pct = cfg$pivka_trigger_pct,
      metastasis_map    = meta_v,
      lymph_node_map    = ln_v
    )
    pat_tmp <- bl_data %>%
      select(patient_id, os_months, death_status, sw_iptw, actual_treatment_day) %>%
      left_join(trig_t, by = "patient_id") %>%
      mutate(
        triggered = as.integer(!is.na(trigger_day)),
        trigger_stage = coalesce(trigger_stage, "never"),
        eff_trigger_day = ifelse(is.na(trigger_day), 9999L, trigger_day),
        eff_triggered = triggered
      )

    early_s <- pat_tmp %>% mutate(
      arm = "early_combo",
      nc = actual_treatment_day > EARLY_GRACE_DAYS,
      censored = as.integer(nc), os_m = os_months, event = death_status)
    cm <- EARLY_GRACE_DAYS / 30.44
    cn <- early_s$nc & early_s$os_months > cm
    early_s$os_m[cn] <- cm; early_s$event[cn] <- 0

    dyn_s <- pat_tmp %>% mutate(
      arm = "dynamic", censored = 0L, os_m = os_months, event = death_status)
    dc1 <- dyn_s$eff_triggered == 1 & dyn_s$actual_treatment_day > (dyn_s$eff_trigger_day + grace_d)
    dc3 <- dyn_s$eff_triggered == 0 & dyn_s$actual_treatment_day < 9999
    dc2 <- dyn_s$eff_triggered == 1 & dyn_s$actual_treatment_day < dyn_s$eff_trigger_day
    dc  <- dc1 | dc3 | dc2
    cd1 <- (dyn_s$eff_trigger_day + grace_d) / 30.44
    cd3 <- dyn_s$actual_treatment_day / 30.44
    cd2 <- dyn_s$actual_treatment_day / 30.44
    cd_d <- ifelse(dc1, cd1, ifelse(dc2, cd2, ifelse(dc3, cd3, 9999)))
    cn_d <- dc & dyn_s$os_months > cd_d
    dyn_s$os_m[cn_d] <- cd_d[cn_d]; dyn_s$event[cn_d] <- 0; dyn_s$censored[dc] <- 1L

    cd_all <- bind_rows(
      dyn_s   %>% select(patient_id, arm, os_m, event, sw_iptw, censored),
      early_s %>% select(patient_id, arm, os_m, event, sw_iptw, censored)
    ) %>% filter(os_m > 0)

    cd_all <- cd_all %>% left_join(bl_data %>% select(patient_id, all_of(ps_vars_list)), by = "patient_id")
    cd_all$sw_ipcw_s <- NA_real_
    for (a in c("dynamic", "early_combo")) {
      mask <- cd_all$arm == a
      cd_all$sw_ipcw_s[mask] <- fit_stabilized_ipcw(cd_all[mask, ], a, ps_vars_list)
    }
    cd_all$sw_s <- cd_all$sw_ipcw_s
    for (a in c("dynamic", "early_combo")) {
      mask <- cd_all$arm == a
      pt <- quantile(cd_all$sw_s[mask], IPCW_TRUNCATION, na.rm = TRUE)
      cd_all$sw_s[mask] <- pmin(cd_all$sw_s[mask], pt)
    }
    cd_all$A <- as.integer(cd_all$arm == "dynamic")
    cd_all <- cd_all[cd_all$os_m > 0, ]
    fit <- tryCatch(
      coxph(Surv(os_m, event) ~ A, data = cd_all, weights = sw_s, robust = TRUE),
      error = function(e) NULL)
    if (is.null(fit)) return(data.frame(analysis=label, HR=NA, CI_lo=NA, CI_hi=NA, P=NA))
    s <- summary(fit)
    data.frame(
      analysis = label,
      HR = round(exp(coef(fit)["A"]), 4),
      CI_lo = round(exp(confint(fit)["A", 1]), 4),
      CI_hi = round(exp(confint(fit)["A", 2]), 4),
      P = round(s$coefficients["A", "Pr(>|z|)"], 4)
    )
  }

  sens_rows <- list()
  e_val    <- compute_evalue(hr)
  e_val_lo <- if (!is.na(ci[1]) && ci[1] > 1) compute_evalue(ci[1]) else 1.0
  sens_rows[[1]] <- data.frame(
    analysis = sprintf("Primary (trunc=%.0fth pctl)", IPCW_TRUNCATION * 100),
    HR = round(hr, 4), CI_lo = round(ci[1], 4), CI_hi = round(ci[2], 4),
    P = round(pv, 4), E_value = round(e_val, 2), E_value_CI = round(e_val_lo, 2))
  cat("  Truncation 95th...\n")
  sens_rows[[2]] <- run_sens_cox(clone_df, 0.95,  "Truncation: 95th pctl")
  cat("  Truncation 99.5th...\n")
  sens_rows[[3]] <- run_sens_cox(clone_df, 0.995, "Truncation: 99.5th pctl")

  cat("  Unweighted...\n")
  cox_uw <- coxph(Surv(os_m, event) ~ A, data = clone_df[clone_df$os_m > 0, ])
  s_uw <- summary(cox_uw)
  sens_rows[[4]] <- data.frame(
    analysis = "Unweighted (no IPCW)",
    HR = round(exp(coef(cox_uw)["A"]), 4),
    CI_lo = round(exp(confint(cox_uw)["A", 1]), 4),
    CI_hi = round(exp(confint(cox_uw)["A", 2]), 4),
    P = round(s_uw$coefficients["A", "Pr(>|z|)"], 4),
    E_value = NA, E_value_CI = NA)

  cat("  IPTW x IPCW...\n")
  cox_ii <- coxph(Surv(os_m, event) ~ A, data = clone_df[clone_df$os_m > 0, ],
                  weights = sw_iptw_ipcw, robust = TRUE)
  s_ii <- summary(cox_ii)
  sens_rows[[5]] <- data.frame(
    analysis = "IPTW x IPCW (non-standard for CCW)",
    HR = round(exp(coef(cox_ii)["A"]), 4),
    CI_lo = round(exp(confint(cox_ii)["A", 1]), 4),
    CI_hi = round(exp(confint(cox_ii)["A", 2]), 4),
    P = round(s_ii$coefficients["A", "Pr(>|z|)"], 4),
    E_value = NA, E_value_CI = NA)

  idx <- 6
  for (gd_info in list(c(90, "90 days"), c(120, "120 days"), c(240, "240 days"))) {
    gd <- as.integer(gd_info[1]); glabel <- gd_info[2]
    cat(sprintf("  Grace period: %s...\n", glabel))
    r <- run_grace_sens(bl, lon, post_haic_df, bl_afp_vec, ps_vars, gd,
                        paste("Grace period:", glabel),
                        meta_v = meta_vec, ln_v = ln_vec)
    r$E_value <- NA; r$E_value_CI <- NA
    sens_rows[[idx]] <- r; idx <- idx + 1
  }

  sens_df <- bind_rows(sens_rows)
  write.csv(sens_df, file.path(out_dir, "R_Table4_sensitivity.csv"), row.names = FALSE)
  cat("  Sensitivity results:\n"); print(sens_df)

  # ── Phase 10: Save main results, RMST, clone, risk, baseline ──────────────
  cat("\n--- Phase 10: Saving results for visualization ---\n")

  main_res <- data.frame(
    HR = hr, HR_lo = ci[1], HR_hi = ci[2], HR_p = pv,
    E_value = e_val, E_value_CI = e_val_lo,
    n_eligible = n_eligible,
    n_dyn = nrow(dyn_sub), n_early = nrow(early_sub),
    events_dyn = sum(dyn_sub$event), events_early = sum(early_sub$event),
    censored_dyn = sum(dyn_sub$censored), censored_early = sum(early_sub$censored),
    median_dyn = median(dyn_sub$os_m), median_early = median(early_sub$os_m),
    ess_dyn = ess_dyn, ess_early = ess_early,
    ph_test_p = ph_p
  )
  write.csv(main_res, file.path(out_dir, "R_main_results.csv"), row.names = FALSE)

  rmst_df <- do.call(rbind, lapply(names(rmst_results), function(tau) {
    r <- rmst_results[[tau]]
    data.frame(tau = as.integer(tau),
               rmst_dyn = round(r$rmst_dyn, 3), rmst_early = round(r$rmst_early, 3),
               delta = round(r$delta, 3), ci_lo = round(r$ci_lo, 3),
               ci_hi = round(r$ci_hi, 3), p = round(r$p, 4))
  }))
  write.csv(rmst_df, file.path(out_dir, "R_Table3_RMST.csv"), row.names = FALSE)

  write.csv(
    clone_df %>% select(patient_id, arm, os_m, event, sw, sw_ipcw, sw_iptw, sw_iptw_ipcw, censored),
    file.path(out_dir, "R_clone_dataset.csv"), row.names = FALSE)

  write.csv(pat_df %>% select(patient_id, trigger_day, trigger_stage, eff_triggered),
            file.path(out_dir, "R_trigger_table.csv"), row.names = FALSE)

  risk_data <- list()
  for (a in c("dynamic", "early_combo")) {
    sub <- clone_df[clone_df$arm == a, ]
    risk_times <- c(0, 6, 12, 18, 24, 30, 36)
    counts <- sapply(risk_times, function(tp) sum(sub$os_m >= tp))
    risk_data[[a]] <- data.frame(time = risk_times, n_at_risk = counts, arm = a)
  }
  write.csv(bind_rows(risk_data), file.path(out_dir, "R_risk_table.csv"), row.names = FALSE)

  bl_table_vars <- list(
    c("age", "Age (years)", "continuous"),
    c("sex_binary", "Male sex", "binary"),
    c("alt", "ALT (U/L)", "continuous"),
    c("ast", "AST (U/L)", "continuous"),
    c("tbil", "Total bilirubin", "continuous"),
    c("alb", "Albumin (g/L)", "continuous"),
    c("inr", "INR", "continuous"),
    c("plt", "Platelet count", "continuous"),
    c("creatinine", "Creatinine", "continuous"),
    c("tumor_max_diameter_cm", "Max tumor size (cm)", "continuous"),
    c("tumor_large", "Tumor >10 cm", "binary"),
    c("pvtt_grade", "PVTT grade (0-2)", "continuous"),
    c("pvtt_advanced", "PVTT Vp3/4", "binary"),
    c("hvtt_binary", "HVTT present", "binary"),
    c("ivc_ra_binary", "IVC/RA thrombus", "binary"),
    c("metastasis_binary", "Extrahepatic metastasis", "binary"),
    c("lymph_node_meta_binary", "Lymph node metastasis", "binary"),
    c("ascites_score_enc", "Ascites score (0-2)", "continuous"),
    c("log_afp_bl", "log(AFP+1)", "continuous"),
    c("afp_high", "AFP >400 ng/mL", "binary"),
    c("log_pivka_bl", "log(PIVKA+1)", "continuous"),
    c("pivka_high", "PIVKA >8000", "binary"),
    c("albi_bl", "ALBI score", "continuous"),
    c("albi_grade_enc", "ALBI grade", "continuous"),
    c("neut_bl", "Neutrophil count", "continuous"),
    c("lymph_bl", "Lymphocyte count", "continuous"),
    c("nlr_bl", "NLR", "continuous"),
    c("plr_bl", "PLR", "continuous"),
    c("piv_bl", "PIV (baseline)", "continuous"),
    c("mono_bl", "Monocyte count (baseline)", "continuous")
  )
  bl_rows <- lapply(bl_table_vars, function(info) {
    v <- info[1]; label <- info[2]; vtype <- info[3]
    if (!v %in% names(bl)) return(NULL)
    x_t <- bl[[v]][treated]; x_c <- bl[[v]][!treated]
    if (vtype == "continuous") {
      t_str <- sprintf("%.1f +/- %.1f", mean(x_t, na.rm = TRUE), sd(x_t, na.rm = TRUE))
      c_str <- sprintf("%.1f +/- %.1f", mean(x_c, na.rm = TRUE), sd(x_c, na.rm = TRUE))
    } else {
      t_str <- sprintf("%.0f (%.1f%%)", sum(x_t, na.rm = TRUE), mean(x_t, na.rm = TRUE) * 100)
      c_str <- sprintf("%.0f (%.1f%%)", sum(x_c, na.rm = TRUE), mean(x_c, na.rm = TRUE) * 100)
    }
    smd_b <- smd_before[match(v, ps_vars)]
    smd_a <- smd_after[match(v, ps_vars)]
    data.frame(Variable = label,
               treatment_added = t_str, HAIC_only = c_str,
               SMD_unadjusted = round(abs(smd_b), 3),
               SMD_after_IPTW = round(abs(smd_a), 3))
  })
  write.csv(bind_rows(bl_rows), file.path(out_dir, "R_Table2_baseline.csv"), row.names = FALSE)

  cat(sprintf("\n[%s] DONE — HR=%.3f (%.3f-%.3f), P=%.4f, ESS(%s)=%.0f, ESS(%s)=%.0f\n",
              cfg$out_subdir, hr, ci[1], ci[2], pv,
              DISPLAY_DYN, ess_dyn, DISPLAY_EARLY, ess_early))
}

# =============================================================================
# Driver: load data once, run both cohorts
# =============================================================================

cat("======================================================================\n")
cat("TTE IT-Rules on TWO cohorts (Adaptive On Demand vs Early Combination)\n")
cat("======================================================================\n")
cat(sprintf("data_dir     : %s\n", data_dir))
cat(sprintf("project_root : %s\n", project_root))
cat(sprintf("base output  : %s\n", base_out_dir))

baseline_all     <- read.csv(file.path(data_dir, "HAIC_NO_TACE_4_TIDY_baseline.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
longitudinal_all <- read.csv(file.path(data_dir, "HAIC_NO_TACE_4_TIDY_longitudinal.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
if (!"days_haic_to_immune_y" %in% names(baseline_all) &&
     "days_haic_to_immune" %in% names(baseline_all)) {
  baseline_all$days_haic_to_immune_y <- baseline_all$days_haic_to_immune
}
if (!"albi_score_calculated" %in% names(baseline_all))
  baseline_all$albi_score_calculated <- baseline_all$albi_score
if (!"albi_grade_calculated" %in% names(baseline_all))
  baseline_all$albi_grade_calculated <- baseline_all$albi_grade
cat(sprintf("Baseline: %d patients, %d cols\n", nrow(baseline_all), ncol(baseline_all)))
cat(sprintf("Longitudinal: %d rows, %d cols\n", nrow(longitudinal_all), ncol(longitudinal_all)))

psm_dir <- file.path(
  normalizePath(file.path(data_dir, "..", ".."), winslash = "/", mustWork = TRUE),
  "HAIC_NO_TACE_4_TIDY", "update_group_7", "results", "psm_balance_tables_complete"
)
if (!dir.exists(psm_dir)) {
  # Fallback: try HAIC_Immunotherapy_Decision_TTE/HAIC_NO_TACE_4_TIDY/... (sibling to data_dir)
  psm_dir_alt <- file.path(data_dir, "HAIC_NO_TACE_4_TIDY", "update_group_7",
                           "results", "psm_balance_tables_complete")
  if (dir.exists(psm_dir_alt)) psm_dir <- psm_dir_alt
}
cat(sprintf("PSM dir      : %s\n", psm_dir))

ar_path <- file.path(
  normalizePath(file.path(data_dir, "..", ".."), winslash = "/", mustWork = TRUE),
  "HAIC_NO_TACE_4_TIDY", "update_group_7", "data", "analysis_ready.csv"
)
if (!file.exists(ar_path)) {
  ar_path_alt <- file.path(data_dir, "HAIC_NO_TACE_4_TIDY", "update_group_7",
                           "data", "analysis_ready.csv")
  if (file.exists(ar_path_alt)) ar_path <- ar_path_alt
}
ar_df <- if (file.exists(ar_path)) {
  cat(sprintf("analysis_ready: %s\n", ar_path))
  read.csv(ar_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  cat("analysis_ready.csv not found — will rely on baseline.main_group\n")
  NULL
}

for (cfg_name in names(CONFIGS)) {
  cfg <- CONFIGS[[cfg_name]]
  set.seed(RANDOM_SEED)
  run_cohort(cfg, baseline_all, longitudinal_all, psm_dir, project_root, ar_df)
}

cat("\n======================================================================\n")
cat("All cohorts complete. Results under:\n")
cat(sprintf("  %s\n", base_out_dir))
cat("======================================================================\n")
