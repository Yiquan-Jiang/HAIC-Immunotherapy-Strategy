#!/usr/bin/env Rscript
# =============================================================================
# Generate TTE Patient Flow Draw.io Diagram — IT_RULES_R_two_cohorts (parameterized)
# =============================================================================
# This script generates the Dynamic vs Early-Combo CCW flow diagram for
# the two cohorts produced by tte_IT_R_two_cohorts.R:
#   - cohort_7group_psm02   (HAIC then I on-demand)
#   - cohort_3matched       (HAIC then I+T on-demand)
#
# It is based on the two user-provided templates
#   scripts/tte_core/generate_tte_flow_drawio.R
#   scripts/tte_core/generate_tte_flow_drawio_IT_cohort3matched.R
# and uses the same vtx()/edg()/S$* style vocabulary, but the numbers are
# pulled directly from the R result CSVs, so the diagram auto-updates when
# the analysis re-runs.
#
# Usage:
#   Rscript generate_tte_flow_drawio_two_cohorts.R <cohort>
#   cohort ∈ {cohort_7group_psm02, cohort_3matched, all}
#
# Data sources (per cohort folder under
#   output/step3_tte/IT_RULES_R_two_cohorts/<cohort>/):
#   - R_main_results.csv
#   - R_Table3_RMST.csv
#   - R_trigger_table.csv          patient_id, trigger_day, trigger_stage, eff_triggered
#   - R_clone_dataset.csv          patient_id, arm, os_m, event, sw, censored
#   - R_cox_zph_test.csv           PH-assumption test
#   - pathway_classification_table.csv  actual_day, trigger_day, arm_allegiance
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# -----------------------------------------------------------------------------
# 1. Resolve paths and arguments
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Resolve script directory across both Rscript and interactive modes
resolve_script_dir <- function() {
  args_raw <- commandArgs(trailingOnly = FALSE)
  fa <- args_raw[grepl("^--file=", args_raw)]
  if (length(fa)) return(dirname(normalizePath(sub("^--file=", "", fa[1]))))
  f <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(f) && nzchar(f)) return(dirname(normalizePath(f)))
  getwd()
}
script_dir <- resolve_script_dir()
project_root <- normalizePath(file.path(script_dir, "../.."),
                              winslash = "/", mustWork = FALSE)

COHORTS_ALL <- c("cohort_7group_psm02", "cohort_3matched")
if (length(args) >= 1 && nzchar(args[1])) {
  sel <- args[1]
  cohorts <- if (sel == "all") COHORTS_ALL else sel
} else {
  cohorts <- COHORTS_ALL
}

cohorts <- intersect(cohorts, c(COHORTS_ALL, "cohort_7group_psm02", "cohort_3matched"))
if (!length(cohorts)) stop("No valid cohort selected. Valid: ",
                           paste(COHORTS_ALL, collapse = ", "), ", or 'all'.")

# Cohort-specific display labels
COHORT_META <- list(
  cohort_7group_psm02 = list(
    title = "HAIC then I On-Demand",
    rule_ab = "I (immune checkpoint inhibitor)",
    early_abbr = "I",
    rule_text = paste0(
      "<b>Trigger Rules</b> (evaluated at each pre-HAIC visit, cycle \u2265 3)<br>",
      "<font style=\"font-size:10px\">Rule 1: AFP change from baseline &gt; -46.7% (drop &lt; 46.7%)<br>",
      "Rule 2: NLR / PLR / PIVKA / LN / Meta flags<br>",
      "First visit meeting any rule \u2192 Trigger</font>"),
    cohort_note = paste0(
      "<b>Cohort:</b> HAIC alone + HAIC_then_I matched (PSM 1:1) ",
      "+ all HAIC+I_concurrent (unmatched)")
  ),
  cohort_3matched = list(
    title = "HAIC then I+T On-Demand",
    rule_ab = "I+T (immune + antiangiogenic therapy)",
    early_abbr = "I+T",
    rule_text = paste0(
      "<b>Trigger Rules</b> (evaluated at each pre-HAIC visit, cycle \u2265 3)<br>",
      "<font style=\"font-size:10px\">Rule 1: AFP change from baseline &gt; -46.7% (drop &lt; 46.7%)<br>",
      "Rule 2: SII / PLR / PIVKA / LN / Meta flags<br>",
      "First visit meeting any rule \u2192 Trigger</font>"),
    cohort_note = paste0(
      "<b>Cohort:</b> matched_06 (HAIC_alone vs HAIC_then_I+T, PSM 1:1) ",
      "+ all HAIC+I+T_concurrent (unmatched)")
  )
)

# -----------------------------------------------------------------------------
# 2. XML helpers (identical API to the two user-supplied templates)
# -----------------------------------------------------------------------------
esc <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;", s, fixed = TRUE)
  s <- gsub(">", "&gt;", s, fixed = TRUE)
  s <- gsub('"', "&quot;", s, fixed = TRUE)
  s
}
vtx <- function(id, value, style, x, y, w, h) {
  sprintf('        <mxCell id="%s" value="%s" style="%s" parent="1" vertex="1">\n          <mxGeometry x="%d" y="%d" width="%d" height="%d" as="geometry" />\n        </mxCell>\n',
    id, esc(value), style, x, y, w, h)
}
edg <- function(id, style, src, tgt, value = "") {
  val_attr <- if (nzchar(value)) sprintf(' value="%s"', esc(value)) else ""
  sprintf('        <mxCell id="%s"%s style="%s" parent="1" source="%s" target="%s" edge="1">\n          <mxGeometry relative="1" as="geometry" />\n        </mxCell>\n',
    id, val_attr, style, src, tgt)
}

# Reuse the template's style vocabulary (NPG-adjacent colours).
S <- list(
  start  = "rounded=1;whiteSpace=wrap;html=1;fillColor=#2C3E50;strokeColor=#2C3E50;fontColor=#FFFFFF;fontSize=13;",
  stage  = "rounded=1;whiteSpace=wrap;html=1;fillColor=#2980B9;strokeColor=#2471A3;fontColor=#FFFFFF;fontSize=11;",
  trig   = "rounded=1;whiteSpace=wrap;html=1;fillColor=#C0392B;strokeColor=#A93226;fontColor=#FFFFFF;fontSize=11;",
  exempt = "rounded=1;whiteSpace=wrap;html=1;fillColor=#27AE60;strokeColor=#1E8449;fontColor=#FFFFFF;fontSize=11;",
  grey   = "rounded=1;whiteSpace=wrap;html=1;fillColor=#7F8C8D;strokeColor=#616A6B;fontColor=#FFFFFF;fontSize=11;",
  never  = "rounded=1;whiteSpace=wrap;html=1;fillColor=#6C757D;strokeColor=#566573;fontColor=#FFFFFF;fontSize=11;",
  layer2 = "rounded=1;whiteSpace=wrap;html=1;fillColor=#8E44AD;strokeColor=#7D3C98;fontColor=#FFFFFF;fontSize=12;",
  censor = "rounded=1;whiteSpace=wrap;html=1;fillColor=#FFF3E0;strokeColor=#E67E22;strokeWidth=2;fontSize=10;",
  case3  = "rounded=1;whiteSpace=wrap;html=1;fillColor=#FCE4EC;strokeColor=#AD1457;strokeWidth=2;fontSize=10;",
  arrow  = "edgeStyle=orthogonalEdgeStyle;html=1;strokeColor=#444444;strokeWidth=2;",
  arrow2 = "edgeStyle=orthogonalEdgeStyle;html=1;strokeColor=#444444;strokeWidth=1.5;",
  feed   = "edgeStyle=orthogonalEdgeStyle;html=1;strokeColor=#9B59B6;strokeWidth=2;fontColor=#9B59B6;fontSize=10;fontStyle=1;labelBackgroundColor=#FAFBFC;",
  dash_o = "edgeStyle=orthogonalEdgeStyle;html=1;strokeColor=#E67E22;strokeWidth=1;dashed=1;",
  dash_p = "edgeStyle=orthogonalEdgeStyle;html=1;strokeColor=#AD1457;strokeWidth=1;dashed=1;",
  result = "rounded=1;whiteSpace=wrap;html=1;fillColor=#16A085;strokeColor=#138D75;fontColor=#FFFFFF;fontSize=12;",
  orange = "rounded=1;whiteSpace=wrap;html=1;fillColor=#E67E22;strokeColor=#CA6F1E;fontColor=#FFFFFF;fontSize=12;",
  pink   = "rounded=1;whiteSpace=wrap;html=1;fillColor=#AD1457;strokeColor=#880E4F;fontColor=#FFFFFF;fontSize=12;",
  text   = "text;strokeColor=none;fillColor=none;html=1;fontSize=14;align=center;whiteSpace=wrap;",
  note   = "text;strokeColor=none;fillColor=none;html=1;fontSize=10;align=left;fontColor=#888888;fontStyle=2;whiteSpace=wrap;",
  rule_box = "rounded=1;whiteSpace=wrap;html=1;fillColor=#8E44AD;strokeColor=#7D3C98;fontColor=#FFFFFF;fontSize=11;",
  comply_detail = "rounded=1;whiteSpace=wrap;html=1;fillColor=#E8F8F5;strokeColor=#16A085;strokeWidth=2;fontSize=10;fontColor=#333333;",
  ph_bad = "rounded=1;whiteSpace=wrap;html=1;fillColor=#FDECEA;strokeColor=#C0392B;strokeWidth=2;fontSize=10;fontColor=#922B21;",
  ph_ok  = "rounded=1;whiteSpace=wrap;html=1;fillColor=#EAF7EF;strokeColor=#16A085;strokeWidth=2;fontSize=10;fontColor=#117864;"
)

pct_n <- function(n, tot) sprintf("%.1f%%", n / max(tot, 1) * 100)

fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("&lt; 0.001")
  return(sprintf("= %.3f", p))
}

# -----------------------------------------------------------------------------
# 3. Per-cohort data aggregation
# -----------------------------------------------------------------------------
DYNAMIC_GRACE_DAYS <- 90

aggregate_cohort <- function(cohort_dir) {
  req_files <- c("R_main_results.csv", "R_Table3_RMST.csv", "R_trigger_table.csv",
                 "R_clone_dataset.csv", "R_cox_zph_test.csv",
                 "pathway_classification_table.csv")
  for (f in req_files) {
    p <- file.path(cohort_dir, f)
    if (!file.exists(p)) stop("Missing required file: ", p)
  }
  main  <- read.csv(file.path(cohort_dir, "R_main_results.csv"), stringsAsFactors = FALSE)
  rmst  <- read.csv(file.path(cohort_dir, "R_Table3_RMST.csv"), stringsAsFactors = FALSE)
  trig  <- read.csv(file.path(cohort_dir, "R_trigger_table.csv"), stringsAsFactors = FALSE)
  clone <- read.csv(file.path(cohort_dir, "R_clone_dataset.csv"), stringsAsFactors = FALSE)
  ph    <- read.csv(file.path(cohort_dir, "R_cox_zph_test.csv"), stringsAsFactors = FALSE)
  pway  <- read.csv(file.path(cohort_dir, "pathway_classification_table.csv"),
                    stringsAsFactors = FALSE)

  # ── Top-line numbers
  out <- list(
    N = as.integer(main$n_eligible[1]),
    n_dyn = as.integer(main$n_dyn[1]),
    n_early = as.integer(main$n_early[1]),
    events_dyn = as.integer(main$events_dyn[1]),
    events_early = as.integer(main$events_early[1]),
    censored_dyn = as.integer(main$censored_dyn[1]),
    censored_early = as.integer(main$censored_early[1]),
    hr = as.numeric(main$HR[1]),
    hr_lo = as.numeric(main$HR_lo[1]),
    hr_hi = as.numeric(main$HR_hi[1]),
    hr_p = as.numeric(main$HR_p[1]),
    evalue = as.numeric(main$E_value[1]),
    ess_dyn = as.numeric(main$ess_dyn[1]),
    ess_early = as.numeric(main$ess_early[1]),
    ph_p = as.numeric(main$ph_test_p[1]),
    ph_fail = !is.na(as.numeric(main$ph_test_p[1])) && as.numeric(main$ph_test_p[1]) < 0.05,
    rmst = rmst
  )

  # ── Layer / cycle buckets from trigger_stage
  trig <- trig %>%
    mutate(
      trigger_stage = as.character(trigger_stage),
      layer = case_when(
        str_detect(trigger_stage, "^layer1") ~ "L1",
        str_detect(trigger_stage, "^layer2") ~ "L2",
        TRUE ~ "Never"
      ),
      cycle = suppressWarnings(as.integer(str_replace(
        str_extract(trigger_stage, "cycle\\d+"),
        "cycle", ""))),
      l2_rule = ifelse(layer == "L2",
                       str_replace(trigger_stage, "^layer2_post_haic_", ""),
                       NA_character_)
    )

  # Layer 1 by cycle (collapse cycle 4+ into one bucket for the diagram)
  l1 <- trig %>% filter(layer == "L1")
  out$l1_total <- nrow(l1)
  out$l1_cycle3  <- sum(l1$cycle == 3, na.rm = TRUE)
  out$l1_cycle4  <- sum(l1$cycle == 4, na.rm = TRUE)
  out$l1_cycle5p <- sum(l1$cycle >= 5, na.rm = TRUE)

  # Layer 1 "dominant rule" breakdown — parse rule tokens after cycleN_
  rule_tokens <- l1 %>%
    mutate(rules = str_replace(trigger_stage, "^layer1_cycle\\d+_", "")) %>%
    pull(rules)
  token_split <- lapply(strsplit(rule_tokens, "\\+", perl = TRUE), unique)
  token_flat <- unlist(token_split)
  top_tokens <- sort(table(token_flat), decreasing = TRUE)
  out$l1_top_rules <- top_tokens

  # Layer 2
  l2 <- trig %>% filter(layer == "L2")
  out$l2_total <- nrow(l2)
  out$l2_afp  <- sum(l2$l2_rule == "afp", na.rm = TRUE)
  out$l2_pivka <- sum(l2$l2_rule == "pivka", na.rm = TRUE)
  out$l2_afp_pivka <- sum(l2$l2_rule == "afp+pivka", na.rm = TRUE)

  # Never
  out$never_total <- sum(trig$layer == "Never")

  # ── Join with pathway_classification_table to compute Case 1/2/3
  pway <- pway %>% transmute(
    patient_id,
    actual_day = as.integer(actual_day),
    pw_trigger_day = as.integer(trigger_day),
    pw_eff_triggered = as.integer(eff_triggered),
    arm_allegiance = as.character(arm_allegiance)
  )
  joined <- trig %>% left_join(pway, by = "patient_id") %>%
    mutate(
      actual_day = ifelse(is.na(actual_day), 9999L, actual_day),
      # Prefer trigger_day from trigger_table (authoritative) but fall back
      trig_day = ifelse(is.na(trigger_day), pw_trigger_day, trigger_day)
    )

  # For each patient compute their case class under dynamic arm
  joined <- joined %>% mutate(
    case_cls = case_when(
      eff_triggered == 1 & actual_day < trig_day                               ~ "Case3",
      eff_triggered == 1 & actual_day > (trig_day + DYNAMIC_GRACE_DAYS)        ~ "Case1",
      eff_triggered == 1 & actual_day >= trig_day & actual_day <= (trig_day + DYNAMIC_GRACE_DAYS) ~ "Comply",
      eff_triggered == 1 & actual_day == 9999                                  ~ "Comply_noimm",
      eff_triggered == 0 & actual_day < 9999                                   ~ "Case2",
      eff_triggered == 0 & actual_day == 9999                                  ~ "Never_truly",
      TRUE                                                                     ~ "Other"
    )
  )

  # Tally by layer × case
  cc <- joined %>%
    mutate(layer = ifelse(layer %in% c("L1", "L2"), layer, "Never")) %>%
    group_by(layer, case_cls) %>%
    summarise(n = n(), .groups = "drop")
  get_cc <- function(layer_val, case_val) {
    v <- cc %>% filter(layer == layer_val, case_cls == case_val) %>% pull(n)
    if (length(v) == 0) 0L else as.integer(sum(v))
  }
  out$l1_comply <- get_cc("L1", "Comply") + get_cc("L1", "Comply_noimm")
  out$l1_case1  <- get_cc("L1", "Case1")
  out$l1_case3  <- get_cc("L1", "Case3")
  out$l2_comply <- get_cc("L2", "Comply") + get_cc("L2", "Comply_noimm")
  out$l2_case1  <- get_cc("L2", "Case1")
  out$l2_case3  <- get_cc("L2", "Case3")
  out$never_case2 <- get_cc("Never", "Case2")
  out$never_truly <- get_cc("Never", "Never_truly")

  # Aggregate totals
  out$total_triggered <- out$l1_total + out$l2_total
  out$total_case1 <- out$l1_case1 + out$l2_case1
  out$total_case3 <- out$l1_case3 + out$l2_case3
  out$total_case2 <- out$never_case2
  out$total_cens_pre <- out$total_case1 + out$total_case2 + out$total_case3

  # Early-combo arm breakdown: from clone data
  ec_clone <- clone %>% filter(arm == "early_combo")
  out$ec_n <- nrow(ec_clone)
  out$ec_cens <- as.integer(sum(ec_clone$censored))
  out$ec_comply <- out$ec_n - out$ec_cens
  out$ec_events <- as.integer(sum(ec_clone$event))

  out
}

# -----------------------------------------------------------------------------
# 4. Build Page 1 (Dynamic/AoD arm) for one cohort
# -----------------------------------------------------------------------------
build_page1 <- function(d, meta, cohort_key) {
  pct <- function(n) pct_n(n, d$N)
  page <- '    <mxGraphModel dx="1400" dy="1200" grid="1" gridSize="10" guides="1" tooltips="1" connect="1" arrows="1" fold="1" page="1" pageScale="1" pageWidth="1800" pageHeight="1700" math="0" shadow="0">
      <root>
        <mxCell id="0" />
        <mxCell id="1" parent="0" />
'

  # ── Title
  page <- paste0(page, vtx("title",
    sprintf(paste0("<b>TTE Dynamic (Adaptive On-Demand) Arm \u2014 %s</b><br>",
                   "<b>Two-Layer Trigger Architecture + Case 1/2/3 Censoring</b><br>",
                   "Clone-Censor-Weight Framework | %s"),
            meta$title, cohort_key),
    S$text, 260, 10, 1200, 60))

  # ── Total cohort
  page <- paste0(page, vtx("total",
    sprintf("<b>Total Cohort (Cloned into Dynamic Arm)</b><br>N = %d", d$N),
    S$start, 620, 90, 460, 55))

  # ── Layer 1 header
  page <- paste0(page, vtx("l1_header",
    "<b>LAYER 1: Pre-HAIC Trigger Evaluation (cycles &ge; 3)</b>",
    "text;strokeColor=none;fillColor=none;html=1;fontSize=13;align=center;whiteSpace=wrap;fontStyle=1;fontColor=#2980B9;",
    380, 175, 500, 25))

  # ── Rule box (cohort-specific text)
  page <- paste0(page, vtx("rule_box", meta$rule_text,
    S$rule_box, 320, 205, 420, 80))
  page <- paste0(page, edg("e_total_rule", S$arrow, "total", "rule_box"))

  # ── L1 triggered + not triggered
  page <- paste0(page, vtx("l1_triggered",
    sprintf("<b>Layer 1 Triggered</b><br>n = %d (%s)<br><font style=\"font-size:9px;color:#FADBD8\">Cycle 3: %d | Cycle 4: %d | Cycle 5+: %d</font>",
      d$l1_total, pct(d$l1_total), d$l1_cycle3, d$l1_cycle4, d$l1_cycle5p),
    S$trig, 150, 325, 290, 65))

  page <- paste0(page, vtx("l1_not_trig",
    "<b>Layer 1 Not Triggered</b><br><font style=\"font-size:10px\">Rules never met OR insufficient cycle \u2265 3 data<br>\u2192 Forward to Layer 2</font>",
    S$exempt, 500, 325, 260, 65))
  page <- paste0(page, edg("e_rule_trig", S$arrow2, "rule_box", "l1_triggered"))
  page <- paste0(page, edg("e_rule_notrig", S$arrow2, "rule_box", "l1_not_trig"))

  # Top-3 rule tokens annotation
  tok <- d$l1_top_rules
  tok_keep <- head(tok, 6)
  tok_str <- paste(sprintf("%s: %d (%.0f%%)", names(tok_keep), as.integer(tok_keep),
                           as.integer(tok_keep) / max(d$l1_total, 1) * 100),
                   collapse = " | ")
  page <- paste0(page, vtx("l1_rules_note",
    sprintf("<font style=\"font-size:10px\"><b>Dominant L1 rule tokens (multi-rule patients counted in each token):</b><br>%s</font>",
            tok_str),
    S$note, 100, 400, 600, 45))

  # ── L1 outcomes (Comply / Case 1 / Case 3)
  page <- paste0(page, vtx("l1_comply_box",
    sprintf("<font color=\"#16A085\"><b>COMPLY</b></font><br>%s within [trigger, +%dd]<br>n = %d (%d%%)<br><font style=\"font-size:9px;color:#888\">Remain in Dynamic arm</font>",
      meta$early_abbr, DYNAMIC_GRACE_DAYS, d$l1_comply,
      round(d$l1_comply / max(d$l1_total, 1) * 100)),
    S$comply_detail, 20, 460, 180, 75))

  page <- paste0(page, vtx("l1_case1_box",
    sprintf("<b>CASE 1: Late</b><br>%s added &gt; trigger+%dd<br>n = %d (%d%%)<br><font style=\"font-size:9px;color:#888\">Censor at trigger+%dd</font>",
      meta$early_abbr, DYNAMIC_GRACE_DAYS, d$l1_case1,
      round(d$l1_case1 / max(d$l1_total, 1) * 100),
      DYNAMIC_GRACE_DAYS),
    S$censor, 210, 460, 180, 75))

  page <- paste0(page, vtx("l1_case3_box",
    sprintf("<b>CASE 3: Early</b><br>%s added BEFORE trigger<br>n = %d (%d%%)<br><font style=\"font-size:9px;color:#FFF\">Censor at %s start</font>",
      meta$early_abbr, d$l1_case3,
      round(d$l1_case3 / max(d$l1_total, 1) * 100),
      meta$early_abbr),
    S$pink, 400, 460, 180, 75))

  page <- paste0(page, edg("e_l1_comply", S$arrow2, "l1_triggered", "l1_comply_box"))
  page <- paste0(page, edg("e_l1_case1",  S$dash_o, "l1_triggered", "l1_case1_box"))
  page <- paste0(page, edg("e_l1_case3",  S$dash_p, "l1_triggered", "l1_case3_box"))

  # ── Layer 2
  page <- paste0(page, vtx("l2_header",
    "<b>LAYER 2: Post-HAIC AFP/PIVKA Continuous Monitoring (Fallback)</b>",
    "text;strokeColor=none;fillColor=none;html=1;fontSize=13;align=center;whiteSpace=wrap;fontStyle=1;fontColor=#8E44AD;",
    860, 175, 500, 25))

  page <- paste0(page, vtx("l2_box",
    "<b>Layer 2: Post-HAIC Monitoring</b><br><font style=\"font-size:9px\">Condition A: AFP &gt; 20 ng/mL<br>Condition B: AFP nadir &lt; 20 AND rise &gt; 1.3 ng/mL<br>Condition C: PIVKA abnormal</font>",
    S$layer2, 880, 205, 360, 80))
  page <- paste0(page, edg("e_l1nt_l2", S$feed, "l1_not_trig", "l2_box"))

  page <- paste0(page, vtx("l2_triggered",
    sprintf("<b>Layer 2 Triggered</b><br>n = %d (%s)<br><font style=\"font-size:9px;color:#FADBD8\">AFP only: %d | PIVKA only: %d | AFP+PIVKA: %d</font>",
      d$l2_total, pct(d$l2_total), d$l2_afp, d$l2_pivka, d$l2_afp_pivka),
    S$trig, 820, 325, 260, 65))

  page <- paste0(page, vtx("l2_never",
    sprintf("<b>Never Triggered</b><br>n = %d (%s)<br><font style=\"font-size:9px;color:#E8F8F5\">AFP/PIVKA normal or missing</font>",
      d$never_total, pct(d$never_total)),
    S$never, 1120, 325, 230, 65))

  page <- paste0(page, edg("e_l2_trig", S$arrow2, "l2_box", "l2_triggered"))
  page <- paste0(page, edg("e_l2_never", S$arrow2, "l2_box", "l2_never"))

  # ── L2 outcomes
  page <- paste0(page, vtx("l2_comply_box",
    sprintf("<font color=\"#16A085\"><b>COMPLY</b></font><br>%s within 90d<br>n = %d (%d%%)",
      meta$early_abbr, d$l2_comply,
      round(d$l2_comply / max(d$l2_total, 1) * 100)),
    S$comply_detail, 700, 460, 160, 70))

  page <- paste0(page, vtx("l2_case1_box",
    sprintf("<b>CASE 1: Late</b><br>&gt; trigger+90d<br>n = %d (%d%%)",
      d$l2_case1, round(d$l2_case1 / max(d$l2_total, 1) * 100)),
    S$censor, 870, 460, 140, 70))

  page <- paste0(page, vtx("l2_case3_box",
    sprintf("<b>CASE 3: Early</b><br>Added BEFORE trigger<br>n = %d (%d%%)",
      d$l2_case3, round(d$l2_case3 / max(d$l2_total, 1) * 100)),
    S$pink, 1020, 460, 160, 70))

  page <- paste0(page, edg("e_l2_comply", S$arrow2, "l2_triggered", "l2_comply_box"))
  page <- paste0(page, edg("e_l2_case1",  S$dash_o, "l2_triggered", "l2_case1_box"))
  page <- paste0(page, edg("e_l2_case3",  S$dash_p, "l2_triggered", "l2_case3_box"))

  # ── Never outcomes (Case 2 vs truly never)
  page <- paste0(page, vtx("never_case2",
    sprintf("<b>CASE 2: Got %s</b><br>Protocol deviation<br>n = %d (%d%%)<br><font style=\"font-size:9px;color:#888\">Censor at %s start</font>",
      meta$early_abbr, d$never_case2,
      round(d$never_case2 / max(d$never_total, 1) * 100),
      meta$early_abbr),
    S$censor, 1200, 460, 160, 75))

  page <- paste0(page, vtx("never_truly_box",
    sprintf("<font color=\"#16A085\"><b>COMPLY</b></font><br>Truly never got %s<br>n = %d (%d%%)<br><font style=\"font-size:9px;color:#888\">Remain uncensored</font>",
      meta$early_abbr, d$never_truly,
      round(d$never_truly / max(d$never_total, 1) * 100)),
    S$comply_detail, 1370, 460, 170, 75))

  page <- paste0(page, edg("e_nev_case2", S$dash_o, "l2_never", "never_case2"))
  page <- paste0(page, edg("e_nev_truly", S$arrow2, "l2_never", "never_truly_box"))

  # ── Summaries
  page <- paste0(page, vtx("sum_trig",
    sprintf("<b>Total Triggered: %d (%s)</b><br>Layer 1: %d | Layer 2 (fallback): %d",
      d$total_triggered, pct(d$total_triggered), d$l1_total, d$l2_total),
    S$trig, 100, 580, 460, 50))

  page <- paste0(page, vtx("sum_case1",
    sprintf("<b>Case 1 Total: %d</b><br>Late (&gt; trigger+90d)", d$total_case1),
    S$orange, 600, 580, 200, 50))

  page <- paste0(page, vtx("sum_case3",
    sprintf("<b>Case 3 Total: %d</b><br>Early (before trigger day)", d$total_case3),
    S$pink, 820, 580, 220, 50))

  page <- paste0(page, vtx("sum_case2",
    sprintf("<b>Case 2 Total: %d</b><br>Untriggered but got %s",
      d$total_case2, meta$early_abbr),
    S$orange, 1060, 580, 240, 50))

  page <- paste0(page, vtx("sum_cens_pre",
    sprintf("<b>Total Artificially Censored (pre-filter): %d</b><br>Case 1: %d | Case 2: %d | Case 3: %d   (each patient uniquely classified)",
      d$total_cens_pre, d$total_case1, d$total_case2, d$total_case3),
    S$orange, 300, 660, 800, 55))

  # ── Final result (cohort-specific PH banner colour)
  ph_banner <- if (d$ph_fail) {
    sprintf("\u26A0 PH assumption violated: P %s \u2192 RMST is primary estimand",
            fmt_p(d$ph_p))
  } else {
    sprintf("\u2713 PH assumption holds: P %s", fmt_p(d$ph_p))
  }
  ph_style <- if (d$ph_fail) S$ph_bad else S$ph_ok

  page <- paste0(page, vtx("ph_banner", ph_banner, ph_style, 300, 730, 800, 45))

  page <- paste0(page, vtx("sum_final",
    sprintf("<b>Dynamic Arm Final (after os_m &gt; 0 filter): n = %d</b><br>Events (deaths): %d | Retained artificially-censored rows: %d | ESS = %.0f<br><b>HR = %.2f (95%% CI %.2f\u2013%.2f), P %s | E-value = %.2f</b>",
      d$n_dyn, d$events_dyn, d$censored_dyn, d$ess_dyn,
      d$hr, d$hr_lo, d$hr_hi, fmt_p(d$hr_p), d$evalue),
    S$result, 300, 790, 800, 80))

  # ── Notes
  page <- paste0(page, vtx("note_layer2",
    "<font style=\"font-size:10px\"><i><b>Layer 2 rationale:</b> Many patients had &lt; 3 HAIC cycles, so Layer 1 cannot evaluate them. Layer 2 uses post-HAIC AFP/PIVKA monitoring as a fallback.</i></font>",
    S$note, 100, 890, 1100, 40))

  page <- paste0(page, vtx("note_params",
    sprintf("<font style=\"font-size:10px\"><i>Grace period: %d days from trigger day | Time zero: first HAIC date | IPCW truncation: 99th percentile | Bootstrap: 500 paired iterations</i></font>",
            DYNAMIC_GRACE_DAYS),
    S$note, 100, 935, 1100, 35))

  page <- paste0(page, vtx("cohort_note", meta$cohort_note,
    S$note, 100, 975, 1100, 35))

  page <- paste0(page, '      </root>\n    </mxGraphModel>')
  page
}

# -----------------------------------------------------------------------------
# 5. Build Page 2 (Early Combination arm + combined result card)
# -----------------------------------------------------------------------------
build_page2 <- function(d, meta, cohort_key) {
  r12 <- d$rmst %>% filter(tau == 12)
  r24 <- d$rmst %>% filter(tau == 24)
  r36 <- d$rmst %>% filter(tau == 36)
  rmst_str <- function(row) {
    if (nrow(row) == 0) return("NA")
    sprintf("%+.2f mo (%.2f, %.2f), P %s",
            row$delta[1], row$ci_lo[1], row$ci_hi[1], fmt_p(row$p[1]))
  }
  rmst12_str <- rmst_str(r12)
  rmst24_str <- rmst_str(r24)
  rmst36_str <- rmst_str(r36)

  cells <- ""
  cells <- paste0(cells, vtx("ec_title",
    sprintf("<b>TTE Early Combination Arm \u2014 %s</b><br><b>Patient Flow with Censoring</b><br>Clone-Censor-Weight Framework | %s",
            meta$title, cohort_key),
    S$text, 200, 20, 800, 55))

  cells <- paste0(cells, vtx("ec_total",
    sprintf("<b>Total Cohort (Cloned into Early Combo Arm)</b><br>N = %d", d$N),
    S$start, 320, 100, 360, 50))

  cells <- paste0(cells, vtx("ec_rule",
    sprintf("<b>Early Combination Rule</b><br><font style=\"font-size:11px\">Add %s within <b>14 days</b> of first HAIC</font>",
            meta$rule_ab),
    "rounded=1;whiteSpace=wrap;html=1;fillColor=#2980B9;strokeColor=#2471A3;fontColor=#FFFFFF;fontSize=12;",
    320, 210, 360, 55))
  cells <- paste0(cells, edg("ec_e1", S$arrow, "ec_total", "ec_rule"))

  cells <- paste0(cells, vtx("ec_comply",
    sprintf("<b>Compliant: Added %s \u2264 14 days</b><br>n = %d (%.1f%%)",
            meta$early_abbr, d$ec_comply, d$ec_comply / max(d$N, 1) * 100),
    "rounded=1;whiteSpace=wrap;html=1;fillColor=#16A085;strokeColor=#138D75;fontColor=#FFFFFF;fontSize=12;",
    160, 340, 280, 55))

  cells <- paste0(cells, vtx("ec_censored",
    sprintf("<b>Not Compliant: No %s \u2264 14 days</b><br>n = %d (%.1f%%)<br><font style=\"font-size:10px\">\u2192 Artificially Censored at Day 14</font>",
            meta$early_abbr, d$ec_cens, d$ec_cens / max(d$N, 1) * 100),
    "rounded=1;whiteSpace=wrap;html=1;fillColor=#E67E22;strokeColor=#CA6F1E;fontColor=#FFFFFF;fontSize=12;",
    540, 340, 300, 65))

  cells <- paste0(cells, edg("ec_e2", S$arrow, "ec_rule", "ec_comply"))
  cells <- paste0(cells, edg("ec_e3", S$arrow, "ec_rule", "ec_censored"))

  cells <- paste0(cells, vtx("ec_final",
    sprintf("<b>Early Combo Arm Final: n = %d</b><br>Events: %d | Censored: %d | ESS = %.0f",
            d$ec_n, d$events_early, d$ec_cens, d$ess_early),
    S$result, 260, 470, 420, 55))

  # Results card
  ph_str <- if (d$ph_fail) {
    sprintf("\u26A0 PH test P %s (violated)", fmt_p(d$ph_p))
  } else {
    sprintf("\u2713 PH test P %s (holds)", fmt_p(d$ph_p))
  }
  cells <- paste0(cells, vtx("results_box",
    sprintf(paste0("<b>Primary Analysis Results</b><br><br>",
                   "<b>HR (Dynamic vs Early Combo) = %.2f (95%% CI %.2f\u2013%.2f), P %s</b><br>",
                   "E-value = %.2f | %s<br><br>",
                   "Dynamic (AoD): n = %d, events = %d, ESS = %.0f<br>",
                   "Early combo: n = %d, events = %d, ESS = %.0f<br><br>",
                   "&Delta;RMST (12 mo): %s<br>",
                   "&Delta;RMST (24 mo): %s<br>",
                   "&Delta;RMST (36 mo): %s"),
            d$hr, d$hr_lo, d$hr_hi, fmt_p(d$hr_p), d$evalue, ph_str,
            d$n_dyn, d$events_dyn, d$ess_dyn,
            d$ec_n, d$events_early, d$ess_early,
            rmst12_str, rmst24_str, rmst36_str),
    "rounded=1;whiteSpace=wrap;html=1;fillColor=#F0F3F4;strokeColor=#2C3E50;strokeWidth=2;fontSize=11;fontColor=#2C3E50;",
    150, 560, 700, 220))

  cells <- paste0(cells, vtx("cohort_note", meta$cohort_note,
    S$note, 150, 800, 700, 40))

  sprintf('    <mxGraphModel dx="938" dy="900" grid="1" gridSize="10" guides="1" tooltips="1" connect="1" arrows="1" fold="1" page="1" pageScale="1" pageWidth="1200" pageHeight="900" math="0" shadow="0">\n      <root>\n        <mxCell id="0" />\n        <mxCell id="1" parent="0" />\n%s      </root>\n    </mxGraphModel>', cells)
}

# -----------------------------------------------------------------------------
# 6. Assemble .drawio file per cohort
# -----------------------------------------------------------------------------
generate_for_cohort <- function(cohort_key) {
  cohort_dir <- file.path(project_root, "output", "step3_tte",
                          "IT_RULES_R_two_cohorts", cohort_key)
  if (!dir.exists(cohort_dir)) stop("Cohort folder not found: ", cohort_dir)
  meta <- COHORT_META[[cohort_key]]
  if (is.null(meta)) stop("No metadata for cohort: ", cohort_key)

  cat(sprintf("\n=== Building drawio flow for %s ===\n", cohort_key))
  d <- aggregate_cohort(cohort_dir)

  p1 <- build_page1(d, meta, cohort_key)
  p2 <- build_page2(d, meta, cohort_key)
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S.000Z")
  drawio <- paste0(
    '<mxfile host="Electron" modified="', ts, '" version="26.0.0" pages="2">\n',
    sprintf('  <diagram id="tte-aod-flow-%s" name="Dynamic (AoD) Arm">\n', cohort_key),
    p1, '\n',
    '  </diagram>\n',
    sprintf('  <diagram id="tte-early-%s" name="Early Combo + Results">\n', cohort_key),
    p2, '\n',
    '  </diagram>\n',
    '</mxfile>\n')

  out_file <- file.path(cohort_dir, "TTE_AoD_Flow_Censoring.drawio")
  writeLines(drawio, out_file)

  # Console summary
  cat(sprintf("  N = %d | Triggered L1=%d L2=%d | Never=%d\n",
              d$N, d$l1_total, d$l2_total, d$never_total))
  cat(sprintf("  Case counts: Case1=%d, Case2=%d, Case3=%d (pre-filter total=%d)\n",
              d$total_case1, d$total_case2, d$total_case3, d$total_cens_pre))
  cat(sprintf("  Dynamic arm final: n=%d events=%d ESS=%.0f\n",
              d$n_dyn, d$events_dyn, d$ess_dyn))
  cat(sprintf("  Early-combo   final: n=%d events=%d ESS=%.0f\n",
              d$ec_n, d$events_early, d$ess_early))
  cat(sprintf("  HR=%.3f (%.2f-%.2f), P=%.4f, E-value=%.2f, PH P=%.4f %s\n",
              d$hr, d$hr_lo, d$hr_hi, d$hr_p, d$evalue, d$ph_p,
              ifelse(d$ph_fail, "[VIOLATED]", "[OK]")))
  cat(sprintf("  \u2192 Wrote %s\n", out_file))
}

# -----------------------------------------------------------------------------
# 7. Run for selected cohorts
# -----------------------------------------------------------------------------
for (ck in unique(cohorts)) generate_for_cohort(ck)

cat("\nDone. Open in VS Code (draw.io extension) or draw.io web/desktop app.\n")
