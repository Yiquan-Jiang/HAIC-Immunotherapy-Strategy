#!/usr/bin/env Rscript
# =============================================================================
# Generate TTE Patient Flow Draw.io Diagram — IT_RULES v1.2 cohort_3matched
# =============================================================================
# Reproduces patient flow for:
#   tte_IT_R_core_cohort_3matched.R (with Layer 2 fallback + Case 3 censoring)
#   Cohort: matched_06 (HAIC_alone vs HAIC_then_I+T) + all HAIC+I+T_concurrent
#   Strategy A: Dynamic (HAIC then I+T on-demand, IT_RULES_v1.2)
#     - Case 1: triggered, added I+T late (> grace) -> censor at trigger+grace
#     - Case 2: untriggered, but added I+T          -> censor at I+T start
#     - Case 3: triggered, added I+T BEFORE trigger -> censor at I+T start (new)
#   Strategy B: Early Combination (I+T within 14 days)
#
# Page 1: Dynamic Strategy Arm (Two-Layer Architecture, with Case 3)
# Page 2: Early Combination Arm + Results
#
# Usage:
#   Rscript generate_tte_flow_drawio_IT_cohort3matched.R [output_dir]
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
script_dir <- tryCatch(
  dirname(normalizePath(sys.frames()[[1]]$ofile)),
  error = function(e) getwd()
)
project_root <- normalizePath(file.path(script_dir, "../.."), winslash = "/", mustWork = FALSE)

if (length(args) >= 1 && nzchar(args[1])) {
  out_dir <- normalizePath(args[1], winslash = "/", mustWork = FALSE)
} else {
  out_dir <- file.path(project_root, "output", "step3_tte", "IT_RULES_R", "cohort_3matched")
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Flow data — from tte_IT_R_core_cohort_3matched.R output (v1.2)
# =============================================================================

# Total cohort
N <- 1938  # matched_06 (219+219) + all HAIC+I+T_concurrent (1500)

# Cohort composition
n_haic_alone       <- 219
n_haic_then_it     <- 219
n_haic_conc        <- 1500

# ── Dynamic arm: Layer 1 Trigger Evaluation ──
l1_total_trig <- 823      # triggered in Layer 1
l1_cycle3     <- 712
l1_cycle4     <- 80
l1_cycle5p    <- 31

l1_afp_only   <- 131
l1_plr_only   <- 514
l1_afp_plr    <- 178

# Layer 1 outcomes
l1_comply     <- 80       # comply: added I+T within [trigger, trigger+90d]
l1_case1      <- 143      # Case 1: added too LATE (> trigger+90d)
l1_case3      <- 600      # Case 3: added EARLIER than trigger

# ── Dynamic arm: Layer 2 (post-HAIC AFP monitoring fallback) ──
l2_trig       <- 269
l2_comply     <- 9
l2_case1      <- 59
l2_case3      <- 201

# ── Never triggered ──
total_never    <- 846
never_case2    <- 665    # Case 2: got I+T despite never triggering
never_truly    <- 181    # truly never got I+T

total_triggered <- l1_total_trig + l2_trig  # 1092

# ── Censoring totals (pre-filter) ──
total_case1 <- l1_case1 + l2_case1             # 202
total_case2 <- never_case2                      # 665
total_case3 <- l1_case3 + l2_case3             # 801
total_cens_pre <- total_case1 + total_case2 + total_case3  # 1668

# ── Dynamic arm final (post filter os_m>0, from R_main_results.csv) ──
dyn_arm_n      <- 623
dyn_events     <- 129
dyn_censored   <- 353     # after filter (many censored at day 0 dropped)

# ── Early combo arm ──
ec_arm_n      <- 1938
ec_events     <- 541
ec_censored   <- 581

# ── Results ──
hr_val   <- "0.67"
hr_ci    <- "0.55-0.82"
hr_p     <- "<0.001"
e_val    <- "2.33"
ph_p     <- "0.229"
dyn_ess  <- 534
ec_ess   <- 1837
rmst24   <- "+1.97 mo (0.81, 2.53), P <0.001"
rmst36   <- "+3.37 mo (1.24, 4.41), P = 0.004"

# =============================================================================
# Helpers
# =============================================================================
esc <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;", s, fixed = TRUE)
  s <- gsub(">", "&gt;", s, fixed = TRUE)
  s <- gsub('"', "&quot;", s, fixed = TRUE)
  s
}

vtx <- function(id, value, style, x, y, w, h) {
  sprintf('        <mxCell id="%s" value="%s" style="%s" parent="1" vertex="1">
          <mxGeometry x="%d" y="%d" width="%d" height="%d" as="geometry" />
        </mxCell>\n', id, esc(value), style, x, y, w, h)
}

edg <- function(id, style, src, tgt, value = "") {
  val_attr <- if (nzchar(value)) sprintf(' value="%s"', esc(value)) else ""
  pts_xml <- '\n          <mxGeometry relative="1" as="geometry" />'
  sprintf('        <mxCell id="%s"%s style="%s" parent="1" source="%s" target="%s" edge="1">%s\n        </mxCell>\n',
    id, val_attr, style, src, tgt, pts_xml)
}

# Style constants
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
  comply_detail = "rounded=1;whiteSpace=wrap;html=1;fillColor=#E8F8F5;strokeColor=#16A085;strokeWidth=2;fontSize=10;fontColor=#333333;"
)

pct_s <- function(n) sprintf("%.1f%%", n / N * 100)

# =============================================================================
# Page 1: Dynamic Strategy Arm (Two-Layer + Case 3)
# =============================================================================
p1 <- '    <mxGraphModel dx="1400" dy="1200" grid="1" gridSize="10" guides="1" tooltips="1" connect="1" arrows="1" fold="1" page="1" pageScale="1" pageWidth="1800" pageHeight="1700" math="0" shadow="0">
      <root>
        <mxCell id="0" />
        <mxCell id="1" parent="0" />
'

# ── Title ──
p1 <- paste0(p1, vtx("title",
  "<b>TTE IT_RULES v1.2: Dynamic Strategy Arm (HAIC then I+T On-Demand)</b><br><b>Two-Layer Trigger Architecture + Case 3 Censoring | Patient Flow</b><br>Clone-Censor-Weight Framework | cohort_3matched (matched_06 + all HAIC+I+T_conc)",
  S$text, 300, 10, 1100, 55))

# ── Total Cohort ──
p1 <- paste0(p1, vtx("total",
  sprintf("<b>Total Cohort (Cloned into Dynamic Arm)</b><br>N = %d<br><font style=\"font-size:9px;color:#D5D8DC\">HAIC_alone: %d | HAIC+I+T_concurrent: %d | HAIC_then_I+T: %d</font>",
    N, n_haic_alone, n_haic_conc, n_haic_then_it),
  S$start, 620, 90, 460, 65))

# ── Layer 1 Header & Rule Box ──
p1 <- paste0(p1, vtx("l1_header",
  "<b>LAYER 1: Pre-HAIC Trigger Evaluation (from Cycle 3 onwards)</b>",
  "text;strokeColor=none;fillColor=none;html=1;fontSize=13;align=center;whiteSpace=wrap;fontStyle=1;fontColor=#2980B9;",
  380, 185, 500, 25))

p1 <- paste0(p1, vtx("rule_box",
  "<b>Trigger Rules (evaluated at each pre-HAIC visit, cycle &ge; 3)</b><br><font style=\"font-size:10px\">Rule 1: AFP change from baseline &gt; -32.5% (AFP drop &lt; 32.5%)<br>Rule 2: PLR &gt; 102.4<br>First visit meeting Rule 1 OR Rule 2 &rarr; <b>Trigger</b></font>",
  S$rule_box, 320, 215, 420, 75))

p1 <- paste0(p1, edg("e_total_rule", S$arrow, "total", "rule_box"))

# ── Layer 1 Outcomes ──
p1 <- paste0(p1, vtx("l1_triggered",
  sprintf("<b>Layer 1 Triggered</b><br>n = %d (%s)<br><font style=\"font-size:9px;color:#FADBD8\">Cycle 3: %d | Cycle 4: %d | Cycle 5+: %d</font>",
    l1_total_trig, pct_s(l1_total_trig), l1_cycle3, l1_cycle4, l1_cycle5p),
  S$trig, 150, 335, 290, 65))

p1 <- paste0(p1, vtx("l1_not_trig",
  sprintf("<b>Layer 1 Not Triggered</b><br><font style=\"font-size:10px\">Rules never met OR insufficient cycle &ge; 3 data<br>&rarr; go to Layer 2</font>"),
  S$exempt, 500, 335, 250, 65))

p1 <- paste0(p1, edg("e_rule_trig", S$arrow2, "rule_box", "l1_triggered"))
p1 <- paste0(p1, edg("e_rule_notrig", S$arrow2, "rule_box", "l1_not_trig"))

# ── Layer 1 Rule Breakdown ──
p1 <- paste0(p1, vtx("l1_rules_note",
  sprintf("<font style=\"font-size:10px\"><b>L1 Rule Breakdown (n=%d):</b><br>AFP only: %d (%d%%) | PLR only: %d (%d%%) | AFP+PLR: %d (%d%%)</font>",
    l1_total_trig, l1_afp_only, round(l1_afp_only/l1_total_trig*100),
    l1_plr_only, round(l1_plr_only/l1_total_trig*100),
    l1_afp_plr, round(l1_afp_plr/l1_total_trig*100)),
  S$note, 100, 410, 380, 35))

# ── Layer 1 Three Outcomes (Comply / Case 1 / Case 3) ──
p1 <- paste0(p1, vtx("l1_comply_box",
  sprintf("<font color=\"#16A085\"><b>COMPLY</b></font><br>I+T within [trigger, trigger+90d]<br>n = %d (%d%%)<br><font style=\"font-size:9px;color:#888\">Remain in Dynamic arm</font>",
    l1_comply, round(l1_comply/l1_total_trig*100)),
  S$comply_detail, 20, 465, 180, 75))

p1 <- paste0(p1, vtx("l1_case1_box",
  sprintf("<b>CASE 1: Late</b><br>I+T added &gt; trigger+90d<br>n = %d (%d%%)<br><font style=\"font-size:9px;color:#888\">Censor at trigger+90d</font>",
    l1_case1, round(l1_case1/l1_total_trig*100)),
  S$censor, 210, 465, 180, 75))

p1 <- paste0(p1, vtx("l1_case3_box",
  sprintf("<b>CASE 3: Early</b><br>I+T added BEFORE trigger<br>n = %d (%d%%)<br><font style=\"font-size:9px;color:#FFF\">Censor at I+T start (NEW)</font>",
    l1_case3, round(l1_case3/l1_total_trig*100)),
  S$pink, 400, 465, 180, 75))

p1 <- paste0(p1, edg("e_l1_comply", S$arrow2, "l1_triggered", "l1_comply_box"))
p1 <- paste0(p1, edg("e_l1_case1",  S$dash_o, "l1_triggered", "l1_case1_box"))
p1 <- paste0(p1, edg("e_l1_case3",  S$dash_p, "l1_triggered", "l1_case3_box"))

# ── Layer 2 Header & Box ──
p1 <- paste0(p1, vtx("l2_header",
  "<b>LAYER 2: Post-HAIC AFP Continuous Monitoring (Fallback)</b>",
  "text;strokeColor=none;fillColor=none;html=1;fontSize=13;align=center;whiteSpace=wrap;fontStyle=1;fontColor=#8E44AD;",
  860, 185, 500, 25))

p1 <- paste0(p1, vtx("l2_box",
  "<b>Layer 2: Post-HAIC AFP Monitoring</b><br><font style=\"font-size:9px\">Condition A: AFP &gt; 20 ng/mL<br>Condition B: AFP nadir &lt; 20 AND rise from nadir &gt; 1.3 ng/mL<br>First post_haic visit meeting A or B &rarr; Trigger</font>",
  S$layer2, 880, 215, 360, 75))

p1 <- paste0(p1, edg("e_l1nt_l2", S$feed, "l1_not_trig", "l2_box"))

# ── Layer 2 Outcomes ──
p1 <- paste0(p1, vtx("l2_triggered",
  sprintf("<b>Layer 2 Triggered</b><br>(AFP abnormal)<br>n = %d (%s)",
    l2_trig, pct_s(l2_trig)),
  S$trig, 820, 335, 220, 60))

p1 <- paste0(p1, vtx("l2_never",
  sprintf("<b>Never Triggered</b><br>n = %d (%s)<br><font style=\"font-size:9px;color:#E8F8F5\">AFP normal / no post_haic data</font>",
    total_never, pct_s(total_never)),
  S$never, 1100, 335, 230, 60))

p1 <- paste0(p1, edg("e_l2_trig", S$arrow2, "l2_box", "l2_triggered"))
p1 <- paste0(p1, edg("e_l2_never", S$arrow2, "l2_box", "l2_never"))

# ── Layer 2 Three Outcomes (Comply / Case 1 / Case 3) ──
p1 <- paste0(p1, vtx("l2_comply_box",
  sprintf("<font color=\"#16A085\"><b>COMPLY</b></font><br>I+T within 90d of L2 trigger<br>n = %d (%d%%)",
    l2_comply, round(l2_comply/l2_trig*100)),
  S$comply_detail, 710, 465, 170, 70))

p1 <- paste0(p1, vtx("l2_case1_box",
  sprintf("<b>CASE 1: Late</b><br>&gt; trigger+90d<br>n = %d (%d%%)",
    l2_case1, round(l2_case1/l2_trig*100)),
  S$censor, 890, 465, 140, 70))

p1 <- paste0(p1, vtx("l2_case3_box",
  sprintf("<b>CASE 3: Early</b><br>Added BEFORE trigger<br>n = %d (%d%%)",
    l2_case3, round(l2_case3/l2_trig*100)),
  S$pink, 1040, 465, 150, 70))

p1 <- paste0(p1, edg("e_l2_comply", S$arrow2, "l2_triggered", "l2_comply_box"))
p1 <- paste0(p1, edg("e_l2_case1",  S$dash_o, "l2_triggered", "l2_case1_box"))
p1 <- paste0(p1, edg("e_l2_case3",  S$dash_p, "l2_triggered", "l2_case3_box"))

# ── Never Triggered Outcomes (Case 2) ──
p1 <- paste0(p1, vtx("never_case2",
  sprintf("<b>CASE 2: Got I+T</b><br>Protocol deviation<br>n = %d (%d%%)<br><font style=\"font-size:9px;color:#888\">Censor at I+T start</font>",
    never_case2, round(never_case2/total_never*100)),
  S$censor, 1200, 465, 150, 75))

p1 <- paste0(p1, vtx("never_truly_box",
  sprintf("<font color=\"#16A085\"><b>COMPLY</b></font><br>Truly never got I+T<br>n = %d (%d%%)<br><font style=\"font-size:9px;color:#888\">Remain uncensored</font>",
    never_truly, round(never_truly/total_never*100)),
  S$comply_detail, 1360, 465, 170, 75))

p1 <- paste0(p1, edg("e_nev_case2", S$dash_o, "l2_never", "never_case2"))
p1 <- paste0(p1, edg("e_nev_truly", S$arrow2, "l2_never", "never_truly_box"))

# ── Summary Boxes (Aggregate censoring) ──
p1 <- paste0(p1, vtx("sum_trig",
  sprintf("<b>Total Triggered: %d (%s)</b><br>Layer 1: %d (AFP/PLR rules) | Layer 2: %d (AFP monitoring fallback)",
    total_triggered, pct_s(total_triggered), l1_total_trig, l2_trig),
  S$trig, 80, 580, 480, 50))

p1 <- paste0(p1, vtx("sum_case1",
  sprintf("<b>Case 1 Total: %d</b><br>Late (added after trigger+90d)",
    total_case1),
  S$orange, 600, 580, 220, 50))

p1 <- paste0(p1, vtx("sum_case3",
  sprintf("<b>Case 3 Total: %d (NEW)</b><br>Early (added BEFORE trigger day)",
    total_case3),
  S$pink, 840, 580, 240, 50))

p1 <- paste0(p1, vtx("sum_case2",
  sprintf("<b>Case 2 Total: %d</b><br>Untriggered but received I+T",
    total_case2),
  S$orange, 1100, 580, 230, 50))

# ── Aggregate & Final ──
p1 <- paste0(p1, vtx("sum_censor",
  sprintf("<b>Total Artificially Censored (pre-filter): %d</b><br>Case 1: %d | Case 2: %d | Case 3: %d  (each patient uniquely classified)",
    total_cens_pre, total_case1, total_case2, total_case3),
  S$orange, 330, 660, 740, 55))

p1 <- paste0(p1, vtx("sum_final",
  sprintf("<b>Dynamic Arm Final (after os_m &gt; 0 filter): n = %d</b><br>Events (deaths): %d | Artificially censored (retained): %d | ESS = %d<br><b>HR = %s (95%% CI %s), P %s | E-value = %s | PH test P = %s</b>",
    dyn_arm_n, dyn_events, dyn_censored, dyn_ess, hr_val, hr_ci, hr_p, e_val, ph_p),
  S$result, 280, 740, 840, 75))

# ── Notes ──
p1 <- paste0(p1, vtx("note_case3",
  "<font style=\"font-size:11px;color:#AD1457\"><b>NEW — Case 3 censoring (v1.2):</b> Patients who added I+T BEFORE their trigger day<br>are now censored at their I+T start day in the Dynamic arm. This removes<br>801 early-combination adherents that previously polluted the Dynamic arm<br>under v1.1, and is the primary driver of the improved HR.</font>",
  "text;strokeColor=none;fillColor=none;html=1;fontSize=11;align=left;fontColor=#AD1457;whiteSpace=wrap;",
  100, 835, 800, 85))

p1 <- paste0(p1, vtx("note_layer2",
  "<font style=\"font-size:10px\"><i><b>Layer 2 rationale:</b> Many patients have &lt; 3 HAIC cycles, so Layer 1 (pre-HAIC cycle &ge; 3) cannot evaluate them.<br>Layer 2 uses post-HAIC AFP monitoring as a fallback. Layer 2 rescued 269 additional trigger events.</i></font>",
  S$note, 100, 925, 800, 45))

p1 <- paste0(p1, vtx("note_params",
  "<font style=\"font-size:10px\"><i>Grace period: 90 days from trigger day | Time zero: first HAIC date<br>Layer 1: AFP change &gt; -32.5%% OR PLR &gt; 102.4 | Layer 2: AFP &gt; 20 ng/mL</i></font>",
  S$note, 100, 975, 800, 40))

# Close page 1
p1 <- paste0(p1, '      </root>
    </mxGraphModel>')

# =============================================================================
# Page 2: Early Combination Arm + Combined Results
# =============================================================================
p2_cells <- ""

p2_cells <- paste0(p2_cells, vtx("ec_title",
  "<b>TTE IT_RULES v1.2: Early Combination Strategy Arm</b><br><b>Patient Flow with Censoring</b><br>Clone-Censor-Weight Framework | cohort_3matched",
  S$text, 200, 20, 700, 50))

p2_cells <- paste0(p2_cells, vtx("ec_total",
  sprintf("<b>Total Cohort (Cloned into Early Combo Arm)</b><br>N = %d", N),
  S$start, 320, 100, 360, 50))

p2_cells <- paste0(p2_cells, vtx("ec_rule",
  "<b>Early Combination Rule</b><br><font style=\"font-size:11px\">Add I+T (immunotherapy + targeted therapy) within <b>14 days</b> of first HAIC</font>",
  "rounded=1;whiteSpace=wrap;html=1;fillColor=#2980B9;strokeColor=#2471A3;fontColor=#FFFFFF;fontSize=12;",
  320, 210, 360, 55))
p2_cells <- paste0(p2_cells, edg("ec_e1", S$arrow, "ec_total", "ec_rule"))

ec_comply <- ec_arm_n - ec_censored

p2_cells <- paste0(p2_cells, vtx("ec_comply",
  sprintf("<b>Compliant: Added I+T &le; 14 days</b><br>n = %d (%.1f%%)",
    ec_comply, ec_comply/N*100),
  "rounded=1;whiteSpace=wrap;html=1;fillColor=#16A085;strokeColor=#138D75;fontColor=#FFFFFF;fontSize=12;",
  160, 340, 280, 55))

p2_cells <- paste0(p2_cells, vtx("ec_censored",
  sprintf("<b>Not Compliant: No I+T &le; 14 days</b><br>n = %d (%.1f%%)<br><font style=\"font-size:10px\">&rarr; Artificially Censored at Day 14</font>",
    ec_censored, ec_censored/N*100),
  "rounded=1;whiteSpace=wrap;html=1;fillColor=#E67E22;strokeColor=#CA6F1E;fontColor=#FFFFFF;fontSize=12;",
  530, 340, 300, 65))

p2_cells <- paste0(p2_cells, edg("ec_e2", S$arrow, "ec_rule", "ec_comply"))
p2_cells <- paste0(p2_cells, edg("ec_e3", S$arrow, "ec_rule", "ec_censored"))

p2_cells <- paste0(p2_cells, vtx("ec_final",
  sprintf("<b>Early Combo Arm Final: n = %d</b><br>Events: %d | Censored: %d | ESS = %d",
    ec_arm_n, ec_events, ec_censored, ec_ess),
  S$result, 250, 470, 400, 55))

# Combined results box
p2_cells <- paste0(p2_cells, vtx("results_box",
  sprintf("<b>Primary Analysis Results (v1.2)</b><br><br><b>HR (Dynamic vs Early Combo) = %s (95%% CI %s), P %s</b><br>E-value = %s | PH test P = %s<br><br>Dynamic: n=%d, events=%d, ESS=%d<br>Early combo: n=%d, events=%d, ESS=%d<br><br>&Delta;RMST (24 mo): %s<br>&Delta;RMST (36 mo): %s",
    hr_val, hr_ci, hr_p, e_val, ph_p,
    dyn_arm_n, dyn_events, dyn_ess,
    ec_arm_n, ec_events, ec_ess,
    rmst24, rmst36),
  "rounded=1;whiteSpace=wrap;html=1;fillColor=#F0F3F4;strokeColor=#2C3E50;strokeWidth=2;fontSize=11;fontColor=#2C3E50;",
  150, 560, 600, 200))

# Cohort note
p2_cells <- paste0(p2_cells, vtx("cohort_note",
  sprintf("<font style=\"font-size:10px\"><i><b>Cohort composition (N = %d):</b><br>&bull; matched_06 (HAIC_alone vs HAIC_then_I+T, 1:1 PSM): %d<br>&bull; all HAIC+I+T_concurrent (unmatched): %d<br><br><b>Rationale:</b> delayed-vs-alone PSM backbone + full concurrent pool as counterfactual early-combo.</i></font>",
    N, n_haic_alone + n_haic_then_it, n_haic_conc),
  S$note, 150, 780, 600, 100))

p2 <- sprintf('    <mxGraphModel dx="938" dy="900" grid="1" gridSize="10" guides="1" tooltips="1" connect="1" arrows="1" fold="1" page="1" pageScale="1" pageWidth="1200" pageHeight="1000" math="0" shadow="0">
      <root>
        <mxCell id="0" />
        <mxCell id="1" parent="0" />
%s      </root>
    </mxGraphModel>', p2_cells)

# =============================================================================
# Assemble and write
# =============================================================================
ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S.000Z")
drawio <- paste0('<mxfile host="Electron" modified="', ts, '" version="26.0.0" pages="2">\n',
  '  <diagram id="tte-it-flow-v1.2" name="Dynamic Strategy (Two-Layer + Case 3)">\n',
  p1, '\n',
  '  </diagram>\n',
  '  <diagram id="early-combo-results" name="Early Combo + Results">\n',
  p2, '\n',
  '  </diagram>\n',
  '</mxfile>\n')

out_file <- file.path(out_dir, "TTE_IT_RULES_Flow_Censoring_cohort3matched.drawio")
writeLines(drawio, out_file)
cat(sprintf("Draw.io file written: %s\n", out_file))
cat("Open in VS Code (draw.io extension) or draw.io web/desktop app.\n")
cat(sprintf("\nKey flow numbers (v1.2):\n"))
cat(sprintf("  Total cohort: N = %d (HAIC_alone=%d + HAIC_then_I+T=%d + HAIC+I+T_conc=%d)\n",
    N, n_haic_alone, n_haic_then_it, n_haic_conc))
cat(sprintf("  Layer 1 triggered: %d | Layer 2 triggered: %d | Never: %d\n",
    l1_total_trig, l2_trig, total_never))
cat(sprintf("  Censoring: Case 1 (late)=%d, Case 2 (untrig got I+T)=%d, Case 3 (early, NEW)=%d\n",
    total_case1, total_case2, total_case3))
cat(sprintf("  Dynamic arm final: n=%d, events=%d, censored=%d, ESS=%d\n",
    dyn_arm_n, dyn_events, dyn_censored, dyn_ess))
cat(sprintf("  Early combo:       n=%d, events=%d, censored=%d, ESS=%d\n",
    ec_arm_n, ec_events, ec_censored, ec_ess))
cat(sprintf("  Primary: HR = %s (%s), P %s | E-value = %s\n", hr_val, hr_ci, hr_p, e_val))
