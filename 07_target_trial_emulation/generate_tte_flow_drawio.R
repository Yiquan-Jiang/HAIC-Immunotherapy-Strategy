#!/usr/bin/env Rscript
# =============================================================================
# Generate TTE Patient Flow Draw.io Diagram
# =============================================================================
# Reproduces: TTE_Dynamic_Flow_Censoring.drawio
#   Page 1: Dynamic Strategy Arm + Early Combo (side by side on same page)
#   Page 2: Early Combination Arm (standalone view)
#
# Usage:
#   Rscript generate_tte_flow_drawio.R [output_dir]
#   Default output: ../output/step3_tte/PIV_BASED_RULES_R/cohort_7group_psm02/
#
# Data source: R_clone_dataset.csv + R_trigger_table.csv from TTE analysis
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
  out_dir <- file.path(project_root, "output", "step3_tte", "PIV_BASED_RULES_R", "cohort_7group_psm02")
}
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Flow data — update these if TTE rules or cohort change
# =============================================================================
# Dynamic arm flow
N            <- 574
no_cycle3    <- 250    # no cycle3 data → Layer 2
enter_s1     <- 324    # entered Stage 1
s1_trig      <- 136    # S1 triggered + not exempt → add immune
s1_to_s2     <- 188    # S1 exempt (100) + no trigger (88) → Stage 2
s1_exempt    <- 100
s1_no_trig   <- 88
s2_trig      <- 19     # S2 triggered → add immune
s2_pass      <- 169    # S2 pass → Stage 3
s3_trig      <- 111    # S3 triggered → add immune
s3_exempt    <- 34     # S3 exempt → Layer 2
s3_nodata    <- 24     # S3 no data → Layer 2
l2_total     <- 308    # Layer 2 total
l2_trig      <- 65     # L2 AFP triggered
l2_never     <- 243    # L2 never triggered
total_trig   <- s1_trig + s2_trig + s3_trig + l2_trig  # 331

# Dynamic arm censoring (from clone dataset)
s1_comply <- 104; s1_censor <- 69
s2_comply <- 12;  s2_censor <- 6
s3_comply <- 23;  s3_censor <- 13
l2_comply <- 47;  l2_censor <- 37
never_excluded     <- 132  # censored at immune start (removed from dynamic arm)
never_in_dyn       <- 131  # entered dynamic arm
never_censor_later <- 33   # later got immune → censored
never_remain       <- 98   # truly never got immune
total_art_censor   <- s1_censor + s2_censor + s3_censor + l2_censor + never_censor_later  # 158
dyn_arm_n          <- 442
dyn_events         <- 145

# Early combo arm
ec_comply    <- 234    # added immune within 14 days
ec_censored  <- 340    # censored at day 14
ec_events    <- 89     # deaths in comply group
ec_alive     <- 145    # alive in comply group
ec_no_immune <- 200    # never received immune
ec_after_14d <- 140    # immune after day 14
ec_ess       <- 548

# HR results
hr_val   <- "0.770"
hr_ci    <- "0.582-1.019"
hr_p     <- "0.068"
e_val    <- "1.92"
ph_p     <- "0.004"
dyn_ess  <- 350

# =============================================================================
# Helper: XML escape
# =============================================================================
esc <- function(s) {
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;", s, fixed = TRUE)
  s <- gsub(">", "&gt;", s, fixed = TRUE)
  s <- gsub('"', "&quot;", s, fixed = TRUE)
  s
}

# =============================================================================
# Build Page 1: Dynamic Strategy + Early Combo side by side
# =============================================================================
p1 <- '    <mxGraphModel dx="1340" dy="1004" grid="1" gridSize="10" guides="1" tooltips="1" connect="1" arrows="1" fold="1" page="1" pageScale="1" pageWidth="2600" pageHeight="1400" math="0" shadow="0">
      <root>
        <mxCell id="0" />
        <mxCell id="1" parent="0" />
'

# --- Vertex helper ---
vtx <- function(id, value, style, x, y, w, h) {
  sprintf('        <mxCell id="%s" value="%s" style="%s" parent="1" vertex="1">
          <mxGeometry x="%d" y="%d" width="%d" height="%d" as="geometry" />
        </mxCell>\n', id, value, style, x, y, w, h)
}

# --- Edge helper ---
edg <- function(id, style, src, tgt, value = "", pts = NULL) {
  val_attr <- if (nzchar(value)) sprintf(' value="%s"', value) else ""
  pts_xml <- ""
  if (!is.null(pts)) {
    pts_xml <- sprintf('\n          <mxGeometry relative="1" as="geometry">\n            <Array as="points">\n%s\n            </Array>\n          </mxGeometry>',
      paste(sprintf('              <mxPoint x="%d" y="%d" />', pts[,1], pts[,2]), collapse = "\n"))
  } else {
    pts_xml <- '\n          <mxGeometry relative="1" as="geometry" />'
  }
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
  layer2 = "rounded=1;whiteSpace=wrap;html=1;fillColor=#8E44AD;strokeColor=#7D3C98;fontColor=#FFFFFF;fontSize=12;",
  censor = "rounded=1;whiteSpace=wrap;html=1;fillColor=#FFF3E0;strokeColor=#E67E22;strokeWidth=2;fontSize=10;",
  never  = "rounded=1;whiteSpace=wrap;html=1;fillColor=#6C757D;strokeColor=#566573;fontColor=#FFFFFF;fontSize=11;",
  arrow  = "edgeStyle=orthogonalEdgeStyle;html=1;strokeColor=#444444;strokeWidth=2;",
  arrow2 = "edgeStyle=orthogonalEdgeStyle;html=1;strokeColor=#444444;strokeWidth=1.5;",
  feed   = "edgeStyle=orthogonalEdgeStyle;html=1;strokeColor=#9B59B6;strokeWidth=2;fontColor=#9B59B6;fontSize=10;fontStyle=1;labelBackgroundColor=#FAFBFC;",
  dash   = "edgeStyle=orthogonalEdgeStyle;html=1;strokeColor=#E67E22;strokeWidth=1;dashed=1;",
  result = "rounded=1;whiteSpace=wrap;html=1;fillColor=#16A085;strokeColor=#138D75;fontColor=#FFFFFF;fontSize=12;",
  orange = "rounded=1;whiteSpace=wrap;html=1;fillColor=#E67E22;strokeColor=#CA6F1E;fontColor=#FFFFFF;fontSize=12;",
  text   = "text;strokeColor=none;fillColor=none;html=1;fontSize=14;align=center;whiteSpace=wrap;",
  note   = "text;strokeColor=none;fillColor=none;html=1;fontSize=10;align=left;fontColor=#888888;fontStyle=2;whiteSpace=wrap;",
  l2src  = "rounded=1;whiteSpace=wrap;html=1;fillColor=#F3E5F5;strokeColor=#9B59B6;strokeWidth=1;dashed=1;fontSize=9;fontColor=#555555;",
  cendet = "rounded=1;whiteSpace=wrap;html=1;fillColor=#FFF3E0;strokeColor=#E67E22;strokeWidth=2;fontSize=11;fontColor=#333333;align=left;spacingLeft=15;",
  comply_detail = "rounded=1;whiteSpace=wrap;html=1;fillColor=#FFF3E0;strokeColor=#16A085;strokeWidth=2;fontSize=11;fontColor=#333333;",
  line   = "shape=line;strokeColor=#CCCCCC;strokeWidth=1;dashed=1;html=1;"
)

pct_s <- function(n) sprintf("%.1f%%", n / N * 100)

# ── Dynamic arm (left half, x offset = 0) ──

# Title
p1 <- paste0(p1, vtx("title",
  "<b>TTE HAIC THEN I On Demand Dynamic Strategy Arm</b><br><b>Patient Flow with Censoring</b><br>Clone-Censor-Weight Framework",
  S$text, 300, 20, 700, 50))

# Total
p1 <- paste0(p1, vtx("total",
  sprintf("<b>Total Cohort</b><br>N = %d", N),
  S$start, 520, 90, 200, 50))

# Row 1: Enter S1 / No Cycle3
p1 <- paste0(p1, vtx("enter_s1",
  sprintf("<b>Enter Stage 1</b><br>n = %d (%s)", enter_s1, pct_s(enter_s1)),
  S$stage, 330, 200, 180, 50))
p1 <- paste0(p1, vtx("no_c3",
  sprintf("<b>No Cycle 3 Data</b><br>n = %d (%s)", no_cycle3, pct_s(no_cycle3)),
  S$grey, 730, 200, 180, 50))
p1 <- paste0(p1, edg("e_total_s1", S$arrow, "total", "enter_s1"))
p1 <- paste0(p1, edg("e_total_noc3", S$arrow, "total", "no_c3"))

# Row 2: Stage 1
p1 <- paste0(p1, vtx("s1_box",
  sprintf("<b>Stage 1: pre-HAIC-3 Assessment</b><br><font style=\"font-size:10px\">Trigger: PVTT/Meta | PIV_bl &gt; 149.5 | PIV_pre3 &gt; 272.4<br>Exempt: AFP drop &gt; 42.9%%</font>"),
  S$stage, 280, 310, 280, 60))
p1 <- paste0(p1, edg("e_enter_s1", S$arrow2, "enter_s1", "s1_box"))

# Row 3: S1 outcomes
p1 <- paste0(p1, vtx("s1_trig",
  sprintf("<b>Triggered + Not Exempt</b><br>n = %d (%s)", s1_trig, pct_s(s1_trig)),
  S$trig, 100, 430, 200, 50))
p1 <- paste0(p1, vtx("s1_pass",
  sprintf("<b>Not Triggered / Exempt → S2</b><br>n = %d (%s)<br><font style=\"font-size:9px;color:#E8F8F5\">(Exempt: %d + No Trigger: %d)</font>",
    s1_to_s2, pct_s(s1_to_s2), s1_exempt, s1_no_trig),
  S$exempt, 430, 425, 220, 60))
p1 <- paste0(p1, edg("e_s1_trig", S$arrow2, "s1_box", "s1_trig"))
p1 <- paste0(p1, edg("e_s1_pass", S$arrow2, "s1_box", "s1_pass"))

# S1 censoring
p1 <- paste0(p1, vtx("s1_censor",
  sprintf("<font color=\"#16A085\"><b>✔ Comply: %d (%d%%)</b></font><br><font color=\"#E67E22\"><b>✘ Censor: %d (%d%%)</b></font><br><font style=\"font-size:8px;color:#999\">(no immune within 90d grace)</font>",
    s1_comply, round(s1_comply/(s1_comply+s1_censor)*100), s1_censor, round(s1_censor/(s1_comply+s1_censor)*100)),
  S$censor, 100, 500, 200, 55))
p1 <- paste0(p1, edg("e_s1_cen", S$dash, "s1_trig", "s1_censor"))

# Row 4: Stage 2
p1 <- paste0(p1, vtx("s2_box",
  sprintf("<b>Stage 2: pre-HAIC-5 Assessment</b><br><font style=\"font-size:10px\">AFP change &gt; -46.7%% (drop &lt; 46.7%%)?<br>Entered: %d</font>", s1_to_s2),
  S$stage, 390, 555, 280, 55))
p1 <- paste0(p1, edg("e_s1pass_s2", S$arrow2, "s1_pass", "s2_box"))

# Row 5: S2 outcomes
p1 <- paste0(p1, vtx("s2_trig",
  sprintf("<b>Triggered (AFP↓ &lt; 46.7%%)</b><br>n = %d (%s)", s2_trig, pct_s(s2_trig)),
  S$trig, 280, 670, 210, 50))
p1 <- paste0(p1, vtx("s2_pass",
  sprintf("<b>Pass (AFP↓ ≥ 46.7%%) → S3</b><br>n = %d (%s)", s2_pass, pct_s(s2_pass)),
  S$exempt, 570, 670, 210, 50))
p1 <- paste0(p1, edg("e_s2_trig", S$arrow2, "s2_box", "s2_trig"))
p1 <- paste0(p1, edg("e_s2_pass", S$arrow2, "s2_box", "s2_pass"))

# S2 censoring
p1 <- paste0(p1, vtx("s2_censor",
  sprintf("<font color=\"#16A085\"><b>✔ Comply: %d (%d%%)</b></font><br><font color=\"#E67E22\"><b>✘ Censor: %d (%d%%)</b></font><br><font style=\"font-size:8px;color:#999\">(no immune within 90d grace)</font>",
    s2_comply, round(s2_comply/(s2_comply+s2_censor)*100), s2_censor, round(s2_censor/(s2_comply+s2_censor)*100)),
  S$censor, 280, 740, 210, 55))
p1 <- paste0(p1, edg("e_s2_cen", S$dash, "s2_trig", "s2_censor"))

# Row 6: Stage 3
p1 <- paste0(p1, vtx("s3_box",
  sprintf("<b>Stage 3: Post 6th HAIC Assessment</b><br><font style=\"font-size:10px\">AFP↓ &gt; 87.9%% AND AFP &lt; 20?<br>Entered: %d</font>", s2_pass),
  S$stage, 500, 790, 280, 55))
p1 <- paste0(p1, edg("e_s2pass_s3", S$arrow2, "s2_pass", "s3_box"))

# Row 7: S3 outcomes
p1 <- paste0(p1, vtx("s3_trig",
  sprintf("<b>Not Exempt → Add Immune</b><br>n = %d (%s)", s3_trig, pct_s(s3_trig)),
  S$trig, 380, 910, 200, 50))
p1 <- paste0(p1, vtx("s3_exempt",
  sprintf("<b>Exempt (criteria met)</b><br>n = %d (%s)", s3_exempt, pct_s(s3_exempt)),
  S$exempt, 640, 910, 180, 50))
p1 <- paste0(p1, vtx("s3_nodata",
  sprintf("<b>No Stage 3 Data</b><br>n = %d (%s)", s3_nodata, pct_s(s3_nodata)),
  S$grey, 880, 910, 170, 50))
p1 <- paste0(p1, edg("e_s3_trig", S$arrow2, "s3_box", "s3_trig"))
p1 <- paste0(p1, edg("e_s3_exempt", S$arrow2, "s3_box", "s3_exempt"))
p1 <- paste0(p1, edg("e_s3_nodata", S$arrow2, "s3_box", "s3_nodata"))

# S3 censoring
p1 <- paste0(p1, vtx("s3_censor",
  sprintf("<font color=\"#16A085\"><b>✔ Comply: %d (%d%%)</b></font><br><font color=\"#E67E22\"><b>✘ Censor: %d (%d%%)</b></font><br><font style=\"font-size:8px;color:#999\">(no immune within 90d grace)</font>",
    s3_comply, round(s3_comply/(s3_comply+s3_censor)*100), s3_censor, round(s3_censor/(s3_comply+s3_censor)*100)),
  S$censor, 380, 980, 200, 55))
p1 <- paste0(p1, edg("e_s3_cen", S$dash, "s3_trig", "s3_censor"))

# Layer 2
p1 <- paste0(p1, vtx("l2_box",
  sprintf("<b>Layer 2: Post-HAIC AFP Monitoring</b><br>n = %d (%s)<br><font style=\"font-size:9px\">Condition A: AFP &gt; 20 | Condition B: nadir &lt; 20 &amp; rise &gt; 1.3</font>", l2_total, pct_s(l2_total)),
  S$layer2, 1120, 520, 290, 65))

# Feed edges into L2
p1 <- paste0(p1, edg("e_noc3_l2", S$feed, "no_c3", "l2_box", value = as.character(no_cycle3)))
p1 <- paste0(p1, edg("e_s3e_l2",
  paste0(S$feed, "exitX=0.5;exitY=1;exitDx=0;exitDy=0;"),
  "s3_exempt", "l2_box", value = as.character(s3_exempt),
  pts = matrix(c(730, 1030, 1080, 1030, 1080, 553), ncol = 2, byrow = TRUE)))
p1 <- paste0(p1, edg("e_s3n_l2",
  paste0(S$feed, "exitX=1;exitY=0.5;exitDx=0;exitDy=0;"),
  "s3_nodata", "l2_box", value = as.character(s3_nodata),
  pts = matrix(c(1090, 935, 1090, 560), ncol = 2, byrow = TRUE)))

# L2 sources annotation
p1 <- paste0(p1, vtx("l2_sources",
  sprintf("<b>Layer 2 Entry Sources:</b><br>No Cycle3: %d | S3 Exempt: %d | S3 NoData: %d", no_cycle3, s3_exempt, s3_nodata),
  S$l2src, 1140, 460, 250, 40))

# L2 outcomes
p1 <- paste0(p1, vtx("l2_trig",
  sprintf("<b>AFP Abnormal → Add Immune</b><br>n = %d (%s)", l2_trig, pct_s(l2_trig)),
  S$trig, 1080, 650, 210, 50))
p1 <- paste0(p1, vtx("l2_never",
  sprintf("<b>Never Triggered</b><br>n = %d (%s)", l2_never, pct_s(l2_never)),
  S$never, 1340, 650, 180, 50))
p1 <- paste0(p1, edg("e_l2_trig", S$arrow2, "l2_box", "l2_trig"))
p1 <- paste0(p1, edg("e_l2_never", S$arrow2, "l2_box", "l2_never"))

# L2 censoring boxes
p1 <- paste0(p1, vtx("l2_censor",
  sprintf("<font color=\"#16A085\"><b>✔ Comply: %d (%d%%)</b></font><br><font color=\"#E67E22\"><b>✘ Censor: %d (%d%%)</b></font><br><font style=\"font-size:8px;color:#999\">(no immune within 90d grace)</font>",
    l2_comply, round(l2_comply/(l2_comply+l2_censor)*100), l2_censor, round(l2_censor/(l2_comply+l2_censor)*100)),
  S$censor, 1070, 720, 220, 55))
p1 <- paste0(p1, edg("e_l2_cen", S$dash, "l2_trig", "l2_censor"))

p1 <- paste0(p1, vtx("never_censor",
  sprintf("<font color=\"#E67E22\"><b>Excluded (immune early): %d</b></font><br>Enter Dyn arm: %d<br><font color=\"#E67E22\"><b>  ✘ Censor (got immune later): %d</b></font><br><font color=\"#16A085\"><b>  ✔ Remain (truly no immune): %d</b></font>",
    never_excluded, never_in_dyn, never_censor_later, never_remain),
  paste0(S$censor, "align=left;spacingLeft=10;"), 1310, 720, 240, 75))
p1 <- paste0(p1, edg("e_nev_cen", S$dash, "l2_never", "never_censor"))

# Summary boxes
p1 <- paste0(p1, vtx("sum_trig",
  sprintf("<b>Total Triggered: %d (%s)</b><br>S1: %d | S2: %d | S3: %d | L2: %d",
    total_trig, pct_s(total_trig), s1_trig, s2_trig, s3_trig, l2_trig),
  S$trig, 280, 1070, 320, 55))
p1 <- paste0(p1, vtx("sum_censor",
  sprintf("<b>Total Artificially Censored: %d</b><br>Grace non-comply: %d | Never-trig got immune: %d+%d",
    s1_censor+s2_censor+s3_censor+l2_censor+never_excluded+never_censor_later,
    s1_censor+s2_censor+s3_censor+l2_censor, never_excluded, never_censor_later),
  S$orange, 660, 1070, 360, 55))
p1 <- paste0(p1, vtx("sum_final",
  sprintf("<b>Dynamic Arm Final: n = %d</b><br>Events: %d | Art.Censored: %d | Alive uncensored: %d",
    dyn_arm_n, dyn_events, total_art_censor, dyn_arm_n - dyn_events - total_art_censor + total_art_censor),
  S$result, 1080, 1070, 360, 55))

# ── Early Combo arm (right half, x offset = 1600) ──
ox <- 1600

p1 <- paste0(p1, vtx("ec_title",
  "<b>TTE HAIC+I_CONC Early Combination Strategy Arm</b><br><b>Patient Flow with Censoring</b><br>Clone-Censor-Weight Framework",
  S$text, ox + 100, 20, 700, 50))

p1 <- paste0(p1, vtx("ec_total",
  sprintf("<b>Total Cohort (Cloned into Early Combo Arm)</b><br>N = %d", N),
  S$start, ox + 250, 100, 360, 50))
p1 <- paste0(p1, vtx("ec_rule",
  "<b>Early Combination Rule</b><br><font style=\"font-size:11px\">Add immunotherapy within <b>14 days</b> of first HAIC</font>",
  "rounded=1;whiteSpace=wrap;html=1;fillColor=#2980B9;strokeColor=#2471A3;fontColor=#FFFFFF;fontSize=12;",
  ox + 250, 210, 360, 55))
p1 <- paste0(p1, edg("ec_e1", S$arrow, "ec_total", "ec_rule"))

p1 <- paste0(p1, vtx("ec_comply",
  sprintf("<b>Compliant: Added Immune ≤ 14 days</b><br>n = %d (%s)<br><font style=\"font-size:9px;color:#FADBD8\">Median: 0 days (same day as first HAIC)</font>",
    ec_comply, pct_s(ec_comply)),
  "rounded=1;whiteSpace=wrap;html=1;fillColor=#C0392B;strokeColor=#A93226;fontColor=#FFFFFF;fontSize=12;",
  ox + 110, 340, 290, 65))
p1 <- paste0(p1, vtx("ec_censored",
  sprintf("<b>Not Compliant: No Immune ≤ 14 days</b><br>n = %d (%s)<br><font style=\"font-size:10px\">→ Artificially Censored at Day 14</font>",
    ec_censored, pct_s(ec_censored)),
  "rounded=1;whiteSpace=wrap;html=1;fillColor=#E67E22;strokeColor=#CA6F1E;fontColor=#FFFFFF;fontSize=12;",
  ox + 445, 340, 300, 65))
p1 <- paste0(p1, edg("ec_e2", S$arrow, "ec_rule", "ec_comply"))
p1 <- paste0(p1, edg("ec_e3", S$arrow, "ec_rule", "ec_censored"))

p1 <- paste0(p1, vtx("ec_comply_detail",
  sprintf("<font color=\"#16A085\"><b>Remain in Early Combo Arm</b></font><br>n = %d", ec_comply),
  S$comply_detail, ox + 120, 470, 270, 95))
p1 <- paste0(p1, edg("ec_e4", S$arrow2, "ec_comply", "ec_comply_detail"))

p1 <- paste0(p1, vtx("ec_censor_detail",
  sprintf("<b>Censoring Breakdown (n = %d):</b><br><br><font color=\"#6C757D\"><b>No immunotherapy at all: %d (%.1f%%)</b></font><br><font style=\"font-size:9px;color:#888\">HAIC-only patients, never received immune</font><br><br><font color=\"#E67E22\"><b>Immune added after Day 14: %d (%.1f%%)</b></font><br><font style=\"font-size:9px;color:#888\">HAIC-then-Immune patients (delayed addition)</font>",
    ec_censored, ec_no_immune, ec_no_immune/ec_censored*100, ec_after_14d, ec_after_14d/ec_censored*100),
  S$cendet, ox + 445, 460, 300, 130))
p1 <- paste0(p1, edg("ec_e5", "edgeStyle=orthogonalEdgeStyle;html=1;strokeColor=#E67E22;strokeWidth=1.5;dashed=1;",
  "ec_censored", "ec_censor_detail"))

# Close page 1
p1 <- paste0(p1, '      </root>
    </mxGraphModel>')

# =============================================================================
# Build Page 2: Early Combo standalone (simpler view)
# =============================================================================
p2_cells <- paste0(
  vtx("ec2_title",
    "<b>TTE Early Combination Strategy Arm:</b><br><b>Patient Flow with Censoring</b><br>Clone-Censor-Weight Framework",
    S$text, 200, 20, 700, 50),
  vtx("ec2_total",
    sprintf("<b>Total Cohort (Cloned into Early Combo Arm)</b><br>N = %d", N),
    S$start, 350, 100, 360, 50),
  vtx("ec2_rule",
    "<b>Early Combination Rule</b><br><font style=\"font-size:11px\">Add immunotherapy within <b>14 days</b> of first HAIC</font>",
    "rounded=1;whiteSpace=wrap;html=1;fillColor=#2980B9;strokeColor=#2471A3;fontColor=#FFFFFF;fontSize=12;",
    350, 210, 360, 55),
  edg("ec2_e1", S$arrow, "ec2_total", "ec2_rule"),
  vtx("ec2_comply",
    sprintf("<b>Compliant: Added Immune ≤ 14 days</b><br>n = %d (%s)<br><font style=\"font-size:9px;color:#FADBD8\">Median: 0 days (same day as first HAIC)</font>",
      ec_comply, pct_s(ec_comply)),
    "rounded=1;whiteSpace=wrap;html=1;fillColor=#C0392B;strokeColor=#A93226;fontColor=#FFFFFF;fontSize=12;",
    210, 340, 290, 65),
  vtx("ec2_censored",
    sprintf("<b>Not Compliant: No Immune ≤ 14 days</b><br>n = %d (%s)<br><font style=\"font-size:10px\">→ Artificially Censored at Day 14</font>",
      ec_censored, pct_s(ec_censored)),
    "rounded=1;whiteSpace=wrap;html=1;fillColor=#E67E22;strokeColor=#CA6F1E;fontColor=#FFFFFF;fontSize=12;",
    545, 340, 300, 65),
  edg("ec2_e2", S$arrow, "ec2_rule", "ec2_comply"),
  edg("ec2_e3", S$arrow, "ec2_rule", "ec2_censored"),
  vtx("ec2_comply_detail",
    sprintf("<font color=\"#16A085\"><b>Remain in Early Combo Arm</b></font><br>n = %d", ec_comply),
    S$comply_detail, 220, 470, 270, 95),
  edg("ec2_e4", S$arrow2, "ec2_comply", "ec2_comply_detail"),
  vtx("ec2_censor_detail",
    sprintf("<b>Censoring Breakdown (n = %d):</b><br><br><font color=\"#6C757D\"><b>No immunotherapy at all: %d (%.1f%%)</b></font><br><font style=\"font-size:9px;color:#888\">HAIC-only patients, never received immune</font><br><br><font color=\"#E67E22\"><b>Immune added after Day 14: %d (%.1f%%)</b></font><br><font style=\"font-size:9px;color:#888\">HAIC-then-Immune patients (delayed addition)</font>",
      ec_censored, ec_no_immune, ec_no_immune/ec_censored*100, ec_after_14d, ec_after_14d/ec_censored*100),
    S$cendet, 545, 460, 300, 130)
)
p2_cells <- paste0(p2_cells, edg("ec2_e5", "edgeStyle=orthogonalEdgeStyle;html=1;strokeColor=#E67E22;strokeWidth=1.5;dashed=1;",
  "ec2_censored", "ec2_censor_detail"))

p2 <- sprintf('    <mxGraphModel dx="938" dy="703" grid="1" gridSize="10" guides="1" tooltips="1" connect="1" arrows="1" fold="1" page="1" pageScale="1" pageWidth="1200" pageHeight="1000" math="0" shadow="0">
      <root>
        <mxCell id="0" />
        <mxCell id="1" parent="0" />
%s      </root>
    </mxGraphModel>', p2_cells)

# =============================================================================
# Assemble and write
# =============================================================================
drawio <- sprintf('<mxfile host="Electron" modified="%s" version="26.0.0" pages="2">
  <diagram id="tte-flow" name="TTE Dynamic Strategy Flow">
%s
  </diagram>
  <diagram id="early-combo" name="Early Combination Arm">
%s
  </diagram>
</mxfile>\n', format(Sys.time(), "%Y-%m-%dT%H:%M:%S.000Z"), p1, p2)

out_file <- file.path(out_dir, "TTE_Dynamic_Flow_Censoring.drawio")
writeLines(drawio, out_file)
cat(sprintf("Draw.io file written: %s\n", out_file))
cat("Open in VS Code (draw.io extension) or draw.io web/desktop app.\n")
