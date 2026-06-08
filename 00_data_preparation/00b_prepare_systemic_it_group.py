#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the 8th group (Systemic I+T) and append to analysis_ready.csv.

Systemic I+T = systemic immunotherapy + anti-angiogenic / targeted therapy, NO induction
HAIC. Concurrent I+T only (系统类型 == 靶免). Deterministic; no manual edits.

Missing covariates are imputed with the SAME recipe as the main cohort
(HAIC_NO_TACE_4_TIDY/scripts/impute_baseline.py): MICE (IterativeImputer +
RandomForestRegressor, max_iter=10, seed 42) on continuous vars, then derived fields
(log_afp/pivka, afp_high/pivka_high, tumor_size_category, ALBI, Child-Pugh, BCLC)
recomputed from imputed values. Imputation is fit on the new group ALONE so the main
7-group cohort is untouched.
"""
import os
import numpy as np
import pandas as pd
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.ensemble import RandomForestRegressor

BASE = "/Users/yqj/Nutstore Files/我的坚果云/Liver_tumor_big_data/FIRST_LINE_HAIC_2025-12-30"
RAW = os.path.join(BASE, "补充单纯系统治疗患者", "systemic_0519.xlsx")
DATA = os.path.join(BASE, "HAIC_NO_TACE_4_TIDY", "update_group_7", "data")
QCDIR = os.path.join(BASE, "HAIC_NO_TACE_4_TIDY", "update_group_7", "results", "_8group_qc")
RANDOM_STATE = 42
os.makedirs(QCDIR, exist_ok=True)

main = pd.read_csv(os.path.join(DATA, "analysis_ready.csv"))
raw = pd.read_excel(RAW, sheet_name="Sheet1")


def s(x):
    return "" if pd.isna(x) else str(x).strip()


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


# ── cohort selection (waterfall) ───────────────────────────────────────────
wf = []
df = raw[raw["系统类型"].apply(s) == "靶免"].copy()
wf.append(("靶免 I+T concurrent", len(df)))
df = df[df["是否需要删除"].apply(lambda v: s(v) not in ("1", "1.0"))]
wf.append(("- exclude 是否需要删除=1", len(df)))
df = df[df["治疗后HAIC"].apply(lambda v: s(v) not in ("1", "1.0"))]
wf.append(("- exclude 治疗后HAIC=1", len(df)))
df = df[df["是否行介入治疗"].apply(lambda v: s(v) != "是")]
wf.append(("- exclude 是否行介入治疗=是", len(df)))
# invalid follow-up time (death/last-contact before treatment → os_months <= 0)
df = df[df["OS_months"].apply(num) > 0]
wf.append(("- exclude os_months<=0 (invalid)", len(df)))
df = df.reset_index(drop=True)

# ── raw mapping (continuous kept raw; derived left for post-MICE recompute) ─
g = {}
g["patient_id"] = "SYS_" + df["ID"].apply(s)
g["first_treatment_date"] = df["治疗基线时间点_首次目标治疗日期"]
g["first_haic_date"] = np.nan
g["haic_episodes"] = 0
g["immune_episodes"] = np.nan
g["has_immunotherapy"] = 1
g["has_target_therapy"] = 1
g["has_tace"] = 0
g["first_immune_date"] = df["首次免疫开始日期"]
g["first_immune_drug"] = df["首次免疫药物"]
g["days_haic_to_immune"] = np.nan
g["first_target_drug"] = df["首次靶向药药物"]
g["days_haic_to_target"] = np.nan
g["age"] = df["Age"].apply(num)
g["sex"] = df["sex"].map({"男": "Male", "女": "Female"})
_hbv = df["Etiology_HBV"].apply(num) == 1   # column is float (has NaN); compare numerically
g["etiology"] = np.where(_hbv, "HBV", "Other")
g["hbsag"] = np.where(_hbv, "Positive", "Negative")
for tgt, src in [("alt", "ALT.(U/L)"), ("ast", "AST.(U/L)"), ("tbil", "TBIL.(umol/L)"),
                 ("alb", "ALB(g/L)"), ("pt", "PT.(s)"), ("inr", "INR"),
                 ("plt", "PLT.(×10^9/L)"), ("hb", "Hb.(g/L)"), ("wbc", "WBC.(×10^9/L)"),
                 ("neut", "NEUT.(×10^9/L)"), ("lymph", "Lymphocytes.(×10^9/L)"),
                 ("nlr", "NLR"), ("plr", "PLR")]:
    g[tgt] = df[src].apply(num)
for c in ["dbil", "pta", "aptt", "fbg", "creatinine", "mono"]:
    g[c] = np.nan          # not collected in this cohort (not IPTW covariates)
g["afp"] = df["AFP"].apply(num)
g["afp_source"] = "lab"
g["pivka"] = df["PIVKA-II"].apply(num)
g["pivka_source"] = "lab"
# derived (recomputed after MICE): log_afp, log_pivka, afp_high, pivka_high,
# albi_score, albi_grade, child_pugh_score, child_pugh_grade, bclc_stage,
# tumor_size_category — leave NaN here.
for c in ["log_afp", "log_pivka", "afp_high", "pivka_high", "albi_score", "albi_grade",
          "child_pugh_score", "child_pugh_grade", "bclc_stage", "tumor_size_category"]:
    g[c] = np.nan
g["tumor_max_diameter_cm"] = df["Tumor_Diameter"].apply(num)
tn3 = df["Tumor_Number_3"].apply(s)
tn = df["Tumor_Number"].apply(s)


def tcat(a, b):
    if a == "1":
        return ">3"
    if b in ("0", "1"):
        return "Solitary"
    if b in ("2", "3"):
        return "2-3"
    return "Solitary"


g["tumor_count_category"] = [tcat(a, b) for a, b in zip(tn3, tn)]
vpmap = {"Vp3-Vp4": "Vp3/4", "Vp1-Vp2": "Vp1/2", "Absent": "Absent"}


def pvtt(row):
    v = vpmap.get(s(row["Vp_class"]))
    if v is not None:
        return v
    return "Vp3/4" if s(row["PVTT"]) in ("1", "1.0") else "Absent"


g["pvtt_classification"] = [pvtt(r) for _, r in df.iterrows()]
g["hvtt"] = np.where(df["HVTT"].apply(num) == 1, "Yes", "No")
g["ivc_or_ra_thrombus"] = np.where(
    (df["IVCTT"].apply(num) == 1) | (df["RATT"].apply(num) == 1), "Yes", "No")
g["distant_metastasis"] = np.where(df["EHS"].apply(num) == 1, "Yes", "No")
g["lymph_node_metastasis"] = np.where(
    df["EHS_腹部淋巴结转移部位"].apply(lambda v: s(v) not in ("", "None")), "Yes", "No")
asc = df["腹水"].apply(s)
g["ascites"] = asc.map({"无": "Absent", "少量": "Mild", "中量": "Moderate-Severe",
                        "大量": "Moderate-Severe"})
g["ascites_score"] = asc.map({"无": 0, "少量": 1, "中量": 2, "大量": 2})
g["varices"] = np.where(df["varices"].apply(lambda v: s(v) not in ("", "None")), "Yes", "No")
g["hepatic_encephalopathy"] = np.where(df["肝性脑病"].apply(s) == "1", "Yes", "No")
g["main_group"] = "Systemic_I+T"
# match main-cohort encoding ('Yes'/'No' strings)
g["death_status"] = np.where(df["Death_status"].apply(num) == 1, "Yes", "No")
g["death_date"] = df["死亡_死亡日期"]
g["event_date"] = df["死亡_死亡日期"].where(df["Death_status"].apply(s) == "1",
                                          df["死亡_存活判断时间"])
g["os_days"] = df["OS_cal"].apply(num)
g["os_months"] = df["OS_months"].apply(num)

new = pd.DataFrame(g)
miss_pre = {c: int(new[c].isna().sum()) for c in
            ["afp", "pivka", "tumor_max_diameter_cm", "tbil", "alb", "plt", "age", "nlr"]}

# ── MICE imputation (mirror impute_baseline.py; new group alone) ───────────
CONT_AVAIL = ["alt", "ast", "tbil", "alb", "pt", "inr", "plt", "hb", "wbc",
              "neut", "lymph", "nlr", "plr", "tumor_max_diameter_cm", "afp", "pivka"]
AUX_CAT = {
    "etiology": {"HBV": 1, "HCV": 2, "Other": 0},
    "pvtt_classification": {"Absent": 0, "Vp1/2": 1, "Vp3/4": 2},
    "hvtt": {"No": 0, "Yes": 1},
    "ivc_or_ra_thrombus": {"No": 0, "Yes": 1},
    "distant_metastasis": {"No": 0, "Yes": 1},
    "lymph_node_metastasis": {"No": 0, "Yes": 1},
    "ascites": {"Absent": 0, "Mild": 1, "Moderate-Severe": 2},
    "tumor_count_category": {"Solitary": 1, "2-3": 2, ">3": 3},
}
mice_df = new[CONT_AVAIL + ["age", "ascites_score"]].copy()
for col, mapping in AUX_CAT.items():
    mice_df[col + "_enc"] = new[col].map(mapping)
imputer = IterativeImputer(
    estimator=RandomForestRegressor(n_estimators=50, max_depth=8,
                                    random_state=RANDOM_STATE, n_jobs=-1),
    max_iter=10, random_state=RANDOM_STATE, verbose=0)
imp = pd.DataFrame(imputer.fit_transform(mice_df), columns=mice_df.columns, index=new.index)
for col in CONT_AVAIL + ["age"]:
    m = new[col].isna()
    if m.sum():
        new.loc[m, col] = imp.loc[m, col]

# ── recompute derived (exact impute_baseline.py Step 4) ───────────────────
new["log_afp"] = np.log1p(new["afp"].clip(lower=0))
new["log_pivka"] = np.log1p(new["pivka"].clip(lower=0))
new["afp_high"] = new["afp"].apply(lambda x: "Yes" if x > 400 else "No")
new["pivka_high"] = new["pivka"].apply(lambda x: "Yes" if x > 8000 else "No")
new["tumor_size_category"] = new["tumor_max_diameter_cm"].apply(
    lambda d: np.nan if pd.isna(d) else ("≤10cm" if d <= 10 else ">10cm"))


def calc_albi(alb, tbil):
    if pd.isna(alb) or pd.isna(tbil):
        return np.nan
    return round(-0.085 * alb + 0.66 * np.log10(tbil * 0.05848), 3)


def albi_grade_fn(score):
    if pd.isna(score):
        return np.nan
    if score <= -2.60:
        return 1
    if score <= -1.39:
        return 2
    return 3


new["albi_score"] = new.apply(lambda r: calc_albi(r["alb"], r["tbil"]), axis=1)
new["albi_grade"] = new["albi_score"].apply(albi_grade_fn)


def cp_tbil(v):
    return np.nan if pd.isna(v) else (1 if v < 34 else (2 if v <= 51 else 3))


def cp_alb(v):
    return np.nan if pd.isna(v) else (1 if v > 35 else (2 if v >= 28 else 3))


def cp_inr(v):
    return np.nan if pd.isna(v) else (1 if v < 1.7 else (2 if v <= 2.3 else 3))


_cp = (new["tbil"].apply(cp_tbil) + new["alb"].apply(cp_alb) + new["inr"].apply(cp_inr)
       + new["ascites"].map({"Absent": 1, "Mild": 2, "Moderate-Severe": 3}) + 1)
new["child_pugh_score"] = _cp.round(0).astype("Int64")
new["child_pugh_grade"] = new["child_pugh_score"].apply(
    lambda x: np.nan if pd.isna(x) else ("A" if x <= 6 else ("B" if x <= 9 else "C")))


def calc_bclc(row):
    cp, tc, td = row["child_pugh_grade"], row["tumor_count_category"], row["tumor_max_diameter_cm"]
    if cp == "C":
        return "D"
    vascular = (row["pvtt_classification"] in ("Vp1/2", "Vp3/4")) or (row["hvtt"] == "Yes") \
        or (row["ivc_or_ra_thrombus"] == "Yes")
    extrahepatic = (row["distant_metastasis"] == "Yes") or (row["lymph_node_metastasis"] == "Yes")
    if vascular or extrahepatic:
        return "C"
    if tc == "Solitary" and pd.notna(td) and td <= 2.0 and cp == "A":
        return "0"
    if tc == "Solitary":
        return "A"
    if tc == "2-3" and pd.notna(td) and td <= 3.0:
        return "A"
    if tc in ("2-3", ">3"):
        return "B"
    return np.nan


new["bclc_stage"] = new.apply(calc_bclc, axis=1)

# ── decimal rounding (impute_baseline.py Step 5b) ─────────────────────────
for col in ["plt", "hb"]:
    new[col] = new[col].round(0).astype("Int64")
for col in ["alt", "ast", "tbil", "alb", "pt", "neut", "lymph", "wbc",
            "tumor_max_diameter_cm"]:
    new[col] = new[col].round(1)
for col in ["inr", "nlr", "plr", "albi_score", "log_afp", "log_pivka", "os_months"]:
    new[col] = new[col].round(2)
for col in ["afp", "pivka"]:
    new[col] = new[col].round(1)

# ── assemble + save ───────────────────────────────────────────────────────
new = new.reindex(columns=main.columns)
out = pd.concat([main, new], ignore_index=True)
out.to_csv(os.path.join(DATA, "analysis_ready_8group.csv"), index=False)

# ── QC report ─────────────────────────────────────────────────────────────
cov = ["afp", "pivka", "tumor_max_diameter_cm", "tumor_count_category", "pvtt_classification",
       "hvtt", "ivc_or_ra_thrombus", "distant_metastasis", "lymph_node_metastasis", "ascites",
       "varices", "albi_grade", "tbil", "alb", "plt", "age", "nlr"]
lines = ["# 8-group data QC (Systemic I+T)", "", "## Exclusion waterfall"]
for k, v in wf:
    lines.append(f"- {k}: {v}")
lines += ["", "## Pre-imputation missingness (new group)"]
for c, v in miss_pre.items():
    lines.append(f"- {c}: {v}/{len(new)}")
lines += ["", "## Output",
          f"- total rows: {len(out)} (main {len(main)} + new {len(new)})",
          f"- columns: {len(out.columns)}",
          f"- new main_group=='Systemic_I+T': {(out['main_group'] == 'Systemic_I+T').sum()}",
          f"- events in new group: {(new['death_status'] == 'Yes').sum()}",
          "", "## Post-imputation covariate completeness (new group)"]
miss = new[cov].isna().sum()
for c in cov:
    lines.append(f"- {c}: missing {int(miss[c])}/{len(new)}")
lines += ["", f"## IPTW complete-cases (new group): {new[cov].dropna().shape[0]}/{len(new)}",
          "", "## Flags",
          "- sex uniformly 'Male' in source (data-quality artifact; NOT an IPTW covariate)",
          "- varices uniformly recorded absent in source",
          "- dbil/creatinine/mono/pta/aptt/fbg not collected (left missing; not IPTW covariates)",
          "- MICE fit on new group alone (seed 42); main 7-group cohort untouched"]
report = "\n".join(lines)
with open(os.path.join(QCDIR, "data_qc_report.md"), "w") as f:
    f.write(report)
print(report)
print("\nWROTE", os.path.join(DATA, "analysis_ready_8group.csv"))
