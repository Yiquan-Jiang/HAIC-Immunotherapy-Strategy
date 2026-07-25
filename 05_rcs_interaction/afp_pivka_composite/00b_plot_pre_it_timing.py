#!/usr/bin/env python3
"""
绘制 pre-IT 采血时间点相对免疫治疗与 HAIC 的间隔分布（补充材料图）。

回答两个问题：
  (a) pre-IT 采血距免疫治疗开始多少天？
  (b) pre-IT 采血距最近一次 HAIC 多少天？（负值 = HAIC 之前，0 = HAIC 当天输注前）

数据来源：
  - 00_swimmer_plot_events.csv  → 首次免疫/靶向日期、全部 HAIC 治疗日期
  - composite_*_cohort.csv      → 由 00a_extract_pre_it_labs.py 写入的
                                  pre_it_lab_date / pre_it_source

纳入：4 个含免疫治疗的 composite 队列中 trt=1 且成功匹配到 pre-IT 检查的患者，
按 patient_id 去重（N≈2105：序贯 363 + 并发 1742）。

上游：00a_extract_pre_it_labs.py
"""

import os
import warnings

warnings.filterwarnings("ignore")
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, "..", "..", ".."))

SWIMMER_CSV = os.path.join(
    PROJECT_ROOT, "data", "publication_cohorts", "00_swimmer_plot_events.csv"
)
OUT_DIR = os.path.join(SCRIPT_DIR, "output_pre_it_timing")

# composite 队列 → 该队列的治疗组标记列
COHORTS = {
    "composite_THEN_I_cohort.csv": "trt_haic_then_i",
    "composite_THEN_IT_cohort.csv": "trt_haic_then_it",
    "composite_I_CONC_cohort.csv": "trt_haic_i_conc",
    "composite_IT_CONC_cohort.csv": "trt_haic_it_conc",
}
MATCHED_SOURCES = ["pre_it_matched", "baseline_concurrent"]

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",  # SVG 内文字可编辑
        "pdf.fonttype": 42,  # PDF 内文字可编辑
        "font.size": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
    }
)

BLUE, RED, TEAL, INK = "#3a6ea5", "#b2182b", "#1b7a5a", "#222222"


# ── 数据组装 ─────────────────────────────────────────────────────────────────


def load_intervals():
    """返回每位患者的 lab→IT 与 lab→最近一次 HAIC 的间隔（天）。"""
    sw = pd.read_csv(SWIMMER_CSV)
    sw["start_date"] = pd.to_datetime(sw["start_date"], errors="coerce")
    sw = sw[sw["start_date"].notna()]

    haic = sw[sw["treatment_category"].isin(["HAIC", "HAIC+TACE"])]
    haic_dates = haic.groupby("patient_id")["start_date"].apply(
        lambda s: sorted(set(s))
    )

    it = sw[sw["treatment_category"].isin(["Immunotherapy", "Targeted Therapy"])]
    first_it = it.groupby("patient_id")["start_date"].min()

    rows = []
    for csv_name, trt_col in COHORTS.items():
        path = os.path.join(SCRIPT_DIR, csv_name)
        if not os.path.exists(path):
            print(f"  WARNING: {csv_name} not found — skipped")
            continue
        df = pd.read_csv(path, low_memory=False)
        if "first_it_date" not in df.columns:
            df["first_it_date"] = df.get("first_it_date_x")
        keep = (df[trt_col] == 1) & df["pre_it_source"].isin(MATCHED_SOURCES)
        rows.append(df.loc[keep, ["patient_id", "pre_it_lab_date"]])

    p = pd.concat(rows).drop_duplicates("patient_id")
    p["pre_it_lab_date"] = pd.to_datetime(p["pre_it_lab_date"], errors="coerce")
    p = p[p["pre_it_lab_date"].notna()].copy()
    p["first_it_date"] = p["patient_id"].map(first_it)

    # 距免疫/靶向开始的天数
    p["lab_to_it"] = (p["first_it_date"] - p["pre_it_lab_date"]).dt.days

    # 距最近一次 HAIC 的天数：正=HAIC 之后，0=HAIC 当天，负=所有 HAIC 之前
    def rel_haic(pid, lab_date):
        dates = haic_dates.get(pid)
        if not dates:
            return np.nan
        prior = [d for d in dates if d <= lab_date]
        if prior:
            return (lab_date - max(prior)).days
        return -(min(dates) - lab_date).days

    p["d_after_haic"] = [
        rel_haic(r.patient_id, r.pre_it_lab_date) for r in p.itertuples()
    ]
    return p.dropna(subset=["lab_to_it", "d_after_haic"])


# ── 绘图 ─────────────────────────────────────────────────────────────────────


def make_figure(p):
    n = len(p)
    fig = plt.figure(figsize=(172 / 25.4, 70 / 25.4))
    gs = fig.add_gridspec(
        1, 2, width_ratios=[1.0, 1.5], wspace=0.30,
        left=0.068, right=0.985, bottom=0.30, top=0.82,
    )

    # -- a: 采血 → 免疫治疗 --
    ax_a = fig.add_subplot(gs[0, 0])
    days = list(range(0, 11))
    pct = [(p.lab_to_it == d).mean() * 100 for d in days]
    pct.append((p.lab_to_it > 10).mean() * 100)
    xp = np.arange(len(pct))
    ax_a.bar(xp, pct, width=0.78, color=BLUE, edgecolor="white", linewidth=0.4, zorder=3)
    for i in range(2):
        ax_a.text(xp[i], pct[i] + 1.0, f"{pct[i]:.0f}%", ha="center", va="bottom",
                  fontsize=6.3, color=INK, fontweight="bold")
    ax_a.set_xticks(xp)
    ax_a.set_xticklabels([str(d) for d in days] + [">10"], fontsize=6.2)
    ax_a.set_xlabel("Days from blood draw to\nimmunotherapy initiation", labelpad=2)
    ax_a.set_ylabel("Patients (%)")
    ax_a.set_ylim(0, 50)

    cum = np.cumsum(pct)
    ax_cum = ax_a.twinx()
    ax_cum.plot(xp, cum, color=RED, lw=1.3, marker="o", ms=2.4, zorder=4)
    ax_cum.set_ylim(0, 105)
    ax_cum.set_ylabel("Cumulative (%)", color=RED, labelpad=1)
    ax_cum.tick_params(axis="y", colors=RED, pad=1)
    ax_cum.spines["right"].set_visible(True)
    ax_cum.spines["right"].set_color(RED)
    ax_cum.spines["top"].set_visible(False)
    ax_cum.annotate(f"≤3 d: {cum[3]:.0f}%", xy=(3, cum[3]), xytext=(3.6, cum[3] - 26),
                    fontsize=6.1, color=RED,
                    arrowprops=dict(arrowstyle="-", color=RED, lw=0.6))
    ax_a.set_title("Blood draw immediately precedes immunotherapy", fontsize=7.2, pad=4)
    ax_a.text(-0.22, 1.13, "a", transform=ax_a.transAxes, fontsize=11,
              fontweight="bold", va="top")

    # -- b: 采血 → 最近一次 HAIC --
    ax_b = fig.add_subplot(gs[0, 1])
    bins = [
        ("≤−15", -10**9, -15), ("−14 to −8", -14, -8), ("−7", -7, -7), ("−6", -6, -6),
        ("−5", -5, -5), ("−4", -4, -4), ("−3", -3, -3), ("−2", -2, -2), ("−1", -1, -1),
        ("0", 0, 0), ("1 to 7", 1, 7), ("8 to 14", 8, 14), ("15 to 21", 15, 21),
        ("22 to 28", 22, 28), ("29 to 42", 29, 42), ("43 to 90", 43, 90),
        (">90", 91, 10**9),
    ]
    xb = np.arange(len(bins))
    vals = np.array(
        [((p.d_after_haic >= lo) & (p.d_after_haic <= hi)).sum() / n * 100
         for _, lo, hi in bins]
    )
    i0 = [b[0] for b in bins].index("0")
    cols = ["#9fb8ce"] * i0 + [BLUE] + ["#7fb3a0"] * (len(bins) - i0 - 1)
    ax_b.bar(xb, vals, width=0.82, color=cols, edgecolor="white", lw=0.4, zorder=3)
    ax_b.axvline(i0 - 0.5, color="0.35", lw=0.9, ls=(0, (3, 2)), zorder=4)
    ax_b.text(xb[i0 - 1] - 0.1, vals[i0 - 1] + 1.2, f"{vals[i0 - 1]:.0f}%",
              ha="center", fontsize=6.3, color=INK, fontweight="bold")
    ax_b.text(xb[i0] + 0.25, vals[i0] + 1.2, f"{vals[i0]:.0f}%",
              ha="center", fontsize=6.3, color=INK, fontweight="bold")
    ax_b.set_xticks(xb)
    ax_b.set_xticklabels([b[0] for b in bins], fontsize=5.8, rotation=45, ha="right")
    ax_b.set_xlabel(
        "Days from nearest HAIC session to blood draw\n"
        "(negative = drawn before HAIC · 0 = HAIC day, pre-infusion)", labelpad=1,
    )
    ax_b.set_ylabel("Patients (%)")
    ax_b.set_ylim(0, 52)

    pre, post = vals[:i0].sum(), vals[i0 + 1:].sum()
    ax_b.annotate("", xy=(i0 - 0.58, 24.0), xytext=(-0.45, 24.0),
                  arrowprops=dict(arrowstyle="-", color="0.5", lw=0.9))
    ax_b.text(4.0, 24.8, f"before HAIC — {pre:.0f}%", ha="center", va="bottom",
              fontsize=6.0, color="0.42")
    ax_b.annotate("", xy=(len(bins) - 0.55, 13.5), xytext=(i0 + 0.6, 13.5),
                  arrowprops=dict(arrowstyle="-", color=TEAL, lw=0.9))
    ax_b.text(13.0, 14.3, f"after HAIC — {post:.0f}%", ha="center", va="bottom",
              fontsize=6.0, color=TEAL)
    ax_b.set_title(
        f"Relative to HAIC, the draw clusters at the HAIC visit\n(all patients, n={n:,})",
        fontsize=7.2, pad=4,
    )
    ax_b.text(-0.145, 1.13, "b", transform=ax_b.transAxes, fontsize=11,
              fontweight="bold", va="top")

    return fig, pct, cum, vals, i0


def main():
    print("=" * 70)
    print("Pre-IT blood-draw timing figure")
    print("=" * 70)

    p = load_intervals()
    print(f"\n患者数: {len(p)}")

    fig, pct, cum, vals, i0 = make_figure(p)

    os.makedirs(OUT_DIR, exist_ok=True)
    stem = os.path.join(OUT_DIR, "FigS_preIT_timing")
    for ext in ["png", "pdf", "svg", "tiff"]:
        fig.savefig(f"{stem}.{ext}",
                    dpi=600 if ext in ("png", "tiff") else None,
                    bbox_inches="tight")
    print(f"已保存: {stem}.{{png,pdf,svg,tiff}}  (png/tiff 600 dpi)")

    within_week = ((p.d_after_haic >= -7) & (p.d_after_haic <= 0)).mean() * 100
    print("\n关键数值:")
    print(f"  (a) 采血当天开始免疫治疗 {pct[0]:.0f}%, ≤3 天 {cum[3]:.0f}%, "
          f"≤7 天 {cum[7]:.0f}%, 中位 {p.lab_to_it.median():.0f} 天")
    print(f"  (b) HAIC 之前 {vals[:i0].sum():.0f}%, HAIC 当天 {vals[i0]:.0f}%, "
          f"HAIC 之后 {vals[i0 + 1:].sum():.0f}%")
    print(f"      HAIC 前 7 天至 HAIC 当天合计 {within_week:.0f}%")


if __name__ == "__main__":
    main()
