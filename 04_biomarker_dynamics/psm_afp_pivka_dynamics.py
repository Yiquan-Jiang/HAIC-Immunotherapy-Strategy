"""
PSM-matched pairwise biomarker longitudinal dynamics (v4).

Variables: AFP, PIVKA-II, ALBI score, NLR, SII, PIV
Timepoints plotted: Baseline, Pre-HAIC 2-4, Post-6mo, Post-1yr, Post-2yr

Statistical comparisons (two-sample, independent groups):
  - Performed at TWO timepoints only:
      (1) Baseline  — first timepoint
      (2) Last available post-HAIC timepoint per comparison
          (post_2yr → post_1yr → post_6mo, whichever has ≥3 per group)
  - Test: Mann-Whitney U (two-sided); appropriate for skewed, non-normal
    biomarker distributions in independent PSM-matched groups.
  - Multiple-comparison correction: Benjamini-Hochberg FDR across all
    (comparison × variable × 2 timepoints) tests.

Rationale for Mann-Whitney U:
  Biomarkers (AFP, PIVKA, SII, PIV) are heavily right-skewed even after
  log-transformation at the group level. PSM produces matched pairs by
  group, not within-patient pairs, so paired tests (Wilcoxon signed-rank)
  are not appropriate. Mann-Whitney U is the standard non-parametric
  two-sample test for such data.

Output: figures/psm_biomarker_dynamics/remake/
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import glob
import os
import warnings
warnings.filterwarnings('ignore')

# ── rcParams ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':           7,
    'axes.labelsize':      7,
    'axes.titlesize':      8,
    'xtick.labelsize':     6,
    'ytick.labelsize':     6,
    'legend.fontsize':     6,
    'axes.linewidth':      0.5,
    'axes.spines.top':     False,
    'axes.spines.right':   False,
    'xtick.major.width':   0.5,
    'ytick.major.width':   0.5,
    'xtick.major.size':    2.5,
    'ytick.major.size':    2.5,
    'lines.linewidth':     1.2,
    'lines.markersize':    4,
    'legend.frameon':      False,
    'legend.borderpad':    0.2,
    'legend.handlelength': 1.5,
    'figure.dpi':         150,
    'savefig.dpi':        300,
    'savefig.bbox':       'tight',
    'savefig.pad_inches':  0.05,
    'pdf.fonttype':        42,
    'ps.fonttype':         42,
    'axes.grid':           False,
})

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE_DIR, 'data')
PSM_DIR  = os.path.join(BASE_DIR, 'results', 'psm_balance_tables_complete')
FIG_DIR  = os.path.join(BASE_DIR, 'figures', 'psm_biomarker_dynamics', 'remake')
os.makedirs(FIG_DIR, exist_ok=True)

# ── Constants ─────────────────────────────────────────────────────────────────
TIMEPOINTS = [
    'baseline', 'pre_haic_2', 'pre_haic_3', 'pre_haic_4',
    'post_6mo', 'post_1yr', 'post_2yr',
]
TP_LABELS = ['BL', 'Pre-2', 'Pre-3', 'Pre-4', '6mo', '1yr', '2yr']
N_TP      = len(TIMEPOINTS)
ALPHA     = 0.05

# Post-HAIC landmark windows (days from baseline)
POST_LANDMARKS = {
    'post_6mo': {'target_days': 183, 'lo':  90, 'hi': 269},
    'post_1yr': {'target_days': 365, 'lo': 270, 'hi': 450},
    'post_2yr': {'target_days': 730, 'lo': 630, 'hi': 810},
}
# Ordered from latest to earliest for "last available" fallback
POST_PRIORITY = ['post_2yr', 'post_1yr', 'post_6mo']

COLOR_G1 = '#4DBBD5'
COLOR_G2 = '#E64B35'

GROUP_PRETTY = {
    'HAIC_alone':          'HAIC alone',
    'HAIC+I_concurrent':   'HAIC+I (conc.)',
    'HAIC_then_I':         'HAIC\u2192I (deferred)',
    'HAIC+T_concurrent':   'HAIC+T (conc.)',
    'HAIC_then_T':         'HAIC\u2192T (deferred)',
    'HAIC+I+T_concurrent': 'HAIC+I+T (conc.)',
    'HAIC_then_I+T':       'HAIC\u2192I+T (deferred)',
}

STATS_BOX = dict(boxstyle='round,pad=0.3', facecolor='white',
                 edgecolor='#CCCCCC', alpha=0.95, linewidth=0.4)

# (col_name, display_label, scale_type)
VARIABLES = [
    ('afp',        'AFP (ng/mL)',        'log'),
    ('pivka',      'PIVKA-II (mAU/mL)', 'log'),
    ('albi_score', 'ALBI Score',         'linear'),
    ('nlr',        'NLR',               'linear'),
    ('sii',        'SII',               'log'),
    ('piv',        'PIV',               'log'),
]
N_VAR = len(VARIABLES)


# ── Load and preprocess longitudinal data ────────────────────────────────────
print("Loading longitudinal data...")
df_raw = pd.read_csv(os.path.join(DATA_DIR, 'HAIC_NO_TACE_4_TIDY_longitudinal.csv'))
df_raw['haic_date'] = pd.to_datetime(df_raw['haic_date'])

df_raw['sii'] = df_raw['plt'] * df_raw['neut'] / df_raw['lymph']
df_raw['piv'] = df_raw['plt'] * df_raw['mono'] * df_raw['neut'] / df_raw['lymph']
df_raw.loc[df_raw['lymph'] <= 0, ['sii', 'piv']] = np.nan

bl_dates = (df_raw[df_raw['timepoint_type'] == 'baseline']
            [['patient_id', 'haic_date']]
            .rename(columns={'haic_date': 'baseline_date'}))

pre_tps = ['baseline', 'pre_haic_2', 'pre_haic_3', 'pre_haic_4']
df_pre = df_raw[df_raw['timepoint_type'].isin(pre_tps)].copy()

df_post_raw = df_raw[df_raw['timepoint_type'] == 'post_haic'].copy()
df_post_raw = df_post_raw.merge(bl_dates, on='patient_id', how='left')
df_post_raw['days_from_bl'] = (
    df_post_raw['haic_date'] - df_post_raw['baseline_date']).dt.days

post_frames = []
for tp_name, spec in POST_LANDMARKS.items():
    window = df_post_raw[
        (df_post_raw['days_from_bl'] >= spec['lo']) &
        (df_post_raw['days_from_bl'] <= spec['hi'])
    ].copy()
    window['dist_to_target'] = (window['days_from_bl'] - spec['target_days']).abs()
    closest = (window.sort_values('dist_to_target')
               .groupby('patient_id').first().reset_index())
    closest['timepoint_type'] = tp_name
    post_frames.append(closest)

df_post = pd.concat(post_frames, ignore_index=True)
drop_cols = [c for c in ['baseline_date', 'days_from_bl', 'dist_to_target']
             if c in df_post.columns]
df_post.drop(columns=drop_cols, inplace=True)

df_long = pd.concat([df_pre, df_post], ignore_index=True)
df_long['timepoint_type'] = pd.Categorical(
    df_long['timepoint_type'], categories=TIMEPOINTS, ordered=True)

print(f"  Pre-HAIC rows : {len(df_pre)}")
for tp in POST_PRIORITY:
    n = len(df_post[df_post['timepoint_type'] == tp])
    print(f"  {tp} patients : {n}")
print(f"  Total rows    : {len(df_long)}")


# ── Discover PSM matched-ID files ─────────────────────────────────────────────
matched_files = sorted(glob.glob(os.path.join(PSM_DIR, 'matched_ids_*.csv')))
print(f"\n  Found {len(matched_files)} PSM comparison files.\n")


# ── Statistics helpers ─────────────────────────────────────────────────────────
def compute_group_stats(data, var):
    vals = data[var].dropna()
    n = len(vals)
    if n < 3:
        return np.nan, np.nan, np.nan, n
    return np.median(vals), np.percentile(vals, 25), np.percentile(vals, 75), n


def mw_test(data, var, group_col, g1_name, g2_name):
    """Mann-Whitney U, two-sided. Returns (U, p)."""
    a = data.loc[data[group_col] == g1_name, var].dropna()
    b = data.loc[data[group_col] == g2_name, var].dropna()
    if len(a) < 3 or len(b) < 3:
        return np.nan, np.nan
    stat, p = mannwhitneyu(a, b, alternative='two-sided')
    return stat, p


def fmt_p(p, fdr_p=None):
    """Format raw p and optionally FDR-corrected p."""
    if np.isnan(p):
        return 'NA'
    raw = 'P < 0.001' if p < 0.001 else f'P = {p:.3f}'
    if fdr_p is not None and not np.isnan(fdr_p):
        fdr = 'q < 0.001' if fdr_p < 0.001 else f'q = {fdr_p:.3f}'
        return f'{raw}\n{fdr}'
    return raw


def last_available_tp(df_comp, g1_name, g2_name, var, min_n=3):
    """Return the latest post-HAIC timepoint with ≥min_n per group."""
    for tp in POST_PRIORITY:
        tp_data = df_comp[df_comp['timepoint_type'] == tp]
        n1 = tp_data.loc[tp_data['group_label'] == g1_name, var].dropna().shape[0]
        n2 = tp_data.loc[tp_data['group_label'] == g2_name, var].dropna().shape[0]
        if n1 >= min_n and n2 >= min_n:
            return tp
    return None


# ── Pass 1: collect all raw p-values for FDR correction ──────────────────────
print("Pass 1: collecting raw p-values for FDR correction...")
raw_records = []   # list of dicts: comparison, var, tp_role, g1, g2, U, raw_p

for fpath in matched_files:
    fname    = os.path.basename(fpath)
    comp_key = fname.replace('matched_ids_', '').replace('.csv', '')
    df_ids   = pd.read_csv(fpath)
    g1_name  = df_ids['group1'].iloc[0]
    g2_name  = df_ids['group2'].iloc[0]

    df_comp = df_long[df_long['patient_id'].isin(df_ids['patient_id'])].copy()
    df_comp = df_comp.merge(
        df_ids[['patient_id', 'group_label']].drop_duplicates(),
        on='patient_id', how='left')

    for var, _, _ in VARIABLES:
        last_tp = last_available_tp(df_comp, g1_name, g2_name, var)

        for tp_role, tp in [('baseline', 'baseline'), ('last', last_tp)]:
            if tp is None:
                raw_records.append(dict(
                    comparison=comp_key, variable=var, tp_role=tp_role,
                    tp=tp, group1=g1_name, group2=g2_name,
                    U=np.nan, raw_p=np.nan))
                continue
            tp_data = df_comp[df_comp['timepoint_type'] == tp]
            U, p = mw_test(tp_data, var, 'group_label', g1_name, g2_name)
            raw_records.append(dict(
                comparison=comp_key, variable=var, tp_role=tp_role,
                tp=tp, group1=g1_name, group2=g2_name,
                U=U, raw_p=p))

pval_df = pd.DataFrame(raw_records)

# BH-FDR across all valid tests
valid_mask = ~pval_df['raw_p'].isna()
if valid_mask.sum() > 0:
    _, fdr_ps, _, _ = multipletests(
        pval_df.loc[valid_mask, 'raw_p'].values, alpha=0.05, method='fdr_bh')
    pval_df.loc[valid_mask, 'fdr_p'] = fdr_ps
else:
    pval_df['fdr_p'] = np.nan

print(f"  Total tests: {valid_mask.sum()} valid / {len(pval_df)} total")
sig_raw = (pval_df['raw_p'] < 0.05).sum()
sig_fdr = (pval_df['fdr_p'] < 0.05).sum()
print(f"  Significant: {sig_raw} raw P<0.05, {sig_fdr} FDR q<0.05\n")


def get_pvals(comp_key, var, tp_role):
    row = pval_df[
        (pval_df['comparison'] == comp_key) &
        (pval_df['variable']   == var) &
        (pval_df['tp_role']    == tp_role)]
    if len(row) == 0:
        return np.nan, np.nan, None
    r = row.iloc[0]
    return r['raw_p'], r['fdr_p'], r['tp']


# ── Pass 2: plot ──────────────────────────────────────────────────────────────
print("Pass 2: generating figures...")
all_stats = []

for fpath in matched_files:
    fname    = os.path.basename(fpath)
    comp_key = fname.replace('matched_ids_', '').replace('.csv', '')
    comp_num = comp_key[:2]
    print(f"  [{comp_num}] {comp_key}")

    df_ids       = pd.read_csv(fpath)
    g1_name      = df_ids['group1'].iloc[0]
    g2_name      = df_ids['group2'].iloc[0]
    n_per_group  = len(df_ids) // 2
    g1_pretty    = GROUP_PRETTY.get(g1_name, g1_name)
    g2_pretty    = GROUP_PRETTY.get(g2_name, g2_name)

    df_comp = df_long[df_long['patient_id'].isin(df_ids['patient_id'])].copy()
    df_comp = df_comp.merge(
        df_ids[['patient_id', 'group_label']].drop_duplicates(),
        on='patient_id', how='left')

    # ── Multi-panel figure (3 rows × 2 cols) ──────────────────────────────
    n_rows, n_cols = 3, 2
    fig, axes_flat = plt.subplots(n_rows, n_cols, figsize=(6, 6))
    axes_all = axes_flat.flatten()

    for vi, (var, var_label, scale) in enumerate(VARIABLES):
        ax = axes_all[vi]

        # Collect per-timepoint stats
        xs = np.arange(N_TP)
        jitter = 0.08
        meds1, lo1, hi1, ns1 = [], [], [], []
        meds2, lo2, hi2, ns2 = [], [], [], []

        for tp in TIMEPOINTS:
            tp_data = df_comp[df_comp['timepoint_type'] == tp]
            g1_data = tp_data[tp_data['group_label'] == g1_name]
            g2_data = tp_data[tp_data['group_label'] == g2_name]
            med1, q25_1, q75_1, n1 = compute_group_stats(g1_data, var)
            med2, q25_2, q75_2, n2 = compute_group_stats(g2_data, var)
            meds1.append(med1); lo1.append(q25_1); hi1.append(q75_1); ns1.append(n1)
            meds2.append(med2); lo2.append(q25_2); hi2.append(q75_2); ns2.append(n2)

            all_stats.append({
                'comparison': comp_key, 'group1': g1_name, 'group2': g2_name,
                'variable': var, 'timepoint': tp,
                'n1': n1, 'n2': n2,
                'median_g1': med1,
                'iqr_g1': f'{q25_1:.2f}-{q75_1:.2f}' if not np.isnan(q25_1) else 'NA',
                'median_g2': med2,
                'iqr_g2': f'{q25_2:.2f}-{q75_2:.2f}' if not np.isnan(q25_2) else 'NA',
            })

        meds1 = np.array(meds1, dtype=float); meds2 = np.array(meds2, dtype=float)
        lo1   = np.array(lo1,   dtype=float); lo2   = np.array(lo2,   dtype=float)
        hi1   = np.array(hi1,   dtype=float); hi2   = np.array(hi2,   dtype=float)

        err1_lo = np.clip(meds1 - lo1, 0, None); err1_hi = np.clip(hi1 - meds1, 0, None)
        err2_lo = np.clip(meds2 - lo2, 0, None); err2_hi = np.clip(hi2 - meds2, 0, None)
        valid1  = ~np.isnan(meds1);               valid2  = ~np.isnan(meds2)

        show_legend = (vi == 0)

        if valid1.sum() >= 2:
            ax.errorbar(xs[valid1] - jitter, meds1[valid1],
                        yerr=[err1_lo[valid1], err1_hi[valid1]],
                        fmt='none', ecolor=COLOR_G1, elinewidth=0.6,
                        capsize=2.0, capthick=0.6, alpha=0.5, zorder=2)
            ax.plot(xs[valid1] - jitter, meds1[valid1],
                    color=COLOR_G1, marker='o', markersize=4, linewidth=1.2,
                    label=g1_pretty if show_legend else None,
                    markeredgecolor='white', markeredgewidth=0.3, zorder=3)

        if valid2.sum() >= 2:
            ax.errorbar(xs[valid2] + jitter, meds2[valid2],
                        yerr=[err2_lo[valid2], err2_hi[valid2]],
                        fmt='none', ecolor=COLOR_G2, elinewidth=0.6,
                        capsize=2.0, capthick=0.6, alpha=0.5, zorder=2)
            ax.plot(xs[valid2] + jitter, meds2[valid2],
                    color=COLOR_G2, marker='s', markersize=4, linewidth=1.2,
                    linestyle='--',
                    label=g2_pretty if show_legend else None,
                    markeredgecolor='white', markeredgewidth=0.3, zorder=3)

        # Y-axis scale
        if scale == 'log':
            all_hi_v = np.concatenate([hi1[~np.isnan(hi1)], hi2[~np.isnan(hi2)]])
            all_lo_v = np.concatenate([lo1[~np.isnan(lo1)], lo2[~np.isnan(lo2)]])
            if len(all_hi_v) > 0 and len(all_lo_v) > 0:
                y_max = np.max(all_hi_v) * 2.0
                y_min = (np.min(all_lo_v[all_lo_v > 0]) * 0.4
                         if np.any(all_lo_v > 0) else 1)
                ax.set_ylim(y_min, y_max)
            ax.set_yscale('log')
            ax.yaxis.set_major_formatter(mticker.FuncFormatter(
                lambda x, _: f'{x:g}' if x >= 1 else f'{x:.2g}'))

        # ── Annotate p-values at baseline and last timepoint only ─────────
        raw_bl, fdr_bl, _   = get_pvals(comp_key, var, 'baseline')
        raw_last, fdr_last, last_tp = get_pvals(comp_key, var, 'last')

        tp_idx_bl   = TIMEPOINTS.index('baseline')
        tp_idx_last = TIMEPOINTS.index(last_tp) if last_tp else None

        for tp_idx, raw_p, fdr_p in [
            (tp_idx_bl,   raw_bl,   fdr_bl),
            (tp_idx_last, raw_last, fdr_last),
        ]:
            if tp_idx is None or np.isnan(raw_p):
                continue
            # Star marker at top of column
            if fdr_p < 0.001:
                star = '***'
            elif fdr_p < 0.01:
                star = '**'
            elif fdr_p < 0.05:
                star = '*'
            elif raw_p < 0.05:
                star = '†'   # raw sig but FDR not
            else:
                star = ''

            if star:
                ax.text(xs[tp_idx], 0.96, star,
                        ha='center', va='top', fontsize=9,
                        fontweight='bold',
                        color='#333333' if fdr_p < 0.05 else '#999999',
                        transform=ax.get_xaxis_transform())

            # Small p-value box below x-axis (show raw p)
            p_txt = ('P<.001' if raw_p < 0.001 else f'P={raw_p:.3f}')
            ax.text(xs[tp_idx], -0.14, p_txt,
                    ha='center', va='top', fontsize=5.5,
                    color='#333333' if fdr_p < 0.05 else '#888888',
                    transform=ax.get_xaxis_transform(),
                    style='italic')

        ax.set_xticks(xs)
        ax.set_xticklabels(TP_LABELS, fontsize=7)
        ax.set_ylabel(var_label, fontsize=9)
        ax.set_xlim(-0.5, N_TP - 0.5)

        if show_legend:
            ax.legend(loc='upper right', fontsize=7, handlelength=1.8,
                      frameon=True, framealpha=0.95,
                      edgecolor='#CCCCCC', fancybox=False)

        panel_letter = chr(65 + vi)
        ax.text(-0.08, 1.05, panel_letter, transform=ax.transAxes,
                fontsize=11, fontweight='bold', va='top')

    # Hide unused subplot(s)
    for idx in range(N_VAR, n_rows * n_cols):
        axes_all[idx].set_visible(False)

    fig.suptitle(
        f'{g1_pretty}  vs  {g2_pretty}   (PSM matched, n = {n_per_group} per group)',
        fontsize=12, fontweight='bold', y=0.995)

    # Legend for significance notation
    fig.text(0.5, -0.01,
             '*/†: BH-FDR q<0.05/raw P<0.05 only; **: q<0.01; ***: q<0.001  '
             '(Mann-Whitney U, annotated at baseline & last available timepoint)',
             ha='center', fontsize=5.5, color='#555555', style='italic')

    plt.tight_layout(rect=[0, 0.01, 1, 0.98], h_pad=1.4, w_pad=1.0)

    out_stem = os.path.join(FIG_DIR, f'comp_{comp_num}_{g1_name}_vs_{g2_name}')
    fig.savefig(out_stem + '.pdf', bbox_inches='tight', pad_inches=0.05)
    fig.savefig(out_stem + '.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(f"    Saved: {out_stem}.pdf/.png")


# ── Export stats CSV ──────────────────────────────────────────────────────────
stats_df = pd.DataFrame(all_stats)
# Merge FDR p-values for baseline and last timepoint
for role in ['baseline', 'last']:
    sub = pval_df[pval_df['tp_role'] == role][
        ['comparison', 'variable', 'tp', 'U', 'raw_p', 'fdr_p']
    ].rename(columns={
        'tp': f'tp_{role}', 'U': f'U_{role}',
        'raw_p': f'raw_p_{role}', 'fdr_p': f'fdr_p_{role}'})
    stats_df = stats_df.merge(
        sub, on=['comparison', 'variable'], how='left')

csv_path = os.path.join(FIG_DIR, 'psm_biomarker_stats_summary.csv')
stats_df.to_csv(csv_path, index=False)
print(f"\nSaved stats CSV: {csv_path}")

pval_path = os.path.join(FIG_DIR, 'mw_pvalues_fdr.csv')
pval_df.to_csv(pval_path, index=False)
print(f"Saved p-value table: {pval_path}")


# ── P-value heatmap (2 columns: baseline & last timepoint) ───────────────────
print("\nGenerating P-value heatmaps (baseline vs last timepoint)...")

comparisons_all = pval_df['comparison'].unique()
n_comp = len(comparisons_all)

for var, var_label, _ in VARIABLES:
    var_pv = pval_df[pval_df['variable'] == var].copy()
    comp_labels = []
    pval_matrix = np.full((n_comp, 2), np.nan)   # col0=baseline, col1=last
    fdr_matrix  = np.full((n_comp, 2), np.nan)

    for ci, comp in enumerate(comparisons_all):
        sub = var_pv[var_pv['comparison'] == comp]
        g1  = GROUP_PRETTY.get(sub['group1'].iloc[0], sub['group1'].iloc[0])
        g2  = GROUP_PRETTY.get(sub['group2'].iloc[0], sub['group2'].iloc[0])
        comp_labels.append(f'{comp[:2]}: {g1} vs {g2}')

        for col_i, role in enumerate(['baseline', 'last']):
            row = sub[sub['tp_role'] == role]
            if len(row) > 0:
                pval_matrix[ci, col_i] = row['raw_p'].values[0]
                fdr_matrix[ci,  col_i] = row['fdr_p'].values[0]

    log_pvals = -np.log10(np.clip(pval_matrix, 1e-10, 1.0))

    fig, ax = plt.subplots(figsize=(4.5, max(4.0, n_comp * 0.38 + 1.2)))
    im = ax.imshow(log_pvals, aspect='auto', cmap='RdYlGn_r',
                   vmin=0, vmax=4, interpolation='nearest')

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Baseline', 'Last timepoint'], fontsize=8)
    ax.set_yticks(range(n_comp))
    ax.set_yticklabels(comp_labels, fontsize=7)

    for ci in range(n_comp):
        for col_i in range(2):
            p   = pval_matrix[ci, col_i]
            q   = fdr_matrix[ci,  col_i]
            if np.isnan(p):
                txt, color, fw = 'NA', '#999999', 'normal'
            else:
                p_str = '<.001' if p < 0.001 else f'.{int(round(p * 1000)):03d}'
                txt   = f'P={p_str}'
                color = 'white' if log_pvals[ci, col_i] > 2.5 else '#333333'
                fw    = 'bold' if (not np.isnan(q) and q < 0.05) else 'normal'
            ax.text(col_i, ci, txt, ha='center', va='center',
                    fontsize=5.5, color=color, fontweight=fw,
                    multialignment='center')

    for ci in range(n_comp):
        ax.axhline(ci + 0.5, color='white', linewidth=0.3)
    ax.axhline(-0.5, color='white', linewidth=0.5)

    cbar = fig.colorbar(im, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label('$-\\log_{10}(P_{raw})$', fontsize=8)
    cbar.ax.axhline(-np.log10(0.05), color='black', linewidth=1, linestyle='--')
    cbar.ax.text(1.5, -np.log10(0.05), 'P=0.05', fontsize=6, va='center')

    ax.set_title(
        f'Mann-Whitney U: {var_label}\n'
        f'({n_comp} PSM pairs; BH-FDR corrected; bold = q<0.05)',
        fontsize=9, fontweight='bold', pad=8)

    plt.tight_layout()
    out_stem = os.path.join(FIG_DIR, f'pvalue_heatmap_{var}')
    fig.savefig(out_stem + '.pdf', bbox_inches='tight', pad_inches=0.05)
    fig.savefig(out_stem + '.png', dpi=300, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    print(f"  Saved: {out_stem}.pdf/.png")


print("\n" + "=" * 60)
print("METHODOLOGY NOTES:")
print(f"  Variables : {N_VAR} (AFP, PIVKA, ALBI, NLR, SII, PIV)")
print(f"  Plotted   : all {N_TP} timepoints (median ± IQR)")
print(f"  Tested    : baseline  AND  last available post-HAIC timepoint")
print(f"              (post_2yr → post_1yr → post_6mo, ≥3 per group)")
print(f"  Test      : Mann-Whitney U (two-sided, independent samples)")
print(f"  Correction: Benjamini-Hochberg FDR across all valid tests")
print(f"  SII = PLT × Neut / Lymph")
print(f"  PIV = PLT × Mono × Neut / Lymph")
print("=" * 60)
print(f"\nAll done! Output: {FIG_DIR}")
