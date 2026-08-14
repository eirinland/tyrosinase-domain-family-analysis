import csv, glob, os, statistics
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASEDIR = os.path.dirname(os.path.abspath(__file__))

# --- Load classification (per-site His identity) ---
structs = {}
with open(os.path.join(BASEDIR, 'helix_and_gap_filtered_structures.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        acc = r['accession']
        cua_his = [r.get(f'CuA_His{i}', '---') for i in [1,2,3]]
        cub_his = [r.get(f'CuB_His{i}', '---') for i in [1,2,3]]
        cua_nhis = sum(1 for h in cua_his if h == 'HIS')
        cub_nhis = sum(1 for h in cub_his if h == 'HIS')
        n_his = int(r.get('n_coord_his', 0))
        cls = r.get('classification', '')
        structs[acc] = {
            'cua_nhis': cua_nhis, 'cub_nhis': cub_nhis,
            'n_his': n_his, 'cls': cls,
        }

# --- Load new AllMetal3D NC results (per canonical position) ---
m3d = {}
with open(os.path.join(BASEDIR, 'allmetal3d', 'metal3d_nc_canonical_cu.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r.get('status') != 'ok':
            continue
        acc = r['accession']
        cua_p = float(r['CuA_closest_cu_prob']) if r.get('CuA_closest_cu_prob') else None
        cub_p = float(r['CuB_closest_cu_prob']) if r.get('CuB_closest_cu_prob') else None
        m3d[acc] = {'cua_prob': cua_p, 'cub_prob': cub_p}

# --- Load canonical 200 baseline (old per-AF3-Cu format) ---
can200_probs = []
for fn in sorted(glob.glob(os.path.join(BASEDIR, 'metal3d', 'results', 'metal3d_canonical200_*.tsv'))):
    with open(fn) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r.get('status') != 'ok':
                continue
            prob = float(r['closest_cu_prob']) if r.get('closest_cu_prob') else None
            if prob is not None:
                can200_probs.append(prob)

# --- Load PinMyMetal results ---
pmm = {}
with open(os.path.join(BASEDIR, 'pinmymetal', 'pmm_nc_results.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        pmm[r['accession']] = {
            'n_cu': int(r['n_cu']),
            'n_sites': int(r['n_sites']),
        }

# --- Panel A: AllMetal3D Cu probability by per-site His count ---
site_3his = []
site_2his = []
site_01his = []

for acc, s in structs.items():
    if acc not in m3d:
        continue
    for site, nhis in [('cua', s['cua_nhis']), ('cub', s['cub_nhis'])]:
        prob = m3d[acc][f'{site}_prob']
        if prob is None:
            continue
        if nhis == 3:
            site_3his.append(prob)
        elif nhis == 2:
            site_2his.append(prob)
        else:
            site_01his.append(prob)

# --- Panel B: Paired mononuclear (intact vs degraded site) ---
paired_intact = []
paired_degraded = []

for acc, s in structs.items():
    if s['cls'] != 'mononuclear' or acc not in m3d:
        continue
    intact_p = degraded_p = None
    if s['cua_nhis'] >= 2 and s['cub_nhis'] <= 1:
        intact_p = m3d[acc]['cua_prob']
        degraded_p = m3d[acc]['cub_prob']
    elif s['cub_nhis'] >= 2 and s['cua_nhis'] <= 1:
        intact_p = m3d[acc]['cub_prob']
        degraded_p = m3d[acc]['cua_prob']
    if intact_p is not None and degraded_p is not None:
        paired_intact.append(intact_p)
        paired_degraded.append(degraded_p)

# --- Panel C: Cu prediction rate by His count (3 methods) ---
his_counts = sorted(set(s['n_his'] for s in structs.values()))
af3_rates = []
m3d_rates = []
pmm_rates = []
ns = []

for nh in his_counts:
    accs = [a for a, s in structs.items() if s['n_his'] == nh]
    n = len(accs)
    ns.append(n)
    af3_rates.append(100.0)
    m3d_any = sum(1 for a in accs if a in m3d and (
        m3d[a]['cua_prob'] is not None or m3d[a]['cub_prob'] is not None))
    m3d_rates.append(100 * m3d_any / n if n else 0)
    pmm_any = sum(1 for a in accs if a in pmm and pmm[a]['n_cu'] > 0)
    pmm_rates.append(100 * pmm_any / n if n else 0)

# ---- Figure: 3 panels ----
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(13.5, 4.5),
    gridspec_kw={'width_ratios': [1.3, 1, 1.3]})

C_CAN200 = '#08519C'
C_3HIS = '#2166AC'
C_2HIS = '#67A9CF'
C_01HIS = '#7B2D8E'
C_AF3 = '#D94801'
C_M3D = '#2166AC'
C_PMM = '#41AB5D'

# --- Panel A ---
datasets_a = [can200_probs, site_3his, site_2his, site_01his]
positions_a = [1, 2.3, 3.3, 4.3]
labels_a = ['Canonical\npool', '3 His', '2 His', '0–1 His']
colors_a = [C_CAN200, C_3HIS, C_2HIS, C_01HIS]

parts = ax1.violinplot(datasets_a, positions=positions_a, showmeans=False,
                       showmedians=False, showextrema=False)
for pc, color in zip(parts['bodies'], colors_a):
    pc.set_facecolor(color); pc.set_alpha(0.35)
    pc.set_edgecolor(color); pc.set_linewidth(0.8)

bp = ax1.boxplot(datasets_a, positions=positions_a, widths=0.18, patch_artist=True,
                 showfliers=False, zorder=3,
                 medianprops=dict(color='white', linewidth=1.5),
                 whiskerprops=dict(color='#333333', linewidth=0.8),
                 capprops=dict(color='#333333', linewidth=0.8))
for patch, color in zip(bp['boxes'], colors_a):
    patch.set_facecolor(color); patch.set_alpha(0.85)
    patch.set_edgecolor('#333333'); patch.set_linewidth(0.8)

for d, pos in zip(datasets_a, positions_a):
    ax1.text(pos, -0.08, f'n={len(d)}', ha='center', va='top', fontsize=7.5, color='#555555')

ax1.axvline(1.65, color='#CCCCCC', linewidth=0.7, linestyle=':', zorder=0)
ax1.set_xticks(positions_a)
ax1.set_xticklabels(labels_a, fontsize=8.5)
ax1.set_ylabel('AllMetal3D Cu probability', fontsize=10)
ax1.set_title('Per-site histidine retention', fontsize=10.5, fontweight='bold', pad=8)
ax1.set_ylim(-0.15, 1.05)
ax1.set_xlim(0.3, 5.0)
ax1.axhline(0.5, color='#999999', linewidth=0.7, linestyle='--', zorder=1)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# --- Panel B: paired mononuclear ---
np.random.seed(42)
jitter = np.random.normal(0, 0.04, len(paired_intact))

for g, d, j in zip(paired_intact, paired_degraded, jitter):
    ax2.plot([1+j, 2+j], [g, d], color='#AAAAAA', linewidth=0.25, alpha=0.4, zorder=1)

ax2.scatter(np.ones(len(paired_intact))+jitter, paired_intact,
            s=6, color=C_3HIS, alpha=0.45, zorder=2, edgecolors='none')
ax2.scatter(np.ones(len(paired_degraded))*2+jitter, paired_degraded,
            s=6, color=C_01HIS, alpha=0.45, zorder=2, edgecolors='none')

bp2 = ax2.boxplot([paired_intact, paired_degraded], positions=[1, 2], widths=0.25,
                  patch_artist=True, showfliers=False, zorder=4,
                  medianprops=dict(color='white', linewidth=1.8),
                  whiskerprops=dict(color='#333333', linewidth=0.8),
                  capprops=dict(color='#333333', linewidth=0.8))
for patch, color in zip(bp2['boxes'], [C_3HIS, C_01HIS]):
    patch.set_facecolor(color); patch.set_alpha(0.7)
    patch.set_edgecolor('#333333'); patch.set_linewidth(0.8)

ax2.set_xticks([1, 2])
ax2.set_xticklabels(['Intact\nsite', 'Degraded\nsite'], fontsize=9)
ax2.set_title('Mononuclear: paired sites', fontsize=10.5, fontweight='bold', pad=8)
ax2.set_ylim(-0.15, 1.05)
ax2.axhline(0.5, color='#999999', linewidth=0.7, linestyle='--', zorder=1)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

from scipy.stats import wilcoxon
if len(paired_intact) > 10:
    stat, p = wilcoxon(paired_intact, paired_degraded)
    ptext = f'p = {p:.1e}' if p < 0.001 else f'p = {p:.3f}'
    y_bar = 1.01
    ax2.plot([1, 1, 2, 2], [y_bar-0.02, y_bar, y_bar, y_bar-0.02], color='#333333', linewidth=0.8)
    ax2.text(1.5, y_bar+0.005, ptext, ha='center', va='bottom', fontsize=7.5, style='italic')

ax2.text(1, -0.08, f'n={len(paired_intact)}', ha='center', va='top', fontsize=7.5, color='#555555')
ax2.text(2, -0.08, f'n={len(paired_degraded)}', ha='center', va='top', fontsize=7.5, color='#555555')

# --- Panel C: Cu prediction rate by His count, 3 methods ---
x = np.arange(len(his_counts))
w = 0.25

bars_af3 = ax3.bar(x - w, af3_rates, w, color=C_AF3, alpha=0.85, label='AF3', edgecolor='white', linewidth=0.5)
bars_m3d = ax3.bar(x, m3d_rates, w, color=C_M3D, alpha=0.85, label='AllMetal3D', edgecolor='white', linewidth=0.5)
bars_pmm = ax3.bar(x + w, pmm_rates, w, color=C_PMM, alpha=0.85, label='PinMyMetal', edgecolor='white', linewidth=0.5)

ax3.set_xticks(x)
ax3.set_xticklabels([f'{h}\n(n={n})' for h, n in zip(his_counts, ns)], fontsize=8)
ax3.set_xlabel('Coordinating histidines', fontsize=10)
ax3.set_ylabel('Structures with Cu prediction (%)', fontsize=10)
ax3.set_title('Cu prediction by method', fontsize=10.5, fontweight='bold', pad=8)
ax3.set_ylim(0, 115)
ax3.legend(fontsize=8.5, loc='upper left', framealpha=0.9)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)

for bars in [bars_m3d, bars_pmm]:
    for bar in bars:
        h = bar.get_height()
        if h > 0:
            ax3.text(bar.get_x() + bar.get_width()/2, h + 1.5, f'{h:.0f}',
                     ha='center', va='bottom', fontsize=6.5, color='#555555')

# Panel labels
ax1.text(-0.08, 1.06, 'A', transform=ax1.transAxes, fontsize=14, fontweight='bold', va='top')
ax2.text(-0.08, 1.06, 'B', transform=ax2.transAxes, fontsize=14, fontweight='bold', va='top')
ax3.text(-0.08, 1.06, 'C', transform=ax3.transAxes, fontsize=14, fontweight='bold', va='top')

plt.tight_layout(w_pad=2.5)
outpath = os.path.join(BASEDIR, 'metal3d_confidence_vs_classification.pdf')
plt.savefig(outpath, dpi=300, bbox_inches='tight')
plt.savefig(outpath.replace('.pdf', '.png'), dpi=300, bbox_inches='tight')
print(f'Saved: {outpath}')
print(f'Saved: {outpath.replace(".pdf", ".png")}')

# Print stats
print(f'\nPanel A: can200={len(can200_probs)}, 3His={len(site_3his)}, 2His={len(site_2his)}, 0-1His={len(site_01his)}')
print(f'Panel B: {len(paired_intact)} paired mononuclear structures')
print(f'Panel C: His counts {his_counts}, ns={ns}')
print(f'  AF3:  {[f"{r:.0f}" for r in af3_rates]}')
print(f'  M3D:  {[f"{r:.0f}" for r in m3d_rates]}')
print(f'  PMM:  {[f"{r:.0f}" for r in pmm_rates]}')
