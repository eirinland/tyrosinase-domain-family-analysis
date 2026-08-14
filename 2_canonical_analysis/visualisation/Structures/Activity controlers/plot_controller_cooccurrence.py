"""Co-occurrence heatmap of Asn205 x Arg209 activity controller substitutions,
with pie markers for characterized PPOs (matching crystal_annotation_heatmap.py)."""

import csv
from pathlib import Path
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

VECTOR_CSV = Path(__file__).resolve().parent.parent.parent.parent / 'position_vectors.csv'
OUTDIR = Path(__file__).resolve().parent

CHAR_XLSX = Path('/cluster/work/projects/nn1003k/eirin/bioinf/characterized_PPOs.xlsx')
if not CHAR_XLSX.exists():
    for p in [Path(__file__).resolve().parent.parent.parent.parent.parent / 'characterized_PPOs.xlsx',
              Path(__file__).resolve().parent.parent.parent.parent / 'characterized_PPOs.xlsx']:
        if p.exists():
            CHAR_XLSX = p
            break

ACTIVITY_COLORS = {
    'TYR':        '#58A0DF',
    'CaOx':       '#9BD0F0',
    'AUS':        '#B5EAE4',
    'oAPO':       '#BEF0BE',
    'oMP':        '#7BC486',
    'DHICA ox':   '#FFDBBB',
    'DCT':        '#F9BE9E',
    'hemocyanin': '#D0B5E4',
}
ACTIVITY_ORDER = ['TYR', 'CaOx', 'AUS', 'oAPO', 'oMP', 'DHICA ox', 'DCT', 'hemocyanin']

with open(VECTOR_CSV) as f:
    rows = list(csv.DictReader(f))

pairs = [(r['Asn205'], r['Arg209']) for r in rows if r['Asn205'] not in ('?', '~') and r['Arg209'] not in ('?', '~')]
counts = Counter(pairs)

n205_totals = Counter(p[0] for p in pairs)
r209_totals = Counter(p[1] for p in pairs)
n205_order = [r for r, _ in n205_totals.most_common()]
r209_order = [r for r, _ in r209_totals.most_common()]

mat = np.zeros((len(n205_order), len(r209_order)))
for i, n in enumerate(n205_order):
    for j, r in enumerate(r209_order):
        mat[i, j] = counts.get((n, r), 0)

# Load characterized PPOs
import pandas as pd
char_df = pd.read_excel(CHAR_XLSX)
vec_by_acc = {r['accession'].lower(): r for r in rows}

crystal_map = {}
for _, row in char_df.iterrows():
    acc = str(row['Accession']).strip().lower()
    act = str(row['Activity']).strip()
    v = vec_by_acc.get(acc)
    if v is None or v['Asn205'] in ('?', '~') or v['Arg209'] in ('?', '~'):
        continue
    n205, r209 = v['Asn205'], v['Arg209']
    if n205 in n205_order and r209 in r209_order:
        i = n205_order.index(n205)
        j = r209_order.index(r209)
        crystal_map.setdefault((i, j), {})
        crystal_map[(i, j)][act] = crystal_map[(i, j)].get(act, 0) + 1

fig, ax = plt.subplots(figsize=(7, 8))

cmap = mcolors.LinearSegmentedColormap.from_list(
    'steel_blue',
    ['#EEF2F7', '#C5D5E4', '#8BAFC8', '#5B8FAE', '#3A7198', '#1D5580', '#0A3A5E'])
cmap.set_under('white')
norm = mcolors.LogNorm(vmin=1, vmax=mat.max())

im = ax.imshow(mat, cmap=cmap, norm=norm, aspect='auto')


ax.set_xticks(range(len(r209_order)))
ax.set_xticklabels(r209_order, fontsize=12, fontfamily='monospace')
ax.set_yticks(range(len(n205_order)))
ax.set_yticklabels(n205_order, fontsize=12, fontfamily='monospace')

ax.set_xlabel('Arg209 substitution', fontsize=13)
ax.set_ylabel('Asn205 substitution', fontsize=13)

n205_ref = n205_order.index('N') if 'N' in n205_order else None
r209_ref = r209_order.index('R') if 'R' in r209_order else None
if n205_ref is not None:
    ax.get_yticklabels()[n205_ref].set_fontweight('bold')
if r209_ref is not None:
    ax.get_xticklabels()[r209_ref].set_fontweight('bold')

# Pie markers for characterized PPOs
pie_frac = 0.7
xlim = ax.get_xlim()
ylim = ax.get_ylim()
x_span = xlim[1] - xlim[0]
y_span = ylim[0] - ylim[1]
sz = min(pie_frac / x_span, pie_frac / y_span)

fig.canvas.draw()
for (i, j), act_counts in crystal_map.items():
    acts = [a for a in ACTIVITY_ORDER if a in act_counts]
    if not acts:
        continue
    colors = [ACTIVITY_COLORS[a] for a in acts]
    x_ax = (j - xlim[0]) / x_span
    y_ax = (ylim[0] - i) / y_span
    inax = ax.inset_axes([x_ax - sz / 2, y_ax - sz / 2, sz, sz],
                         transform=ax.transAxes, zorder=5)
    inax.pie([1] * len(acts), colors=colors, startangle=90,
             wedgeprops={'edgecolor': 'black', 'linewidth': 0.4})
    inax.set_aspect('equal')
    inax.patch.set_alpha(0)

handles = [plt.Line2D([0], [0], marker='o', color='w',
                       markerfacecolor=ACTIVITY_COLORS[act],
                       markeredgecolor='black', markeredgewidth=0.5,
                       markersize=9, label=('Hc' if act == 'hemocyanin' else act))
           for act in ACTIVITY_ORDER]
leg = ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(0.0, -0.16),
          fontsize=9, frameon=False, ncol=4, handletextpad=0.3,
          borderpad=0.3, labelspacing=0.4, columnspacing=0.8)

cbar = fig.colorbar(im, ax=ax, orientation='horizontal', shrink=0.3,
                    pad=0.12, aspect=15, anchor=(1.0, 1.0))
cbar.set_label('Count', fontsize=9)
cbar.ax.tick_params(labelsize=8)

plt.tight_layout()
for ext in ('pdf', 'png'):
    fig.savefig(OUTDIR / f'controller_cooccurrence.{ext}', dpi=300, bbox_inches='tight')
print(f'Saved to {OUTDIR}/controller_cooccurrence.pdf/.png')
print(f'Total pairs: {len(pairs)}, unique combos: {len(counts)}')
print(f'Characterized combos with pies: {len(crystal_map)}')
print(f'Top 10 combos:')
for (n, r), c in counts.most_common(10):
    print(f'  N205={n} + R209={r}: {c}')
