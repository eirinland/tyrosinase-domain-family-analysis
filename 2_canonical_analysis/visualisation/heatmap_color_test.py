#!/usr/bin/env python3
"""Generate 3 heatmap color variants for comparison. Same data, different cmaps."""
import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe

VECTOR_CSV = Path(__file__).parent.parent / 'position_vectors.csv'
CHAR_XLSX = Path(__file__).parent.parent.parent.parent / 'characterized_PPOs.xlsx'
OUT = Path(__file__).parent / 'heatmap_color_comparison.png'

POSITION_LABELS = [
    'Gly46', 'Phe65', 'Trp68', 'Glu195', 'Asn205',
    'Arg209', 'Val218', 'Ala221', 'Phe227', 'His230', 'thioether',
]
DISPLAY_LABELS = [
    'Gly\n46', 'Phe\n65', 'Trp\n68', 'Glu\n195', 'Asn\n205',
    'Arg\n209', 'Val\n218', 'Ala\n221', 'Phe\n227', 'His\n230', 'Cys',
]
AA_ORDER = list('GAVLIFWYCMSTPNQDEHKR') + ['-', '~']

ACTIVITY_COLORS = {
    'TYR': '#58A0DF', 'CaOx': '#9BD0F0', 'AUS': '#B5EAE4',
    'oAPO': '#BEF0BE', 'oMP': '#7BC486', 'DHICA ox': '#FFDBBB',
    'DCT': '#F9BE9E', 'hemocyanin': '#D0B5E4',
}
ACTIVITY_ORDER = ['TYR', 'CaOx', 'AUS', 'oAPO', 'oMP', 'DHICA ox', 'DCT', 'hemocyanin']

# Three candidate colormaps
CMAPS = {
    'A: Steel blue': mcolors.LinearSegmentedColormap.from_list('steel',
        ['#EEF2F7', '#C5D5E4', '#8BAFC8', '#5B8FAE', '#3A7198', '#1D5580', '#0A3A5E']),
    'B: Teal': mcolors.LinearSegmentedColormap.from_list('teal',
        ['#EDF5F4', '#BDD9D5', '#7FBCB5', '#4FA092', '#2E8273', '#1A6557', '#0B4A3D']),
    'C: Blue-purple': mcolors.LinearSegmentedColormap.from_list('bluepurp',
        ['#EEF0F7', '#C5CADE', '#9198C1', '#6B6FAA', '#564E95', '#43337D', '#2E1D5E']),
}

for fname in ["Helvetica Neue", "Helvetica", "Arial"]:
    matches = matplotlib.font_manager.findSystemFonts()
    if any(fname.lower().replace(" ", "") in f.lower().replace(" ", "") for f in matches):
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = [fname]
        break


def _cell_value(r, pos):
    val = r.get(pos, '?')
    if val == '?':
        return val
    if pos == 'thioether':
        return 'C' if val in ('C', 'C*') else '-'
    if pos == 'Gly46' and r.get('Gly46_ss') == 'c':
        return '~'
    return val


def load_protein_activities(path):
    import pandas as pd
    df = pd.read_excel(path)
    return {str(a).strip(): [str(b).strip()] for a, b in zip(df['Accession'], df['Activity'])}


with open(VECTOR_CSV) as f:
    rows = list(csv.DictReader(f))
vec_by_acc = {r['accession'].lower(): r for r in rows}
protein_activities = load_protein_activities(CHAR_XLSX)

n_pos = len(POSITION_LABELS)
n_aa = len(AA_ORDER)
mat = np.zeros((n_aa, n_pos))
for j, pos in enumerate(POSITION_LABELS):
    counts = {}
    for r in rows:
        key = _cell_value(r, pos)
        if key == '?':
            continue
        counts[key] = counts.get(key, 0) + 1
    total = sum(counts.values())
    for i, aa in enumerate(AA_ORDER):
        mat[i, j] = counts.get(aa, 0) / total if total > 0 else 0

crystal_map = {}
for acc, activities in protein_activities.items():
    r = vec_by_acc.get(acc.lower())
    if not r:
        continue
    for j, pos in enumerate(POSITION_LABELS):
        val = _cell_value(r, pos)
        if val == '?' or val not in AA_ORDER:
            continue
        i = AA_ORDER.index(val)
        crystal_map.setdefault((j, i), {})
        for act in activities:
            crystal_map[(j, i)][act] = crystal_map[(j, i)].get(act, 0) + 1

mat_plot = mat.copy()
mat_plot[mat_plot == 0] = np.nan
norm = mcolors.PowerNorm(gamma=0.4, vmin=0, vmax=1)

fig, axes = plt.subplots(1, 3, figsize=(18, 7.5))
fig.subplots_adjust(wspace=0.30)

for idx, (title, cmap) in enumerate(CMAPS.items()):
    ax = axes[idx]
    cmap.set_bad('#F5F5F5')
    im = ax.imshow(mat_plot, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')

    ax.set_yticks(range(n_aa))
    ax.set_yticklabels(AA_ORDER if idx == 0 else [], fontsize=9, fontfamily='monospace')
    ax.xaxis.set_ticks_position('top')
    ax.set_xticks(range(n_pos))
    ax.set_xticklabels(DISPLAY_LABELS, fontsize=9, rotation=0)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=45)

    # Column header colors — muted versions that work with all cmaps
    col_colors = [
        '#C9B860', '#9E84B5', '#5CB8A8', '#D17580', '#7BB562',
        '#3EAD75', '#D4A050', '#CD8A5F', '#9E84B5', '#5590D0', '#C9B860',
    ]
    stroke = pe.withStroke(linewidth=1.5, foreground='black')
    for i, label in enumerate(ax.get_xticklabels()):
        label.set_color(col_colors[i])
        label.set_fontweight('bold')
        label.set_path_effects([stroke])

    # Pie charts for characterised structures
    xlim, ylim = ax.get_xlim(), ax.get_ylim()
    x_span, y_span = xlim[1] - xlim[0], ylim[0] - ylim[1]
    sz = min(0.8 / x_span, 0.8 / y_span)
    fig.canvas.draw()

    for (j, i), act_counts in crystal_map.items():
        acts = [a for a in ACTIVITY_ORDER if a in act_counts]
        if not acts:
            continue
        colors = [ACTIVITY_COLORS[a] for a in acts]
        x_ax = (j - xlim[0]) / x_span
        y_ax = (ylim[0] - i) / y_span
        inax = ax.inset_axes([x_ax - sz/2, y_ax - sz/2, sz, sz],
                             transform=ax.transAxes, zorder=5)
        inax.pie([1]*len(acts), colors=colors, startangle=90,
                 wedgeprops={'edgecolor': 'black', 'linewidth': 0.4})
        inax.set_aspect('equal')
        inax.patch.set_alpha(0)

    cb = fig.colorbar(im, ax=ax, orientation='horizontal', shrink=0.5,
                      pad=0.08, aspect=15)
    cb.set_label('Residue frequency', fontsize=8)

# Activity legend under first panel
handles = [plt.Line2D([0], [0], marker='o', color='w',
                      markerfacecolor=ACTIVITY_COLORS[act],
                      markeredgecolor='black', markeredgewidth=0.3,
                      markersize=6, label=('Hc' if act == 'hemocyanin' else act))
           for act in ACTIVITY_ORDER]
fig.legend(handles=handles, loc='lower center', fontsize=8, frameon=False,
           ncol=8, handletextpad=0.3, columnspacing=1.0,
           bbox_to_anchor=(0.5, -0.02))

fig.savefig(OUT, dpi=200, bbox_inches='tight')
print(f'Saved {OUT}')
