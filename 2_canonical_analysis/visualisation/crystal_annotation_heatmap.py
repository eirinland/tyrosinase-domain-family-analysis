"""
Heatmap of residue frequencies at active-site positions across the canonical
PPO-fold structures, with per-cell pie markers for the characterised entries.
The characterised set (84 entries, 8 activity classes) is loaded directly from
characterized_PPOs.xlsx (single source of truth), so it tracks the current file.
One pie per cell; wedges coloured by activity (one wedge per activity present).
Two figures: crystal_annotation_heatmap (frequency) + crystal_logcount (counts).
"""

import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe

VECTOR_CSV = Path(__file__).parent.parent / 'position_vectors.csv'
CHAR_XLSX = Path('/cluster/work/projects/nn1003k/eirin/bioinf/characterized_PPOs.xlsx')
if not CHAR_XLSX.exists():
    for p in [Path(__file__).parent.parent.parent.parent / 'characterized_PPOs.xlsx',
              Path(__file__).parent.parent.parent / 'characterized_PPOs.xlsx']:
        if p.exists():
            CHAR_XLSX = p
            break

POSITION_LABELS = [
    'Gly46', 'Phe65', 'Trp68', 'Glu195', 'Asn205',
    'Arg209', 'Val218', 'Ala221', 'Phe227', 'His230', 'thioether',
]

DISPLAY_LABELS = [
    'Gly\n46', 'Phe\n65', 'Trp\n68', 'Glu\n195', 'Asn\n205',
    'Arg\n209', 'Val\n218', 'Ala\n221', 'Phe\n227', 'His\n230', 'Cys',
]

COL_COLORS = [
    (1.0, 0.95294118, 0.69019610),
    (0.80392158, 0.70588237, 0.85882354),
    (0.72156864, 0.94901961, 0.90196079),
    (1.0, 0.71764706, 0.75686275),
    (0.73725490, 0.88627451, 0.68235294),
    (0.52941176, 0.78039216, 0.64705882),
    (1.0, 0.82352941, 0.65098039),
    (0.961, 0.678, 0.506),
    (0.80392158, 0.70588237, 0.85882354),
    (0.55686277, 0.77254903, 0.98823529),
    (1.0, 0.95294118, 0.69019610),
]

AA_ORDER = list('GAVLIFWYCMSTPNQDEHKR') + ['-', '~']

# Pastel palette, consistent with the taxonomy and vector-space figures
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


def load_protein_activities(path):
    """Load the characterised set from the Excel: {accession: [activity]}."""
    import pandas as pd
    df = pd.read_excel(path)
    return {str(a).strip(): [str(b).strip()] for a, b in zip(df['Accession'], df['Activity'])}


CMAP = mcolors.LinearSegmentedColormap.from_list(
    'steel_blue',
    ['#EEF2F7', '#C5D5E4', '#8BAFC8', '#5B8FAE', '#3A7198', '#1D5580', '#0A3A5E'])
CMAP.set_bad('#F5F5F5')


def _cell_value(r, pos):
    val = r.get(pos, '?')
    if val == '?':
        return val
    if pos == 'thioether':
        return 'C' if val in ('C', 'C*') else '-'
    if pos == 'Gly46' and r.get('Gly46_ss') == 'c':
        return '~'
    return val


def _draw(rows, vec_by_acc, protein_activities, mode='frequency'):
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
            if mode == 'frequency':
                mat[i, j] = counts.get(aa, 0) / total if total > 0 else 0
            else:
                mat[i, j] = counts.get(aa, 0)

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

    fig, ax = plt.subplots(figsize=(6, 7.5))
    if mode == 'frequency':
        norm = mcolors.PowerNorm(gamma=0.4, vmin=0, vmax=1)
    else:
        norm = mcolors.LogNorm(vmin=1, vmax=max(len(rows), 100))

    mat[mat == 0] = np.nan
    im = ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm, interpolation='nearest')

    ax.set_yticks(range(n_aa))
    ax.set_yticklabels(AA_ORDER, fontsize=10, fontfamily='monospace')
    ax.set_ylabel('')
    ax.xaxis.set_ticks_position('top')
    ax.set_xticks(range(n_pos))
    ax.set_xticklabels(DISPLAY_LABELS, fontsize=11, rotation=0)

    stroke = pe.withStroke(linewidth=1.5, foreground='black')
    for i, label in enumerate(ax.get_xticklabels()):
        label.set_color(COL_COLORS[i])
        label.set_fontweight('bold')
        label.set_path_effects([stroke])

    ax.set_title('')

    pie_frac = 0.8
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_span = xlim[1] - xlim[0]
    y_span = ylim[0] - ylim[1]
    sz = min(pie_frac / x_span, pie_frac / y_span)

    fig.canvas.draw()
    for (j, i), act_counts in crystal_map.items():
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

    if mode == 'frequency':
        fig.text(0.5, 0.08, 'Residue frequency', ha='center', fontsize=12)
    else:
        fig.text(0.5, 0.08, 'Structures per cell', ha='center', fontsize=12)

    cb = fig.colorbar(im, ax=ax, orientation='horizontal', shrink=0.35,
                      pad=0.10, aspect=20, anchor=(0.0, 1.0))
    cb.set_label('')

    handles = [plt.Line2D([0], [0], marker='o', color='w',
                          markerfacecolor=ACTIVITY_COLORS[act],
                          markeredgecolor='black', markeredgewidth=0.3,
                          markersize=6, label=('Hc' if act == 'hemocyanin' else act))
               for act in ACTIVITY_ORDER]
    ax.legend(handles=handles, loc='upper right', bbox_to_anchor=(1.0, -0.14),
              fontsize=8, frameon=False, ncol=3, handletextpad=0.3,
              borderpad=0.5, labelspacing=0.4, columnspacing=1.0)

    return fig


def main():
    with open(VECTOR_CSV) as f:
        rows = list(csv.DictReader(f))
    vec_by_acc = {r['accession'].lower(): r for r in rows}
    protein_activities = load_protein_activities(CHAR_XLSX)

    base = Path(__file__).parent

    fig = _draw(rows, vec_by_acc, protein_activities, mode='frequency')
    fig.savefig(base / 'crystal_annotation_heatmap.png', dpi=1200, bbox_inches='tight')
    fig.savefig(base / 'crystal_annotation_heatmap.pdf', bbox_inches='tight')
    plt.close()
    print('Saved: crystal_annotation_heatmap.pdf / .png')

    fig = _draw(rows, vec_by_acc, protein_activities, mode='logcount')
    fig.savefig(base / 'crystal_logcount.png', dpi=200, bbox_inches='tight')
    fig.savefig(base / 'crystal_logcount.pdf', bbox_inches='tight')
    plt.close()
    print('Saved: crystal_logcount.pdf / .png')


if __name__ == '__main__':
    main()
