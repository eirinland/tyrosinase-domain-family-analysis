"""
Taxonomy-stratified frequency heatmaps: one panel per kingdom showing
residue frequencies at the 11 variable active-site positions.
"""

import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D

VECTORS = Path(__file__).parent.parent / 'position_vectors.csv'
TAXONOMY = Path(__file__).parent / 'taxonomy_lookup.csv'
CHAR_XLSX = Path(__file__).parent.parent / 'characterized_PPOs.xlsx'

POSITIONS = [
    'Gly46', 'Phe65', 'Trp68', 'Glu195', 'Asn205',
    'Arg209', 'Val218', 'Ala221', 'Phe227', 'His230', 'thioether',
]
DISPLAY = [
    'Gly46', 'Phe65', 'Trp68', 'Glu195', 'Asn205',
    'Arg209', 'Val218', 'Ala221', 'Phe227', 'His230', 'Cys',
]
AA_ORDER = list('GAVLIFWYCMSTPNQDEHKR') + ['-', '~']
KINGDOMS = ['Fungi', 'Animals', 'Plants', 'Bacteria', 'Oomycota']
KINGDOM_LABELS = {'Fungi': 'Fungi', 'Animals': 'Metazoa', 'Plants': 'Viridiplantae',
                  'Bacteria': 'Bacteria', 'Oomycota': 'Oomycota'}

ACTIVITY_COLORS = {
    'TYR': '#58A0DF', 'CaOx': '#9BD0F0', 'AUS': '#B5EAE4', 'oAPO': '#BEF0BE',
    'oMP': '#7BC486', 'DHICA ox': '#FFDBBB', 'DCT': '#F9BE9E', 'hemocyanin': '#D0B5E4',
}
ACTIVITY_ORDER = ['TYR', 'CaOx', 'AUS', 'oAPO', 'oMP', 'DHICA ox', 'DCT', 'hemocyanin']

# Warm occupancy palette (matches the crystal heatmap): wheat -> deep brick,
# with a faint tint for empty cells so nothing reads as stark white.
CMAP = mcolors.LinearSegmentedColormap.from_list(
    'warm_soft',
    ['#fbe3c2', '#f9d2a3', '#f4b67f', '#ec9559', '#dd6b49', '#c0563b', '#94381f', '#6e2a1a'])
CMAP.set_bad('#fdf3e6')


def load_protein_activities(path):
    import pandas as pd
    df = pd.read_excel(path)
    return {str(a).strip(): [str(b).strip()] for a, b in zip(df['Accession'], df['Activity'])}


def draw_activity_pies(ax, cmap):
    """Overlay per-cell activity pies (same method as crystal_annotation_heatmap.py)."""
    xlim = ax.get_xlim(); ylim = ax.get_ylim()
    xs = xlim[1] - xlim[0]; ys = ylim[0] - ylim[1]  # y axis is inverted
    sz = min(1.35 / xs, 1.35 / ys)
    for (j, i), ac in cmap.items():
        acts = [a for a in ACTIVITY_ORDER if a in ac]
        if not acts:
            continue
        cols = [ACTIVITY_COLORS[a] for a in acts]
        x_ax = (j - xlim[0]) / xs
        y_ax = (ylim[0] - i) / ys
        inax = ax.inset_axes([x_ax - sz / 2, y_ax - sz / 2, sz, sz],
                             transform=ax.transAxes, zorder=10)
        inax.pie([1] * len(acts), colors=cols, startangle=90,
                 wedgeprops={'edgecolor': 'black', 'linewidth': 0.3})
        inax.set_aspect('equal')
        inax.patch.set_alpha(0)


def _collect_data():
    tax = {}
    with open(TAXONOMY) as f:
        for row in csv.DictReader(f):
            tax[row['accession']] = row['kingdom']

    protein_activities = load_protein_activities(CHAR_XLSX)
    vec_by_acc = {}
    rows_by_kingdom = {k: [] for k in KINGDOMS}
    with open(VECTORS) as f:
        for row in csv.DictReader(f):
            vec_by_acc[row['accession']] = row
            k = tax.get(row['accession'], '?')
            if k in rows_by_kingdom:
                rows_by_kingdom[k].append(row)
    return tax, protein_activities, vec_by_acc, rows_by_kingdom


def _cell_value(r, pos):
    val = r.get(pos, '?')
    if pos == 'thioether':
        val = 'C' if val in ('C', 'C*') else '-'
    elif pos == 'Gly46' and r.get('Gly46_ss') == 'c' and val != '?':
        val = '~'
    return val


def _build_pie_map(protein_activities, tax, vec_by_acc, kingdom):
    cmap = {}
    for acc, acts in protein_activities.items():
        if acc not in vec_by_acc:
            continue
        if tax.get(acc) != kingdom:
            continue
        r = vec_by_acc[acc]
        for j, pos in enumerate(POSITIONS):
            val = _cell_value(r, pos)
            if val == '?' or val not in AA_ORDER:
                continue
            i = AA_ORDER.index(val)
            cmap.setdefault((j, i), {})
            for a in acts:
                cmap[(j, i)][a] = cmap[(j, i)].get(a, 0) + 1
    return cmap


def _draw_figure(rows_by_kingdom, tax, protein_activities, vec_by_acc,
                 mode='frequency'):
    fig, axes = plt.subplots(1, len(KINGDOMS), figsize=(3.2 * len(KINGDOMS), 7),
                              sharey=True)

    if mode == 'frequency':
        norm = mcolors.PowerNorm(gamma=0.4, vmin=0, vmax=1)
    else:
        norm = mcolors.LogNorm(vmin=1, vmax=12000)

    pie_data = []
    for ax, kingdom in zip(axes, KINGDOMS):
        krows = rows_by_kingdom[kingdom]
        n_pos = len(POSITIONS)
        n_aa = len(AA_ORDER)
        mat = np.zeros((n_aa, n_pos))

        for j, pos in enumerate(POSITIONS):
            counts = {}
            for r in krows:
                val = _cell_value(r, pos)
                if val == '?':
                    continue
                counts[val] = counts.get(val, 0) + 1
            total = sum(counts.values())
            for i, aa in enumerate(AA_ORDER):
                if mode == 'frequency':
                    mat[i, j] = counts.get(aa, 0) / total if total > 0 else 0
                else:
                    mat[i, j] = counts.get(aa, 0)

        mat[mat == 0] = np.nan
        ax.imshow(mat, aspect='auto', cmap=CMAP, norm=norm,
                  interpolation='nearest')

        pie_data.append((ax, _build_pie_map(protein_activities, tax, vec_by_acc, kingdom)))

        ax.set_title(f'{KINGDOM_LABELS[kingdom]}\n(n={len(krows):,})', fontsize=10, fontweight='bold')
        ax.xaxis.set_ticks_position('top')
        ax.set_xticks(range(n_pos))
        ax.set_xticklabels(DISPLAY, fontsize=7, rotation=45, ha='left')

        if ax == axes[0]:
            ax.set_yticks(range(n_aa))
            ax.set_yticklabels(AA_ORDER, fontsize=8, fontfamily='monospace')
        else:
            ax.tick_params(left=False)

    fig.subplots_adjust(wspace=0.05, top=0.88, bottom=0.12)

    cbar_ax = fig.add_axes([0.14, 0.06, 0.26, 0.018])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP),
                      cax=cbar_ax, orientation='horizontal')
    cb.ax.tick_params(labelsize=8)
    if mode == 'frequency':
        cb.ax.set_title('Frequency', fontsize=10)
    else:
        cb.ax.set_title('Structures per cell', fontsize=10)

    handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=ACTIVITY_COLORS[a],
                      markeredgecolor='black', markeredgewidth=0.3, markersize=7,
                      label=('Hc' if a == 'hemocyanin' else a))
               for a in ACTIVITY_ORDER]
    fig.legend(handles=handles, loc='center', bbox_to_anchor=(0.70, 0.07),
               ncol=4, fontsize=8, frameon=False, title='Characterised activity',
               handletextpad=0.3, columnspacing=1.2)

    fig.canvas.draw()
    for pax, pcmap in pie_data:
        draw_activity_pies(pax, pcmap)

    return fig


def main():
    tax, protein_activities, vec_by_acc, rows_by_kingdom = _collect_data()

    fig = _draw_figure(rows_by_kingdom, tax, protein_activities, vec_by_acc,
                       mode='frequency')
    out = Path(__file__).parent / 'taxonomy_heatmaps'
    fig.savefig(f'{out}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{out}.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}.pdf / .png")

    fig = _draw_figure(rows_by_kingdom, tax, protein_activities, vec_by_acc,
                       mode='logcount')
    out = Path(__file__).parent / 'taxonomy_heatmaps_logcount'
    fig.savefig(f'{out}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{out}.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}.pdf / .png")


if __name__ == '__main__':
    main()
