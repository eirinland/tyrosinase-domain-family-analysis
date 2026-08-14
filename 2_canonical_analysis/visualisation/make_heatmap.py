"""
Generate active-site position frequency heatmap from position_vectors.csv.
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
PDB_MAPPINGS = Path('/cluster/work/projects/nn1003k/eirin/bioinf/trimming_test/pdb_mappings.tsv')

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

ACTIVITY_COLORS = {
    'TYR': '#e41a1c',
    'CaOx': '#377eb8',
    'AUS': '#ff7f00',
    'DHICA ox': '#984ea3',
    'oAPO': '#4daf4a',
    'oMP': '#a65628',
}

PROTEIN_ACTIVITIES = {
    # Crystal structures
    'A0A075DN54': ['AUS'],
    'A0A0K1ZP03': ['TYR'],
    'A0A1S9DK56': ['TYR'],
    'A0A261GRE4': ['TYR'],
    'A0A8D3X086': ['TYR'],
    'A0AAJ6N653': ['TYR'],
    'B2ZB02': ['TYR'],
    'C0LU17': ['TYR'],
    'C7FF04': ['TYR'],
    'C7FF05': ['TYR'],
    'D6RTB9': ['oAPO'],
    'G2QLD3': ['CaOx', 'oMP'],
    'P17643': ['DHICA ox'],
    'P43309': ['CaOx'],
    'P43311': ['CaOx'],
    'Q08303': ['CaOx'],
    'Q2T7K1': ['TYR'],
    'Q2UNF9': ['CaOx', 'oMP'],
    'Q83WS2': ['TYR'],
    'Q93HL2': ['TYR'],
    'Q9ZP19': ['CaOx'],
    # Experimentally characterized (non-crystal)
    'P14679': ['TYR'],
    'P11344': ['TYR'],
    'Q0MVP0': ['TYR'],
    'Q00234': ['TYR'],
    'A0A261GVB1': ['TYR'],
    'B8NM74': ['TYR'],
    'O81103': ['CaOx'],
    'Q6UIL3': ['CaOx'],
    'Q9FRX6': ['AUS'],
    'P07147': ['DHICA ox'],
    'B1VTI5': ['oAPO'],
    # oMP (+ CaOx) — published activity not yet in UniProt
    'G2QC95': ['CaOx', 'oMP'],
    'A0A9P1ME48': ['CaOx', 'oMP'],
    'Q2H7I7': ['CaOx', 'oMP'],
    'Q2GZJ4': ['CaOx', 'oMP'],
}


def load_crystal_counts(rows):
    pdb_accs = set()
    with open(PDB_MAPPINGS) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                parts += [''] * (5 - len(parts))
            if parts[0] == 'uniprot_acc' or 'FALSE_POSITIVE' in parts[4]:
                continue
            pdb_accs.add(parts[0])
    crystal_rows = [r for r in rows if r['accession'] in pdb_accs]
    counts = {}
    for pos in POSITION_LABELS:
        c = {}
        for r in crystal_rows:
            val = r.get(pos, '?')
            if pos == 'thioether' and val == 'C*':
                val = 'C'
            elif pos == 'Gly46' and r.get('Gly46_ss') == 'c':
                val = '~'
            if val in ('?', ''):
                continue
            c[val] = c.get(val, 0) + 1
        counts[pos] = c
    return counts


def main():
    with open(VECTOR_CSV) as f:
        rows = list(csv.DictReader(f))

    crystal_counts = load_crystal_counts(rows)

    freqs = {}
    raw_counts = {}
    for pos in POSITION_LABELS:
        c = {}
        for r in rows:
            val = r.get(pos, '?')
            if val == '?':
                continue
            if pos == 'thioether':
                key = 'C' if val in ('C', 'C*') else '-'
            elif pos == 'Gly46' and r.get('Gly46_ss') == 'c':
                key = '~'
            else:
                key = val
            c[key] = c.get(key, 0) + 1
        total = sum(c.values())
        freqs[pos] = {aa: c.get(aa, 0) / total if total > 0 else 0 for aa in AA_ORDER}
        raw_counts[pos] = c

    n_pos = len(POSITION_LABELS)
    n_aa = len(AA_ORDER)
    mat = np.zeros((n_aa, n_pos))
    for j, pos in enumerate(POSITION_LABELS):
        for i, aa in enumerate(AA_ORDER):
            mat[i, j] = freqs[pos].get(aa, 0)

    fig, ax = plt.subplots(figsize=(6, 7))

    norm = mcolors.PowerNorm(gamma=0.4, vmin=0, vmax=1)
    im = ax.imshow(mat, aspect='auto', cmap='YlOrRd', norm=norm, interpolation='nearest')

    ax.set_yticks(range(n_aa))
    ax.set_yticklabels(AA_ORDER, fontsize=10, fontfamily='monospace')
    ax.set_ylabel('')

    ax.xaxis.set_ticks_position('top')
    ax.set_xticks(range(n_pos))
    ax.set_xticklabels(DISPLAY_LABELS, fontsize=11, rotation=0)

    stroke = pe.withStroke(linewidth=0.8, foreground='black')
    for i, label in enumerate(ax.get_xticklabels()):
        label.set_color(COL_COLORS[i])
        label.set_fontweight('bold')
        label.set_path_effects([stroke])

    ax.set_title('')

    for i in range(n_aa):
        for j in range(n_pos):
            val = mat[i, j]
            if val >= 0.05:
                color = 'white' if val > 0.5 else 'black'
                ax.text(j, i, f'{val:.0%}', ha='center', va='center',
                        fontsize=7, color=color, fontweight='bold')
            aa = AA_ORDER[i]
            pos = POSITION_LABELS[j]
            n_af3 = raw_counts.get(pos, {}).get(aa, 0)
            n_xtal = crystal_counts.get(pos, {}).get(aa, 0)
            if n_af3 > 3 and n_xtal == 0:
                ax.text(j, i, '★', ha='left', va='top', fontsize=6,
                        color='white',
                        path_effects=[pe.withStroke(linewidth=0.6, foreground='black')])

    for sep in [5.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 19.5, 20.5]:
        ax.axhline(sep, color='grey', linewidth=0.5, alpha=0.5)

    # Activity circles overlay
    acc_to_row = {r['accession']: r for r in rows}
    for acc, activities in PROTEIN_ACTIVITIES.items():
        if acc not in acc_to_row:
            continue
        r = acc_to_row[acc]
        for j, pos in enumerate(POSITION_LABELS):
            val = r.get(pos, '?')
            if pos == 'thioether':
                val = 'C' if val in ('C', 'C*') else '-'
            elif pos == 'Gly46' and r.get('Gly46_ss') == 'c':
                val = '~'
            if val == '?' or val not in AA_ORDER:
                continue
            i = AA_ORDER.index(val)

            # Gly46: open circle if loop context
            filled = True
            if pos == 'Gly46' and r.get('Gly46_ss') == 'c':
                filled = False

            for k, act in enumerate(activities):
                color = ACTIVITY_COLORS[act]
                offset = (k - (len(activities) - 1) / 2) * 0.15
                if filled:
                    ax.plot(j + offset, i, 'o', color=color, markersize=4,
                            markeredgecolor='black', markeredgewidth=0.3)
                else:
                    ax.plot(j + offset, i, 'o', color='none', markersize=4,
                            markeredgecolor=color, markeredgewidth=1.0)

    fig.colorbar(im, ax=ax, orientation='horizontal', shrink=0.6, pad=0.06,
                 aspect=30, label='Frequency')

    out = Path(__file__).parent / 'position_frequency_heatmap.png'
    fig.savefig(out, dpi=600, bbox_inches='tight')
    print(f'Saved: {out}')

    out_pdf = Path(__file__).parent / 'position_frequency_heatmap.pdf'
    fig.savefig(out_pdf, bbox_inches='tight')
    print(f'Saved: {out_pdf}')


if __name__ == '__main__':
    main()
