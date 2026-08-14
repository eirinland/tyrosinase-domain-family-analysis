"""
UpSet-style plot showing co-occurrence of non-canonical residues.
For each structure, flag which of the 11 positions deviate from canonical,
then show the most common deviation combinations.
"""

import csv
from pathlib import Path
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

VECTORS = Path(__file__).parent.parent / 'position_vectors.csv'

POSITIONS = [
    'Gly46', 'Phe65', 'Trp68', 'Glu195', 'Asn205',
    'Arg209', 'Val218', 'Ala221', 'Phe227', 'His230', 'thioether',
]
DISPLAY = [
    'Gly46', 'Phe65', 'Trp68', 'Glu195', 'Asn205',
    'Arg209', 'Val218', 'Ala221', 'Phe227', 'His230', 'Cys',
]

CANONICAL = {
    'Gly46': 'G', 'Phe65': 'F', 'Trp68': 'W', 'Glu195': 'E',
    'Asn205': 'N', 'Arg209': 'R', 'Val218': 'V', 'Ala221': 'A',
    'Phe227': 'F', 'His230': 'H', 'thioether': '-',
}

TOP_N = 30


def main():
    combos = Counter()
    total = 0

    with open(VECTORS) as f:
        for row in csv.DictReader(f):
            total += 1
            devs = []
            for pos in POSITIONS:
                val = row.get(pos, '?')
                if pos == 'thioether':
                    val = 'C' if val in ('C', 'C*') else '-'
                if val != CANONICAL[pos] and val != '?':
                    devs.append(pos)
            key = tuple(devs)
            combos[key] += 1

    top = combos.most_common(TOP_N)

    # --- Build UpSet plot ---
    n_sets = len(POSITIONS)
    n_combos = len(top)

    fig = plt.figure(figsize=(12, 8))

    # Grid: top = bar chart, bottom = dot matrix
    gs = fig.add_gridspec(2, 2, height_ratios=[2, 1.5], width_ratios=[1, 5],
                          hspace=0.05, wspace=0.02)

    # Bar chart (top right)
    ax_bar = fig.add_subplot(gs[0, 1])
    counts = [c for _, c in top]
    x = np.arange(n_combos)
    bars = ax_bar.bar(x, counts, color='#2166ac', width=0.6)
    ax_bar.set_ylabel('Structures', fontsize=10)
    ax_bar.set_xlim(-0.5, n_combos - 0.5)
    ax_bar.set_xticks([])
    ax_bar.spines['top'].set_visible(False)
    ax_bar.spines['right'].set_visible(False)

    for i, c in enumerate(counts):
        ax_bar.text(i, c + total * 0.003, f'{c:,}', ha='center', va='bottom',
                    fontsize=6, rotation=45)

    # Dot matrix (bottom right)
    ax_dot = fig.add_subplot(gs[1, 1])
    ax_dot.set_xlim(-0.5, n_combos - 0.5)
    ax_dot.set_ylim(-0.5, n_sets - 0.5)
    ax_dot.set_xticks([])
    ax_dot.set_yticks(range(n_sets))
    ax_dot.set_yticklabels([DISPLAY[i] for i in range(n_sets)], fontsize=9)
    ax_dot.invert_yaxis()
    ax_dot.spines['top'].set_visible(False)
    ax_dot.spines['right'].set_visible(False)
    ax_dot.spines['bottom'].set_visible(False)

    for j in range(n_sets):
        ax_dot.axhline(j, color='#eeeeee', linewidth=0.5, zorder=0)

    for i, (combo, _) in enumerate(top):
        active = [POSITIONS.index(p) for p in combo]

        # Draw gray dots for inactive
        for j in range(n_sets):
            if j not in active:
                ax_dot.plot(i, j, 'o', color='#dddddd', markersize=6, zorder=1)

        # Draw dark dots for active, connected by line
        if len(active) > 0:
            for j in active:
                ax_dot.plot(i, j, 'o', color='#333333', markersize=7, zorder=3)
            if len(active) > 1:
                ax_dot.plot([i, i], [min(active), max(active)], '-',
                            color='#333333', linewidth=1.5, zorder=2)

    # Set sizes (bottom left) — how many structures have each individual deviation
    ax_size = fig.add_subplot(gs[1, 0])

    with open(VECTORS) as f:
        rows = list(csv.DictReader(f))

    set_sizes = []
    for pos in POSITIONS:
        canon = CANONICAL[pos]
        n = sum(1 for r in rows
                if (r.get(pos, '?') if pos != 'thioether'
                    else ('C' if r.get(pos, '?') in ('C', 'C*') else '-')) != canon
                and r.get(pos, '?') != '?')
        set_sizes.append(n)

    ax_size.barh(range(n_sets), set_sizes, color='#999999', height=0.6)
    ax_size.set_ylim(-0.5, n_sets - 0.5)
    ax_size.invert_yaxis()
    ax_size.invert_xaxis()
    ax_size.set_yticks([])
    ax_size.set_xlabel('Set size', fontsize=9)
    ax_size.spines['top'].set_visible(False)
    ax_size.spines['right'].set_visible(False)

    for i, s in enumerate(set_sizes):
        ax_size.text(s + total * 0.01, i, f'{s:,}', ha='left', va='center', fontsize=7)

    # Hide top-left
    ax_empty = fig.add_subplot(gs[0, 0])
    ax_empty.axis('off')

    fig.suptitle('Co-occurrence of non-canonical active-site residues (top 30 combinations)',
                  fontsize=12, fontweight='bold', y=0.95)

    out = Path(__file__).parent / 'upset_substitutions'
    fig.savefig(f'{out}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{out}.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}.pdf / .png")

    # Print top combos
    print(f"\nTop {TOP_N} deviation combinations:")
    for combo, c in top:
        names = [DISPLAY[POSITIONS.index(p)] for p in combo] if combo else ['(all canonical)']
        print(f"  {c:>5}  {' + '.join(names)}")


if __name__ == '__main__':
    main()
