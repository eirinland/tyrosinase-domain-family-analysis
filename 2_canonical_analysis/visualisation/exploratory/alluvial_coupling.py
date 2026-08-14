"""
Alluvial/Sankey diagrams for the most coupled position pairs (by Cramér's V).
Shows how residue identities at one position flow into residue identities at another.
Uses plotly for Sankey rendering, exports as static images.
"""

import csv
from pathlib import Path
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path as MPath

VECTORS = Path(__file__).parent.parent / 'position_vectors.csv'

POSITIONS = [
    'Gly46', 'Phe65', 'Trp68', 'Glu195', 'Asn205',
    'Arg209', 'Val218', 'Ala221', 'Phe227', 'His230', 'thioether',
]
DISPLAY = {
    'Gly46': 'Gly46', 'Phe65': 'Phe65', 'Trp68': 'Trp68',
    'Glu195': 'Glu195', 'Asn205': 'Asn205', 'Arg209': 'Arg209',
    'Val218': 'Val218', 'Ala221': 'Ala221', 'Phe227': 'Phe227',
    'His230': 'His230', 'thioether': 'Cys (thioether)',
}

AA_COLORS = {
    'G': '#66c2a5', 'A': '#8da0cb', 'V': '#e78ac3', 'L': '#a6d854',
    'I': '#ffd92f', 'F': '#e5c494', 'W': '#b3b3b3', 'Y': '#fc8d62',
    'C': '#fb9a99', 'M': '#cab2d6', 'S': '#80b1d3', 'T': '#fdb462',
    'P': '#bc80bd', 'N': '#d9d9d9', 'Q': '#ccebc5', 'D': '#ffed6f',
    'E': '#f46d43', 'H': '#74add1', 'K': '#fdae61', 'R': '#abd9e9',
    '-': '#cccccc',
}
DEFAULT_COLOR = '#999999'

MIN_FLOW = 0.01  # minimum fraction to show


def cramers_v(x, y):
    cats_x = sorted(set(x))
    cats_y = sorted(set(y))
    ix = {c: i for i, c in enumerate(cats_x)}
    iy = {c: i for i, c in enumerate(cats_y)}
    table = np.zeros((len(cats_x), len(cats_y)))
    for xi, yi in zip(x, y):
        table[ix[xi], iy[yi]] += 1
    n = table.sum()
    if n == 0:
        return 0.0
    row_sums = table.sum(axis=1, keepdims=True)
    col_sums = table.sum(axis=0, keepdims=True)
    expected = row_sums * col_sums / n
    mask = expected > 0
    chi2 = np.sum((table[mask] - expected[mask])**2 / expected[mask])
    r, k = table.shape
    denom = n * (min(r, k) - 1)
    if denom == 0:
        return 0.0
    return np.sqrt(chi2 / denom)


def draw_alluvial(ax, left_vals, right_vals, left_name, right_name):
    n = len(left_vals)
    joint = Counter(zip(left_vals, right_vals))
    left_counts = Counter(left_vals)
    right_counts = Counter(right_vals)

    left_order = sorted(left_counts.keys(), key=lambda x: -left_counts[x])
    right_order = sorted(right_counts.keys(), key=lambda x: -right_counts[x])

    # Positions
    left_y = {}
    y = 0
    gap = n * 0.005
    for aa in left_order:
        left_y[aa] = (y, y + left_counts[aa])
        y += left_counts[aa] + gap

    right_y = {}
    y = 0
    for aa in right_order:
        right_y[aa] = (y, y + right_counts[aa])
        y += right_counts[aa] + gap

    x_left = 0.0
    x_right = 1.0
    bar_w = 0.06

    # Draw flows
    left_cursor = {aa: left_y[aa][0] for aa in left_order}
    right_cursor = {aa: right_y[aa][0] for aa in right_order}

    for l_aa in left_order:
        for r_aa in right_order:
            count = joint.get((l_aa, r_aa), 0)
            if count / n < MIN_FLOW:
                continue

            y0_l = left_cursor[l_aa]
            y1_l = y0_l + count
            left_cursor[l_aa] = y1_l

            y0_r = right_cursor[r_aa]
            y1_r = y0_r + count
            right_cursor[r_aa] = y1_r

            color = AA_COLORS.get(l_aa, DEFAULT_COLOR)

            verts = [
                (x_left + bar_w, y0_l),
                (0.35, y0_l),
                (0.65, y0_r),
                (x_right - bar_w, y0_r),
                (x_right - bar_w, y1_r),
                (0.65, y1_r),
                (0.35, y1_l),
                (x_left + bar_w, y1_l),
                (x_left + bar_w, y0_l),
            ]
            codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.LINETO,
                     MPath.LINETO, MPath.CURVE4, MPath.CURVE4, MPath.LINETO,
                     MPath.CLOSEPOLY]
            path = MPath(verts, codes)
            patch = mpatches.PathPatch(path, facecolor=color, alpha=0.4,
                                        edgecolor='none')
            ax.add_patch(patch)

    # Draw bars
    for aa in left_order:
        y0, y1 = left_y[aa]
        color = AA_COLORS.get(aa, DEFAULT_COLOR)
        ax.add_patch(plt.Rectangle((x_left, y0), bar_w, y1 - y0,
                                     facecolor=color, edgecolor='white', linewidth=0.5))
        if (y1 - y0) / n > 0.03:
            ax.text(x_left - 0.02, (y0 + y1) / 2, f'{aa} ({left_counts[aa]:,})',
                    ha='right', va='center', fontsize=8, fontweight='bold')

    for aa in right_order:
        y0, y1 = right_y[aa]
        color = AA_COLORS.get(aa, DEFAULT_COLOR)
        ax.add_patch(plt.Rectangle((x_right - bar_w, y0), bar_w, y1 - y0,
                                     facecolor=color, edgecolor='white', linewidth=0.5))
        if (y1 - y0) / n > 0.03:
            ax.text(x_right + 0.02, (y0 + y1) / 2, f'{aa} ({right_counts[aa]:,})',
                    ha='left', va='center', fontsize=8, fontweight='bold')

    ax.set_xlim(-0.25, 1.25)
    max_y = max(max(v[1] for v in left_y.values()), max(v[1] for v in right_y.values()))
    ax.set_ylim(-n * 0.02, max_y + n * 0.02)
    ax.set_xticks([bar_w / 2, 1 - bar_w / 2])
    ax.set_xticklabels([left_name, right_name], fontsize=11, fontweight='bold')
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)


def main():
    with open(VECTORS) as f:
        rows = list(csv.DictReader(f))

    data = {}
    for pos in POSITIONS:
        vals = []
        for r in rows:
            v = r.get(pos, '?')
            if pos == 'thioether':
                v = 'C' if v in ('C', 'C*') else '-'
            vals.append(v)
        data[pos] = vals

    # Find top coupled pairs
    pairs = []
    for i, p1 in enumerate(POSITIONS):
        for j, p2 in enumerate(POSITIONS):
            if j <= i:
                continue
            v = cramers_v(data[p1], data[p2])
            pairs.append((p1, p2, v))
    pairs.sort(key=lambda x: -x[2])

    top4 = pairs[:4]
    print("Top 4 coupled pairs:")
    for p1, p2, v in top4:
        print(f"  {DISPLAY[p1]} × {DISPLAY[p2]}: V = {v:.3f}")

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()

    for idx, (p1, p2, v) in enumerate(top4):
        ax = axes[idx]
        draw_alluvial(ax, data[p1], data[p2], DISPLAY[p1], DISPLAY[p2])
        ax.set_title(f"{DISPLAY[p1]}  ↔  {DISPLAY[p2]}\nCramér's V = {v:.3f}",
                      fontsize=11, fontweight='bold')

    fig.suptitle('Residue co-variation at the most coupled active-site positions',
                  fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    out = Path(__file__).parent / 'alluvial_coupling'
    fig.savefig(f'{out}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{out}.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}.pdf / .png")


if __name__ == '__main__':
    main()
