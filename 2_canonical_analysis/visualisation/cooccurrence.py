"""
Co-occurrence analysis of active-site position vectors.

1. Cramér's V heatmap across all 14 positions
2. Top enriched/depleted residue pair combinations for strongly coupled positions
"""

import csv
from pathlib import Path
from itertools import combinations

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

VECTOR_CSV = Path(__file__).parent.parent / 'position_vectors.csv'

POSITIONS = [
    'Gly46', 'Phe65', 'Trp68', 'Glu195', 'Asn205',
    'Arg209', 'Val218', 'Ala221', 'Phe227', 'His230',
    'thioether',
]

SHORT_LABELS = [
    'Gly46', 'Phe65', 'Trp68', 'Glu195', 'Asn205',
    'Arg209', 'Val218', 'Ala221', 'Phe227', 'His230',
    'Thioether',
]


def cramers_v(x, y):
    """Cramér's V from two categorical arrays."""
    from collections import Counter
    # Build contingency table
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
    # Avoid division by zero
    mask = expected > 0
    chi2 = np.sum((table[mask] - expected[mask])**2 / expected[mask])
    r, k = table.shape
    denom = n * (min(r, k) - 1)
    if denom == 0:
        return 0.0
    return np.sqrt(chi2 / denom)


def log_odds_enrichment(x, y, min_count=20):
    """Compute log2 odds ratios for each (residue_x, residue_y) pair."""
    from collections import Counter
    n = len(x)
    joint = Counter(zip(x, y))
    mx = Counter(x)
    my = Counter(y)
    results = []
    for (a, b), obs in joint.items():
        exp = mx[a] * my[b] / n
        if obs < min_count or exp < 5:
            continue
        ratio = obs / exp
        log2r = np.log2(ratio)
        results.append((a, b, obs, exp, ratio, log2r))
    results.sort(key=lambda r: -abs(r[5]))
    return results


def main():
    with open(VECTOR_CSV) as f:
        rows = list(csv.DictReader(f))

    # Extract column data, map thioether C/C* -> C
    data = {}
    for pos in POSITIONS:
        vals = []
        for r in rows:
            v = r.get(pos, '?')
            if pos == 'thioether' and v == 'C*':
                v = 'C'
            elif pos == 'Gly46' and r.get('Gly46_ss') == 'c' and v != '?':
                v = '~'
            vals.append(v)
        data[pos] = vals

    n = len(POSITIONS)

    # --- 1. Cramér's V matrix ---
    vmat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            v = cramers_v(data[POSITIONS[i]], data[POSITIONS[j]])
            vmat[i, j] = v
            vmat[j, i] = v
        vmat[i, i] = 1.0

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(vmat, cmap='YlOrRd', vmin=0, vmax=0.8, interpolation='nearest')
    ax.set_xticks(range(n))
    ax.set_xticklabels(SHORT_LABELS, rotation=45, ha='right', fontsize=9)
    ax.set_yticks(range(n))
    ax.set_yticklabels(SHORT_LABELS, fontsize=9)

    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            val = vmat[i, j]
            if val >= 0.05:
                color = 'white' if val > 0.5 else 'black'
                ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                        fontsize=7, color=color, fontweight='bold')

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Cramér's V", fontsize=11)
    ax.set_title("Pairwise association between active-site positions (Cramér's V)", fontsize=12)

    plt.tight_layout()
    out = Path(__file__).parent / 'cramers_v_heatmap.png'
    fig.savefig(out, dpi=300, bbox_inches='tight')
    print(f'Saved: {out}')
    out_pdf = Path(__file__).parent / 'cramers_v_heatmap.pdf'
    fig.savefig(out_pdf, bbox_inches='tight')
    print(f'Saved: {out_pdf}')
    plt.close()

    # --- 2. Top coupled pairs: enrichment analysis ---
    # Find top pairs by Cramér's V
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            pairs.append((i, j, vmat[i, j]))
    pairs.sort(key=lambda x: -x[2])

    print(f'\n{"="*70}')
    print('Top 15 coupled position pairs (by Cramér\'s V)')
    print(f'{"="*70}')

    for rank, (i, j, v) in enumerate(pairs[:15], 1):
        pi, pj = POSITIONS[i], POSITIONS[j]
        print(f'\n--- #{rank}: {SHORT_LABELS[i]} × {SHORT_LABELS[j]}  (V = {v:.3f}) ---')
        enrichments = log_odds_enrichment(data[pi], data[pj], min_count=20)
        if not enrichments:
            print('  (no combinations with count ≥ 20)')
            continue
        print(f'  {"AA_1":>5} {"AA_2":>5} {"obs":>7} {"exp":>7} {"ratio":>7} {"log2OR":>7}')
        for a, b, obs, exp, ratio, log2r in enrichments[:10]:
            tag = '+++' if log2r > 1 else '++' if log2r > 0.5 else '--' if log2r < -1 else '-' if log2r < -0.5 else ''
            print(f'  {a:>5} {b:>5} {obs:>7.0f} {exp:>7.1f} {ratio:>7.2f} {log2r:>+7.2f} {tag}')

    # --- 3. Thioether-specific enrichment table ---
    print(f'\n{"="*70}')
    print('Thioether co-occurrence enrichment (top per position)')
    print(f'{"="*70}')
    thio_idx = POSITIONS.index('thioether')
    for i in range(n):
        if i == thio_idx:
            continue
        v = vmat[i, thio_idx]
        if v < 0.1:
            continue
        pi = POSITIONS[i]
        print(f'\n{SHORT_LABELS[i]} (V = {v:.3f}):')
        enrichments = log_odds_enrichment(data[pi], data[POSITIONS[thio_idx]], min_count=20)
        for a, b, obs, exp, ratio, log2r in enrichments[:5]:
            direction = 'enriched' if log2r > 0 else 'depleted'
            print(f'  {a} + thio={b}: {obs:.0f} obs vs {exp:.1f} exp ({ratio:.2f}x, {direction})')


if __name__ == '__main__':
    main()
