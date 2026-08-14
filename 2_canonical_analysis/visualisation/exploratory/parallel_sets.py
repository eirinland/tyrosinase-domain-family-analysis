"""
Parallel sets diagram: 11 active-site positions as vertical axes,
flows between adjacent positions showing residue co-occurrence.
Characterized vectors colored by activity, uncharacterized in grey.
"""

import csv
from pathlib import Path
from collections import Counter, defaultdict

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
DISPLAY = [
    'Gly46', 'Phe65', 'Trp68', 'Glu195', 'Asn205',
    'Arg209', 'Val218', 'Ala221', 'Phe227', 'His230', 'Cys',
]

CHARACTERIZED = {
    'A0A075DN54': 'AUS',    'A0A0K1ZP03': 'TYR',   'A0A1S9DK56': 'TYR',
    'A0A261GRE4': 'TYR',   'A0A8D3X086': 'TYR',   'A0AAJ6N653': 'TYR',
    'B2ZB02': 'TYR',       'C0LU17': 'TYR',       'C7FF04': 'TYR',
    'C7FF05': 'TYR',       'D6RTB9': 'oAPO',      'G2QLD3': 'oMP',
    'P17643': 'DHICA ox',  'P43309': 'CaOx',      'P43311': 'CaOx',
    'Q08303': 'CaOx',      'Q2T7K1': 'TYR',       'Q2UNF9': 'oMP',
    'Q83WS2': 'TYR',       'Q93HL2': 'TYR',       'Q9ZP19': 'CaOx',
    'P14679': 'TYR',       'P11344': 'TYR',       'Q0MVP0': 'TYR',
    'Q00234': 'TYR',       'A0A261GVB1': 'TYR',   'B8NM74': 'TYR',
    'O81103': 'CaOx',      'Q6UIL3': 'CaOx',      'Q9FRX6': 'AUS',
    'P07147': 'DHICA ox',  'B1VTI5': 'oAPO',      'G2QC95': 'oMP',
    'A0A9P1ME48': 'oMP',   'Q2H7I7': 'oMP',       'Q2GZJ4': 'oMP',
    'G2Q526': 'oMP',
}

ACTIVITY_COLOR = {
    'TYR': '#2166ac',
    'CaOx': '#b2182b',
    'oMP': '#4daf4a',
    'oAPO': '#984ea3',
    'AUS': '#ff7f00',
    'DHICA ox': '#e7298a',
}

UNCHAR_COLOR = '#d0d0d0'
MIN_FRAC = 0.025  # group residues below this as "other"
CHAR_MIN_HEIGHT = 0.012  # minimum visual height for characterized flows (fraction of total)


def get_residue(row, pos):
    val = row.get(pos, '?')
    if pos == 'thioether':
        return 'C' if val in ('C', 'C*') else '-'
    return val


def draw_flow(ax, x0, x1, y0_start, y0_end, y1_start, y1_end, color, alpha):
    verts = [
        (x0, y0_start),
        ((x0 + x1) / 2, y0_start),
        ((x0 + x1) / 2, y1_start),
        (x1, y1_start),
        (x1, y1_end),
        ((x0 + x1) / 2, y1_end),
        ((x0 + x1) / 2, y0_end),
        (x0, y0_end),
        (x0, y0_start),
    ]
    codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.LINETO,
             MPath.LINETO, MPath.CURVE4, MPath.CURVE4, MPath.LINETO,
             MPath.CLOSEPOLY]
    path = MPath(verts, codes)
    patch = mpatches.PathPatch(path, facecolor=color, alpha=alpha,
                                edgecolor='none', zorder=1)
    ax.add_patch(patch)


def main():
    rows = []
    with open(VECTORS) as f:
        for row in csv.DictReader(f):
            rows.append(row)
    total = len(rows)

    # Build per-structure activity label
    acc_activity = {}
    for row in rows:
        acc = row['accession']
        if acc in CHARACTERIZED:
            acc_activity[acc] = CHARACTERIZED[acc]

    # For each position, get residue values and determine "other" grouping
    pos_vals = []  # list of lists, one per position
    pos_categories = []  # ordered category list per position
    for pos in POSITIONS:
        vals = [get_residue(r, pos) for r in rows]
        pos_vals.append(vals)
        counts = Counter(vals)
        major = {aa for aa, c in counts.items() if c / total >= MIN_FRAC}
        cats = sorted(major, key=lambda x: -counts[x])
        has_other = any(c / total < MIN_FRAC for c in counts.values())
        if has_other:
            cats.append('other')
        pos_categories.append(cats)

    # Map each structure to its category at each position
    struct_cats = []  # list of tuples
    for i, row in enumerate(rows):
        cats = []
        for j, pos in enumerate(POSITIONS):
            val = pos_vals[j][i]
            if val in pos_categories[j] and val != 'other':
                cats.append(val)
            else:
                cats.append('other')
        struct_cats.append(tuple(cats))

    # Count flows between adjacent positions, split by activity
    # activity = None for uncharacterized
    activities_order = ['TYR', 'CaOx', 'oMP', 'oAPO', 'AUS', 'DHICA ox', None]

    # For each pair of adjacent positions, count (cat_left, cat_right, activity) -> count
    pair_flows = []
    for j in range(len(POSITIONS) - 1):
        flow_counts = Counter()
        for i, row in enumerate(rows):
            acc = row['accession']
            act = acc_activity.get(acc, None)
            cl = struct_cats[i][j]
            cr = struct_cats[i][j + 1]
            flow_counts[(cl, cr, act)] += 1
        pair_flows.append(flow_counts)

    # Layout
    n_pos = len(POSITIONS)
    fig_width = 2.0 * n_pos
    fig_height = 10
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    bar_width = 0.3
    gap_frac = 0.02  # gap between categories as fraction of total
    x_positions = np.arange(n_pos) * 2.0

    # Compute y-ranges for each category at each position
    pos_y_ranges = []  # list of dicts: {cat: (y_start, y_end)}
    for j in range(n_pos):
        cats = pos_categories[j]
        counts = Counter(struct_cats[i][j] for i in range(total))
        total_gap = gap_frac * total * (len(cats) - 1)
        scale = (total - total_gap) / total if total > 0 else 1

        y_ranges = {}
        y = 0
        for cat in cats:
            height = counts.get(cat, 0) * scale
            y_ranges[cat] = (y, y + height)
            y += height + gap_frac * total
        pos_y_ranges.append(y_ranges)

    # Draw flows in two passes: uncharacterized (grey, background), then characterized on top
    # Pass 1: uncharacterized only — compute positions for all flows
    min_char_h = CHAR_MIN_HEIGHT * total  # minimum pixel height for characterized flows

    for j in range(n_pos - 1):
        x0 = x_positions[j] + bar_width / 2
        x1 = x_positions[j + 1] - bar_width / 2
        flow = pair_flows[j]

        total_gap_l = gap_frac * total * (len(pos_categories[j]) - 1)
        scale_l = (total - total_gap_l) / total
        total_gap_r = gap_frac * total * (len(pos_categories[j+1]) - 1)
        scale_r = (total - total_gap_r) / total

        left_cursors = {cat: pos_y_ranges[j][cat][0] for cat in pos_categories[j]}
        right_cursors = {cat: pos_y_ranges[j+1][cat][0] for cat in pos_categories[j+1]}

        # Draw uncharacterized first
        for cl in pos_categories[j]:
            for cr in pos_categories[j + 1]:
                count = flow.get((cl, cr, None), 0)
                if count == 0:
                    continue
                h_left = count * scale_l
                h_right = count * scale_r
                y0_s = left_cursors[cl]
                y0_e = y0_s + h_left
                left_cursors[cl] = y0_e
                y1_s = right_cursors[cr]
                y1_e = y1_s + h_right
                right_cursors[cr] = y1_e
                draw_flow(ax, x0, x1, y0_s, y0_e, y1_s, y1_e, UNCHAR_COLOR, 0.2)

        # Draw characterized on top with boosted visual height
        char_flows = []
        for act in activities_order:
            if act is None:
                continue
            for cl in pos_categories[j]:
                for cr in pos_categories[j + 1]:
                    count = flow.get((cl, cr, act), 0)
                    if count == 0:
                        continue
                    char_flows.append((cl, cr, act, count))

        for cl, cr, act, count in char_flows:
            h_left_real = count * scale_l
            h_right_real = count * scale_r
            h_left = max(h_left_real, min_char_h)
            h_right = max(h_right_real, min_char_h)

            y0_s = left_cursors[cl]
            y0_e = y0_s + h_left
            left_cursors[cl] = y0_e
            y1_s = right_cursors[cr]
            y1_e = y1_s + h_right
            right_cursors[cr] = y1_e

            color = ACTIVITY_COLOR.get(act, '#333333')
            draw_flow(ax, x0, x1, y0_s, y0_e, y1_s, y1_e, color, 0.9)

    # Draw bars at each position
    for j in range(n_pos):
        x = x_positions[j]
        cats = pos_categories[j]
        counts = Counter(struct_cats[i][j] for i in range(total))

        for cat in cats:
            y0, y1 = pos_y_ranges[j][cat]
            # Color bar by most common activity in this category, or grey
            cat_acts = Counter()
            for i in range(total):
                if struct_cats[i][j] == cat:
                    acc = rows[i]['accession']
                    act = acc_activity.get(acc, None)
                    if act is not None:
                        cat_acts[act] += 1

            bar_color = '#888888' if not cat_acts else '#666666'

            ax.add_patch(plt.Rectangle(
                (x - bar_width / 2, y0), bar_width, y1 - y0,
                facecolor=bar_color, edgecolor='white', linewidth=0.5, zorder=2))

            # Label
            height = y1 - y0
            if height / total > 0.03:
                label = cat if cat != 'other' else '...'
                pct = counts.get(cat, 0) / total * 100
                ax.text(x, (y0 + y1) / 2, f'{label}\n{pct:.0f}%',
                        ha='center', va='center', fontsize=7,
                        fontweight='bold', color='white', zorder=3)

    # Position labels at top
    for j in range(n_pos):
        max_y = max(v[1] for v in pos_y_ranges[j].values())
        ax.text(x_positions[j], max_y + total * 0.02, DISPLAY[j],
                ha='center', va='bottom', fontsize=9, fontweight='bold',
                rotation=0)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        mpatches.Patch(facecolor=UNCHAR_COLOR, alpha=0.4,
                       label=f'Uncharacterized (n={total - len(acc_activity):,})'),
    ]
    for act in ['TYR', 'CaOx', 'oMP', 'oAPO', 'AUS', 'DHICA ox']:
        n_act = sum(1 for a in acc_activity.values() if a == act)
        if n_act > 0:
            legend_elements.append(
                mpatches.Patch(facecolor=ACTIVITY_COLOR[act], alpha=0.85,
                               label=f'{act} (n={n_act})'))

    ax.legend(handles=legend_elements, loc='upper right', fontsize=8,
              framealpha=0.9, title='Activity', title_fontsize=9)

    ax.set_xlim(x_positions[0] - 1, x_positions[-1] + 1)
    ax.set_ylim(-total * 0.01, total * 1.08)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')

    ax.set_title('Active-site residue paths across 21,928 PPO structures',
                  fontsize=13, fontweight='bold', pad=20)

    plt.tight_layout()
    out = Path(__file__).parent / 'parallel_sets'
    fig.savefig(f'{out}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{out}.png', dpi=250, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}.pdf / .png")

    # Print stats
    n_char = len(acc_activity)
    print(f"\nCharacterized: {n_char} structures across {len(set(acc_activity.values()))} activities")
    print(f"Uncharacterized: {total - n_char}")
    for j, pos in enumerate(POSITIONS):
        cats = pos_categories[j]
        counts = Counter(struct_cats[i][j] for i in range(total))
        print(f"  {DISPLAY[j]}: {', '.join(f'{c}={counts[c]}' for c in cats)}")


if __name__ == '__main__':
    main()
