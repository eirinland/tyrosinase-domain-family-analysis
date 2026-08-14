"""
Network graph of active-site vectors. Nodes = unique vectors (sized by count),
edges = Hamming distance <= threshold. Colored by characterized activity.
"""

import csv
from pathlib import Path
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx

VECTORS = Path(__file__).parent.parent / 'position_vectors.csv'

POSITIONS = [
    'Gly46', 'Phe65', 'Trp68', 'Glu195', 'Asn205',
    'Arg209', 'Val218', 'Ala221', 'Phe227', 'His230', 'thioether',
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

MIN_SIZE = 3
HAMMING_THRESHOLD = 2


def get_vector_parts(row):
    parts = []
    for pos in POSITIONS:
        val = row.get(pos, '?')
        if pos == 'thioether':
            val = 'C' if val in ('C', 'C*') else '-'
        elif pos == 'Gly46' and row.get('Gly46_ss') == 'c' and val != '?':
            val = '~'
        parts.append(val)
    return tuple(parts)


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def main():
    vec_count = Counter()
    vec_parts = {}
    vec_activities = {}

    with open(VECTORS) as f:
        for row in csv.DictReader(f):
            v = row.get('vector', '')
            parts = get_vector_parts(row)
            vec_count[v] += 1
            vec_parts[v] = parts

            acc = row['accession']
            if acc in CHARACTERIZED:
                if v not in vec_activities:
                    vec_activities[v] = set()
                vec_activities[v].add(CHARACTERIZED[acc])

    keep = {v: c for v, c in vec_count.items() if c >= MIN_SIZE}
    vectors = sorted(keep.keys())
    print(f"Vectors with >= {MIN_SIZE} members: {len(vectors)}")

    G = nx.Graph()
    for v in vectors:
        G.add_node(v, count=keep[v])

    n = len(vectors)
    edge_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            d = hamming(vec_parts[vectors[i]], vec_parts[vectors[j]])
            if d <= HAMMING_THRESHOLD:
                G.add_edge(vectors[i], vectors[j], hamming=d,
                           weight=4.0 if d == 1 else 1.0)
                edge_count += 1
    edges_d1 = [(u, v) for u, v, d in G.edges(data='hamming') if d == 1]
    edges_d2 = [(u, v) for u, v, d in G.edges(data='hamming') if d == 2]
    print(f"Edges (Hamming <= {HAMMING_THRESHOLD}): {edge_count} ({len(edges_d1)} d=1, {len(edges_d2)} d=2)")

    components = list(nx.connected_components(G))
    print(f"Connected components: {len(components)}")
    print(f"Largest component: {max(len(c) for c in components)} nodes")

    pos = nx.spring_layout(G, k=1.5, iterations=120, seed=42, weight='weight')

    sizes = np.array([keep[v] for v in vectors])
    log_sizes = np.log10(np.clip(sizes, 1, None))
    pt_sizes = 15 + (log_sizes - log_sizes.min()) / (max(log_sizes.max() - log_sizes.min(), 1)) * 250

    colors = []
    markers = []
    for v in vectors:
        if v in vec_activities:
            act = sorted(vec_activities[v])[0]
            colors.append(ACTIVITY_COLOR.get(act, '#333333'))
            markers.append('D')
        else:
            colors.append('#cccccc')
            markers.append('o')

    fig, ax = plt.subplots(figsize=(12, 10))

    nx.draw_networkx_edges(G, pos, edgelist=edges_d2, ax=ax,
                           edge_color='#eeeeee', width=0.3, alpha=0.4, style='dashed')
    nx.draw_networkx_edges(G, pos, edgelist=edges_d1, ax=ax,
                           edge_color='#bbbbbb', width=0.7, alpha=0.7)

    # Uncharacterized
    unchar_idx = [i for i, v in enumerate(vectors) if v not in vec_activities]
    if unchar_idx:
        xy = np.array([pos[vectors[i]] for i in unchar_idx])
        ax.scatter(xy[:, 0], xy[:, 1],
                   s=[pt_sizes[i] for i in unchar_idx],
                   c=[colors[i] for i in unchar_idx],
                   edgecolors='#999999', linewidths=0.3, alpha=0.6, zorder=2)

    # Characterized
    char_idx = [i for i, v in enumerate(vectors) if v in vec_activities]
    if char_idx:
        xy = np.array([pos[vectors[i]] for i in char_idx])
        ax.scatter(xy[:, 0], xy[:, 1],
                   s=[pt_sizes[i] * 2 for i in char_idx],
                   c=[colors[i] for i in char_idx],
                   edgecolors='black', linewidths=1.2, alpha=1.0, zorder=3,
                   marker='D')

    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#cccccc',
               markersize=8, label=f'Uncharacterized (n={len(unchar_idx)})'),
    ]
    for act in sorted(ACTIVITY_COLOR):
        if any(act in vec_activities.get(v, set()) for v in vectors):
            legend_elements.append(
                Line2D([0], [0], marker='D', color='w',
                       markerfacecolor=ACTIVITY_COLOR[act],
                       markeredgecolor='black', markersize=8, label=act))

    for s_example in [5, 50, 500]:
        if s_example <= sizes.max():
            ls = np.log10(s_example)
            ps = 15 + (ls - log_sizes.min()) / (max(log_sizes.max() - log_sizes.min(), 1)) * 250
            legend_elements.append(
                Line2D([0], [0], marker='o', color='w', markerfacecolor='#aaaaaa',
                       markersize=np.sqrt(ps) * 0.6, label=f'n={s_example}'))

    legend_elements.append(
        Line2D([0], [0], color='#bbbbbb', linewidth=1.2, label='Hamming = 1'))
    legend_elements.append(
        Line2D([0], [0], color='#cccccc', linewidth=0.8, linestyle='dashed',
               label='Hamming = 2'))

    ax.legend(handles=legend_elements, loc='upper left', fontsize=9,
              framealpha=0.9, title='Activity / Size / Edge')

    ax.set_title(f'Vector network (Hamming distance ≤ {HAMMING_THRESHOLD}, '
                  f'min size ≥ {MIN_SIZE})\n'
                  f'{len(vectors)} vectors, {edge_count} edges, '
                  f'{len(components)} components',
                  fontsize=12)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')

    out = Path(__file__).parent / 'network_graph'
    fig.savefig(f'{out}.pdf', dpi=300, bbox_inches='tight')
    fig.savefig(f'{out}.png', dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out}.pdf / .png")


if __name__ == '__main__':
    main()
