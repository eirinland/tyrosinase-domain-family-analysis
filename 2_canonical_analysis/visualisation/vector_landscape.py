"""
Visualize active-site residue vector landscape using MDS and t-SNE on
Hamming distance. Highlights the 37 experimentally characterized PPOs
against the uncharacterized majority.

Input:  ../position_vectors.csv
Output: vector_landscape_mds.pdf, vector_landscape_tsne.pdf
"""

import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import TSNE, MDS
from collections import Counter

CHARACTERIZED = {
    'A0A075DN54': 'CaOx',   'A0A0K1ZP03': 'TYR',   'A0A1S9DK56': 'TYR',
    'A0A261GRE4': 'TYR',    'A0A8D3X086': 'TYR',    'A0AAJ6N653': 'TYR',
    'B2ZB02':     'TYR',    'C0LU17':     'TYR',    'C7FF04':     'TYR',
    'C7FF05':     'TYR',    'D6RTB9':     'oAPO',   'G2QLD3':     'oMP',
    'P17643':     'DHICA ox','P43309':     'CaOx',   'P43311':     'CaOx',
    'Q08303':     'CaOx',   'Q2T7K1':     'TYR',    'Q2UNF9':     'oMP',
    'Q83WS2':     'TYR',    'Q93HL2':     'TYR',    'Q9ZP19':     'CaOx',
    'P14679':     'TYR',    'P11344':     'TYR',    'Q0MVP0':     'TYR',
    'Q00234':     'TYR',    'A0A261GVB1': 'TYR',    'B8NM74':     'TYR',
    'O81103':     'CaOx',   'Q6UIL3':     'CaOx',   'Q9FRX6':     'AUS',
    'P07147':     'DHICA ox','B1VTI5':     'oAPO',   'G2QC95':     'oMP',
    'A0A9P1ME48': 'oMP',    'Q2H7I7':     'oMP',    'Q2GZJ4':     'oMP',
    'G2Q526':     'oMP',
}

ACTIVITY_COLOR = {
    'TYR':      '#2166ac',
    'CaOx':     '#b2182b',
    'oMP':      '#4daf4a',
    'oAPO':     '#984ea3',
    'AUS':      '#ff7f00',
    'DHICA ox': '#e7298a',
}

MIN_SIZE = 5

def parse_vector(vec_str):
    parts = vec_str.split('-')
    cleaned = []
    for p in parts:
        if p == '':
            cleaned.append('_gap')
        else:
            cleaned.append(p)
    return cleaned

def main():
    acc_vec = {}
    with open('../position_vectors.csv') as f:
        reader = csv.reader(f)
        header = next(reader)
        vi = header.index('vector')
        for row in reader:
            acc_vec[row[0]] = row[vi]

    vec_counts = Counter(acc_vec.values())
    char_vectors = {}
    for acc, act in CHARACTERIZED.items():
        if acc in acc_vec:
            v = acc_vec[acc]
            if v not in char_vectors:
                char_vectors[v] = set()
            char_vectors[v].add(act)

    keep = {v: c for v, c in vec_counts.items() if c >= MIN_SIZE}
    print(f"Vectors >= {MIN_SIZE} members: {len(keep)}")
    print(f"Characterized among these: {sum(1 for v in keep if v in char_vectors)}")

    vectors = sorted(keep.keys())
    vec_idx = {v: i for i, v in enumerate(vectors)}
    n = len(vectors)

    parsed = [parse_vector(v) for v in vectors]
    max_len = max(len(p) for p in parsed)
    for p in parsed:
        while len(p) < max_len:
            p.append('_gap')

    hamming = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = sum(1 for a, b in zip(parsed[i], parsed[j]) if a != b)
            hamming[i, j] = d
            hamming[j, i] = d

    sizes = np.array([keep[v] for v in vectors])
    log_sizes = np.log10(sizes)
    pt_sizes = 8 + (log_sizes - log_sizes.min()) / (log_sizes.max() - log_sizes.min()) * 80

    is_char = np.array([v in char_vectors for v in vectors])
    has_thioether = np.array(['C' in v.split('-')[-1] or v.endswith('-C') or v.endswith('-C*')
                               for v in vectors])

    colors = []
    for v in vectors:
        if v in char_vectors:
            acts = char_vectors[v]
            act = sorted(acts)[0]
            colors.append(ACTIVITY_COLOR.get(act, '#333333'))
        else:
            colors.append('#cccccc')
    colors = np.array(colors)

    edgecolors = []
    for v in vectors:
        last = v.split('-')[-1]
        if last in ('C', 'C*'):
            edgecolors.append('#444444')
        else:
            edgecolors.append('none')
    edgecolors = np.array(edgecolors)

    for method_name, embedder in [
        ('MDS',  MDS(n_components=2, dissimilarity='precomputed', random_state=42,
                     normalized_stress='auto')),
        ('t-SNE', TSNE(n_components=2, metric='precomputed', random_state=42,
                       perplexity=min(30, n - 1), init='random')),
    ]:
        print(f"Running {method_name}...")
        coords = embedder.fit_transform(hamming)

        fig, ax = plt.subplots(figsize=(8, 7))

        order = np.argsort(is_char.astype(int))

        ax.scatter(coords[order][~is_char[order], 0],
                   coords[order][~is_char[order], 1],
                   s=pt_sizes[order][~is_char[order]],
                   c=colors[order][~is_char[order]],
                   edgecolors=edgecolors[order][~is_char[order]],
                   linewidths=0.4, alpha=0.5, zorder=2)

        char_mask = is_char[order]
        ax.scatter(coords[order][char_mask, 0],
                   coords[order][char_mask, 1],
                   s=pt_sizes[order][char_mask] * 1.5,
                   c=colors[order][char_mask],
                   edgecolors='black', linewidths=1.0, alpha=1.0, zorder=3,
                   marker='D')

        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='#cccccc',
                   markersize=8, label=f'Uncharacterized (n={sum(~is_char)})'),
        ]
        for act in sorted(ACTIVITY_COLOR):
            if any(act in char_vectors.get(v, set()) for v in vectors if v in char_vectors):
                legend_elements.append(
                    Line2D([0], [0], marker='D', color='w',
                           markerfacecolor=ACTIVITY_COLOR[act],
                           markeredgecolor='black', markersize=8, label=act))

        size_examples = [5, 50, 500]
        for s in size_examples:
            if s <= sizes.max():
                ls = np.log10(s)
                ps = 8 + (ls - log_sizes.min()) / (log_sizes.max() - log_sizes.min()) * 80
                legend_elements.append(
                    Line2D([0], [0], marker='o', color='w', markerfacecolor='#aaaaaa',
                           markersize=np.sqrt(ps), label=f'n={s}'))

        ax.legend(handles=legend_elements, loc='upper left', fontsize=8,
                  framealpha=0.9, title='Activity / Size', title_fontsize=9)

        n_char = sum(is_char)
        n_unchar = sum(~is_char)
        ax.set_title(f'Active-site residue vector landscape ({method_name})\n'
                     f'{n} vectors (>={MIN_SIZE} members) | '
                     f'{n_char} characterized, {n_unchar} uncharacterized',
                     fontsize=11)
        ax.set_xlabel(f'{method_name} 1')
        ax.set_ylabel(f'{method_name} 2')
        ax.set_xticks([])
        ax.set_yticks([])

        plt.tight_layout()
        out = f'vector_landscape_{method_name.lower().replace("-","")}.pdf'
        fig.savefig(out, dpi=300)
        fig.savefig(out.replace('.pdf', '.png'), dpi=200)
        plt.close(fig)
        print(f"  -> {out}")

if __name__ == '__main__':
    main()
