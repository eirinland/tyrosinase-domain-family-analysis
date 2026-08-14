"""
Analyze the network graph structure: what drives clustering,
what's in disconnected components, and thioether distribution.
"""

import csv
from pathlib import Path
from collections import Counter, defaultdict

import numpy as np
import networkx as nx

VECTORS = Path(__file__).parent.parent / 'position_vectors.csv'
TAXONOMY = Path(__file__).parent / 'taxonomy_lookup.csv'

POSITIONS = [
    'Gly46', 'Phe65', 'Trp68', 'Glu195', 'Asn205',
    'Arg209', 'Val218', 'Ala221', 'Phe227', 'His230', 'thioether',
]
DISPLAY = {p: p for p in POSITIONS}
DISPLAY['thioether'] = 'Cys'

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

MIN_SIZE = 3
HAMMING_THRESHOLD = 2


def get_parts(row):
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
    tax = {}
    with open(TAXONOMY) as f:
        for row in csv.DictReader(f):
            tax[row['accession']] = row['kingdom']

    vec_count = Counter()
    vec_parts = {}
    vec_activities = defaultdict(set)
    vec_kingdoms = defaultdict(Counter)

    with open(VECTORS) as f:
        for row in csv.DictReader(f):
            v = row.get('vector', '')
            parts = get_parts(row)
            vec_count[v] += 1
            vec_parts[v] = parts
            acc = row['accession']
            if acc in CHARACTERIZED:
                vec_activities[v].add(CHARACTERIZED[acc])
            k = tax.get(acc, '?')
            vec_kingdoms[v][k] += 1

    keep = {v: c for v, c in vec_count.items() if c >= MIN_SIZE}
    vectors = sorted(keep.keys())

    G = nx.Graph()
    for v in vectors:
        G.add_node(v, count=keep[v])
    for i in range(len(vectors)):
        for j in range(i + 1, len(vectors)):
            d = hamming(vec_parts[vectors[i]], vec_parts[vectors[j]])
            if d <= HAMMING_THRESHOLD:
                G.add_edge(vectors[i], vectors[j])

    components = sorted(nx.connected_components(G), key=len, reverse=True)

    print(f"Total vectors (>= {MIN_SIZE}): {len(vectors)}")
    print(f"Total structures represented: {sum(keep[v] for v in vectors)}")
    print(f"Components: {len(components)}")
    print()

    # Analyze each component
    for ci, comp in enumerate(components):
        comp_size = len(comp)
        comp_structures = sum(keep[v] for v in comp)
        comp_acts = set()
        for v in comp:
            comp_acts.update(vec_activities.get(v, set()))

        # Thioether breakdown
        thio_yes = sum(keep[v] for v in comp if vec_parts[v][POSITIONS.index('thioether')] == 'C')
        thio_no = comp_structures - thio_yes

        # Kingdom breakdown
        kingdoms = Counter()
        for v in comp:
            for k, c in vec_kingdoms[v].items():
                kingdoms[k] += c

        # Common residues at each position
        pos_consensus = []
        for pi, pos in enumerate(POSITIONS):
            res_counts = Counter()
            for v in comp:
                res_counts[vec_parts[v][pi]] += keep[v]
            top = res_counts.most_common(3)
            top_str = ', '.join(f'{r}:{c}' for r, c in top)
            pos_consensus.append((pos, top_str))

        if comp_size >= 5 or comp_acts:
            label = f"Component {ci+1}"
            if comp_acts:
                label += f" [contains: {', '.join(sorted(comp_acts))}]"

            print(f"{'='*70}")
            print(f"{label}")
            print(f"  Vectors: {comp_size}, Structures: {comp_structures:,}")
            print(f"  Thioether: C={thio_yes:,} ({thio_yes/comp_structures*100:.0f}%), no={thio_no:,} ({thio_no/comp_structures*100:.0f}%)")
            top_k = kingdoms.most_common(4)
            print(f"  Kingdoms: {', '.join(f'{k}={c:,}' for k, c in top_k)}")
            print(f"  Position signatures:")
            for pos, top_str in pos_consensus:
                print(f"    {DISPLAY[pos]:>8}: {top_str}")
            print()

        elif comp_size >= 3:
            # Brief summary for small components
            thio_pct = thio_yes / comp_structures * 100 if comp_structures > 0 else 0
            top_k = kingdoms.most_common(1)[0] if kingdoms else ('?', 0)
            # Find what makes this component unique - which positions differ from main component
            if ci > 0:
                pass  # only detail larger ones

    # Disconnected components summary
    print(f"\n{'='*70}")
    print("DISCONNECTED COMPONENTS (not in largest component)")
    print(f"{'='*70}")
    main_comp = components[0]
    for ci, comp in enumerate(components[1:], 2):
        comp_structures = sum(keep[v] for v in comp)
        comp_acts = set()
        for v in comp:
            comp_acts.update(vec_activities.get(v, set()))

        thio_yes = sum(keep[v] for v in comp if vec_parts[v][POSITIONS.index('thioether')] == 'C')

        kingdoms = Counter()
        for v in comp:
            for k, c in vec_kingdoms[v].items():
                kingdoms[k] += c
        top_k = kingdoms.most_common(2)

        # Characteristic positions: what's different from typical?
        example_v = max(comp, key=lambda v: keep[v])
        parts = vec_parts[example_v]

        act_str = f" [{', '.join(sorted(comp_acts))}]" if comp_acts else ""
        thio_str = "C" if thio_yes > comp_structures / 2 else "-"
        king_str = ', '.join(f'{k}={c}' for k, c in top_k)

        # Find minimum hamming to main component
        min_ham = 11
        for v in comp:
            for mv in main_comp:
                d = hamming(vec_parts[v], vec_parts[mv])
                if d < min_ham:
                    min_ham = d

        # Find the distinctive positions (what's unusual about this cluster)
        # Compare most common residue at each position vs the overall most common
        overall_common = {}
        for pi, pos in enumerate(POSITIONS):
            rc = Counter()
            for v in vectors:
                rc[vec_parts[v][pi]] += keep[v]
            overall_common[pi] = rc.most_common(1)[0][0]

        distinctive = []
        for pi, pos in enumerate(POSITIONS):
            rc = Counter()
            for v in comp:
                rc[vec_parts[v][pi]] += keep[v]
            top_res = rc.most_common(1)[0][0]
            if top_res != overall_common[pi]:
                distinctive.append(f"{DISPLAY[pos]}={top_res}")

        print(f"\n  Comp {ci}: {len(comp)} vectors, {comp_structures:,} structures, "
              f"thioether={thio_str}, min_hamming_to_main={min_ham}{act_str}")
        print(f"    Kingdoms: {king_str}")
        if distinctive:
            print(f"    Distinctive: {', '.join(distinctive)}")
        print(f"    Largest vector: {'-'.join(parts)} (n={keep[example_v]})")


if __name__ == '__main__':
    main()
