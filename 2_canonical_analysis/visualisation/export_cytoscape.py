"""
Export vector network for Cytoscape.
Outputs:
  - network_nodes.csv: one row per vector (node attributes)
  - network_edges.csv: all edges up to MAX_HAMMING (edge attributes)
"""

import csv
from pathlib import Path
from collections import Counter, defaultdict

VECTORS = Path(__file__).parent.parent / 'position_vectors.csv'
TAXONOMY = Path(__file__).parent / 'taxonomy_lookup.csv'

POSITIONS = [
    'Gly46', 'Phe65', 'Trp68', 'Glu195', 'Asn205',
    'Arg209', 'Val218', 'Ala221', 'Phe227', 'His230', 'thioether',
]

TAXONOMY_KINGDOMS = ['Fungi', 'Animals', 'Plants', 'Bacteria', 'Oomycota']
KINGDOM_COLORS = {
    'Fungi': '#EAD2B9', 'Animals': '#83B3F4', 'Plants': '#5A9163',
    'Bacteria': '#A9D8E8', 'Oomycota': '#BAA3C5',
}

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

ACTIVITY_COLORS = {
    'TYR': '#58A0DF', 'CaOx': '#9BD0F0', 'AUS': '#B5EAE4', 'oAPO': '#BEF0BE',
    'oMP': '#7BC486', 'DHICA ox': '#FFDBBB', 'DCT': '#F9BE9E', 'hemocyanin': '#D0B5E4',
}

MIN_SIZE = 3
MAX_HAMMING = 5
HD_MAX_COLOR = 7


def _hex_to_rgb(h):
    h = h.lstrip('#')
    return tuple(int(h[i:i+2], 16) / 255 for i in (0, 2, 4))


def _rgb_to_hex(r, g, b):
    return f'#{int(r*255):02X}{int(g*255):02X}{int(b*255):02X}'


def _lighten(hex_color, fraction):
    r, g, b = _hex_to_rgb(hex_color)
    r = r + (1.0 - r) * fraction
    g = g + (1.0 - g) * fraction
    b = b + (1.0 - b) * fraction
    return _rgb_to_hex(r, g, b)


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


def differing_positions(a, b):
    diffs = []
    for pos, x, y in zip(POSITIONS, a, b):
        if x != y:
            diffs.append(f"{pos}:{x}>{y}")
    return '; '.join(diffs)


def main():
    tax = {}
    if TAXONOMY.exists():
        with open(TAXONOMY) as f:
            for row in csv.DictReader(f):
                tax[row['accession']] = {
                    'kingdom': row.get('kingdom', '?'),
                    'phylum': row.get('phylum', '?'),
                    'genus': row.get('genus', '?'),
                }

    vec_count = Counter()
    vec_parts = {}
    vec_activities = defaultdict(set)
    vec_kingdoms = defaultdict(Counter)
    vec_genera = defaultdict(Counter)
    vec_accessions = defaultdict(list)

    with open(VECTORS) as f:
        for row in csv.DictReader(f):
            v = row.get('vector', '')
            parts = get_parts(row)
            vec_count[v] += 1
            vec_parts[v] = parts
            acc = row['accession']
            vec_accessions[v].append(acc)
            if acc in CHARACTERIZED:
                vec_activities[v].add(CHARACTERIZED[acc])
            t = tax.get(acc, {})
            vec_kingdoms[v][t.get('kingdom', '?')] += 1
            vec_genera[v][t.get('genus', '?')] += 1

    keep = {v: c for v, c in vec_count.items() if c >= MIN_SIZE}
    vectors = sorted(keep.keys())
    print(f"Vectors with >= {MIN_SIZE} members: {len(vectors)}")

    char_vectors = {v: sorted(vec_activities[v]) for v in vec_activities if v in keep}
    print(f"Characterized vectors in network: {len(char_vectors)}")

    def nearest_characterized(parts):
        best_hd, best_act = 999, 'none'
        for cv, acts in char_vectors.items():
            cv_parts = vec_parts[cv]
            d = hamming(parts, cv_parts)
            if d < best_hd:
                best_hd = d
                best_act = acts[0]
        return best_hd, best_act

    # --- Nodes ---
    outdir = Path(__file__).parent
    node_file = outdir / 'network_nodes.csv'
    frac_cols = [f'frac_{k}' for k in TAXONOMY_KINGDOMS]
    node_fields = ['vector', 'count', 'log_count', 'activity', 'characterized',
                   'has_thioether', 'nearest_activity', 'min_hd', 'hd_color',
                   ] + frac_cols + POSITIONS + [
                   'top_kingdom', 'kingdom_breakdown', 'top_genus', 'example_accessions']

    import math
    with open(node_file, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=node_fields)
        w.writeheader()
        for v in vectors:
            parts = vec_parts[v]
            acts = sorted(vec_activities.get(v, set()))
            kingdoms = vec_kingdoms[v]
            genera = vec_genera[v]
            top_k = kingdoms.most_common(1)[0][0] if kingdoms else '?'
            k_str = '; '.join(f'{k}={c}' for k, c in kingdoms.most_common(5))
            top_g = genera.most_common(1)[0][0] if genera else '?'
            examples = '; '.join(vec_accessions[v][:5])

            if acts:
                hd_val = 0
                near_act = acts[0]
                color = ACTIVITY_COLORS.get(near_act, '#CCCCCC')
            else:
                hd_val, near_act = nearest_characterized(parts)
                if hd_val >= 4:
                    color = '#CCCCCC'
                else:
                    frac = (hd_val / 3.0) * 0.75
                    base = ACTIVITY_COLORS.get(near_act, '#CCCCCC')
                    color = _lighten(base, frac)

            row = {
                'vector': v,
                'count': keep[v],
                'log_count': f'{math.log10(max(keep[v], 1)):.2f}',
                'activity': '|'.join(acts) if acts else 'none',
                'characterized': 'yes' if acts else 'no',
                'has_thioether': 'yes' if parts[POSITIONS.index('thioether')] == 'C' else 'no',
                'nearest_activity': near_act,
                'min_hd': hd_val,
                'hd_color': color,
                'top_kingdom': top_k,
                'kingdom_breakdown': k_str,
                'top_genus': top_g,
                'example_accessions': examples,
            }
            ktotal = sum(kingdoms.values())
            for k in TAXONOMY_KINGDOMS:
                row[f'frac_{k}'] = f'{kingdoms.get(k, 0) / ktotal:.3f}' if ktotal else '0'
            for pos, val in zip(POSITIONS, parts):
                row[pos] = val
            w.writerow(row)
    print(f"Nodes: {node_file} ({len(vectors)} rows)")

    # --- Edges ---
    edge_file = outdir / 'network_edges.csv'
    edge_count = 0
    with open(edge_file, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['source', 'target', 'hamming', 'differences'])
        for i in range(len(vectors)):
            for j in range(i + 1, len(vectors)):
                d = hamming(vec_parts[vectors[i]], vec_parts[vectors[j]])
                if d <= MAX_HAMMING:
                    diffs = differing_positions(vec_parts[vectors[i]], vec_parts[vectors[j]])
                    w.writerow([vectors[i], vectors[j], d, diffs])
                    edge_count += 1
    print(f"Edges: {edge_file} ({edge_count:,} rows, Hamming <= {MAX_HAMMING})")

    # Edge count per threshold
    for t in range(1, MAX_HAMMING + 1):
        print(f"  Hamming <= {t}: filter in Cytoscape with 'hamming <= {t}'")


if __name__ == '__main__':
    main()
