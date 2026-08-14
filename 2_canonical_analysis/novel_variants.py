"""
Screen active-site vectors against crystallized structures to identify
novel variants — residues at positions never observed in any crystal structure.

Reads pdb_mappings.tsv (genuine PDB accessions) and position_vectors.csv,
builds per-position "seen in crystal" sets, then flags structures with
residues not represented in any crystal.
"""

import csv
from collections import Counter, defaultdict
from pathlib import Path

PDB_MAPPINGS = Path('/cluster/work/projects/nn1003k/eirin/bioinf/trimming_test/pdb_mappings.tsv')
VECTOR_CSV = Path(__file__).parent / 'position_vectors.csv'

POSITIONS = [
    'HisCuA1+1', 'pos46', 'pos59', 'HisCuA2+1',
    'Phe_cons', 'Trp_cons', 'HisCuA3+1',
    'waterkeeper', 'HisCuB1+1', 'HisCuB2+1', 'gatekeeper', 'pos227', 'H230', 'HisCuB3+1',
    'thioether',
]

SHORT = {
    'HisCuA1+1': 'CuA His1+1', 'pos46': 'pos46', 'pos59': 'pos59',
    'HisCuA2+1': 'CuA His2+1', 'Phe_cons': 'Phe (cons)', 'Trp_cons': 'Trp (cons)',
    'HisCuA3+1': 'CuA His3+1', 'waterkeeper': 'Waterkeeper',
    'HisCuB1+1': 'CuB His1+1', 'HisCuB2+1': 'CuB His2+1',
    'gatekeeper': 'Gatekeeper', 'pos227': 'pos227', 'H230': 'H230',
    'HisCuB3+1': 'CuB His3+1', 'thioether': 'Thioether',
}


def main():
    # Load genuine PDB accessions
    pdb_accs = set()
    with open(PDB_MAPPINGS) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                parts += [''] * (5 - len(parts))
            if parts[0] == 'uniprot_acc':
                continue
            if 'FALSE_POSITIVE' in parts[4]:
                continue
            pdb_accs.add(parts[0])
    print(f"Genuine PDB accessions: {len(pdb_accs)}")

    # Load all vectors
    with open(VECTOR_CSV) as f:
        rows = list(csv.DictReader(f))
    print(f"Total structures: {len(rows)}")

    # Build per-position crystallized residue sets
    crystal_rows = [r for r in rows if r['accession'] in pdb_accs]
    print(f"Structures with crystal coverage: {len(crystal_rows)}")

    crystal_residues = {}
    for pos in POSITIONS:
        vals = set()
        for r in crystal_rows:
            v = r.get(pos, '?')
            if pos == 'thioether' and v == 'C*':
                v = 'C'
            if v not in ('?', ''):
                vals.add(v)
        crystal_residues[pos] = vals

    print(f"\nCrystallized residues per position:")
    print(f"{'Position':<16} {'Residues'}")
    print('-' * 50)
    for pos in POSITIONS:
        res = sorted(crystal_residues[pos])
        print(f"{SHORT[pos]:<16} {', '.join(res)}")

    # Screen all structures for novel residues
    novel_at_pos = defaultdict(Counter)  # pos -> {residue: count}
    novel_structures = []  # rows with at least one novel residue
    novel_count_per_struct = Counter()  # n_novel_positions -> count

    for r in rows:
        novel_positions = []
        for pos in POSITIONS:
            v = r.get(pos, '?')
            if pos == 'thioether' and v == 'C*':
                v = 'C'
            if v in ('?', ''):
                continue
            if v not in crystal_residues[pos]:
                novel_positions.append((pos, v))
                novel_at_pos[pos][v] += 1

        novel_count_per_struct[len(novel_positions)] += 1
        if novel_positions:
            novel_structures.append({
                'accession': r['accession'],
                'vector': r.get('vector', ''),
                'novel': novel_positions,
                'n_novel': len(novel_positions),
            })

    # Summary
    total_novel = len(novel_structures)
    print(f"\n{'='*60}")
    print(f"NOVEL VARIANT ANALYSIS")
    print(f"{'='*60}")
    print(f"\nStructures with >= 1 novel residue: {total_novel} ({100*total_novel/len(rows):.1f}%)")
    print(f"Structures fully covered by crystal: {len(rows)-total_novel} ({100*(len(rows)-total_novel)/len(rows):.1f}%)")

    print(f"\nNovel positions per structure:")
    for n in sorted(novel_count_per_struct):
        c = novel_count_per_struct[n]
        print(f"  {n:>2} novel positions: {c:>6} structures ({100*c/len(rows):5.1f}%)")

    print(f"\nNovel residues per position:")
    print(f"{'Position':<16} {'Crystallized':<25} {'Novel residues (count)'}")
    print('-' * 80)
    for pos in POSITIONS:
        cry = sorted(crystal_residues[pos])
        nov = novel_at_pos.get(pos, {})
        if nov:
            nov_str = ', '.join(f"{aa}:{n}" for aa, n in sorted(nov.items(), key=lambda x: -x[1]))
            print(f"{SHORT[pos]:<16} {','.join(cry):<25} {nov_str}")
        else:
            print(f"{SHORT[pos]:<16} {','.join(cry):<25} (none)")

    # Write detailed output
    out = Path(__file__).parent / 'novel_variants.csv'
    with open(out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['accession', 'vector', 'n_novel_positions', 'novel_positions', 'novel_residues'])
        for s in sorted(novel_structures, key=lambda x: -x['n_novel']):
            pos_str = ';'.join(p for p, _ in s['novel'])
            res_str = ';'.join(f"{p}={v}" for p, v in s['novel'])
            w.writerow([s['accession'], s['vector'], s['n_novel'], pos_str, res_str])
    print(f"\nDetailed output: {out}")

    # Top novel vectors by group size
    novel_vectors = Counter()
    for s in novel_structures:
        novel_vectors[s['vector']] += 1
    print(f"\nTop 20 novel vectors (by group size):")
    for v, c in novel_vectors.most_common(20):
        # Find which positions are novel
        parts = v.split('-')
        var_positions = [p for p in POSITIONS if p != 'thioether']  # thioether is last
        novel_pos = []
        for s in novel_structures:
            if s['vector'] == v:
                novel_pos = [f"{p}={res}" for p, res in s['novel']]
                break
        print(f"  {c:>5}  {v}  novel: {', '.join(novel_pos)}")


if __name__ == '__main__':
    main()
