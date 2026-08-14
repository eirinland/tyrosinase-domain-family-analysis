"""
Per-structure outlier detection within vector groups.

For each vector group (min size ≥ 2), superpose all members onto PmTYR,
extract transformed Cα coordinates at variable positions, compute the group
centroid (mean coordinates), then flag outliers where any single position
deviates > 1.0 Å from the group centroid at that position.
"""

import csv
import argparse
from pathlib import Path
from multiprocessing import Pool
from collections import defaultdict

import numpy as np

from extract_position_vectors import (
    parse_atoms, get_cu, get_coord_his_ne2, get_ca_map,
    kabsch, load_reference, VARIABLE_PMTYR, PMTYR_LABELS,
)

_pmtyr_ref = None

POS_LABELS = [PMTYR_LABELS[p] for p in VARIABLE_PMTYR]


def init_worker(pmtyr_path):
    global _pmtyr_ref
    _pmtyr_ref = load_reference(pmtyr_path)


def process_one(cif_path):
    accession = Path(cif_path).name.split('_taxID_')[0]
    try:
        atoms = parse_atoms(cif_path)
    except Exception:
        return None

    cu1, cu2 = get_cu(atoms)
    if cu1 is None:
        return None

    coord_his = get_coord_his_ne2(atoms, cu1, cu2)
    if len(coord_his) < 6:
        return None

    qry_ca = get_ca_map(atoms)
    qry_ne2_coords = np.array([[a['x'], a['y'], a['z']] for a in coord_his[:6]])
    R, t = kabsch(qry_ne2_coords, _pmtyr_ref['ne2_coords'])

    qry_ca_t = {seq: R @ coord + t for seq, coord in qry_ca.items()}

    coords = []
    for ref_pos in VARIABLE_PMTYR:
        if ref_pos not in _pmtyr_ref['ca']:
            return None
        ref_coord = _pmtyr_ref['ca'][ref_pos]
        best_dist, best_coord = 999.0, None
        for qseq, qcoord in qry_ca_t.items():
            d = float(np.linalg.norm(qcoord - ref_coord))
            if d < best_dist:
                best_dist = d
                best_coord = qcoord
        coords.append(best_coord)

    return {'accession': accession, 'coords': np.array(coords)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cifs', required=True)
    parser.add_argument('--pmtyr', required=True)
    parser.add_argument('--vectors-csv', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--min-group', type=int, default=2)
    parser.add_argument('--multi-cutoff', type=float, default=1.25,
                        help='Per-position cutoff for multi-position outliers (default 1.25 Å)')
    parser.add_argument('--multi-min-pos', type=int, default=2,
                        help='Min positions exceeding multi-cutoff (default 2)')
    parser.add_argument('--single-cutoff', type=float, default=2.0,
                        help='Per-position cutoff for single dramatic outliers (default 2.0 Å)')
    parser.add_argument('--workers', type=int, default=16)
    args = parser.parse_args()

    acc_vector = {}
    with open(args.vectors_csv) as f:
        for row in csv.DictReader(f):
            if row.get('error'):
                continue
            acc_vector[row['accession']] = row['vector']

    groups = defaultdict(list)
    for acc, vec in acc_vector.items():
        groups[vec].append(acc)
    valid_groups = {v: accs for v, accs in groups.items() if len(accs) >= args.min_group}
    valid_accs = set()
    for accs in valid_groups.values():
        valid_accs.update(accs)
    print(f"Vector groups ≥{args.min_group}: {len(valid_groups)}", flush=True)
    print(f"Structures to process: {len(valid_accs)}", flush=True)

    cif_dir = Path(args.cifs)
    cif_files = [
        str(p) for p in sorted(cif_dir.glob('*.cif'))
        if p.name.split('_taxID_')[0] in valid_accs
    ]
    print(f"CIF files matched: {len(cif_files)}", flush=True)

    acc_coords = {}
    done = 0
    with Pool(processes=args.workers,
              initializer=init_worker,
              initargs=(args.pmtyr,)) as pool:
        for result in pool.imap_unordered(process_one, cif_files, chunksize=50):
            if result is not None:
                acc_coords[result['accession']] = result['coords']
            done += 1
            if done % 2000 == 0:
                print(f"  {done}/{len(cif_files)}...", flush=True)

    print(f"Coordinates extracted: {len(acc_coords)}", flush=True)

    multi_cutoff = args.multi_cutoff
    multi_min = args.multi_min_pos
    single_cutoff = args.single_cutoff

    rows = []
    counts = {'multi_position': 0, 'single_position': 0, 'ok': 0}
    for vec, accs in sorted(valid_groups.items(), key=lambda x: -len(x[1])):
        member_data = [(a, acc_coords[a]) for a in accs if a in acc_coords]
        if len(member_data) < 2:
            continue

        stacked = np.array([c for _, c in member_data])  # (N, 10, 3)
        centroid = stacked.mean(axis=0)                   # (10, 3)

        for i, (acc, coords) in enumerate(member_data):
            pos_dists = np.sqrt(np.sum((coords - centroid)**2, axis=1))  # (10,)
            rmsd = float(np.sqrt(np.mean(pos_dists**2)))

            flagged = []
            for j, d in enumerate(pos_dists):
                if d > multi_cutoff:
                    flagged.append(f'{POS_LABELS[j]}={d:.2f}')

            max_dist = float(pos_dists.max())
            if len(flagged) >= multi_min:
                outlier_type = 'multi_position'
            elif max_dist > single_cutoff:
                outlier_type = 'single_position'
            else:
                outlier_type = 'ok'
            counts[outlier_type] += 1

            row = {
                'accession': acc,
                'vector': vec,
                'group_size': len(member_data),
                'rmsd_to_centroid': f'{rmsd:.4f}',
                'max_pos_dist': f'{max_dist:.3f}',
                'n_flagged_positions': len(flagged),
                'flagged_positions': ';'.join(flagged) if flagged else '',
                'outlier_type': outlier_type,
            }
            for j, label in enumerate(POS_LABELS):
                row[f'{label}_dist'] = f'{pos_dists[j]:.3f}'

            rows.append(row)

    rows.sort(key=lambda r: r['accession'])

    fieldnames = [
        'accession', 'vector', 'group_size',
        'rmsd_to_centroid', 'max_pos_dist',
        'n_flagged_positions', 'flagged_positions', 'outlier_type',
    ]
    for label in POS_LABELS:
        fieldnames.append(f'{label}_dist')

    with open(args.output, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)

    total = len(rows)
    print(f"\nDone. {total} structures in {len(valid_groups)} groups.", flush=True)
    print(f"Multi-position (≥{multi_min} pos > {multi_cutoff} Å): {counts['multi_position']:>5} ({100*counts['multi_position']/total:.1f}%)")
    print(f"Single-position (max pos > {single_cutoff} Å):    {counts['single_position']:>5} ({100*counts['single_position']/total:.1f}%)")
    print(f"Not outlier:                            {counts['ok']:>5} ({100*counts['ok']/total:.1f}%)")
    print(f"Output: {args.output}")


if __name__ == '__main__':
    main()
