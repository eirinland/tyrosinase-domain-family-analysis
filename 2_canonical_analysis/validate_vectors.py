"""
Validate position vectors: for each vector group, superpose members onto PmTYR,
extract transformed Cα coordinates at variable positions, compute mean coords,
then report mean RMSD of each member to the group centroid.
"""

import csv
import argparse
import math
from pathlib import Path
from multiprocessing import Pool
from collections import defaultdict

import numpy as np

from extract_position_vectors import (
    parse_atoms, get_cu, get_coord_his_ne2, get_ca_map,
    kabsch, load_reference, PMTYR_POSITIONS, PMTYR_ANCHORS,
    VARIABLE_PMTYR, PMTYR_LABELS,
)

_pmtyr_ref = None


def init_worker(pmtyr_path):
    global _pmtyr_ref
    _pmtyr_ref = load_reference(pmtyr_path)


def process_one(cif_path):
    """Superpose onto PmTYR, return transformed Cα coords at variable positions."""
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

    # For each variable position, find nearest transformed Cα
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

    return {
        'accession': accession,
        'coords': np.array(coords),  # (12, 3)
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cifs', required=True)
    parser.add_argument('--pmtyr', required=True)
    parser.add_argument('--vectors-csv', required=True, help='position_vectors.csv')
    parser.add_argument('--output', required=True)
    parser.add_argument('--min-group', type=int, default=5)
    parser.add_argument('--workers', type=int, default=16)
    args = parser.parse_args()

    # Read vector assignments
    acc_vector = {}
    with open(args.vectors_csv) as f:
        for row in csv.DictReader(f):
            if row.get('error'):
                continue
            acc_vector[row['accession']] = row['vector']

    # Find groups above min size
    groups = defaultdict(list)
    for acc, vec in acc_vector.items():
        groups[vec].append(acc)
    valid_groups = {v: accs for v, accs in groups.items() if len(accs) >= args.min_group}
    valid_accs = set()
    for accs in valid_groups.values():
        valid_accs.update(accs)
    print(f"Vector groups ≥{args.min_group}: {len(valid_groups)}", flush=True)
    print(f"Structures to process: {len(valid_accs)}", flush=True)

    # Find CIF files
    cif_dir = Path(args.cifs)
    cif_files = [
        str(p) for p in sorted(cif_dir.glob('*.cif'))
        if p.name.split('_taxID_')[0] in valid_accs
    ]
    print(f"CIF files matched: {len(cif_files)}", flush=True)

    # Extract coordinates
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

    # Compute per-group stats
    results = []
    for vec, accs in sorted(valid_groups.items(), key=lambda x: -len(x[1])):
        member_coords = [acc_coords[a] for a in accs if a in acc_coords]
        if len(member_coords) < 2:
            continue

        stacked = np.array(member_coords)  # (N, 12, 3)
        centroid = stacked.mean(axis=0)     # (12, 3)

        rmsds = []
        for c in member_coords:
            rmsd = float(np.sqrt(np.mean(np.sum((c - centroid)**2, axis=1))))
            rmsds.append(rmsd)

        results.append({
            'vector': vec,
            'n_structures': len(accs),
            'n_with_coords': len(member_coords),
            'mean_rmsd': f'{np.mean(rmsds):.3f}',
            'median_rmsd': f'{np.median(rmsds):.3f}',
            'max_rmsd': f'{np.max(rmsds):.3f}',
            'std_rmsd': f'{np.std(rmsds):.3f}',
        })

    with open(args.output, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=[
            'vector', 'n_structures', 'n_with_coords',
            'mean_rmsd', 'median_rmsd', 'max_rmsd', 'std_rmsd'])
        w.writeheader()
        w.writerows(results)

    print(f"\nGroups validated: {len(results)}", flush=True)
    print(f"Output: {args.output}", flush=True)

    # Summary
    rmsds_all = [float(r['mean_rmsd']) for r in results]
    print(f"\nMean RMSD across groups: {np.mean(rmsds_all):.3f} Å")
    print(f"Median: {np.median(rmsds_all):.3f} Å")

    print("\nTop 10 largest groups:")
    for r in results[:10]:
        print(f"  n={r['n_structures']:<6} mean={r['mean_rmsd']}  max={r['max_rmsd']}  {r['vector']}")

    print("\n10 worst (highest mean RMSD):")
    worst = sorted(results, key=lambda r: -float(r['mean_rmsd']))[:10]
    for r in worst:
        print(f"  n={r['n_structures']:<6} mean={r['mean_rmsd']}  max={r['max_rmsd']}  {r['vector']}")


if __name__ == '__main__':
    main()
