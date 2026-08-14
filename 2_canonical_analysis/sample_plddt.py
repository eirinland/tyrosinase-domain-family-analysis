"""Sample per-residue pLDDT at matched variable positions."""

import csv
import random
import argparse
from pathlib import Path
from collections import defaultdict

import numpy as np

from extract_position_vectors import (
    parse_atoms, get_cu, get_coord_his_ne2, get_ca_map, get_aa_map,
    kabsch, load_reference, VARIABLE_PMTYR, PMTYR_LABELS, AA3,
)

def get_plddt_map(atoms):
    """Extract per-residue pLDDT from B-factor column (AF3 convention)."""
    plddt = {}
    for a in atoms:
        if a['atom'] == 'CA' and a['seq'].isdigit():
            seq = int(a['seq'])
            if seq not in plddt:
                plddt[seq] = a.get('bfactor', None)
    return plddt


def parse_atoms_with_bfactor(path):
    atoms = []
    col_names = []
    in_atom = False
    collecting = False
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if line == 'loop_':
                collecting = False; in_atom = False; col_names = []
                continue
            if line.startswith('_atom_site.'):
                collecting = True; col_names.append(line)
                continue
            if collecting and col_names:
                collecting = False; in_atom = True
            if in_atom:
                if line.startswith('_') or line == '#' or not line:
                    break
                parts = line.split()
                if len(parts) != len(col_names):
                    continue
                row = dict(zip(col_names, parts))
                try:
                    atoms.append({
                        'elem': row.get('_atom_site.type_symbol', '').upper(),
                        'atom': row.get('_atom_site.label_atom_id', '').upper(),
                        'resn': row.get('_atom_site.label_comp_id', ''),
                        'chain': row.get('_atom_site.label_asym_id', 'A'),
                        'seq': row.get('_atom_site.label_seq_id', ''),
                        'x': float(row['_atom_site.Cartn_x']),
                        'y': float(row['_atom_site.Cartn_y']),
                        'z': float(row['_atom_site.Cartn_z']),
                        'bfactor': float(row.get('_atom_site.B_iso_or_equiv', '0')),
                    })
                except (KeyError, ValueError):
                    continue
    return atoms


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cifs', required=True)
    parser.add_argument('--pmtyr', required=True)
    parser.add_argument('--vectors-csv', required=True)
    parser.add_argument('--n-sample', type=int, default=0, help='0 = all')
    parser.add_argument('--out-csv', default=None, help='Save raw per-position pLDDT values')
    args = parser.parse_args()

    ref = load_reference(args.pmtyr)

    # Get accessions from vectors CSV
    accs = []
    with open(args.vectors_csv) as f:
        for row in csv.DictReader(f):
            if not row.get('error'):
                accs.append(row['accession'])

    cif_dir = Path(args.cifs)
    cif_map = {}
    for p in cif_dir.glob('*.cif'):
        acc = p.name.split('_taxID_')[0]
        cif_map[acc] = str(p)

    if args.n_sample > 0:
        sample_accs = random.sample(accs, min(args.n_sample, len(accs)))
    else:
        sample_accs = accs
    print(f"Processing {len(sample_accs)} structures...", flush=True)

    # Collect pLDDT per position
    pos_plddts = defaultdict(list)
    pos_below70 = defaultdict(int)
    pos_total = defaultdict(int)
    done = 0

    for acc in sample_accs:
        if acc not in cif_map:
            continue
        atoms = parse_atoms_with_bfactor(cif_map[acc])
        cu1, cu2 = get_cu(atoms)
        if cu1 is None:
            continue
        coord_his = get_coord_his_ne2(atoms, cu1, cu2)
        if len(coord_his) < 6:
            continue

        qry_ca = get_ca_map(atoms)
        qry_ne2_coords = np.array([[a['x'], a['y'], a['z']] for a in coord_his[:6]])
        R, t = kabsch(qry_ne2_coords, ref['ne2_coords'])
        qry_ca_t = {seq: R @ coord + t for seq, coord in qry_ca.items()}

        # pLDDT map from bfactor
        plddt = {}
        for a in atoms:
            if a['atom'] == 'CA' and a['seq'].isdigit():
                seq = int(a['seq'])
                if seq not in plddt:
                    plddt[seq] = a['bfactor']

        for ref_pos in VARIABLE_PMTYR:
            if ref_pos not in ref['ca']:
                continue
            ref_coord = ref['ca'][ref_pos]
            best_dist, best_seq = 999.0, None
            for qseq, qcoord in qry_ca_t.items():
                d = float(np.linalg.norm(qcoord - ref_coord))
                if d < best_dist:
                    best_dist = d
                    best_seq = qseq
            if best_seq is not None and best_seq in plddt:
                label = PMTYR_LABELS[ref_pos]
                val = plddt[best_seq]
                pos_plddts[label].append(val)
                pos_total[label] += 1
                if val < 70:
                    pos_below70[label] += 1

        done += 1
        if done % 500 == 0:
            print(f"  {done}...", flush=True)

    print(f"\nProcessed {done} structures\n")
    print(f"{'Position':<15} {'mean':>6} {'med':>6} {'<70':>6} {'%<70':>7} {'min':>6}")
    print('-' * 52)
    for ref_pos in VARIABLE_PMTYR:
        label = PMTYR_LABELS[ref_pos]
        vals = pos_plddts.get(label, [])
        if not vals:
            continue
        arr = np.array(vals)
        n_below = pos_below70[label]
        pct = 100 * n_below / len(vals)
        print(f"{label:<15} {arr.mean():6.1f} {np.median(arr):6.1f} {n_below:6d} {pct:6.1f}% {arr.min():6.1f}")

    if args.out_csv:
        with open(args.out_csv, 'w') as fout:
            fout.write("position,plddt\n")
            for ref_pos in VARIABLE_PMTYR:
                label = PMTYR_LABELS[ref_pos]
                for val in pos_plddts.get(label, []):
                    fout.write(f"{label},{val:.2f}\n")
        print(f"\nRaw values saved to {args.out_csv}")


if __name__ == '__main__':
    main()
