#!/usr/bin/env python3
"""Detailed Cu coordination analysis for the non-canonical pool (1,036 structures).

For each structure with ≥2 Cu atoms:
  - Find all residues with any heavy atom within 3.0 Å of each Cu
  - Report the closest atom per residue, its distance, and the residue number

Output: cu_coordination_3A.tsv — one row per Cu-residue contact.
Columns: accession, cu_id (Cu1/Cu2), resname, resnum, closest_atom, distance, n_coord_his
"""

import argparse
import csv
import os
import sys
import numpy as np

COORD_CUTOFF = 3.0


def get_cu_atoms(cif_path):
    cus = []
    cols = []
    in_atom_site = False
    in_data = False
    x_col = y_col = z_col = elem_col = bfac_col = None

    with open(cif_path) as f:
        for line in f:
            if line.startswith('_atom_site.'):
                if not in_atom_site:
                    in_atom_site = True
                    cols = []
                cols.append(line.strip().split('.')[1].strip())
                continue
            if in_atom_site and not line.startswith('_atom_site.'):
                in_atom_site = False
                in_data = True
                x_col = cols.index('Cartn_x')
                y_col = cols.index('Cartn_y')
                z_col = cols.index('Cartn_z')
                elem_col = cols.index('type_symbol')
                bfac_col = cols.index('B_iso_or_equiv') if 'B_iso_or_equiv' in cols else None
            if in_data:
                if line.startswith('#') or line.startswith('loop_') or line.startswith('_'):
                    break
                parts = line.split()
                if len(parts) > max(x_col, y_col, z_col, elem_col):
                    if parts[elem_col] == 'CU':
                        cus.append({
                            'coord': np.array([float(parts[x_col]), float(parts[y_col]), float(parts[z_col])]),
                            'bfactor': float(parts[bfac_col]) if bfac_col is not None else 0,
                        })
    return cus


def get_contacts(cif_path, cu_coord, cutoff):
    cols = []
    in_atom_site = False
    in_data = False
    contacts_by_res = {}

    with open(cif_path) as f:
        for line in f:
            if line.startswith('_atom_site.'):
                if not in_atom_site:
                    in_atom_site = True
                    cols = []
                cols.append(line.strip().split('.')[1].strip())
                continue
            if in_atom_site and not line.startswith('_atom_site.'):
                in_atom_site = False
                in_data = True
                x_col = cols.index('Cartn_x')
                y_col = cols.index('Cartn_y')
                z_col = cols.index('Cartn_z')
                elem_col = cols.index('type_symbol')
                resname_col = cols.index('label_comp_id')
                resnum_col = cols.index('label_seq_id')
                atomname_col = cols.index('label_atom_id')
                group_col = cols.index('group_PDB')
            if in_data:
                if line.startswith('#') or line.startswith('loop_') or line.startswith('_'):
                    break
                parts = line.split()
                if len(parts) <= max(x_col, y_col, z_col):
                    continue
                if parts[group_col] != 'ATOM':
                    continue
                if parts[elem_col] == 'H':
                    continue
                coord = np.array([float(parts[x_col]), float(parts[y_col]), float(parts[z_col])])
                d = float(np.linalg.norm(coord - cu_coord))
                resname = parts[resname_col]
                resnum = parts[resnum_col]
                atomname = parts[atomname_col]
                key = (resname, resnum)
                if key not in contacts_by_res or d < contacts_by_res[key]['distance']:
                    contacts_by_res[key] = {
                        'resname': resname,
                        'resnum': resnum,
                        'closest_atom': atomname,
                        'distance': round(d, 3),
                    }

    return [c for c in contacts_by_res.values() if c['distance'] <= cutoff]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cif-dir', required=True)
    ap.add_argument('--pool', required=True, help='nc_pool.csv')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    pool = []
    with open(args.pool) as f:
        for r in csv.DictReader(f):
            pool.append(r['accession'])
    print(f"Pool: {len(pool)} accessions", flush=True)

    cif_files = {}
    for fn in os.listdir(args.cif_dir):
        if fn.endswith('.cif'):
            acc = fn.split('_taxID_')[0]
            cif_files[acc] = fn

    out = open(args.out, 'w', newline='')
    w = csv.writer(out, delimiter='\t')
    w.writerow(['accession', 'cu_id', 'resname', 'resnum',
                'closest_atom', 'distance', 'n_coord_his'])

    n_ok = n_skip = 0
    for i, acc in enumerate(pool):
        if acc not in cif_files:
            n_skip += 1
            continue

        cif_path = os.path.join(args.cif_dir, cif_files[acc])
        try:
            cus = get_cu_atoms(cif_path)
            if len(cus) < 2:
                n_skip += 1
                continue

            all_contacts = []
            for ci, cu in enumerate(cus[:2]):
                contacts = get_contacts(cif_path, cu['coord'], COORD_CUTOFF)
                for c in contacts:
                    c['cu_id'] = f'Cu{ci+1}'
                all_contacts.extend(contacts)

            n_his = sum(1 for c in all_contacts
                        if c['resname'] == 'HIS'
                        and c['closest_atom'] in ('NE2', 'ND1')
                        and c['distance'] <= 3.5)

            for c in all_contacts:
                w.writerow([acc, c['cu_id'], c['resname'], c['resnum'],
                            c['closest_atom'], c['distance'], n_his])

            n_ok += 1
        except Exception as e:
            print(f"  {acc}: {e}", file=sys.stderr, flush=True)
            n_skip += 1

        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{len(pool)} processed ({n_ok} ok, {n_skip} skipped)", flush=True)

    out.close()
    print(f"\nDone: {n_ok} ok, {n_skip} skipped", flush=True)


if __name__ == '__main__':
    main()
