#!/usr/bin/env python3
"""Extract Cu/His active-site geometry from all AF3 seed/sample CIFs.
Stdlib only (no gemmi) -- mirrors check_canonical_criteria.py parser."""

import csv
import math
import re
import sys
from pathlib import Path


def parse_atoms(path):
    atoms = []
    col_names = []
    in_atom = False
    collecting = False
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if line == "loop_":
                collecting = False; in_atom = False; col_names = []
                continue
            if line.startswith("_atom_site."):
                collecting = True; col_names.append(line)
                continue
            if collecting and col_names:
                collecting = False; in_atom = True
            if in_atom:
                if line.startswith("_") or line == "#" or not line:
                    break
                parts = line.split()
                if len(parts) != len(col_names):
                    continue
                row = dict(zip(col_names, parts))
                try:
                    atoms.append({
                        "elem": row.get("_atom_site.type_symbol", "").upper(),
                        "atom": row.get("_atom_site.label_atom_id", "").upper(),
                        "resn": row.get("_atom_site.label_comp_id", ""),
                        "seq": row.get("_atom_site.label_seq_id", ""),
                        "x": float(row["_atom_site.Cartn_x"]),
                        "y": float(row["_atom_site.Cartn_y"]),
                        "z": float(row["_atom_site.Cartn_z"]),
                        "bfactor": float(row.get("_atom_site.B_iso_or_equiv", "0")),
                    })
                except (KeyError, ValueError):
                    continue
    return atoms


def dist(a, b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)


def parse_cif_name(path):
    name = re.sub(r'_model\.cif$', '', path.name)
    m = re.match(r'^([A-Za-z0-9]+)_taxID_\d+', name)
    if not m:
        return None, None, None
    accession = m.group(1)
    seed_m = re.search(r'_seed-(\d+)', name)
    sample_m = re.search(r'_sample-(\d+)', name)
    seed = int(seed_m.group(1)) if seed_m else 0
    sample = int(sample_m.group(1)) if sample_m else -1
    return accession, seed, sample


def extract_geometry(cif_path):
    try:
        atoms = parse_atoms(str(cif_path))
    except Exception as e:
        return {'error': str(e)}

    cu = [a for a in atoms if a["elem"] == "CU"]
    cu.sort(key=lambda a: a.get("seq", ""))

    n_cu = len(cu)
    cu_cu_dist = None
    cu1_plddt = cu[0]["bfactor"] if n_cu >= 1 else None
    cu2_plddt = cu[1]["bfactor"] if n_cu >= 2 else None
    if n_cu >= 2:
        cu_cu_dist = dist(cu[0], cu[1])

    his_ne2 = [a for a in atoms if a["resn"] == "HIS" and a["atom"] == "NE2"]
    his_ca = {a["seq"]: a for a in atoms if a["resn"] == "HIS" and a["atom"] == "CA"}

    seen_seq = set()
    n_coordinating = 0
    coord_his_ca_plddts = []
    for h in his_ne2:
        if h["seq"] in seen_seq:
            continue
        if cu:
            min_cu_dist = min(dist(h, c) for c in cu)
            if min_cu_dist <= 3.5:
                seen_seq.add(h["seq"])
                n_coordinating += 1
                ca = his_ca.get(h["seq"])
                if ca:
                    coord_his_ca_plddts.append(ca["bfactor"])

    n_his_total = len({a["seq"] for a in atoms if a["resn"] == "HIS" and a["atom"] == "CA"})
    min_coord_his_ca_plddt = min(coord_his_ca_plddts) if coord_his_ca_plddts else None

    canonical = (
        n_cu >= 2
        and cu_cu_dist is not None
        and 2.8 <= cu_cu_dist <= 5.5
        and n_coordinating >= 6
        and min_coord_his_ca_plddt is not None
        and min_coord_his_ca_plddt >= 70.0
    )

    return {
        'n_cu': n_cu,
        'cu_cu_dist': cu_cu_dist,
        'cu1_plddt': cu1_plddt,
        'cu2_plddt': cu2_plddt,
        'n_his_total': n_his_total,
        'n_coordinating_his': n_coordinating,
        'min_coord_his_ca_plddt': min_coord_his_ca_plddt,
        'canonical': canonical,
        'error': None,
    }


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <mount_dir> <output.tsv>", file=sys.stderr)
        sys.exit(1)

    mount_dir = Path(sys.argv[1])
    output_path = Path(sys.argv[2])

    cifs = sorted(mount_dir.rglob('*_model.cif'))
    print(f"Found {len(cifs)} CIF files in {mount_dir}", file=sys.stderr)

    fieldnames = [
        'accession', 'seed', 'sample', 'n_cu', 'cu_cu_dist',
        'cu1_plddt', 'cu2_plddt', 'n_his_total', 'n_coordinating_his',
        'min_coord_his_ca_plddt', 'canonical', 'error',
    ]

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()

        for i, cif in enumerate(cifs):
            if (i + 1) % 1000 == 0:
                print(f"  {i+1}/{len(cifs)}", file=sys.stderr)

            accession, seed, sample = parse_cif_name(cif)
            if accession is None:
                continue

            geom = extract_geometry(cif)
            row = {
                'accession': accession,
                'seed': seed,
                'sample': 'top' if sample < 0 else sample,
            }
            for k, v in geom.items():
                row[k] = f'{v:.4f}' if isinstance(v, float) else v
            writer.writerow(row)

    print(f"Wrote {output_path}", file=sys.stderr)


if __name__ == '__main__':
    main()
