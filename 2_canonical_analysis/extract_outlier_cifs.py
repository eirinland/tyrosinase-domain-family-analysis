"""
Extract active-site residues from vector group members into mini-CIF files
for visual inspection in PyMOL. Outputs one directory per vector group,
containing one CIF per structure with only the mapped residues + Cu atoms.
"""

import csv
import argparse
from pathlib import Path
from multiprocessing import Pool

import numpy as np

from extract_position_vectors import (
    parse_atoms, get_cu, get_coord_his_ne2, get_ca_map, get_aa_map,
    kabsch, load_reference, PMTYR_POSITIONS, VARIABLE_PMTYR,
)

_pmtyr_ref = None


def init_worker(pmtyr_path):
    global _pmtyr_ref
    _pmtyr_ref = load_reference(pmtyr_path)


def get_mapped_residues(atoms, ref):
    """Return set of query residue seq numbers that map to reference positions."""
    cu1, cu2 = get_cu(atoms)
    if cu1 is None:
        return None, None, None

    coord_his = get_coord_his_ne2(atoms, cu1, cu2)
    if len(coord_his) < 6:
        return None, None, None

    qry_ca = get_ca_map(atoms)
    qry_ne2_coords = np.array([[a['x'], a['y'], a['z']] for a in coord_his[:6]])
    R, t = kabsch(qry_ne2_coords, ref['ne2_coords'])

    qry_ca_t = {seq: R @ coord + t for seq, coord in qry_ca.items()}

    mapped_seqs = set()
    for ref_pos in PMTYR_POSITIONS:
        if ref_pos not in ref['ca']:
            continue
        ref_coord = ref['ca'][ref_pos]
        best_dist, best_seq = 999.0, None
        for qseq, qcoord in qry_ca_t.items():
            d = float(np.linalg.norm(qcoord - ref_coord))
            if d < best_dist:
                best_dist = d
                best_seq = qseq
        if best_seq is not None:
            mapped_seqs.add(best_seq)

    return mapped_seqs, R, t


def extract_active_site_cif(cif_path, output_path):
    """Read CIF, keep only residues mapped to reference positions + Cu,
    apply Kabsch transformation so output is in PmTYR reference frame."""
    with open(cif_path) as f:
        lines = f.readlines()

    atoms_raw = parse_atoms(cif_path)
    mapped_seqs, R, t = get_mapped_residues(atoms_raw, _pmtyr_ref)
    if mapped_seqs is None:
        return False

    # Find coordinate column indices and transform
    col_names = []
    header_lines = []
    atom_lines = []
    in_atom = False
    collecting = False
    past_atom = False

    x_idx = y_idx = z_idx = None

    for raw in lines:
        line = raw.strip()
        if past_atom:
            break
        if line == 'loop_' and not in_atom:
            header_lines.append(raw)
            collecting = False
            col_names = []
            continue
        if line.startswith('_atom_site.'):
            collecting = True
            col_names.append(line)
            header_lines.append(raw)
            continue
        if collecting and col_names and not line.startswith('_atom_site.'):
            collecting = False
            in_atom = True
            x_idx = col_names.index('_atom_site.Cartn_x')
            y_idx = col_names.index('_atom_site.Cartn_y')
            z_idx = col_names.index('_atom_site.Cartn_z')

        if in_atom:
            if line.startswith('_') or line == '#' or not line:
                past_atom = True
                break
            parts = line.split()
            if len(parts) != len(col_names):
                continue
            row = dict(zip(col_names, parts))
            seq = row.get('_atom_site.label_seq_id', '')
            elem = row.get('_atom_site.type_symbol', '').upper()
            if elem == 'CU' or (seq.isdigit() and int(seq) in mapped_seqs):
                # Transform coordinates into PmTYR reference frame
                xyz = np.array([float(parts[x_idx]), float(parts[y_idx]), float(parts[z_idx])])
                xyz_t = R @ xyz + t
                parts[x_idx] = f'{xyz_t[0]:.3f}'
                parts[y_idx] = f'{xyz_t[1]:.3f}'
                parts[z_idx] = f'{xyz_t[2]:.3f}'
                atom_lines.append(' '.join(parts) + '\n')
        else:
            header_lines.append(raw)

    with open(output_path, 'w') as f:
        for line in header_lines:
            f.write(line)
        for line in atom_lines:
            f.write(line)
        f.write('#\n')

    return True


def process_one(args):
    cif_path, output_path = args
    try:
        return extract_active_site_cif(cif_path, output_path)
    except Exception:
        return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cifs', required=True)
    parser.add_argument('--pmtyr', required=True)
    parser.add_argument('--outliers-csv', required=True)
    parser.add_argument('--vectors', nargs='+', required=True,
                        help='Vector strings to extract (quote each)')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--workers', type=int, default=16)
    args = parser.parse_args()

    target_vectors = set(args.vectors)

    # Read outlier data to identify members and outlier status
    members = {}  # vector -> [(accession, is_outlier, rmsd)]
    with open(args.outliers_csv) as f:
        for row in csv.DictReader(f):
            vec = row['vector']
            if vec in target_vectors:
                if vec not in members:
                    members[vec] = []
                members[vec].append((
                    row['accession'],
                    row['is_outlier'] == 'True',
                    float(row['rmsd_to_centroid']),
                ))

    cif_dir = Path(args.cifs)
    out_base = Path(args.output_dir)
    out_base.mkdir(parents=True, exist_ok=True)

    tasks = []
    for vec, mems in members.items():
        safe_name = vec.replace('*', 'x').replace('-', '_')
        vec_dir = out_base / safe_name
        vec_dir.mkdir(exist_ok=True)

        for acc, is_outlier, rmsd in mems:
            cif_matches = list(cif_dir.glob(f'{acc}_taxID_*_model.cif'))
            if not cif_matches:
                continue
            prefix = 'OUTLIER_' if is_outlier else ''
            out_name = f'{safe_name}__{prefix}{acc}_rmsd{rmsd:.3f}.cif'
            tasks.append((str(cif_matches[0]), str(vec_dir / out_name)))

    print(f"Vectors requested: {len(target_vectors)}", flush=True)
    print(f"Structures to extract: {len(tasks)}", flush=True)

    done = 0
    success = 0
    with Pool(processes=args.workers,
              initializer=init_worker,
              initargs=(args.pmtyr,)) as pool:
        for result in pool.imap_unordered(process_one, tasks, chunksize=10):
            done += 1
            if result:
                success += 1
            if done % 200 == 0:
                print(f"  {done}/{len(tasks)}...", flush=True)

    print(f"\nDone. Extracted {success}/{len(tasks)} structures.", flush=True)
    print(f"Output: {args.output_dir}")

    for vec in target_vectors:
        safe_name = vec.replace('*', 'x').replace('-', '_')
        vec_dir = out_base / safe_name
        n_files = len(list(vec_dir.glob('*.cif')))
        n_outliers = len(list(vec_dir.glob('OUTLIER_*.cif')))
        print(f"  {vec}: {n_files} CIFs ({n_outliers} outliers)")


if __name__ == '__main__':
    main()
