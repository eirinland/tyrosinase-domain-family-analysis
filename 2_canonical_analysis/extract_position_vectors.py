"""
Active-site position vector extraction with geometry-based thioether detection.

Superpose each structure onto PmTYR (B2ZB02) using Kabsch SVD on 6 coordinating
His NE2 atoms. For each reference position, find the nearest query Cα after
superposition and record the residue identity.

Thioether Cys is detected geometrically (no reference needed):
  - Primary (3D):  Cys SG within 3.5 Å of CuA His2 imidazole ring
  - Fallback (seq): Cys between CuA His1 and CuA His2 in sequence
  CuA = first three coordinating His by sequence number.

Output: CSV with one row per structure, columns for each reference position
        plus RMSD, thioether detection, and a combined vector string.
"""

import math
import csv
import argparse
from pathlib import Path
from multiprocessing import Pool
from collections import Counter

import numpy as np
import biotite.structure.io.pdbx as bpdbx
import biotite.structure as bstruc

AA3 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
}

IMIDAZOLE_ATOMS = {'ND1', 'NE2', 'CD2', 'CE1', 'CG'}

# PmTYR reference positions (B2ZB02 AF3 numbering)
# His 42,60,69 = CuA (first sequential triplet); His 204,208,231 = CuB
PMTYR_POSITIONS = [42, 46, 60, 65, 68, 69, 195, 204, 205, 208, 209, 218, 221, 227, 230, 231]
PMTYR_ANCHORS = {42, 60, 69, 204, 208, 231}

VARIABLE_PMTYR = [p for p in PMTYR_POSITIONS if p not in PMTYR_ANCHORS]

PMTYR_LABELS = {
    42: 'His42', 46: 'Gly46',
    60: 'His60', 65: 'Phe65', 68: 'Trp68',
    69: 'His69', 195: 'Glu195',
    204: 'His204', 205: 'Asn205',
    208: 'His208', 209: 'Arg209',
    218: 'Val218', 221: 'Ala221', 227: 'Phe227', 230: 'His230',
    231: 'His231',
}

# Reference secondary structure from PmTYR (biotite: a=helix, b=sheet, c=coil)
PMTYR_SS = {
    42: 'a', 46: 'a',
    60: 'c', 65: 'a', 68: 'a',
    69: 'a', 195: 'a',
    204: 'a', 205: 'a',
    208: 'a', 209: 'a',
    218: 'c', 221: 'c', 227: 'a', 230: 'a',
    231: 'a',
}


def parse_atoms(path):
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
                    })
                except (KeyError, ValueError):
                    continue
    return atoms


def dist3d(a, b):
    return math.sqrt((a['x']-b['x'])**2 + (a['y']-b['y'])**2 + (a['z']-b['z'])**2)


def get_cu(atoms):
    cu = [a for a in atoms if a['elem'] == 'CU']
    cu.sort(key=lambda a: a['chain'])
    if len(cu) < 2:
        return None, None
    return cu[0], cu[1]


def get_coord_his_ne2(atoms, cu1, cu2, cutoff=3.0):
    """Get 6 coordinating His NE2 atoms sorted by sequence number."""
    his_ne2 = [a for a in atoms if a['resn'] == 'HIS' and a['atom'] == 'NE2']
    coord = []
    for a in his_ne2:
        if min(dist3d(a, cu1), dist3d(a, cu2)) <= cutoff:
            coord.append(a)
    coord.sort(key=lambda a: int(a['seq']))
    return coord


def get_ca_map(atoms):
    ca = {}
    for a in atoms:
        if a['atom'] == 'CA' and a['seq'].isdigit():
            ca[int(a['seq'])] = np.array([a['x'], a['y'], a['z']])
    return ca


def get_aa_map(atoms):
    aa = {}
    for a in atoms:
        if a['seq'].isdigit() and a['resn'] in AA3:
            seq = int(a['seq'])
            if seq not in aa:
                aa[seq] = AA3[a['resn']]
    return aa


def get_ss_map(path):
    """Per-residue secondary structure via biotite (a=helix, b=sheet, c=coil)."""
    cif = bpdbx.CIFFile.read(path)
    atoms = bpdbx.get_structure(cif, model=1)
    protein = atoms[bstruc.filter_amino_acids(atoms)]
    sse = bstruc.annotate_sse(protein)
    ca = protein[protein.atom_name == 'CA']
    return {int(rid): ss for rid, ss in zip(ca.res_id, sse)}


def kabsch(P, Q):
    cP, cQ = P.mean(0), Q.mean(0)
    H = (P - cP).T @ (Q - cQ)
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = cQ - R @ cP
    return R, t


def load_reference(path):
    atoms = parse_atoms(path)
    cu1, cu2 = get_cu(atoms)
    coord_his = get_coord_his_ne2(atoms, cu1, cu2)
    return {
        'ne2_order': coord_his,
        'ne2_coords': np.array([[a['x'], a['y'], a['z']] for a in coord_his]),
        'ca': get_ca_map(atoms),
        'aa': get_aa_map(atoms),
    }


def superpose_and_map(ref, qry_atoms, positions):
    """Superpose query onto reference and map positions. Returns (rmsd, results_dict).

    Pairing: sort both reference and query coordinating His NE2 by sequence
    number, then pair 1:1 (His1↔His1, His2↔His2, ... His6↔His6).
    """
    cu1, cu2 = get_cu(qry_atoms)
    if cu1 is None:
        return None, {p: ('?', None, 999) for p in positions}

    coord_his = get_coord_his_ne2(qry_atoms, cu1, cu2)
    if len(coord_his) < 6:
        return None, {p: ('?', None, 999) for p in positions}

    qry_ca = get_ca_map(qry_atoms)
    qry_aa = get_aa_map(qry_atoms)

    # Pair by sequence order: ref ne2_order[i] ↔ coord_his[i]
    qry_ne2_coords = np.array([[a['x'], a['y'], a['z']] for a in coord_his[:6]])
    R, t = kabsch(qry_ne2_coords, ref['ne2_coords'])

    ta = (R @ qry_ne2_coords.T).T + t
    rmsd = float(np.sqrt(np.mean(np.sum((ta - ref['ne2_coords'])**2, axis=1))))

    qry_ca_t = {seq: R @ coord + t for seq, coord in qry_ca.items()}

    results = {}
    for ref_pos in positions:
        if ref_pos not in ref['ca']:
            results[ref_pos] = ('?', None, 999)
            continue
        ref_coord = ref['ca'][ref_pos]
        best_dist, best_seq = 999.0, None
        for qseq, qcoord in qry_ca_t.items():
            d = float(np.linalg.norm(qcoord - ref_coord))
            if d < best_dist:
                best_dist = d
                best_seq = qseq
        results[ref_pos] = (qry_aa.get(best_seq, '?'), best_seq, best_dist)

    return rmsd, results


def detect_thioether(atoms, sg_cutoff=3.5):
    """Detect thioether Cys by geometry: Cys SG near CuA His2 imidazole ring.

    CuA = first 3 coordinating His by sequence number.
    Returns (label, sg_dist):
      label:   'C'  if SG within sg_cutoff of CuA His2 ring (3D positive)
               'C*' if Cys between CuA His1-His2 in sequence only (seq positive)
               '-'  if no thioether detected
      sg_dist: distance from nearest Cys SG to CuA His2 ring (float or None)
    """
    cu1, cu2 = get_cu(atoms)
    if cu1 is None:
        return '-', None

    his_ne2 = [a for a in atoms if a['resn'] == 'HIS' and a['atom'] == 'NE2']
    coord_his = []
    for a in his_ne2:
        d1, d2 = dist3d(a, cu1), dist3d(a, cu2)
        if min(d1, d2) <= 3.0:
            coord_his.append(a)
    coord_his.sort(key=lambda a: int(a['seq']))
    if len(coord_his) < 6:
        return '-', None

    cua_his1_seq = int(coord_his[0]['seq'])
    cua_his2_seq = int(coord_his[1]['seq'])

    # 3D method: Cys SG within cutoff of CuA His2 imidazole ring
    his2_ring = [a for a in atoms if a['resn'] == 'HIS'
                 and int(a['seq']) == cua_his2_seq and a['atom'] in IMIDAZOLE_ATOMS]
    cys_sg = [a for a in atoms if a['resn'] == 'CYS' and a['atom'] == 'SG']

    best_sg_dist = 999.0
    for sg in cys_sg:
        for ha in his2_ring:
            d = dist3d(sg, ha)
            if d < best_sg_dist:
                best_sg_dist = d

    if best_sg_dist <= sg_cutoff:
        return 'C', best_sg_dist

    # Sequence fallback: Cys between CuA His1 and CuA His2
    aa_map = get_aa_map(atoms)
    for seq in range(cua_his1_seq + 1, cua_his2_seq):
        if aa_map.get(seq) == 'C':
            return 'C*', best_sg_dist if best_sg_dist < 900 else None

    return '-', best_sg_dist if best_sg_dist < 900 else None


# Global reference (set in worker init)
_pmtyr_ref = None


def init_worker(pmtyr_path):
    global _pmtyr_ref
    _pmtyr_ref = load_reference(pmtyr_path)


def process_one(cif_path):
    accession = Path(cif_path).name.split('_taxID_')[0]
    try:
        qry_atoms = parse_atoms(cif_path)
    except Exception as e:
        return {'accession': accession, 'error': str(e)}

    rmsd_pm, res_pm = superpose_and_map(_pmtyr_ref, qry_atoms, PMTYR_POSITIONS)
    thio_label, thio_sg_dist = detect_thioether(qry_atoms)

    try:
        ss_map = get_ss_map(cif_path)
    except Exception:
        ss_map = {}

    row = {'accession': accession}
    row['pmtyr_rmsd'] = f'{rmsd_pm:.3f}' if rmsd_pm is not None else 'NA'

    for pos in PMTYR_POSITIONS:
        aa, qseq, ca_d = res_pm[pos]
        label = PMTYR_LABELS[pos]
        row[label] = aa
        row[f'{label}_cadist'] = f'{ca_d:.2f}' if ca_d < 900 else 'NA'
        qry_ss = ss_map.get(qseq) if qseq is not None else None
        ref_ss = PMTYR_SS.get(pos)
        if qry_ss is not None and ref_ss is not None:
            row[f'{label}_ss'] = qry_ss
            row[f'{label}_ss_match'] = str(qry_ss == ref_ss)
        else:
            row[f'{label}_ss'] = 'NA'
            row[f'{label}_ss_match'] = 'NA'

    row['thioether'] = thio_label
    row['thioether_sg_dist'] = f'{thio_sg_dist:.2f}' if thio_sg_dist is not None else 'NA'

    # Build vector string from variable positions + thioether
    var_parts = []
    for pos in VARIABLE_PMTYR:
        aa, qseq, ca_d = res_pm[pos]
        if pos == 46 and ss_map.get(qseq) == 'c':
            var_parts.append('~')
        else:
            var_parts.append(aa if ca_d < 3.0 else f'{aa}*')
    var_parts.append(thio_label)
    row['vector'] = '-'.join(var_parts)

    row['error'] = ''
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cifs', required=True, help='CIF directory (squashFS mount)')
    parser.add_argument('--passed', required=True, help='after_his_coordination_filter.csv')
    parser.add_argument('--pmtyr', required=True, help='PmTYR CIF (B2ZB02)')
    parser.add_argument('--output', required=True, help='Output CSV')
    parser.add_argument('--workers', type=int, default=16)
    args = parser.parse_args()

    passed_acc = set()
    with open(args.passed) as f:
        for row in csv.DictReader(f):
            if row.get('passes_filter', 'True') != 'True':
                continue
            sid = row.get('accession', row.get('sequence_id', ''))
            bare = sid.split('_taxID_')[0] if '_taxID_' in sid else sid
            passed_acc.add(bare)
    print(f"Accessions: {len(passed_acc)}", flush=True)

    cif_dir = Path(args.cifs)
    cif_files = [
        str(p) for p in sorted(cif_dir.glob('*.cif'))
        if p.name.split('_taxID_')[0] in passed_acc
    ]
    print(f"CIF files matched: {len(cif_files)}", flush=True)

    fieldnames = ['accession', 'pmtyr_rmsd']
    for pos in PMTYR_POSITIONS:
        label = PMTYR_LABELS[pos]
        fieldnames.extend([label, f'{label}_cadist', f'{label}_ss', f'{label}_ss_match'])
    fieldnames.extend(['thioether', 'thioether_sg_dist', 'vector', 'error'])

    results = []
    done = 0
    with Pool(processes=args.workers,
              initializer=init_worker,
              initargs=(args.pmtyr,)) as pool:
        for row in pool.imap_unordered(process_one, cif_files, chunksize=50):
            results.append(row)
            done += 1
            if done % 2000 == 0:
                print(f"  {done}/{len(cif_files)}...", flush=True)

    results.sort(key=lambda r: r['accession'])

    with open(args.output, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        w.writeheader()
        w.writerows(results)

    errors = sum(1 for r in results if r.get('error'))
    print(f"\nDone. {len(results)} structures, {errors} errors.", flush=True)
    print(f"Output: {args.output}", flush=True)

    # Quick summary
    vectors = Counter(r.get('vector', '') for r in results if not r.get('error'))
    print(f"\nUnique vectors: {len(vectors)}")
    print("Top 10:")
    for v, n in vectors.most_common(10):
        print(f"  {n:>6}  {v}")

    # Thioether summary
    thio = Counter(r.get('thioether', '') for r in results if not r.get('error'))
    print(f"\nThioether detection:")
    for label, n in thio.most_common():
        print(f"  {label:>3}: {n}")


if __name__ == '__main__':
    main()
