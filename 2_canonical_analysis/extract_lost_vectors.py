"""
Extract position vectors for the 229 'lost canonical' structures using
multi-reference alignment.

These structures pass the alignment-free Cu/His/pLDDT check but fail the
pipeline's alignment-based filters, so PmTYR alone may give poor superpositions.
We try 5 diverse references and keep the best RMSD per structure.

References:
  B2ZB02  - PmTYR (fungal TYR)
  P14679  - HsTYR (animal TYR)
  Q83WS2  - BtTYR (bacterial TYR)
  G2QLD3  - AboMP (fungal oMP)
  Q9ZP19  - LcCaOx (plant CaOx)

Usage:
  python extract_lost_vectors.py --cifs /mnt/models \
      --lost alignmentfree_test.tsv --output lost_vectors.csv
"""

import csv, math, argparse
from pathlib import Path
from multiprocessing import Pool
from collections import Counter

import numpy as np
import biotite.structure.io.pdbx as bpdbx
import biotite.structure as bstruc

from extract_position_vectors import (
    parse_atoms, dist3d, get_cu, get_coord_his_ne2, get_ca_map, get_aa_map,
    get_ss_map, kabsch, detect_thioether, AA3,
    PMTYR_POSITIONS, PMTYR_ANCHORS, PMTYR_LABELS, PMTYR_SS,
    VARIABLE_PMTYR, IMIDAZOLE_ATOMS,
)

REFERENCES = ['B2ZB02', 'P14679', 'Q83WS2', 'G2QLD3', 'Q9ZP19']


def load_pmtyr_ref(cif_path):
    atoms = parse_atoms(cif_path)
    cu1, cu2 = get_cu(atoms)
    coord_his = get_coord_his_ne2(atoms, cu1, cu2)
    return {
        'name': 'B2ZB02',
        'ne2_coords': np.array([[a['x'], a['y'], a['z']] for a in coord_his[:6]]),
        'ca': get_ca_map(atoms),
        'aa': get_aa_map(atoms),
        'pos_map': {p: p for p in PMTYR_POSITIONS},
    }


def map_reference(pmtyr_ref, ref_path, ref_name):
    """Superpose a non-PmTYR reference onto PmTYR, find equivalent positions."""
    atoms = parse_atoms(ref_path)
    cu1, cu2 = get_cu(atoms)
    if cu1 is None:
        return None
    coord_his = get_coord_his_ne2(atoms, cu1, cu2)
    if len(coord_his) < 6:
        return None

    ref_ne2 = np.array([[a['x'], a['y'], a['z']] for a in coord_his[:6]])
    R, t = kabsch(ref_ne2, pmtyr_ref['ne2_coords'])

    ref_ca = get_ca_map(atoms)
    ref_ca_t = {seq: R @ coord + t for seq, coord in ref_ca.items()}

    pos_map = {}
    for pmtyr_pos in PMTYR_POSITIONS:
        target = pmtyr_ref['ca'][pmtyr_pos]
        best_d, best_seq = 999.0, None
        for seq, coord in ref_ca_t.items():
            d = float(np.linalg.norm(coord - target))
            if d < best_d:
                best_d = d
                best_seq = seq
        if best_d < 5.0:
            pos_map[pmtyr_pos] = best_seq
        else:
            pos_map[pmtyr_pos] = None

    rmsd = float(np.sqrt(np.mean(np.sum(
        ((R @ ref_ne2.T).T + t - pmtyr_ref['ne2_coords'])**2, axis=1))))

    return {
        'name': ref_name,
        'ne2_coords': ref_ne2,
        'ca': ref_ca,
        'aa': get_aa_map(atoms),
        'pos_map': pos_map,
        'rmsd_to_pmtyr': rmsd,
    }


def superpose_query(ref, qry_atoms):
    """Superpose query onto reference, extract mapped positions. Returns (rmsd, results)."""
    cu1, cu2 = get_cu(qry_atoms)
    if cu1 is None:
        return None, {}

    coord_his = get_coord_his_ne2(qry_atoms, cu1, cu2)
    if len(coord_his) < 6:
        return None, {}

    qry_ne2 = np.array([[a['x'], a['y'], a['z']] for a in coord_his[:6]])
    R, t = kabsch(qry_ne2, ref['ne2_coords'])

    rmsd = float(np.sqrt(np.mean(np.sum(
        ((R @ qry_ne2.T).T + t - ref['ne2_coords'])**2, axis=1))))

    qry_ca = get_ca_map(qry_atoms)
    qry_aa = get_aa_map(qry_atoms)
    qry_ca_t = {seq: R @ coord + t for seq, coord in qry_ca.items()}

    results = {}
    for pmtyr_pos in PMTYR_POSITIONS:
        ref_pos = ref['pos_map'].get(pmtyr_pos)
        if ref_pos is None:
            results[pmtyr_pos] = ('?', None, 999)
            continue
        ref_coord = ref['ca'][ref_pos]
        best_d, best_seq = 999.0, None
        for qseq, qcoord in qry_ca_t.items():
            d = float(np.linalg.norm(qcoord - ref_coord))
            if d < best_d:
                best_d = d
                best_seq = qseq
        results[pmtyr_pos] = (qry_aa.get(best_seq, '?'), best_seq, best_d)

    return rmsd, results


_refs = None


def init_worker(refs):
    global _refs
    _refs = refs


def process_one(cif_path):
    accession = Path(cif_path).name.split('_taxID_')[0]
    try:
        qry_atoms = parse_atoms(cif_path)
    except Exception as e:
        return {'accession': accession, 'error': str(e)}

    best_rmsd, best_res, best_ref = 999.0, None, None
    for ref in _refs:
        rmsd, res = superpose_query(ref, qry_atoms)
        if rmsd is not None and rmsd < best_rmsd:
            best_rmsd = rmsd
            best_res = res
            best_ref = ref['name']

    if best_res is None:
        row = {'accession': accession, 'error': 'no_alignment'}
        for pos in PMTYR_POSITIONS:
            label = PMTYR_LABELS[pos]
            row[label] = '?'
            row[f'{label}_cadist'] = 'NA'
            row[f'{label}_ss'] = 'NA'
            row[f'{label}_ss_match'] = 'NA'
        row.update({'thioether': '-', 'thioether_sg_dist': 'NA',
                    'vector': '', 'best_ref': '', 'best_rmsd': 'NA'})
        return row

    thio_label, thio_sg_dist = detect_thioether(qry_atoms)
    try:
        ss_map = get_ss_map(cif_path)
    except Exception:
        ss_map = {}

    row = {'accession': accession, 'best_ref': best_ref,
           'best_rmsd': f'{best_rmsd:.3f}'}

    for pos in PMTYR_POSITIONS:
        aa, qseq, ca_d = best_res[pos]
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
    row['thioether_sg_dist'] = f'{thio_sg_dist:.2f}' if thio_sg_dist else 'NA'

    var_parts = []
    for pos in VARIABLE_PMTYR:
        aa, qseq, ca_d = best_res[pos]
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
    parser.add_argument('--cifs', required=True)
    parser.add_argument('--lost', required=True, help='alignmentfree_test.tsv')
    parser.add_argument('--output', required=True)
    parser.add_argument('--workers', type=int, default=8)
    args = parser.parse_args()

    lost_acc = set()
    with open(args.lost) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row.get('passes_all') == 'True':
                lost_acc.add(row['accession'])
    print(f"Lost canonical accessions: {len(lost_acc)}", flush=True)

    cif_dir = Path(args.cifs)
    cif_map = {}
    for p in cif_dir.glob('*.cif'):
        acc = p.name.split('_taxID_')[0]
        if acc in lost_acc:
            cif_map[acc] = str(p)
    print(f"CIF files matched: {len(cif_map)}", flush=True)

    # Load references
    ref_cifs = {}
    for acc in REFERENCES:
        matches = list(cif_dir.glob(f'{acc}_taxID_*.cif'))
        if matches:
            ref_cifs[acc] = str(matches[0])
        else:
            print(f"  WARNING: reference {acc} not found in CIF dir")
    print(f"References found: {list(ref_cifs.keys())}", flush=True)

    pmtyr_ref = load_pmtyr_ref(ref_cifs['B2ZB02'])
    refs = [pmtyr_ref]
    for acc in REFERENCES[1:]:
        if acc in ref_cifs:
            mapped = map_reference(pmtyr_ref, ref_cifs[acc], acc)
            if mapped:
                n_mapped = sum(1 for v in mapped['pos_map'].values() if v is not None)
                print(f"  {acc}: RMSD to PmTYR = {mapped['rmsd_to_pmtyr']:.3f}, "
                      f"positions mapped = {n_mapped}/{len(PMTYR_POSITIONS)}")
                refs.append(mapped)

    fieldnames = ['accession', 'best_ref', 'best_rmsd']
    for pos in PMTYR_POSITIONS:
        label = PMTYR_LABELS[pos]
        fieldnames.extend([label, f'{label}_cadist', f'{label}_ss', f'{label}_ss_match'])
    fieldnames.extend(['thioether', 'thioether_sg_dist', 'vector', 'error'])

    cif_files = sorted(cif_map.values())
    results = []
    done = 0
    with Pool(processes=args.workers, initializer=init_worker, initargs=(refs,)) as pool:
        for row in pool.imap_unordered(process_one, cif_files, chunksize=10):
            results.append(row)
            done += 1
            if done % 50 == 0:
                print(f"  {done}/{len(cif_files)}...", flush=True)

    results.sort(key=lambda r: r['accession'])

    with open(args.output, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        w.writeheader()
        w.writerows(results)

    errors = sum(1 for r in results if r.get('error'))
    print(f"\nDone. {len(results)} structures, {errors} errors.", flush=True)

    # Per-reference usage
    ref_usage = Counter(r.get('best_ref', '') for r in results if not r.get('error'))
    print(f"\nBest reference usage:")
    for ref, n in ref_usage.most_common():
        print(f"  {ref}: {n}")

    # RMSD distribution
    rmsds = [float(r['best_rmsd']) for r in results if r.get('best_rmsd', 'NA') != 'NA']
    if rmsds:
        print(f"\nRMSD: min={min(rmsds):.3f} median={sorted(rmsds)[len(rmsds)//2]:.3f} "
              f"max={max(rmsds):.3f}")

    vectors = Counter(r.get('vector', '') for r in results if not r.get('error'))
    print(f"\nUnique vectors: {len(vectors)}")
    print("Top 10:")
    for v, n in vectors.most_common(10):
        print(f"  {n:>4}  {v}")

    thio = Counter(r.get('thioether', '') for r in results if not r.get('error'))
    print(f"\nThioether detection:")
    for label, n in thio.most_common():
        print(f"  {label:>3}: {n}")


if __name__ == '__main__':
    main()
