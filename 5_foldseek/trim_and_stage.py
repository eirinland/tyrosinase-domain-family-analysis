"""
Trim low-pLDDT terminal regions from AlphaFold3 CIF files and stage to output dir.

Identical to the parent directory version but reads 'sequence_id' column
from the accession CSV (matching the filtering pipeline output format).
"""

import argparse
import csv
import os
from multiprocessing import Pool
from pathlib import Path


PLDDT_COL = '_atom_site.B_iso_or_equiv'


def parse_cif(path):
    header_lines = []
    col_names = []
    rows = []
    in_atom_loop = False
    collecting_cols = False

    with open(path) as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            stripped = line.strip()

            if in_atom_loop:
                if stripped.startswith('_') or stripped in ('#', 'loop_') or not stripped:
                    break
                parts = stripped.split()
                if len(parts) == len(col_names):
                    rows.append(parts)
                continue

            if stripped == 'loop_':
                collecting_cols = False
                col_names = []
                header_lines.append(raw)
                continue

            if stripped.startswith('_atom_site.'):
                collecting_cols = True
                col_names.append(stripped)
                header_lines.append(raw)
                continue

            if collecting_cols and col_names and not stripped.startswith('_atom_site.'):
                collecting_cols = False
                in_atom_loop = True
                parts = stripped.split()
                if len(parts) == len(col_names):
                    rows.append(parts)
                continue

            header_lines.append(raw)

    return col_names, rows, header_lines


def get_residue_plddt(col_names, rows):
    seq_col  = col_names.index('_atom_site.label_seq_id')
    atom_col = col_names.index('_atom_site.label_atom_id')
    b_col    = col_names.index(PLDDT_COL)

    res_ca    = {}
    res_all   = {}
    res_order = []

    for row in rows:
        seq_id  = row[seq_col]
        atom    = row[atom_col].upper()
        try:
            b = float(row[b_col])
        except ValueError:
            continue

        if seq_id not in res_all:
            res_all[seq_id] = []
            res_order.append(seq_id)
        res_all[seq_id].append(b)

        if atom == 'CA':
            res_ca[seq_id] = b

    result = []
    for seq_id in res_order:
        if seq_id in res_ca:
            result.append((seq_id, res_ca[seq_id]))
        elif res_all[seq_id]:
            result.append((seq_id, sum(res_all[seq_id]) / len(res_all[seq_id])))

    return result


def find_trim_bounds(res_plddt, cutoff, window):
    n = len(res_plddt)
    if n < window:
        return None

    start_idx = None
    for i in range(n - window + 1):
        if all(res_plddt[i + j][1] >= cutoff for j in range(window)):
            start_idx = i
            break

    end_idx = None
    for i in range(n - 1, window - 2, -1):
        if all(res_plddt[i - j][1] >= cutoff for j in range(window)):
            end_idx = i
            break

    if start_idx is None or end_idx is None or start_idx >= end_idx:
        return None

    kept_ids = {res_plddt[i][0] for i in range(start_idx, end_idx + 1)}
    n_trimmed_n = start_idx
    n_trimmed_c = n - 1 - end_idx
    return kept_ids, n_trimmed_n, n_trimmed_c


def write_trimmed_cif(col_names, rows, header_lines, kept_ids, out_path):
    seq_col = col_names.index('_atom_site.label_seq_id')
    kept_rows = [r for r in rows if r[seq_col] in kept_ids]
    if not kept_rows:
        return False
    with open(out_path, 'w') as fh:
        for line in header_lines:
            fh.write(line if line.endswith('\n') else line + '\n')
        for row in kept_rows:
            fh.write(' '.join(row) + '\n')
        fh.write('#\n')
    return True


def process_one(args):
    acc, cif_dir, out_dir, cutoff, window = args
    candidates = [
        os.path.join(cif_dir, f'{acc}.cif'),
        os.path.join(cif_dir, f'{acc}_model_A.cif'),
    ]
    cif_path = None
    for c in candidates:
        if os.path.exists(c):
            cif_path = c
            break
    if cif_path is None:
        for fname in os.listdir(cif_dir):
            if fname.startswith(acc) and fname.endswith('.cif'):
                cif_path = os.path.join(cif_dir, fname)
                break

    if cif_path is None:
        return acc, 'not_found', 0, 0

    try:
        col_names, rows, header_lines = parse_cif(cif_path)
        if PLDDT_COL not in col_names:
            return acc, 'no_plddt_col', 0, 0

        res_plddt = get_residue_plddt(col_names, rows)
        result = find_trim_bounds(res_plddt, cutoff, window)

        if result is None:
            return acc, 'no_valid_window', 0, 0

        kept_ids, n_trim_n, n_trim_c = result
        out_path = os.path.join(out_dir, os.path.basename(cif_path))
        ok = write_trimmed_cif(col_names, rows, header_lines, kept_ids, out_path)
        if not ok:
            return acc, 'empty_after_trim', 0, 0

        return acc, 'ok', n_trim_n, n_trim_c

    except Exception as e:
        return acc, f'error:{e}', 0, 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cif-dir',    required=True)
    ap.add_argument('--accessions', required=True)
    ap.add_argument('--output',     required=True)
    ap.add_argument('--cutoff',     type=float, default=70.0)
    ap.add_argument('--window',     type=int,   default=5)
    ap.add_argument('--workers',    type=int,   default=16)
    args = ap.parse_args()

    os.makedirs(args.output, exist_ok=True)

    accessions = []
    with open(args.accessions) as f:
        for row in csv.DictReader(f):
            accessions.append(row['sequence_id'])

    print(f'Accessions to process : {len(accessions)}')
    print(f'pLDDT cutoff          : {args.cutoff}')
    print(f'Window size           : {args.window}')
    print(f'Workers               : {args.workers}')
    print(flush=True)

    tasks = [(acc, args.cif_dir, args.output, args.cutoff, args.window)
             for acc in accessions]

    counts = {'ok': 0, 'not_found': 0, 'no_valid_window': 0,
              'no_plddt_col': 0, 'empty_after_trim': 0, 'error': 0}
    trim_n_total = trim_c_total = 0

    with Pool(args.workers) as pool:
        for i, (acc, status, tn, tc) in enumerate(pool.imap_unordered(process_one, tasks), 1):
            key = status if status in counts else 'error'
            counts[key] += 1
            trim_n_total += tn
            trim_c_total += tc
            if i % 2000 == 0:
                print(f'  {i}/{len(tasks)} processed ...', flush=True)

    n_ok = counts['ok']
    print(f'\nDone.')
    print(f'  Written (ok)          : {n_ok}')
    print(f'  Not found in mount    : {counts["not_found"]}')
    print(f'  Skipped (no window)   : {counts["no_valid_window"]}')
    print(f'  Skipped (other)       : {counts["no_plddt_col"] + counts["empty_after_trim"] + counts["error"]}')
    if n_ok:
        print(f'  Mean N-term trimmed   : {trim_n_total/n_ok:.1f} residues')
        print(f'  Mean C-term trimmed   : {trim_c_total/n_ok:.1f} residues')


if __name__ == '__main__':
    main()
