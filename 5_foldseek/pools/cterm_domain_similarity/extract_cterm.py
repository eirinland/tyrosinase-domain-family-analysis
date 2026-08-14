"""
Extract each cluster member's C-terminal non-PPO domain into a standalone CIF
for within-cluster structural similarity (Foldseek all-vs-all).

Reads untrimmed AF3 CIFs from a squashfuse mount, keeps only atom rows whose
label_seq_id falls in a C-terminal non-PPO segment (segment start > PPO-core end).
Output CIF is named by bare accession so Foldseek tags one uniform chain (_A).
Stdlib only -> runs inside python_tools.sif.
"""
import argparse, csv, os, re
from multiprocessing import Pool


def parse_cif(path):
    col_names, rows, header_lines = [], [], []
    in_atom_loop = collecting = False
    with open(path) as fh:
        for raw in fh:
            s = raw.rstrip('\n').strip()
            if in_atom_loop:
                if s.startswith('_') or s in ('#', 'loop_') or not s:
                    break
                parts = s.split()
                if len(parts) == len(col_names):
                    rows.append(parts)
                continue
            if s == 'loop_':
                collecting = False; col_names = []; header_lines.append(raw); continue
            if s.startswith('_atom_site.'):
                collecting = True; col_names.append(s); header_lines.append(raw); continue
            if collecting and col_names and not s.startswith('_atom_site.'):
                collecting = False; in_atom_loop = True
                parts = s.split()
                if len(parts) == len(col_names):
                    rows.append(parts)
                continue
            header_lines.append(raw)
    return col_names, rows, header_lines


def parse_ranges(s):
    out = []
    if not s or s in ('-', 'nan', ''):
        return out
    for seg in s.split(','):
        m = re.match(r'(\d+)-(\d+)', seg.strip())
        if m:
            out.append((int(m.group(1)), int(m.group(2))))
    return out


def cterm_segments(ppo, nonppo):
    if not ppo or not nonppo:
        return []
    pe = max(e for _, e in ppo)
    return [(s, e) for s, e in nonppo if s > pe]


def write_cif(col_names, rows, header_lines, segs, out_path):
    seq_col = col_names.index('_atom_site.label_seq_id')

    def keep(v):
        try:
            i = int(v)
        except ValueError:
            return False
        return any(lo <= i <= hi for lo, hi in segs)

    kept = [r for r in rows if keep(r[seq_col])]
    if not kept:
        return False
    with open(out_path, 'w') as fh:
        for line in header_lines:
            fh.write(line if line.endswith('\n') else line + '\n')
        for r in kept:
            fh.write(' '.join(r) + '\n')
        fh.write('#\n')
    return True


def process_one(args):
    member_raw, cif_dir, out_dir, segs = args
    # mount CIFs are uniformly '..._model.cif'; some cluster-TSV member names
    # carry a spurious '_A' chain tag -> strip it before lookup.
    lookup = member_raw[:-2] if member_raw.endswith('_A') else member_raw
    cif = os.path.join(cif_dir, lookup + '.cif')
    if not os.path.exists(cif):
        return member_raw, 'not_found'
    try:
        cn, rows, hdr = parse_cif(cif)
        if '_atom_site.label_seq_id' not in cn:
            return member_raw, 'no_seq_col'
        acc = member_raw.split('_taxID_')[0]
        ok = write_cif(cn, rows, hdr, segs, os.path.join(out_dir, acc + '.cif'))
        return member_raw, 'ok' if ok else 'empty'
    except Exception as e:
        return member_raw, f'error:{e}'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cif-dir', required=True)
    ap.add_argument('--cluster-tsv', required=True)
    ap.add_argument('--ppo-csv', required=True)
    ap.add_argument('--rep', required=True, help='bare accession of cluster rep')
    ap.add_argument('--output', required=True)
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--min-len', type=int, default=40)
    ap.add_argument('--workers', type=int, default=16)
    args = ap.parse_args()
    os.makedirs(args.output, exist_ok=True)

    ppo = {}
    with open(args.ppo_csv) as f:
        for r in csv.DictReader(f):
            ppo[r['accession']] = (parse_ranges(r['ppo_range']),
                                   parse_ranges(r['non_ppo_range']))

    members = []
    with open(args.cluster_tsv) as f:
        for line in f:
            rep, mem = line.rstrip('\n').split('\t')
            if rep.split('_taxID_')[0] == args.rep:
                members.append(mem)

    tasks, meta, skipped = [], {}, []
    for mem in members:
        acc = mem.split('_taxID_')[0]
        d = ppo.get(acc)
        if not d:
            skipped.append((mem, acc, '', '', 'no_ppo_row')); continue
        segs = cterm_segments(*d)
        clen = sum(e - s + 1 for s, e in segs)
        if not segs or clen < args.min_len:
            skipped.append((mem, acc, clen, len(segs), 'no_cterm')); continue
        tasks.append((mem, args.cif_dir, args.output, segs))
        meta[mem] = (acc, clen, len(segs))

    print(f'rep {args.rep}: {len(members)} members, '
          f'{len(tasks)} with C-term >= {args.min_len} res, {len(skipped)} skipped',
          flush=True)

    results = {}
    counts = {}
    with Pool(args.workers) as p:
        for mem, st in p.imap_unordered(process_one, tasks):
            results[mem] = st
            k = st.split(':')[0]
            counts[k] = counts.get(k, 0) + 1

    with open(args.manifest, 'w', newline='') as man:
        mw = csv.writer(man)
        mw.writerow(['member', 'accession', 'cterm_len', 'n_seg', 'status'])
        for mem in members:
            if mem in meta:
                acc, clen, nseg = meta[mem]
                mw.writerow([mem, acc, clen, nseg, results.get(mem, 'missing')])
        for row in skipped:
            mw.writerow(list(row))

    print('extract status:', counts, flush=True)


if __name__ == '__main__':
    main()
