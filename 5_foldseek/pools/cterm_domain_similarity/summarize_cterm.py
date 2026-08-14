"""
Summarize within-cluster C-terminal domain structural similarity.

Inputs (per cluster):
  --allvsall   foldseek easy-search TSV: query,target,qtmscore,ttmscore,alnlen,qlen,tlen
  --manifest   extract_cterm.py manifest (member,accession,cterm_len,n_seg,status)
  --clust05 / --clust07   foldseek easy-cluster cluster_cluster.tsv (rep<TAB>member)
Reports: pairwise min(qTM,tTM) distribution, nearest-neighbour TM, sub-cluster
structure (n clusters, largest fraction) and median C-term length per sub-cluster
(to test whether a length bimodality maps onto distinct structural sub-populations).
Stdlib only.
"""
import argparse, csv, statistics as st
from collections import defaultdict


def strip(name):
    return name[:-2] if name.endswith('_A') else name


def pct(sorted_vals, p):
    if not sorted_vals:
        return float('nan')
    return sorted_vals[min(len(sorted_vals) - 1, int(p * len(sorted_vals)))]


def load_lengths(manifest):
    L = {}
    with open(manifest) as f:
        for r in csv.DictReader(f):
            if r['status'] == 'ok':
                try:
                    L[r['accession']] = int(r['cterm_len'])
                except (ValueError, KeyError):
                    pass
    return L


def cluster_report(path, lengths):
    members = defaultdict(list)
    with open(path) as f:
        for line in f:
            rep, mem = line.rstrip('\n').split('\t')
            members[strip(rep)].append(strip(mem))
    sizes = sorted(members.items(), key=lambda kv: len(kv[1]), reverse=True)
    n_struct = sum(len(v) for v in members.values())
    return members, sizes, n_struct


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--label', required=True)
    ap.add_argument('--allvsall', required=True)
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--clust05', required=True)
    ap.add_argument('--clust07', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    lengths = load_lengths(args.manifest)

    vals = []
    nn = {}
    n_struct_seen = set()
    with open(args.allvsall) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 4:
                continue
            q, t = strip(p[0]), strip(p[1])
            n_struct_seen.add(q)
            if q == t:
                continue
            try:
                m = min(float(p[2]), float(p[3]))
            except ValueError:
                continue
            vals.append(m)
            if m > nn.get(q, -1.0):
                nn[q] = m

    vals.sort()
    nnv = sorted(nn.values())
    n = len(vals)

    lines = []
    def out(s=''):
        lines.append(s)
        print(s, flush=True)

    out(f'=== {args.label} : within-cluster C-terminal domain similarity ===')
    out(f'structures with extracted C-term domain : {len(n_struct_seen)}')
    out(f'C-term length (extracted)               : '
        f'median={pct(sorted(lengths.values()),0.5):.0f} '
        f'IQR=[{pct(sorted(lengths.values()),0.25):.0f},'
        f'{pct(sorted(lengths.values()),0.75):.0f}] '
        f'n={len(lengths)}')
    out('')
    out('-- pairwise min(qTM,tTM), all ordered non-self hits --')
    if n:
        out(f'  pairs={n}  median={pct(vals,0.5):.3f}  '
            f'IQR=[{pct(vals,0.25):.3f},{pct(vals,0.75):.3f}]  '
            f'mean={st.mean(vals):.3f} sd={st.pstdev(vals):.3f}  '
            f'min={vals[0]:.3f} max={vals[-1]:.3f}')
        for thr in (0.5, 0.7, 0.8, 0.9):
            f = sum(1 for v in vals if v >= thr) / n
            out(f'    fraction >= {thr:.1f} : {100*f:5.1f}%')
    out('')
    out('-- nearest-neighbour min(qTM,tTM) per structure --')
    if nnv:
        out(f'  n={len(nnv)}  median={pct(nnv,0.5):.3f}  '
            f'IQR=[{pct(nnv,0.25):.3f},{pct(nnv,0.75):.3f}]  '
            f'min={nnv[0]:.3f}')
        for thr in (0.5, 0.7, 0.8):
            f = sum(1 for v in nnv if v >= thr) / len(nnv)
            out(f'    fraction with a neighbour >= {thr:.1f} : {100*f:5.1f}%')
    out('')

    for tag, path in (('TM0.5', args.clust05), ('TM0.7', args.clust07)):
        members, sizes, ns = cluster_report(path, lengths)
        singletons = sum(1 for _, v in sizes if len(v) == 1)
        out(f'-- foldseek easy-cluster {tag} (connected-component) --')
        out(f'  {len(members)} sub-clusters over {ns} domains '
            f'({singletons} singletons); largest={len(sizes[0][1])} '
            f'({100*len(sizes[0][1])/ns:.1f}%)')
        out(f'  {"rep":<14}{"size":>6}{"%":>7}  median_len  IQR_len')
        for rep, mem in sizes[:8]:
            ls = sorted(lengths[m] for m in mem if m in lengths)
            if ls:
                out(f'  {rep:<14}{len(mem):>6}{100*len(mem)/ns:>6.1f}%'
                    f'  {pct(ls,0.5):>9.0f}  [{pct(ls,0.25):.0f},{pct(ls,0.75):.0f}]')
            else:
                out(f'  {rep:<14}{len(mem):>6}{100*len(mem)/ns:>6.1f}%        n/a')
        out('')

    with open(args.out, 'w') as f:
        f.write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
