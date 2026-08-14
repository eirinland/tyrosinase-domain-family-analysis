#!/usr/bin/env python3
"""Re-score BGC context with kingdom-expanded markers.

Uses existing pfam_map.tsv + neighbourhoods.tsv — no new API calls.

Original markers (fungal-only TFs):
  TF: PF00172, PF04082, PF11951 (all fungal Zn2-Cys6 / fungal-specific TF)
  transporter: PF07690, PF00083, PF00005, PF00664
  P450: PF00067
  OMT: PF00891

Expanded markers:
  TF_fungal: PF00172, PF04082, PF11951 (unchanged)
  TF_bacterial: PF00126 (LysR), PF00440 (TetR), PF01047 (MarR/Lrp),
                PF00165 (AraC), PF00392 (GntR)
  transporter: unchanged (already kingdom-agnostic)
  P450: unchanged
  OMT: unchanged
  biosynthetic: PF00501 (AMP-binding/NRPS), PF00109 (Beta-ketoacyl synthase/PKS),
                PF02801 (PKS C-term), PF00668 (Condensation domain/NRPS)
  melC1: PF06236 (Tyrosinase co-factor MelC1) — bacterial operon marker

Scoring:
  UstYa_cluster: ≥1 UstYa (PF11807) — unchanged
  melC1_operon: ≥1 MelC1 (PF06236) — NEW, bacterial-specific
  BGC_accessory: ≥2 of {TF_fungal, TF_bacterial, transporter, P450, OMT, biosynthetic}
  none: everything else
"""

import csv
from collections import Counter, defaultdict
from pathlib import Path

GROUPS_DIR = Path(__file__).parent / 'groups'

USTYA = 'PF11807'
MELC1 = 'PF06236'

MARKERS = {
    'TF_fungal':     {'PF00172', 'PF04082', 'PF11951'},
    'TF_bacterial':  {'PF00126', 'PF00440', 'PF01047', 'PF00165', 'PF00392'},
    'transporter':   {'PF07690', 'PF00083', 'PF00005', 'PF00664'},
    'P450':          {'PF00067'},
    'OMT':           {'PF00891'},
    'biosynthetic':  {'PF00501', 'PF00109', 'PF02801', 'PF00668'},
}

# Also keep the old fungal-only scoring for comparison
MARKERS_OLD = {
    'TF':          {'PF00172', 'PF04082', 'PF11951'},
    'transporter': {'PF07690', 'PF00083', 'PF00005', 'PF00664'},
    'P450':        {'PF00067'},
    'OMT':         {'PF00891'},
}


def score_group(group_dir):
    nb_file = group_dir / 'neighbourhoods.tsv'
    pfam_file = group_dir / 'pfam_map.tsv'
    if not nb_file.exists() or not pfam_file.exists():
        return None

    pid2pfam = {}
    with open(pfam_file) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            fams = set()
            for entry in r['pfam'].split(';'):
                fid = entry.split('|')[0].strip()
                if fid.startswith('PF'):
                    fams.add(fid)
            if fams:
                pid2pfam[r['protein_id']] = fams

    nb_rows = list(csv.DictReader(open(nb_file), delimiter='\t'))
    queries = sorted({r['query_accession'] for r in nb_rows})
    neigh = [r for r in nb_rows if r['is_target'] != '1' and r.get('protein_id')]

    nb_by_q = defaultdict(list)
    for r in neigh:
        nb_by_q[r['query_accession']].append((int(r['offset']), r['protein_id'].split('.')[0]))

    rows = []
    for q in queries:
        fams_all = set()
        ust_offsets = []
        melc1_offsets = []
        for off, pid in nb_by_q.get(q, []):
            fams = pid2pfam.get(pid, set())
            fams_all |= fams
            if USTYA in fams:
                ust_offsets.append(abs(off))
            if MELC1 in fams:
                melc1_offsets.append(abs(off))

        n_ust = len(ust_offsets)
        n_melc1 = len(melc1_offsets)

        marker_hits = {k: bool(fams_all & v) for k, v in MARKERS.items()}
        n_marker_cats = sum(marker_hits.values())
        old_marker_cats = sum(bool(fams_all & v) for v in MARKERS_OLD.values())

        if n_ust >= 1:
            bgc_new = 'UstYa_cluster'
        elif n_melc1 >= 1:
            bgc_new = 'melC1_operon'
        elif n_marker_cats >= 2:
            bgc_new = 'BGC_accessory'
        else:
            bgc_new = 'none'

        if n_ust >= 1:
            bgc_old = 'UstYa_cluster'
        elif old_marker_cats >= 2:
            bgc_old = 'BGC_accessory'
        else:
            bgc_old = 'none'

        rows.append({
            'accession': q,
            'n_neighbours': len(nb_by_q.get(q, [])),
            'n_ustya': n_ust,
            'closest_ustya': min(ust_offsets) if ust_offsets else '',
            'n_melc1': n_melc1,
            'has_TF_fungal': int(marker_hits['TF_fungal']),
            'has_TF_bacterial': int(marker_hits['TF_bacterial']),
            'has_transporter': int(marker_hits['transporter']),
            'has_P450': int(marker_hits['P450']),
            'has_OMT': int(marker_hits['OMT']),
            'has_biosynthetic': int(marker_hits['biosynthetic']),
            'bgc_context': bgc_new,
            'bgc_old': bgc_old,
        })

    fields = ['accession', 'n_neighbours', 'n_ustya', 'closest_ustya', 'n_melc1',
              'has_TF_fungal', 'has_TF_bacterial', 'has_transporter',
              'has_P450', 'has_OMT', 'has_biosynthetic', 'bgc_context', 'bgc_old']
    out_file = group_dir / 'bgc_scoring_expanded.tsv'
    with open(out_file, 'w', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    n_ust = sum(1 for r in rows if r['bgc_context'] == 'UstYa_cluster')
    n_melc1 = sum(1 for r in rows if r['bgc_context'] == 'melC1_operon')
    n_bgc = sum(1 for r in rows if r['bgc_context'] == 'BGC_accessory')
    n_none = sum(1 for r in rows if r['bgc_context'] == 'none')
    n_ust_old = sum(1 for r in rows if r['bgc_old'] == 'UstYa_cluster')
    n_bgc_old = sum(1 for r in rows if r['bgc_old'] == 'BGC_accessory')
    n_none_old = sum(1 for r in rows if r['bgc_old'] == 'none')

    return {
        'n': len(rows), 'ustya': n_ust, 'melc1': n_melc1, 'bgc_acc': n_bgc, 'none': n_none,
        'ustya_old': n_ust_old, 'bgc_old': n_bgc_old, 'none_old': n_none_old,
    }


def main():
    print(f'{"Group":<28s} {"N":>4s}  {"UstYa":>5s} {"melC1":>5s} {"BGCac":>5s} {"none":>5s}  |  {"UstYa":>5s} {"BGCac":>5s} {"none":>5s}')
    print(f'{"":28s} {"":>4s}  {"--- expanded ---":^22s}  |  {"--- old ---":^17s}')
    print('-' * 90)

    for group_dir in sorted(GROUPS_DIR.iterdir()):
        if not group_dir.is_dir():
            continue
        result = score_group(group_dir)
        if result is None:
            print(f'{group_dir.name:<28s}  SKIPPED (missing files)')
            continue
        r = result
        bgc_pct = (r['ustya'] + r['melc1'] + r['bgc_acc']) / r['n'] * 100 if r['n'] else 0
        bgc_old_pct = (r['ustya_old'] + r['bgc_old']) / r['n'] * 100 if r['n'] else 0
        print(f'{group_dir.name:<28s} {r["n"]:>4d}  {r["ustya"]:>5d} {r["melc1"]:>5d} {r["bgc_acc"]:>5d} {r["none"]:>5d} ({bgc_pct:4.1f}%)  |  {r["ustya_old"]:>5d} {r["bgc_old"]:>5d} {r["none_old"]:>5d} ({bgc_old_pct:4.1f}%)')


if __name__ == '__main__':
    main()
