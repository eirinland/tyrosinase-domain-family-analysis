#!/usr/bin/env python3
"""Re-triage the three pools using the copper-anchored helicity core test.

Replaces the 2026-06-14 (M1 + cealign + pLDDT-floor) core decision for the
<6-His failed-canonical scope ONLY. Canonical and the six-His fails are untouched.

New rule (32,069 structures):
  canonical (canonical==True)                       -> canonical   (unchanged)
  failed-canonical, n_his>=6                         -> discarded   (six-His, out of scope)
  failed-canonical, n_his<6  (the 9,898 scope):
        core_ok (new copper-anchored test) == True  -> non-canonical
        else                                        -> discarded

Reports the delta vs the old pools: structures rescued (old discarded ->
non-canonical) and removed (old non-canonical -> discarded), each broken down by
the old discard/keep reason, plus His-spread of the new non-canonical pool.
Stdlib only.
"""
import argparse, csv, os
from collections import Counter, defaultdict


def nhis(v):
    try:
        return int(float(v))
    except (ValueError, TypeError):
        return -1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--assignment', required=True)  # three_pool_assignment_final.csv
    ap.add_argument('--core', required=True)         # core_helix_results.tsv
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args()

    old = {r['accession']: r for r in csv.DictReader(open(a.assignment))}
    core = {}
    with open(a.core) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            core[r['accession']] = r

    rows = []
    cnt = Counter()
    transitions = Counter()
    rescued, removed = [], []
    nc_his = Counter()
    no_eval = []  # in-scope but no core result row (should be 0)

    for acc, o in old.items():
        canon = o['canonical'] == 'True'
        nh = nhis(o['n_his'])
        old_pool = o['pool']
        cr = core.get(acc)
        ok = bool(cr) and cr.get('core_ok') == 'True'

        if canon:
            new_pool = 'canonical'
        elif nh >= 6:
            new_pool = 'discarded'         # six-His fail, out of scope
        else:
            new_pool = 'non-canonical' if ok else 'discarded'
            if not cr:
                no_eval.append(acc)

        cnt[new_pool] += 1
        transitions[(old_pool, new_pool)] += 1
        if new_pool == 'non-canonical':
            nc_his[nh] += 1

        row = dict(accession=acc, new_pool=new_pool, old_pool=old_pool,
                   canonical=o['canonical'], n_his=o['n_his'],
                   core_ok_new=('True' if ok else ('' if not cr else 'False')),
                   best_ref=(cr or {}).get('best_ref', ''),
                   best_qtm=(cr or {}).get('best_qtm', ''),
                   n_helix_ok=(cr or {}).get('n_helix_ok', ''),
                   old_reason=o['reason'])
        # carry per-helix detail for inspection samples
        for h in (1, 2, 3, 4):
            row[f'a{h}_dist'] = (cr or {}).get(f'a{h}_dist', '')
            row[f'a{h}_plddt'] = (cr or {}).get(f'a{h}_plddt', '')
            row[f'a{h}_helical'] = (cr or {}).get(f'a{h}_helical', '')
            row[f'a{h}_sat'] = (cr or {}).get(f'a{h}_sat', '')
        rows.append(row)

        if old_pool == 'discarded' and new_pool == 'non-canonical':
            rescued.append(row)
        elif old_pool == 'noncanonical' and new_pool == 'discarded':
            removed.append(row)

    # ---- write full new assignment ----
    fields = ['accession', 'new_pool', 'old_pool', 'canonical', 'n_his',
              'core_ok_new', 'best_ref', 'best_qtm', 'n_helix_ok', 'old_reason']
    for h in (1, 2, 3, 4):
        fields += [f'a{h}_dist', f'a{h}_plddt', f'a{h}_helical', f'a{h}_sat']
    out_assign = os.path.join(a.outdir, 'three_pool_assignment_corehelix.csv')
    with open(out_assign, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

    # ---- accession lists ----
    for pool, fn in (('non-canonical', 'noncanonical_accessions_corehelix.csv'),
                     ('discarded', 'discarded_accessions_corehelix.csv'),
                     ('canonical', 'canonical_accessions_corehelix.csv')):
        with open(os.path.join(a.outdir, fn), 'w', newline='') as f:
            w = csv.writer(f); w.writerow(['accession'])
            for r in rows:
                if r['new_pool'] == pool:
                    w.writerow([r['accession']])

    # ---- transition tables (sorted, for inspection) ----
    def sort_for_view(lst):
        # most interesting first: lowest best_qtm last? show by n_helix_ok then qtm
        def key(r):
            try:
                q = float(r['best_qtm']) if r['best_qtm'] != '' else -1
            except ValueError:
                q = -1
            return (-q,)
        return sorted(lst, key=key)

    for lst, fn in ((sort_for_view(rescued), 'rescued_from_discarded.csv'),
                    (sort_for_view(removed), 'removed_from_noncanonical.csv')):
        with open(os.path.join(a.outdir, fn), 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction='ignore')
            w.writeheader(); w.writerows(lst)

    # ---- summary ----
    lines = []
    P = lines.append
    P("=== copper-anchored core-test re-triage ===")
    P(f"input assignment : {a.assignment}")
    P(f"core results     : {a.core}  ({len(core)} rows)")
    P("")
    P("NEW pools (32,069):")
    for k in ('canonical', 'non-canonical', 'discarded'):
        P(f"  {k:14}: {cnt[k]}")
    P("")
    P("OLD pools (2026-06-14 final): canonical 21,893 / non-canonical 1,106 / discarded 9,070")
    P("")
    P(f"  non-canonical delta : {cnt['non-canonical']} - 1106 = {cnt['non-canonical']-1106:+d}")
    P(f"  discarded delta     : {cnt['discarded']} - 9070 = {cnt['discarded']-9070:+d}")
    P("")
    P("transitions (old_pool -> new_pool):")
    for (op, np_), n in sorted(transitions.items(), key=lambda kv: -kv[1]):
        flag = ''
        if op == 'discarded' and np_ == 'non-canonical':
            flag = '   <-- RESCUED'
        elif op == 'noncanonical' and np_ == 'discarded':
            flag = '   <-- REMOVED'
        P(f"  {op:14} -> {np_:14}: {n}{flag}")
    P("")
    P(f"RESCUED (old discarded -> non-canonical): {len(rescued)}")
    rc = Counter()
    for r in rescued:
        rc[r['old_reason'].split('+floor')[0]] += 1  # collapse floorNN
    for reason, n in rc.most_common():
        P(f"    {reason:24}: {n}")
    P(f"  rescued His-spread: {dict(sorted(Counter(nhis(r['n_his']) for r in rescued).items()))}")
    P("")
    P(f"REMOVED (old non-canonical -> discarded): {len(removed)}")
    P(f"  removed His-spread: {dict(sorted(Counter(nhis(r['n_his']) for r in removed).items()))}")
    P("")
    P(f"NEW non-canonical His-spread: {dict(sorted(nc_his.items()))}")
    if no_eval:
        P("")
        P(f"WARNING: {len(no_eval)} in-scope accessions had NO core result row "
          f"(treated as discard): {no_eval[:10]}{' ...' if len(no_eval)>10 else ''}")
    P("")
    P("written: three_pool_assignment_corehelix.csv, {noncanonical,discarded,canonical}_accessions_corehelix.csv,")
    P("         rescued_from_discarded.csv, removed_from_noncanonical.csv, retriage_summary.txt")

    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(a.outdir, 'retriage_summary.txt'), 'w').write(txt + "\n")


if __name__ == '__main__':
    main()
