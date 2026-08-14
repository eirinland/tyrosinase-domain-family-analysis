#!/usr/bin/env python3
"""Score the copper-anchored core test (M7) against the hand-labelled eye set and
compare head-to-head with the six methods already in benchmark_results.tsv.

The benchmark's hard decision is **core present** (belongs in canonical OR
non-canonical) vs **discard**. Ground truth = `your_label`:
    present  := your_label in {canonical, noncanonical}
    absent   := your_label == discard
Only hand-labelled rows are scored. Stdlib only.
"""
import argparse, csv
from collections import Counter, defaultdict


def fnum(v, d=None):
    try:
        return float(v)
    except (ValueError, TypeError):
        return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--results', required=True)   # benchmark_results.tsv (has your_label + m1..m6)
    ap.add_argument('--m7', required=True)        # core_helix_bench.tsv (new filter)
    ap.add_argument('--out-prefix', default='benchmark_core')
    a = ap.parse_args()

    res = list(csv.DictReader(open(a.results), delimiter='\t'))
    m7 = {}
    with open(a.m7) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            m7[r['accession']] = r
    def m7ok(acc):
        r = m7.get(acc)
        return bool(r) and r.get('core_ok') == 'True'

    # ----- method "core-present" calls -----
    methods = {
        'M1 helix(current)':   lambda r: r['m1_core_ok'] == 'True',
        'M2 refTM>=0.50':      lambda r: (fnum(r['m2_best_ttm'], -1) >= 0.50),
        'M4 cealign<=6.0':     lambda r: (fnum(r['m4_ce_cov'], -1) >= 0.55 and fnum(r['m4_ce_rmsd'], 1e9) <= 6.0),
        'M4 cealign<=5.0':     lambda r: (fnum(r['m4_ce_cov'], -1) >= 0.55 and fnum(r['m4_ce_rmsd'], 1e9) <= 5.0),
        'M6 bundle&identity':  lambda r: r['m6_combined'] == 'True',
        'M4<=5.0 AND M1':      lambda r: (fnum(r['m4_ce_cov'], -1) >= 0.55 and fnum(r['m4_ce_rmsd'], 1e9) <= 5.0)
                                          and r['m1_core_ok'] == 'True',
        'M7 Cu-anchored':      lambda r: m7ok(r['accession']),
    }

    labeled = [r for r in res if r['your_label'] in ('canonical', 'noncanonical', 'discard')]
    def truth_present(r): return r['your_label'] in ('canonical', 'noncanonical')

    # ----- confusion per method -----
    print(f"Eye-labelled structures scored: {len(labeled)} "
          f"(present={sum(truth_present(r) for r in labeled)}, "
          f"discard={sum(not truth_present(r) for r in labeled)})\n")
    hdr = f"{'method':20} {'acc%':>6} {'TP':>4} {'FP':>4} {'FN':>4} {'TN':>4} {'prec':>6} {'recall':>7} {'F1':>6}"
    print(hdr); print('-' * len(hdr))
    scoreboard = []
    for name, fn in methods.items():
        tp = fp = fn_ = tn = 0
        for r in labeled:
            call = fn(r); tru = truth_present(r)
            if call and tru: tp += 1
            elif call and not tru: fp += 1
            elif not call and tru: fn_ += 1
            else: tn += 1
        n = len(labeled); acc = 100 * (tp + tn) / n
        prec = tp / (tp + fp) if (tp + fp) else 0.0
        rec = tp / (tp + fn_) if (tp + fn_) else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
        scoreboard.append((name, acc, tp, fp, fn_, tn, prec, rec, f1))
        print(f"{name:20} {acc:6.1f} {tp:4} {fp:4} {fn_:4} {tn:4} {prec:6.2f} {rec:7.2f} {f1:6.2f}")

    best = max(scoreboard, key=lambda s: s[1])
    print(f"\nbest accuracy: {best[0]}  ({best[1]:.1f}%)\n")

    # ----- per-stratum: M7 present-calls vs truth -----
    print("PER-STRATUM (n / truth-present / M7-present / M1-present / M4<=5&M1-present)")
    print("-" * 78)
    by = defaultdict(list)
    for r in labeled:
        by[r['stratum']].append(r)
    m1f = methods['M1 helix(current)']; m41 = methods['M4<=5.0 AND M1']; m7f = methods['M7 Cu-anchored']
    for strat in sorted(by):
        rs = by[strat]; n = len(rs)
        tp = sum(truth_present(r) for r in rs)
        print(f"{strat:20} {n:4} {tp:6}   M7={sum(m7f(r) for r in rs):3}  "
              f"M1={sum(m1f(r) for r in rs):3}  M4&M1={sum(m41(r) for r in rs):3}")

    # ----- the 3 discriminating strata, structure-level M7 errors -----
    print("\nM7 ERRORS on the discriminating strata (acc | your_label | M7 call):")
    print("-" * 78)
    for strat in ('canon_disagree', 'noncanon_lowqtm', 'noncanon_rand', 'discard_fulllen'):
        errs = [r for r in by.get(strat, []) if m7f(r) != truth_present(r)]
        print(f"  {strat}: {len(errs)} M7 errors / {len(by.get(strat, []))}")
        for r in errs[:12]:
            mr = m7.get(r['accession'], {})
            dists = ",".join(mr.get(f'a{h}_dist', '-') for h in (1, 2, 3, 4))
            plds = ",".join(mr.get(f'a{h}_plddt', '-') for h in (1, 2, 3, 4))
            print(f"      {r['accession']:13} truth={r['your_label']:12} "
                  f"M7={'present' if m7f(r) else 'absent':7} nhel={mr.get('n_helix_ok','-')} "
                  f"qtm={mr.get('best_qtm','-')} d=[{dists}] pl=[{plds}]")

    # ----- write per-structure compare -----
    out_tsv = f"{a.out_prefix}_compare.tsv"
    cols = ['accession', 'stratum', 'your_label', 'truth_present',
            'm1', 'm4_50', 'm4_60', 'm6', 'm4m1', 'm7',
            'm7_best_ref', 'm7_qtm', 'm7_nhel',
            'a1_dist', 'a1_plddt', 'a2_dist', 'a2_plddt',
            'a3_dist', 'a3_plddt', 'a4_dist', 'a4_plddt']
    with open(out_tsv, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t'); w.writerow(cols)
        for r in labeled:
            mr = m7.get(r['accession'], {})
            w.writerow([r['accession'], r['stratum'], r['your_label'], int(truth_present(r)),
                        int(methods['M1 helix(current)'](r)),
                        int(methods['M4 cealign<=5.0'](r)),
                        int(methods['M4 cealign<=6.0'](r)),
                        int(methods['M6 bundle&identity'](r)),
                        int(m41(r)), int(m7f(r)),
                        mr.get('best_ref', ''), mr.get('best_qtm', ''), mr.get('n_helix_ok', ''),
                        mr.get('a1_dist', ''), mr.get('a1_plddt', ''),
                        mr.get('a2_dist', ''), mr.get('a2_plddt', ''),
                        mr.get('a3_dist', ''), mr.get('a3_plddt', ''),
                        mr.get('a4_dist', ''), mr.get('a4_plddt', '')])
    print(f"\nwritten {out_tsv}")


if __name__ == '__main__':
    main()
