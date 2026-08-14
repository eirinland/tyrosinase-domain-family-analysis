#!/usr/bin/env python3
"""Pick a labelled inspection sample from the M7 re-triage for visual PyMOL review.

Four labels = the four old->new pool transitions (mutually exclusive):
  new_rescue          discarded   -> non-canonical   (flipped up; all 4/4 core)
  new_discard         noncanonical-> discarded       (flipped down; all <4/4, incl. SYG4)
  borderline_kept     noncanonical-> non-canonical   (stable keep; tightest 4/4 pass)
  borderline_discarded discarded  -> discarded        (stable discard; 3/4, missing helix just failed)

Writes borderlines_manifest.tsv (accession, label, out_filename, metrics, why).
Stdlib only.
"""
import argparse, csv

PER = 6


def f(v, d=None):
    try:
        return float(v)
    except (ValueError, TypeError):
        return d


def truthy(v):
    return str(v) in ('1', 'True', 'true')


def helices(r):
    return [(f(r[f'a{h}_dist'], 99.0), f(r[f'a{h}_plddt'], 0.0),
             truthy(r[f'a{h}_helical']), truthy(r[f'a{h}_sat'])) for h in (1, 2, 3, 4)]


def qtm(r):
    return f(r['best_qtm'], -1.0)


def nhel(r):
    try:
        return int(float(r['n_helix_ok']))
    except (ValueError, TypeError):
        return -1


def tightness(r):
    """For a 4/4 keep: how close the weakest helix is to failing (smaller = more borderline)."""
    hs = helices(r)
    d_room = 4.0 - max(h[0] for h in hs)        # distance headroom to dmax=4.0
    p_room = min(h[1] for h in hs) - 70.0       # pLDDT headroom to pmin=70
    return min(d_room / 4.0, p_room / 30.0)


def deficit(r):
    """For a 3/4 discard: how badly the single failing (but structurally helical) helix missed."""
    fails = [h for h in helices(r) if not h[3] and h[2]]   # sat=False but helical=True
    if not fails:
        return None
    h = min(fails, key=lambda x: max(0.0, x[0] - 4.0) + max(0.0, 70.0 - x[1]) / 50.0)
    return max(0.0, h[0] - 4.0) + max(0.0, 70.0 - h[1]) / 50.0


def hsumm(r):
    return " ".join(f"a{i+1}(d={h[0]:.2f},pl={h[1]:.0f},{'H' if h[2] else '-'}{'+' if h[3] else 'x'})"
                    for i, h in enumerate(helices(r)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--assign', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.assign)))
    scope = [r for r in rows if r['core_ok_new'] in ('True', 'False')]
    op, npp = (lambda r: r['old_pool']), (lambda r: r['new_pool'])

    rescue = [r for r in scope if op(r) == 'discarded' and npp(r) == 'non-canonical']
    discard = [r for r in scope if op(r) == 'noncanonical' and npp(r) == 'discarded']
    keep = [r for r in scope if op(r) == 'noncanonical' and npp(r) == 'non-canonical']
    stable_disc = [r for r in scope if op(r) == 'discarded' and npp(r) == 'discarded']

    picks = []  # (label, row, why)

    # new_rescue: 4 lowest-qTM (riskiest) + 2 highest (clear-good controls)
    rs = sorted(rescue, key=qtm)
    chosen = rs[:4] + rs[-2:]
    seen = set()
    for r in chosen:
        if r['accession'] in seen:
            continue
        seen.add(r['accession'])
        tag = 'low-qTM rescue (riskiest)' if qtm(r) <= qtm(rs[3]) else 'high-qTM rescue (control)'
        picks.append(('new_rescue', r, f"{tag}; full 4/4 core M7 validated, cealign gate had cut it"))

    # new_discard: SYG4 reference + 5 highest-qTM 3/4 cuts (most contested)
    d3 = sorted([r for r in discard if nhel(r) == 3], key=qtm, reverse=True)
    syg = next((r for r in discard if r['accession'] == 'A0A8H8SYG4'), None)
    chosen = ([syg] if syg else []) + [r for r in d3 if r['accession'] != 'A0A8H8SYG4'][:5]
    for r in chosen[:PER]:
        why = "the motivating half-bundle (2/4)" if r['accession'] == 'A0A8H8SYG4' \
              else f"best-aligned cut (qTM {qtm(r):.2f}) but only {nhel(r)}/4 helices"
        picks.append(('new_discard', r, why))

    # borderline_kept: tightest 4/4 pass among stable keeps
    bk = sorted([r for r in keep if nhel(r) == 4], key=tightness)[:PER]
    for r in bk:
        picks.append(('borderline_kept', r, f"weakest helix near threshold (tightness {tightness(r):.3f})"))

    # borderline_discarded: 3/4 stable discards whose missing helix barely failed
    cand = [r for r in stable_disc if nhel(r) == 3 and deficit(r) is not None]
    bd = sorted(cand, key=deficit)[:PER]
    for r in bd:
        picks.append(('borderline_discarded', r, f"missing helix just failed (deficit {deficit(r):.3f})"))

    cols = ['accession', 'label', 'out_filename', 'old_pool', 'new_pool',
            'best_qtm', 'n_helix_ok', 'helices', 'why']
    with open(a.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(cols)
        for label, r, why in picks:
            acc, q, nh = r['accession'], qtm(r), nhel(r)
            outfn = f"{label}__{acc}__qtm{q:.3f}_nhel{nh}.cif"
            w.writerow([acc, label, outfn, r['old_pool'], r['new_pool'],
                        f"{q:.3f}", nh, hsumm(r), why])

    # readable summary
    print(f"selected {len(picks)} structures "
          f"(rescue {sum(p[0]=='new_rescue' for p in picks)}, "
          f"discard {sum(p[0]=='new_discard' for p in picks)}, "
          f"bord-kept {sum(p[0]=='borderline_kept' for p in picks)}, "
          f"bord-disc {sum(p[0]=='borderline_discarded' for p in picks)})\n")
    last = None
    for label, r, why in picks:
        if label != last:
            print(f"== {label} ==")
            last = label
        print(f"  {r['accession']:13} qtm={qtm(r):.3f} nhel={nhel(r)}  {hsumm(r)}  | {why}")


if __name__ == '__main__':
    main()
