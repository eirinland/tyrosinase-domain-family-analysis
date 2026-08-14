"""
Test whether Phe65 and Trp68 non-aromatic substitutions are structurally
compensated by another aromatic ring occupying the vacated cavity.

Same Cu-anchored descriptor method as test_phe227_compensation.py:
each aromatic ring is described by (dmin, dmax) to the two Cu atoms.
We compare the aromatic environment of substituted structures against
controls (Phe65=F, Trp68=W) to see if the lost ring is replaced.

Usage:
    python test_aromatic_compensation.py /path/to/cifs

Input:
    - position_vectors.csv (same directory as this script)
    - CIF directory (first positional arg, or /mnt/models)

Creates accession lists automatically from position_vectors.csv.
"""
import glob, os, math, statistics, csv, sys
from collections import Counter
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'electrostatics'))
from extract_features import parse_cif, dist

VECTORS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'position_vectors.csv')

RING = {
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TRP': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
}


def cu_pair(atoms):
    cus = [(a['x'], a['y'], a['z']) for a in atoms if a['elem'] == 'CU']
    if len(cus) < 2:
        return None
    best = None
    for i in range(len(cus)):
        for j in range(i + 1, len(cus)):
            d = dist(cus[i], cus[j])
            if best is None or d < best[0]:
                best = (d, cus[i], cus[j])
    return best[1], best[2]


def aromatics(atoms):
    groups = {}
    for a in atoms:
        if a['res'] in RING and a['name'] in RING[a['res']]:
            groups.setdefault((a['chain'], a['seq'], a['res']), []).append(a)
    out = []
    for (ch, seq, res), ats in groups.items():
        if len(ats) < 4:
            continue
        c = (sum(x['x'] for x in ats) / len(ats),
             sum(x['y'] for x in ats) / len(ats),
             sum(x['z'] for x in ats) / len(ats))
        out.append((res, c))
    return out


def descr(path):
    atoms = parse_cif(path)
    cp = cu_pair(atoms)
    if not cp:
        return None
    cu1, cu2 = cp
    mid = tuple((cu1[k] + cu2[k]) / 2 for k in range(3))
    res = []
    for restype, c in aromatics(atoms):
        if dist(mid, c) <= 12:
            d1, d2 = dist(cu1, c), dist(cu2, c)
            res.append((restype, min(d1, d2), max(d1, d2)))
    return res


def count_near(desclist, slot, tol=1.6, types=('PHE', 'TYR', 'TRP')):
    s_dmin, s_dmax = slot
    hits = 0
    for acc, d in desclist:
        ok = any(r[0] in types and abs(r[1] - s_dmin) <= tol and abs(r[2] - s_dmax) <= tol
                 for r in d)
        if ok:
            hits += 1
    return 100 * hits / len(desclist) if desclist else 0


def near_counts(desclist, rmax=9.0):
    nphe = []; narom = []
    for acc, d in desclist:
        nphe.append(sum(1 for r in d if r[0] == 'PHE' and r[2] <= rmax))
        narom.append(sum(1 for r in d if r[2] <= rmax))
    return statistics.mean(nphe), statistics.mean(narom)


def build_groups(cifdir):
    with open(VECTORS) as f:
        rows = list(csv.DictReader(f))

    groups = {
        'Phe65_L': set(), 'Phe65_F': set(),
        'Trp68_A': set(), 'Trp68_W': set(),
    }
    for r in rows:
        acc = r['accession']
        if r['Phe65'] == 'L': groups['Phe65_L'].add(acc)
        if r['Phe65'] == 'F': groups['Phe65_F'].add(acc)
        if r['Trp68'] == 'A': groups['Trp68_A'].add(acc)
        if r['Trp68'] == 'W': groups['Trp68_W'].add(acc)

    cifs = glob.glob(os.path.join(cifdir, '*.cif'))
    acc_to_cif = {}
    for p in cifs:
        acc = os.path.basename(p).split('_taxID_')[0]
        acc_to_cif[acc] = p

    print(f"CIF files found: {len(cifs)}")
    for g, accs in groups.items():
        avail = sum(1 for a in accs if a in acc_to_cif)
        print(f"  {g}: {len(accs)} in vectors, {avail} with CIF")

    return groups, acc_to_cif


def parse_group(accs, acc_to_cif, max_ctrl=500):
    """Parse CIFs for a group. Cap control sets to max_ctrl for speed."""
    import random
    todo = [a for a in accs if a in acc_to_cif]
    if len(todo) > max_ctrl:
        random.seed(42)
        todo = random.sample(todo, max_ctrl)
    results = []
    for acc in todo:
        d = descr(acc_to_cif[acc])
        if d is not None:
            results.append((acc, d))
    return results


def analyse(name, test_data, ctrl_data, test_label, ctrl_label):
    print(f"\n{'=' * 60}")
    print(f"  {name}")
    print(f"{'=' * 60}")
    print(f"Parsed: {test_label}={len(test_data)}  {ctrl_label}={len(ctrl_data)}")

    if not test_data or not ctrl_data:
        print("  ** Not enough structures, skipping **")
        return

    # Find conserved aromatic slots from control set
    all_arom = [(r[1], r[2]) for _, d in ctrl_data for r in d]
    grid = Counter((round(x[0]), round(x[1])) for x in all_arom)
    print(f"\nTop aromatic slots in {ctrl_label} (dmin, dmax) [rounded]:")
    for slot, cnt in grid.most_common(6):
        print(f"  {slot}: {cnt} hits")

    print(f"\nPer-slot compensation: % of structures with an aromatic there")
    print(f"  {'slot':>12s}  PHE({ctrl_label})  PHE({test_label})  ANY({ctrl_label})  ANY({test_label})")
    for slot, _ in grid.most_common(6):
        fc = count_near(ctrl_data, slot, types=('PHE',))
        ft = count_near(test_data, slot, types=('PHE',))
        ac = count_near(ctrl_data, slot, types=('PHE', 'TYR', 'TRP'))
        at = count_near(test_data, slot, types=('PHE', 'TYR', 'TRP'))
        print(f"  {str(slot):>12s}    {fc:5.0f}%      {ft:5.0f}%      {ac:5.0f}%      {at:5.0f}%")

    fp_c, fa_c = near_counts(ctrl_data)
    fp_t, fa_t = near_counts(test_data)
    print(f"\nMean aromatic rings within 9 A of Cu midpoint:")
    print(f"  PHE:  {ctrl_label}={fp_c:.2f}  {test_label}={fp_t:.2f}  (diff={fp_t - fp_c:+.2f})")
    print(f"  ALL:  {ctrl_label}={fa_c:.2f}  {test_label}={fa_t:.2f}  (diff={fa_t - fa_c:+.2f})")


def main():
    cifdir = sys.argv[1] if len(sys.argv) > 1 else '/mnt/models'
    groups, acc_to_cif = build_groups(cifdir)

    print("\nParsing Phe65 groups...")
    phe65_L = parse_group(groups['Phe65_L'], acc_to_cif, max_ctrl=500)
    phe65_F = parse_group(groups['Phe65_F'], acc_to_cif, max_ctrl=500)
    analyse("Phe65: Leu vs Phe (control)", phe65_L, phe65_F, "Phe65=L", "Phe65=F")

    print("\nParsing Trp68 groups...")
    trp68_A = parse_group(groups['Trp68_A'], acc_to_cif, max_ctrl=500)
    trp68_W = parse_group(groups['Trp68_W'], acc_to_cif, max_ctrl=500)
    analyse("Trp68: Ala vs Trp (control)", trp68_A, trp68_W, "Trp68=A", "Trp68=W")


if __name__ == '__main__':
    main()
