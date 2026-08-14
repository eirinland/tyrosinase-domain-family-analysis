"""Metal3D confidence vs chemical-plausibility classification."""
import csv, glob, os, statistics
from collections import defaultdict

BASEDIR = os.path.dirname(os.path.abspath(__file__))
PLAUSIBLE = {'GLU', 'ASP', 'CYS', 'TYR', 'MET'}

def tier(h1, h2, h3):
    res = [h1, h2, h3]
    n = res.count('HIS')
    if n == 3: return 'canonical'
    if n == 2 and [r for r in res if r != 'HIS'][0] in PLAUSIBLE: return 'plausible'
    return 'divergent'

def overall(a, b):
    if a != 'divergent' and b != 'divergent': return 'binuclear'
    if a != 'divergent' or b != 'divergent': return 'mononuclear'
    return 'no_cu'

# Load classified structures
structs = {}
with open(os.path.join(BASEDIR, 'noncanonical_analysis.tsv')) as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if any(r[p] == '---' for p in ['CuA_His1','CuA_His2','CuA_His3','CuB_His1','CuB_His2','CuB_His3']):
            continue
        acc = r['accession'].split('_taxID_')[0]
        cua_t = tier(r['CuA_His1'], r['CuA_His2'], r['CuA_His3'])
        cub_t = tier(r['CuB_His1'], r['CuB_His2'], r['CuB_His3'])
        structs[acc] = {
            'cls': overall(cua_t, cub_t),
            'cua_tier': cua_t, 'cub_tier': cub_t,
            'cu1_assign': r.get('cu1_assignment', ''),
            'cu2_assign': r.get('cu2_assignment', ''),
        }

# Load Metal3D main format (per AF3 Cu atom)
main_by_acc = defaultdict(dict)
for fn in sorted(glob.glob(os.path.join(BASEDIR, 'metal3d/results/metal3d_*.tsv'))):
    with open(fn) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r.get('status') != 'ok': continue
            acc = r['accession'].split('_taxID_')[0]
            idx = r.get('af3_cu_index', '')
            main_by_acc[acc][idx] = {
                'dist': float(r['closest_cu_dist']) if r.get('closest_cu_dist') else None,
                'prob': float(r['closest_cu_prob']) if r.get('closest_cu_prob') else None,
            }

# Load Metal3D fallback format (per-site CuA/CuB)
fallback = {}
for fn in sorted(glob.glob(os.path.join(BASEDIR, 'metal3d/results/metal3d_neither_*.tsv'))) +           [os.path.join(BASEDIR, 'metal3d/results/combined_mono.tsv'),
           os.path.join(BASEDIR, 'metal3d/results/combined_all.tsv')]:
    if not os.path.exists(fn): continue
    with open(fn) as f:
        reader = csv.DictReader(f, delimiter='\t')
        if 'CuA_closest_cu_prob' not in (reader.fieldnames or []): continue
        for r in reader:
            if r.get('status') != 'ok': continue
            acc = r['accession'].split('_taxID_')[0]
            fallback[acc] = {
                'CuA_prob': float(r['CuA_closest_cu_prob']) if r.get('CuA_closest_cu_prob') else None,
                'CuA_dist': float(r['CuA_closest_cu_dist']) if r.get('CuA_closest_cu_dist') else None,
                'CuB_prob': float(r['CuB_closest_cu_prob']) if r.get('CuB_closest_cu_prob') else None,
                'CuB_dist': float(r['CuB_closest_cu_dist']) if r.get('CuB_closest_cu_dist') else None,
            }

def get_site_probs(acc):
    """Return (CuA_prob, CuB_prob) from Metal3D, or (None, None)."""
    s = structs[acc]
    cua_prob = cub_prob = None

    if acc in main_by_acc:
        for idx, vals in main_by_acc[acc].items():
            assign_key = f'cu{idx}_assign'
            site = s.get(assign_key, '')
            if site == 'CuA' and vals['prob'] is not None:
                cua_prob = vals['prob']
            elif site == 'CuB' and vals['prob'] is not None:
                cub_prob = vals['prob']

    if acc in fallback:
        fb = fallback[acc]
        if cua_prob is None and fb['CuA_prob'] is not None:
            cua_prob = fb['CuA_prob']
        if cub_prob is None and fb['CuB_prob'] is not None:
            cub_prob = fb['CuB_prob']

    return cua_prob, cub_prob

# Collect per-site probabilities grouped by tier and overall class
site_data = defaultdict(list)  # (site, tier) -> [probs]
class_data = defaultdict(lambda: {'cua': [], 'cub': []})

n_with_data = 0
for acc, s in structs.items():
    cua_p, cub_p = get_site_probs(acc)
    if cua_p is None and cub_p is None:
        continue
    n_with_data += 1

    if cua_p is not None:
        site_data[('CuA', s['cua_tier'])].append(cua_p)
        class_data[s['cls']]['cua'].append(cua_p)
    if cub_p is not None:
        site_data[('CuB', s['cub_tier'])].append(cub_p)
        class_data[s['cls']]['cub'].append(cub_p)

def stats(vals):
    if not vals: return 'n/a'
    m = statistics.median(vals)
    mu = statistics.mean(vals)
    q1 = sorted(vals)[len(vals)//4]
    q3 = sorted(vals)[3*len(vals)//4]
    hi = sum(1 for v in vals if v >= 0.5)
    return f'n={len(vals):>4}  mean={mu:.3f}  median={m:.3f}  Q1={q1:.3f}  Q3={q3:.3f}  >=0.5: {hi:>4} ({100*hi/len(vals):.0f}%)'

print(f'Structures with Metal3D data: {n_with_data}/{len(structs)}')
print()
print('=' * 80)
print(' PER-SITE TIER: Metal3D Cu probability at canonical His positions')
print('=' * 80)
for site in ['CuA', 'CuB']:
    print(f'\n  {site}:')
    for t in ['canonical', 'plausible', 'divergent']:
        vals = site_data.get((site, t), [])
        print(f'    {t:12s}  {stats(vals)}')

print()
print('=' * 80)
print(' PER-CLASS: Metal3D Cu probability by overall classification')
print('=' * 80)
for cls in ['binuclear', 'mononuclear', 'no_cu']:
    d = class_data[cls]
    print(f'\n  {cls}:')
    print(f'    CuA  {stats(d["cua"])}')
    print(f'    CuB  {stats(d["cub"])}')

# For mononuclear: split into good-site vs divergent-site
print()
print('=' * 80)
print(' MONONUCLEAR detail: good site vs divergent site')
print('=' * 80)
good_probs = []
div_probs = []
for acc, s in structs.items():
    if s['cls'] != 'mononuclear': continue
    cua_p, cub_p = get_site_probs(acc)
    if s['cua_tier'] != 'divergent' and cua_p is not None:
        good_probs.append(cua_p)
    elif s['cua_tier'] == 'divergent' and cua_p is not None:
        div_probs.append(cua_p)
    if s['cub_tier'] != 'divergent' and cub_p is not None:
        good_probs.append(cub_p)
    elif s['cub_tier'] == 'divergent' and cub_p is not None:
        div_probs.append(cub_p)
print(f'  Good site (canonical/plausible):  {stats(good_probs)}')
print(f'  Divergent site:                   {stats(div_probs)}')

# Probability distribution bins
print()
print('=' * 80)
print(' PROBABILITY DISTRIBUTION by classification')
print('=' * 80)
bins = [(0, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5),
        (0.5, 0.6), (0.6, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 1.01)]
header = f'{cls:14s} {site:5s}'
for lo, hi in bins:
    header += f' {lo:.1f}-{hi:.1f}'
print(header)
print('-' * len(header))

for cls in ['binuclear', 'mononuclear', 'no_cu']:
    for site_label in ['CuA', 'CuB']:
        vals = class_data[cls][site_label.lower()[:3]]
        row = f'{cls:14s} {site_label:5s}'
        for lo, hi in bins:
            c = sum(1 for v in vals if lo <= v < hi)
            row += f' {c:>7d}'
        print(row)
