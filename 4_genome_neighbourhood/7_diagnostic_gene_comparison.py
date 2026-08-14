"""
Step 7: Check if pool-diagnostic gene categories appear near characterized exemplars.

Instead of exact product matching, uses functional keyword categories
to compare characterized PPO neighbourhoods against pool signatures.
"""

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

WORK = Path(__file__).parent
TARGETS = WORK / 'target_accessions.tsv'
NEIGHBOURHOODS = WORK / 'neighbourhoods.tsv'

FUNCTIONAL_CATEGORIES = {
    'copper_protein': r'copper|di-copper|tyrosinase|laccase|multicopper|polyphenol.oxid|catechol.oxid',
    'copper_chaperone': r'copper.chap|chaperon.*copper|melc1|melc2|tyrosinase.co.?factor',
    'transporter_MFS': r'major.facilitator|mfs.transporter|mfs.general',
    'transporter_other': r'transporter|permease|efflux|abc.transporter',
    'cytochrome_P450': r'cytochrome.p450|p450.monooxygenase',
    'methyltransferase': r'methyltransferase',
    'HET_incompatibility': r'heterokaryon|incompatibility|het-',
    'transcription_factor': r'transcription|regulator|activator|repressor',
    'oxidoreductase': r'oxidoreductase|dehydrogenase|reductase|oxygenase|oxidase',
    'hydrolase': r'hydrolase|peptidase|protease|esterase',
    'kinase': r'kinase',
    'biosynthetic_cluster': r'polyketide|nonribosomal|nrps|synthase',
    'melanin_pathway': r'dopachrome|tyrp|trp-1|trp-2|melanin|pigment',
}


def categorize(product: str) -> set[str]:
    p = product.lower()
    cats = set()
    for cat, pattern in FUNCTIONAL_CATEGORIES.items():
        if re.search(pattern, p):
            cats.add(cat)
    return cats


def main():
    with open(TARGETS) as f:
        target_rows = list(csv.DictReader(f, delimiter='\t'))

    acc_info = {}
    for r in target_rows:
        acc = r['accession']
        if acc not in acc_info:
            acc_info[acc] = {'groups': set(), 'subgroups': set()}
        acc_info[acc]['groups'].add(r['group'])
        acc_info[acc]['subgroups'].add(r['subgroup'])

    with open(NEIGHBOURHOODS) as f:
        nb_rows = list(csv.DictReader(f, delimiter='\t'))

    per_query = defaultdict(lambda: {'products': [], 'categories': Counter()})
    for r in nb_rows:
        if r['is_target'] == '1':
            continue
        acc = r['query_accession']
        product = r.get('product', '')
        per_query[acc]['products'].append(product)
        for cat in categorize(product):
            per_query[acc]['categories'][cat] += 1

    pool_subgroups = {'oMP', 'oAPO', 'DCT_DHICA', 'hemocyanin'}
    pool_cats = defaultdict(Counter)
    pool_n = Counter()
    pool_has_cat = defaultdict(lambda: Counter())

    for acc, data in per_query.items():
        for sg in acc_info.get(acc, {}).get('subgroups', set()):
            if sg in pool_subgroups:
                pool_n[sg] += 1
                for cat, count in data['categories'].items():
                    pool_cats[sg][cat] += count
                    pool_has_cat[sg][cat] += 1

    print('=' * 75)
    print('POOL FUNCTIONAL CATEGORY PROFILES (% of queries with ≥1 gene in category)')
    print('=' * 75)
    for pool in sorted(pool_subgroups):
        n = pool_n[pool]
        if n == 0:
            continue
        print(f'\n  {pool} (n={n}):')
        for cat in sorted(FUNCTIONAL_CATEGORIES.keys()):
            has = pool_has_cat[pool][cat]
            total = pool_cats[pool][cat]
            if has > 0:
                print(f'    {cat:30s}  {has:4d}/{n} queries ({has/n*100:5.1f}%)  '
                      f'{total} total genes')

    print('\n' + '=' * 75)
    print('CHARACTERIZED EXEMPLARS: functional categories in their neighbourhoods')
    print('=' * 75)

    for char_sg in ['oMP', 'oAPO', 'TYR', 'biosynthetic']:
        char_accs = [a for a in per_query
                     if 'characterized' in acc_info.get(a, {}).get('groups', set())
                     and char_sg in acc_info.get(a, {}).get('subgroups', set())]
        if not char_accs:
            continue

        print(f'\n--- Characterized {char_sg} ({len(char_accs)} with flanking data) ---')

        for acc in sorted(char_accs):
            data = per_query[acc]
            n_flanking = len(data['products'])
            cats = data['categories']
            print(f'\n  {acc} ({n_flanking} flanking genes):')
            if not cats:
                non_hyp = [p for p in data['products']
                           if 'hypothetical' not in p.lower()
                           and 'unknown' not in p.lower()
                           and p.strip()]
                if non_hyp:
                    print(f'    No category matches. Products: '
                          f'{"; ".join(non_hyp[:5])}')
                else:
                    print(f'    No informative products')
                continue
            for cat in sorted(cats.keys()):
                target_pool = {'oMP': 'oMP', 'oAPO': 'oAPO'}.get(char_sg)
                pool_pct = ''
                if target_pool and pool_n[target_pool] > 0:
                    pct = pool_has_cat[target_pool][cat] / pool_n[target_pool] * 100
                    pool_pct = f'  (pool: {pct:.1f}%)'
                print(f'    {cat:30s}  {cats[cat]:2d} genes{pool_pct}')

    print('\n' + '=' * 75)
    print('SIGNATURE MATCH SUMMARY')
    print('=' * 75)

    oMP_signature = {'transporter_MFS', 'HET_incompatibility', 'cytochrome_P450',
                     'methyltransferase'}
    oAPO_signature = {'copper_protein'}
    DCT_signature = {'melanin_pathway', 'copper_protein'}

    for label, sig, pool in [
        ('oMP BGC', oMP_signature, 'oMP'),
        ('oAPO copper-tandem', oAPO_signature, 'oAPO'),
    ]:
        char_accs = [a for a in per_query
                     if 'characterized' in acc_info.get(a, {}).get('groups', set())
                     and pool in acc_info.get(a, {}).get('subgroups', set())]
        pool_accs = [a for a in per_query
                     if pool in acc_info.get(a, {}).get('subgroups', set())
                     and 'candidate_pool' in acc_info.get(a, {}).get('groups', set())]

        if not char_accs and not pool_accs:
            continue

        print(f'\n  {label} signature ({", ".join(sorted(sig))}):')

        for acc in char_accs:
            cats = set(per_query[acc]['categories'].keys())
            matches = cats & sig
            print(f'    Characterized {acc}: '
                  f'{len(matches)}/{len(sig)} signature categories '
                  f'({", ".join(sorted(matches)) if matches else "none"})')

        if pool_accs:
            n_with_any = sum(1 for a in pool_accs
                             if set(per_query[a]['categories'].keys()) & sig)
            n_with_all = sum(1 for a in pool_accs
                             if sig <= set(per_query[a]['categories'].keys()))
            print(f'    Pool {pool} (n={len(pool_accs)}): '
                  f'{n_with_any} ({n_with_any/len(pool_accs)*100:.1f}%) have ≥1 sig. category, '
                  f'{n_with_all} ({n_with_all/len(pool_accs)*100:.1f}%) have all')


if __name__ == '__main__':
    main()
