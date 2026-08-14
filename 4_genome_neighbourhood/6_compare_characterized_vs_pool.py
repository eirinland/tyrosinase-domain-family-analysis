"""
Step 6: Compare characterized PPO neighbourhoods to their candidate pool.

For each characterized protein with flanking data, checks whether its
neighbourhood resembles the aggregate pattern of its corresponding pool.

Produces:
  - Per-protein neighbourhood profile vs pool profile comparison
  - Jaccard similarity of product sets
  - Enrichment of pool-characteristic genes near characterized exemplars
"""

import csv
from collections import Counter, defaultdict
from pathlib import Path

WORK = Path(__file__).parent
TARGETS = WORK / 'target_accessions.tsv'
NEIGHBOURHOODS = WORK / 'neighbourhoods.tsv'

POOL_MAP = {
    'oMP': 'oMP',
    'oAPO': 'oAPO',
    'DCT': 'DCT_DHICA',
    'DHICA_ox': 'DCT_DHICA',
    'hemocyanin': 'hemocyanin',
    'TYR': None,
    'CaOx': None,
    'AUS': None,
    'biosynthetic': None,
}


def normalize(product: str) -> str:
    p = product.strip().lower()
    if p in ('hypothetical protein', 'uncharacterized protein', ''):
        return ''
    return p


def main():
    with open(TARGETS) as f:
        target_rows = list(csv.DictReader(f, delimiter='\t'))

    acc_to_subgroups = defaultdict(set)
    acc_to_groups = defaultdict(set)
    for r in target_rows:
        acc_to_subgroups[r['accession']].add(r['subgroup'])
        acc_to_groups[r['accession']].add(r['group'])

    characterized = {r['accession'] for r in target_rows if r['group'] == 'characterized'}
    pool_members = {r['accession'] for r in target_rows if r['group'] == 'candidate_pool'}

    with open(NEIGHBOURHOODS) as f:
        nb_rows = list(csv.DictReader(f, delimiter='\t'))

    per_query_products = defaultdict(list)
    for r in nb_rows:
        if r['is_target'] == '1':
            continue
        norm = normalize(r.get('product', ''))
        if norm:
            per_query_products[r['query_accession']].append(norm)

    pool_subgroups = {'oMP', 'oAPO', 'DCT_DHICA', 'hemocyanin'}

    pool_product_counts = defaultdict(Counter)
    pool_query_counts = Counter()
    for acc, products in per_query_products.items():
        for sg in acc_to_subgroups.get(acc, set()):
            if sg in pool_subgroups:
                pool_product_counts[sg].update(products)
                pool_query_counts[sg] += 1

    print('=' * 70)
    print('POOL AGGREGATE PROFILES')
    print('=' * 70)
    for pool in sorted(pool_subgroups):
        n = pool_query_counts[pool]
        top = pool_product_counts[pool].most_common(10)
        print(f'\n{pool} (n={n} queries with flanking data):')
        for prod, c in top:
            print(f'  {c:5d} ({c/n*100:4.1f}%)  {prod}')

    print('\n' + '=' * 70)
    print('CHARACTERIZED EXEMPLAR NEIGHBOURHOODS vs POOL PROFILES')
    print('=' * 70)

    for char_sg in ['oMP', 'oAPO', 'TYR', 'biosynthetic']:
        char_accs = [a for a in characterized
                     if char_sg in acc_to_subgroups.get(a, set())
                     and a in per_query_products]
        if not char_accs:
            continue

        target_pool = POOL_MAP.get(char_sg)
        print(f'\n--- Characterized {char_sg} ({len(char_accs)} with data) ---')
        if target_pool:
            print(f'    Corresponding pool: {target_pool} '
                  f'(n={pool_query_counts[target_pool]})')

        for acc in sorted(char_accs):
            products = per_query_products[acc]
            prod_set = set(products)
            prod_counts = Counter(products)

            print(f'\n  {acc} ({len(products)} flanking genes):')
            for prod, c in prod_counts.most_common(8):
                print(f'    {c}x  {prod}')

            if target_pool and pool_product_counts[target_pool]:
                pool_prods = set(pool_product_counts[target_pool].keys())
                pool_top20 = {p for p, _ in
                              pool_product_counts[target_pool].most_common(20)}

                shared_any = prod_set & pool_prods
                shared_top = prod_set & pool_top20
                jaccard = (len(prod_set & pool_prods) /
                           len(prod_set | pool_prods)) if prod_set | pool_prods else 0

                print(f'    Pool overlap: {len(shared_any)}/{len(prod_set)} '
                      f'products also in pool '
                      f'(Jaccard={jaccard:.3f})')
                print(f'    Matches pool top-20: {len(shared_top)} '
                      f'({", ".join(sorted(shared_top)[:5])}{"..." if len(shared_top) > 5 else ""})')

    print('\n' + '=' * 70)
    print('CROSS-POOL COMPARISON: Characterized oMP/oAPO vs ALL pools')
    print('=' * 70)

    for char_sg in ['oMP', 'oAPO']:
        char_accs = [a for a in characterized
                     if char_sg in acc_to_subgroups.get(a, set())
                     and a in per_query_products]
        if not char_accs:
            continue

        char_products = set()
        for acc in char_accs:
            char_products.update(per_query_products[acc])

        print(f'\n  Characterized {char_sg} '
              f'({len(char_accs)} proteins, {len(char_products)} unique products):')
        for pool in sorted(pool_subgroups):
            pool_prods = set(pool_product_counts[pool].keys())
            if not pool_prods:
                continue
            shared = char_products & pool_prods
            jaccard = (len(shared) / len(char_products | pool_prods)
                       if char_products | pool_prods else 0)
            print(f'    vs {pool:12s}: {len(shared):3d} shared products, '
                  f'Jaccard={jaccard:.3f}')


if __name__ == '__main__':
    main()
