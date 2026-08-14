"""
Step 5: Summarize genome neighbourhoods by group.

Reads neighbourhoods.tsv + target_accessions.tsv, produces:
  - Per-group product frequency tables (which gene products co-occur near PPOs)
  - Per-group full product string frequencies (complete gene names)
  - Per-group neighbourhood composition summaries

Output: summary_by_group.tsv, product_frequencies.tsv, full_products_by_group.tsv
"""

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

WORK = Path(__file__).parent
TARGETS = WORK / 'target_accessions.tsv'
NEIGHBOURHOODS = WORK / 'neighbourhoods.tsv'
SUMMARY_OUT = WORK / 'summary_by_group.tsv'
PRODUCTS_OUT = WORK / 'product_frequencies.tsv'
FULL_PRODUCTS_OUT = WORK / 'full_products_by_group.tsv'

STOP_WORDS = {
    'protein', 'hypothetical', 'putative', 'probable', 'uncharacterized',
    'predicted', 'domain-containing', 'like', 'family', 'containing',
    'related', 'with', 'type', 'that', 'this', 'from', 'have', 'been',
}


def normalize_product(product: str) -> str:
    p = product.strip().lower()
    p = re.sub(r'\s+', ' ', p)
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

    with open(NEIGHBOURHOODS) as f:
        nb_rows = list(csv.DictReader(f, delimiter='\t'))

    print(f'Loaded {len(nb_rows)} neighbourhood entries')

    group_products = defaultdict(Counter)
    group_full_products = defaultdict(Counter)
    group_counts = Counter()
    group_flanking = Counter()

    for r in nb_rows:
        qacc = r['query_accession']
        subgroups = acc_to_subgroups.get(qacc, {'unknown'})

        if r['is_target'] == '1':
            for g in subgroups:
                group_counts[g] += 1
            continue

        product = r.get('product', '')
        norm = normalize_product(product)

        for g in subgroups:
            group_flanking[g] += 1

            if norm:
                group_full_products[g][norm] += 1

            for word in product.lower().split():
                if len(word) > 3 and word not in STOP_WORDS:
                    group_products[g][word] += 1

    # Write summary per group
    with open(SUMMARY_OUT, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['subgroup', 'n_queries', 'n_flanking_genes',
                     'avg_flanking', 'top_keywords', 'top_products'])
        for g in sorted(group_counts.keys()):
            n_q = group_counts[g]
            n_f = group_flanking[g]
            avg = f'{n_f / n_q:.1f}' if n_q > 0 else '0'
            top_kw = '; '.join(f'{w}({c})' for w, c in group_products[g].most_common(15))
            top_prod = '; '.join(
                f'{p}({c})' for p, c in group_full_products[g].most_common(10))
            w.writerow([g, n_q, n_f, avg, top_kw, top_prod])

    # Write full product frequencies across all groups
    all_products = set()
    for counts in group_full_products.values():
        for prod, c in counts.items():
            if c >= 3:
                all_products.add(prod)

    groups = sorted(group_counts.keys())
    with open(FULL_PRODUCTS_OUT, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['product'] + groups)
        for prod in sorted(all_products):
            row = [prod] + [str(group_full_products[g].get(prod, 0)) for g in groups]
            w.writerow(row)

    # Write keyword frequencies across all groups
    all_keywords = set()
    for counts in group_products.values():
        for kw, c in counts.items():
            if c >= 5:
                all_keywords.add(kw)

    with open(PRODUCTS_OUT, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['keyword'] + groups)
        for kw in sorted(all_keywords):
            row = [kw] + [str(group_products[g].get(kw, 0)) for g in groups]
            w.writerow(row)

    print(f'\nGroups with data: {len(group_counts)}')
    for g in sorted(group_counts.keys()):
        print(f'  {g}: {group_counts[g]} queries, '
              f'{group_flanking[g]} flanking genes')
    print(f'\nUnique annotated products (count>=3): {len(all_products)}')
    print(f'Unique keywords (count>=5): {len(all_keywords)}')
    print(f'Output: {SUMMARY_OUT.name}, {PRODUCTS_OUT.name}, {FULL_PRODUCTS_OUT.name}')


if __name__ == '__main__':
    main()
