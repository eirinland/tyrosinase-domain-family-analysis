"""
Step 9: EFI-GNT-style neighbour co-occurrence for the His->Pro mononuclear group.

EFI-GNT defines a neighbour family's "co-occurrence frequency" as the FRACTION OF
QUERY sequences in the cluster that have >=1 neighbour of that family (not the raw
neighbour count). This script computes that per-query fraction for the his5pro
group, over GenBank product strings (family proxy) and over keyword tokens.

Inputs : target_accessions.tsv (subgroup==his5pro), neighbourhoods.tsv
Outputs: his5pro_cooccurrence_products.tsv, his5pro_cooccurrence_keywords.tsv
Pfam    : if his5pro_pfam_map.tsv exists (protein_id -> pfam_id;pfam_name from
          10_pfam_via_uniprot.py) it is also rolled up -> his5pro_cooccurrence_pfam.tsv
"""

import csv
import re
from collections import Counter, defaultdict
from pathlib import Path

WORK = Path(__file__).parent
TARGETS = WORK / 'target_accessions.tsv'
NEIGHBOURHOODS = WORK / 'neighbourhoods.tsv'
PFAM_MAP = WORK / 'his5pro_pfam_map.tsv'
OUT_PROD = WORK / 'his5pro_cooccurrence_products.tsv'
OUT_KW = WORK / 'his5pro_cooccurrence_keywords.tsv'
OUT_PFAM = WORK / 'his5pro_cooccurrence_pfam.tsv'

SUBGROUP = 'his5pro'

STOP_WORDS = {
    'protein', 'hypothetical', 'putative', 'probable', 'uncharacterized',
    'predicted', 'domain-containing', 'like', 'family', 'containing',
    'related', 'with', 'type', 'that', 'this', 'from', 'have', 'been',
    'subunit', 'partial', 'unnamed', 'product', 'and', 'the', 'gene',
}


def norm_product(p: str) -> str:
    p = re.sub(r'\s+', ' ', p.strip().lower())
    if p in ('hypothetical protein', 'uncharacterized protein',
             'unnamed protein product', 'putative protein', ''):
        return ''
    return p


def write_table(path, header, rows):
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(header)
        w.writerows(rows)


def main():
    # his5pro query set
    his5 = {r['accession'] for r in csv.DictReader(open(TARGETS), delimiter='\t')
            if r['subgroup'] == SUBGROUP}
    n_total = len(his5)

    nb = [r for r in csv.DictReader(open(NEIGHBOURHOODS), delimiter='\t')
          if r['query_accession'] in his5]

    # queries that actually got fetched (>=1 row, incl. their own target row)
    queries_with_rows = {r['query_accession'] for r in nb}
    # neighbour rows (exclude the query's own CDS)
    neigh = [r for r in nb if r['is_target'] != '1']
    queries_with_neighbour = {r['query_accession'] for r in neigh}

    # ---- product-level co-occurrence (per-query presence) ----
    per_query_prod = defaultdict(set)        # query -> {normalized product}
    prod_total = Counter()                   # raw neighbour occurrences
    prod_min_dist = {}                       # product -> closest |offset| seen
    for r in neigh:
        np_ = norm_product(r.get('product', ''))
        if not np_:
            continue
        q = r['query_accession']
        per_query_prod[q].add(np_)
        prod_total[np_] += 1
        d = abs(int(r['offset']))
        prod_min_dist[np_] = min(prod_min_dist.get(np_, 99), d)

    prod_qcount = Counter()                  # product -> #distinct queries
    for q, prods in per_query_prod.items():
        for p in prods:
            prod_qcount[p] += 1

    n_with_data = len(queries_with_rows)
    n_with_neigh = len(queries_with_neighbour)
    n_with_info = len(per_query_prod)        # >=1 informative (non-hyp) neighbour

    prod_rows = []
    for p, qc in prod_qcount.most_common():
        prod_rows.append([
            p, qc, f'{qc / n_total:.3f}', f'{qc / n_with_info:.3f}',
            prod_total[p], prod_min_dist[p],
        ])
    write_table(OUT_PROD,
                ['product', 'n_queries', 'freq_of_all_queries',
                 'freq_of_queries_with_neighbours', 'total_occurrences',
                 'closest_offset'],
                prod_rows)

    # ---- keyword-level co-occurrence (per-query presence of a token) ----
    per_query_kw = defaultdict(set)
    kw_total = Counter()
    for r in neigh:
        prod = r.get('product', '').lower()
        if norm_product(prod) == '':
            continue
        q = r['query_accession']
        toks = {w for w in re.split(r'[^a-z0-9]+', prod)
                if len(w) > 3 and w not in STOP_WORDS}
        for w in toks:
            per_query_kw[q].add(w)
        for w in toks:
            kw_total[w] += 1
    kw_qcount = Counter()
    for q, kws in per_query_kw.items():
        for w in kws:
            kw_qcount[w] += 1
    kw_rows = [[w, qc, f'{qc / n_total:.3f}', kw_total[w]]
               for w, qc in kw_qcount.most_common() if qc >= 2]
    write_table(OUT_KW,
                ['keyword', 'n_queries', 'freq_of_all_queries', 'total_occurrences'],
                kw_rows)

    # ---- optional Pfam rollup ----
    pfam_note = 'no his5pro_pfam_map.tsv (run 10_pfam_via_uniprot.py first)'
    if PFAM_MAP.exists():
        pid2pfam = {}
        for r in csv.DictReader(open(PFAM_MAP), delimiter='\t'):
            fams = [f for f in r.get('pfam', '').split(';') if f]
            pid2pfam[r['protein_id'].split('.')[0]] = fams
        per_query_pfam = defaultdict(set)
        pfam_total = Counter()
        pfam_name = {}
        for r in neigh:
            pid = r.get('protein_id', '').split('.')[0]
            for fam in pid2pfam.get(pid, []):
                pid_, _, nm = fam.partition('|')
                per_query_pfam[r['query_accession']].add(pid_)
                pfam_total[pid_] += 1
                if nm:
                    pfam_name[pid_] = nm
        pfam_qcount = Counter()
        for q, fams in per_query_pfam.items():
            for fam in fams:
                pfam_qcount[fam] += 1
        n_with_pfam = len(per_query_pfam)
        pfam_rows = [[fam, pfam_name.get(fam, ''), qc,
                      f'{qc / n_total:.3f}',
                      f'{qc / n_with_pfam:.3f}' if n_with_pfam else '0',
                      pfam_total[fam]]
                     for fam, qc in pfam_qcount.most_common()]
        write_table(OUT_PFAM,
                    ['pfam_id', 'pfam_name', 'n_queries', 'freq_of_all_queries',
                     'freq_of_queries_with_pfam', 'total_occurrences'],
                    pfam_rows)
        pfam_note = (f'{len(pfam_rows)} Pfam families '
                     f'({n_with_pfam} queries with >=1 Pfam-annotated neighbour) '
                     f'-> {OUT_PFAM.name}')

    # ---- console summary ----
    print(f'his5pro queries (target list)          : {n_total}')
    print(f'  with fetched neighbourhood record    : {n_with_data}')
    print(f'  with >=1 neighbour gene              : {n_with_neigh}')
    print(f'  with >=1 informative (non-hyp) neigh : {n_with_info}')
    print(f'distinct informative products          : {len(prod_qcount)}')
    print()
    print(f'TOP co-occurring products (by # queries / {n_total}):')
    for p, qc, fa, fi, tot, d in prod_rows[:25]:
        print(f'  {qc:3d} ({float(fa)*100:4.1f}%)  off>={d}  {p[:62]}')
    print()
    print(f'TOP co-occurring keywords (by # queries / {n_total}):')
    for w, qc, fa, tot in kw_rows[:20]:
        print(f'  {qc:3d} ({float(fa)*100:4.1f}%)  {w}')
    print()
    print('Pfam:', pfam_note)
    print(f'Outputs: {OUT_PROD.name}, {OUT_KW.name}')


if __name__ == '__main__':
    main()
