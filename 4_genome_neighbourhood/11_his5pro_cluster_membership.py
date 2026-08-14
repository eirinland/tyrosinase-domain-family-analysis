"""
Step 11: Per-query BGC-context evidence table for the His->Pro group.

For each his5pro query with neighbourhood data, scores the fungal-RiPP/mycotoxin
biosynthetic-gene-cluster (BGC) context: number of UstYa/DUF3328 paralogs in the
window, closest UstYa distance, and presence of canonical BGC accessory Pfams
(pathway-specific TF, MFS/ABC transporter, P450, O-methyltransferase).

Output: his5pro_cluster_membership.tsv
"""

import csv
from collections import defaultdict
from pathlib import Path

WORK = Path(__file__).parent
STAGE3 = Path('/Users/eirinlandsem/Library/Mobile Documents/com~apple~CloudDocs/'
              'Proteinkjemi_PhD/Skriving/Thesis/manuscripts:articles/Manuscript_TyrY/'
              'New_bioinf/bioinf_redo/Super_reference_pipeline/3_noncanonical_analysis/'
              'helix_and_gap_filtered_structures.tsv')
OUT = WORK / 'his5pro_cluster_membership.tsv'

MARKERS = {
    'TF':          {'PF00172', 'PF04082', 'PF11951'},          # Zn2Cys6 / fungal-TF
    'transporter': {'PF07690', 'PF00083', 'PF00005', 'PF00664'},  # MFS / sugar / ABC
    'P450':        {'PF00067'},
    'OMT':         {'PF00891'},                                 # O-methyltransferase
}
USTYA = 'PF11807'


def main():
    pid2fam = {}
    for r in csv.DictReader(open(WORK / 'his5pro_pfam_map.tsv'), delimiter='\t'):
        pid2fam[r['protein_id'].split('.')[0]] = {
            f.split('|')[0] for f in r['pfam'].split(';') if f}

    his5 = {r['accession'] for r in csv.DictReader(open(WORK / 'target_accessions.tsv'),
            delimiter='\t') if r['subgroup'] == 'his5pro'}

    sp, gen = {}, {}
    for r in csv.DictReader(open(STAGE3), delimiter='\t'):
        sp[r['accession']] = r.get('species', '')
        gen[r['accession']] = r.get('genus', '')

    # gather neighbours per query
    nb = defaultdict(list)   # query -> list of (offset, pid)
    for r in csv.DictReader(open(WORK / 'neighbourhoods.tsv'), delimiter='\t'):
        if r['query_accession'] in his5 and r['is_target'] != '1' and r['protein_id']:
            nb[r['query_accession']].append((int(r['offset']),
                                             r['protein_id'].split('.')[0]))

    rows = []
    for q, neigh in nb.items():
        fams_all = set()
        ust_offsets = []
        for off, pid in neigh:
            fams = pid2fam.get(pid, set())
            fams_all |= fams
            if USTYA in fams:
                ust_offsets.append(abs(off))
        n_ust = len(ust_offsets)
        rows.append({
            'accession': q,
            'genus': gen.get(q, ''),
            'species': sp.get(q, ''),
            'n_neighbours': len(neigh),
            'n_ustya': n_ust,
            'closest_ustya_offset': min(ust_offsets) if ust_offsets else '',
            'has_TF': int(bool(fams_all & MARKERS['TF'])),
            'has_transporter': int(bool(fams_all & MARKERS['transporter'])),
            'has_P450': int(bool(fams_all & MARKERS['P450'])),
            'has_OMT': int(bool(fams_all & MARKERS['OMT'])),
            'bgc_context': '',
        })
    for r in rows:
        markers = r['has_TF'] + r['has_transporter'] + r['has_P450'] + r['has_OMT']
        r['bgc_context'] = ('UstYa_cluster' if r['n_ustya'] >= 1
                            else ('BGC_accessory' if markers >= 2 else 'none'))

    rows.sort(key=lambda r: (-r['n_ustya'],
                             r['closest_ustya_offset'] if r['closest_ustya_offset'] != '' else 99))
    fields = ['accession', 'genus', 'species', 'n_neighbours', 'n_ustya',
              'closest_ustya_offset', 'has_TF', 'has_transporter', 'has_P450',
              'has_OMT', 'bgc_context']
    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    n = len(rows)
    n_ust = sum(1 for r in rows if r['n_ustya'] >= 1)
    n_bgc = sum(1 for r in rows if r['bgc_context'] != 'none')
    print(f'his5pro queries with neighbourhood data: {n}')
    print(f'  in a UstYa cluster (>=1 UstYa neighbour): {n_ust} ({n_ust/178*100:.1f}% of 178)')
    print(f'  any BGC context (UstYa or >=2 accessory): {n_bgc} ({n_bgc/178*100:.1f}% of 178)')
    print(f'  no BGC context: {n - n_bgc}')
    print(f'Output: {OUT.name}')


if __name__ == '__main__':
    main()
