"""
Map H6Gln neighbour proteins (GenBank CDS IDs) -> Pfam via UniProt ID-mapping.
Then re-run co-occurrence + cluster membership scoring.

Uses the same approach as 10_pfam_via_uniprot.py (idmapping service).
"""

import csv
import json
import re
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.parse import urlencode

WORK = Path(__file__).parent
NB_FILE = WORK / 'h6gln_neighbourhoods.tsv'
OUT_PFAM = WORK / 'h6gln_pfam_map.tsv'
OUT_NAMES = WORK / 'h6gln_pfam_names.json'
OUT_COOC_PFAM = WORK / 'h6gln_cooccurrence_pfam.tsv'
OUT_CLUSTER = WORK / 'h6gln_cluster_membership.tsv'

STAGE3 = WORK.parent / '3_noncanonical_analysis' / 'helix_and_gap_filtered_structures.tsv'

API = 'https://rest.uniprot.org'
FROM_DB = 'EMBL-GenBank-DDBJ_CDS'

MARKERS = {
    'TF':          {'PF00172', 'PF04082', 'PF11951'},
    'transporter': {'PF07690', 'PF00083', 'PF00005', 'PF00664'},
    'P450':        {'PF00067'},
    'OMT':         {'PF00891'},
}
USTYA = 'PF11807'


def submit(ids):
    data = urlencode({'from': FROM_DB, 'to': 'UniProtKB', 'ids': ','.join(ids)}).encode()
    req = Request(f'{API}/idmapping/run', data=data)
    return json.loads(urlopen(req, timeout=120).read())['jobId']


def wait(job):
    url = f'{API}/idmapping/status/{job}'
    for _ in range(120):
        req = Request(url, headers={'Accept': 'application/json'})
        r = json.loads(urlopen(req, timeout=60).read())
        st = r.get('jobStatus')
        if st in ('FINISHED', None) or 'results' in r:
            return
        if st in ('RUNNING', 'NEW'):
            time.sleep(3)
            continue
        raise RuntimeError(f'idmapping job {job}: {r}')
    raise TimeoutError(f'idmapping job {job} did not finish')


def fetch_results(job):
    url = (f'{API}/idmapping/uniprotkb/results/{job}'
           f'?fields=accession,xref_pfam&format=tsv&size=500')
    rows = []
    while url:
        req = Request(url, headers={'Accept': 'text/plain'})
        with urlopen(req, timeout=120) as resp:
            text = resp.read().decode()
            link = resp.headers.get('Link', '')
        lines = text.strip().split('\n')
        for ln in lines[1:]:
            parts = ln.split('\t')
            if len(parts) >= 3:
                rows.append((parts[0], parts[1], parts[2]))
            elif len(parts) == 2:
                rows.append((parts[0], parts[1], ''))
        m = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = m.group(1) if m else ''
        time.sleep(0.2)
    return rows


def pfam_name(pid, cache):
    if pid in cache:
        return cache[pid]
    url = f'https://www.ebi.ac.uk/interpro/api/entry/pfam/{pid}'
    try:
        req = Request(url, headers={'Accept': 'application/json'})
        with urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read())
        nm = data.get('metadata', {}).get('name', {})
        nm = nm.get('name', '') if isinstance(nm, dict) else str(nm)
    except Exception:
        nm = ''
    cache[pid] = nm
    time.sleep(0.15)
    return nm


def main():
    # Collect neighbour protein IDs
    nb_rows = list(csv.DictReader(open(NB_FILE), delimiter='\t'))
    query_accs = {r['query_accession'] for r in nb_rows}
    neigh = [r for r in nb_rows if r['is_target'] != '1' and r['protein_id']]
    pids = sorted({r['protein_id'].split('.')[0] for r in neigh})
    print(f'{len(pids)} unique neighbour protein IDs to map')

    # ID-mapping
    job = submit(pids)
    print(f'  idmapping job {job} submitted, waiting...')
    wait(job)
    rows = fetch_results(job)
    print(f'  {len(rows)} UniProt hits returned')

    pid2pfam = defaultdict(set)
    for frm, entry, pfam_col in rows:
        for fam in pfam_col.split(';'):
            fam = fam.strip()
            if fam.startswith('PF'):
                pid2pfam[frm.split('.')[0]].add(fam)

    mapped = {k: v for k, v in pid2pfam.items() if v}
    print(f'  {len(mapped)}/{len(pids)} protein IDs have >=1 Pfam')

    # Resolve names
    cache = json.loads(OUT_NAMES.read_text()) if OUT_NAMES.exists() else {}
    all_fams = sorted({f for fams in mapped.values() for f in fams})
    print(f'  resolving {len(all_fams)} unique Pfam names...')
    for f in all_fams:
        pfam_name(f, cache)
    OUT_NAMES.write_text(json.dumps(cache, indent=1))

    # Write pfam map
    with open(OUT_PFAM, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['protein_id', 'pfam'])
        for pid in pids:
            fams = sorted(mapped.get(pid, []))
            val = ';'.join(f'{f}|{cache.get(f, "")}' for f in fams)
            w.writerow([pid, val])
    print(f'  -> {OUT_PFAM.name} ({len(mapped)} annotated)')

    # --- Pfam co-occurrence ---
    pid2fam_list = {}
    for pid, fams in mapped.items():
        pid2fam_list[pid] = [f'{f}|{cache.get(f,"")}' for f in sorted(fams)]

    per_query_pfam = defaultdict(set)
    pfam_total = Counter()
    pfam_name_map = {}
    for r in neigh:
        pid = r['protein_id'].split('.')[0]
        for fam_str in pid2fam_list.get(pid, []):
            fid, _, nm = fam_str.partition('|')
            per_query_pfam[r['query_accession']].add(fid)
            pfam_total[fid] += 1
            if nm:
                pfam_name_map[fid] = nm

    pfam_qcount = Counter()
    for q, fams in per_query_pfam.items():
        for f in fams:
            pfam_qcount[f] += 1

    n_total = len(query_accs)
    n_with_pfam = len(per_query_pfam)
    with open(OUT_COOC_PFAM, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['pfam_id', 'pfam_name', 'n_queries', 'freq_of_all_queries',
                     'freq_of_queries_with_pfam', 'total_occurrences'])
        for fid, qc in pfam_qcount.most_common():
            w.writerow([fid, pfam_name_map.get(fid, ''), qc, f'{qc/n_total:.3f}',
                        f'{qc/n_with_pfam:.3f}' if n_with_pfam else '0', pfam_total[fid]])
    print(f'  -> {OUT_COOC_PFAM.name}')

    # --- BGC cluster membership ---
    sp_map, gen_map = {}, {}
    with open(STAGE3) as f:
        raw = f.read().replace('\r', '')
    for r in csv.DictReader(raw.splitlines(), delimiter='\t'):
        sp_map[r['accession']] = r.get('species', '')
        gen_map[r['accession']] = r.get('genus', '')

    nb_by_q = defaultdict(list)
    for r in neigh:
        nb_by_q[r['query_accession']].append((int(r['offset']), r['protein_id'].split('.')[0]))

    clust_rows = []
    for q in sorted(query_accs):
        fams_all = set()
        ust_offsets = []
        for off, pid in nb_by_q.get(q, []):
            fams = pid2pfam.get(pid, set())
            fams_all |= fams
            if USTYA in fams:
                ust_offsets.append(abs(off))
        n_ust = len(ust_offsets)
        markers = sum([
            bool(fams_all & MARKERS['TF']),
            bool(fams_all & MARKERS['transporter']),
            bool(fams_all & MARKERS['P450']),
            bool(fams_all & MARKERS['OMT']),
        ])
        bgc = 'UstYa_cluster' if n_ust >= 1 else ('BGC_accessory' if markers >= 2 else 'none')
        clust_rows.append({
            'accession': q, 'genus': gen_map.get(q, ''), 'species': sp_map.get(q, ''),
            'n_neighbours': len(nb_by_q.get(q, [])),
            'n_ustya': n_ust, 'closest_ustya_offset': min(ust_offsets) if ust_offsets else '',
            'has_TF': int(bool(fams_all & MARKERS['TF'])),
            'has_transporter': int(bool(fams_all & MARKERS['transporter'])),
            'has_P450': int(bool(fams_all & MARKERS['P450'])),
            'has_OMT': int(bool(fams_all & MARKERS['OMT'])),
            'bgc_context': bgc,
        })

    clust_rows.sort(key=lambda r: (-r['n_ustya'],
                                    r['closest_ustya_offset'] if r['closest_ustya_offset'] != '' else 99))
    fields = ['accession', 'genus', 'species', 'n_neighbours', 'n_ustya',
              'closest_ustya_offset', 'has_TF', 'has_transporter', 'has_P450', 'has_OMT', 'bgc_context']
    with open(OUT_CLUSTER, 'w', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fields)
        w.writeheader()
        w.writerows(clust_rows)

    # Summary
    n_data = len(query_accs)
    n_ust = sum(1 for r in clust_rows if r['bgc_context'] == 'UstYa_cluster')
    n_bgc_acc = sum(1 for r in clust_rows if r['bgc_context'] == 'BGC_accessory')
    n_none = sum(1 for r in clust_rows if r['bgc_context'] == 'none')
    print(f'\n=== BGC Summary ===')
    print(f'Queries with data: {n_data}/{n_total}')
    print(f'  UstYa cluster: {n_ust}')
    print(f'  BGC accessory: {n_bgc_acc}')
    print(f'  no BGC context: {n_none}')
    print()
    print(f'Top Pfam co-occurrences:')
    for fid, qc in pfam_qcount.most_common(20):
        print(f'  {qc:3d} ({qc/n_total*100:4.1f}%)  {fid} {pfam_name_map.get(fid,"")}')


if __name__ == '__main__':
    main()
