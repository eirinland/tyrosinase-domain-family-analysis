"""
Step 10: Map his5pro neighbour proteins (GenBank CDS IDs) -> Pfam families.

GenBank CDS protein IDs -> UniProtKB (UniProt ID-mapping service) retrieving the
Pfam cross-references (xref_pfam), then resolve Pfam accessions to names via the
InterPro API. Output feeds the Pfam rollup in 9_his5pro_cooccurrence.py.

Output: his5pro_pfam_map.tsv  (protein_id <tab> PFxxxxx|Name;PFyyyyy|Name)
        his5pro_pfam_names.json (cache of pfam_id -> name)
"""

import csv
import json
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.parse import urlencode

WORK = Path(__file__).parent
NEIGHBOURHOODS = WORK / 'neighbourhoods.tsv'
TARGETS = WORK / 'target_accessions.tsv'
OUT = WORK / 'his5pro_pfam_map.tsv'
NAME_CACHE = WORK / 'his5pro_pfam_names.json'

API = 'https://rest.uniprot.org'
FROM_DB = 'EMBL-GenBank-DDBJ_CDS'


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
    """Stream TSV results: columns From, Entry, Pfam. Follow Link rel=next."""
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
    his5 = {r['accession'] for r in csv.DictReader(open(TARGETS), delimiter='\t')
            if r['subgroup'] == 'his5pro'}
    pids = sorted({r['protein_id'].split('.')[0]
                   for r in csv.DictReader(open(NEIGHBOURHOODS), delimiter='\t')
                   if r['query_accession'] in his5 and r['is_target'] != '1'
                   and r['protein_id']})
    print(f'{len(pids)} unique neighbour protein IDs to map')

    # ID-map in one job (limit is 100k)
    job = submit(pids)
    print(f'  idmapping job {job} submitted, waiting...')
    wait(job)
    rows = fetch_results(job)
    print(f'  {len(rows)} UniProt hits returned')

    # protein_id -> set of pfam accessions
    pid2pfam = defaultdict(set)
    for frm, entry, pfam in rows:
        for fam in pfam.split(';'):
            fam = fam.strip()
            if fam.startswith('PF'):
                pid2pfam[frm.split('.')[0]].add(fam)

    mapped = {k: v for k, v in pid2pfam.items() if v}
    print(f'  {len(mapped)}/{len(pids)} protein IDs have >=1 Pfam '
          f'({len(mapped)/len(pids)*100:.0f}% coverage)')

    # resolve names
    cache = json.loads(NAME_CACHE.read_text()) if NAME_CACHE.exists() else {}
    all_fams = sorted({f for fams in mapped.values() for f in fams})
    print(f'  resolving {len(all_fams)} unique Pfam names...')
    for f in all_fams:
        pfam_name(f, cache)
    NAME_CACHE.write_text(json.dumps(cache))

    with open(OUT, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['protein_id', 'pfam'])
        for pid in pids:
            fams = sorted(mapped.get(pid, []))
            val = ';'.join(f'{f}|{cache.get(f, "")}' for f in fams)
            w.writerow([pid, val])
    print(f'Output: {OUT.name} ({len(mapped)} annotated of {len(pids)})')


if __name__ == '__main__':
    main()
