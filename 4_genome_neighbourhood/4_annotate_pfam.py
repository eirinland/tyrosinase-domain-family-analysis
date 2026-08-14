"""
Step 4: Annotate flanking gene products with Pfam domains via InterPro API.

Reads neighbourhoods.tsv, collects unique protein IDs, queries InterPro
for domain annotations, and writes an enriched output.

Output: neighbourhoods_annotated.tsv (adds pfam_ids, pfam_names columns)
"""

import csv
import json
import sys
import time
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError

WORK = Path(__file__).parent
INPUT = WORK / 'neighbourhoods.tsv'
OUTPUT = WORK / 'neighbourhoods_annotated.tsv'
DOMAIN_CACHE = WORK / 'domain_cache.json'

INTERPRO_API = 'https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot'
BATCH_SIZE = 10


def load_cache() -> dict:
    if DOMAIN_CACHE.exists():
        return json.loads(DOMAIN_CACHE.read_text())
    return {}


def save_cache(cache: dict):
    DOMAIN_CACHE.write_text(json.dumps(cache))


def fetch_interpro_protein(protein_acc: str) -> list[dict]:
    """Fetch Pfam domain annotations for a single protein from InterPro."""
    # Try by UniProt accession first, then by EMBL protein ID
    url = f'https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{protein_acc}?format=json'
    try:
        req = Request(url, headers={'Accept': 'application/json'})
        with urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read())
        domains = []
        for result in data.get('results', []):
            meta = result.get('metadata', {})
            domains.append({
                'pfam_id': meta.get('accession', ''),
                'pfam_name': meta.get('name', ''),
            })
        return domains
    except (HTTPError, URLError):
        return []


def main():
    with open(INPUT) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))

    print(f'Loaded {len(rows)} neighbourhood entries')

    # Collect unique protein IDs
    protein_ids = set()
    for r in rows:
        pid = r['protein_id'].split('.')[0]
        if pid:
            protein_ids.add(pid)
    print(f'Unique protein IDs to annotate: {len(protein_ids)}')

    cache = load_cache()
    print(f'Domain cache: {len(cache)} entries')

    todo = [p for p in protein_ids if p not in cache]
    print(f'To fetch: {len(todo)}')

    for i, pid in enumerate(todo):
        if (i + 1) % 100 == 1:
            print(f'  Fetching {i + 1}/{len(todo)}...')
        domains = fetch_interpro_protein(pid)
        cache[pid] = domains
        if (i + 1) % 500 == 0:
            save_cache(cache)
        time.sleep(0.2)

    save_cache(cache)
    print(f'Domain cache now: {len(cache)} entries')

    with open(OUTPUT, 'w', newline='') as f:
        fieldnames = list(rows[0].keys()) + ['pfam_ids', 'pfam_names']
        writer = csv.DictWriter(f, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            pid = r['protein_id'].split('.')[0]
            domains = cache.get(pid, [])
            r['pfam_ids'] = ';'.join(d['pfam_id'] for d in domains)
            r['pfam_names'] = ';'.join(d['pfam_name'] for d in domains)
            writer.writerow(r)

    print(f'Output: {OUTPUT.name}')


if __name__ == '__main__':
    main()
