"""
Build accession -> taxonomy mapping from FASTA headers + UniProt taxonomy API.
Outputs taxonomy_lookup.csv with columns: accession, taxid, kingdom, phylum, genus
"""

import csv
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

FASTA = Path('/cluster/work/projects/nn1003k/eirin/bioinf/HMMsearch/query.fasta')
VECTORS = Path(__file__).parent.parent / 'position_vectors.csv'
OUT = Path(__file__).parent / 'taxonomy_lookup.csv'

PROXY = 'http://10.63.2.48:3128/'
BATCH = 80


def parse_fasta_taxids(fasta_path):
    acc_tax = {}
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                acc = line[1:].split('|')[0]
                m = re.search(r'taxID:(\d+)', line)
                if m:
                    acc_tax[acc] = int(m.group(1))
    return acc_tax


def classify(entry):
    lineage_names = [l.get('scientificName', '') for l in entry.get('lineage', [])]
    sci = entry.get('scientificName', '')
    genus = sci.split()[0] if sci else '?'

    kingdom = '?'
    phylum = '?'
    for l in entry.get('lineage', []):
        rank = l.get('rank', '')
        name = l.get('scientificName', '')
        if rank == 'phylum' and phylum == '?':
            phylum = name

    if 'Fungi' in lineage_names:
        kingdom = 'Fungi'
    elif any(n in lineage_names for n in ('Viridiplantae', 'Streptophyta', 'Embryophyta')):
        kingdom = 'Plants'
    elif any(n in lineage_names for n in ('Oomycota', 'Stramenopiles')):
        kingdom = 'Oomycota'
    elif 'Metazoa' in lineage_names:
        kingdom = 'Animals'
    elif 'Bacteria' in lineage_names:
        kingdom = 'Bacteria'
    elif 'Archaea' in lineage_names:
        kingdom = 'Archaea'

    return kingdom, phylum, genus


def fetch_batch(taxids):
    handler = urllib.request.ProxyHandler({'http': PROXY, 'https': PROXY})
    opener = urllib.request.build_opener(handler)

    query = '+OR+'.join(f'tax_id:{t}' for t in taxids)
    url = f'https://rest.uniprot.org/taxonomy/search?query={query}&fields=id,scientific_name,lineage&format=json&size=500'

    results = {}
    try:
        req = urllib.request.Request(url)
        req.add_header('User-Agent', 'Python taxonomy-fetch/1.0')
        resp = opener.open(req, timeout=60)
        data = json.loads(resp.read().decode())
        for entry in data.get('results', []):
            tid = entry.get('taxonId')
            results[tid] = classify(entry)
    except Exception as e:
        print(f"  error: {e}", file=sys.stderr)
    return results


def main():
    vec_accs = set()
    with open(VECTORS) as f:
        for row in csv.DictReader(f):
            vec_accs.add(row['accession'])
    print(f"Accessions in vectors: {len(vec_accs)}")

    acc_tax = parse_fasta_taxids(FASTA)
    acc_tax = {a: t for a, t in acc_tax.items() if a in vec_accs}
    print(f"Accessions with taxID: {len(acc_tax)}")

    unique_taxids = sorted(set(acc_tax.values()))
    print(f"Unique taxIDs to look up: {len(unique_taxids)}")

    tax_info = {}
    for i in range(0, len(unique_taxids), BATCH):
        batch = unique_taxids[i:i+BATCH]
        print(f"  {i+1}-{i+len(batch)} / {len(unique_taxids)}...", end='', flush=True)
        info = fetch_batch(batch)
        tax_info.update(info)
        print(f" {len(info)} ok")
        if i + BATCH < len(unique_taxids):
            time.sleep(0.3)

    print(f"Resolved {len(tax_info)} / {len(unique_taxids)} taxIDs")

    with open(OUT, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['accession', 'taxid', 'kingdom', 'phylum', 'genus'])
        for acc in sorted(vec_accs):
            tid = acc_tax.get(acc, 0)
            if tid in tax_info:
                kingdom, phylum, genus = tax_info[tid]
            else:
                kingdom, phylum, genus = '?', '?', '?'
            w.writerow([acc, tid, kingdom, phylum, genus])

    # Summary
    kingdoms = {}
    with open(OUT) as f:
        for row in csv.DictReader(f):
            k = row['kingdom']
            kingdoms[k] = kingdoms.get(k, 0) + 1
    for k, c in sorted(kingdoms.items(), key=lambda x: -x[1]):
        print(f"  {k}: {c}")

    print(f"Written: {OUT}")


if __name__ == '__main__':
    main()
