"""
Step 2: Fetch genome cross-references from UniProt for all target accessions.

For each UniProt accession, retrieves:
  - Source nucleotide accession(s) (EMBL/GenBank)
  - Gene coordinates on the source sequence
  - Proteome / genome assembly ID (if available)

Uses UniProt REST API with batched requests.
Output: genome_crossrefs.tsv
"""

import csv
import json
import time
import sys
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError
from urllib.parse import quote

WORK = Path(__file__).parent
INPUT = WORK / 'target_accessions.tsv'
OUTPUT = WORK / 'genome_crossrefs.tsv'
FAILED = WORK / 'fetch_failed.txt'

UNIPROT_API = 'https://rest.uniprot.org'
BATCH_SIZE = 100
FIELDS = 'accession,xref_embl,xref_geneid,organism_id,lineage'


def fetch_uniprot_batch(accessions: list[str]) -> dict:
    """Query UniProt for a batch of accessions, return parsed JSON entries."""
    query = ' OR '.join(f'accession:{acc}' for acc in accessions)
    url = (
        f'{UNIPROT_API}/uniprotkb/search?'
        f'query={quote(query)}'
        f'&fields={FIELDS}'
        f'&format=json'
        f'&size={len(accessions)}'
    )
    req = Request(url, headers={'Accept': 'application/json'})
    try:
        with urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read())
        return {r['primaryAccession']: r for r in data.get('results', [])}
    except (HTTPError, URLError, json.JSONDecodeError) as e:
        print(f'  Error fetching batch: {e}', file=sys.stderr)
        return {}


def extract_embl_refs(entry: dict) -> list[dict]:
    """Extract EMBL nucleotide cross-references from a UniProt entry."""
    refs = []
    for xref in entry.get('uniProtKBCrossReferences', []):
        if xref.get('database') != 'EMBL':
            continue
        props = {p['key']: p['value'] for p in xref.get('properties', [])}
        nuc_id = xref.get('id', '')
        protein_id = props.get('ProteinId', '')
        mol_type = props.get('MoleculeType', '')
        # Skip entries without protein ID (non-coding, etc.)
        if protein_id in ('-', ''):
            continue
        refs.append({
            'nucleotide_acc': nuc_id,
            'protein_id': protein_id,
            'molecule_type': mol_type,
        })
    return refs


def main():
    with open(INPUT) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))

    unique_accs = sorted(set(r['accession'] for r in rows))
    print(f'Fetching genome cross-references for {len(unique_accs)} unique accessions')

    # Load any existing progress
    done = set()
    if OUTPUT.exists():
        with open(OUTPUT) as f:
            for r in csv.DictReader(f, delimiter='\t'):
                done.add(r['accession'])
        print(f'  Resuming: {len(done)} already fetched')

    todo = [a for a in unique_accs if a not in done]
    if not todo:
        print('All accessions already fetched.')
        return

    print(f'  Remaining: {len(todo)}')

    # Open output for appending
    write_header = not OUTPUT.exists() or len(done) == 0
    outf = open(OUTPUT, 'a', newline='')
    writer = csv.DictWriter(outf, delimiter='\t', fieldnames=[
        'accession', 'organism_id', 'lineage', 'nucleotide_acc',
        'protein_id', 'molecule_type', 'n_embl_refs',
    ])
    if write_header:
        writer.writeheader()

    failed = []
    n_with_refs = 0
    n_no_refs = 0

    for i in range(0, len(todo), BATCH_SIZE):
        batch = todo[i:i + BATCH_SIZE]
        batch_num = i // BATCH_SIZE + 1
        total_batches = (len(todo) + BATCH_SIZE - 1) // BATCH_SIZE

        if batch_num % 10 == 1 or batch_num == total_batches:
            print(f'  Batch {batch_num}/{total_batches} '
                  f'({n_with_refs} with refs, {n_no_refs} without, {len(failed)} failed)')

        results = fetch_uniprot_batch(batch)

        for acc in batch:
            entry = results.get(acc)
            if entry is None:
                failed.append(acc)
                continue

            org_id = str(entry.get('organism', {}).get('taxonId', ''))
            lineage_list = entry.get('organism', {}).get('lineage', [])
            lineage = ';'.join(lineage_list[:5]) if lineage_list else ''

            embl_refs = extract_embl_refs(entry)
            if not embl_refs:
                writer.writerow({
                    'accession': acc,
                    'organism_id': org_id,
                    'lineage': lineage,
                    'nucleotide_acc': '',
                    'protein_id': '',
                    'molecule_type': '',
                    'n_embl_refs': '0',
                })
                n_no_refs += 1
                continue

            n_with_refs += 1
            # Write the first genomic ref (prefer Genomic_DNA over mRNA)
            genomic = [r for r in embl_refs if r['molecule_type'] == 'Genomic_DNA']
            best = genomic[0] if genomic else embl_refs[0]
            writer.writerow({
                'accession': acc,
                'organism_id': org_id,
                'lineage': lineage,
                'nucleotide_acc': best['nucleotide_acc'],
                'protein_id': best['protein_id'],
                'molecule_type': best['molecule_type'],
                'n_embl_refs': str(len(embl_refs)),
            })

        outf.flush()
        time.sleep(0.5)  # respect rate limits

    outf.close()

    if failed:
        with open(FAILED, 'w') as f:
            f.write('\n'.join(failed) + '\n')

    print(f'\nDone. {n_with_refs} with genome refs, {n_no_refs} without, {len(failed)} failed.')
    print(f'Output: {OUTPUT.name}')


if __name__ == '__main__':
    main()
