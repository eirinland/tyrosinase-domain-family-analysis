"""
Step 3: Fetch gene neighbourhoods from NCBI/ENA.

For each protein with a nucleotide cross-reference:
  1. Fetch the nucleotide record (GenBank flat file) from NCBI
  2. Find the CDS feature matching our protein
  3. Extract ±N flanking genes
  4. Record gene product annotations

Uses NCBI Entrez E-utilities with concurrent fetching (3 req/sec).
Output: neighbourhoods.tsv (one row per flanking gene)

Usage: python 3_fetch_neighbourhoods.py [--flank 10] [--workers 8]
"""

import argparse
import os
import csv
import re
import sys
import time
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode

WORK = Path(__file__).parent
CROSSREFS = WORK / 'genome_crossrefs.tsv'
OUTPUT = WORK / 'neighbourhoods.tsv'
PROGRESS = WORK / 'neighbourhood_progress.txt'

NCBI_EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'


class RateLimiter:
    def __init__(self, max_per_second=3):
        self.min_interval = 1.0 / max_per_second
        self.lock = threading.Lock()
        self.last = 0.0

    def acquire(self):
        with self.lock:
            now = time.time()
            wait = self.min_interval - (now - self.last)
            if wait > 0:
                time.sleep(wait)
            self.last = time.time()


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--flank', type=int, default=10,
                   help='Number of flanking genes on each side (default: 10)')
    p.add_argument('--email', default=os.environ.get('NCBI_EMAIL', ''),
                   help='Email for NCBI Entrez. Defaults to $NCBI_EMAIL; '
                        'NCBI asks that requests identify the caller.')
    p.add_argument('--workers', type=int, default=8,
                   help='Number of concurrent download threads (default: 8)')
    return p.parse_args()


def fetch_genbank(nuc_acc: str, email: str) -> str | None:
    params = urlencode({
        'db': 'nuccore',
        'id': nuc_acc,
        'rettype': 'gb',
        'retmode': 'text',
        'email': email,
        'tool': 'ppo_gna',
    })
    url = f'{NCBI_EFETCH}?{params}'
    for attempt in range(3):
        try:
            with urlopen(url, timeout=120) as resp:
                return resp.read().decode('utf-8', errors='replace')
        except (HTTPError, URLError):
            if attempt < 2:
                time.sleep(2 ** attempt)
    return None


def parse_cds_features(gb_text: str) -> list[dict]:
    """Minimal GenBank CDS parser. Extracts location, product, protein_id, locus_tag."""
    features = []
    in_cds = False
    current = {}
    current_qualifier = None
    current_value = ''

    for line in gb_text.split('\n'):
        if line.startswith('     CDS '):
            if in_cds and current:
                features.append(current)
            in_cds = True
            current = {'location_raw': line[21:].strip(), 'qualifiers': {}}
            current_qualifier = None
            current_value = ''
        elif in_cds and line.startswith('                     /'):
            if current_qualifier:
                current['qualifiers'][current_qualifier] = current_value.strip('"')
            match = re.match(r'\s+/(\w+)(?:="?(.*))?', line)
            if match:
                current_qualifier = match.group(1)
                current_value = match.group(2) or ''
                if current_value.endswith('"'):
                    current_value = current_value[:-1]
            else:
                current_qualifier = None
                current_value = ''
        elif in_cds and line.startswith('                     ') and current_qualifier:
            current_value += ' ' + line.strip().strip('"')
        elif in_cds and not line.startswith('                     '):
            if current_qualifier:
                current['qualifiers'][current_qualifier] = current_value.strip('"')
            features.append(current)
            in_cds = False
            current = {}
            current_qualifier = None

    if in_cds and current:
        if current_qualifier:
            current['qualifiers'][current_qualifier] = current_value.strip('"')
        features.append(current)

    result = []
    for feat in features:
        q = feat['qualifiers']
        loc = feat['location_raw']
        complement = 'complement' in loc
        numbers = re.findall(r'(\d+)', loc)
        if len(numbers) >= 2:
            start, end = int(numbers[0]), int(numbers[-1])
        else:
            continue

        result.append({
            'start': start,
            'end': end,
            'strand': '-' if complement else '+',
            'protein_id': q.get('protein_id', ''),
            'locus_tag': q.get('locus_tag', q.get('gene', '')),
            'product': q.get('product', 'hypothetical protein'),
            'pseudo': 'pseudo' in q or 'pseudogene' in q,
        })

    result.sort(key=lambda x: x['start'])
    return result


def find_target_and_neighbours(features: list[dict], protein_id: str, flank: int):
    """Find target CDS by protein_id and return it + flanking genes."""
    target_idx = None
    for i, feat in enumerate(features):
        pid = feat['protein_id'].split('.')[0]
        if pid == protein_id.split('.')[0]:
            target_idx = i
            break

    if target_idx is None:
        return None, []

    start = max(0, target_idx - flank)
    end = min(len(features), target_idx + flank + 1)
    neighbours = []
    for i in range(start, end):
        offset = i - target_idx
        neighbours.append({**features[i], 'offset': offset, 'is_target': offset == 0})

    return features[target_idx], neighbours


def main():
    args = parse_args()
    limiter = RateLimiter(3)

    with open(CROSSREFS) as f:
        crossrefs = list(csv.DictReader(f, delimiter='\t'))

    targets = [r for r in crossrefs if r['nucleotide_acc']]
    print(f'{len(targets)} accessions with nucleotide cross-references')

    done = set()
    if PROGRESS.exists():
        done = set(PROGRESS.read_text().strip().split('\n'))
        print(f'  Resuming: {len(done)} already processed')

    todo = [r for r in targets if r['accession'] not in done]
    print(f'  Remaining: {len(todo)}')

    by_nuc = {}
    for r in todo:
        by_nuc.setdefault(r['nucleotide_acc'], []).append(r)
    nuc_list = list(by_nuc.keys())
    print(f'  Unique nucleotide accessions: {len(nuc_list)}')

    write_header = not OUTPUT.exists() or len(done) == 0
    outf = open(OUTPUT, 'a', newline='')
    writer = csv.DictWriter(outf, delimiter='\t', fieldnames=[
        'query_accession', 'nucleotide_acc', 'offset', 'is_target',
        'gene_start', 'gene_end', 'strand', 'protein_id',
        'locus_tag', 'product',
    ])
    if write_header:
        writer.writeheader()

    progf = open(PROGRESS, 'a')
    n_success = 0
    n_fail = 0
    n_fetched = 0
    t_start = time.time()

    def rate_limited_fetch(nuc_acc):
        limiter.acquire()
        return nuc_acc, fetch_genbank(nuc_acc, args.email)

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(rate_limited_fetch, nuc): nuc for nuc in nuc_list}

        for future in as_completed(futures):
            n_fetched += 1
            try:
                nuc_acc, gb_text = future.result()
            except Exception:
                nuc_acc = futures[future]
                gb_text = None

            entries = by_nuc[nuc_acc]

            if n_fetched % 200 == 0 or n_fetched == len(nuc_list):
                elapsed = time.time() - t_start
                rate = n_fetched / elapsed if elapsed > 0 else 0
                eta = (len(nuc_list) - n_fetched) / rate if rate > 0 else 0
                print(f'  {n_fetched}/{len(nuc_list)} records '
                      f'({n_success} ok, {n_fail} fail, '
                      f'{rate:.1f}/s, ETA {eta / 60:.0f}m)')

            if gb_text is None:
                n_fail += len(entries)
                for e in entries:
                    progf.write(e['accession'] + '\n')
                continue

            features = parse_cds_features(gb_text)
            if not features:
                n_fail += len(entries)
                for e in entries:
                    progf.write(e['accession'] + '\n')
                continue

            for entry in entries:
                target, neighbours = find_target_and_neighbours(
                    features, entry['protein_id'], args.flank)

                if target is None:
                    n_fail += 1
                else:
                    n_success += 1
                    for nb in neighbours:
                        writer.writerow({
                            'query_accession': entry['accession'],
                            'nucleotide_acc': nuc_acc,
                            'offset': nb['offset'],
                            'is_target': '1' if nb['is_target'] else '0',
                            'gene_start': nb['start'],
                            'gene_end': nb['end'],
                            'strand': nb['strand'],
                            'protein_id': nb['protein_id'],
                            'locus_tag': nb['locus_tag'],
                            'product': nb['product'],
                        })

                progf.write(entry['accession'] + '\n')

            if n_fetched % 50 == 0:
                outf.flush()
                progf.flush()

    outf.close()
    progf.close()
    print(f'\nDone. {n_success} neighbourhoods extracted, {n_fail} failed.')
    print(f'Output: {OUTPUT.name}')


if __name__ == '__main__':
    main()
