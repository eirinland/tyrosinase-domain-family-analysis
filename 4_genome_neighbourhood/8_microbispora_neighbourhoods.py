"""
Step 8: Fetch and compare neighbourhoods for all Microbispora-group proteins.

Fetches genome neighbourhoods for the 6 Microbispora proteins not in the main
pipeline, combines with the 1 already fetched, and compares their flanking genes.
"""

import csv
import re
import time
from collections import Counter
from pathlib import Path
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode

WORK = Path(__file__).parent
NEIGHBOURHOODS = WORK / 'neighbourhoods.tsv'
NCBI_EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
EMAIL = 'eirinlandsem1@gmail.com'
FLANK = 10

MICROBISPORA = {
    'A0A8H9LF69': ('BMMN01000009', 'MBD3146478.1', 'MET', 'M. bryophytorum'),
    'A0ABR8L1I5': ('JACXRZ010000011', 'MBD3144832.1', 'MET', 'M. bryophytorum'),
    'A0A1N6ZK30': ('FTNI01000007', 'SIR27185.1', 'ILE', 'M. rosea'),
    'A0A544XUX3': ('VIRM01000096', 'TQS08294.1', 'ILE', 'M. hainanensis'),
    'A0A5N6BWN4': ('VDMA02000006', 'KAB8184936.1', 'ILE', 'M. catharanthi'),
    'A0ABV0AFJ4': ('JBDJAW010000002', 'MEN3534112.1', 'ILE', 'M. maris'),
    'A0ABZ1SYY9': ('CP108085', 'WUP78280.1', 'ILE', 'M. hainanensis'),
}


def fetch_genbank(nuc_acc: str) -> str | None:
    params = urlencode({
        'db': 'nuccore', 'id': nuc_acc,
        'rettype': 'gb', 'retmode': 'text',
        'email': EMAIL, 'tool': 'ppo_gna',
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
            'start': start, 'end': end,
            'strand': '-' if complement else '+',
            'protein_id': q.get('protein_id', ''),
            'locus_tag': q.get('locus_tag', q.get('gene', '')),
            'product': q.get('product', 'hypothetical protein'),
        })
    result.sort(key=lambda x: x['start'])
    return result


def find_neighbours(features, protein_id, flank=FLANK):
    target_idx = None
    for i, feat in enumerate(features):
        if feat['protein_id'].split('.')[0] == protein_id.split('.')[0]:
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
    all_neighbourhoods = {}

    # Load existing A0A8H9LF69 from pipeline
    with open(NEIGHBOURHOODS) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['query_accession'] == 'A0A8H9LF69':
                acc = r['query_accession']
                all_neighbourhoods.setdefault(acc, []).append(r)

    print(f'Loaded A0A8H9LF69 from pipeline: {len(all_neighbourhoods.get("A0A8H9LF69", []))} entries')

    # Fetch the other 6
    to_fetch = {acc: info for acc, info in MICROBISPORA.items() if acc != 'A0A8H9LF69'}
    for acc, (nuc_acc, prot_id, variant, species) in to_fetch.items():
        print(f'Fetching {acc} ({species}) from {nuc_acc}...')
        time.sleep(0.4)
        gb_text = fetch_genbank(nuc_acc)
        if gb_text is None:
            print(f'  FAILED to fetch {nuc_acc}')
            continue

        features = parse_cds_features(gb_text)
        if not features:
            print(f'  No CDS features parsed from {nuc_acc}')
            continue

        target, neighbours = find_neighbours(features, prot_id)
        if target is None:
            print(f'  Could not find {prot_id} in {nuc_acc} ({len(features)} CDS)')
            continue

        print(f'  Found: {len(neighbours)} genes (±{FLANK})')
        for nb in neighbours:
            all_neighbourhoods.setdefault(acc, []).append({
                'query_accession': acc,
                'nucleotide_acc': nuc_acc,
                'offset': str(nb['offset']),
                'is_target': '1' if nb['is_target'] else '0',
                'product': nb['product'],
                'protein_id': nb['protein_id'],
                'locus_tag': nb['locus_tag'],
            })

    # Compare neighbourhoods
    print('\n' + '=' * 75)
    print('MICROBISPORA GROUP: FLANKING GENE COMPARISON')
    print('=' * 75)

    for acc in sorted(all_neighbourhoods.keys()):
        nuc, prot, variant, species = MICROBISPORA[acc]
        entries = all_neighbourhoods[acc]
        flanking = [e for e in entries if e.get('is_target') != '1']
        products = [e['product'].lower().strip() for e in flanking
                    if e['product'].lower().strip() not in
                    ('hypothetical protein', 'uncharacterized protein', '')]

        print(f'\n  {acc} ({species}, CuB_His2={variant})')
        print(f'    Contig: {nuc}, {len(flanking)} flanking genes')
        print(f'    Products:')
        for e in sorted(flanking, key=lambda x: int(x['offset'])):
            marker = '  ' if e.get('is_target') != '1' else '>>'
            print(f'      [{int(e["offset"]):+3d}] {e["product"]}')

    # Cross-compare: shared products
    print('\n' + '=' * 75)
    print('SHARED PRODUCTS BETWEEN MICROBISPORA NEIGHBOURS')
    print('=' * 75)

    acc_products = {}
    for acc in sorted(all_neighbourhoods.keys()):
        flanking = [e for e in all_neighbourhoods[acc] if e.get('is_target') != '1']
        products = set()
        for e in flanking:
            p = e['product'].lower().strip()
            if p not in ('hypothetical protein', 'uncharacterized protein', ''):
                products.add(p)
        acc_products[acc] = products

    all_accs = sorted(acc_products.keys())
    for i, a1 in enumerate(all_accs):
        for a2 in all_accs[i+1:]:
            shared = acc_products[a1] & acc_products[a2]
            union = acc_products[a1] | acc_products[a2]
            jaccard = len(shared) / len(union) if union else 0
            s1 = MICROBISPORA[a1][3]
            s2 = MICROBISPORA[a2][3]
            print(f'  {a1} ({s1}) vs {a2} ({s2}):')
            print(f'    Shared: {len(shared)}/{len(union)} (Jaccard={jaccard:.3f})')
            if shared:
                print(f'    Products: {"; ".join(sorted(shared)[:8])}')

    # Aggregate: most common flanking products across all 7
    print('\n' + '=' * 75)
    print('MOST COMMON FLANKING PRODUCTS ACROSS ALL MICROBISPORA')
    print('=' * 75)

    all_products = Counter()
    presence = Counter()
    for acc, entries in all_neighbourhoods.items():
        local = set()
        for e in entries:
            if e.get('is_target') == '1':
                continue
            p = e['product'].lower().strip()
            if p not in ('hypothetical protein', 'uncharacterized protein', ''):
                all_products[p] += 1
                local.add(p)
        for p in local:
            presence[p] += 1

    print(f'\n  Products found near ≥2 of 7 Microbispora proteins:')
    for prod, n_prot in presence.most_common():
        if n_prot >= 2:
            print(f'    {n_prot}/7 proteins: {prod} ({all_products[prod]} total)')


if __name__ == '__main__':
    main()
