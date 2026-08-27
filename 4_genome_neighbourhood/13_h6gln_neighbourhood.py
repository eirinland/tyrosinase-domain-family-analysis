"""
Genome neighbourhood analysis for the fungal H6Gln mononuclear group (46 structures).

Self-contained pipeline:
  1. Fetch genome cross-references from UniProt REST API
  2. Fetch gene neighbourhoods from NCBI Entrez (±10 flanking genes)
  3. Fetch Pfam annotations for neighbour proteins via UniProt
  4. Compute co-occurrence statistics (product, keyword, Pfam)
  5. Score BGC context (UstYa cluster / BGC accessory / none)

All outputs go to 4_genome_neighbourhood/h6gln_*.tsv

Usage: python 13_h6gln_neighbourhood.py
"""

import os
import csv
import json
import re
import sys
import time
import threading
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, quote

WORK = Path(__file__).parent
STAGE3 = WORK.parent / 'Super_reference_pipeline' / '3_noncanonical_analysis' / 'helix_and_gap_filtered_structures.tsv'

OUT_XREFS = WORK / 'h6gln_genome_crossrefs.tsv'
OUT_NB = WORK / 'h6gln_neighbourhoods.tsv'
OUT_PFAM = WORK / 'h6gln_pfam_map.tsv'
OUT_COOC_PROD = WORK / 'h6gln_cooccurrence_products.tsv'
OUT_COOC_KW = WORK / 'h6gln_cooccurrence_keywords.tsv'
OUT_COOC_PFAM = WORK / 'h6gln_cooccurrence_pfam.tsv'
OUT_CLUSTER = WORK / 'h6gln_cluster_membership.tsv'

NCBI_EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
UNIPROT_API = 'https://rest.uniprot.org'
# NCBI asks that Entrez requests identify the caller. Set NCBI_EMAIL to your
# own address before running; it is deliberately not hardcoded so this script
# can be shared without carrying a personal address.
EMAIL = os.environ.get('NCBI_EMAIL', '')
if not EMAIL:
    raise SystemExit('Set NCBI_EMAIL to the address NCBI Entrez should see, e.g. export NCBI_EMAIL=you@example.org')
FLANK = 10

POSS = ['CuA_His1','CuA_His2','CuA_His3','CuB_His1','CuB_His2','CuB_His3']

MARKERS = {
    'TF':          {'PF00172', 'PF04082', 'PF11951'},
    'transporter': {'PF07690', 'PF00083', 'PF00005', 'PF00664'},
    'P450':        {'PF00067'},
    'OMT':         {'PF00891'},
}
USTYA = 'PF11807'

STOP_WORDS = {
    'protein', 'hypothetical', 'putative', 'probable', 'uncharacterized',
    'predicted', 'domain-containing', 'like', 'family', 'containing',
    'related', 'with', 'type', 'that', 'this', 'from', 'have', 'been',
    'subunit', 'partial', 'unnamed', 'product', 'and', 'the', 'gene',
}


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


def get_h6gln_fungal_accessions():
    with open(STAGE3) as f:
        raw = f.read().replace('\r', '')
    rows = list(csv.DictReader(raw.splitlines(), delimiter='\t'))
    return [r for r in rows
            if r['CuB_His3'] == 'GLN'
            and r.get('classification', '') == 'mononuclear'
            and r.get('phylum', '') in ('Ascomycota', 'Basidiobolomycota')]


# --- Step 1: UniProt cross-refs ---

def fetch_uniprot_batch(accessions):
    query = ' OR '.join(f'accession:{acc}' for acc in accessions)
    url = (f'{UNIPROT_API}/uniprotkb/search?'
           f'query={quote(query)}&fields=accession,xref_embl,organism_id,lineage'
           f'&format=json&size={len(accessions)}')
    req = Request(url, headers={'Accept': 'application/json'})
    try:
        with urlopen(req, timeout=60) as resp:
            data = json.loads(resp.read())
        return {r['primaryAccession']: r for r in data.get('results', [])}
    except (HTTPError, URLError, json.JSONDecodeError) as e:
        print(f'  UniProt batch error: {e}', file=sys.stderr)
        return {}


def extract_embl_refs(entry):
    refs = []
    for xref in entry.get('uniProtKBCrossReferences', []):
        if xref.get('database') != 'EMBL':
            continue
        props = {p['key']: p['value'] for p in xref.get('properties', [])}
        protein_id = props.get('ProteinId', '')
        if protein_id in ('-', ''):
            continue
        refs.append({
            'nucleotide_acc': xref.get('id', ''),
            'protein_id': protein_id,
            'molecule_type': props.get('MoleculeType', ''),
        })
    return refs


def step1_crossrefs(accessions):
    print(f'\n=== Step 1: Fetching UniProt cross-refs for {len(accessions)} accessions ===')
    results = []
    for i in range(0, len(accessions), 100):
        batch = accessions[i:i+100]
        data = fetch_uniprot_batch(batch)
        for acc in batch:
            entry = data.get(acc)
            if entry is None:
                continue
            org_id = str(entry.get('organism', {}).get('taxonId', ''))
            lineage = ';'.join(entry.get('organism', {}).get('lineage', [])[:5])
            embl = extract_embl_refs(entry)
            if not embl:
                results.append({'accession': acc, 'organism_id': org_id, 'lineage': lineage,
                                'nucleotide_acc': '', 'protein_id': '', 'molecule_type': '', 'n_embl_refs': '0'})
            else:
                genomic = [r for r in embl if r['molecule_type'] == 'Genomic_DNA']
                best = genomic[0] if genomic else embl[0]
                results.append({'accession': acc, 'organism_id': org_id, 'lineage': lineage,
                                **best, 'n_embl_refs': str(len(embl))})
        time.sleep(0.5)

    fields = ['accession', 'organism_id', 'lineage', 'nucleotide_acc', 'protein_id', 'molecule_type', 'n_embl_refs']
    with open(OUT_XREFS, 'w', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    n_with = sum(1 for r in results if r['nucleotide_acc'])
    print(f'  {n_with}/{len(results)} have nucleotide cross-refs -> {OUT_XREFS.name}')
    return [r for r in results if r['nucleotide_acc']]


# --- Step 2: NCBI neighbourhoods ---

def fetch_genbank(nuc_acc):
    params = urlencode({'db': 'nuccore', 'id': nuc_acc, 'rettype': 'gb',
                        'retmode': 'text', 'email': EMAIL, 'tool': 'ppo_gna'})
    for attempt in range(3):
        try:
            with urlopen(f'{NCBI_EFETCH}?{params}', timeout=120) as resp:
                return resp.read().decode('utf-8', errors='replace')
        except (HTTPError, URLError):
            if attempt < 2:
                time.sleep(2 ** attempt)
    return None


def parse_cds_features(gb_text):
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
            'start': start, 'end': end, 'strand': '-' if complement else '+',
            'protein_id': q.get('protein_id', ''), 'locus_tag': q.get('locus_tag', q.get('gene', '')),
            'product': q.get('product', 'hypothetical protein'),
        })
    result.sort(key=lambda x: x['start'])
    return result


def step2_neighbourhoods(xrefs):
    print(f'\n=== Step 2: Fetching NCBI neighbourhoods for {len(xrefs)} entries ===')
    limiter = RateLimiter(3)
    by_nuc = {}
    for r in xrefs:
        by_nuc.setdefault(r['nucleotide_acc'], []).append(r)
    nuc_list = list(by_nuc.keys())
    print(f'  {len(nuc_list)} unique nucleotide accessions')

    fields = ['query_accession', 'nucleotide_acc', 'offset', 'is_target',
              'gene_start', 'gene_end', 'strand', 'protein_id', 'locus_tag', 'product']
    outf = open(OUT_NB, 'w', newline='')
    writer = csv.DictWriter(outf, delimiter='\t', fieldnames=fields)
    writer.writeheader()

    n_ok = n_fail = 0

    def rate_fetch(nuc):
        limiter.acquire()
        return nuc, fetch_genbank(nuc)

    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = {pool.submit(rate_fetch, nuc): nuc for nuc in nuc_list}
        for i, future in enumerate(as_completed(futures), 1):
            nuc_acc, gb_text = future.result()
            entries = by_nuc[nuc_acc]
            if i % 50 == 0 or i == len(nuc_list):
                print(f'  {i}/{len(nuc_list)} ({n_ok} ok, {n_fail} fail)')
            if gb_text is None:
                n_fail += len(entries)
                continue
            features = parse_cds_features(gb_text)
            if not features:
                n_fail += len(entries)
                continue
            for entry in entries:
                pid = entry['protein_id'].split('.')[0]
                target_idx = None
                for j, feat in enumerate(features):
                    if feat['protein_id'].split('.')[0] == pid:
                        target_idx = j
                        break
                if target_idx is None:
                    n_fail += 1
                    continue
                n_ok += 1
                start = max(0, target_idx - FLANK)
                end = min(len(features), target_idx + FLANK + 1)
                for j in range(start, end):
                    offset = j - target_idx
                    writer.writerow({
                        'query_accession': entry['accession'], 'nucleotide_acc': nuc_acc,
                        'offset': offset, 'is_target': '1' if offset == 0 else '0',
                        'gene_start': features[j]['start'], 'gene_end': features[j]['end'],
                        'strand': features[j]['strand'], 'protein_id': features[j]['protein_id'],
                        'locus_tag': features[j]['locus_tag'], 'product': features[j]['product'],
                    })

    outf.close()
    print(f'  {n_ok} neighbourhoods extracted, {n_fail} failed -> {OUT_NB.name}')
    return n_ok


# --- Step 3: Pfam via UniProt ---

def step3_pfam():
    print(f'\n=== Step 3: Fetching Pfam annotations for neighbour proteins ===')
    pids = set()
    with open(OUT_NB) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['is_target'] != '1' and r['protein_id']:
                pids.add(r['protein_id'].split('.')[0])
    print(f'  {len(pids)} unique neighbour protein IDs')

    pid2pfam = {}
    pids_list = sorted(pids)
    for i in range(0, len(pids_list), 100):
        batch = pids_list[i:i+100]
        query = ' OR '.join(f'xref:embl-{pid}' for pid in batch)
        url = (f'{UNIPROT_API}/uniprotkb/search?query={quote(query)}'
               f'&fields=accession,xref_pfam&format=json&size=500')
        req = Request(url, headers={'Accept': 'application/json'})
        try:
            with urlopen(req, timeout=60) as resp:
                data = json.loads(resp.read())
            for entry in data.get('results', []):
                acc = entry['primaryAccession']
                pfams = []
                for xref in entry.get('uniProtKBCrossReferences', []):
                    if xref.get('database') == 'Pfam':
                        pfam_id = xref.get('id', '')
                        props = {p['key']: p['value'] for p in xref.get('properties', [])}
                        name = props.get('EntryName', '')
                        pfams.append(f'{pfam_id}|{name}')
                if pfams:
                    for xref2 in entry.get('uniProtKBCrossReferences', []):
                        if xref2.get('database') == 'EMBL':
                            embl_props = {p['key']: p['value'] for p in xref2.get('properties', [])}
                            epid = embl_props.get('ProteinId', '').split('.')[0]
                            if epid in pids:
                                pid2pfam[epid] = ';'.join(pfams)
        except (HTTPError, URLError, json.JSONDecodeError) as e:
            print(f'  Pfam batch error: {e}', file=sys.stderr)
        if (i // 100) % 5 == 0:
            print(f'  batch {i//100 + 1}/{(len(pids_list)+99)//100}')
        time.sleep(0.5)

    with open(OUT_PFAM, 'w', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=['protein_id', 'pfam'])
        w.writeheader()
        for pid, pfam in sorted(pid2pfam.items()):
            w.writerow({'protein_id': pid, 'pfam': pfam})

    print(f'  {len(pid2pfam)} proteins with Pfam annotations -> {OUT_PFAM.name}')


# --- Step 4: Co-occurrence + cluster membership ---

def norm_product(p):
    p = re.sub(r'\s+', ' ', p.strip().lower())
    if p in ('hypothetical protein', 'uncharacterized protein',
             'unnamed protein product', 'putative protein', ''):
        return ''
    return p


def step4_analysis(target_accs):
    print(f'\n=== Step 4: Co-occurrence analysis ===')

    targets = set(target_accs)
    nb = [r for r in csv.DictReader(open(OUT_NB), delimiter='\t')
          if r['query_accession'] in targets]
    neigh = [r for r in nb if r['is_target'] != '1']
    queries_with_data = {r['query_accession'] for r in nb}
    n_total = len(targets)

    # Product co-occurrence
    per_query_prod = defaultdict(set)
    prod_total = Counter()
    for r in neigh:
        np_ = norm_product(r.get('product', ''))
        if not np_:
            continue
        per_query_prod[r['query_accession']].add(np_)
        prod_total[np_] += 1

    prod_qcount = Counter()
    for q, prods in per_query_prod.items():
        for p in prods:
            prod_qcount[p] += 1

    n_with_info = len(per_query_prod)
    prod_rows = []
    for p, qc in prod_qcount.most_common():
        prod_rows.append([p, qc, f'{qc/n_total:.3f}', f'{qc/n_with_info:.3f}' if n_with_info else '0', prod_total[p]])
    with open(OUT_COOC_PROD, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['product', 'n_queries', 'freq_of_all_queries', 'freq_of_queries_with_neighbours', 'total_occurrences'])
        w.writerows(prod_rows)

    # Keyword co-occurrence
    per_query_kw = defaultdict(set)
    kw_total = Counter()
    for r in neigh:
        prod = r.get('product', '').lower()
        if norm_product(prod) == '':
            continue
        toks = {w for w in re.split(r'[^a-z0-9]+', prod) if len(w) > 3 and w not in STOP_WORDS}
        for w in toks:
            per_query_kw[r['query_accession']].add(w)
            kw_total[w] += 1
    kw_qcount = Counter()
    for q, kws in per_query_kw.items():
        for w in kws:
            kw_qcount[w] += 1
    with open(OUT_COOC_KW, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['keyword', 'n_queries', 'freq_of_all_queries', 'total_occurrences'])
        for kw, qc in kw_qcount.most_common():
            if qc >= 2:
                w.writerow([kw, qc, f'{qc/n_total:.3f}', kw_total[kw]])

    # Pfam co-occurrence
    pid2pfam = {}
    if OUT_PFAM.exists():
        for r in csv.DictReader(open(OUT_PFAM), delimiter='\t'):
            fams = [f for f in r.get('pfam', '').split(';') if f]
            pid2pfam[r['protein_id'].split('.')[0]] = fams

    per_query_pfam = defaultdict(set)
    pfam_total = Counter()
    pfam_name = {}
    for r in neigh:
        pid = r.get('protein_id', '').split('.')[0]
        for fam in pid2pfam.get(pid, []):
            fid, _, nm = fam.partition('|')
            per_query_pfam[r['query_accession']].add(fid)
            pfam_total[fid] += 1
            if nm:
                pfam_name[fid] = nm

    pfam_qcount = Counter()
    for q, fams in per_query_pfam.items():
        for f in fams:
            pfam_qcount[f] += 1

    n_with_pfam = len(per_query_pfam)
    with open(OUT_COOC_PFAM, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['pfam_id', 'pfam_name', 'n_queries', 'freq_of_all_queries',
                     'freq_of_queries_with_pfam', 'total_occurrences'])
        for fid, qc in pfam_qcount.most_common():
            w.writerow([fid, pfam_name.get(fid, ''), qc, f'{qc/n_total:.3f}',
                        f'{qc/n_with_pfam:.3f}' if n_with_pfam else '0', pfam_total[fid]])

    # Cluster membership (BGC scoring)
    sp_map, gen_map = {}, {}
    with open(STAGE3) as f:
        raw = f.read().replace('\r', '')
    for r in csv.DictReader(raw.splitlines(), delimiter='\t'):
        sp_map[r['accession']] = r.get('species', '')
        gen_map[r['accession']] = r.get('genus', '')

    nb_by_q = defaultdict(list)
    for r in neigh:
        if r['protein_id']:
            nb_by_q[r['query_accession']].append((int(r['offset']), r['protein_id'].split('.')[0]))

    clust_rows = []
    for q in sorted(queries_with_data):
        fams_all = set()
        ust_offsets = []
        for off, pid in nb_by_q.get(q, []):
            fams = set()
            for fam_str in pid2pfam.get(pid, []):
                fams.add(fam_str.split('|')[0])
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
    n_data = len(queries_with_data)
    n_ust = sum(1 for r in clust_rows if r['bgc_context'] == 'UstYa_cluster')
    n_bgc = sum(1 for r in clust_rows if r['bgc_context'] != 'none')
    n_none = sum(1 for r in clust_rows if r['bgc_context'] == 'none')
    print(f'\n=== Results ===')
    print(f'H6Gln fungal targets: {n_total}')
    print(f'  with neighbourhood data: {n_data}')
    print(f'  UstYa cluster: {n_ust}')
    print(f'  BGC accessory: {n_bgc - n_ust}')
    print(f'  no BGC context: {n_none}')
    print()
    print(f'Top co-occurring products:')
    for row in prod_rows[:15]:
        print(f'  {row[1]:3d} ({float(row[2])*100:4.1f}%)  {row[0][:70]}')
    print()
    print(f'Top co-occurring Pfams:')
    for fid, qc in pfam_qcount.most_common(15):
        print(f'  {qc:3d} ({qc/n_total*100:4.1f}%)  {fid} {pfam_name.get(fid,"")}')
    print()
    print(f'Outputs: {OUT_COOC_PROD.name}, {OUT_COOC_KW.name}, {OUT_COOC_PFAM.name}, {OUT_CLUSTER.name}')


def main():
    targets = get_h6gln_fungal_accessions()
    accs = [r['accession'] for r in targets]
    print(f'Fungal H6Gln mononuclear: {len(accs)} accessions')

    xrefs = step1_crossrefs(accs)
    if not xrefs:
        print('No cross-refs found, stopping.')
        return

    n_ok = step2_neighbourhoods(xrefs)
    if n_ok == 0:
        print('No neighbourhoods extracted, stopping.')
        return

    step3_pfam()
    step4_analysis(accs)


if __name__ == '__main__':
    main()
