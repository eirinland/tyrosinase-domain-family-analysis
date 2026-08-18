"""
Genome neighbourhood analysis for all manuscript groups not yet tested.

Canonical groups (from position_vectors.csv):
  E195R (136), H230Y (206), G46N (1047), G46I (164), N205K (111),
  N205R (143), F227W (84), V218R (104)

Non-canonical groups (from helix_and_gap_filtered_structures.tsv):
  Phytophthora ALN|HHH (20), Trichinella HHH|YTQ (14),
  Microbispora no_cu (8), Clonorchis/Opisthorchis (6), Oomycota no_cu (52)

Runs: UniProt cross-refs → NCBI neighbourhoods → Pfam (ID-mapping) → co-occurrence + BGC scoring.
All outputs go to 4_genome_neighbourhood/groups/<group_name>/

Usage: python 15_all_groups_neighbourhood.py [--group GROUP_NAME]
"""

import argparse
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
BASE = WORK.parent
PV_FILE = BASE / '2_canonical_analysis' / 'position_vectors.csv'
NC_FILE = BASE / '3_noncanonical_analysis' / 'helix_and_gap_filtered_structures.tsv'

NCBI_EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
UNIPROT_API = 'https://rest.uniprot.org'
EMAIL = 'eirinlandsem1@gmail.com'
FLANK = 10

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


def define_groups():
    """Return dict of group_name -> list of accessions."""
    groups = {}

    # Canonical groups from position vectors
    pv = {}
    with open(PV_FILE) as f:
        for r in csv.DictReader(f):
            pv[r['accession']] = r

    groups['E195R'] = [a for a, r in pv.items() if r['Glu195'] == 'R']
    groups['H230Y'] = [a for a, r in pv.items() if r['His230'] == 'Y']
    groups['G46N'] = [a for a, r in pv.items() if r['Gly46'] == 'N']
    groups['G46I'] = [a for a, r in pv.items() if r['Gly46'] == 'I']
    groups['N205K'] = [a for a, r in pv.items() if r['Asn205'] == 'K']
    groups['N205R'] = [a for a, r in pv.items() if r['Asn205'] == 'R']
    groups['F227W'] = [a for a, r in pv.items() if r['Phe227'] == 'W']
    groups['V218R'] = [a for a, r in pv.items() if r['Val218'] == 'R']

    # Non-canonical groups
    with open(NC_FILE) as f:
        raw = f.read().replace('\r', '')
    nc = list(csv.DictReader(raw.splitlines(), delimiter='\t'))

    groups['Phytophthora_ALN_HHH'] = [
        r['accession'] for r in nc
        if r['CuA_His1'] == 'ALA' and r['CuA_His2'] == 'LEU' and r['CuA_His3'] == 'ASN'
        and r['CuB_His1'] == 'HIS' and r['CuB_His2'] == 'HIS' and r['CuB_His3'] == 'HIS']

    groups['Trichinella_HHH_YTQ'] = [
        r['accession'] for r in nc
        if r['CuA_His1'] == 'HIS' and r['CuA_His2'] == 'HIS' and r['CuA_His3'] == 'HIS'
        and r['CuB_His1'] == 'TYR' and r['CuB_His2'] == 'THR' and r['CuB_His3'] == 'GLN']

    groups['Microbispora_nocu'] = [
        r['accession'] for r in nc
        if r.get('classification', '') == 'no_cu' and 'Microbispora' in r.get('genus', '')]

    groups['Clonorchis_Opisthorchis'] = [
        r['accession'] for r in nc
        if r.get('classification', '') == 'no_cu'
        and r.get('genus', '') in ('Clonorchis', 'Opisthorchis')]

    groups['Oomycota_nocu'] = [
        r['accession'] for r in nc
        if r.get('classification', '') == 'no_cu' and r.get('phylum', '') == 'Oomycota']

    return groups


# --- UniProt cross-refs ---

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
        print(f'    UniProt error: {e}', file=sys.stderr)
        return {}


def extract_embl_refs(entry):
    refs = []
    for xref in entry.get('uniProtKBCrossReferences', []):
        if xref.get('database') != 'EMBL':
            continue
        props = {p['key']: p['value'] for p in xref.get('properties', [])}
        pid = props.get('ProteinId', '')
        if pid in ('-', ''):
            continue
        refs.append({'nucleotide_acc': xref.get('id', ''), 'protein_id': pid,
                     'molecule_type': props.get('MoleculeType', '')})
    return refs


def step_crossrefs(accessions, out_file):
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
    with open(out_file, 'w', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fields)
        w.writeheader()
        w.writerows(results)

    with_refs = [r for r in results if r['nucleotide_acc']]
    print(f'    {len(with_refs)}/{len(results)} have nucleotide cross-refs')
    return with_refs


# --- NCBI neighbourhoods ---

def fetch_genbank(nuc_acc):
    params = urlencode({'db': 'nuccore', 'id': nuc_acc, 'rettype': 'gb',
                        'retmode': 'text', 'email': EMAIL, 'tool': 'ppo_gna'})
    for attempt in range(3):
        try:
            with urlopen(f'{NCBI_EFETCH}?{params}', timeout=120) as resp:
                return resp.read().decode('utf-8', errors='replace')
        except (HTTPError, URLError, Exception):
            if attempt < 2:
                time.sleep(2 ** attempt)
    return None


def parse_cds_features(gb_text):
    features = []
    in_cds = False
    current = {}
    cq = None
    cv = ''
    for line in gb_text.split('\n'):
        if line.startswith('     CDS '):
            if in_cds and current:
                features.append(current)
            in_cds = True
            current = {'location_raw': line[21:].strip(), 'qualifiers': {}}
            cq = None; cv = ''
        elif in_cds and line.startswith('                     /'):
            if cq:
                current['qualifiers'][cq] = cv.strip('"')
            match = re.match(r'\s+/(\w+)(?:="?(.*))?', line)
            if match:
                cq = match.group(1)
                cv = match.group(2) or ''
                if cv.endswith('"'):
                    cv = cv[:-1]
            else:
                cq = None; cv = ''
        elif in_cds and line.startswith('                     ') and cq:
            cv += ' ' + line.strip().strip('"')
        elif in_cds and not line.startswith('                     '):
            if cq:
                current['qualifiers'][cq] = cv.strip('"')
            features.append(current)
            in_cds = False
            current = {}; cq = None
    if in_cds and current:
        if cq:
            current['qualifiers'][cq] = cv.strip('"')
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
            'protein_id': q.get('protein_id', ''),
            'locus_tag': q.get('locus_tag', q.get('gene', '')),
            'product': q.get('product', 'hypothetical protein'),
        })
    result.sort(key=lambda x: x['start'])
    return result


def step_neighbourhoods(xrefs, out_file):
    limiter = RateLimiter(3)
    by_nuc = {}
    for r in xrefs:
        by_nuc.setdefault(r['nucleotide_acc'], []).append(r)
    nuc_list = list(by_nuc.keys())

    fields = ['query_accession', 'nucleotide_acc', 'offset', 'is_target',
              'gene_start', 'gene_end', 'strand', 'protein_id', 'locus_tag', 'product']
    outf = open(out_file, 'w', newline='')
    writer = csv.DictWriter(outf, delimiter='\t', fieldnames=fields)
    writer.writeheader()
    n_ok = n_fail = 0

    def rf(nuc):
        limiter.acquire()
        return nuc, fetch_genbank(nuc)

    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = {pool.submit(rf, nuc): nuc for nuc in nuc_list}
        for i, future in enumerate(as_completed(futures), 1):
            try:
                nuc_acc, gb_text = future.result()
            except Exception:
                nuc_acc = futures[future]
                gb_text = None
            entries = by_nuc[nuc_acc]
            if i % 100 == 0 or i == len(nuc_list):
                print(f'    {i}/{len(nuc_list)} ({n_ok} ok, {n_fail} fail)')
            if gb_text is None:
                n_fail += len(entries); continue
            features = parse_cds_features(gb_text)
            if not features:
                n_fail += len(entries); continue
            for entry in entries:
                pid = entry['protein_id'].split('.')[0]
                target_idx = None
                for j, feat in enumerate(features):
                    if feat['protein_id'].split('.')[0] == pid:
                        target_idx = j; break
                if target_idx is None:
                    n_fail += 1; continue
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
    print(f'    {n_ok} neighbourhoods, {n_fail} failed')
    return n_ok


# --- Pfam via ID-mapping ---

def idmap_submit(ids):
    data = urlencode({'from': 'EMBL-GenBank-DDBJ_CDS', 'to': 'UniProtKB', 'ids': ','.join(ids)}).encode()
    req = Request(f'{UNIPROT_API}/idmapping/run', data=data)
    return json.loads(urlopen(req, timeout=120).read())['jobId']


def idmap_wait(job):
    url = f'{UNIPROT_API}/idmapping/status/{job}'
    for _ in range(120):
        req = Request(url, headers={'Accept': 'application/json'})
        r = json.loads(urlopen(req, timeout=60).read())
        st = r.get('jobStatus')
        if st in ('FINISHED', None) or 'results' in r:
            return
        if st in ('RUNNING', 'NEW'):
            time.sleep(3); continue
        raise RuntimeError(f'idmapping {job}: {r}')
    raise TimeoutError(f'idmapping {job} timeout')


def idmap_fetch(job):
    url = (f'{UNIPROT_API}/idmapping/uniprotkb/results/{job}'
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
                rows.append((parts[0], parts[2]))
            elif len(parts) == 2:
                rows.append((parts[0], ''))
        m = re.search(r'<([^>]+)>;\s*rel="next"', link)
        url = m.group(1) if m else ''
        time.sleep(0.2)
    return rows


def pfam_name_lookup(pid, cache):
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


def step_pfam(nb_file, out_pfam, name_cache_file):
    nb_rows = list(csv.DictReader(open(nb_file), delimiter='\t'))
    neigh = [r for r in nb_rows if r['is_target'] != '1' and r['protein_id']]
    pids = sorted({r['protein_id'].split('.')[0] for r in neigh})
    if not pids:
        print('    No neighbour proteins to map')
        return {}
    print(f'    {len(pids)} neighbour proteins to ID-map')

    pid2pfam = defaultdict(set)
    for chunk_start in range(0, len(pids), 5000):
        chunk = pids[chunk_start:chunk_start + 5000]
        try:
            job = idmap_submit(chunk)
            print(f'    idmapping job submitted ({len(chunk)} IDs), waiting...')
            idmap_wait(job)
            rows = idmap_fetch(job)
            print(f'    {len(rows)} hits returned')
            for frm, pfam_col in rows:
                for fam in pfam_col.split(';'):
                    fam = fam.strip()
                    if fam.startswith('PF'):
                        pid2pfam[frm.split('.')[0]].add(fam)
        except Exception as e:
            print(f'    idmapping error: {e}', file=sys.stderr)

    mapped = {k: v for k, v in pid2pfam.items() if v}
    print(f'    {len(mapped)}/{len(pids)} have Pfam')

    cache = json.loads(name_cache_file.read_text()) if name_cache_file.exists() else {}
    all_fams = sorted({f for fams in mapped.values() for f in fams})
    for f in all_fams:
        pfam_name_lookup(f, cache)
    name_cache_file.write_text(json.dumps(cache, indent=1))

    with open(out_pfam, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['protein_id', 'pfam'])
        for pid in pids:
            fams = sorted(mapped.get(pid, []))
            w.writerow([pid, ';'.join(f'{f}|{cache.get(f,"")}' for f in fams)])

    return {pid: [f'{f}|{cache.get(f,"")}' for f in sorted(fams)] for pid, fams in mapped.items()}


# --- Analysis ---

def norm_product(p):
    p = re.sub(r'\s+', ' ', p.strip().lower())
    if p in ('hypothetical protein', 'uncharacterized protein',
             'unnamed protein product', 'putative protein', ''):
        return ''
    return p


def step_analysis(nb_file, pid2pfam_list, target_accs, out_dir):
    nb_rows = list(csv.DictReader(open(nb_file), delimiter='\t'))
    neigh = [r for r in nb_rows if r['is_target'] != '1']
    queries = {r['query_accession'] for r in nb_rows}
    n_total = len(target_accs)

    # Product co-occurrence
    per_q_prod = defaultdict(set)
    prod_total = Counter()
    for r in neigh:
        np_ = norm_product(r.get('product', ''))
        if not np_:
            continue
        per_q_prod[r['query_accession']].add(np_)
        prod_total[np_] += 1
    prod_qc = Counter()
    for q, prods in per_q_prod.items():
        for p in prods:
            prod_qc[p] += 1
    n_info = len(per_q_prod)

    with open(out_dir / 'cooccurrence_products.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['product', 'n_queries', 'freq_all', 'freq_with_data', 'total'])
        for p, qc in prod_qc.most_common():
            w.writerow([p, qc, f'{qc/n_total:.3f}', f'{qc/n_info:.3f}' if n_info else '0', prod_total[p]])

    # Keyword co-occurrence
    per_q_kw = defaultdict(set)
    kw_total = Counter()
    for r in neigh:
        prod = r.get('product', '').lower()
        if norm_product(prod) == '':
            continue
        toks = {w for w in re.split(r'[^a-z0-9]+', prod) if len(w) > 3 and w not in STOP_WORDS}
        for w in toks:
            per_q_kw[r['query_accession']].add(w)
            kw_total[w] += 1
    kw_qc = Counter()
    for q, kws in per_q_kw.items():
        for w in kws:
            kw_qc[w] += 1
    with open(out_dir / 'cooccurrence_keywords.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['keyword', 'n_queries', 'freq_all', 'total'])
        for kw, qc in kw_qc.most_common():
            if qc >= 2:
                w.writerow([kw, qc, f'{qc/n_total:.3f}', kw_total[kw]])

    # Pfam co-occurrence
    pid2pfam = {}
    for pid, fam_list in pid2pfam_list.items():
        pid2pfam[pid] = {f.split('|')[0] for f in fam_list}

    per_q_pfam = defaultdict(set)
    pfam_total = Counter()
    pfam_name = {}
    for r in neigh:
        pid = r.get('protein_id', '').split('.')[0]
        for fam_str in pid2pfam_list.get(pid, []):
            fid, _, nm = fam_str.partition('|')
            per_q_pfam[r['query_accession']].add(fid)
            pfam_total[fid] += 1
            if nm:
                pfam_name[fid] = nm
    pfam_qc = Counter()
    for q, fams in per_q_pfam.items():
        for f in fams:
            pfam_qc[f] += 1
    n_pfam = len(per_q_pfam)

    with open(out_dir / 'cooccurrence_pfam.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['pfam_id', 'pfam_name', 'n_queries', 'freq_all', 'freq_pfam', 'total'])
        for fid, qc in pfam_qc.most_common():
            w.writerow([fid, pfam_name.get(fid, ''), qc, f'{qc/n_total:.3f}',
                        f'{qc/n_pfam:.3f}' if n_pfam else '0', pfam_total[fid]])

    # BGC scoring
    nb_by_q = defaultdict(list)
    for r in neigh:
        if r['protein_id']:
            nb_by_q[r['query_accession']].append((int(r['offset']), r['protein_id'].split('.')[0]))

    bgc_rows = []
    for q in sorted(queries):
        fams_all = set()
        ust_offsets = []
        for off, pid in nb_by_q.get(q, []):
            fams = pid2pfam.get(pid, set())
            fams_all |= fams
            if USTYA in fams:
                ust_offsets.append(abs(off))
        n_ust = len(ust_offsets)
        markers = sum([bool(fams_all & MARKERS[k]) for k in MARKERS])
        bgc = 'UstYa_cluster' if n_ust >= 1 else ('BGC_accessory' if markers >= 2 else 'none')
        bgc_rows.append({
            'accession': q, 'n_neighbours': len(nb_by_q.get(q, [])),
            'n_ustya': n_ust, 'closest_ustya': min(ust_offsets) if ust_offsets else '',
            'has_TF': int(bool(fams_all & MARKERS['TF'])),
            'has_transporter': int(bool(fams_all & MARKERS['transporter'])),
            'has_P450': int(bool(fams_all & MARKERS['P450'])),
            'has_OMT': int(bool(fams_all & MARKERS['OMT'])),
            'bgc_context': bgc,
        })

    fields = ['accession', 'n_neighbours', 'n_ustya', 'closest_ustya',
              'has_TF', 'has_transporter', 'has_P450', 'has_OMT', 'bgc_context']
    with open(out_dir / 'bgc_scoring.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fields)
        w.writeheader()
        w.writerows(bgc_rows)

    n_ust = sum(1 for r in bgc_rows if r['bgc_context'] == 'UstYa_cluster')
    n_bgc = sum(1 for r in bgc_rows if r['bgc_context'] == 'BGC_accessory')
    n_none = sum(1 for r in bgc_rows if r['bgc_context'] == 'none')

    return {
        'n_targets': n_total,
        'n_with_data': len(queries),
        'n_ustya': n_ust,
        'n_bgc_accessory': n_bgc,
        'n_none': n_none,
        'top_pfam': [(fid, pfam_name.get(fid, ''), qc, f'{qc/n_total*100:.1f}%')
                     for fid, qc in pfam_qc.most_common(10)],
        'top_products': [(p, qc, f'{qc/n_total*100:.1f}%')
                         for p, qc in prod_qc.most_common(10)],
    }


def run_group(name, accessions):
    out_dir = WORK / 'groups' / name
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f'\n{"="*60}')
    print(f'  {name}: {len(accessions)} accessions')
    print(f'{"="*60}')

    # Step 1: cross-refs
    print('  Step 1: UniProt cross-refs')
    xrefs = step_crossrefs(accessions, out_dir / 'crossrefs.tsv')
    if not xrefs:
        print('    No cross-refs found, skipping.')
        return None

    # Step 2: neighbourhoods
    print('  Step 2: NCBI neighbourhoods')
    n_ok = step_neighbourhoods(xrefs, out_dir / 'neighbourhoods.tsv')
    if n_ok == 0:
        print('    No neighbourhoods extracted, skipping.')
        return None

    # Step 3: Pfam
    print('  Step 3: Pfam ID-mapping')
    pid2pfam = step_pfam(out_dir / 'neighbourhoods.tsv', out_dir / 'pfam_map.tsv',
                         out_dir / 'pfam_names.json')

    # Step 4: analysis
    print('  Step 4: Co-occurrence + BGC scoring')
    result = step_analysis(out_dir / 'neighbourhoods.tsv', pid2pfam, accessions, out_dir)

    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--group', help='Run only this group (by name)')
    args = ap.parse_args()

    groups = define_groups()

    if args.group:
        if args.group not in groups:
            print(f'Unknown group: {args.group}')
            print(f'Available: {", ".join(sorted(groups.keys()))}')
            sys.exit(1)
        groups = {args.group: groups[args.group]}

    results = {}
    for name, accs in sorted(groups.items()):
        results[name] = run_group(name, accs)

    # Summary table
    print(f'\n{"="*60}')
    print(f'  SUMMARY')
    print(f'{"="*60}')
    print(f'{"Group":<25s} {"N":>5s} {"Data":>5s} {"UstYa":>6s} {"BGC":>5s} {"None":>5s}')
    print('-' * 55)
    for name in sorted(results.keys()):
        r = results[name]
        if r is None:
            print(f'{name:<25s} {"?":>5s} {"0":>5s} {"-":>6s} {"-":>5s} {"-":>5s}')
        else:
            print(f'{name:<25s} {r["n_targets"]:>5d} {r["n_with_data"]:>5d} '
                  f'{r["n_ustya"]:>6d} {r["n_bgc_accessory"]:>5d} {r["n_none"]:>5d}')

    # Top Pfam per group
    print(f'\nTop Pfam families per group:')
    for name in sorted(results.keys()):
        r = results[name]
        if r is None or not r['top_pfam']:
            continue
        print(f'\n  {name}:')
        for fid, nm, qc, pct in r['top_pfam'][:5]:
            print(f'    {qc:3d} ({pct})  {fid} {nm}')


if __name__ == '__main__':
    main()
