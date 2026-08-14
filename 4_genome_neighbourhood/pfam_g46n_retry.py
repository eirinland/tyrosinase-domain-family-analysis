#!/usr/bin/env python3
"""Standalone Pfam + analysis for G46N — 1000-ID chunks with retries."""

import csv, json, re, sys, time
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.parse import urlencode

UNIPROT_API = 'https://rest.uniprot.org'
WORK = Path(__file__).parent / 'groups' / 'G46N'

USTYA = 'PF11807'
MARKERS = {
    'TF': {'PF00172', 'PF04082', 'PF11951'},
    'transporter': {'PF07690', 'PF00083', 'PF00005', 'PF00664'},
    'P450': {'PF00067'},
    'OMT': {'PF00891'},
}
STOP_WORDS = {'protein','domain','family','like','type','containing','with','terminal',
              'repeat','superfamily','related','probable','putative','predicted','involved',
              'associated','dependent','specific','that','this','from','have','more','also',
              'other','some','into','over','such','than','only','very','same','these',
              'hypothetical','uncharacterized'}


def idmap_submit(ids):
    data = urlencode({'from': 'EMBL-GenBank-DDBJ_CDS', 'to': 'UniProtKB', 'ids': ','.join(ids)}).encode()
    req = Request(f'{UNIPROT_API}/idmapping/run', data=data)
    return json.loads(urlopen(req, timeout=120).read())['jobId']


def idmap_wait(job, max_polls=120):
    url = f'{UNIPROT_API}/idmapping/status/{job}'
    for i in range(max_polls):
        req = Request(url, headers={'Accept': 'application/json'})
        r = json.loads(urlopen(req, timeout=60).read())
        st = r.get('jobStatus')
        if st in ('FINISHED', None) or 'results' in r:
            return
        if st in ('RUNNING', 'NEW'):
            time.sleep(3)
            continue
        raise RuntimeError(f'idmapping {job}: {r}')
    raise TimeoutError(f'idmapping {job} timeout after {max_polls} polls')


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


def norm_product(p):
    p = re.sub(r'\s+', ' ', p.strip().lower())
    if p in ('hypothetical protein', 'uncharacterized protein',
             'unnamed protein product', 'putative protein', ''):
        return ''
    return p


def main():
    nb_file = WORK / 'neighbourhoods.tsv'
    pfam_file = WORK / 'pfam_map.tsv'
    name_cache = WORK / 'pfam_names.json'

    nb_rows = list(csv.DictReader(open(nb_file), delimiter='\t'))
    neigh = [r for r in nb_rows if r['is_target'] != '1' and r['protein_id']]
    pids = sorted({r['protein_id'].split('.')[0] for r in neigh})
    print(f'{len(pids)} neighbour proteins to ID-map')

    CHUNK = 1000
    pid2pfam = defaultdict(set)
    n_chunks = (len(pids) + CHUNK - 1) // CHUNK

    for ci in range(n_chunks):
        chunk = pids[ci * CHUNK:(ci + 1) * CHUNK]
        for attempt in range(3):
            try:
                job = idmap_submit(chunk)
                print(f'  Chunk {ci+1}/{n_chunks} ({len(chunk)} IDs) submitted, job={job[:8]}...')
                idmap_wait(job)
                rows = idmap_fetch(job)
                print(f'  → {len(rows)} hits')
                for frm, pfam_col in rows:
                    for fam in pfam_col.split(';'):
                        fam = fam.strip()
                        if fam.startswith('PF'):
                            pid2pfam[frm.split('.')[0]].add(fam)
                break
            except Exception as e:
                print(f'  Chunk {ci+1} attempt {attempt+1} failed: {e}', file=sys.stderr)
                if attempt < 2:
                    time.sleep(10)
        else:
            print(f'  Chunk {ci+1} FAILED after 3 attempts', file=sys.stderr)

    mapped = {k: v for k, v in pid2pfam.items() if v}
    print(f'{len(mapped)}/{len(pids)} have Pfam')

    cache = json.loads(name_cache.read_text()) if name_cache.exists() else {}
    all_fams = sorted({f for fams in mapped.values() for f in fams})
    for f in all_fams:
        pfam_name_lookup(f, cache)
    name_cache.write_text(json.dumps(cache, indent=1))

    pid2pfam_list = {pid: [f'{f}|{cache.get(f,"")}' for f in sorted(fams)] for pid, fams in mapped.items()}

    with open(pfam_file, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['protein_id', 'pfam'])
        for pid in pids:
            fams = sorted(mapped.get(pid, []))
            w.writerow([pid, ';'.join(f'{f}|{cache.get(f,"")}' for f in fams)])

    # --- Analysis ---
    print('\nRunning analysis...')

    # Get target accessions
    queries = {r['query_accession'] for r in nb_rows}
    # Load canonical position vectors to get G46N group size
    pv_path = Path(__file__).parent.parent / '2_canonical_analysis' / 'position_vectors.csv'
    pv = list(csv.DictReader(open(pv_path)))
    target_accs = [r['accession'] for r in pv if r.get('Gly46') == 'N']
    n_total = len(target_accs)
    print(f'  G46N group: {n_total} accessions, {len(queries)} with neighbourhood data')

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

    with open(WORK / 'cooccurrence_products.tsv', 'w', newline='') as f:
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
    with open(WORK / 'cooccurrence_keywords.tsv', 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['keyword', 'n_queries', 'freq_all', 'total'])
        for kw, qc in kw_qc.most_common():
            if qc >= 2:
                w.writerow([kw, qc, f'{qc/n_total:.3f}', kw_total[kw]])

    # Pfam co-occurrence
    pid2pfam_sets = {}
    for pid, fam_list in pid2pfam_list.items():
        pid2pfam_sets[pid] = {f.split('|')[0] for f in fam_list}

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

    with open(WORK / 'cooccurrence_pfam.tsv', 'w', newline='') as f:
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
            fams = pid2pfam_sets.get(pid, set())
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
    with open(WORK / 'bgc_scoring.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fields)
        w.writeheader()
        w.writerows(bgc_rows)

    n_ust = sum(1 for r in bgc_rows if r['bgc_context'] == 'UstYa_cluster')
    n_bgc = sum(1 for r in bgc_rows if r['bgc_context'] == 'BGC_accessory')
    n_none = sum(1 for r in bgc_rows if r['bgc_context'] == 'none')

    print(f'\n=== G46N Results ===')
    print(f'  Total: {n_total}, With data: {len(queries)}')
    print(f'  UstYa: {n_ust}, BGC accessory: {n_bgc}, None: {n_none}')
    print(f'  BGC%: {(n_ust + n_bgc) / len(queries) * 100:.1f}% (of those with data)')
    print(f'\n  Top 10 Pfam:')
    for fid, qc in pfam_qc.most_common(10):
        print(f'    {fid} {pfam_name.get(fid, "")}: {qc} ({qc/n_total*100:.1f}%)')
    print(f'\n  Top 10 Products:')
    for p, qc in prod_qc.most_common(10):
        print(f'    {p}: {qc} ({qc/n_total*100:.1f}%)')


if __name__ == '__main__':
    main()
