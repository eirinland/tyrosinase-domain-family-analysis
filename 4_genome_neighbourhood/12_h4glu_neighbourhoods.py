"""
Step 12: Genome neighbourhoods for the four H4Glu (His4->Glu) mononuclear structures.

The 4 H4Glu accessions are oomycete (Phytophthora fragariae x3) + plant
(Sesamum latifolium x1) -- a different taxonomic story from the Ascomycota-heavy
H5Pro group. Three are DELETED from UniProt ("not part of a reference proteome");
their genomic context is recovered through UniParc -> source EMBLWGS GenBank CDS
protein IDs (resolved in the prior step). This script:

  1. efetch each GenBank protein (GenPept) -> parse /coded_by -> WGS nucleotide
     scaffold accession + CDS location
  2. efetch the scaffold GenBank record
  3. find the target CDS by protein_id and extract +/-10 flanking genes
     (reuses parse_cds_features / find_target_and_neighbours from
      3_fetch_neighbourhoods.py)

For accessions with several candidate protein IDs (same UniParc sequence in
multiple assemblies) the best-annotated scaffold (most flanking CDS) is used;
alternates are reported.

Output: h4glu_neighbourhoods.tsv  (same columns as neighbourhoods.tsv)
"""

import csv
import re
import sys
import time
from pathlib import Path
from urllib.request import urlopen
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode

sys.path.insert(0, str(Path(__file__).parent))
from importlib import import_module
_m = import_module('3_fetch_neighbourhoods')
parse_cds_features = _m.parse_cds_features
find_target_and_neighbours = _m.find_target_and_neighbours

WORK = Path(__file__).parent
OUTPUT = WORK / 'h4glu_neighbourhoods.tsv'
EMAIL = 'eirinlandsem1@gmail.com'
NCBI_EFETCH = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
FLANK = 10

# H4Glu set: accession -> taxonomy + recovered GenBank CDS protein IDs (UniParc EMBLWGS)
H4GLU = [
    {'accession': 'A0A6A3Q4B5', 'organism': 'Phytophthora fragariae',
     'phylum': 'Oomycota', 'protein_ids': ['KAE9169894.1']},
    {'accession': 'A0A6A3HLB0', 'organism': 'Phytophthora fragariae',
     'phylum': 'Oomycota', 'protein_ids': ['KAE8969384', 'KAE9176695', 'KAE9067838']},
    {'accession': 'A0A6G0QFF2', 'organism': 'Phytophthora fragariae',
     'phylum': 'Oomycota', 'protein_ids': ['KAE9285141']},
    {'accession': 'A0AAW2WSM8', 'organism': 'Sesamum latifolium',
     'phylum': 'Streptophyta', 'protein_ids': ['KAL0443686']},
]


def efetch(db, acc, rettype):
    params = urlencode({'db': db, 'id': acc, 'rettype': rettype,
                        'retmode': 'text', 'email': EMAIL, 'tool': 'ppo_gna'})
    url = f'{NCBI_EFETCH}?{params}'
    for attempt in range(3):
        try:
            with urlopen(url, timeout=120) as resp:
                return resp.read().decode('utf-8', errors='replace')
        except (HTTPError, URLError):
            if attempt < 2:
                time.sleep(2 ** attempt)
    return None


def extract_qualifier(text, name):
    """Value of /name="..." in a GenBank/GenPept block, collapsing line wraps."""
    idx = text.find(f'/{name}="')
    if idx == -1:
        return ''
    start = idx + len(name) + 3
    end = text.find('"', start)
    return re.sub(r'\s+', '', text[start:end]) if end != -1 else ''


def parse_coded_by(gp_text):
    """/coded_by -> (nuc_acc, strand, start, end). Handles complement()/join()."""
    cb = extract_qualifier(gp_text, 'coded_by')
    if not cb:
        return None
    strand = '-' if cb.startswith('complement') else '+'
    m = re.search(r'([A-Za-z]+\d+(?:\.\d+)?):', cb)
    if not m:
        return None
    nuc = m.group(1)
    nums = re.findall(r'\d+', cb.split(':', 1)[1])
    if len(nums) < 2:
        return None
    return nuc, strand, int(nums[0]), int(nums[-1])


def main():
    rows = []
    summary = []
    for rec in H4GLU:
        acc = rec['accession']
        print(f'\n=== {acc}  ({rec["organism"]}, {rec["phylum"]}) ===')
        candidates = []
        for pid in rec['protein_ids']:
            gp = efetch('protein', pid, 'gp')
            time.sleep(0.4)
            if not gp:
                print(f'  {pid}: GenPept fetch failed')
                continue
            cb = parse_coded_by(gp)
            if not cb:
                print(f'  {pid}: no /coded_by')
                continue
            nuc, strand, start, end = cb
            gb = efetch('nuccore', nuc, 'gb')
            time.sleep(0.4)
            if not gb:
                print(f'  {pid}: scaffold {nuc} fetch failed')
                continue
            feats = parse_cds_features(gb)
            target, neigh = find_target_and_neighbours(feats, pid, FLANK)
            n_n = len(neigh) - 1 if target else 0
            print(f'  {pid}: scaffold {nuc} ({len(feats)} CDS) '
                  f'-> target {"found" if target else "NOT found"}, {n_n} neighbours')
            if target:
                candidates.append((n_n, pid, nuc, neigh))

        if not candidates:
            print(f'  -> no resolvable neighbourhood for {acc}')
            summary.append((acc, rec['organism'], '', '', 0, []))
            continue

        candidates.sort(reverse=True)
        n_n, pid, nuc, neigh = candidates[0]
        alts = [c[1] for c in candidates[1:]]
        if alts:
            print(f'  -> using {pid} (scaffold {nuc}); alternates: {", ".join(alts)}')

        prods = []
        for nb in neigh:
            rows.append({
                'query_accession': acc, 'nucleotide_acc': nuc,
                'offset': nb['offset'], 'is_target': '1' if nb['is_target'] else '0',
                'gene_start': nb['start'], 'gene_end': nb['end'], 'strand': nb['strand'],
                'protein_id': nb['protein_id'], 'locus_tag': nb['locus_tag'],
                'product': nb['product'],
            })
            tag = '>>>' if nb['is_target'] else f'{nb["offset"]:+d}'
            print(f'     {tag:>4}  {nb["strand"]}  {nb["protein_id"]:<16}  {nb["product"]}')
            if not nb['is_target']:
                prods.append(nb['product'])
        summary.append((acc, rec['organism'], pid, nuc, n_n, prods))

    fields = ['query_accession', 'nucleotide_acc', 'offset', 'is_target',
              'gene_start', 'gene_end', 'strand', 'protein_id', 'locus_tag', 'product']
    with open(OUTPUT, 'w', newline='') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    print('\n' + '=' * 70)
    print('SUMMARY: H4Glu genome neighbourhoods')
    for acc, org, pid, nuc, n_n, prods in summary:
        print(f'  {acc} ({org}): {n_n} neighbours on {nuc or "?"}')
    print(f'\nOutput: {OUTPUT.name} ({len(rows)} rows)')


if __name__ == '__main__':
    main()
