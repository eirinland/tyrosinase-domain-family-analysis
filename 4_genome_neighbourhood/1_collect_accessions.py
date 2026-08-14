"""
Collect all target accessions for genome neighbourhood analysis.
Groups:
  - 4 candidate pools from canonical analysis (oMP, oAPO, DCT/DHICA, hemocyanin)
  - 85 characterized PPOs
  - 12 showcase structures from noncanonical analysis
Output: target_accessions.tsv (accession, group, subgroup)
"""

import csv
from pathlib import Path

BASE = Path(__file__).parent.parent
OUT = Path(__file__).parent / 'target_accessions.tsv'

# --- Candidate pools from position_vectors.csv ---
with open(BASE / '2_canonical_analysis' / 'position_vectors.csv') as f:
    pv_rows = list(csv.DictReader(f))

entries = {}

for r in pv_rows:
    acc = r['accession']
    if r.get('Arg209') == 'Y':
        entries.setdefault(acc, []).append(('candidate_pool', 'oMP'))
    if r.get('Gly46') == 'N':
        entries.setdefault(acc, []).append(('candidate_pool', 'oAPO'))
    if r.get('His230') == 'L':
        entries.setdefault(acc, []).append(('candidate_pool', 'DCT_DHICA'))
    if r.get('Gly46') == 'E':
        entries.setdefault(acc, []).append(('candidate_pool', 'hemocyanin'))

# --- 85 characterized PPOs ---
CHARACTERIZED = {
    'A0A0K1ZP03': 'TYR', 'A0A1S9DK56': 'TYR', 'A0A261GRE4': 'TYR',
    'A0A261GVB1': 'TYR', 'A0A8D3X086': 'TYR', 'A0AAJ6N653': 'TYR',
    'B2ZB02': 'TYR', 'B8NM74': 'TYR', 'C0LU17': 'TYR', 'C7FF04': 'TYR',
    'C7FF05': 'TYR', 'O42713': 'TYR', 'P00440': 'TYR', 'P06845': 'TYR',
    'P07524': 'TYR', 'P11344': 'TYR', 'P14679': 'TYR', 'P33180': 'TYR',
    'P54834': 'TYR', 'P55023': 'TYR', 'P55024': 'TYR', 'P55025': 'TYR',
    'P55033': 'TYR', 'Q00024': 'TYR', 'Q00234': 'TYR', 'Q04604': 'TYR',
    'Q08303': 'TYR', 'Q0MVP0': 'TYR', 'Q2T7K1': 'TYR', 'Q83WS2': 'TYR',
    'Q8MIU0': 'TYR', 'Q92396': 'TYR', 'Q93HL2': 'TYR', 'Q9BDE0': 'TYR',
    'A0A8E6HRF2': 'TYR', 'A0A0K0NPU9': 'TYR', 'A0A238GSS3': 'TYR',
    'A0A238GT27': 'TYR', 'L0D705': 'TYR', 'A0A9D5H103': 'TYR',
    'B3WFP2': 'TYR', 'K0B2W2': 'TYR', 'Q5VM57': 'TYR', 'P43309': 'TYR',
    'Q6UIL3': 'TYR',
    'O81103': 'CaOx', 'P43311': 'CaOx', 'Q06355': 'CaOx',
    'Q9MB14': 'CaOx', 'Q9ZP19': 'CaOx', 'I7HUF2': 'CaOx',
    'A0A075DN54': 'AUS', 'Q9FRX6': 'AUS', 'A0A2R6XDI4': 'AUS',
    'O57405': 'DHICA_ox', 'P07147': 'DHICA_ox', 'P17643': 'DHICA_ox',
    'P55027': 'DHICA_ox', 'P55028': 'DHICA_ox', 'Q2VPW6': 'DHICA_ox',
    'Q8WN57': 'DHICA_ox',
    'B1VTI5': 'oAPO', 'D6RTB9': 'oAPO',
    'A0A9P1ME48': 'oMP', 'G2Q526': 'oMP', 'G2QC95': 'oMP',
    'G2QLD3': 'oMP', 'Q2GZJ4': 'oMP', 'Q2H7I7': 'oMP', 'Q2UNF9': 'oMP',
    'O93505': 'DCT', 'P29812': 'DCT', 'P40126': 'DCT',
    'Q4R1H1': 'DCT', 'Q95119': 'DCT',
    'A0A336U966': 'biosynthetic', 'A0A8J9RRY2': 'biosynthetic',
    'P0DUQ0': 'biosynthetic', 'Q0CRX0': 'biosynthetic',
    'Q5AUW8': 'biosynthetic', 'Q5BGU9': 'biosynthetic',
    'P80960': 'hemocyanin', 'P83040': 'hemocyanin',
    'P56823': 'hemocyanin', 'P56826': 'hemocyanin',
}

for acc, act in CHARACTERIZED.items():
    entries.setdefault(acc, []).append(('characterized', act))

# --- Showcase structures ---
SHOWCASE = {
    'A0A0L0VSW2': 'binuclear_TyrCuB',
    'A0AA39K704': 'binuclear_TyrCuA',
    'A0A9Q3GXD1': 'binuclear_AspCuB',
    'A0A168EXD0': 'mononuclear_ProCuB',
    'A0AAE0D8U4': 'mononuclear_GlnCuB',
    'A0A6A3MUA3': 'mononuclear_CuA_loss',
    'A0AAN9YIB4': 'mononuclear_LeuCuA',
    'A0A0V1G4S3': 'mononuclear_CuB_loss',
    'A0ABV9I266': 'mononuclear_AspCuB_Zn',
    'A0A8H9LF69': 'neither_Microbispora',
    'H3GEM4':     'neither_oomycete',
    'H2KPL1':     'neither_liverfluke',
}

for acc, label in SHOWCASE.items():
    entries.setdefault(acc, []).append(('showcase', label))

# --- Write output ---
rows_out = []
for acc, groups in sorted(entries.items()):
    for group, subgroup in groups:
        rows_out.append({'accession': acc, 'group': group, 'subgroup': subgroup})

with open(OUT, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['accession', 'group', 'subgroup'], delimiter='\t')
    w.writeheader()
    w.writerows(rows_out)

unique_accs = len(entries)
print(f'Written {len(rows_out)} rows ({unique_accs} unique accessions) to {OUT.name}')
for g in ['candidate_pool', 'characterized', 'showcase']:
    accs_in_group = {a for a, gs in entries.items() if any(gg == g for gg, _ in gs)}
    print(f'  {g}: {len(accs_in_group)}')
