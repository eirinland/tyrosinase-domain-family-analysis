import csv, math, sys
from collections import Counter, defaultdict
from pathlib import Path

if len(sys.argv) < 2:
    print("Usage: python3 make_node_table.py <tm_label>  (e.g. tm05, tm07, tm08, tm09)")
    sys.exit(1)

tm_label = sys.argv[1]

BASE = Path('/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/foldseek/pools')
OUTDIR = BASE / ('all_vs_all_' + tm_label)

if tm_label == 'tm08':
    CLUSTER_TSV = BASE / 'results' / 'cluster_cluster.tsv'
else:
    CLUSTER_TSV = BASE / tm_label / 'results' / 'cluster_cluster.tsv'

REP_CSV = OUTDIR / 'rep_accessions.csv'
TAXONOMY = Path('/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/taxonomy_lookup.csv')
AGGER_XLSX = Path('/cluster/work/projects/nn1003k/eirin/bioinf/Agger_sequences_and_groups.xlsx')
THIOETHER_TSV = BASE / 'thioether_check.tsv'

TYPE1_GROUPS = set('abcdefgh')
TYPE2_GROUPS = set('ijkl')
KINGDOM_COLS = ['Fungi', 'Metazoa', 'Bacteria', 'Viridiplantae', 'Other_Eukaryota', 'Unknown']
KINGDOM_MAP = {
    'Animals': 'Metazoa', 'Plants': 'Viridiplantae', 'Oomycota': 'Other_Eukaryota',
    'Archaea': 'Bacteria', 'Fungi': 'Fungi', 'Metazoa': 'Metazoa',
    'Bacteria': 'Bacteria', 'Viridiplantae': 'Viridiplantae',
}

def map_kingdom(k):
    if k in KINGDOM_MAP: return KINGDOM_MAP[k]
    if k in ('?', '', 'Unknown', None): return 'Unknown'
    return 'Other_Eukaryota'

def acc_from_member(m):
    s = m.replace('.cif', '').replace('_model_A', '_model').replace('_model', '')
    return s.split('_taxID_')[0] if '_taxID_' in s else s

# --- Load taxonomy ---
tax = {}
with open(str(TAXONOMY)) as f:
    for row in csv.DictReader(f):
        tax[row['accession']] = {'kingdom': row.get('kingdom', 'Unknown'), 'phylum': row.get('phylum', 'Unknown')}
print("Taxonomy loaded: {} accessions".format(len(tax)))

# --- Load Agger groups ---
agger = {}
try:
    import openpyxl
    wb = openpyxl.load_workbook(str(AGGER_XLSX), read_only=True)
    ws = wb.active
    for row in ws.iter_rows(min_row=2, values_only=True):
        acc, group = row[0], row[1]
        if acc and group: agger[acc] = group.strip().lower()
    wb.close()
    print("Agger groups loaded: {} sequences".format(len(agger)))
except Exception as e:
    print("Warning: could not load Agger Excel: {}".format(e))

# --- Load thioether ---
thioether = {}
with open(str(THIOETHER_TSV)) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        te = row.get('thioether', '')
        thioether[row['accession']] = 1 if te == 'C' else 0
print("Thioether data loaded: {} accessions".format(len(thioether)))

# --- Load clusters at this TM ---
clusters = defaultdict(list)
with open(str(CLUSTER_TSV)) as f:
    for line in f:
        rep, mem = line.strip().split('\t')
        clusters[rep].append(mem)
print("{} clusters: {}".format(tm_label, len(clusters)))

# --- Load rep list (network nodes) ---
reps = set()
with open(str(REP_CSV)) as f:
    for row in csv.DictReader(f):
        reps.add(row['sequence_id'])
print("Network nodes ({} reps): {}".format(tm_label, len(reps)))

# For each rep, get its cluster members
# Reps ARE cluster reps at this TM, so members = clusters[rep]
rep_cluster_members = {}
for rep in reps:
    cif_name = rep + '.cif' if not rep.endswith('.cif') else rep
    bare = rep.replace('.cif', '')
    if cif_name in clusters:
        rep_cluster_members[bare] = clusters[cif_name]
    elif bare in clusters:
        rep_cluster_members[bare] = clusters[bare]
    else:
        rep_cluster_members[bare] = [bare + '.cif']

# --- Build tables ---
node_rows, tax_rows, tax_cluster_rows = [], [], []

for rep in sorted(rep_cluster_members.keys(), key=lambda r: -len(rep_cluster_members.get(r, []))):
    members = rep_cluster_members.get(rep, [])
    size = len(members)
    rep_acc = acc_from_member(rep)

    kingdoms = Counter()
    cys_count = cys_total = 0
    agger_counts = Counter()

    for m in members:
        acc = acc_from_member(m)
        k = map_kingdom(tax.get(acc, {}).get('kingdom', 'Unknown'))
        kingdoms[k] += 1
        if acc in thioether:
            cys_total += 1
            cys_count += thioether[acc]
        if acc in agger:
            agger_counts[agger[acc]] += 1

    present_groups = set(agger_counts.keys())
    if present_groups & TYPE2_GROUPS: atype = 'type2'
    elif present_groups & TYPE1_GROUPS: atype = 'type1'
    elif sum(agger_counts.values()) > 0: atype = 'mixed'
    else: atype = 'unassigned'

    cys_frac = cys_count / cys_total if cys_total > 0 else 0.0

    node_row = {
        'cluster_rep': rep, 'cluster_size': size,
        'log_cluster_size': '{:.4f}'.format(math.log10(max(size, 1))),
        'agger_type': atype, 'cys_in_shell': '{:.4f}'.format(cys_frac),
        'n_agger_members': sum(agger_counts.values()),
    }
    for g in 'abcdefghijkl':
        node_row['group_{}'.format(g)] = agger_counts.get(g, 0)
    node_rows.append(node_row)

    rep_tax = tax.get(rep_acc, {})
    rep_k = rep_tax.get('kingdom', 'Unknown') or 'Unknown'
    tax_rows.append({'name': rep, 'kingdom': rep_k, 'superkingdom': rep_k,
                     'phylum': rep_tax.get('phylum', 'Unknown') or 'Unknown'})

    total = sum(kingdoms.values()) or 1
    majority = kingdoms.most_common(1)[0][0] if kingdoms else 'Unknown'
    tc_row = {'name': rep, 'majority_kingdom': majority}
    for col in KINGDOM_COLS:
        tc_row['frac_{}'.format(col)] = '{:.4f}'.format(kingdoms.get(col, 0) / total)
    tax_cluster_rows.append(tc_row)

node_fields = ['cluster_rep', 'cluster_size', 'log_cluster_size', 'agger_type',
               'cys_in_shell', 'n_agger_members'] + ['group_{}'.format(g) for g in 'abcdefghijkl']
with open(str(OUTDIR / 'node_table.tsv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=node_fields, delimiter='\t')
    w.writeheader()
    w.writerows(node_rows)
print("Written {} rows to node_table.tsv".format(len(node_rows)))

with open(str(OUTDIR / 'node_taxonomy.tsv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['name', 'kingdom', 'superkingdom', 'phylum'], delimiter='\t')
    w.writeheader()
    w.writerows(tax_rows)
print("Written {} rows to node_taxonomy.tsv".format(len(tax_rows)))

tc_fields = ['name', 'majority_kingdom'] + ['frac_{}'.format(c) for c in KINGDOM_COLS]
with open(str(OUTDIR / 'node_taxonomy_cluster.tsv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=tc_fields, delimiter='\t')
    w.writeheader()
    w.writerows(tax_cluster_rows)
print("Written {} rows to node_taxonomy_cluster.tsv".format(len(tax_cluster_rows)))
