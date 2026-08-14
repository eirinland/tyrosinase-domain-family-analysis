#!/bin/bash
#SBATCH --account=nn1003k
#SBATCH --job-name=ppo_validation
#SBATCH --ntasks-per-node=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=validation_%j.log

set -eu
SUBMITDIR=$(pwd)
cd "$SUBMITDIR"

module load NRIS/CPU

SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MNTDIR=/tmp/sqsh_validation_$$
REFDIR=~/Super_reference_pipeline/1_filtering/foldseek
DOMAIN_CSV=~/Super_reference_pipeline/1_filtering/foldseek/ppo_domain_assignment_multiref.csv
CANON_SCRIPT=~/Super_reference_pipeline/check_canonical_criteria.py
FOLDSEEK=/cluster/projects/nn1003k/prog/foldseek/bin/foldseek

ACCS="A0A075DN54 A0A0K1ZP03 A0A1S9DK56 A0A261GRE4 A0A8D3X086 C7FF05 G2QLD3 P43309 P43311 Q08303 Q2T7K1 Q2UNF9 Q83WS2 Q93HL2 Q2UP46 P83040"

# 1. Extract CIFs for validation set
echo "=== Extracting CIFs ==="
mkdir -p cifs_full cifs_trimmed "$MNTDIR"
squashfuse "$SQSH" "$MNTDIR"
for acc in $ACCS; do
    cp "$MNTDIR"/${acc}_*.cif cifs_full/ 2>/dev/null || echo "WARN: $acc not found"
done
echo "Extracted: $(ls cifs_full/ | wc -l) CIFs"

# 2. Trim to PPO core
echo "=== Trimming to PPO core ==="
ACCS_CSV="A0A075DN54,A0A0K1ZP03,A0A1S9DK56,A0A261GRE4,A0A8D3X086,C7FF05,G2QLD3,P43309,P43311,Q08303,Q2T7K1,Q2UNF9,Q83WS2,Q93HL2,Q2UP46,P83040"
python3 trim_to_ppo_core.py \
    --domain-csv "$DOMAIN_CSV" \
    --cif-dir cifs_full \
    --out-dir cifs_trimmed \
    --accessions "$ACCS_CSV"

# 3. Foldseek pass 2 on trimmed structures
echo "=== Foldseek pass 2 (trimmed) ==="
mkdir -p foldseek_tmp
$FOLDSEEK createdb cifs_trimmed foldseek_tmp/queryDB
$FOLDSEEK createdb "$REFDIR"/ref_*.pdb foldseek_tmp/targetDB
$FOLDSEEK search foldseek_tmp/queryDB foldseek_tmp/targetDB foldseek_tmp/alnDB foldseek_tmp/tmp \
    --max-seqs 9 --alignment-type 2 -e 100 -a
$FOLDSEEK convertalis foldseek_tmp/queryDB foldseek_tmp/targetDB foldseek_tmp/alnDB \
    foldseek_pass2.tsv \
    --format-output "query,target,qstart,qend,qlen,tstart,tend,tlen,alntmscore,qtmscore,lddt,rmsd,alnlen,evalue"
echo "Foldseek pass 2 hits: $(wc -l < foldseek_pass2.tsv)"

# 4. Canonical criteria on trimmed structures
echo "=== Canonical criteria (trimmed) ==="
python3 "$CANON_SCRIPT" --cifs cifs_trimmed --output canonical_trimmed.csv --workers 4

# 5. Summary
echo "=== Building summary ==="
python3 - << PYEOF
import csv

accs = "A0A075DN54 A0A0K1ZP03 A0A1S9DK56 A0A261GRE4 A0A8D3X086 C7FF05 G2QLD3 P43309 P43311 Q08303 Q2T7K1 Q2UNF9 Q83WS2 Q93HL2 Q2UP46 P83040".split()

fs2 = {}
with open("foldseek_pass2.tsv") as f:
    for line in f:
        if line.startswith("query"):
            continue
        parts = line.strip().split("\t")
        q = parts[0].split("_taxID_")[0]
        qtm = float(parts[9])
        if q not in fs2 or qtm > float(fs2[q]["qtmscore"]):
            fs2[q] = {"qtmscore": parts[9], "alnlen": parts[12], "target": parts[1], "alntmscore": parts[8]}

canon_t = {}
with open("canonical_trimmed.csv") as f:
    for r in csv.DictReader(f):
        acc = r["accession"].split("_taxID_")[0]
        canon_t[acc] = r

canon_f = {}
with open("../../canonical_criteria_all.csv") as f:
    for r in csv.DictReader(f):
        if r["accession"] in accs:
            canon_f[r["accession"]] = r

dom = {}
with open("../../1_filtering/foldseek/ppo_domain_assignment_multiref.csv") as f:
    for r in csv.DictReader(f):
        if r["accession"] in accs:
            dom[r["accession"]] = r

with open("validation_summary.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["accession", "ppo_range", "full_canon", "full_CuCu", "full_nHis", "full_pLDDT",
                "trim_canon", "trim_CuCu", "trim_nHis", "trim_pLDDT",
                "fs2_alnlen", "fs2_qtmscore", "fs2_best_ref"])
    for acc in accs:
        cf = canon_f.get(acc, {})
        ct = canon_t.get(acc, {})
        f2 = fs2.get(acc, {})
        d = dom.get(acc, {})
        w.writerow([
            acc, d.get("ppo_range", ""),
            cf.get("canonical", ""), cf.get("cu_dist", ""), cf.get("n_coord_his", ""), cf.get("min_plddt", ""),
            ct.get("canonical", ""), ct.get("cu_dist", ""), ct.get("n_coord_his", ""), ct.get("min_plddt", ""),
            f2.get("alnlen", ""), f2.get("qtmscore", ""), f2.get("target", ""),
        ])

print("Wrote validation_summary.csv")
print()
print("%-14s %-10s  %5s %6s  %5s %6s  %6s %8s %s" % (
    "accession", "ppo_range", "fCanon", "fCuCu", "tCanon", "tCuCu", "fs2_al", "fs2_qtm", "fs2_ref"))
print("-" * 100)
for acc in accs:
    cf = canon_f.get(acc, {})
    ct = canon_t.get(acc, {})
    f2 = fs2.get(acc, {})
    d = dom.get(acc, {})
    print("%-14s %-10s  %5s %6s  %5s %6s  %6s %8s %s" % (
        acc, d.get("ppo_range", ""),
        cf.get("canonical", "?"), cf.get("cu_dist", "?"),
        ct.get("canonical", "?"), ct.get("cu_dist", "?"),
        f2.get("alnlen", "?"), f2.get("qtmscore", "?"), f2.get("target", "?"),
    ))
PYEOF

fusermount -u "$MNTDIR" && rmdir "$MNTDIR"
echo "=== Done ==="
