#!/bin/bash
#SBATCH --job-name=fs_multiref_bt
#SBATCH --account=nn1003k
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=foldseek_multiref_bt_%j.log

set -euo pipefail

SUBMITDIR=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline
REFDIR=$SUBMITDIR/1_filtering/foldseek
VALDIR=$SUBMITDIR/1_filtering/validation
FOLDSEEK=/cluster/projects/nn1003k/prog/foldseek/bin/foldseek
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MNTDIR=/tmp/sqsh_foldseek_$$
TMPDIR=/tmp/foldseek_tmp_$$

module load NRIS/CPU

mkdir -p "$MNTDIR" "$TMPDIR" "$TMPDIR/refs"
squashfuse "$SQSH" "$MNTDIR"

cleanup() {
    fusermount -u "$MNTDIR" 2>/dev/null || true
    rm -rf "$TMPDIR" "$MNTDIR" 2>/dev/null || true
}
trap cleanup EXIT

cp "$REFDIR"/ref_*.pdb "$TMPDIR/refs/"
echo "References: $(ls "$TMPDIR/refs/" | wc -l)"

echo "Creating target DB..."
$FOLDSEEK createdb "$TMPDIR/refs" "$TMPDIR/targetdb" --threads 1

echo "Creating query DB from $(ls "$MNTDIR" | wc -l) structures..."
$FOLDSEEK createdb "$MNTDIR" "$TMPDIR/querydb" --threads 16

echo "Running search with backtrace..."
$FOLDSEEK search "$TMPDIR/querydb" "$TMPDIR/targetdb" "$TMPDIR/result" "$TMPDIR/tmp" \
    --threads 16 \
    -a \
    -e 10 \
    --exhaustive-search 1 \
    --max-seqs 9

echo "Converting all hits (with cigar)..."
$FOLDSEEK convertalis "$TMPDIR/querydb" "$TMPDIR/targetdb" "$TMPDIR/result" \
    "$REFDIR/foldseek_multiref_all_bt.tsv" \
    --format-mode 4 \
    --format-output query,target,qstart,qend,qlen,tstart,tend,tlen,alntmscore,qtmscore,lddt,rmsd,alnlen,evalue,cigar

echo "Total hits: $(wc -l < "$REFDIR/foldseek_multiref_all_bt.tsv")"

echo "Selecting best hit per query..."
python3 << 'PYEOF'
import csv

REFDIR = "/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/1_filtering/foldseek"
best = {}
with open(f"{REFDIR}/foldseek_multiref_all_bt.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    fields = reader.fieldnames
    for r in reader:
        q = r["query"]
        score = float(r["alntmscore"])
        if q not in best or score > float(best[q]["alntmscore"]):
            best[q] = r

with open(f"{REFDIR}/foldseek_multiref_best_bt.tsv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for q in sorted(best):
        writer.writerow(best[q])

print(f"Best hits: {len(best)} queries")
from collections import Counter
ref_counts = Counter(r["target"] for r in best.values())
for ref, cnt in ref_counts.most_common():
    print(f"  {ref}: {cnt}")
PYEOF

echo ""
echo "Running helix coverage check on all structures..."
python3 "$VALDIR/check_core_helices.py" \
    "$REFDIR/foldseek_multiref_all_bt.tsv" \
    -o "$SUBMITDIR/1_filtering/helix_coverage.csv"

echo ""
echo "Running helix pLDDT check on full structures..."
python3 "$VALDIR/check_helix_plddt.py" \
    "$REFDIR/foldseek_multiref_all_bt.tsv" \
    --cif-dir "$MNTDIR" \
    -o "$SUBMITDIR/1_filtering/helix_plddt.csv"

echo ""
echo "Done at $(date)"
