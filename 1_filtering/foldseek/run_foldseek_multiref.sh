#!/bin/bash
#SBATCH --job-name=foldseek_multiref
#SBATCH --account=nn1003k
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=foldseek_multiref_%j.log

SUBMITDIR=/cluster/home/eirinlandsem/Super_reference_pipeline
REFDIR=$SUBMITDIR/1_filtering/foldseek
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

# Copy reference PDBs to a single directory for createdb
cp "$REFDIR"/ref_*.pdb "$TMPDIR/refs/"
echo "References:"
ls "$TMPDIR/refs/"

echo "Creating target DB from 9 references..."
$FOLDSEEK createdb "$TMPDIR/refs" "$TMPDIR/targetdb" --threads 1

echo "Creating query DB from all structures..."
$FOLDSEEK createdb "$MNTDIR" "$TMPDIR/querydb" --threads 16

echo "Running search (all queries vs all references)..."
$FOLDSEEK search "$TMPDIR/querydb" "$TMPDIR/targetdb" "$TMPDIR/result" "$TMPDIR/tmp" \
    --threads 16 \
    -a \
    -e 10 \
    --exhaustive-search 1 \
    --max-seqs 9

echo "Converting all hits..."
$FOLDSEEK convertalis "$TMPDIR/querydb" "$TMPDIR/targetdb" "$TMPDIR/result" \
    "$REFDIR/foldseek_multiref_all.tsv" \
    --format-mode 4 \
    --format-output query,target,qstart,qend,qlen,tstart,tend,tlen,alntmscore,qtmscore,lddt,rmsd,alnlen,evalue

echo "Total hits:"
wc -l "$REFDIR/foldseek_multiref_all.tsv"

echo "Selecting best hit per query (by alntmscore)..."
python3 -c "
import csv, sys

best = {}
with open('$REFDIR/foldseek_multiref_all.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for r in reader:
        q = r['query']
        score = float(r['alntmscore'])
        if q not in best or score > float(best[q]['alntmscore']):
            best[q] = r

with open('$REFDIR/foldseek_multiref_best.tsv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=reader.fieldnames, delimiter='\t')
    writer.writeheader()
    for q in sorted(best):
        writer.writerow(best[q])

print(f'Best hits: {len(best)} queries')

# Summary: how many queries matched each reference
from collections import Counter
ref_counts = Counter(r['target'] for r in best.values())
for ref, cnt in ref_counts.most_common():
    print(f'  {ref}: {cnt}')
"

echo "Done."
