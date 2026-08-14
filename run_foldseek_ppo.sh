#!/bin/bash
#SBATCH --job-name=foldseek_ppo
#SBATCH --account=nn1003k
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=foldseek_ppo_%j.log

SUBMITDIR=$(pwd)
FOLDSEEK=/cluster/projects/nn1003k/prog/foldseek/bin/foldseek
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MNTDIR=/tmp/sqsh_foldseek_$$
TMPDIR=/tmp/foldseek_tmp_$$

module load NRIS/CPU

mkdir -p "$MNTDIR" "$TMPDIR"
squashfuse "$SQSH" "$MNTDIR"

# Extract PmTYR as target
PMTYR=$(ls "$MNTDIR"/B2ZB02_taxID_*.cif | head -1)
echo "PmTYR reference: $PMTYR"

# Create foldseek databases
echo "Creating target DB..."
$FOLDSEEK createdb "$PMTYR" "$TMPDIR/targetdb" --threads 1

echo "Creating query DB..."
$FOLDSEEK createdb "$MNTDIR" "$TMPDIR/querydb" --threads 16

echo "Running search..."
$FOLDSEEK search "$TMPDIR/querydb" "$TMPDIR/targetdb" "$TMPDIR/result" "$TMPDIR/tmp" \
    --threads 16 \
    -a \
    -e 10 \
    --exhaustive-search 1 \
    --max-seqs 1

echo "Converting results..."
$FOLDSEEK convertalis "$TMPDIR/querydb" "$TMPDIR/targetdb" "$TMPDIR/result" \
    "$SUBMITDIR/foldseek_ppo_alignment.tsv" \
    --format-mode 4 \
    --format-output query,target,qstart,qend,qlen,tstart,tend,tlen,alntmscore,qtmscore,lddt,rmsd,alnlen,evalue

echo "Lines in output:"
wc -l "$SUBMITDIR/foldseek_ppo_alignment.tsv"

echo "Sample:"
head -5 "$SUBMITDIR/foldseek_ppo_alignment.tsv"

fusermount -u "$MNTDIR" 2>/dev/null || true
rm -rf "$TMPDIR" "$MNTDIR" 2>/dev/null || true
