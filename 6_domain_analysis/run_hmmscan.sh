#!/bin/bash
#SBATCH --job-name=hmmscan
#SBATCH --account=nn1003k
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=hmmscan_%A_%a.log
#SBATCH --array=0-64

set -eu

WORKDIR=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/6_domain_analysis
HMMSCAN=/cluster/work/projects/nn1003k/eirin/tools/hmmer3/bin/hmmscan
PFAM=/cluster/work/projects/nn1003k/eirin/databases/pfam/Pfam-A.hmm

CHUNK=$(printf "chunk_%04d.fasta" $SLURM_ARRAY_TASK_ID)
INPUT="$WORKDIR/sequences/$CHUNK"

if [ ! -f "$INPUT" ]; then
    echo "No such chunk: $CHUNK — skipping"
    exit 0
fi

OUTDIR="$WORKDIR/hmmscan_results"
mkdir -p "$OUTDIR"

DOMTBL="$OUTDIR/domtbl_$(printf '%04d' $SLURM_ARRAY_TASK_ID).tsv"

echo "Running hmmscan on $CHUNK ($(grep -c '^>' "$INPUT") sequences)..."
$HMMSCAN --cpu 4 --domtblout "$DOMTBL" --noali -E 1e-5 --domE 1e-3 \
    "$PFAM" "$INPUT" > /dev/null

echo "Done: $(wc -l < "$DOMTBL") lines in $DOMTBL"
