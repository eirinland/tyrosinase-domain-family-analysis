#!/bin/bash
#SBATCH --job-name=extract_seqs
#SBATCH --account=nn1003k
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=extract_seqs_%j.log

set -eu

WORKDIR=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/6_domain_analysis
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MOUNT=/tmp/${USER}_cifmnt_$$

mkdir -p "$MOUNT"
squashfuse "$SQSH" "$MOUNT"
echo "Mounted squashfs at $MOUNT"
echo "Files in mount: $(ls "$MOUNT" | wc -l)"

python3 "$WORKDIR/extract_sequences.py" "$MOUNT" 500

fusermount -u "$MOUNT"
rmdir "$MOUNT"
echo "Unmounted squashfs"

# Also write a single combined FASTA for reference
cat "$WORKDIR/sequences/chunk_"*.fasta > "$WORKDIR/sequences/all_sequences.fasta"
echo "Combined FASTA: $(grep -c '^>' "$WORKDIR/sequences/all_sequences.fasta") sequences"
