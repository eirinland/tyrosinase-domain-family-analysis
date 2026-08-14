#!/bin/bash
#SBATCH --job-name=canon_check
#SBATCH --account=nn1003k
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=canonical_check_%j.log

SUBMITDIR=$(pwd)
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MNTDIR=/tmp/sqsh_canon_$$

module load NRIS/CPU
module load Python/3.11.5-GCCcore-13.2.0

mkdir -p "$MNTDIR"
squashfuse "$SQSH" "$MNTDIR"

python3 "$SUBMITDIR/check_canonical_criteria.py" \
    --cifs "$MNTDIR" \
    --output "$SUBMITDIR/canonical_criteria_all.csv" \
    --workers 16

fusermount -u "$MNTDIR" 2>/dev/null || true
rmdir "$MNTDIR" 2>/dev/null || true
