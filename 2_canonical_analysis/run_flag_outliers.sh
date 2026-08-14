#!/bin/bash
#SBATCH --job-name=flag_outliers
#SBATCH --account=nn1003k
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=flag_outliers_%j.log

set -euo pipefail

SUBMITDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
PYTHON_SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MOUNT=/tmp/outlier_mount_$$
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "$MOUNT"
cleanup() { fusermount -u "$MOUNT" 2>/dev/null || true; rmdir "$MOUNT" 2>/dev/null || true; }
trap cleanup EXIT

echo "Mounting squashFS..."
squashfuse "$SQSH" "$MOUNT"
echo "Mounted. Files: $(ls "$MOUNT" | wc -l)"

PMTYR=$(ls "$MOUNT"/B2ZB02_taxID_*_model.cif)
echo "PmTYR ref: $PMTYR"

echo "Flagging outliers..."
apptainer exec --cleanenv \
    --bind "$SUBMITDIR:$SUBMITDIR" \
    --bind "$MOUNT:$MOUNT:ro" \
    --bind /cluster/work/projects/nn1003k/eirin:/cluster/work/projects/nn1003k/eirin \
    "$PYTHON_SIF" \
    python3 "$SUBMITDIR/flag_outliers.py" \
        --cifs "$MOUNT" \
        --pmtyr "$PMTYR" \
        --vectors-csv "$SUBMITDIR/position_vectors.csv" \
        --output "$SUBMITDIR/structure_outliers.csv" \
        --min-group 2 \
        --workers "$THREADS"

echo ""
echo "Done: $(date)"
