#!/bin/bash
#SBATCH --job-name=sample_plddt
#SBATCH --account=nn1003k
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=sample_plddt_%j.log

set -euo pipefail

SUBMITDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
PYTHON_SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MOUNT=/tmp/plddt_mount_$$

mkdir -p "$MOUNT"
cleanup() { fusermount -u "$MOUNT" 2>/dev/null || true; rmdir "$MOUNT" 2>/dev/null || true; }
trap cleanup EXIT

squashfuse "$SQSH" "$MOUNT"
PMTYR=$(ls "$MOUNT"/B2ZB02_taxID_*_model.cif)

apptainer exec --cleanenv \
    --bind "$SUBMITDIR:$SUBMITDIR" \
    --bind "$MOUNT:$MOUNT:ro" \
    "$PYTHON_SIF" \
    python3 "$SUBMITDIR/sample_plddt.py" \
        --cifs "$MOUNT" \
        --pmtyr "$PMTYR" \
        --vectors-csv "$SUBMITDIR/position_vectors.csv" \
        --n-sample 0 \
        --out-csv "$SUBMITDIR/plddt_per_position.csv"

echo "Done: $(date)"
