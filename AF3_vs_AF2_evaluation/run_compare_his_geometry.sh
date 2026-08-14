#!/bin/bash
#SBATCH --job-name=his_geometry
#SBATCH --account=nn1003k
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=his_geometry_%j.log

set -euo pipefail

SUBMITDIR=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline
EVALDIR=$SUBMITDIR/AF3_vs_AF2_evaluation
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MNTDIR=/tmp/sqsh_hisgeo_$$
AFDBDIR=/tmp/afdb_cifs_$$

module load NRIS/CPU

mkdir -p "$MNTDIR" "$AFDBDIR"
squashfuse "$SQSH" "$MNTDIR"

cleanup() {
    fusermount -u "$MNTDIR" 2>/dev/null || true
    rm -rf "$AFDBDIR" "$MNTDIR" 2>/dev/null || true
}
trap cleanup EXIT

echo "Mounted AF3 structures: $(ls "$MNTDIR" | wc -l)"
echo "Starting His geometry comparison at $(date)"

python3 "$EVALDIR/compare_his_geometry.py" \
    --canonical-csv "$SUBMITDIR/canonical_criteria_all_ca.csv" \
    --af3-dir "$MNTDIR" \
    --afdb-dir "$AFDBDIR" \
    --output "$EVALDIR/af3_vs_afdb_his_geometry.csv" \
    --workers 8

echo ""
echo "Cleaning up AFDB downloads..."
echo "Done at $(date)"
