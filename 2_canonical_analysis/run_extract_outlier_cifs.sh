#!/bin/bash
#SBATCH --job-name=extract_outlier_cifs
#SBATCH --account=nn1003k
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=extract_outlier_cifs_%j.log

set -euo pipefail

SUBMITDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
PYTHON_SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MOUNT=/tmp/extract_cifs_mount_$$
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "$MOUNT"
cleanup() { fusermount -u "$MOUNT" 2>/dev/null || true; rmdir "$MOUNT" 2>/dev/null || true; }
trap cleanup EXIT

echo "Mounting squashFS..."
squashfuse "$SQSH" "$MOUNT"
echo "Mounted. Files: $(ls "$MOUNT" | wc -l)"

PMTYR=$(ls "$MOUNT"/B2ZB02_taxID_*_model.cif)
echo "PmTYR ref: $PMTYR"

# Select groups for inspection:
# 1. V-F-W-E-N-L-T-S-F-L--      largest group (n=1041), 27 outliers
# 2. S-F-W-E-N-I-V-S-F-H--      2nd largest with outliers (n=607), 35 outliers
# 3. V-F-W-E-N-S-P-A-F-L--      high worst RMSD (1.659), 15 outliers
# 4. C-F-F-E-G-I-F-A-F-H-C      thioether group (n=203), 7 outliers
# 5. Y-F-W-E-N-R-G-S-F-H--      moderate size (n=286), 19 outliers, worst=1.107

echo "Extracting active-site CIFs for selected groups..."
apptainer exec --cleanenv \
    --bind "$SUBMITDIR:$SUBMITDIR" \
    --bind "$MOUNT:$MOUNT:ro" \
    --bind /cluster/work/projects/nn1003k/eirin:/cluster/work/projects/nn1003k/eirin \
    "$PYTHON_SIF" \
    python3 "$SUBMITDIR/extract_outlier_cifs.py" \
        --cifs "$MOUNT" \
        --pmtyr "$PMTYR" \
        --outliers-csv "$SUBMITDIR/structure_outliers.csv" \
        --output-dir "$SUBMITDIR/outlier_inspection" \
        --workers "$THREADS" \
        --vectors \
            "V-F-W-E-N-L-T-S-F-L--" \
            "S-F-W-E-N-I-V-S-F-H--" \
            "V-F-W-E-N-S-P-A-F-L--" \
            "C-F-F-E-G-I-F-A-F-H-C" \
            "Y-F-W-E-N-R-G-S-F-H--"

echo ""
echo "Done: $(date)"
