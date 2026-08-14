#!/usr/bin/env bash
#SBATCH --job-name=metal3d_dct
#SBATCH --account=nn1003k
#SBATCH --partition=accel
#SBATCH --gpus=1
#SBATCH --mem-per-gpu=120G
#SBATCH --time=02:00:00
#SBATCH --output=/cluster/work/projects/nn1003k/eirin/bioinf/bioinf_redo/2_canonical_analysis/supplementary/metal3d_dct_dhica/metal3d_dct_remaining_%j.log

set -euo pipefail

module load NRIS/GPU

SUBMITDIR=/cluster/work/projects/nn1003k/eirin/bioinf/bioinf_redo/2_canonical_analysis/supplementary/metal3d_dct_dhica
SIF=/cluster/projects/nn1003k/prog/allmetal3d/allmetal3d_gpu.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
RUNNER=/cluster/work/projects/nn1003k/eirin/bioinf/bioinf_redo/3_noncanonical_analysis/metal3d/run_metal3d.py

MOUNT=/tmp/${SLURM_JOB_ID}_sqsh
STAGING=/tmp/${SLURM_JOB_ID}_cifs

mkdir -p "$MOUNT" "$STAGING"

cleanup() {
    if mountpoint -q "$MOUNT" 2>/dev/null; then
        fusermount -u "$MOUNT" || true
    fi
    rm -rf "$STAGING"
    rmdir "$MOUNT" 2>/dev/null || true
}
trap cleanup EXIT

squashfuse "$SQSH" "$MOUNT"

while read -r name; do
    src="$MOUNT/${name}.cif"
    [ -f "$src" ] && cp "$src" "$STAGING/"
done < "$SUBMITDIR/dct_remaining_list.txt"
echo "Staged $(ls "$STAGING" | wc -l) CIFs"

apptainer exec --nv --cleanenv \
    --bind "$STAGING:$STAGING:ro" \
    --bind "$SUBMITDIR:$SUBMITDIR" \
    --bind "$(dirname "$RUNNER"):$(dirname "$RUNNER"):ro" \
    "$SIF" \
    python3 "$RUNNER" \
        --cif-dir "$STAGING" \
        --output "$SUBMITDIR/metal3d_dct_remaining_results.tsv"

echo "Done: $(date)"
