#!/usr/bin/env bash
#SBATCH --job-name=m3d_nc_rem
#SBATCH --account=nn1003k
#SBATCH --partition=accel
#SBATCH --gpus=1
#SBATCH --mem-per-gpu=120G
#SBATCH --time=06:00:00
#SBATCH --array=0-9%5
#SBATCH --output=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/metal3d_nc_remaining/logs/metal3d_%A_%a.log

set -euo pipefail

module load NRIS/GPU

SUBMITDIR=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/metal3d_nc_remaining
SIF=/cluster/projects/nn1003k/prog/allmetal3d/allmetal3d_gpu.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh

MOUNT=/tmp/${SLURM_JOB_ID}_sqsh
STAGING=/tmp/${SLURM_JOB_ID}_cifs

mkdir -p "$MOUNT" "$STAGING" "$SUBMITDIR/logs" "$SUBMITDIR/results"

cleanup() {
    if mountpoint -q "$MOUNT" 2>/dev/null; then
        fusermount -u "$MOUNT" || true
    fi
    rm -rf "$STAGING"
    rmdir "$MOUNT" 2>/dev/null || true
}
trap cleanup EXIT

TOTAL=$(wc -l < "$SUBMITDIR/structure_list.txt")
N_TASKS=10
PER_TASK=$(( (TOTAL + N_TASKS - 1) / N_TASKS ))
START=$(( SLURM_ARRAY_TASK_ID * PER_TASK ))
END=$(( START + PER_TASK ))
if [ "$END" -gt "$TOTAL" ]; then END=$TOTAL; fi

echo "Task $SLURM_ARRAY_TASK_ID: structures $START-$((END-1)) of $TOTAL"
echo "Started: $(date)"

squashfuse "$SQSH" "$MOUNT"

CHUNK_LIST=$(sed -n "$((START+1)),$((END))p" "$SUBMITDIR/structure_list.txt")
for name in $CHUNK_LIST; do
    src="$MOUNT/${name}.cif"
    [ -f "$src" ] && cp "$src" "$STAGING/"
done
echo "Staged $(ls "$STAGING" | wc -l) CIFs"

apptainer exec --nv --cleanenv \
    --bind "$STAGING:$STAGING:ro" \
    --bind "$SUBMITDIR:$SUBMITDIR" \
    "$SIF" \
    python3 "$SUBMITDIR/run_metal3d.py" \
        --cif-dir "$STAGING" \
        --output "$SUBMITDIR/results/metal3d_remaining_${SLURM_ARRAY_TASK_ID}.tsv"

echo "Done: $(date)"
