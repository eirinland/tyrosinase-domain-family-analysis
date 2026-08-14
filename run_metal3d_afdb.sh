#!/bin/bash
#SBATCH --job-name=metal3d_afdb
#SBATCH --account=nn1003k
#SBATCH --partition=accel
#SBATCH --gpus=1
#SBATCH --mem-per-gpu=120G
#SBATCH --time=06:00:00
#SBATCH --array=0-3
#SBATCH --output=metal3d_afdb_%A_%a.log

set -euo pipefail

module load NRIS/GPU

SUBMITDIR=/cluster/home/eirinlandsem/Super_reference_pipeline
SIF=/cluster/projects/nn1003k/prog/allmetal3d/allmetal3d_gpu.sif
AFDB_DIR=$SUBMITDIR/afdb_cifs

TOTAL=136
OFFSET=31
N_TASKS=4
PER_TASK=$(( (TOTAL + N_TASKS - 1) / N_TASKS ))
START=$(( OFFSET + SLURM_ARRAY_TASK_ID * PER_TASK ))
END=$(( START + PER_TASK ))
if [ "$END" -gt "$(( OFFSET + TOTAL ))" ]; then END=$(( OFFSET + TOTAL )); fi

echo "Task $SLURM_ARRAY_TASK_ID: structures $START-$((END-1))"

apptainer exec --nv --cleanenv \
    --bind "$AFDB_DIR:$AFDB_DIR:ro" \
    --bind "$SUBMITDIR:$SUBMITDIR" \
    "$SIF" \
    python3 "$SUBMITDIR/run_metal3d_afdb.py" \
        --afdb-dir "$AFDB_DIR" \
        --af3-cu-csv "$SUBMITDIR/canonical_criteria_all.csv" \
        --output "$SUBMITDIR/metal3d_afdb_${SLURM_ARRAY_TASK_ID}.csv" \
        --start "$START" \
        --end "$END"
