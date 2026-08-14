#!/bin/bash
#SBATCH --job-name=seed_var
#SBATCH --account=nn1003k
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --array=0-8

set -eu

SUBMITDIR="/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/1_filtering/seed_variability"
SQSH_DIR="/cluster/work/projects/nn1003k/eirin/tmp/af3_seeds"
OUTDIR="${SUBMITDIR}/results"
SCRIPT="${SUBMITDIR}/extract_seed_geometry.py"
SIF="/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif"

module load NRIS/CPU

mkdir -p "$OUTDIR"

IDX=$(printf '%03d' $SLURM_ARRAY_TASK_ID)
SQSH="${SQSH_DIR}/results_${IDX}.sqsh"
MOUNT="/tmp/${USER}/seed_var_${IDX}"

mkdir -p "$MOUNT"
squashfuse "$SQSH" "$MOUNT"

trap "fusermount -u $MOUNT 2>/dev/null; rmdir $MOUNT 2>/dev/null" EXIT

apptainer exec --cleanenv \
    --bind "$MOUNT:$MOUNT:ro" \
    --bind "$OUTDIR:$OUTDIR" \
    --bind "$SCRIPT:$SCRIPT:ro" \
    "$SIF" \
    python3 "$SCRIPT" "$MOUNT" "${OUTDIR}/seed_geometry_${IDX}.tsv"

echo "Done: task $SLURM_ARRAY_TASK_ID"
