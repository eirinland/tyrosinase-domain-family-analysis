#!/bin/bash
#SBATCH --job-name=metal3d_showcase
#SBATCH --account=nn1003k
#SBATCH --partition=accel
#SBATCH --gpus=1
#SBATCH --mem-per-gpu=120G
#SBATCH --time=02:00:00
#SBATCH --output=metal3d_showcase_%j.log

set -euo pipefail

module load NRIS/GPU

SUBMITDIR=/cluster/home/eirinlandsem/Super_reference_pipeline/AF3_vs_AF2_evaluation
SIF=/cluster/projects/nn1003k/prog/allmetal3d/allmetal3d_gpu.sif

apptainer exec --nv --cleanenv \
    --bind "$SUBMITDIR:$SUBMITDIR" \
    "$SIF" \
    python3 "$SUBMITDIR/run_metal3d_showcase.py"
