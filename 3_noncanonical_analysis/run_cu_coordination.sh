#!/bin/bash
#SBATCH --account=nn1003k
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --job-name=cu_coord_3A
#SBATCH --output=cu_coord_3A_%j.log

set -eu

SUBMITDIR=$(pwd)
WORKDIR=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/3_noncanonical_analysis
SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MOUNT=/tmp/cifs_$$

module load NRIS/CPU

mkdir -p "$MOUNT"
squashfuse "$SQSH" "$MOUNT"

cd "$WORKDIR"

apptainer exec --cleanenv \
  --bind "$MOUNT":/data/cifs \
  --bind "$WORKDIR":/work \
  "$SIF" \
  python3 /work/cu_coordination_detail.py \
    --cif-dir /data/cifs \
    --pool /work/nc_pool.csv \
    --out /work/cu_coordination_3A.tsv

fusermount -u "$MOUNT" 2>/dev/null || true
rmdir "$MOUNT" 2>/dev/null || true

echo "Done. Output: $WORKDIR/cu_coordination_3A.tsv"
