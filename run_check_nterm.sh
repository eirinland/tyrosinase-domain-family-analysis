#!/bin/bash
#SBATCH --job-name=nterm_check
#SBATCH --account=nn1003k
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=nterm_check_%j.log

set -euo pipefail
SUBMITDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MNTDIR=/tmp/sqsh_$$
SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
POOLS=$SUBMITDIR/1_filtering/final_pools

module load NRIS/CPU

mkdir -p $MNTDIR
squashfuse $SQSH $MNTDIR

cd $SUBMITDIR
apptainer exec --cleanenv --bind /cluster/work --bind $MNTDIR:/mnt/models $SIF \
    python3 check_nterm.py \
        --cifs /mnt/models \
        --canonical $POOLS/canonical_accessions.csv \
        --noncanonical $POOLS/noncanonical_accessions.csv \
        --criteria $SUBMITDIR/canonical_criteria_all_ca.csv \
        --workers 16

fusermount -u $MNTDIR
rmdir $MNTDIR
