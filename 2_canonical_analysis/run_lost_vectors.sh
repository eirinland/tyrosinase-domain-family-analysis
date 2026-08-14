#!/bin/bash
#SBATCH --job-name=lost_vec
#SBATCH --account=nn1003k
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=lost_vec_%j.log

SUBMITDIR=$(pwd)
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MNTDIR=/tmp/sqsh_$$
SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif

module load NRIS/CPU

mkdir -p $MNTDIR
squashfuse $SQSH $MNTDIR

cd $SUBMITDIR
apptainer exec --cleanenv --bind /cluster/work --bind $MNTDIR:/mnt/models $SIF \
    python3 extract_lost_vectors.py \
        --cifs /mnt/models \
        --lost alignmentfree_test.tsv \
        --output lost_vectors.csv \
        --workers 8

fusermount -u $MNTDIR
rmdir $MNTDIR
