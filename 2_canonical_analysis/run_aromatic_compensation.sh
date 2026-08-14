#!/bin/bash
#SBATCH --job-name=arom_comp
#SBATCH --account=nn1003k
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=arom_comp_%j.log

SUBMITDIR=$(pwd)
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MNTDIR=/tmp/sqsh_$$
SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif

module load NRIS/CPU

mkdir -p $MNTDIR
squashfuse $SQSH $MNTDIR

cd $SUBMITDIR
apptainer exec --cleanenv --bind /cluster/work --bind $MNTDIR:/mnt/models $SIF python3 test_aromatic_compensation.py /mnt/models

fusermount -u $MNTDIR
rmdir $MNTDIR
