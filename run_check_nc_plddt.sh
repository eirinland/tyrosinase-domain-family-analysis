#!/bin/bash
#SBATCH --job-name=nc_plddt
#SBATCH --account=nn1003k
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=00:10:00
#SBATCH --output=/cluster/home/eirinlandsem/Super_reference_pipeline/nc_plddt_%j.log

SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MNTDIR=/tmp/sqsh_$$
SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif

module load NRIS/CPU
mkdir -p $MNTDIR && squashfuse $SQSH $MNTDIR
apptainer exec --cleanenv --bind /cluster/work --bind $MNTDIR:/mnt/models $SIF \
    python3 /cluster/home/eirinlandsem/Super_reference_pipeline/check_nc_plddt.py
fusermount -u $MNTDIR; rmdir $MNTDIR
