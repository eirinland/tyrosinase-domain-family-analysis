#!/bin/bash
#SBATCH --job-name=chainsaw
#SBATCH --account=nn1003k
#SBATCH --array=0-7
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=chainsaw_%A_%a.log

SUBMITDIR=$(pwd)
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MNTDIR=/tmp/sqsh_chainsaw_${SLURM_ARRAY_TASK_ID}_$$

module load NRIS/CPU
module load Python/3.11.5-GCCcore-13.2.0
source /cluster/work/projects/nn1003k/eirin/bioinf/chainsaw_env/bin/activate

CHUNK=$(printf "acc_chunk_%02d.txt" $SLURM_ARRAY_TASK_ID)

mkdir -p $MNTDIR
squashfuse $SQSH $MNTDIR

cd $SUBMITDIR
python3 run_chainsaw_batch.py \
    --cifs $MNTDIR \
    --accessions $CHUNK \
    --output chainsaw_results_${SLURM_ARRAY_TASK_ID}.csv \
    --scratch /tmp/chainsaw_scratch_${SLURM_ARRAY_TASK_ID}_$$

fusermount -u $MNTDIR
rmdir $MNTDIR
