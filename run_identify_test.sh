#!/bin/bash
#SBATCH --job-name=ppo_id_test
#SBATCH --account=nn1003k
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:10:00
#SBATCH --output=ppo_id_test_%j.log

SUBMITDIR=$(pwd)
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MNTDIR=/tmp/sqsh_ppoid_$$

module load NRIS/CPU
module load Python/3.11.5-GCCcore-13.2.0
source /cluster/work/projects/nn1003k/eirin/bioinf/chainsaw_env/bin/activate

mkdir -p "$MNTDIR"
squashfuse "$SQSH" "$MNTDIR"

python3 "$SUBMITDIR/identify_ppo_domain.py" \
    --cifs "$MNTDIR" \
    --chainsaw "$SUBMITDIR/chainsaw_results_all.csv" \
    --reference "$SUBMITDIR/consensus_reference.json" \
    --output "$SUBMITDIR/ppo_domain_test.csv" \
    --workers 4 \
    --limit 200

fusermount -u "$MNTDIR" 2>/dev/null || true
rmdir "$MNTDIR" 2>/dev/null || true
