#!/bin/bash
#SBATCH --job-name=af3_vs_afdb
#SBATCH --account=nn1003k
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=af3_vs_afdb_%j.log

SUBMITDIR=$(pwd)
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
MNTDIR=/tmp/sqsh_compare_$$
AFDB_DIR=$SUBMITDIR/afdb_cifs

export http_proxy=http://10.63.2.48:3128/
export https_proxy=http://10.63.2.48:3128/

module load NRIS/CPU
module load Python/3.11.5-GCCcore-13.2.0

mkdir -p "$MNTDIR" "$AFDB_DIR"
squashfuse "$SQSH" "$MNTDIR"

python3 "$SUBMITDIR/compare_af3_afdb.py" \
    --cifs "$MNTDIR" \
    --canonical-csv "$SUBMITDIR/canonical_criteria_all.csv" \
    --output "$SUBMITDIR/af3_vs_afdb_comparison.csv" \
    --afdb-dir "$AFDB_DIR" \
    --sample 200

fusermount -u "$MNTDIR" 2>/dev/null || true
rmdir "$MNTDIR" 2>/dev/null || true
