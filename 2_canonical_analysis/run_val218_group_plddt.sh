#!/usr/bin/env bash
#SBATCH --job-name=val218_grp
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH --output=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/2_canonical_analysis/val218_grp_%j.log
# Per-structure Val218 Cα pLDDT for S/Y/K/E/H/R carriers (figure-rep selection).
set -euo pipefail
SUBMITDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
M=$(mktemp -d /tmp/val218grp_mnt_XXXX)
cleanup(){ fusermount -u "$M" 2>/dev/null || true; rmdir "$M" 2>/dev/null || true; }
trap cleanup EXIT
cd "$SUBMITDIR"
squashfuse "$SQSH" "$M"
apptainer exec --bind /cluster/work/projects/nn1003k:/cluster/work/projects/nn1003k --bind "$M:$M:ro" \
  "$SIF" python3 "$SUBMITDIR/val218_group_plddt.py" "$M"
