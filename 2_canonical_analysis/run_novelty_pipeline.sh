#!/usr/bin/env bash
#SBATCH --job-name=novelty
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/2_canonical_analysis/novelty_%j.log
# One entry point: mounts the structures, runs the full consolidated pipeline
# (vector stages A-D,F + copper-distance stage E + HMM stage G).
set -euo pipefail
SUBMITDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
AFA="$SUBMITDIR/hmm/all_hmmalign.afa"   # produced by hmm/run_hmmalign.sh; enables stages G,I
M=$(mktemp -d /tmp/novelty_mnt_XXXX)
cleanup(){ fusermount -u "$M" 2>/dev/null || true; rmdir "$M" 2>/dev/null || true; }
trap cleanup EXIT
cd "$SUBMITDIR"
squashfuse "$SQSH" "$M"
if [ -f "$AFA" ]; then AFA_ARG="--afa $AFA"
else echo "NOTE: $AFA missing -> stages G,I skipped (run hmm/run_hmmalign.sh to enable)"; AFA_ARG=""; fi
apptainer exec --bind /cluster/work/projects/nn1003k:/cluster/work/projects/nn1003k --bind "$M:$M:ro" \
  "$SIF" python3 "$SUBMITDIR/novelty_pipeline.py" --cifs "$M" $AFA_ARG
