#!/usr/bin/env bash
#SBATCH --job-name=hmmalign_ppo
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=hmmalign_%j.log
# Regenerate all_hmmalign.afa in-folder: PF00264 HMM alignment of the query seqs.
# Enables novelty_pipeline.py stages G,I. Inputs (PF00264.hmm, query.fasta) are
# vendored alongside this script, so no reach into the old bioinf_redo pipeline.
set -euo pipefail
SUBMITDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
SIF=/cluster/projects/nn1003k/prog/af3/af3_cpu_amd64.sif
module load NRIS/CPU
apptainer exec --cleanenv --bind "$SUBMITDIR:$SUBMITDIR" "$SIF" \
    /hmmer/bin/hmmalign --outformat afa \
        "$SUBMITDIR/PF00264.hmm" "$SUBMITDIR/query.fasta" > "$SUBMITDIR/all_hmmalign.afa"
echo "wrote $SUBMITDIR/all_hmmalign.afa ($(grep -c '^>' "$SUBMITDIR/all_hmmalign.afa") seqs)"
