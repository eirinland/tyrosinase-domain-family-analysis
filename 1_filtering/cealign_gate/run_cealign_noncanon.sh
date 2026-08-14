#!/usr/bin/env bash
#SBATCH --job-name=cealign_nc
#SBATCH --account=nn1003k
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=4:00:00
#SBATCH --output=cealign_nc_%j.log
set -euo pipefail

SUBMITDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
B=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline
SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
cd "$SUBMITDIR"

M=/tmp/${SLURM_JOB_ID:-$$}_mnt
mkdir -p "$M"
cleanup(){ fusermount -u "$M" 2>/dev/null||true; rmdir "$M" 2>/dev/null||true; rm -rf nc_cifs refdir; }
trap cleanup EXIT
squashfuse "$SQSH" "$M"
echo "mounted: $(ls "$M"|wc -l) files"

# 1. copy non-canonical candidate CIFs (strip CR: accession list is CRLF)
rm -rf nc_cifs; mkdir -p nc_cifs
while read -r acc; do
    [ -z "$acc" ] && continue
    f=$(ls "$M/${acc}_taxID_"*_model.cif 2>/dev/null | head -1 || true)
    [ -n "$f" ] && cp "$f" nc_cifs/ || true
done < <(tail -n +2 "$B/1_filtering/final_pools/noncanonical_accessions.csv" | tr -d "\r")
echo "copied: $(ls nc_cifs|wc -l) cifs"

# 2. references (9 PPO refs)
rm -rf refdir; mkdir -p refdir
cp "$B/1_filtering/foldseek/"ref_*.pdb refdir/

# 3. PyMOL cealign vs 9 refs
apptainer exec --bind "$B:$B" --bind "$SUBMITDIR:$SUBMITDIR" "$SIF" \
    python3 "$SUBMITDIR/cealign_noncanon.py" \
        --cif-dir nc_cifs --ref-dir refdir \
        --output "$SUBMITDIR/cealign_noncanon_results.csv"
echo "results: $(($(wc -l < cealign_noncanon_results.csv)-1)) rows"
echo "Done: $(date)"
