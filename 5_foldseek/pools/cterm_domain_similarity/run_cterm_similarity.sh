#!/usr/bin/env bash
#SBATCH --job-name=cterm_sim
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/foldseek/pools/cterm_domain_similarity/cterm_sim_%j.log

set -euo pipefail

PIPE=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline
SUBMITDIR=$PIPE/foldseek/pools/cterm_domain_similarity

FOLDSEEK=/cluster/projects/nn1003k/prog/foldseek/bin/foldseek
PYTHON_SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
CLUSTER_TSV=$PIPE/foldseek/pools/results/cluster_cluster.tsv      # TM0.8 membership
PPO_CSV=$PIPE/1_filtering/foldseek/ppo_domain_assignment_multiref.csv

MOUNT=/tmp/${SLURM_JOB_ID}_sqsh_mount
THREADS=16
MIN_LEN=40

mkdir -p "$MOUNT"

echo "============================================"
echo " Within-cluster C-terminal domain similarity"
echo "============================================"
echo "Job ID  : $SLURM_JOB_ID"
echo "Started : $(date)"
echo ""

cleanup() {
    echo ""
    echo "Cleaning up at $(date)..."
    if mountpoint -q "$MOUNT" 2>/dev/null; then
        fusermount -u "$MOUNT" || true
    fi
    for d in /tmp/${SLURM_JOB_ID}_*; do
        [ -d "$d" ] && rm -rf "$d"
    done
}
trap cleanup EXIT

echo "Step 1: Mounting SquashFS..."
squashfuse "$SQSH" "$MOUNT"
echo "Mounted — $(ls "$MOUNT" | wc -l) CIF files available"
echo ""

# clusters by TM0.8 rank: c3 = Fungi (A0AAJ0MJ41), c4 = plant (A0AA88AFA7)
for SPEC in "c3_fungi:A0AAJ0MJ41" "c4_plant:A0AA88AFA7"; do
    LABEL=${SPEC%%:*}
    REP=${SPEC##*:}
    OUTDIR=$SUBMITDIR/$LABEL
    STAGING=/tmp/${SLURM_JOB_ID}_staging_${LABEL}
    mkdir -p "$OUTDIR" "$STAGING" "$OUTDIR/sub05" "$OUTDIR/sub07"

    echo "============================================"
    echo " $LABEL  (rep $REP)"
    echo "============================================"

    echo "  Extracting C-terminal domains..."
    apptainer exec \
        --bind "$MOUNT:$MOUNT:ro" \
        --bind "$STAGING:$STAGING" \
        --bind "$PIPE:$PIPE" \
        "$PYTHON_SIF" \
        python3 "$SUBMITDIR/extract_cterm.py" \
            --cif-dir     "$MOUNT" \
            --cluster-tsv "$CLUSTER_TSV" \
            --ppo-csv     "$PPO_CSV" \
            --rep         "$REP" \
            --output      "$STAGING" \
            --manifest    "$OUTDIR/manifest.csv" \
            --min-len     "$MIN_LEN" \
            --workers     "$THREADS"
    N_STAGED=$(ls "$STAGING" | wc -l)
    echo "  Staged C-term CIFs: $N_STAGED"
    echo ""

    echo "  Foldseek easy-search all-vs-all..."
    FSTMP=/tmp/${SLURM_JOB_ID}_fstmp_${LABEL}; mkdir -p "$FSTMP"
    $FOLDSEEK easy-search \
        "$STAGING" "$STAGING" \
        "$OUTDIR/allvsall.tsv" "$FSTMP" \
        --format-output "query,target,qtmscore,ttmscore,alnlen,qlen,tlen" \
        -s 4.0 --threads "$THREADS" -e inf --max-seqs 5000 -v 3
    echo "  Raw hits: $(wc -l < "$OUTDIR/allvsall.tsv")"
    echo ""

    echo "  Foldseek easy-cluster TM0.5 / TM0.7 (sub-populations)..."
    FSC5=/tmp/${SLURM_JOB_ID}_fsc5_${LABEL}; mkdir -p "$FSC5"
    $FOLDSEEK easy-cluster "$STAGING" "$OUTDIR/sub05/clu" "$FSC5" \
        --tmscore-threshold 0.50 --cluster-mode 1 --threads "$THREADS" -v 3
    FSC7=/tmp/${SLURM_JOB_ID}_fsc7_${LABEL}; mkdir -p "$FSC7"
    $FOLDSEEK easy-cluster "$STAGING" "$OUTDIR/sub07/clu" "$FSC7" \
        --tmscore-threshold 0.70 --cluster-mode 1 --threads "$THREADS" -v 3
    echo ""

    echo "  Summarizing..."
    apptainer exec --bind "$PIPE:$PIPE" "$PYTHON_SIF" \
        python3 "$SUBMITDIR/summarize_cterm.py" \
            --label    "$LABEL" \
            --allvsall "$OUTDIR/allvsall.tsv" \
            --manifest "$OUTDIR/manifest.csv" \
            --clust05  "$OUTDIR/sub05/clu_cluster.tsv" \
            --clust07  "$OUTDIR/sub07/clu_cluster.tsv" \
            --out      "$OUTDIR/summary.txt"
    echo ""
done

echo "Done at $(date)"
