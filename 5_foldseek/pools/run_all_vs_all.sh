#!/usr/bin/env bash
#SBATCH --job-name=avsa_pools
#SBATCH --account=nn1003k
#SBATCH --partition=small
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/foldseek/pools/allvsall_%j.log

set -euo pipefail

SUBMITDIR=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/foldseek/pools

FOLDSEEK=/cluster/projects/nn1003k/prog/foldseek/bin/foldseek
PYTHON_SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh

MOUNT=/tmp/${SLURM_JOB_ID}_sqsh_mount
THREADS=16
PLDDT_CUTOFF=70
WINDOW=7

mkdir -p "$MOUNT"

echo "============================================"
echo " All-vs-All TM-score: pools cluster reps"
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
echo "Mounted — $(ls $MOUNT | wc -l) CIF files available"
echo ""

for TM_LABEL in tm05 tm07 tm08 tm09; do
    OUTDIR=$SUBMITDIR/all_vs_all_${TM_LABEL}
    STAGING=/tmp/${SLURM_JOB_ID}_staging_${TM_LABEL}
    FSTMP=/tmp/${SLURM_JOB_ID}_fstmp_${TM_LABEL}
    mkdir -p "$STAGING" "$FSTMP" "$OUTDIR"
    if [ "$TM_LABEL" = "tm08" ]; then
        CLUSTER_TSV=$SUBMITDIR/results/cluster_cluster.tsv
    else
        CLUSTER_TSV=$SUBMITDIR/$TM_LABEL/results/cluster_cluster.tsv
    fi
    echo "  Building rep_accessions.csv from $CLUSTER_TSV ..."
    { echo "sequence_id,filter_status"; cut -f1 "$CLUSTER_TSV" | sort -u | sed "s/\$/,unknown/"; } > "$OUTDIR/rep_accessions.csv"

    echo "============================================"
    echo " Processing $TM_LABEL"
    echo "============================================"
    N_REPS=$(tail -n +2 "$OUTDIR/rep_accessions.csv" | wc -l)
    echo "Reps: $N_REPS"
    echo ""

    echo "  Creating lookup CSV..."
    apptainer exec --bind "$OUTDIR:$OUTDIR" "$PYTHON_SIF" python3 -c "
import csv
with open('$OUTDIR/rep_accessions.csv') as fin, \
     open('$OUTDIR/rep_accessions_lookup.csv', 'w', newline='') as fout:
    reader = csv.DictReader(fin)
    writer = csv.DictWriter(fout, fieldnames=['sequence_id', 'filter_status'])
    writer.writeheader()
    for row in reader:
        sid = row['sequence_id']
        lookup = sid[:-2] if sid.endswith('_A') else sid
        writer.writerow({'sequence_id': lookup, 'filter_status': row['filter_status']})
"

    echo "  Trimming low-pLDDT termini..."
    apptainer exec \
        --bind "$MOUNT:$MOUNT:ro" \
        --bind "$STAGING:$STAGING" \
        --bind "$SUBMITDIR:$SUBMITDIR" \
        --bind "$OUTDIR:$OUTDIR" \
        "$PYTHON_SIF" \
        python3 "$SUBMITDIR/trim_and_stage.py" \
            --cif-dir    "$MOUNT" \
            --accessions "$OUTDIR/rep_accessions_lookup.csv" \
            --output     "$STAGING" \
            --cutoff     "$PLDDT_CUTOFF" \
            --window     "$WINDOW" \
            --workers    "$THREADS"

    echo "  Restoring _A suffix where needed..."
    apptainer exec --bind "$STAGING:$STAGING" --bind "$OUTDIR:$OUTDIR" "$PYTHON_SIF" python3 -c "
import csv, os
with open('$OUTDIR/rep_accessions.csv') as f:
    for row in csv.DictReader(f):
        orig = row['sequence_id']
        if not orig.endswith('_A'): continue
        src = os.path.join('$STAGING', orig[:-2] + '.cif')
        dst = os.path.join('$STAGING', orig + '.cif')
        if os.path.exists(src): os.rename(src, dst)
"
    N_STAGED=$(ls "$STAGING" | wc -l)
    echo "  Staged: $N_STAGED CIFs"
    echo ""

    echo "  Running foldseek easy-search all-vs-all..."
    $FOLDSEEK easy-search \
        "$STAGING" "$STAGING" \
        "$OUTDIR/allvsall_raw.tsv" "$FSTMP" \
        --format-output "query,target,qtmscore,ttmscore,alnlen,qlen,tlen" \
        -s 4.0 --threads "$THREADS" -e inf --max-seqs 5000 -v 3
    echo "  Search done at $(date)"
    echo "  Raw hits: $(wc -l < "$OUTDIR/allvsall_raw.tsv")"
    echo ""

    echo "  Building edge table..."
    apptainer exec --bind "$OUTDIR:$OUTDIR" "$PYTHON_SIF" python3 -c "
import csv
edges = {}
with open('$OUTDIR/allvsall_raw.tsv') as f:
    for line in f:
        parts = line.strip().split('\t')
        q, t = parts[0], parts[1]
        if q == t: continue
        min_tm = min(float(parts[2]), float(parts[3]))
        # foldseek appends a chain tag (_A) to staged ..._model_A.cif files,
        # producing a doubled ..._model_A_A; collapse it so edge names match
        # the cluster-rep node names (..._model_A).
        def _nm(x):
            x = x.replace('.cif','')
            return x[:-2] if x.endswith('_A_A') else x
        q_name, t_name = _nm(q), _nm(t)
        pair = tuple(sorted([q_name, t_name]))
        if pair not in edges or min_tm > edges[pair]:
            edges[pair] = min_tm
with open('$OUTDIR/edge_table.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['rep1','rep2','min_tm'])
    for (r1,r2), tm in sorted(edges.items(), key=lambda x: -x[1]):
        w.writerow([r1, r2, '{:.4f}'.format(tm)])
print('Total unique pairs: {}'.format(len(edges)))

with open('$OUTDIR/edge_table_tm04.tsv', 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['rep1','rep2','min_tm'])
    n = 0
    for (r1,r2), tm in sorted(edges.items(), key=lambda x: -x[1]):
        if tm >= 0.4:
            w.writerow([r1, r2, '{:.4f}'.format(tm)])
            n += 1
print('Pairs at TM >= 0.4: {}'.format(n))
"
    echo ""

    echo "  Building node tables..."
    apptainer exec \
        --bind /cluster/work/projects/nn1003k:/cluster/work/projects/nn1003k \
        --bind /cluster/projects/nn1003k:/cluster/projects/nn1003k:ro \
        "$PYTHON_SIF" \
        python3 "$SUBMITDIR/make_node_table.py" "$TM_LABEL"
    echo ""

    rm -rf "$STAGING" "$FSTMP"
    echo "  $TM_LABEL complete at $(date)"
    echo ""
done

echo "============================================"
echo " All done: $(date)"
echo "============================================"
