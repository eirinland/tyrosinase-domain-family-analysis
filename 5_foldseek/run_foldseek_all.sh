#!/usr/bin/env bash
#SBATCH --job-name=fs_all
#SBATCH --account=nn1003k
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=foldseek_all_%j.log

set -euo pipefail

SUBMITDIR=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/foldseek/all_structures
SCRIPTDIR=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/foldseek

FOLDSEEK=/cluster/projects/nn1003k/prog/foldseek/bin/foldseek
PYTHON_SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh

ACCESSION_CSV=$SUBMITDIR/accessions.csv

MOUNT=/tmp/${SLURM_JOB_ID}_sqsh_mount
STAGING=/tmp/${SLURM_JOB_ID}_trimmed_cifs
FSTMP05=/tmp/${SLURM_JOB_ID}_fstmp05
FSTMP07=/tmp/${SLURM_JOB_ID}_fstmp07
FSTMP08=/tmp/${SLURM_JOB_ID}_fstmp08
FSTMP09=/tmp/${SLURM_JOB_ID}_fstmp09

THREADS=16
CLUSTER_MODE=1
PLDDT_CUTOFF=70
WINDOW=7

mkdir -p "$MOUNT" "$STAGING" "$FSTMP05" "$FSTMP07" "$FSTMP08" "$FSTMP09"
mkdir -p "$SUBMITDIR/results" "$SUBMITDIR/tm05/results" "$SUBMITDIR/tm07/results" "$SUBMITDIR/tm09/results"

echo "============================================"
echo " FoldSeek Clustering — ALL STRUCTURES"
echo "============================================"
echo "Job ID        : $SLURM_JOB_ID"
echo "Started       : $(date)"
echo "Accession CSV : $ACCESSION_CSV"
echo "pLDDT cutoff  : ${PLDDT_CUTOFF} (window: ${WINDOW} residues)"
echo "Cluster mode  : $CLUSTER_MODE (connected component)"
echo "Threads       : $THREADS"
echo ""

cleanup() {
    echo ""
    echo "Cleaning up at $(date)..."
    if mountpoint -q "$MOUNT" 2>/dev/null; then
        fusermount -u "$MOUNT" && echo "Unmounted $MOUNT" || echo "Warning: unmount failed"
    fi
    rm -rf "$STAGING" "$FSTMP05" "$FSTMP07" "$FSTMP08" "$FSTMP09"
    rmdir "$MOUNT" 2>/dev/null || true
}
trap cleanup EXIT

echo "Step 1: Mounting SquashFS..."
squashfuse "$SQSH" "$MOUNT"
echo "Mounted at $MOUNT — $(ls $MOUNT | wc -l) CIF files available"
echo ""

echo "Step 2: Trimming low-pLDDT termini and staging to /tmp..."
apptainer exec \
    --bind "$MOUNT:$MOUNT:ro" \
    --bind "$STAGING:$STAGING" \
    --bind "$SCRIPTDIR:$SCRIPTDIR" \
    "$PYTHON_SIF" \
    python3 "$SCRIPTDIR/trim_and_stage.py" \
        --cif-dir    "$MOUNT" \
        --accessions "$ACCESSION_CSV" \
        --output     "$STAGING" \
        --cutoff     "$PLDDT_CUTOFF" \
        --window     "$WINDOW" \
        --workers    "$THREADS"
echo ""

N_STAGED=$(ls "$STAGING" | wc -l)
echo "Trimmed CIFs staged: $N_STAGED"
echo ""

echo "Step 3a: Clustering at TM 0.5..."
$FOLDSEEK easy-cluster \
    "$STAGING" \
    "$SUBMITDIR/tm05/results/cluster" \
    "$FSTMP05" \
    --tmscore-threshold 0.50 \
    --cluster-mode "$CLUSTER_MODE" \
    --threads "$THREADS" \
    -v 3
echo "TM 0.5 done at $(date)"
echo ""

echo "Step 3b: Clustering at TM 0.7..."
$FOLDSEEK easy-cluster \
    "$STAGING" \
    "$SUBMITDIR/tm07/results/cluster" \
    "$FSTMP07" \
    --tmscore-threshold 0.70 \
    --cluster-mode "$CLUSTER_MODE" \
    --threads "$THREADS" \
    -v 3
echo "TM 0.7 done at $(date)"
echo ""

echo "Step 3c: Clustering at TM 0.8 (primary)..."
$FOLDSEEK easy-cluster \
    "$STAGING" \
    "$SUBMITDIR/results/cluster" \
    "$FSTMP08" \
    --tmscore-threshold 0.80 \
    --cluster-mode "$CLUSTER_MODE" \
    --threads "$THREADS" \
    -v 3
echo "TM 0.8 done at $(date)"
echo ""

echo "Step 3d: Clustering at TM 0.9..."
$FOLDSEEK easy-cluster \
    "$STAGING" \
    "$SUBMITDIR/tm09/results/cluster" \
    "$FSTMP09" \
    --tmscore-threshold 0.90 \
    --cluster-mode "$CLUSTER_MODE" \
    --threads "$THREADS" \
    -v 3
echo "TM 0.9 done at $(date)"
echo ""

echo "Step 4: Thioether detection..."
apptainer exec \
    --bind "$MOUNT:$MOUNT:ro" \
    --bind "$SCRIPTDIR:$SCRIPTDIR" \
    --bind "$SUBMITDIR:$SUBMITDIR" \
    "$PYTHON_SIF" \
    python3 "$SCRIPTDIR/detect_thioether.py" \
        --cifs "$MOUNT" \
        --accessions "$ACCESSION_CSV" \
        --output "$SUBMITDIR/thioether_check.tsv" \
        --workers "$THREADS"
echo ""

echo "Step 5: Comparison summary..."
apptainer exec \
    --bind "$SUBMITDIR:$SUBMITDIR" \
    "$PYTHON_SIF" \
    python3 - <<'PYEOF'
import csv
from collections import defaultdict, Counter

BASE = "/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/foldseek/all_structures"

thresholds = {
    "TM 0.5": f"{BASE}/tm05/results/cluster_cluster.tsv",
    "TM 0.7": f"{BASE}/tm07/results/cluster_cluster.tsv",
    "TM 0.8": f"{BASE}/results/cluster_cluster.tsv",
    "TM 0.9": f"{BASE}/tm09/results/cluster_cluster.tsv",
}

status_map = {}
with open(f"{BASE}/accessions.csv") as f:
    for row in csv.DictReader(f):
        status_map[row['sequence_id']] = row['filter_status']

print(f"{'Threshold':<12} {'Clusters':>10} {'Singletons':>12} {'Largest':>10} {'Top 5 sizes':<40}")
print("-" * 90)

for label, path in thresholds.items():
    clusters = defaultdict(list)
    with open(path) as f:
        for line in f:
            rep, member = line.strip().split('\t')
            clusters[rep].append(member)
    sizes = sorted([len(v) for v in clusters.values()], reverse=True)
    n_sing = sum(1 for s in sizes if s == 1)
    top5 = ", ".join(str(s) for s in sizes[:5])
    print(f"{label:<12} {len(clusters):>10} {n_sing:>12} {sizes[0]:>10} {top5:<40}")

print()
print("=== TM 0.8 top 20 clusters with status breakdown ===")
clusters08 = defaultdict(list)
with open(f"{BASE}/results/cluster_cluster.tsv") as f:
    for line in f:
        rep, member = line.strip().split('\t')
        clusters08[rep].append(member)

sizes08 = sorted([(rep, len(m)) for rep, m in clusters08.items()], key=lambda x: -x[1])
for rep, size in sizes08[:20]:
    members = clusters08[rep]
    statuses = Counter()
    for m in members:
        base = m.replace('.cif', '')
        if base.endswith('_A'):
            base = base[:-2]
        base = base.split('_taxID_')[0]
        statuses[status_map.get(base, 'unknown')] += 1
    status_str = ', '.join(f'{s}:{n}' for s, n in statuses.most_common())
    print(f"  {size:>6d}  {status_str}")
PYEOF

echo ""
echo "============================================"
echo " Pipeline complete: $(date)"
echo "============================================"
