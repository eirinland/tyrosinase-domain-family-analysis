#!/usr/bin/env bash
#SBATCH --job-name=ppo_bench
#SBATCH --account=nn1003k
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --output=ppo_bench_%j.log
set -euo pipefail

SUBMITDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
B=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline
SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
FOLDSEEK=/cluster/projects/nn1003k/prog/foldseek/bin/foldseek
V="$B/1_filtering/validation"
cd "$SUBMITDIR"

M=/tmp/${SLURM_JOB_ID:-$$}_mnt; mkdir -p "$M"
cleanup(){ fusermount -u "$M" 2>/dev/null||true; rmdir "$M" 2>/dev/null||true; rm -rf bench_cifs refdir fs_tmp; }
trap cleanup EXIT
squashfuse "$SQSH" "$M"

# 1. Select the stratified benchmark set
apptainer exec --bind "$B:$B" "$SIF" python3 select_benchmark.py \
    --pool "$V/three_pool_assignment.csv" --core "$V/ppo_core_results.csv" \
    --disagree "$V/canonical_helix_disagreements.csv" \
    --bestbt "$B/1_filtering/foldseek/foldseek_multiref_best_bt.tsv" \
    --output benchmark_set.csv
echo "benchmark set: $(($(wc -l < benchmark_set.csv)-1)) structures"

# 2. Copy benchmark CIFs + the 9 reference PDBs
rm -rf bench_cifs refdir; mkdir -p bench_cifs refdir
tail -n +2 benchmark_set.csv | cut -d, -f1 | while read acc; do
    f=$(ls "$M/${acc}_taxID_"*_model.cif 2>/dev/null | head -1)
    [ -n "$f" ] && cp "$f" bench_cifs/
done
cp "$B/1_filtering/foldseek/"ref_*.pdb refdir/
cp "$B/1_filtering/foldseek/"ref_*.pdb bench_cifs/   # refs are queries too (positive controls)
echo "copied: $(ls bench_cifs|wc -l) queries (incl 9 refs) + $(ls refdir|wc -l) refs"

# 3. Foldseek global TMalign (all metrics: qtm, ttm, alntm, lddt, cigar)
rm -rf fs_tmp; mkdir -p fs_tmp
"$FOLDSEEK" easy-search bench_cifs refdir fs_global.tsv fs_tmp/tmp \
    --alignment-type 1 \
    --format-output "query,target,qstart,qend,qlen,tstart,tend,qtmscore,ttmscore,alntmscore,lddt,alnlen,cigar" \
    -s 9.5 --threads 16 -e inf --max-seqs 20 2>&1 | tail -2
echo "foldseek hits: $(wc -l < fs_global.tsv)"

# 4. Score all six methods (foldseek + pymol super/cealign + biotite intrinsic)
apptainer exec --bind "$B:$B" "$SIF" python3 run_methods.py \
    --cif-dir bench_cifs --ref-dir refdir --fs-tsv fs_global.tsv \
    --bench benchmark_set.csv --criteria "$B/canonical_criteria_all_ca.csv" \
    --pool "$V/three_pool_assignment.csv" --output benchmark_methods_raw.tsv

# 5. Summarize + disagreement queue
apptainer exec --bind "$B:$B" "$SIF" python3 summarize.py \
    benchmark_methods_raw.tsv benchmark_summary.txt benchmark_disagreements.tsv

echo "Done: $(date)"
