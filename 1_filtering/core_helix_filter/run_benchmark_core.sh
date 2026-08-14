#!/usr/bin/env bash
#SBATCH --job-name=corehelix_bench
#SBATCH --account=nn1003k
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=0:20:00
#SBATCH --output=corehelix_bench_%j.log
# Score the copper-anchored core test (M7) on the 241-structure eye benchmark and
# compare to M1..M6. Structures are already extracted (no squashfs mount needed).
set -euo pipefail

SUBMITDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
B=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline
SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
FOLDSEEK=/cluster/projects/nn1003k/prog/foldseek/bin/foldseek
REFDIR=$B/1_filtering/foldseek
BENCH=$B/1_filtering/benchmark
PM=$BENCH/structures/B2ZB02_taxID_1404_model.cif      # B2ZB02 AF3 model == ref_PmTYR
cd "$SUBMITDIR"
cleanup(){ rm -rf bench_panel fs_bench_tmp; }
trap cleanup EXIT

# 1. the 251 benchmark accessions
awk -F'\t' 'NR>1{print $1}' "$BENCH/benchmark_results.tsv" > bench_accs.txt
echo "bench accs: $(wc -l < bench_accs.txt)"

# 2. Cu-bearing reference panel (same 5 refs as production)
rm -rf bench_panel; mkdir -p bench_panel
cp "$PM" bench_panel/ref_PmTYR.cif
for r in ref_2Y9W_Abisporus ref_5CE9_Jregia ref_1BT3_Ibatatas ref_1JS8_squid; do
    cp "$REFDIR/$r.pdb" bench_panel/
done

# 3. foldseek easy-search benchmark structures vs panel
rm -rf fs_bench_tmp; mkdir -p fs_bench_tmp
"$FOLDSEEK" easy-search "$BENCH/structures" bench_panel fs_bench_vs_panel.tsv fs_bench_tmp/tmp \
    --format-output "query,target,qstart,qend,qlen,tstart,tend,qtmscore,cigar" \
    -s 9.5 --threads 8 -e inf --max-seqs 50 2>&1 | tail -2
echo "foldseek hits: $(wc -l < fs_bench_vs_panel.tsv)"

# 4. copper-anchored core test on the benchmark
apptainer exec --bind "$B:$B" --bind "$SUBMITDIR:$SUBMITDIR" "$SIF" \
    python3 "$SUBMITDIR/core_helix_check.py" \
        --cif-dir "$BENCH/structures" --acc-list bench_accs.txt \
        --ref-dir "$REFDIR" --pmtyr "$PM" \
        --fs-tsv fs_bench_vs_panel.tsv \
        --output core_helix_bench.tsv --dmax 4.0 --pmin 70 --need 2

# 5. head-to-head vs the eye labels + M1..M6
apptainer exec --bind "$B:$B" --bind "$SUBMITDIR:$SUBMITDIR" "$SIF" \
    python3 "$SUBMITDIR/benchmark_eval.py" \
        --results "$BENCH/benchmark_results.tsv" \
        --m7 core_helix_bench.tsv \
        --out-prefix benchmark_core | tee benchmark_core_summary.txt
echo "Done: $(date)"
