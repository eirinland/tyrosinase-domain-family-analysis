#!/usr/bin/env bash
#SBATCH --job-name=nc3_extract
#SBATCH --account=nn1003k
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --output=stage3_extract_%j.log
set -euo pipefail

SUBMITDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
B=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline
SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
FOLDSEEK=/cluster/projects/nn1003k/prog/foldseek/bin/foldseek
PROXY=http://10.63.2.48:3128/

MOUNT=/tmp/${SLURM_JOB_ID:-$$}_mnt
TAX=/tmp/${SLURM_JOB_ID:-$$}_tax
cd "$SUBMITDIR"
cleanup(){ fusermount -u "$MOUNT" 2>/dev/null||true; rmdir "$MOUNT" 2>/dev/null||true; rm -rf "$TAX" nc_cifs fs_tmp; }
trap cleanup EXIT
mkdir -p "$MOUNT" "$TAX"

echo "=== taxdump ==="
http_proxy="$PROXY" https_proxy="$PROXY" wget -q -O "$TAX/taxdump.tar.gz" \
    https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -xzf "$TAX/taxdump.tar.gz" -C "$TAX" nodes.dmp names.dmp merged.dmp

squashfuse "$SQSH" "$MOUNT"
echo "mounted: $(ls "$MOUNT"|wc -l) files"
PM=$(ls "$MOUNT"/B2ZB02_taxID_*_model.cif)

# 1. Build non-canonical pool CSV from final_pools (2026-06-14 rule), enrich from criteria
apptainer exec --bind "$B:$B" "$SIF" python3 - "$B" "$SUBMITDIR/nc_pool.csv" <<'PY'
import sys,csv
B,out=sys.argv[1],sys.argv[2]
pool=set(r['accession'] for r in csv.DictReader(open(B+'/1_filtering/final_pools/noncanonical_accessions.csv')))
w=csv.writer(open(out,'w')); w.writerow(['accession','failed_step','n_his','cu_dist'])
n=0
for r in csv.DictReader(open(B+'/canonical_criteria_all_ca.csv')):
    a=r['accession']
    if a not in pool: continue
    nh=int(r['n_coord_his']) if r['n_coord_his'].strip() else 0
    his_ok=r['his_ok'].strip().lower()=='true'
    step='step5_his_coordination' if not his_ok else 'step4_cu_distance'
    w.writerow([a,step,nh,r.get('cu_dist','')]); n+=1
print(f'pool: {n}',file=sys.stderr)
PY
echo "pool: $(($(wc -l < nc_pool.csv)-1)) structures"

# 2. Copy non-canonical CIFs from mount
rm -rf nc_cifs; mkdir -p nc_cifs
tail -n +2 nc_pool.csv | cut -d, -f1 | while read acc; do
    f=$(ls "$MOUNT/${acc}_taxID_"*_model.cif 2>/dev/null | head -1)
    [ -n "$f" ] && cp "$f" nc_cifs/
done
echo "copied: $(ls nc_cifs|wc -l) cifs"

# 3. Foldseek easy-search non-canonical vs B2ZB02 (fresh, AF3-model numbering)
rm -rf fs_tmp; mkdir -p fs_tmp
cp "$PM" fs_tmp/B2ZB02_ref.cif
"$FOLDSEEK" easy-search nc_cifs fs_tmp/B2ZB02_ref.cif fs_vs_b2zb02.tsv fs_tmp/tmp \
    --format-output "query,target,qstart,qend,qlen,tstart,tend,qtmscore,cigar" \
    -s 9.5 --threads 8 -e inf --max-seqs 5 2>&1 | tail -2
echo "foldseek hits: $(wc -l < fs_vs_b2zb02.tsv)"

# 4. Extract + align
apptainer exec --bind "$MOUNT:$MOUNT:ro" --bind "$B:$B" --bind "$TAX:$TAX:ro" \
    --bind "$SUBMITDIR:$SUBMITDIR" "$SIF" \
    python3 "$SUBMITDIR/stage3_extract_align.py" \
        --cif-dir nc_cifs --pool-csv nc_pool.csv --pmtyr "$PM" \
        --fs-tsv fs_vs_b2zb02.tsv --taxdump-dir "$TAX" \
        --output "$SUBMITDIR/noncanonical_analysis.tsv"
echo "Done: $(date)"
