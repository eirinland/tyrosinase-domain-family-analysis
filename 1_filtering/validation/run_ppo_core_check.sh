#!/usr/bin/env bash
#SBATCH --job-name=ppo_core
#SBATCH --account=nn1003k
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --output=ppo_core_%j.log
set -euo pipefail

SUBMITDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
B=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline
SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
FOLDSEEK=/cluster/projects/nn1003k/prog/foldseek/bin/foldseek
cd "$SUBMITDIR"

M=/tmp/${SLURM_JOB_ID:-$$}_mnt
mkdir -p "$M"
cleanup(){ fusermount -u "$M" 2>/dev/null||true; rmdir "$M" 2>/dev/null||true; rm -rf nc_cand_cifs refdir fs_tmp; }
trap cleanup EXIT
squashfuse "$SQSH" "$M"
echo "mounted: $(ls "$M"|wc -l) files"

# 1. NOT-canonical candidate accessions (canonical != True)
apptainer exec --bind "$B:$B" "$SIF" python3 - "$B" "$SUBMITDIR/not_canonical_accessions.csv" <<'PY'
import sys,csv
B,out=sys.argv[1],sys.argv[2]
acc=[r["accession"] for r in csv.DictReader(open(B+"/canonical_criteria_all_ca.csv"))
     if r["canonical"].strip().lower()!="true"]
open(out,"w").write("accession\n"+"\n".join(acc)+"\n")
print("not-canonical candidates:",len(acc),file=sys.stderr)
PY
echo "candidates: $(($(wc -l < not_canonical_accessions.csv)-1))"

# 2. Copy candidate CIFs from mount
rm -rf nc_cand_cifs; mkdir -p nc_cand_cifs
tail -n +2 not_canonical_accessions.csv | while read acc; do
    f=$(ls "$M/${acc}_taxID_"*_model.cif 2>/dev/null | head -1)
    [ -n "$f" ] && cp "$f" nc_cand_cifs/
done
echo "copied: $(ls nc_cand_cifs|wc -l) cifs"

# 3. Reference dir (9 PPO refs)
rm -rf refdir; mkdir -p refdir
cp "$B/1_filtering/foldseek/"ref_*.pdb refdir/

# 4. Foldseek global TMalign (--alignment-type 1) vs 9 refs
rm -rf fs_tmp; mkdir -p fs_tmp
"$FOLDSEEK" easy-search nc_cand_cifs refdir fs_global.tsv fs_tmp/tmp \
    --alignment-type 1 \
    --format-output "query,target,qstart,qend,qlen,tstart,tend,qtmscore,cigar" \
    -s 9.5 --threads 16 -e inf --max-seqs 20 2>&1 | tail -3
echo "foldseek hits: $(wc -l < fs_global.tsv)"

# 5. Query-structure 4-helix core check
apptainer exec --bind "$B:$B" --bind "$SUBMITDIR:$SUBMITDIR" "$SIF" \
    python3 "$B/1_filtering/validation/ppo_core_check.py" \
        --cif-dir nc_cand_cifs --ref-dir refdir --fs-tsv fs_global.tsv \
        --output "$SUBMITDIR/ppo_core_results.csv"

# 6. Three-pool assignment
apptainer exec --bind "$B:$B" --bind "$SUBMITDIR:$SUBMITDIR" "$SIF" \
    python3 - "$B" "$SUBMITDIR" <<'PY'
import sys,csv
B,S=sys.argv[1],sys.argv[2]
canon=set(); allacc=set()
for r in csv.DictReader(open(B+"/canonical_criteria_all_ca.csv")):
    allacc.add(r["accession"])
    if r["canonical"].strip().lower()=="true": canon.add(r["accession"])
core_ok={}
for r in csv.DictReader(open(S+"/ppo_core_results.csv")):
    core_ok[r["accession"]] = r["core_ok"].strip().lower()=="true"
noncanon=set(a for a,ok in core_ok.items() if ok and a not in canon)
discarded=allacc - canon - noncanon
with open(S+"/three_pool_assignment.csv","w",newline="\n") as f:
    w=csv.writer(f); w.writerow(["accession","pool"])
    for a in sorted(allacc):
        p="canonical" if a in canon else "non_canonical" if a in noncanon else "discarded"
        w.writerow([a,p])
print(f"TOTAL={len(allacc)}  canonical={len(canon)}  non_canonical={len(noncanon)}  discarded={len(discarded)}",file=sys.stderr)
# how many candidates had no foldseek hit at all
scored=set(core_ok); cand=allacc-canon
print(f"candidates={len(cand)}  scored_by_foldseek={len(scored&cand)}  no_hit={len(cand-scored)}",file=sys.stderr)
PY
echo "Done: $(date)"
