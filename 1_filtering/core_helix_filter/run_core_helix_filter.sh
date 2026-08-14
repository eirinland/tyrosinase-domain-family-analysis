#!/usr/bin/env bash
#SBATCH --job-name=core_helix_filter
#SBATCH --account=nn1003k
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --output=core_helix_filter_%j.log
# Copper-anchored 4-helix core test over ALL <6-His failed-canonical candidates
# (9,898), then re-triage the pools. Run on the small (amd64) partition.
set -euo pipefail

SUBMITDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
B=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline
SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
FOLDSEEK=/cluster/projects/nn1003k/prog/foldseek/bin/foldseek
REFDIR=$B/1_filtering/foldseek
ASSIGN=$B/1_filtering/final_pools/three_pool_assignment_final.csv
DMAX=4.0; PMIN=70; NEED=1

MOUNT=/tmp/${SLURM_JOB_ID:-$$}_mnt
cd "$SUBMITDIR"
cleanup(){ fusermount -u "$MOUNT" 2>/dev/null||true; rmdir "$MOUNT" 2>/dev/null||true;
           rm -rf scope_cifs scope_panel fs_tmp mount_index.txt scope_files.txt; }
trap cleanup EXIT
mkdir -p "$MOUNT"
squashfuse "$SQSH" "$MOUNT"
echo "mounted: $(ls "$MOUNT"|wc -l) files"
PM=$(ls "$MOUNT"/B2ZB02_taxID_*_model.cif)

# 1. scope = failed-canonical AND <6 coordinating His (the re-triable candidates)
awk -F, 'NR>1 && $3=="False" && ($4+0)<6 {print $1}' "$ASSIGN" > scope_accessions.txt
echo "scope accessions: $(wc -l < scope_accessions.txt)   (expect 9898)"

# 2. copy those CIFs from the mount (single-pass index join, no per-acc ls)
ls "$MOUNT" | grep '_model\.cif$' > mount_index.txt
awk -F'_taxID_' 'NR==FNR{w[$1]=1;next} ($1 in w){print}' scope_accessions.txt mount_index.txt > scope_files.txt
rm -rf scope_cifs; mkdir -p scope_cifs
while read -r f; do cp "$MOUNT/$f" scope_cifs/; done < scope_files.txt
echo "copied: $(ls scope_cifs|wc -l) cifs"

# 3. Cu-bearing reference panel: B2ZB02 AF3 model as ref_PmTYR + 4 crystal refs
#    (the 5 refs with a di-Cu site and hard-coded core-helix ranges in the script)
rm -rf scope_panel; mkdir -p scope_panel
cp "$PM" scope_panel/ref_PmTYR.cif
for r in ref_2Y9W_Abisporus ref_5CE9_Jregia ref_1BT3_Ibatatas ref_1JS8_squid; do
    cp "$REFDIR/$r.pdb" scope_panel/
done
echo "panel: $(ls scope_panel)"

# 4. foldseek easy-search scope vs panel (local align, all hits, with cigar)
rm -rf fs_tmp; mkdir -p fs_tmp
"$FOLDSEEK" easy-search scope_cifs scope_panel fs_scope_vs_panel.tsv fs_tmp/tmp \
    --format-output "query,target,qstart,qend,qlen,tstart,tend,qtmscore,cigar" \
    -s 9.5 --threads 8 -e inf --max-seqs 50 2>&1 | tail -2
echo "foldseek hits: $(wc -l < fs_scope_vs_panel.tsv)"

# 5. copper-anchored helicity core test
apptainer exec --bind "$MOUNT:$MOUNT:ro" --bind "$B:$B" --bind "$SUBMITDIR:$SUBMITDIR" "$SIF" \
    python3 "$SUBMITDIR/core_helix_check.py" \
        --cif-dir scope_cifs --acc-list scope_accessions.txt \
        --ref-dir "$REFDIR" --pmtyr "$PM" \
        --fs-tsv fs_scope_vs_panel.tsv \
        --output core_helix_results.tsv \
        --dmax "$DMAX" --pmin "$PMIN" --need "$NEED"

# 6. re-triage the three pools with the new core decision
apptainer exec --bind "$B:$B" --bind "$SUBMITDIR:$SUBMITDIR" "$SIF" \
    python3 "$SUBMITDIR/retriage.py" \
        --assignment "$ASSIGN" \
        --core core_helix_results.tsv \
        --outdir "$SUBMITDIR"
echo "Done: $(date)"
