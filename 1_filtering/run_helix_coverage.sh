#!/usr/bin/env bash
# Unified structural metric extraction: CEAlign onto PmTYR, extract helix coverage,
# pLDDT, residue identity at His positions, Cu ref probe, Cu-Cu distance.
# Then apply presence-based helix filter (cov >= 0.40 on all four helices).
set -euo pipefail
SUBMITDIR="$(cd "$(dirname "$0")" && pwd)"
SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
PMTYR=/cluster/work/projects/nn1003k/eirin/fold_analysis/pmtyr.pdb
HELIX_SS="$SUBMITDIR/../3_noncanonical_analysis/helix_ss.py"

M=$(mktemp -d /tmp/helixcov_XXXX)
cleanup(){ fusermount -u "$M" 2>/dev/null || true; rmdir "$M" 2>/dev/null || true; }
trap cleanup EXIT
squashfuse "$SQSH" "$M"

# CEAlign + extract all metrics
apptainer exec --cleanenv --bind /cluster/work/projects/nn1003k:/cluster/work/projects/nn1003k \
  --bind "$M:$M:ro" --bind "$(dirname $PMTYR):$(dirname $PMTYR):ro" "$SIF" \
  python3 "$HELIX_SS" "$M" "$PMTYR" "$SUBMITDIR/all_accessions.txt" "$SUBMITDIR/helix_perhelix.tsv" "${SLURM_CPUS_PER_TASK:-1}"

# Apply presence rule (cov >= 0.40) and write helix_coverage.tsv
apptainer exec --cleanenv --bind /cluster/work/projects/nn1003k:/cluster/work/projects/nn1003k "$SIF" \
  python3 -c "
import csv
COV_THRESH = 0.40
rows=[r for r in csv.DictReader(open('$SUBMITDIR/helix_perhelix.tsv'),delimiter='\t') if r['status']=='ok']
for r in rows:
    r['passes_helix'] = str(all(float(r[f'helix{h}_cov']) >= COV_THRESH for h in range(1, 5)))
cols = list(rows[0].keys()) if rows else []
with open('$SUBMITDIR/helix_coverage.tsv','w',newline='') as o:
    w=csv.DictWriter(o,fieldnames=cols,delimiter='\t',extrasaction='ignore'); w.writeheader()
    for r in rows: w.writerow(r)
print('helix_coverage.tsv:',len(rows),'rows; passes_helix True=',sum(r['passes_helix']=='True' for r in rows))
"

# Apply three-tier filters
apptainer exec --cleanenv --bind /cluster/work/projects/nn1003k:/cluster/work/projects/nn1003k "$SIF" \
  python3 "$SUBMITDIR/apply_filters.py"

echo "Done: $SUBMITDIR/helix_coverage.tsv + $SUBMITDIR/filter_results.csv"
