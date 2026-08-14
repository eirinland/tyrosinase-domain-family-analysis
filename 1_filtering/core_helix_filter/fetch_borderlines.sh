#!/usr/bin/env bash
# Re-run the selector on the authoritative Olivia assignment, then pull each
# labelled CIF from the squashfs mount into borderlines/ for PyMOL review.
set -euo pipefail
SUBMITDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
cd "$SUBMITDIR"

python3 select_borderlines.py --assign three_pool_assignment_corehelix.csv --out borderlines_manifest.tsv

MOUNT=$(mktemp -d)
cleanup(){ fusermount -u "$MOUNT" 2>/dev/null || true; rmdir "$MOUNT" 2>/dev/null || true; }
trap cleanup EXIT
squashfuse "$SQSH" "$MOUNT"

rm -rf borderlines; mkdir -p borderlines
miss=0
while IFS=$'\t' read -r acc label outfn rest; do
    src=$(ls "$MOUNT/${acc}_taxID_"*_model.cif 2>/dev/null | head -1 || true)
    if [ -n "$src" ]; then cp "$src" "borderlines/$outfn"; else echo "MISSING: $acc"; miss=$((miss+1)); fi
done < <(tail -n +2 borderlines_manifest.tsv)
cp borderlines_manifest.tsv borderlines/
echo "staged $(ls borderlines/*.cif 2>/dev/null | wc -l) cifs, $miss missing"
