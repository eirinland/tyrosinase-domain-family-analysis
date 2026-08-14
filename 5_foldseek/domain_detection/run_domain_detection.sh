#!/bin/bash
#SBATCH --job-name=domain_detect
#SBATCH --account=nn1003k
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=domain_detection_%j.log

set -eu

SUBMITDIR=$(pwd)
BASE=/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline
SIF=/cluster/projects/nn1003k/prog/olivia_scripts/python_tools/python_tools.sif
SQSH=/cluster/work/projects/nn1003k/ark/predicion_projects/eirin/tyr/af3_run/results/all_models.sqsh
CLUSTER_TSV=$BASE/5_foldseek/pools/results/cluster_cluster.tsv
WORKDIR=$BASE/5_foldseek/domain_detection
SCRIPT=$WORKDIR/domain_detection.py

# ── Cluster reps (v3 TM0.8 pool clusters) ───────────────────────────────
NTD_REP="A0AAV4B738_taxID_259542_model"   # 249 members (old cluster 8, mollusc N-terminal cap)
CTD_REP="A0ACC0WCW2_taxID_230839_model"   # 581 members (old cluster 5, oomycete C-terminal block)

# ── NTD: A0AA88XQR2, capping domain residues 27-63 (37 CAs) ────────────
NTD_REF_ACC="A0AA88XQR2_taxID_66713_model"

# ── CTD: A0A485L799, five helices, anchor = h3+h4 (blocking motif) ──────
# Helices: h1=349-368, h2=371-388, h3=400-406, h4=413-422, h5=433-441
# Anchor: h3+h4 (indices 2,3) = 17 CAs; total = 64 CAs
# Blocking residue: TRP406
CTD_REF_ACC="A0A485L799_taxID_120398_model"

module load NRIS/CPU

# Mount squashfs
MNTDIR=$(mktemp -d)
trap "fusermount -u $MNTDIR 2>/dev/null; rmdir $MNTDIR 2>/dev/null" EXIT
squashfuse $SQSH $MNTDIR

# ── Extract cluster members ──────────────────────────────────────────────
echo "Extracting cluster member lists..."
NTD_MEMBERS=$WORKDIR/ntd_members.txt
CTD_MEMBERS=$WORKDIR/ctd_members.txt
grep "^${NTD_REP}" "$CLUSTER_TSV" | cut -f2 > "$NTD_MEMBERS"
grep "^${CTD_REP}" "$CLUSTER_TSV" | cut -f2 > "$CTD_MEMBERS"
echo "  NTD cluster: $(wc -l < "$NTD_MEMBERS") members"
echo "  CTD cluster: $(wc -l < "$CTD_MEMBERS") members"

# ── Stage CIFs ───────────────────────────────────────────────────────────
stage_cifs() {
    local members_file=$1
    local out_dir=$2
    mkdir -p "$out_dir"
    local n=0
    while IFS= read -r member; do
        local name="$member"
        case "$name" in *_A|*_B) name="${name%_?}" ;; esac
        local src="$MNTDIR/${name}.cif"
        if [ -f "$src" ]; then
            cp "$src" "$out_dir/${member}.cif"
            n=$((n + 1))
        fi
    done < "$members_file"
    echo "  Staged $n CIFs to $out_dir"
}

NTD_CIFS=$WORKDIR/ntd_cifs
CTD_CIFS=$WORKDIR/ctd_cifs
stage_cifs "$NTD_MEMBERS" "$NTD_CIFS"
stage_cifs "$CTD_MEMBERS" "$CTD_CIFS"

# Stage reference CIFs
for ref_acc in "$NTD_REF_ACC" "$CTD_REF_ACC"; do
    ref_name="$ref_acc"
    case "$ref_name" in *_A|*_B) ref_name="${ref_name%_?}" ;; esac
    cp "$MNTDIR/${ref_name}.cif" "$WORKDIR/${ref_acc}.cif"
done

echo ""

# ── NTD detection (N-terminal capping domain, 37 CAs) ───────────────────
echo "=== NTD (N-terminal capping domain, 249 targets) ==="
echo "  Reference: A0AA88XQR2, residues 27-63 (37 CAs)"
apptainer exec --cleanenv \
    --bind "$WORKDIR:$WORKDIR" \
    "$SIF" python3 "$SCRIPT" ntd \
        --ref-cif "$WORKDIR/${NTD_REF_ACC}.cif" \
        --domain-start 27 \
        --domain-end 63 \
        --targets "$NTD_CIFS" \
        --output "$WORKDIR/ntd_results.tsv"

echo ""

# ── CTD detection (C-terminal blocking domain, 5 helices / 64 CAs) ──────
echo "=== CTD (C-terminal blocking domain, 581 targets) ==="
echo "  Reference: A0A485L799, anchor h3(400-406)+h4(413-422), gatekeeper TRP406"
apptainer exec --cleanenv \
    --bind "$WORKDIR:$WORKDIR" \
    "$SIF" python3 "$SCRIPT" ctd \
        --ref-cif "$WORKDIR/${CTD_REF_ACC}.cif" \
        --helices '[{"start":349,"end":368},{"start":371,"end":388},{"start":400,"end":406},{"start":413,"end":422},{"start":433,"end":441}]' \
        --anchor-idx "2,3" \
        --ref-resid 406 \
        --search-radius 15 \
        --targets "$CTD_CIFS" \
        --output "$WORKDIR/ctd_results.tsv"

echo ""

# ── Cleanup staged CIFs ─────────────────────────────────────────────────
rm -rf "$NTD_CIFS" "$CTD_CIFS"

echo "Done. Results:"
echo "  $WORKDIR/ntd_results.tsv"
echo "  $WORKDIR/ctd_results.tsv"
