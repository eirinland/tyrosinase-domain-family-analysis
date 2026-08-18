#!/usr/bin/env python3
"""
Repoint hardcoded absolute paths at the repo's own location, so the folder is
portable (works wherever it is cloned) rather than tied to one machine.

Each edit is an exact string replacement, verified to apply exactly once.
parents[N] was checked individually per file against that file's real depth.

Deliberately NOT changed:
  - /cluster/... paths that refer to things which are not repo content:
    external tool installs (chainsaw), raw reference downloads (refs/),
    the AF3 CIF directory. These are a record of the HPC run.
  - SLURM .sh submit scripts: documentary, never intended to run locally.

Run:  python3 03_fix_paths.py [--dry-run]
"""
import sys
from pathlib import Path

REPO = Path("/Users/eirinlandsem/Library/Mobile Documents/com~apple~CloudDocs/"
            "Proteinkjemi_PhD/Skriving/Thesis/manuscripts:articles/Manuscript_TyrY/"
            "New_bioinf/bioinf_redo/ppo-family-structural-analysis")
DRY = "--dry-run" in sys.argv

MAC = ('"/Users/eirinlandsem/Library/Mobile Documents/com~apple~CloudDocs/Proteinkjemi_PhD/'
       'Skriving/Thesis/manuscripts:articles/Manuscript_TyrY/New_bioinf/bioinf_redo/'
       'Super_reference_pipeline"')
CLU = '"/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline"'

PLATFORM_BLOCK = f"""import platform
if platform.system() == "Darwin":
    BASE = {MAC}
else:
    BASE = {CLU}"""

PLATFORM_FIX = """from pathlib import Path as _P
BASE = str(_P(__file__).resolve().parents[2])   # repo root"""

PLATFORM_BLOCK2 = f"""import os, platform

if platform.system() == "Darwin":
    BASE = {MAC}
else:
    BASE = {CLU}"""

PLATFORM_FIX2 = """import os
from pathlib import Path as _P
BASE = str(_P(__file__).resolve().parents[2])   # repo root"""

# (relative path, old, new)  -- parents[N] verified per file
EDITS = [
    # --- A: platform-branch BASE -> resolve from the file's own location -------
    ("1_filtering/pool_summary/plot_pool_summary.py", PLATFORM_BLOCK, PLATFORM_FIX),
    ("1_filtering/pool_summary/plot_discarded_characterization.py", PLATFORM_BLOCK2, PLATFORM_FIX2),
    ("1_filtering/core_helix_filter/plot_threshold_grid.py", PLATFORM_BLOCK2, PLATFORM_FIX2),

    # --- B: PyMOL fallback candidate pointing at the sibling copy --------------
    ("2_canonical_analysis/visualisation/pymol_active_site.py",
     '''PMTYR_PATH = ("/Users/eirinlandsem/Library/Mobile Documents/com~apple~CloudDocs/"
              "Proteinkjemi_PhD/Skriving/Thesis/manuscripts:articles/Manuscript_TyrY/"
              "New_bioinf/bioinf_redo/Super_reference_pipeline/1_filtering/"
              "B2ZB02_taxID_1404_model.cif")''',
     '''import os
_HERE = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
PMTYR_PATH = os.path.join(_HERE, '..', '..', '1_filtering', 'B2ZB02_taxID_1404_model.cif')'''),

    ("2_canonical_analysis/visualisation/pymol_val218_active_site.py",
     """_CANDIDATES.append('/Users/eirinlandsem/Library/Mobile Documents/com~apple~CloudDocs/'
                   'Proteinkjemi_PhD/Skriving/Thesis/manuscripts:articles/Manuscript_TyrY/'
                   'New_bioinf/bioinf_redo/Super_reference_pipeline/2_canonical_analysis/'
                   'visualisation')""",
     """# (no absolute fallback: the folder is portable, so __file__/cwd are the anchors)"""),

    ("2_canonical_analysis/visualisation/pymol_variants_active_site.py",
     """_CANDIDATES.append('/Users/eirinlandsem/Library/Mobile Documents/com~apple~CloudDocs/'
                   'Proteinkjemi_PhD/Skriving/Thesis/manuscripts:articles/Manuscript_TyrY/'
                   'New_bioinf/bioinf_redo/Super_reference_pipeline/2_canonical_analysis/'
                   'visualisation')""",
     """# (no absolute fallback: the folder is portable, so __file__/cwd are the anchors)"""),

    # --- D: stage 4 referring to sibling stages via the old parent layout -------
    # These scripts used to sit one level ABOVE Super_reference_pipeline, so
    # BASE/'Super_reference_pipeline'/X was right. Inside the repo, BASE already
    # IS the root, so the extra segment must go.
    ("4_genome_neighbourhood/13_h6gln_neighbourhood.py",
     "WORK.parent / 'Super_reference_pipeline' / '3_noncanonical_analysis'",
     "WORK.parent / '3_noncanonical_analysis'"),
    ("4_genome_neighbourhood/14_h6gln_pfam.py",
     "WORK.parent / 'Super_reference_pipeline' / '3_noncanonical_analysis'",
     "WORK.parent / '3_noncanonical_analysis'"),
    ("4_genome_neighbourhood/15_all_groups_neighbourhood.py",
     "BASE / 'Super_reference_pipeline' / '2_canonical_analysis'",
     "BASE / '2_canonical_analysis'"),
    ("4_genome_neighbourhood/11_his5pro_cluster_membership.py",
     """STAGE3 = Path('/Users/eirinlandsem/Library/Mobile Documents/com~apple~CloudDocs/'
              'Proteinkjemi_PhD/Skriving/Thesis/manuscripts:articles/Manuscript_TyrY/'
              'New_bioinf/bioinf_redo/Super_reference_pipeline/3_noncanonical_analysis/'
              'helix_and_gap_filtered_structures.tsv')""",
     "STAGE3 = WORK.parent / '3_noncanonical_analysis' / 'helix_and_gap_filtered_structures.tsv'"),

    # --- C: cluster paths whose targets ARE repo content -----------------------
    # note foldseek/pools on the cluster is 5_foldseek/pools in the repo
    ("5_foldseek/pools/make_node_table.py",
     "BASE = Path('/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/foldseek/pools')",
     "BASE = Path(__file__).resolve().parent          # 5_foldseek/pools"),
    ("5_foldseek/pools/make_node_table.py",
     "TAXONOMY = Path('/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/taxonomy_lookup.csv')",
     "TAXONOMY = BASE.parents[1] / 'taxonomy_lookup.csv'"),
    ("5_foldseek/pools/make_node_table.py",
     "AGGER_XLSX = Path('/cluster/work/projects/nn1003k/eirin/bioinf/Agger_sequences_and_groups.xlsx')",
     "AGGER_XLSX = BASE.parents[1] / 'data' / 'Agger_sequences_and_groups.xlsx'"),

    ("5_foldseek/pools/cterm_domain_similarity/tm09_placement.py",
     'P = "/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/foldseek/pools"',
     'import os\nP = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # 5_foldseek/pools'),

    ("5_foldseek/pools/cterm_domain_similarity/characterized_overlap.py",
     'PIPE = "/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline"\n'
     'CTERM = f"{PIPE}/foldseek/pools/cterm_domain_similarity/c3_fungi"\n'
     'CLUST08 = f"{PIPE}/foldseek/pools/results/cluster_cluster.tsv"',
     'import os\n'
     'PIPE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repo root\n'
     'CTERM = f"{PIPE}/5_foldseek/pools/cterm_domain_similarity/c3_fungi"\n'
     'CLUST08 = f"{PIPE}/5_foldseek/pools/results/cluster_cluster.tsv"'),
    ("5_foldseek/pools/cterm_domain_similarity/characterized_overlap.py",
     'TAX = f"{PIPE}/taxonomy_lookup.csv"',
     'TAX = f"{PIPE}/taxonomy_lookup.csv"'),  # already correct once PIPE is the root
]

ok = fail = 0
for rel, old, new in EDITS:
    p = REPO / rel
    if not p.exists():
        print(f"MISSING FILE  {rel}"); fail += 1; continue
    t = p.read_text()
    n = t.count(old)
    if n != 1:
        if old == new and n >= 1:
            print(f"no-op ok      {rel}"); ok += 1; continue
        print(f"MATCH={n} (need 1)  {rel}\n    looking for: {old[:70]}..."); fail += 1; continue
    if not DRY:
        p.write_text(t.replace(old, new, 1))
    print(f"{'would patch' if DRY else 'patched'}   {rel}")
    ok += 1

print(f"\n{'-'*50}\napplied: {ok}   failed: {fail}")
sys.exit(1 if fail else 0)
