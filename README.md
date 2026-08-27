# ppo-family-structural-analysis

Structural bioinformatics pipeline for the family-scale analysis of the polyphenol
oxidase (PPO) fold, Pfam **PF00264**, across **32,069 AlphaFold 3 predicted structures**.

This repository contains the code and result tables behind the manuscript
*Family-scale structural analysis of the PPO fold*. It covers the full path from raw
AF3 models to the three structural pools, the canonical active-site vector analysis,
the non-canonical (pseudoenzyme) classification, structural clustering, genome
neighbourhood analysis, and Pfam domain annotation.

## Headline results

| Pool | Structures | Definition |
|---|---:|---|
| **Canonical** | **21,893** (68.3%) | 2 Cu, Cu–Cu 2.8–5.5 Å, 6 His NE2 ≤3.5 Å, per-His Cα pLDDT ≥70 |
| **Non-canonical** | **1,060** (3.3%) | <6 coordinating His **and** passes the M7 copper-anchored core-helix test |
| **Discarded** | **9,116** (28.4%) | six-His but failing canonical geometry (278), or <6 His failing M7 (8,838) |

Non-canonical site tiers (chemical plausibility, `classify_sites.py` v2):
**117 binuclear / 667 mononuclear / 276 no-Cu (divergent)**.

Structural clustering of the 22,953 pooled structures (22,952 clustered):
**1,075 clusters** at TM 0.8; the four largest account for 82.8% of the dataset.

The M7 core-helix test was validated on a 241-structure hand-labelled benchmark:
**TP 139, FP 0, FN 23, TN 79, accuracy 90.5%, precision 1.00, F1 0.92**. A
400-combination threshold sweep confirmed **FP = 0 across the entire grid**, so the
precision is architectural rather than threshold-tuned.

## Layout

```
1_filtering/                    pools: Foldseek multi-ref, canonical criteria, M7 core-helix test
  benchmark/                    seven-method comparison on the labelled benchmark (Table S1)
  core_helix_filter/            the M7 test, threshold sweep, retriage
  final_pools/                  THE authoritative pool assignment (32,069 rows)
  pool_summary/                 Figure S1 panels, discarded-pool characterisation
  seed_variability/             AF3 seed reproducibility
  separate_helix_check/         sensitivity test: must the four M7 anchors sit on four helices?
  validation/                   query-structure core check, QC of the canonical pool
2_canonical_analysis/           active-site position vectors, novelty enumeration (Table S4)
  supplementary/                curated sub-tables
  visualisation/                heatmaps (Figure 5), co-occurrence (Figure S3), taxonomy
3_noncanonical_analysis/        helix-anchored alignment, site classification (Figures 6–8)
  allmetal3d/                   AF3 vs AllMetal3D vs PinMyMetal comparison (Figure 10)
  hmm_vector_check/             per-structure HMM vs structure at the six His positions
4_genome_neighbourhood/         genome neighbourhood analysis (one stage, both runs)
  (top level)                   non-canonical run: H5Pro, H6Gln (Figure 7, Table S8) + the 17 scripts
  canonical/                    canonical run over the 67 novelty groups
  groups/                       per-group neighbourhood tables
5_foldseek/                     clustering + all-vs-all network (Figure 3)
6_domain_analysis/              Pfam/hmmscan domain annotation, Chainsaw segmentation
data/                           characterised PPOs, Meitil clade mapping, taxonomy lookup
figures/                        rendered manuscript figures (network, heatmap, pipeline, crystal)
supplementary_tables/           the curated supplementary tables (S3, S5, GNA highlights)
docs/                           provenance audit
```

The repository root additionally holds the earliest pipeline steps, which ran before the
stage folders existed: `check_canonical_criteria.py` and `run_canonical_check.sh` (the
canonical geometry test, producing `canonical_criteria_all_ca.csv`), `identify_ppo_domain.py`,
`run_chainsaw_*` (domain segmentation, producing `chainsaw_results_all.csv`), and the
AF3-versus-AFDB comparison scripts. `canonical_criteria_all_ca.csv` is the input to
`1_filtering/assign_pools.py` and therefore to the entire classification.

Stage numbering follows pipeline order. Stage 4 covers two analyses over disjoint group
sets, both produced by the same 17 numbered scripts: the **non-canonical** run (site and
function classes such as `mononuclear_ProCuB` = H5Pro, `hemocyanin`, `oAPO`) sits at the
folder's top level alongside the scripts, and the **canonical** run over the 67 novelty
substitution groups (`Gly46_N`, `Arg209_E`, …) is in `canonical/`.

## Reproducing

Scripts are self-contained and stdlib-first; see `requirements.txt`. Most were written
to run on an HPC cluster (SLURM), so `run_*.sh` wrappers carry the job configuration
and the `*.py` files hold the logic. Paths in the wrappers point at the original cluster
layout and will need adjusting.

The pipeline's input is the set of 32,069 AF3 model CIFs, named
`<accession>_taxID_<taxid>_model.cif`. These are **not** in this repository (2.55 GB);
see the Zenodo deposition listed in `ZENODO_MANIFEST.md`.

Rough order:

```
1_filtering/foldseek/run_foldseek_multiref_bt.sh    # align vs the 9 references
1_filtering/run_canonical_check.sh                  # alignment-free Cu/His/pLDDT criteria
1_filtering/core_helix_filter/run_core_helix_filter.sh   # M7 core-helix test
1_filtering/assign_pools.py                         # -> final_pools/
2_canonical_analysis/run_extract.sh                 # position vectors
2_canonical_analysis/hmm/build_alignment.py         # PF00264 alignment -> match-state table
2_canonical_analysis/run_novelty_pipeline.sh        # novelty enumeration
3_noncanonical_analysis/run_stage3_extract.sh       # align + classify non-canonical sites
3_noncanonical_analysis/hmm_vector_check/           # HMM vs structure at the six His positions
5_foldseek/run_foldseek_pools.sh                    # cluster the pooled structures
```

### The profile-HMM comparisons

Both HMM-versus-structure results (Table S4's `hmm_agreement.tsv` for the canonical
second-sphere positions, and the six-His comparison for the non-canonical pool) run from
files in this repository, with no cluster access and no HMMER installation:

```bash
python3 2_canonical_analysis/hmm/build_alignment.py          # ~10 s, needs pyhmmer
python3 2_canonical_analysis/novelty_pipeline.py \
        --hmm-cols 2_canonical_analysis/hmm/hmm_match_columns.tsv.gz   # stages G + I
python3 3_noncanonical_analysis/hmm_vector_check/hmm_vector_check.py   # pandas only
```

`build_alignment.py` aligns the deposited `hmm/query.fasta.gz` (32,753 sequences) to
`PF00264.hmm` with pyhmmer, which embeds the same HMMER3 code as `hmmalign`, and writes
`hmm_match_columns.tsv.gz` (one row per sequence, one character per match state, 2 MB)
plus `reference_map.tsv` (PmTYR residue number → match state). That table, not the
782 MB `all_hmmalign.afa`, is what the analyses consume; `--afa` still works if you would
rather generate the full alignment with `run_hmmalign.sh`. Regenerating and re-running
stage I reproduces the deposited `supplementary/hmm_agreement.tsv` byte for byte.

Note that stage E and stage H of `novelty_pipeline.py` need `--cifs` (the AF3 model
archive on Zenodo). Without it, stages A–D, F, G and I run, but `flagged_groups.tsv` will
be missing the six groups that stage E rescues on copper distance; `hmm_agreement.tsv`
is unaffected.

## Data availability

- **AF3 models (32,069 CIFs)**, PyMOL sessions, and bulk intermediates: Zenodo, see `ZENODO_MANIFEST.md`.
- **Crystal structure** of A0A8H9LF69 (*Microbispora bryophytorum*): PDB **32AE**.
- **Reference structures**: 9 chain-A PDBs in `1_filtering/foldseek/ref_*.pdb`, spanning
  bacterial, fungal, plant, human, molluscan, archaeal and oomycete PPO diversity.

## Verifying the manuscript numbers

`verify_manuscript_numbers.py` regenerates the manuscript values that no other
script in this repository produces, printing each published value beside the
recomputed one together with the convention used. It reads only tracked files.

```bash
python3 verify_manuscript_numbers.py                # local blocks only
python3 verify_manuscript_numbers.py --identity     # adds a network block
python3 verify_manuscript_numbers.py --tsv out.tsv  # machine-readable summary
```

All claims it checks agree with the manuscript. Add `--identity` for a seventh
that needs network access and `biopython`.

Paths are resolved by `repo_paths.py`, which finds the analysis root by walking
up from the script, so a fresh clone works with no configuration. Set `PPO_BASE`
to point at a working tree elsewhere.

## Provenance notes

This repository is an honest record, including where it is incomplete.

**Stale thresholds in `core_helix_check.py`, corrected.** The file carried the
pre-calibration helicity windows `HELIX_D3 (4.8, 6.4)` and `HELIX_D4 (5.4, 7.4)`
rather than the final widened windows `(4.0, 6.4)` and `(4.8, 8.2)` reported in
the manuscript Methods. Running it as deposited gave `core_ok` = 1,036 and an M7
accuracy of 74.7% instead of the published 1,060 and 90.5%. The windows are now
the widened values, which reproduce 1,060 and agree with 100% of the `core_ok`
calls in `1_filtering/final_pools/three_pool_assignment_final.csv`. See
`1_filtering/si_discarded_examples/widening_check.txt`. The strict windows in
`benchmark/run_methods.py` and `validation/ppo_core_check.py` are **correct** and
unchanged, since those implement methods M1 to M6, which were benchmarked under
the pre-widening definition.

**A thioether statistic was corrected during review.** An earlier draft reported
that 90% of canonical structures bearing the thioether carried a non-core
C-terminal domain (5,516/6,151). That numerator could not be obtained under any
tested definition and has been replaced by the reproducible measurement, 96%
carrying an additional non-PPO domain (5,916/6,151), with the companion figure
corrected from 3,769 to 3,766 so that both come from the same computation. The
description was corrected alongside the numbers, because the rule that
reproduces the second half of the sentence counts any non-PPO domain rather than
specifically a C-terminal one. Restricting to domains C-terminal to the PPO
range gives 4,950/6,151 (80.5%) and breaks the second half as well
(2,026/14,383 = 14.1% against 26%). Of the 5,916 structures that carry an extra
domain, it is C-terminal to the tyrosinase domain in 4,950 (84%), which is why
the Discussion retains that observation as a qualifier. Note that Chainsaw
merges the two tyrosinase lobes into a single PPO domain, so lobe segmentation
does not affect any of these counts. Run `verify_manuscript_numbers.py` to see
every tested definition.

Two generating scripts were lost when the cluster scratch area was purged and are **not**
recoverable, although their outputs are present and were used in the manuscript:

- `detect_thioether.py` — produced `5_foldseek/pools/thioether_check.tsv`, the
  `cys_in_shell` (Type 1 vs Type 2) node attribute in the Figure 3 network.
- `extract_position_vectors_super.py` — produced `2_canonical_analysis/position_vectors_super.csv`.
  The non-`_super` variant, `extract_position_vectors.py`, is present.

The AF3-to-AllMetal3D alignment and merge step behind Figure 10 was also lost; the
resulting comparison tables and all plotting scripts are present, and closely related
v2 implementations survive.

Several supplementary tables in `supplementary_tables/` (`table_s3_*.tsv`, `table_gna_*.tsv`,
`table_bgc_pfam_markers.tsv`, `novelty_enumeration_compact.tsv`) were assembled manually
from the named pipeline outputs rather than by a script.

`table_s5_hmm_nc_agreement.tsv` (the six-His HMM comparison for the non-canonical pool)
was also assembled by hand and its generator is lost. It has since been rebuilt from the
deposited inputs by `3_noncanonical_analysis/hmm_vector_check/`, which writes the
per-structure table `hmm_vs_structure_nc.tsv` that the aggregate table summarises. The
rebuild reproduces the published per-position counts exactly at five of the six positions;
at CuB_His1 it gives 233 agree / 44 disagree against the published 232 / 45. The joint
"complete six-position vectors correct" count, quoted in the manuscript as 504 of 1,060,
comes out as **502**: every one of the six positions recovered as the same residue at the
same residue number, with both a deletion and an unmapped structural position counted as
failures. Two looser figures are reported alongside it so the difference is auditable —
505 if only the residue type must match, and 512 if the two rendered vector strings are
compared, which wrongly scores an unmapped structural position against an HMM deletion as
a match in 8 position-instances across 7 structures. The remaining gap to 504 is two
structures, and the lost script's exact convention cannot be recovered; the per-structure
table lets any reader recompute the figure under a stated definition.

Files superseded during development (`*.preM7.bak`, `*.prewidened.bak`, and similar) were
excluded so that no one recomputes retired numbers. They are in the Zenodo archive.

## Citation

Please cite the manuscript and, for the structures, the Zenodo DOI.

## Licence

See `LICENSE`.
