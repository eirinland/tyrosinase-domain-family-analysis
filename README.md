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
  validation/                   query-structure core check, QC of the canonical pool
2_canonical_analysis/           active-site position vectors, novelty enumeration (Table S4)
  supplementary/                curated sub-tables
  visualisation/                heatmaps (Figure 5), co-occurrence (Figure S3), taxonomy
3_noncanonical_analysis/        helix-anchored alignment, site classification (Figures 6–8)
  allmetal3d/                   AF3 vs AllMetal3D vs PinMyMetal comparison (Figure 10)
4_genome_neighbourhood/         non-canonical GNA: H5Pro, H6Gln (Figure 7, Table S8)
4_genome_neighbourhood_canonical/  canonical GNA over the 67 novelty groups
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

Stage numbering follows pipeline order. `4_genome_neighbourhood` (non-canonical) and
`4_genome_neighbourhood_canonical` are two separate analyses over disjoint group sets;
the 17 numbered scripts in the former generated **both**.

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
2_canonical_analysis/run_novelty_pipeline.sh        # novelty enumeration
3_noncanonical_analysis/run_stage3_extract.sh       # align + classify non-canonical sites
5_foldseek/run_foldseek_pools.sh                    # cluster the pooled structures
```

## Data availability

- **AF3 models (32,069 CIFs)**, PyMOL sessions, and bulk intermediates: Zenodo, see `ZENODO_MANIFEST.md`.
- **Crystal structure** of A0A8H9LF69 (*Microbispora bryophytorum*): PDB **32AE**.
- **Reference structures**: 9 chain-A PDBs in `1_filtering/foldseek/ref_*.pdb`, spanning
  bacterial, fungal, plant, human, molluscan, archaeal and oomycete PPO diversity.

## Provenance notes

This repository is an honest record, including where it is incomplete.

Two generating scripts were lost when the cluster scratch area was purged and are **not**
recoverable, although their outputs are present and were used in the manuscript:

- `detect_thioether.py` — produced `5_foldseek/pools/thioether_check.tsv`, the
  `cys_in_shell` (Type 1 vs Type 2) node attribute in the Figure 3 network.
- `extract_position_vectors_super.py` — produced `2_canonical_analysis/position_vectors_super.csv`.
  The non-`_super` variant, `extract_position_vectors.py`, is present.

The AF3-to-AllMetal3D alignment and merge step behind Figure 10 was also lost; the
resulting comparison tables and all plotting scripts are present, and closely related
v2 implementations survive.

Several supplementary tables at the repository root (`table_s3_*.tsv`, `table_gna_*.tsv`,
`table_bgc_pfam_markers.tsv`, `table_s5_hmm_nc_agreement.tsv`, `novelty_enumeration_compact.tsv`)
were assembled manually from the named pipeline outputs rather than by a script.

Files superseded during development (`*.preM7.bak`, `*.prewidened.bak`, and similar) were
excluded so that no one recomputes retired numbers. They are in the Zenodo archive.

## Citation

Please cite the manuscript and, for the structures, the Zenodo DOI.

## Licence

See `LICENSE`.
