# Manuscript → script provenance audit (independent)

**Goal:** confirm that a fresh repo built from `Super_reference_pipeline/` (SRP) can regenerate every number, figure, and table in the manuscript. Built independently from the pipeline contents, then cross-checked against Claude Code's list.

**Method:** inventoried all 203 scripts under `bioinf_redo/`, split into those inside SRP (101) vs outside SRP but still in `bioinf_redo/` (102); traced detectable output-writes; mapped manuscript figures/tables to generators. Where I could not execute or trace a path with certainty, I say so.

**Headline finding:** SRP is **NOT currently self-contained.** 102 scripts the manuscript depends on live *outside* SRP, in sibling `bioinf_redo/` folders. Several stages have their authoritative/working versions outside SRP, not inside it. A repo built from SRP as-is would be missing the generators for at least Table S8, Figure 7 (non-canonical GNA), most of the Stage-2 novelty figures, the network figure (Fig 3) export code, and the supplementary-table generator.

---

## A. Confirmed GAPS — scripts outside SRP that the manuscript needs (pull these in)

| Manuscript item | Generator(s) — location OUTSIDE SRP | Severity |
|---|---|---|
| **Table S8 + Figure 7** (non-canonical genome-neighbourhood: H5Pro, H6Gln, UstYa/cutinase context) | `4_genome_neighbourhood/` — **17 scripts** (`1_collect_accessions.py` … `15_all_groups_neighbourhood.py`, `rescore_bgc.py`, `pfam_g46n_retry.py`). SRP only contains `4_genome_neighbourhood_canonical/`, which holds **output TSVs but no scripts.** | **HIGH** — entire analysis stage, no code in SRP |
| **Stage-2 position vectors** (Fig 5 vectors, the 11-position encoding) | `2_canonical_analysis/extract_position_vectors.py` + `run_extract.sh` | HIGH |
| **Novel-state enumeration** (67 novel states, Table S4) | `2_canonical_analysis/novel_variants.py`, `find_candidate_groups.py`, `qc_candidate_groups.py` | HIGH |
| **Aromatic-compensation check** (Phe65/Trp68/Phe227) | `2_canonical_analysis/test_aromatic_compensation.py`, `test_phe227_compensation.py`, `run_aromatic_compensation.sh` | MED |
| **Supplementary table generator** (Table S3 and the curated-set sub-tables) | `2_canonical_analysis/supplementary/generate_tables.py` (writes `table_DCT_DHICA.tsv`, `table_oMP_*.tsv`, `table_Fusarium_Pseudomonas_convergence.tsv`, etc.) | HIGH |
| **Figure 3 network** (Foldseek landscape — taxonomy/clade/cys mapping, Cytoscape export) | `2_canonical_analysis/visualisation/` — `network_graph_v3.py`, `export_cytoscape.py`, `fetch_taxonomy.py`, `taxonomy_heatmaps.py`, `cooccurrence.py`, `make_vecfig.py`/`vecfig2.py`. (SRP's `5_foldseek/pools/` has clustering + `make_node_table.py` but not the figure-rendering/export code.) | HIGH |
| **HMM-vs-structure comparison (canonical)** — Fig 5 concordance, Table S5 | `2_canonical_analysis/extract_lost_vectors.py`, `validate_vectors.py`, `make_summary_table.py` | MED |
| **Metal3D on DCT/DHICA** (Intro Zn/Cu ambiguity, supplementary) | `2_canonical_analysis/supplementary/metal3d_dct_dhica/` (4 scripts) | LOW |

## B. Present and authoritative INSIDE SRP (good — keep)

- **Filtering & pools** (Fig 2; the 21,893/1,060/9,116 split): `1_filtering/assign_pools.py`, `canonical_criteria` checks, `validation/ppo_core_check.py`.
- **M7 benchmark** (the authoritative Acc 90.5% / F1 0.92): `1_filtering/core_helix_filter/benchmark_eval.py`, `core_helix_check.py`, `run_benchmark_core.sh`, `plot_threshold_grid.py`; plus `1_filtering/benchmark/` (7-method comparison, Table S1: `run_methods.py`, `summarize.py`, `select_benchmark.py`). **Note:** the stale `benchmark_core_summary.txt`/`_compare.tsv` here are pre-widening — quarantine or regenerate.
- **Seed reproducibility** (Methods SD 0.05/0.14/1.85 Å): `1_filtering/seed_variability/` — complete (scripts + results). ✅
- **Discarded-pool characterization** (answer to ARK comment 13): `1_filtering/pool_summary/plot_discarded_characterization.py`. ✅
- **Stage-3 non-canonical** (117/667/276 tiers; Figs 6, 8): `3_noncanonical_analysis/classify_sites.py`, `stage3_extract_align.py`, `add_coord_dists.py`, `generate_patterns.py`, `run_stage3_extract.sh`. ✅ (SRP copy is the v2 M7-widened version per CLAUDE.md.)
- **Fig 10 metal-prediction plots** (AF3 vs AllMetal3D vs PMM): `3_noncanonical_analysis/allmetal3d/` plotting scripts + `1_filtering/metal3d_noncanonicalpool/run_metal3d.py` are present. ⚠️ but see gap C2.
- **Foldseek clustering** (Fig 3 data layer, TM thresholds): `5_foldseek/` + `5_foldseek/pools/run_all_vs_all.sh`, `make_node_table.py`. ✅
- **Domain analysis** (Pfam/Chainsaw, §Domain annotation): `6_domain_analysis/` + `run_chainsaw_batch.py`. ✅

## C. ORPHANS — outputs in the manuscript with NO generating script anywhere

1. **Loose supplementary tables at SRP root** — `table_s3_characterised_ppos.tsv`, `table_gna_variant_groups.tsv`, `table_gna_highlights.tsv`, `table_bgc_pfam_markers.tsv`, `table_s5_hmm_nc_agreement.tsv`, `novelty_enumeration_compact.tsv`. **None appear by name in any script.** `generate_tables.py` (outside SRP) writes *differently-named* TSVs (`table_oMP_*`, `table_DCT_DHICA`, …), so these six look **hand-assembled or renamed** from other outputs. A reviewer/deposit will ask how they were built. → decide: regenerate via a documented script, or add a README note stating they are manually curated from named sources.
2. **AllMetal3D/AF3 merge step only** (Fig 10 data): Claude Code reports `step1_align.py`, `step2_metal3d.py`, `merge_results.py` were **lost from Olivia and not local** — these produced the AF3-vs-metal3d comparison TSVs (`af3_vs_m3d_comparison.tsv`, `af3_vs_metal3d_comparison.tsv`). **Correction:** the PinMyMetal per-site tables are NOT orphaned — `build_pmm_per_site.py` IS present in `3_noncanonical_analysis/allmetal3d/` and generates `pmm_per_site.tsv`/`pmm_feconi_per_site.tsv`. The metal3d compute driver `1_filtering/metal3d_noncanonicalpool/run_metal3d.py` is also present. So the genuine gap is narrow: only the **AF3↔metal3d alignment/merge step** is missing; PMM generation and the metal3d run are covered. Fig 10 plotting scripts are all present. → recover/rewrite only the align+merge step, or document it as "results provided; upstream merge script not retained."

## D. Cross-check against Claude Code's list

| Claude Code claim | My independent finding | Agree? |
|---|---|---|
| Stage-4 non-canonical GNA outside SRP (biggest gap) | Confirmed — 17 scripts, 46 files in `4_genome_neighbourhood/`; SRP has canonical-only, no scripts | ✅ |
| AllMetal3D rerun scripts lost | Partially — the AF3↔metal3d **align/merge** step (step1/step2/merge_results) is absent; but `build_pmm_per_site.py` and `run_metal3d.py` ARE present, so only the merge step is the true gap | ⚠️ narrowed |
| Loose SRP-root tables have no generating script | Confirmed — none referenced by name in any script | ✅ |
| `detect_thioether.py`, `extract_position_vectors_super.py`, `/tmp/build_network.py` lost | Partially — `extract_position_vectors.py` (non-`_super`) IS outside-SRP and recoverable; `detect_thioether.py` not found in SRP tree (Fig 3 Type1/2 colouring); `build_network.py` not present | ✅ (mostly) |
| Crystallography files (Fig 9, Table S9) | **Superseded by you** — deposited to PDB, not a repo concern | n/a |
| GitHub repo stale/lacking | Confirmed — 1 public repo, last push 2026-02-12; you're rebuilding from SRP instead | ✅ |

## E. Recommended actions before cutting a release

1. **Get a durable backup first** (SRP/`bioinf_redo` is currently the only complete copy). Non-negotiable before any file moves.
2. **Pull section-A scripts into SRP** — especially the whole `4_genome_neighbourhood/` folder and the `2_canonical_analysis/` working scripts + `visualisation/` + `supplementary/generate_tables.py`. This is the bulk of the work.
3. **Resolve the two orphan classes (C1, C2)** — either regenerate with a script or add explicit README provenance notes.
4. **Quarantine superseded files** — `*.prewidened.bak`, `*.preM7.bak`, `3_noncanonical_analysis_old/`, `*_old.*`, `*.preFECONI.*` into a `superseded/` folder or `.gitignore` them, so no one recomputes wrong numbers.
5. **Include a cleaned CLAUDE.md as README/METHODS** — it holds the authoritative provenance (M7 numbers, threshold history) and won't travel otherwise.
6. **Re-run one end-to-end check** on the assembled repo to confirm the headline numbers regenerate before publishing the DOI.

## F. Caveats on this audit
- I traced scripts by inventory + output-write grep + CLAUDE.md, **not by executing the pipeline.** Path-building via variables means some script→output links are inferred, not proven. A definitive check requires running each stage.
- "Authoritative version" calls (SRP vs outside) are based on directory structure and CLAUDE.md notes; where both a SRP and an outside copy exist (e.g. `classify_sites.py`, `crystal_annotation_heatmap.py`), confirm which is newer before copying — do not assume SRP is always the latest.
- The manuscript's **32,753 vs 32,069** count discrepancy (flagged earlier) is a text issue, not a script gap.