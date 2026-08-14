# PPO quality-filtering method benchmark

Evaluates how well different structural methods place each structure into
**canonical / non-canonical / discarded**, so the pool-filtering can be trusted.
Geometry (Cu-Cu, His, pLDDT) separates canonical from non-canonical; the hard,
fragile decision is **core present vs not** (non-canonical vs discarded). This
benchmark stress-tests that decision with six orthogonal methods on a stratified
set spanning the boundary, including known positive/negative controls.

## Methods scored per structure
- **M1** Foldseek global TMalign + query-backbone helix check (current pipeline).
- **M2** Foldseek reference-normalized TM (`ttmscore`) — removes the query-length
  normalization artifact that inflates short half-bundle fragments.
- **M3** PyMOL `super` vs the 9 references — best CA coverage + RMSD.
- **M4** PyMOL `cealign` vs the 9 references — best CA coverage + RMSD.
- **M5** Biotite reference-free SSE — long-helix count + compact-bundle test.
- **M6** Combined: intrinsic bundle (M5) AND a PPO-identity gate (M1 or M2).

## Benchmark strata
Controls with known answers: the 9 reference cores, the Microbispora
non-canonical anchor `A0A8H9LF69` (must stay non-canonical), characterized PPOs.
Boundary cases to inspect: the 64 canonical helix-check disagreements, the 74
non-canonical scored 4/4 at qTM<0.3 (possible false positives / junk), full-length
discarded (possible false negatives / real), short + no-hit discarded (expected junk).

## Run
```
sbatch run_benchmark.sh
```
Outputs: `benchmark_set.csv`, `benchmark_results.tsv` (all methods + geometry
context), `benchmark_summary.txt` (per-stratum call rates, control checks,
pairwise agreement), `benchmark_disagreements.tsv` (the inspection queue —
structures where methods disagree).

## Inspect
Open the disagreement queue in PyMOL, label each canonical / non-canonical / junk,
then pick the method (or threshold) that best reproduces the labels.
