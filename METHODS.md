# Methods and provenance

Authoritative description of the classification rules, thresholds, and the reasoning
behind them. Where an approach was tried and rejected, that is recorded too: several
retired methods produced numbers that appear in earlier drafts, and this file is the
record of why they were replaced.

## Dataset

32,069 AlphaFold 3 models of PF00264 (PPO fold) sequences, predicted with 2 Cu ions.
Naming: `<accession>_taxID_<taxid>_model.cif`.

Insects and arthropods are **not** part of this family and were excluded. O61363 (squid
haemocyanin, 2,896 aa) is in the PF00264 HMM search but was too large for AF3 prediction
and is absent from the 32,069.

Reference set: 9 chain-A structures spanning the family, in `1_filtering/foldseek/`:
PmTYR (bacterial typical, own crystal structure), 8BBR (*V. spinosum*, bacterial atypical),
2Y9W (*A. bisporus*, fungal), 1BT3 (*I. batatas*, plant catechol oxidase), 5CE9
(*J. regia*, plant), 5M8L (human DHICA oxidase), 1JS8 (squid haemocyanin, trimmed),
I3D139 (archaeal, AF3 fill-in), A0A9N8ELP9 (oomycete, AF3 fill-in).

## The classification rule (authoritative)

Applied by `1_filtering/assign_pools.py`, producing `1_filtering/final_pools/`.

1. **Canonical** iff all of: 2 Cu present; Cu–Cu distance 2.8–5.5 Å; 6 distinct His NE2
   within 3.5 Å of a Cu; per-His **Cα** pLDDT ≥ 70. → **21,893**
2. Otherwise, if the structure still has **≥6 coordinating His** → **discarded** (278).
   These have a full six-His site but fail on geometry or confidence. They are
   deliberately *not* rescued: the non-canonical pool is defined biologically as valid
   PPO folds that have **lost** the canonical six-His site, so a structure retaining six
   His has no analytical home there.
3. Otherwise (**<6 His**), **non-canonical** iff it passes the M7 copper-anchored
   core-helix test; else **discarded**. → **1,060** non-canonical, 8,838 discarded.

Canonical is **geometry-only**: no helix filter, no qTM floor. The 4-helix bundle is
implied by the di-Cu/six-His geometry, which is a stronger and less fragile junk filter
than helix coverage. This is justified empirically: an independent run of the helix check
over the canonical pool found 99.7% (21,829/21,893) pass 4/4, and the 64 that disagree are
geometry-perfect (worst-His Cα pLDDT median 94.3, minimum 70.7). Using the helix check as
a *gate* would have wrongly dropped those 64 valid structures.

Per-His Cα pLDDT is used rather than NE2 pLDDT: NE2 penalises side-chain uncertainty
unfairly, and switching to Cα recovered 301 additional genuine canonicals.

### The M7 core-helix test

`1_filtering/core_helix_filter/core_helix_check.py`. For each <6-His structure that failed
canonical: Foldseek-align against the 5 Cu-bearing references, take the best-qTM reference,
perform a helix-anchored Kabsch superposition followed by ICP onto that reference's four
core helices, then **at the copper-anchored position** require the query's own backbone to
be helical at all four core helices:

- Cα i→i+3 distance ∈ **4.0–6.4 Å**
- Cα i→i+4 distance ∈ **4.8–8.2 Å**
- plus proximity to the reference helix and sufficient confidence

`core_ok` = all four helices accepted. Structures with no result are treated as core-fail.

Thresholds were widened from the original (4.8–6.4, 5.4–7.4) after a 400-combination sweep
over the 241-structure benchmark showed **FP = 0 across the whole grid**. The wider windows
admit 3₁₀ turns and mildly distorted helices without admitting a single false positive, and
rescued 24 structures from discard. The 23 false negatives are predominantly low-pLDDT
structures failing the confidence gate rather than the helicity test.

## Superseded approaches

Recorded because they generated numbers that appear in earlier drafts.

- **4-helix bundle coverage on the Foldseek alignment** (original stage 3). Fundamentally
  broken: it tested the *alignment*, not the query structure. Local alignment clipped the
  terminal helices a1/a4 on near-identical structures (false negatives), and M-states could
  be threaded across a cherry-picked poor reference (false positives). It wrongly excluded
  548 genuine canonical PPOs.
- **Helix-pLDDT filter** (average pLDDT over core helices ≥70). Dropped; the per-His Cα
  pLDDT gate carries the confidence requirement.
- **Trimming to the PPO core plus a second Foldseek pass.** Tested and rejected: Foldseek
  local alignment gives identical aligned regions on full and trimmed structures, and
  trimming caused 272 structures to lose canonical status by cutting coordinating His
  residues at Chainsaw domain boundaries.
- **M1 + cealign + pLDDT-floor rule** (the 2026-06-14 rule, giving 1,106 / 9,070). Replaced
  by M7 because the M1 helix test was reference-frame-free (it tested helicity wherever the
  global alignment mapped, not at the actual Cu site) and the sub-gates were ad hoc. M7 asks
  the biologically correct question directly.

## Stage 2, canonical analysis

Kabsch superposition of the 6 coordinating His NE2 onto PmTYR, mapping 16 reference
positions to the nearest query Cα, with geometric thioether-Cys detection.
Outputs `position_vectors.csv` and `position_vectors_super.csv` (21,893 rows each).

The **Kabsch** method is used for residue-level position mapping rather than Foldseek,
which mis-annotates loop residues and returns gaps at many positions. Foldseek is retained
only for coarse helix coverage in stage 1.

Novelty enumeration yields 67 novel states (50 placed, 17 displaced); within-one-substitution
10,991 (50.2%), within-two 17,697 (80.8%).

## Stage 3, non-canonical analysis

Superposition is **helix-anchored, not His-anchored**: a Kabsch seed on the four core-helix
Cα from the Foldseek cigar, then correspondence-free ICP onto the PmTYR core-helix Cα, then
reading the residue nearest each of the six canonical His anchors. This is why the low-His
tail (0–2 His) aligns correctly; the concern that "Kabsch needs ≥3 His" applied only to a
superseded His-anchored method.

`classify_sites.py` v2 scores all 1,060. A canonical His position whose anchor finds no
query Cα within 3.0 Å is scored as a **divergent position, not a drop**. For 2-His sites
with a coordinatable substitution (Glu/Asp/Cys/Tyr/Met), the substituting residue's nearest
coordinating atom must lie within 4.0 Å of the canonical Cu position, or the site is
downgraded from plausible to divergent. Tiers: **117 binuclear, 667 mononuclear, 276 no-Cu**.

## Stage 5, clustering and network

pLDDT terminal trimming, then Foldseek `easy-cluster` in connected-component mode at TM
0.5/0.7/0.8/0.9. Of the 22,953 pooled structures, 22,952 clustered (one dropped: no window).
At TM 0.8: 1,075 clusters, 866 singletons, largest 9,816.

Two bugs found and fixed here, both of which silently corrupted the network figure:

- **Chain-tag doubling.** Foldseek appends `_<chain>` to the structure name, so staging CIFs
  already named `..._model_A.cif` produced edge names `..._A_A` that failed to join against
  node names `..._A`. The four largest clusters imported into Cytoscape as tiny unstyled
  dots. Edge *counts* were never wrong, only the names. Never stage CIFs whose filename
  already encodes the chain unless you normalise afterwards.
- **`cys_in_shell` fallback.** `C*` is the fallback marker meaning a Cys exists in sequence
  between the two CuA His but its SG is >3.5 Å from the His2 ring, i.e. **not** a real
  thioether. It was being counted as present. Only genuine `C` (SG ≤3.5 Å) counts.
  Overall there are 1,359 `C*` versus 6,151 `C`.

## Known limitations

- 74 non-canonical structures with qTM <0.3 scored 4/4 on the core test and remain
  uninspected; these are candidate false positives.
- 32 full-length discarded structures align across the bundle but read as non-helical;
  these are candidate false negatives.
- Confidence in the canonical pool is high (pure geometry, independently corroborated).
  Confidence at the non-canonical boundary is high but rests on the 241-structure benchmark.
