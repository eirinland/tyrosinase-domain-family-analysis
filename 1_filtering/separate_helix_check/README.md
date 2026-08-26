# Separate-helix requirement for the M7 core test

**Question.** M7 tests the four core helices independently: helix *i* is accepted when
one of its reference His anchors finds a query Cα that is close (≤4 Å), locally
helical and confident (pLDDT ≥70). Nothing forbids two anchors being satisfied by the
*same* piece of query backbone. What would the non-canonical pool look like if the
four anchors were also required to land on four **separate** helices?

## Run order

```bash
python3 cache_ca.py             # OPTIONAL: ca_cache.npz is deposited; this rebuilds it
python3 calibrate_segments.py   # choose the Ca-only "same helix" definition
python3 separate_helix_check.py # apply it; writes summary.txt + the two TSVs
```

`ca_cache.npz` (9 MB, Cα coordinates and pLDDT for the 1,256 structures in scope) is
deposited so the check re-runs without the 2.55 GB AF3 archive; `cache_ca.py` regenerates
it from the archive if needed. `calibrate_segments.py` reads the reference structures in
`1_filtering/foldseek/` and `B2ZB02_taxID_1404_model.cif`, both in the repository.

## Calibrating "the same helix" (calibrate_segments.py)

The rule is only meaningful if a textbook PPO passes it, so the segment definition was
picked by reference control rather than assumed. Candidate definitions differ in the
Cα i→i+3 / i→i+4 windows and in how many covering helical turns a residue needs:

| definition | 4 separate core helices on the 5 Cu-bearing refs? |
|---|---|
| widened windows, ≥1 turn | **4/5** – merges PmTYR a3 and a4 through the intervening loop |
| widened windows, ≥3 turns | **5/5**, segment bounds match the hard-coded core-helix ranges |
| pre-widening windows, ≥1 turn | 5/5, but segments are grossly over-extended (PmTYR a1 → 10–53) |
| strict windows, ≥3 or ≥5 turns | 3/5 and 0/5 – anchors fall off the ends of the segments |

**Adopted:** a residue is helical when covered by **≥3** genuine helical turns under the
FINAL widened M7 windows (4.0–6.4 Å / 4.8–8.2 Å); maximal runs of consecutive helical
residues are the query's helices; an anchor joins the run holding its query residue or
a neighbour (±1, M7's own cap fallback). On PmTYR/B2ZB02 this gives 17–47, 66–82,
205–210, 227–244 against the hard-coded 34–46, 65–83, 203–211, 226–244.

The looser "≥1 turn" variant is kept in the output only to show what a naive definition
does: it flags 693/1,060 and fails its own reference control (it would discard PmTYR).

## Result

**9 of the 1,060 non-canonical structures have two anchors on one helix.**
Non-canonical 1,060 → **1,051**, discarded 9,116 → **9,125**; canonical 21,893 is
untouched (it never uses M7). Sensitivity: bridging helical runs across gaps of ≤2
residues gives 11 instead of 9.

Colliding pairs: a2/a4 ×4, a1/a2 ×3, a3/a4 ×1, a1/a2 + a3/a4 ×1. **Four are cross-lobe**
(one helix serving a CuA and a CuB slot) — those are exactly the four already shown as
false positives in SI panel E.

The 9 are worse than the pool on every axis: qTM median 0.405 (pool 0.490), length
median 337 (pool 429), 6/9 classified `no_cu` at stage 3.

Benchmark (241 hand-labelled structures): the rule costs **nothing**.
M7 alone TP=141 FP=0 FN=21 TN=79; M7 + separate-helix identical. No hand-labelled true
core collides, and there were no false positives left to remove (M7 precision is
already 1.00). *(TP=141 vs the published 139 is the recomputation artefact noted in
`si_discarded_examples/recompute_helicity_widened.py`: the stored table keeps one
candidate anchor per helix, so widening cannot promote a different anchor of the same
helix. Both sides of the comparison use the same recomputation, so the rule's delta of
0 is unaffected.)*

## Downstream, if adopted

- stage 3 tiers 117/667/276 → **116/665/270**
- H6Gln GNA group 43 → 42 (loses W7F1D8, *Bipolaris*, `bgc_context = none`; not one of
  the 39 *Colletotrichum*, so the conserved-locus result is unchanged)
- Foldseek pool clustering staged set 22,953 → 22,944
- no canonical-pool, stage-2, or characterized-PPO number changes

## The residue-index proxy was not reliable

SI panel E used `|qres_i − qres_j| ≤ 10` as a stand-in for "same helix". Against the
calibrated definition it gets 7 of the 9 right, raises 21 false alarms (anchors close in
sequence but on genuinely different helices), and misses 2 real collisions where the
anchors are 15 residues apart on one long fused helix (A0A8H5WCB4 *Fusarium* qTM 0.745,
W7F1D8 *Bipolaris* qTM 0.668). If panel E keeps a count in its caption, use 9 (or 28
described explicitly as "anchors within 10 residues", which is a different claim).

## Files

- `separate_helix_nc.tsv` — all 1,060 non-canonical: anchor residues, assigned helical
  segment per definition, collision flags, colliding pair, cross-lobe flag, stage-3 tier,
  taxonomy
- `separate_helix_bench.tsv` — the same for the 241 benchmark structures, plus label
- `summary.txt` — the full printed report, including the reference control table
- `ca_cache.npz` — Cα coordinates + pLDDT for the 1,256 structures in scope
