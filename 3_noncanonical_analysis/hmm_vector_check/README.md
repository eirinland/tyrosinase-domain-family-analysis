# HMM vs structure at the six canonical His positions, per structure

Rebuilds the data behind Table S6 (`3_noncanonical_analysis/hmm_his_comparison.tsv`,
a summary with no surviving generator) so that individual structures and groups can
be looked up, and answers two questions: does the "complete six-position vectors
correct" set cover the non-canonical cases the text highlights, and are the H5Pro
structures assigned correctly by HMM. All counts below use the recommended convention
(502 of 1,060); the named-structure verdicts are identical under every convention.

## Run order

```bash
python3 ../../2_canonical_analysis/hmm/build_alignment.py   # once; needs pyhmmer
python3 hmm_vector_check.py    # per-position comparison  -> hmm_vs_structure_nc.tsv
python3 highlighted_cases.py   # named structures, groups, H5Pro -> highlighted_summary.txt
```

`build_alignment.py` aligns the deposited `query.fasta.gz` (32,753 sequences) to
PF00264.26 and writes `hmm_match_columns.tsv.gz` plus `reference_map.tsv`; it uses
pyhmmer, which embeds the same HMMER3 code as the `hmmalign` binary named in the
Methods. After that step this check needs **only pandas** -- it reads the match-state
table, so the claim is verifiable from the deposit without any alignment software.
`WITH_RESNUMS=1` re-aligns the 1,061 sequences to also recover residue numbers, which
the compact table does not carry (needed only for the stricter variant below).

Other inputs: the structural vectors in `3_noncanonical_analysis/noncanonical_analysis.tsv`
and the pool assignment in `1_filtering/final_pools/`.

B2ZB02 ties the six positions to match states: His42→11, His60→25, His69→34,
His204→174, His208→178, His231→201.

## Reproduction fidelity

Per-position counts reproduce the published Table S6 **exactly at five of six
positions**. CuB_His1 differs by one structure (agree 233 vs published 232).
`n_substituted` matches at all six, which confirms the table counts positions where
the structural vector shows a mapped non-His residue (unmapped `---` excluded).

Complete six-position vectors, published as 504/1,060:

| convention | count | |
|---|---:|---|
| string equality of the two rendered vectors | 512 | **not quotable**, see below |
| every position recovered, residue **type** only | 505 | |
| every position recovered, **same residue** | **502** | **recommended** |
| all structurally mapped positions correct (unmapped ignored) | 547 | |
| all substituted positions correct | 581 | |

The three top rows differ for two unrelated reasons, and separating them is what
settles the number:

- **512 → 505 (7 structures).** Comparing the two rendered vector *strings* scores an
  unmapped structural position (`---`, printed `-`) against an HMM deletion (`-`) as a
  match, in 8 position-instances across 7 structures. Nothing was recovered at those
  positions, and the per-position rule behind Table S6 never lets a deletion agree, so
  the string comparison contradicts the table it is supposed to summarise. It is kept
  in the output only to explain the discrepancy.
- **505 → 502 (3 structures).** Three position-instances where the HMM reports the right
  residue *type* at the wrong residue *number*: ILE151 vs I152, HIS116 vs H121, and
  HIS102 vs H88 — the last 14 residues out. These are genuine misalignments, exactly
  what a position-recovery claim must not count as successes.

**For the manuscript:** quote **502 of 1,060 (47.4%)**, with a position counted as
recovered only when the HMM column reports the same residue *at the same residue
number*, and both a deletion and an unmapped structural position counted as failures.
The per-structure table `hmm_vs_structure_nc.tsv` is the deposited evidence; the
aggregate `hmm_his_comparison.tsv` cannot support a joint claim, because a joint count
needs each structure's six calls together. The published 504 remains unreproducible,
but it sits 2 above the defensible number rather than 8 below the loosest one.

## Q1: does the recovered set cover the highlighted cases? No, and it fails on the ones that matter

Named structures: **7 of 11 covered**.

| structure | where | struct → HMM | in the set? |
|---|---|---|---|
| A0A9N8JDS5 | Fig 6A binuclear | HHHHHY → HHHHHY | yes |
| A0A9P7QPY3 | Fig 6B binuclear | HHHHYH → HHHHYH | yes |
| A0A6G0QFF2 | Fig 6C binuclear | HHH**E**HH → HHH**-**HH | **no**, H4 Glu gapped out |
| G3JPN7 | Fig 7A, H5Pro rep | HHHHPH → HHHHPH | yes |
| H1UWR0 | Fig 7B, H6Gln rep | HHHHHQ → HHHHHQ | yes |
| A0A0V1MUW3 | Fig 7C, Trichinella | HHHYTQ → HHHYTQ | yes |
| A0A8J9R8H5 | phomQ1 | HHHHPH → HHHHPH | yes |
| A0A142I737 | phomQ1' | HHHHPH → HHHHPH | yes |
| H3GEM4 | Fig 8A degenerate | ---Q-Y → --YQHY | **no** |
| A0A8H9LF69 | Fig 8B, PDB 32AE | YHQMQA → **-**HQMQA | **no**, H1 Tyr gapped out |
| H2KPL1 | Fig 8C degenerate | NHHNHQ → **-**HHNHQ | **no**, H1 Asn gapped out |

Groups (all group sizes reproduce the counts stated in the text):

| group | n | full vector correct | defining position(s) correct |
|---|---|---|---|
| H5Pro, dominant mononuclear | 178 | 174 | 175 |
| H6Gln | 43 | 41 | 41 (2 gapped) |
| Trichinella, CuB degenerate | 14 | 14 | 14 |
| binuclear, His at all six | 69 | 58 | – |
| H4 Glu, 2 distant species | 4 | **0** | **0** (all four gapped) |
| Phytophthora, CuA replaced | 20 | **0** | **0** (all 20 gapped at H1–H3) |
| largest degenerate (YVYQHY) | 34 | **0** | **0** (all 34 gapped) |
| Clonorchis / Opisthorchis | 6 | **0** | **0** (all six gapped) |
| Microbispora actinobacteria | 8 | 1 | – |

So the split is clean: the highlighted **binuclear and mononuclear** groups are inside
the 504, and **every highlighted degenerate group is essentially absent from it** —
including the Microbispora protein with the crystal structure and the H4 Glu group.
The failure mode is not a wrong residue, it is a **gap**: hmmalign deletes the match
state rather than mis-assigning it.

The gradient behind the pool-level number:

| site class | complete vectors correct |
|---|---|
| binuclear | 84/117 (72%) |
| mononuclear | 351/667 (53%) |
| degenerate (no_cu) | 67/276 (24%) |

| substituted positions | correct |
|---|---|
| 1 | 333/511 (65%) |
| 2 | 28/107 (26%) |
| 3 | 58/128 (45%) |
| 4 | 20/64 (31%) |
| 5 | 5/107 (5%) |
| 6 | 0/29 (0%) |

## Q2: H5Pro is assigned correctly by HMM

Of the 178 H5Pro structures, **175 read the proline correctly at H5** and none are
gapped there; 174 have the complete six-position vector right. Both annotated
halogenases (phomQ1 A0A8J9R8H5, phomQ1' A0A142I737) are correct at all six.

The four misses:

| accession | struct → HMM | error | species |
|---|---|---|---|
| A0A4Q4NQ20 | HHHHPH → HHHMRH | H4 His→Met, H5 Pro→Arg | *Alternaria alternata* |
| A0A8H7BB98 | HHHHPH → HHHMHH | H4 His→Met, H5 Pro→His | *Alternaria burnsii* |
| A0A8H6FIQ4 | HHHHPH → HHHHSH | H5 Pro→Ser | *Letharia lupina* |
| A0A9W9IRY5 | HHHHPH → -HHHPH | H1 gapped | *Penicillium capsulatum* |

This is expected: H5Pro carries a single substitution in an otherwise canonical
domain, which is where profile alignment works. H5Pro is therefore **not** an example
of HMM failure; the degenerate groups are.

## Noted while doing this

- The results text cites this table as **Table S6** (Table S5 is the canonical-pool
  vector-position table). The "18–81%" range is the `Agreement_pct_incl_gaps` column
  (17.7–81.1), i.e. agreement over all substituted positions including the gapped
  ones, not the `Agreement_pct` column (60.2–96.6).
- **H3GEM4**, the Figure 8A representative, has four structurally unmapped positions
  (`- - - Q - Y`, Cu–Cu 40.4 Å) and is therefore a member of neither exact-vector
  group it is cited alongside: not of the 20 Phytophthora `ALNHHH` structures (which
  number exactly 20 without it) and not of the 34-structure `YVYQHY` group (exactly 34
  without it). Worth re-checking which group the panel is meant to represent.

## Files

- `hmm_vs_structure_nc.tsv` — 1,060 rows: both vectors, per-position residue, HMM and
  structural residue numbers, agree / gapped / substituted flags, and the four
  complete-vector definitions
- `highlighted_cases.tsv`, `h5pro_cases.tsv` — the named structures and the full H5Pro group
- `summary.txt`, `highlighted_summary.txt` — the printed reports
- inputs come from `2_canonical_analysis/hmm/` (`hmm_match_columns.tsv.gz`,
  `reference_map.tsv`), regenerated by `build_alignment.py` from `query.fasta.gz`
