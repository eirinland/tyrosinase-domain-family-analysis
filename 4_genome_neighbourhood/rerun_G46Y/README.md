# Gly46_Y genome-neighbourhood rerun (2026-08-18)

Why this exists: the original bacterial canonical GNA run wrote only three
aggregate tables and never saved its per-locus `neighbourhoods.tsv` or target
list. Those aggregates count *occurrences*, not distinct carriers, so the
published G46Y figure (58% of 31 loci) could not be verified from surviving
data. This rerun regenerates the per-locus records for that one group.

## What was run

Targets were rebuilt from the `vector` column of
`2_canonical_analysis/position_vectors.csv` (Gly46 == 'Y'), the same definition
`novelty_pipeline.py` uses. That yields **68 carriers**, matching
`novelty_enumeration.tsv` n=68 exactly, of which 43 are bacterial.

Then, unchanged from the original pipeline:

```
2_fetch_genome_context.py                 # UniProt REST -> genome_crossrefs.tsv
3_fetch_neighbourhoods.py --flank 10      # NCBI Entrez  -> neighbourhoods.tsv
```

Run in this isolated directory on purpose: both scripts read and write in their
own parent directory, so running them in `4_genome_neighbourhood/` would have
overwritten the non-canonical GNA data. Use `run_rerun.sh`, which copies the
scripts in before running.

## Result

- 68 targets → 50 with UniProt genome cross-references → **50 neighbourhoods**
  (885 flanking genes), of which **35 bacterial** (34 actinobacterial).
- Per-carrier marker frequencies among the 35 bacterial carriers:

| Definition | carriers | % |
|---|---:|---:|
| chaplin | 21/35 | **60%** |
| exopolysaccharide | 18/35 | 51% |
| chaplin OR exopolysaccharide | 26/35 | 74% |
| chaplin AND exopolysaccharide | 13/35 | 37% |

**The published 58% is corroborated.** It corresponds to chaplin carriers over
loci with context: 18/31 = 58% then, 21/35 = 60% now. Scaling the fresh rate
onto the original denominator gives 18.6 of 31 against the published 18, i.e.
the same underlying rate.

The denominator moved from 31 to 35 because more annotated nucleotide records
are available in NCBI now than when the original run was made. That is expected
drift for an API-backed analysis, and is the reason the raw records are kept
here rather than only a summary.

Note the published label reads "Chaplin / exopolysaccharide biosynthesis", but
the 58% tracks chaplin alone; the union of both markers is 74%.
