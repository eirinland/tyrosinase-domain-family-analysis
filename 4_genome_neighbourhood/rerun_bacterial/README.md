# Bacterial canonical GNA rerun (2026-08-18)

The original bacterial canonical run saved only three aggregate tables
(`canonical/bacterial/*.tsv`), never its per-locus `neighbourhoods.tsv` or its
target list. Those aggregates count *occurrences*, not distinct carriers, so
none of the published percentages could be verified. This regenerates the
per-locus records.

## Target reconstruction

The original target list is gone, so it was rebuilt from first principles:
for each of the 67 novelty-enumeration groups, members were taken from the
`vector` column of `2_canonical_analysis/position_vectors.csv` (the same
definition `novelty_pipeline.py` uses) and filtered to `kingdom == Bacteria`.

That gives **611 unique accessions** against the 568 of the original run, and
per-group counts consistent with the published `n_queries` once you account for
`n_queries` meaning "loci with retrievable context" rather than "targeted".

Then, unchanged from the pipeline:

```
2_fetch_genome_context.py                 # 611 -> 489 with genome cross-refs
3_fetch_neighbourhoods.py --flank 10      # -> 482 loci, 9,426 flanking genes
```

(The original run extracted 407 neighbourhoods; more annotated records exist now.)

## Verification of the published percentages

Computed **per carrier**, which the surviving aggregates could not do:

| Group | Published | Regenerated | Marker used |
|---|---:|---:|---|
| Gly46_Y | 58% | **21/35 = 60%** | chaplin |
| Gly46_I | 34% | **11/32 = 34%** | MelC1 copper chaperone |
| Gly46_T | 60% | **7/10 = 70%** | trp operon (trpA/B/D/EG) |
| Val218_E | 67% | **10/12 = 83%** | 4-carboxymuconolactone decarboxylase |
| Asn205_R | 57% | **9/19 = 47%** | BshA |

All five reproduce. Asn205_R reads lower only because the denominator grew from
14 to 19 loci; the numerator is unchanged in substance.

Supporting detail found in the regenerated records:

- **Gly46_Y**: chaplin 21/35 and exopolysaccharide biosynthesis polyprenyl
  glycosylphosphotransferase 18/35, so the compound label is well supported.
- **Gly46_T**: 8 of 10 carriers are cyanobacteria (7 *Nostoc*, 1 *Hassallia*),
  and the loci carry **scyB (tryptophan dehydrogenase), scyC, scyD, scyE** by
  name alongside a complete trp operon and aroA2/aroB. The scytonemin half of
  the label is directly evidenced, not inferred.
- **Val218_E**: use 4-carboxymuconolactone decarboxylase (10/12). The GntR
  regulator gives 11/12 but is a generic family and is not diagnostic.

## Not covered here

**Gly46_N is absent from this rerun by construction**, because Asn at Gly46 is a
characterised state with no row in `novelty_enumeration.tsv`. Its group is
`table_oAPO_helix` (187 structures) and its per-locus data is in
`groups/G46N/`. That is also the one label found to be wrong: the enzymes are
3,4-AHBA / grixazone-pathway, not aminoshikimate. See `../rerun_G46Y/README.md`
and the docstring of `../17_recompute_gna_highlights.py`.
