#!/usr/bin/env python3
"""
17_recompute_gna_highlights.py - regenerate and audit the genome-neighbourhood
highlights table (supplementary_tables/table_gna_highlights.tsv).

That table was originally assembled by hand, so its numbers had no generating
script and could not be checked. This script recomputes every row from the
pipeline outputs and prints the recomputed value beside the published one, so
any drift is visible rather than silent.

WHY A SPEC RATHER THAN ONE UNIFORM RULE
The eight rows are not homogeneous. They differ in three ways that a single
rule cannot capture, and getting any of them wrong changes the answer:

  1. GROUP DEFINITION. Five rows are novelty-enumeration substitution groups.
     G46N is NOT: it is the 187-structure `table_oAPO_helix` set (Asn46 on a
     helix, the characterised o-aminophenol oxidase state). Gly46_N has no row
     in novelty_enumeration.tsv at all, because Asn at Gly46 is a known state
     rather than a novel one, and it appears in no canonical GNA output.
     Defining it the obvious way instead - position_vectors.csv where
     Gly46 == 'N', which is how the surviving pfam_g46n_retry.py does it -
     gives 1047 carriers, 592 with context, and 24-51% instead of 71%.
  2. DENOMINATOR. G46N's frequency is over ACTINOBACTERIAL carriers only.
     The other bacterial rows are over all carriers with genome context.
  3. FREQUENCY BASIS. The two non-canonical rows read a precomputed Pfam
     co-occurrence frequency, and even those disagree with each other:
     H5Pro's 0.404 is over the full group (72/178) while H6Gln's 0.951 is
     over queries with data (39/41).

Each row below therefore carries its own explicit definition.

Run:  python3 17_recompute_gna_highlights.py [--write]
      --write updates supplementary_tables/table_gna_highlights.tsv in place.
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

csv.field_size_limit(10 ** 9)

HERE = Path(__file__).resolve().parent          # 4_genome_neighbourhood
ROOT = HERE.parent                              # repo root
CANON = HERE / 'canonical'
SUPP = ROOT / '2_canonical_analysis' / 'supplementary'
OUT = ROOT / 'supplementary_tables' / 'table_gna_highlights.tsv'


def tsv(path, delim='\t'):
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def load_taxonomy():
    return {r['accession']: r for r in tsv(ROOT / 'taxonomy_lookup.csv', ',')}


def neighbourhood_products(path):
    """query_accession -> set of lowercased flanking product names."""
    per_q = defaultdict(set)
    for r in tsv(path):
        per_q[r['query_accession']].add((r.get('product') or '').lower())
    return per_q


def hits(per_q, accs, markers):
    """Carriers in `accs` with at least one marker substring among flanking products."""
    return {q for q in accs
            if any(any(m in prod for prod in per_q.get(q, ())) for m in markers)}


# ---------------------------------------------------------------- row specs --

# Canonical novelty groups: n from novelty_enumeration, membership from the
# position vectors, context and markers from the bacterial GNA run.
NOVELTY_ROWS = [
    # label   position  residue  neighbourhood            markers
    ('G46Y',  'Gly46',  'Y', 'Chaplin / exopolysaccharide biosynthesis',
     ['chaplin', 'exopolysaccharide'],
     'Chaplin; exopolysaccharide biosynthesis polyprenyl glycosylphosphotransferase'),
    ('G46I',  'Gly46',  'I', 'Melanin copper-chaperone MelC1',
     ['tyrosinase cofactor', 'melc'],
     'Tyrosinase cofactor MelC1; tyrosinase MelC2'),
    ('G46T',  'Gly46',  'T', 'Tryptophan / scytonemin biosynthesis',
     ['anthranilate phosphoribosyltransferase', 'tryptophan synthase'],
     'Anthranilate phosphoribosyltransferase; tryptophan synthase (alpha/beta)'),
    ('V218E', 'Val218', 'E', 'Protocatechuate catabolism',
     ['4-carboxymuconolactone decarboxylase', 'gntr'],
     '4-carboxymuconolactone decarboxylase; GntR family transcriptional regulator'),
    ('N205R', 'Asn205', 'R', 'Bacillithiol biosynthesis',
     ['bsha', 'glucosaminyl l-malate synthase'],
     'N-acetyl-alpha-D-glucosaminyl L-malate synthase BshA; sugar transporter'),
]

PUBLISHED = {
    'G46N':  (187, 97, 71),
    'G46Y':  (68,  31, 58),
    'G46I':  (65,  32, 34),
    'G46T':  (279, 10, 60),
    'V218E': (99,  12, 67),
    'N205R': (143, 14, 57),
    'H5Pro': (178, 95, 40),
    'H6Gln': (43,  39, 95),
}


def main():
    tax = load_taxonomy()
    pv = tsv(SUPP.parent / 'position_vectors.csv', ',')
    novelty = tsv(SUPP / 'novelty_enumeration.tsv')
    rows = []

    # ---- G46N: the oAPO_helix group, actinobacterial denominator ------------
    oapo = {r['accession'] for r in tsv(SUPP / 'table_oAPO_helix.tsv')}
    per_q = neighbourhood_products(HERE / 'groups' / 'G46N' / 'neighbourhoods.tsv')
    ctx = oapo & set(per_q)
    acti = {q for q in ctx if tax.get(q, {}).get('phylum') == 'Actinomycetota'}
    markers = ['2-amino-3,7-dideoxy', '3-dehydroquinate synthase', '3-amino-4-hydroxybenzoic']
    h = hits(per_q, acti, markers)
    rows.append(['G46N', len(oapo), len(ctx), 'Aminoshikimate pathway',
                 round(100 * len(h) / len(acti)) if acti else 0,
                 '2-amino-3,7-dideoxy-D-threo-hept-6-ulosonate synthase; '
                 '3-dehydroquinate synthase II; 3-amino-4-hydroxybenzoic acid synthase',
                 f'{len(h)}/{len(acti)} actinobacterial'])

    # ---- canonical novelty groups ------------------------------------------
    # NOTE: the bacterial run's per-group accession lists are NOT in the repo
    # (canonical/target_accessions.tsv covers only the second-sphere run, whose
    # subgroups are G46N_helix / E195R / H230Y_5pos / F227W / F65L_W68A_*).
    # So per-CARRIER counts cannot be recomputed for these five rows; only the
    # product x group OCCURRENCE matrix survives. For four of the five the two
    # coincide and the published value reproduces exactly. Where a marker's
    # occurrence count exceeds n_context the two demonstrably differ, and the
    # row is reported as not reproducible rather than silently rounded down.
    matrix = tsv(CANON / 'bacterial' / 'full_products_by_group.tsv')
    bact_summary = {r['subgroup']: int(r['n_queries'])
                    for r in tsv(CANON / 'bacterial' / 'summary_by_group.tsv')}
    for label, pos, res, nbh, markers, marker_txt in NOVELTY_ROWS:
        n = next((int(r['n']) for r in novelty
                  if r['position'] == pos and r['residue'] == res), None)
        grp = f'{pos}_{res}'
        n_ctx = bact_summary.get(grp, 0)
        best, best_prod = 0, ''
        for r in matrix:
            if any(m in r['product'].lower() for m in markers):
                v = int(r.get(grp) or 0)
                if v > best:
                    best, best_prod = v, r['product']
        if n_ctx and best > n_ctx:
            freq, basis = None, (f'NOT REPRODUCIBLE: "{best_prod[:28]}" has {best} '
                                 f'occurrences across only {n_ctx} loci; per-carrier '
                                 f'counts not in repo')
        else:
            freq = round(100 * best / n_ctx) if n_ctx else 0
            basis = f'{best}/{n_ctx} loci ("{best_prod[:34]}")'
        rows.append([label, n, n_ctx, nbh, freq, marker_txt, basis])

    # ---- non-canonical rows: precomputed Pfam co-occurrence -----------------
    for label, pfam_file, pfam_id, nbh, marker_txt, n_src in [
        ('H5Pro', 'his5pro_cooccurrence_pfam.tsv', 'PF11807', 'UstYa-type BGC',
         'Mycotoxin biosynthesis protein UstYa (PF11807)', 'PRO'),
        ('H6Gln', 'h6gln_cooccurrence_pfam.tsv', 'PF24864', 'Conserved non-BGC locus',
         'DUF7730 (PF24864); RIO domain (PF01163); zinc carboxypeptidase; sulfotransferase', 'GLN'),
    ]:
        pf = {r['pfam_id']: r for r in tsv(HERE / pfam_file)}
        row = pf.get(pfam_id)
        freq = float(row['freq_of_all_queries']) if row else 0.0
        nq = int(row['n_queries']) if row else 0
        # group size from the stage-3 substitution pattern table
        pat = tsv(ROOT / '3_noncanonical_analysis' / 'nc_patterns_compact.tsv')
        n = next((int(r['n']) for r in pat
                  if r['H1'] == 'HIS' and r['H6' if n_src == 'GLN' else 'H5'] == n_src), None)
        memb = HERE / (f'his5pro_cluster_membership.tsv' if n_src == 'PRO'
                       else 'h6gln_cluster_membership.tsv')
        n_ctx = len(tsv(memb))
        rows.append([label, n, n_ctx, nbh, round(100 * freq), marker_txt,
                     f'{nq} queries, freq_of_all_queries={freq:.3f}'])

    # ---- report -------------------------------------------------------------
    hdr = ['Group', 'n', 'n_context', 'Conserved neighbourhood', 'Freq. (%)', 'Marker proteins']
    print(f'{"Group":7s} {"n (pub)":>14s} {"n_context (pub)":>18s} {"freq% (pub)":>14s}   basis')
    print('-' * 96)
    drift = []
    for r in rows:
        label, n, n_ctx, _, freq, _, basis = r
        pn, pc, pf_ = PUBLISHED[label]
        show = 'n/a' if freq is None else freq
        f = lambda got, pub: f'{got} ({pub}){"" if got == pub else "  <-"}'
        print(f'{label:7s} {f(n, pn):>14s} {f(n_ctx, pc):>18s} {f(show, pf_):>14s}   {basis}')
        for name, got, pub in (('n', n, pn), ('n_context', n_ctx, pc), ('freq', freq, pf_)):
            if got != pub:
                drift.append((label, name, pub, got))

    if drift:
        print('\nDRIFT vs the published table:')
        for label, name, pub, got in drift:
            print(f'  {label:6s} {name:10s} published {pub} -> recomputed {got}')
    else:
        print('\nAll rows reproduce the published table exactly.')

    if '--write' in sys.argv:
        with open(OUT, 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(hdr)
            for r in rows:
                w.writerow(r[:6])
        print(f'\nwrote {OUT.relative_to(ROOT)}')
    else:
        print('\n(dry run; pass --write to update the table)')


if __name__ == '__main__':
    main()
