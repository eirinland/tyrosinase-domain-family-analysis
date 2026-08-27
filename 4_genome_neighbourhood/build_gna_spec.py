#!/usr/bin/env python3
"""
build_gna_spec.py - regenerate the genome-neighbourhood figure spec from the
pipeline's coordinate tables.

For each of the eight figure panels this reads the real +/-10-gene window of a
chosen representative accession, classifies each gene by product keyword, and
writes gna_gene_maps_spec.tsv (consumed by 16_plot_gene_maps.py).

Representative accessions were chosen as carriers of the panel's signature
context gene that also have a retrievable genome neighbourhood. Panel E
(Nostoc, Gly46>Thr tryptophan operon) was recovered via UniProt->NCBI and is
stored in panelE_G46T_nostoc_neighbourhood.tsv.

Usage: python build_gna_spec.py
"""
import csv, os
from collections import defaultdict

WORK = os.path.dirname(os.path.abspath(__file__))

# panel -> (representative accession, species, substitution, locus label,
#           coordinate-source file relative to this folder)
PANELS = {
 "A": ("A0AAX4IJG7", "Colletotrichum destructivum", "CuB His5>Pro", "UstYa-type BGC", "neighbourhoods.tsv"),
 "B": ("A0A066XSE1", "Colletotrichum sublineola", "CuB His6>Gln", "Conserved host-interaction locus", "h6gln_neighbourhoods.tsv"),
 "C": ("A0A0C5G0L4", "Streptomyces sp.", "Gly46>Tyr", "Chaplin / exopolysaccharide", "neighbourhoods.tsv"),
 "D": ("A0ABW8CT20", "Streptomyces sp.", "Gly46>Ile", "MelC1 melanin operon", "groups/G46I/neighbourhoods.tsv"),
 "E": ("A0A2K8SJE4", "Nostoc flagelliforme", "Gly46>Thr", "Tryptophan biosynthesis", "panelE_G46T_nostoc_neighbourhood.tsv"),
 "F": ("A0ABU4T9S4", "Lentzea sp.", "Val218>Glu", "Protocatechuate catabolism", "neighbourhoods.tsv"),
 "G": ("A0A1G7K2Z2", "Cellulophaga sp.", "Asn205>Arg", "Bacillithiol biosynthesis", "groups/N205R/neighbourhoods.tsv"),
 "H": ("A0ABT1P7V4", "Streptantibioticus sp.", "Gly46>Asn", "Grixazone-type (3,4-AHBA)", "groups/G46N/neighbourhoods.tsv"),
}
SECTION = {"A":"NON-CANONICAL","B":"NON-CANONICAL","C":"CANONICAL BACTERIAL",
           "D":"CANONICAL BACTERIAL","E":"CANONICAL BACTERIAL","F":"CANONICAL BACTERIAL",
           "G":"CANONICAL BACTERIAL","H":"CANONICAL BACTERIAL"}

CAT_RULES = [
 ("ppo",         ["tyrosinase","di-copper","dicopper","copper-binding"]),
 ("regulatory",  ["transcriptional regulator","transcription factor","lysr","gntr","tetr","marr","arac","regulator","sigma factor","two-component"]),
 ("transport",   ["transporter","mfs","major facilitator","abc transporter","permease","tonb","outer membrane receptor","solute-binding","substrate-binding","symporter","efflux","carrier"]),
 ("pathway",     ["ustya","melc1","tyrosinase cofactor","tyrosinase co-factor","copper chaperone","copz","bsha","malate synthase","neud","o-acyltransferase","dltb","beta-glucosidase","chaplin","exopolysaccharide","glycosylphosphot","o-antigen","anthranilate","dehydroquinate","dahp","deoxy-7-phospho","tryptophan synthase","chorismate","carboxymuconolactone","protocatechuate","rio","carboxypeptidase","phytanoyl","cutinase","chitodextrinase","ammonia-lyase","decarboxylase","trpa","trpb","trpc","trpd","aro","tyra","prephenate","indole-3-glycerol"]),
 ("biosynthetic",["methyltransferase","p450","cytochrome","oxidoreductase","oxidase","dehydrogenase","hydrolase","glycosyltransferase","glycosyl transferase","synthase","synthetase","monooxygenase","dioxygenase","reductase","transferase","kinase","hydroxylase","deacetylase","amidase","amidohydrolase","ligase","isomerase","mutase","peptidase","sirtuin","cyclase","fad-binding","gnat","acetyltransferase","chitin","glycoside hydrolase","dehydratase","phosphatase","esterase","lyase","oxygenase","duf7730","dur7730","duf","dodecin","fmn-dependent"]),
]
def classify(product, is_target):
    if is_target == "1": return "ppo"
    p = (product or "").lower()
    for cat, kws in CAT_RULES:
        if any(k in p for k in kws): return cat
    return "hypothetical"

LABELS = [
 ("trpa,","TrpA"),("tryptophan synthase alpha","TrpA"),("trpb,","TrpB"),("tryptophan synthase beta","TrpB"),
 ("trpc,","TrpC"),("indole-3-glycerol phosphate synthase","TrpC"),("trpd,","TrpD"),("anthranilate phosphoribosyltransferase","TrpD"),
 ("trpeg","TrpE"),("anthranilate synthase","TrpE"),("arob,","aroB"),("3-dehydroquinate synthase","DHQ syn."),
 ("aroa2","DAHP syn."),("3-deoxy-7-phosphoheptulonate synthase","DAHP syn."),("tyra,","tyrA"),("prephenate dehydrogenase","preph. DH"),
 ("chaplin","chaplin"),("exopolysaccharide biosynthesis polyprenyl","EPS glycPT"),("polysaccharide deacetylase","PS deacet."),
 ("lipopolysaccharide biosynthesis","LPS biosynth."),("glycosyl transferase family 1","glycosTF"),("glycosyl transferase","glycosTF"),("glycosyltransferase","glycosTF"),
 ("o-antigen","O-antigen"),("melc1","MelC1"),("tyrosinase cofactor","MelC1"),("tyrosinase co-factor","MelC1"),
 ("lysr family transcriptional regulator","LysR"),("carboxymuconolactone decarboxylase","4-CML dec."),("gntr","GntR"),("protocatechuate","protocat."),
 ("multicopper oxidase","MCO"),("marr family","MarR"),
 ("n-acetyl-alpha-d-glucosaminyl l-malate synthase bsha","BshA"),("bsha","BshA"),("sugar o-acyltransferase","NeuD"),("neud","NeuD"),
 ("d-alanyl-lipoteichoic acid acyltransferase dltb","DltB"),("dltb","DltB"),("copper chaperone copz","CopZ"),("copz","CopZ"),
 ("outer membrane receptor for ferrienterochelin","TonB rec."),("tonb","TonB rec."),("beta-glucosidase","beta-gluc."),("chitodextrinase","chitodextrin."),
 ("sugar transferase involved in lps","sugarTF"),("mycotoxin biosynthesis protein ustya","UstYa"),("ustya","UstYa"),
 ("s-adenosyl-l-methionine-dependent methyltransfera","SAM MTase"),("cytochrome p450","P450"),("sirtuin","sirtuin"),
 ("chitin-binding","chitinase"),("nodb","chitinase"),("rio1 family","RIO kin."),("rio kinase","RIO kin."),
 ("zinc carboxypeptidase","Zn CPase"),("phytanoyl-coa dioxygenase","PhyH"),("cutinase","cutinase"),("glycoside hydrolase family 2","GH2"),("ergot alkaloid","ergot alk."),
 ("major facilitator superfamily transporter","MFS"),("mfs transporter","MFS"),("major facilitator","MFS"),
]
# Panel H (Gly46>Asn, grixazone) was annotated by hand after the generic rules
# ran, so that the 3,4-AHBA pathway genes the panel exists to show are named.
# That curation is recorded here rather than left implicit, so the spec
# regenerates exactly. Keyed on product string; never applied to the target
# gene, and scoped per panel because the same product is annotated
# differently elsewhere (panel D calls "tyrosinase cofactor" ppo, H calls it
# pathway).
PANEL_OVERRIDES = {
 "H": {
  'methyltransferase domain-containing protein': ('MTase', 'biosynthetic'),
  'cache and HAMP domain-containing protein': ('', 'hypothetical'),
  'poly-gamma-glutamate synthase PgsB': ('', 'pathway'),
  'poly-gamma-glutamate biosynthesis protein PgsC/CapC': ('', 'pathway'),
  'C40 family peptidase': ('', 'biosynthetic'),
  'glycosyltransferase family 39 protein': ('', 'biosynthetic'),
  'hydrolase': ('', 'biosynthetic'),
  '3-dehydroquinate synthase II': ('DHQ-II', 'pathway'),
  '2-amino-3,7-dideoxy-D-threo-hept-6-ulosonate synthase': ('aminoDAHP', 'pathway'),
  'aspartate kinase': ('Asp kin.', 'pathway'),
  'tyrosinase cofactor': ('MelC1', 'pathway'),
  'copper resistance protein CopC': ('CopC', 'transport'),
  'CopD family protein': ('', 'transport'),
  'LysR substrate-binding domain-containing protein': ('', 'regulatory'),
  'FAD/NAD(P)-binding protein': ('', 'biosynthetic'),
  'AfsR/SARP family transcriptional regulator': ('SARP reg.', 'regulatory'),
  '2-oxo acid dehydrogenase subunit E2': ('', 'biosynthetic'),
  'alpha-ketoacid dehydrogenase subunit beta': ('', 'biosynthetic'),
  'pyruvate dehydrogenase (acetyl-transferring) E1 component subunit alpha': ('', 'biosynthetic'),
  'hypothetical protein': ('', 'hypothetical'),
 },
}

def label(product, is_target):
    if is_target == "1": return "TDF"   # tyrosinase-domain family member
    p = (product or "").lower()
    for k, lab in LABELS:
        if k in p: return lab
    return ""

def load(path):
    d = defaultdict(list)
    with open(os.path.join(WORK, path)) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d[r["query_accession"]].append(r)
    return d

def main():
    header = ["panel","section","species","substitution","locus","accession",
              "idx","offset","label","category","strand","labelpos","product"]
    rows = []
    for panel in "ABCDEFGH":
        acc, sp, sub, loc, srcfile = PANELS[panel]
        src = load(srcfile)
        genes = sorted(src[acc], key=lambda g: int(g["offset"]))
        lc = 0
        overrides = PANEL_OVERRIDES.get(panel, {})
        for idx, g in enumerate(genes):
            cat = classify(g["product"], g["is_target"])
            lab = label(g["product"], g["is_target"])
            if lab and cat == "hypothetical": cat = "pathway"
            if g["is_target"] != "1" and g["product"] in overrides:
                lab, cat = overrides[g["product"]]
            if lab:
                labelpos = "above" if lc % 2 == 0 else "below"; lc += 1
            else:
                labelpos = "none"
            rows.append([panel, SECTION[panel], sp, sub, loc, acc, idx, g["offset"],
                         lab, cat, g["strand"], (g["product"] or "")[:90]])
            rows[-1].insert(11, labelpos)  # keep column order
    out = os.path.join(WORK, "gna_gene_maps_spec.tsv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(header); w.writerows(rows)
    print("wrote", out, "-", len(rows), "gene rows")

if __name__ == "__main__":
    main()
