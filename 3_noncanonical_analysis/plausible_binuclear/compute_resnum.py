#!/usr/bin/env python3
"""Recover each substituted His position's query residue number for the 66
plausible-binuclear structures, by replaying stage-3's helix-anchored superposition
(helix_seed -> ICP -> nearest query CA to the B2ZB02 anchor). Reuses the exact
functions/foldseek hits stage3_extract_align.py used, so the residue found here is
the one classify_sites.py scored. Outputs query auth_seq_id (PyMOL `resi`)."""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import stage3_extract_align as s3

CIF_DIR = os.path.join(HERE, "cifs_raw")
PB = os.path.join(HERE, "pb66.tsv")
FS = os.path.join(HERE, "..", "fs_vs_b2zb02.tsv")
PM = os.path.join(CIF_DIR, "B2ZB02_ref.cif")
CA_CUTOFF = 3.0

ref = s3.parse_atoms(PM)
ref_cus, ref_ca, cuA_pos, cuB_pos = s3.build_ref_anchor_sites(ref)
ref_seqs = s3.ordered_seqs(ref)
rcs = s3.ref_core_seqs(ref_ca)
LAB2ANCHOR = {"CuA_His1": cuA_pos[0], "CuA_His2": cuA_pos[1], "CuA_His3": cuA_pos[2],
              "CuB_His1": cuB_pos[0], "CuB_His2": cuB_pos[1], "CuB_His3": cuB_pos[2]}
# sequential index 1..6 by ascending residue number across all 6 anchors
SEQIDX = {p: i + 1 for i, p in enumerate(sorted(s3.ANCHORS))}
print(f"anchors CuA={cuA_pos} CuB={cuB_pos}  seqidx={SEQIDX}", flush=True)


def label_auth_map(path):
    """label_seq_id -> auth_seq_id for CA atoms (PyMOL resi uses auth_seq_id)."""
    cols = []; coll = False; inA = False; m = {}
    with open(path) as f:
        for raw in f:
            l = raw.strip()
            if l == "loop_": coll = False; inA = False; cols = []; continue
            if l.startswith("_atom_site."): coll = True; cols.append(l); continue
            if coll and cols: coll = False; inA = True
            if inA:
                if l.startswith("_") or l == "#" or not l: break
                p = l.split()
                if len(p) != len(cols): continue
                r = dict(zip(cols, p))
                if r["_atom_site.label_atom_id"].upper() == "CA":
                    ls = r["_atom_site.label_seq_id"]
                    if ls.isdigit(): m[int(ls)] = r["_atom_site.auth_seq_id"]
    return m


# foldseek best hit per accession (by qtm) -- same selection as stage3
fs = {}
with open(FS) as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 8: continue
        acc = p[0].split("_taxID_")[0]
        qstart = int(p[2]); tstart = int(p[5]); cigar = p[-1]
        qtm = float(p[7]) if p[7] else 0.0
        if acc not in fs or qtm > fs[acc][3]:
            fs[acc] = (qstart, tstart, cigar, qtm)

acc2cif = {}
for fn in os.listdir(CIF_DIR):
    if fn.endswith(".cif") and fn != "B2ZB02_ref.cif":
        acc2cif[fn.split("_taxID_")[0]] = fn


def nearest_seq(anchor_seq, qca_t, cutoff):
    rc = ref_ca[anchor_seq]; bd = 999.0; bs = None
    for s, c in qca_t.items():
        dd = float(np.linalg.norm(c - rc))
        if dd < bd: bd = dd; bs = s
    return (bs, bd) if (bs is not None and bd <= cutoff) else (None, bd)


out = []
for line in open(PB):
    acc, lab, newres, cadist = line.rstrip("\n").split("\t")
    fn = acc2cif[acc]
    qa = s3.parse_atoms(os.path.join(CIF_DIR, fn))
    qca = s3.ca_map(qa); qres3 = s3.resn3_map(qa); qseqs = s3.ordered_seqs(qa)
    la = label_auth_map(os.path.join(CIF_DIR, fn))
    qstart, tstart, cigar, qtm = fs[acc]
    seed = s3.helix_seed(qca, qseqs, qstart, tstart, cigar, ref_ca, ref_seqs)
    R, t, _ = seed
    R, t, _, _ = s3.icp_refine(qca, R, t, ref_ca, rcs)
    qca_t = {s: R @ c + t for s, c in qca.items()}
    anchor = LAB2ANCHOR[lab]
    bs, bd = nearest_seq(anchor, qca_t, CA_CUTOFF)
    rtype = qres3.get(bs, "UNK")
    auth = la.get(bs, str(bs))
    si = SEQIDX[anchor]
    ok = (rtype == newres)
    out.append((acc, lab, si, newres, rtype, bs, auth, f"{bd:.2f}", ok))
    print(f"{acc:14} {lab:9} seq{si} anchor{anchor} -> label_seq={bs} auth={auth} "
          f"{rtype} (want {newres}) d={bd:.2f} {'OK' if ok else 'MISMATCH'}", flush=True)

nbad = sum(1 for r in out if not r[8])
print(f"\n{len(out)} structures, {nbad} residue-type mismatches", flush=True)
with open(os.path.join(HERE, "pb66_resnum.tsv"), "w") as f:
    f.write("accession\tlabel\tseqidx\tnew_res\tres_type\tlabel_seq\tauth_seq\tca_dist\tmatch\n")
    for r in out:
        f.write("\t".join(str(x) for x in r) + "\n")
print("wrote pb66_resnum.tsv")
