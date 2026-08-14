#!/usr/bin/env python3
"""Stand-alone re-run of novelty_pipeline stage K over the ENTIRE canonical pool
(not the n=2000 sample). Reuses the validated stage-K machinery verbatim; only the
sample size changes. Usage: plddt_fullset.py <mounted_cif_dir>"""
import os, sys
import novelty_pipeline as N

cifs = sys.argv[1]
*_, V, _, _, charv = N.load()
print(f"loaded {len(V)} structures")
seqid = {}
for fn in os.listdir(cifs):
    if fn.endswith('_model.cif'):
        stem = fn[:-4]
        seqid[stem.split('_taxID_')[0]] = stem
N.stage_K(V, cifs, seqid, n_sample=10**9)  # min(n_sample, pool) -> whole pool
