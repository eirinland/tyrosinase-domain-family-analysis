#!/usr/bin/env python3
"""Stream the AF3 model archive once and cache Ca coordinates + pLDDT for the
structures this check needs (1,060 non-canonical + 241 benchmark), so the analysis
itself is re-runnable in seconds. Output: ca_cache.npz
"""
import os
import tarfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
# pipeline root, two levels up: <root>/1_filtering/separate_helix_check
BASE = os.path.dirname(os.path.dirname(HERE))
# The AF3 model archive is a Zenodo item, not part of the repository. Only cache_ca.py
# needs it, and only when ca_cache.npz is missing.
TARBALL = os.environ.get("PPO_AF3_TARBALL", os.path.join(
    os.path.dirname(os.path.dirname(BASE)), "ppo_af3_models.tar.gz"))


def parse_cif_ca(fh):
    cols, inA, coll = [], False, False
    seq, xyz, bf = [], [], []
    for raw in fh:
        l = raw.decode("utf-8", "replace").strip() if isinstance(raw, bytes) else raw.strip()
        if l == "loop_":
            coll, inA, cols = False, False, []
            continue
        if l.startswith("_atom_site."):
            coll = True
            cols.append(l)
            continue
        if coll and cols:
            coll, inA = False, True
        if inA:
            if l.startswith("_") or l == "#" or not l:
                break
            p = l.split()
            if len(p) != len(cols):
                continue
            r = dict(zip(cols, p))
            try:
                if r["_atom_site.label_atom_id"].upper() != "CA":
                    continue
                s = r["_atom_site.label_seq_id"]
                if not s.lstrip("-").isdigit():
                    continue
                seq.append(int(s))
                xyz.append([float(r["_atom_site.Cartn_x"]), float(r["_atom_site.Cartn_y"]),
                            float(r["_atom_site.Cartn_z"])])
                bf.append(float(r.get("_atom_site.B_iso_or_equiv", "0") or 0))
            except (KeyError, ValueError):
                continue
    return np.array(seq, dtype=np.int32), np.array(xyz, dtype=np.float32), \
        np.array(bf, dtype=np.float32)


pools = pd.read_csv(f"{BASE}/1_filtering/final_pools/three_pool_assignment_final.csv")
bench = pd.read_csv(f"{BASE}/1_filtering/benchmark/benchmark_results.tsv", sep="\t")
want = set(pools.loc[pools.pool == "noncanonical", "accession"]) | set(bench.accession)
print(f"caching {len(want):,} structures from {os.path.basename(TARBALL)} ...")

out = {}
with tarfile.open(TARBALL, "r|gz") as tf:
    for m in tf:
        if not m.name.endswith("_model.cif"):
            continue
        acc = os.path.basename(m.name).split("_taxID_")[0]
        if acc not in want or f"{acc}|seq" in out:
            continue
        s, x, b = parse_cif_ca(tf.extractfile(m))
        out[f"{acc}|seq"] = s
        out[f"{acc}|xyz"] = x
        out[f"{acc}|bf"] = b
        if len(out) // 3 % 250 == 0:
            print(f"  {len(out)//3:,}/{len(want):,}")
        if len(out) // 3 == len(want):
            break
np.savez_compressed(f"{HERE}/ca_cache.npz", **out)
print(f"cached {len(out)//3:,}/{len(want):,} -> ca_cache.npz")
