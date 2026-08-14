#!/usr/bin/env python3
"""Per-structure Val218 Cα pLDDT for carriers of selected substitutions.
Reuses novelty_pipeline stage-K Kabsch machinery (6-His superposition onto PmTYR,
nearest-Cα mapping at the Val218 position, read its B_iso = pLDDT).
Usage: val218_group_plddt.py <mounted_cif_dir> [residues]   (default residues: SYKEHR)
Writes val218_group_plddt.csv."""
import os, sys, csv
import numpy as np
from multiprocessing import Pool
import novelty_pipeline as N

TARGET = set(sys.argv[2]) if len(sys.argv) > 2 else set("SYKEHR")
VIDX = N.POS.index("Val218")  # 6


def worker(task):
    path, acc = task
    if not os.path.exists(path):
        return None
    at = N._parse_atoms(path)
    qca = N._superpose(at)
    if not qca:
        return None
    pl = {int(x['seq']): x['bfactor'] for x in at if x['atom'] == 'CA' and x['seq'].isdigit()}
    refc = N._MAP_REFPOS[VIDX]
    bs = min(qca, key=lambda s: np.linalg.norm(qca[s] - refc))
    return (acc, pl.get(bs))


def main():
    cifs = sys.argv[1]
    seqid = {}
    for fn in os.listdir(cifs):
        if fn.endswith('_model.cif'):
            seqid[fn.split('_taxID_')[0]] = fn[:-4]
    sel = {}
    for r in csv.DictReader(open(N.PVEC)):
        if r.get('error'):
            continue
        if r.get('Val218') in TARGET:
            sel[r['accession']] = (r['Val218'], r.get('Val218_cadist'))
    ref = N._ref_variable(cifs)
    if ref is None:
        sys.exit("reference CIF/positions unavailable")
    rne2, refpos = ref
    tasks = [(f"{cifs}/{seqid[a]}.cif", a) for a in sel if a in seqid]
    print(f"target residues {sorted(TARGET)}; selected {len(sel)} carriers; {len(tasks)} with CIF")
    nproc = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
    out = {}
    with Pool(nproc, initializer=N._map_init, initargs=(rne2, refpos)) as pool:
        for r in pool.imap_unordered(worker, tasks, chunksize=50):
            if r:
                out[r[0]] = r[1]
    tax = {}
    for r in csv.DictReader(open(N.TAXf)):
        tax[r['accession']] = (r.get('kingdom', ''), r.get('phylum', ''), r.get('genus', ''))
    with open("val218_group_plddt.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["accession", "val218_res", "val218_plddt", "val218_cadist", "kingdom", "phylum", "genus"])
        for a, (res, cad) in sel.items():
            if a not in out:
                continue
            k = tax.get(a, ("", "", ""))
            pl = out[a]
            w.writerow([a, res, f"{pl:.1f}" if pl is not None else "", cad, *k])
    print(f"wrote val218_group_plddt.csv ({len(out)} rows)")


if __name__ == "__main__":
    main()
