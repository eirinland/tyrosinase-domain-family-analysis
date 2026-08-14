#!/usr/bin/env python3
"""Score every benchmark structure with six orthogonal core-detection methods.

Per structure the question is: does it carry the complete PPO 4-helix di-copper
core (-> belongs in canonical or non-canonical), or not (-> discard)?  Geometry
(Cu-Cu, His, pLDDT) only separates canonical from non-canonical and is recorded
as context; it does NOT decide core presence.

  M1  Foldseek global TMalign + query-backbone helix check (current pipeline).
  M2  Foldseek reference-normalized TM (ttmscore) -- kills the query-length
      normalization artifact that inflates short half-bundle fragments.
  M3  PyMOL super vs the 9 refs: best CA coverage + RMSD.
  M4  PyMOL cealign vs the 9 refs: best CA coverage + RMSD (CE is robust to
      large indels / distant homologs).
  M5  Biotite reference-free SSE: count long helices (>=8 res) and test whether
      the four longest pack into a compact bundle.
  M6  Combined: intrinsic bundle (M5) AND a PPO-identity gate (M1 or M2).
"""
import argparse, csv, math, os, re, sys

# ---- reference core-helix ranges (PDB author numbering) -- same as ppo_core_check.py
CORE_HELICES_PDB = {
    "ref_PmTYR":               [(34, 46), (65, 83), (203, 211), (226, 244)],
    "ref_8BBR_Vspinosum":      [(73, 83), (91, 110), (259, 265), (279, 300)],
    "ref_2Y9W_Abisporus":      [(54, 61), (90, 113), (255, 267), (291, 309)],
    "ref_5CE9_Jregia":         [(80, 91), (113, 131), (241, 246), (269, 283)],
    "ref_1BT3_Ibatatas":       [(81, 92), (114, 133), (240, 247), (269, 287)],
    "ref_5M8L_human":          [(184, 196), (220, 239), (376, 383), (399, 417)],
    "ref_1JS8_squid":          [(2536, 2543), (2567, 2584), (2660, 2679), (2697, 2719)],
    "ref_I3D139_archaea":      [(36, 47), (67, 85), (197, 205), (222, 240)],
    "ref_A0A9N8ELP9_oomycota": [(413, 424), (437, 456), (590, 601), (624, 642)],
}
D3_LO, D3_HI = 4.8, 6.4
D4_LO, D4_HI = 5.4, 7.4


def dist(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


def read_ca(path, chain="A"):
    """CA coords (ordered) and author resnums for chain A from a .cif or .pdb."""
    coords, resn = [], []
    if path.endswith(".pdb"):
        for line in open(path):
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            if line[21] not in (" ", "A"):
                continue
            try:
                resn.append(int(line[22:26]))
                coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except ValueError:
                continue
    else:  # AF3 CIF, whitespace-split
        for line in open(path):
            if not line.startswith("ATOM"):
                continue
            p = line.split()
            if len(p) < 15 or p[3] != "CA" or p[6] != chain:
                continue
            try:
                resn.append(int(p[8]))
                coords.append((float(p[10]), float(p[11]), float(p[12])))
            except ValueError:
                continue
    return coords, resn


def read_ca_pdb_resnums(path):
    out = []
    for line in open(path):
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        if line[21] not in (" ", "A"):
            continue
        try:
            out.append(int(line[22:26]))
        except ValueError:
            continue
    return out


def ref_helix_seq_ranges(pdb, helices):
    resn = read_ca_pdb_resnums(pdb)
    ranges = []
    for a, b in helices:
        idx = [i for i, rn in enumerate(resn, 1) if a <= rn <= b]
        ranges.append((min(idx), max(idx)) if idx else None)
    return ranges, len(resn)


def helical_mask(coords):
    n = len(coords)
    inhel = [False]*n
    for i in range(n-4):
        if D3_LO < dist(coords[i], coords[i+3]) < D3_HI and D4_LO < dist(coords[i], coords[i+4]) < D4_HI:
            for j in range(i, min(i+5, n)):
                inhel[j] = True
    return inhel


def helix_present(coords, seq_idx, min_h=4, min_frac=0.5):
    if not seq_idx:
        return False
    lo, hi = min(seq_idx), max(seq_idx)
    span = list(range(lo, hi+1))
    if len(span) < min_h:
        return False
    inhel = helical_mask(coords)
    cnt = sum(1 for s in span if 1 <= s <= len(coords) and inhel[s-1])
    return cnt >= min_h and cnt/len(span) >= min_frac


def parse_cigar(qstart, tstart, cigar):
    qp, tp, num = qstart, tstart, ""
    for ch in cigar:
        if ch.isdigit():
            num += ch; continue
        n = int(num) if num else 1; num = ""
        if ch == "M":
            for _ in range(n):
                yield qp, tp; qp += 1; tp += 1
        elif ch == "I":
            qp += n
        elif ch == "D":
            tp += n


def load_foldseek(tsv):
    hits = {}
    for line in open(tsv):
        f = line.rstrip("\n").split("\t")
        if len(f) < 13:
            continue
        try:
            d = dict(target=f[1], qstart=int(f[2]), tstart=int(f[5]),
                     qtm=float(f[7]), ttm=float(f[8]), alntm=float(f[9]),
                     lddt=float(f[10]), cigar=f[12])
        except ValueError:
            continue
        hits.setdefault(f[0].split("_taxID_")[0], []).append(d)
    return hits


def biotite_bundle(path):
    """Reference-free: count long helices (>=8 res) and find the largest set of
    mutually-packed helices (min inter-helix CA-CA < 10 A). A 4-helix bundle has
    a connected component of >=4 such helices. Returns (n_long, bundle, comp)."""
    import numpy as np
    from scipy.spatial.distance import cdist
    import biotite.structure as struc
    if path.endswith(".pdb"):
        import biotite.structure.io.pdb as io
        arr = io.PDBFile.read(path).get_structure(model=1)
    else:
        import biotite.structure.io.pdbx as pdbx
        arr = pdbx.get_structure(pdbx.CIFFile.read(path), model=1)
    arr = arr[arr.chain_id == arr.chain_id[0]]
    try:
        sse = "".join(struc.annotate_sse(arr))
    except Exception:
        return 0, False, ""
    ca = arr[arr.atom_name == "CA"]
    segs = [(m.start(), m.end()) for m in re.finditer("a+", sse)
            if m.end()-m.start() >= 8 and m.end() <= len(ca)]
    n_long = len(segs)
    if n_long < 4:
        return n_long, False, str(n_long)
    hc = [ca.coord[s:e] for s, e in segs]
    n = len(hc)
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if cdist(hc[i], hc[j]).min() < 10.0:
                adj[i].append(j); adj[j].append(i)
    seen, best = [False]*n, 0
    for s in range(n):
        if seen[s]:
            continue
        stack, comp = [s], 0
        seen[s] = True
        while stack:
            u = stack.pop(); comp += 1
            for v in adj[u]:
                if not seen[v]:
                    seen[v] = True; stack.append(v)
        best = max(best, comp)
    return n_long, best >= 4, str(best)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif-dir", required=True)
    ap.add_argument("--ref-dir", required=True)
    ap.add_argument("--fs-tsv", required=True)
    ap.add_argument("--bench", required=True)       # benchmark_set.csv
    ap.add_argument("--criteria", required=True)    # canonical_criteria_all_ca.csv
    ap.add_argument("--pool", required=True)        # three_pool_assignment.csv
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    bench = {r["accession"]: r for r in csv.DictReader(open(a.bench))}
    pool = {r["accession"]: r["pool"] for r in csv.DictReader(open(a.pool))}
    geom = {r["accession"]: r for r in csv.DictReader(open(a.criteria))}
    hits = load_foldseek(a.fs_tsv)

    refs, ref_nca = {}, {}
    for ref, hel in CORE_HELICES_PDB.items():
        p = os.path.join(a.ref_dir, ref + ".pdb")
        if os.path.exists(p):
            refs[ref], ref_nca[ref] = ref_helix_seq_ranges(p, hel)

    from pymol import cmd
    cmd.feedback("disable", "all", "everything")
    for ref in refs:
        cmd.load(os.path.join(a.ref_dir, ref + ".pdb"), ref)

    # collect query files (cif-dir holds the benchmark CIFs + the 9 ref PDBs)
    queries = []
    for fn in sorted(os.listdir(a.cif_dir)):
        if fn.endswith(".cif"):
            queries.append((fn.split("_taxID_")[0], os.path.join(a.cif_dir, fn)))
        elif fn.endswith(".pdb"):
            queries.append((fn[:-4], os.path.join(a.cif_dir, fn)))

    cols = ["accession", "stratum", "expected", "pool", "qlen", "cu_dist", "n_his",
            "min_plddt", "canonical",
            "m1_core_ok", "m1_nhel", "m1_best_ref", "m1_qtm",
            "m2_best_ttm", "m2_best_alntm", "m2_ttm_ref",
            "m3_super_cov", "m3_super_rmsd", "m3_super_naln", "m3_super_ref",
            "m4_ce_cov", "m4_ce_rmsd", "m4_ce_naln", "m4_ce_ref",
            "m5_nlonghelix", "m5_bundle", "m5_bundlespan", "m6_combined"]
    w = csv.writer(open(a.output, "w", newline="\n"), delimiter="\t")
    w.writerow(cols)

    for i, (acc, path) in enumerate(queries, 1):
        is_ref = acc in refs
        coords, _ = read_ca(path)
        qlen = len(coords)
        if qlen < 10:
            continue
        b = bench.get(acc, {})
        stratum = b.get("stratum", "ctrl_ref" if is_ref else "NA")
        expected = b.get("expected", "core_present" if is_ref else "NA")
        g = geom.get(acc, {})

        # --- M1 / M2 from foldseek
        m1_core_ok = m1_nhel = m1_ref = ""
        m1_qtm = m2_ttm = m2_alntm = m2_ttm_ref = ""
        hl = [h for h in hits.get(acc, []) if h["target"] in refs]
        if hl:
            best = max(hl, key=lambda h: h["qtm"])
            m1_ref, m1_qtm = best["target"], f"{best['qtm']:.3f}"
            pairs = list(parse_cigar(best["qstart"], best["tstart"], best["cigar"]))
            pres = []
            for hr in refs[best["target"]]:
                if hr is None:
                    pres.append(False); continue
                ts, te = hr
                qidx = [qp for qp, tp in pairs if ts <= tp <= te]
                pres.append(helix_present(coords, qidx))
            m1_nhel = sum(pres); m1_core_ok = (m1_nhel == 4)
            bt = max(hl, key=lambda h: h["ttm"])
            m2_ttm, m2_alntm, m2_ttm_ref = f"{bt['ttm']:.3f}", f"{bt['alntm']:.3f}", bt["target"]

        # --- M3 super / M4 cealign vs 9 refs
        cmd.load(path, "mob")
        best_s = best_c = None
        for ref in refs:
            try:
                rs = cmd.super("mob and name CA", ref + " and name CA")
                cov = rs[1] / ref_nca[ref]
                if best_s is None or cov > best_s[0]:
                    best_s = (cov, rs[0], rs[1], ref)
            except Exception:
                pass
            try:
                rc = cmd.cealign(ref + " and name CA", "mob and name CA")
                cov = rc["alignment_length"] / ref_nca[ref]
                if best_c is None or cov > best_c[0]:
                    best_c = (cov, rc["RMSD"], rc["alignment_length"], ref)
            except Exception:
                pass
        cmd.delete("mob")

        # --- M5 intrinsic bundle
        try:
            n_long, bundle, span = biotite_bundle(path)
        except Exception:
            n_long, bundle, span = "", "", ""

        # --- M6 combined
        m6 = bool(bundle) and (m1_core_ok is True or (m2_ttm != "" and float(m2_ttm) >= 0.5))

        w.writerow([
            acc, stratum, expected, pool.get(acc, "ref" if is_ref else "NA"), qlen,
            g.get("cu_dist", ""), g.get("n_coord_his", ""), g.get("min_plddt", ""),
            g.get("canonical", ""),
            m1_core_ok, m1_nhel, m1_ref, m1_qtm,
            m2_ttm, m2_alntm, m2_ttm_ref,
            f"{best_s[0]:.3f}" if best_s else "", f"{best_s[1]:.2f}" if best_s else "",
            best_s[2] if best_s else "", best_s[3] if best_s else "",
            f"{best_c[0]:.3f}" if best_c else "", f"{best_c[1]:.2f}" if best_c else "",
            best_c[2] if best_c else "", best_c[3] if best_c else "",
            n_long, bundle, span, m6,
        ])
        if i % 25 == 0:
            print(f"  ...{i}/{len(queries)}", file=sys.stderr)
    print(f"scored {len(queries)} structures", file=sys.stderr)


if __name__ == "__main__":
    main()
