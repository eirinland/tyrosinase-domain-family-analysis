#!/usr/bin/env python3
"""Query-structure 4-helix core check (global TMalign superposition).

Replaces the broken alignment-coverage helix filter. For each candidate:
  1. Foldseek --alignment-type 1 (global TMalign) vs the 9 PPO references.
  2. Pick the best-qtmscore reference hit.
  3. Map each of that reference's 4 core-helix residue ranges to query
     residues via the global alignment cigar.
  4. Test the QUERY's OWN backbone secondary structure at each mapped
     window (Ca(i,i+3) 4.8-6.4 A AND Ca(i,i+4) 5.4-7.4 A).
  5. Require all 4 core helices present in the query.

The point: we test the query structure, not the alignment. Global TMalign
does not trim terminal helices the way local alignment does.
"""
import argparse, csv, math, os, sys

# Reference core-helix ranges in each reference's PDB author residue numbering.
CORE_HELICES_PDB = {
    "ref_PmTYR":            [(34, 46),   (65, 83),   (203, 211), (226, 244)],
    "ref_8BBR_Vspinosum":   [(73, 83),   (91, 110),  (259, 265), (279, 300)],
    "ref_2Y9W_Abisporus":   [(54, 61),   (90, 113),  (255, 267), (291, 309)],
    "ref_5CE9_Jregia":      [(80, 91),   (113, 131), (241, 246), (269, 283)],
    "ref_1BT3_Ibatatas":    [(81, 92),   (114, 133), (240, 247), (269, 287)],
    "ref_5M8L_human":       [(184, 196), (220, 239), (376, 383), (399, 417)],
    "ref_1JS8_squid":       [(2536, 2543), (2567, 2584), (2660, 2679), (2697, 2719)],
    "ref_I3D139_archaea":   [(36, 47),   (67, 85),   (197, 205), (222, 240)],
    "ref_A0A9N8ELP9_oomycota": [(413, 424), (437, 456), (590, 601), (624, 642)],
}

# Helix geometry windows (A) on consecutive Ca atoms.
D3_LO, D3_HI = 4.8, 6.4
D4_LO, D4_HI = 5.4, 7.4


def dist(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


def read_ca_cif(path, chain="A"):
    """Ordered list of (resnum, (x,y,z)) for CA atoms of given chain."""
    out = []
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        p = line.split()
        if len(p) < 15:
            continue
        if p[3] != "CA" or p[6] != chain:
            continue
        try:
            out.append((int(p[8]), (float(p[10]), float(p[11]), float(p[12]))))
        except ValueError:
            continue
    return out


def read_ca_pdb(path):
    """Ordered list of PDB author resnums for CA atoms (chain A or blank)."""
    resnums = []
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        if line[12:16].strip() != "CA":
            continue
        ch = line[21]
        if ch not in (" ", "A"):
            continue
        try:
            resnums.append(int(line[22:26]))
        except ValueError:
            continue
    return resnums


def ref_helix_seq_ranges(ref_pdb_path, helices_pdb):
    """Map PDB helix ranges -> 1-based sequential CA-index ranges (Foldseek coords)."""
    resnums = read_ca_pdb(ref_pdb_path)
    idx_of = {}  # pdb resnum -> first seq index (1-based)
    for i, rn in enumerate(resnums, 1):
        idx_of.setdefault(rn, i)
    ranges = []
    for a, b in helices_pdb:
        seq_idx = [i for i, rn in enumerate(resnums, 1) if a <= rn <= b]
        if seq_idx:
            ranges.append((min(seq_idx), max(seq_idx)))
        else:
            ranges.append(None)
    return ranges


def helical_mask(coords):
    """Boolean list: is residue i the start of a helical i..i+4 stretch."""
    n = len(coords)
    mask = [False] * n
    for i in range(n - 4):
        d3 = dist(coords[i], coords[i+3])
        d4 = dist(coords[i], coords[i+4])
        if D3_LO < d3 < D3_HI and D4_LO < d4 < D4_HI:
            mask[i] = True
    return mask


def query_helix_present(q_coords, q_seq_indices, min_helical=4, min_frac=0.5):
    """Given query CA coords (ordered) and the set of 1-based seq indices that
    aligned to a reference helix window, decide whether the query forms a helix."""
    if not q_seq_indices:
        return False, 0, 0
    lo, hi = min(q_seq_indices), max(q_seq_indices)
    span = list(range(lo, hi + 1))  # 1-based inclusive query span
    if len(span) < min_helical:
        return False, len(span), 0
    mask = helical_mask(q_coords)  # 0-based start-of-helix flags
    # A residue is "in a helix" if it participates in any helical i..i+4 window.
    inhelix = [False] * len(q_coords)
    for i, h in enumerate(mask):
        if h:
            for j in range(i, min(i + 5, len(q_coords))):
                inhelix[j] = True
    helical_count = sum(1 for s in span if 1 <= s <= len(q_coords) and inhelix[s-1])
    frac = helical_count / len(span)
    present = helical_count >= min_helical and frac >= min_frac
    return present, len(span), helical_count


def parse_cigar_map(qstart, tstart, cigar):
    """Yield (qpos, tpos) 1-based sequential index pairs for aligned (M) columns."""
    qpos, tpos = qstart, tstart
    num = ""
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        n = int(num) if num else 1
        num = ""
        if ch == "M":
            for _ in range(n):
                yield qpos, tpos
                qpos += 1
                tpos += 1
        elif ch in ("I",):       # insertion in query (gap in target)
            qpos += n
        elif ch in ("D",):       # deletion in query (gap in query) -> advance target
            tpos += n
    return


def load_foldseek(tsv):
    """query -> list of hit dicts. Foldseek query name '<acc>_taxID_<id>_model_A'."""
    hits = {}
    for line in open(tsv):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        q, t = f[0], f[1]
        try:
            d = dict(query=q, target=t,
                     qstart=int(f[2]), qend=int(f[3]), qlen=int(f[4]),
                     tstart=int(f[5]), tend=int(f[6]),
                     qtmscore=float(f[7]), cigar=f[8])
        except ValueError:
            continue
        hits.setdefault(q, []).append(d)
    return hits


def acc_from_query(q):
    # '<acc>_taxID_<id>_model_A' -> acc
    return q.split("_taxID_")[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif-dir", required=True)
    ap.add_argument("--ref-dir", required=True, help="dir with ref_*.pdb")
    ap.add_argument("--fs-tsv", required=True, help="foldseek --alignment-type 1 output")
    ap.add_argument("--output", required=True)
    ap.add_argument("--min-helical", type=int, default=4)
    ap.add_argument("--min-frac", type=float, default=0.5)
    args = ap.parse_args()

    # Precompute reference helix sequential ranges.
    ref_ranges = {}
    for ref, helices in CORE_HELICES_PDB.items():
        pdb = os.path.join(args.ref_dir, ref + ".pdb")
        if not os.path.exists(pdb):
            print(f"WARN missing ref {pdb}", file=sys.stderr)
            continue
        ref_ranges[ref] = ref_helix_seq_ranges(pdb, helices)

    hits = load_foldseek(args.fs_tsv)
    print(f"queries with hits: {len(hits)}", file=sys.stderr)

    cif_for = {}
    for fn in os.listdir(args.cif_dir):
        if fn.endswith(".cif"):
            cif_for[acc_from_query(fn)] = os.path.join(args.cif_dir, fn)

    w = csv.writer(open(args.output, "w"))
    w.writerow(["accession", "best_ref", "best_qtm",
                "a1", "a2", "a3", "a4", "n_helices", "core_ok"])
    n_ok = n_tot = 0
    for q, hlist in hits.items():
        acc = acc_from_query(q)
        cif = cif_for.get(acc)
        if not cif:
            continue
        q_coords = [c for _, c in read_ca_cif(cif)]
        if len(q_coords) < 10:
            continue
        # best-qtm ref that we have helix ranges for
        cand = [h for h in hlist if h["target"] in ref_ranges]
        if not cand:
            continue
        best = max(cand, key=lambda h: h["qtmscore"])
        ref = best["target"]
        # build qpos->set, per helix window collect query seq indices
        pairs = list(parse_cigar_map(best["qstart"], best["tstart"], best["cigar"]))
        present = []
        for hr in ref_ranges[ref]:
            if hr is None:
                present.append(False)
                continue
            ts, te = hr
            qidx = [qp for qp, tp in pairs if ts <= tp <= te]
            ok, _, _ = query_helix_present(q_coords, qidx,
                                           args.min_helical, args.min_frac)
            present.append(ok)
        nhel = sum(present)
        core_ok = nhel == 4
        n_tot += 1
        n_ok += core_ok
        w.writerow([acc, ref, f"{best['qtmscore']:.4f}",
                    *[int(x) for x in present], nhel, core_ok])
    print(f"scored {n_tot}, core_ok {n_ok}", file=sys.stderr)


if __name__ == "__main__":
    main()
