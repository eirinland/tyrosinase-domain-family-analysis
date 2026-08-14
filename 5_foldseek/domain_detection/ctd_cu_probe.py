#!/usr/bin/env python3
"""Detect C-terminal blocking aromatic by proximity to active-site Cu.

For each structure: find Cu atoms, identify the last coordinating His in
sequence, then search for Phe/Trp/Tyr past that His whose sidechain is
within a distance cutoff of either Cu. Reports the closest such aromatic.
"""
import argparse
import csv
import os
import sys

import numpy as np
from Bio.PDB import MMCIFParser

AROMATICS = {"PHE", "TRP", "TYR"}
BACKBONE = {"N", "CA", "C", "O"}


def get_cu_positions(structure):
    cus = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.get_resname().strip() in ("CU", "CU1", "CU2"):
                    for atom in res:
                        cus.append(atom.get_vector().get_array())
    return [np.array(c) for c in cus]


def get_residues(structure):
    residues = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] != " ":
                    continue
                residues.append(res)
        break
    return residues


def sidechain_center(res):
    coords = []
    for atom in res:
        if atom.get_name() not in BACKBONE:
            coords.append(atom.get_vector().get_array())
    if not coords:
        return None
    return np.mean(coords, axis=0)


def check_helix_context(residues, blocker_seqid, window=4):
    """Check if blocker sits on a helix using CA i->i+3 and i->i+4 distances."""
    ca_map = {}
    for res in residues:
        if "CA" in res:
            ca_map[res.id[1]] = res["CA"].get_vector().get_array()

    sid = blocker_seqid
    n_helical = 0
    for offset in range(-window, window + 1):
        r = sid + offset
        if r not in ca_map:
            continue
        d3_ok = False
        d4_ok = False
        if r + 3 in ca_map:
            d3 = np.linalg.norm(np.array(ca_map[r]) - np.array(ca_map[r + 3]))
            d3_ok = 4.0 <= d3 <= 6.4
        if r + 4 in ca_map:
            d4 = np.linalg.norm(np.array(ca_map[r]) - np.array(ca_map[r + 4]))
            d4_ok = 4.8 <= d4 <= 8.2
        if d3_ok and d4_ok:
            n_helical += 1

    return n_helical


def min_sidechain_cu_dist(res, cus):
    best = float("inf")
    for atom in res:
        if atom.get_name() in BACKBONE:
            continue
        coord = atom.get_vector().get_array()
        for cu in cus:
            d = np.linalg.norm(coord - cu)
            if d < best:
                best = d
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", required=True, help="Directory of CIF files")
    parser.add_argument("--output", required=True, help="Output TSV")
    parser.add_argument("--cu-dist", type=float, default=8.0,
                        help="Max distance from aromatic sidechain to Cu (default 8.0)")
    parser.add_argument("--his-cu-dist", type=float, default=3.5,
                        help="Max NE2-Cu distance for coordinating His (default 3.5)")
    args = parser.parse_args()

    cif_parser = MMCIFParser(QUIET=True)
    targets = sorted(f for f in os.listdir(args.targets) if f.endswith(".cif"))
    print("Targets: %d, Cu cutoff: %.1f A" % (len(targets), args.cu_dist))

    rows = []
    for ti, fname in enumerate(targets):
        name = fname.replace(".cif", "")
        path = os.path.join(args.targets, fname)
        try:
            struct = cif_parser.get_structure("t", path)
            cus = get_cu_positions(struct)
            residues = get_residues(struct)

            if len(cus) < 1:
                rows.append({"name": name, "n_cu": 0, "last_his_seqid": "",
                             "blocker_resname": "", "blocker_seqid": "",
                             "blocker_cu_dist": "", "n_helical": "",
                             "n_arom_past_his": 0,
                             "n_arom_near_cu": 0, "has_blocker": False,
                             "all_arom_near": ""})
                continue

            # Find coordinating His
            coord_his = []
            for res in residues:
                if res.get_resname() != "HIS":
                    continue
                ne2 = None
                for atom in res:
                    if atom.get_name() == "NE2":
                        ne2 = atom.get_vector().get_array()
                        break
                if ne2 is None:
                    continue
                for cu in cus:
                    if np.linalg.norm(ne2 - cu) <= args.his_cu_dist:
                        coord_his.append(res.id[1])
                        break

            last_his = max(coord_his) if coord_his else 0

            # Find aromatics past last coordinating His
            arom_past = []
            for res in residues:
                if res.get_resname() not in AROMATICS:
                    continue
                if res.id[1] <= last_his:
                    continue
                d = min_sidechain_cu_dist(res, cus)
                arom_past.append((res.id[1], res.get_resname(), d))

            near_cu = [(sid, rn, d) for sid, rn, d in arom_past if d <= args.cu_dist]
            near_cu.sort(key=lambda x: x[2])

            if near_cu:
                best = near_cu[0]
                all_str = ";".join("%s%d:%.1f" % (rn, sid, d) for sid, rn, d in near_cu)
                n_hel = check_helix_context(residues, best[0])
                rows.append({"name": name, "n_cu": len(cus),
                             "last_his_seqid": last_his,
                             "blocker_resname": best[1], "blocker_seqid": best[0],
                             "blocker_cu_dist": "%.2f" % best[2],
                             "n_helical": n_hel,
                             "n_arom_past_his": len(arom_past),
                             "n_arom_near_cu": len(near_cu),
                             "has_blocker": True,
                             "all_arom_near": all_str})
            else:
                rows.append({"name": name, "n_cu": len(cus),
                             "last_his_seqid": last_his,
                             "blocker_resname": "", "blocker_seqid": "",
                             "blocker_cu_dist": "", "n_helical": "",
                             "n_arom_past_his": len(arom_past),
                             "n_arom_near_cu": 0,
                             "has_blocker": False,
                             "all_arom_near": ""})

        except Exception as e:
            print("  ERROR %s: %s" % (name, e), file=sys.stderr, flush=True)
            rows.append({"name": name, "n_cu": -1, "last_his_seqid": "",
                         "blocker_resname": "", "blocker_seqid": "",
                         "blocker_cu_dist": "", "n_helical": "",
                         "n_arom_past_his": 0,
                         "n_arom_near_cu": 0, "has_blocker": False,
                         "all_arom_near": ""})

        if (ti + 1) % 100 == 0:
            print("  Processed %d/%d..." % (ti + 1, len(targets)), flush=True)

    fields = ["sequence_id", "n_cu", "last_his_seqid", "has_blocker",
              "blocker_resname", "blocker_seqid", "blocker_cu_dist",
              "n_helical", "n_arom_past_his", "n_arom_near_cu", "all_arom_near"]
    with open(args.output, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(fields)
        for r in sorted(rows, key=lambda x: x["name"]):
            w.writerow([r["name"], r["n_cu"], r["last_his_seqid"],
                        r["has_blocker"], r["blocker_resname"],
                        r["blocker_seqid"], r["blocker_cu_dist"],
                        r["n_helical"],
                        r["n_arom_past_his"], r["n_arom_near_cu"],
                        r["all_arom_near"]])

    n_blocker = sum(1 for r in rows if r["has_blocker"])
    print("\n%d/%d have blocking aromatic within %.1f A of Cu" % (
        n_blocker, len(rows), args.cu_dist))

    from collections import Counter
    gate = Counter(r["blocker_resname"] for r in rows if r["has_blocker"])
    for aa, n in gate.most_common():
        print("  %s: %d (%.1f%%)" % (aa, n, 100 * n / n_blocker))

    hels = sorted(int(r["n_helical"]) for r in rows if r["has_blocker"])
    if hels:
        print("Helical context (residues in helix within +/-4 of blocker):")
        for c in [0, 1, 2, 3, 4, 5]:
            n = sum(1 for h in hels if h >= c)
            print("  >= %d helical: %d/%d (%.1f%%)" % (c, n, len(hels), 100*n/len(hels)))

    dists = sorted(float(r["blocker_cu_dist"]) for r in rows if r["has_blocker"])
    if dists:
        print("Distance: median %.1f, min %.1f, max %.1f" % (
            dists[len(dists)//2], dists[0], dists[-1]))


if __name__ == "__main__":
    main()
