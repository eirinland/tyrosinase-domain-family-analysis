#!/usr/bin/env python3
"""Calibrate the Ca-only 'which helix is this residue on' definition against the
five copper-bearing references, whose core-helix ranges are known.

A separate-helix rule is only meaningful if, on a textbook PPO, the four core-helix
anchors fall in four different segments. Any definition that merges a3 and a4 on
PmTYR is wrong, not informative. Prints the segments each candidate definition
gives around the known core-helix ranges, and whether the four anchors separate.
"""
import os
import platform

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(HERE)
FS = f"{BASE}/foldseek"

CORE = {
    "ref_PmTYR": [(34, 46), (65, 83), (203, 211), (226, 244)],
    "ref_2Y9W_Abisporus": [(54, 61), (90, 113), (255, 267), (291, 309)],
    "ref_5CE9_Jregia": [(80, 91), (113, 131), (241, 246), (269, 283)],
    "ref_1BT3_Ibatatas": [(81, 92), (114, 133), (240, 247), (269, 287)],
    "ref_1JS8_squid": [(2536, 2543), (2567, 2584), (2660, 2679), (2697, 2719)],
}
ANCHORS = {
    "ref_PmTYR": [42, 69, 204, 231],
    "ref_2Y9W_Abisporus": [61, 94, 259, 296],
    "ref_5CE9_Jregia": [87, 117, 243, 273],
    "ref_1BT3_Ibatatas": [88, 118, 240, 274],
    "ref_1JS8_squid": [2543, 2571, 2671, 2702],
}

DEFS = [
    ("A widened turns>=1", (4.0, 6.4), (4.8, 8.2), 1),
    ("B widened turns>=3", (4.0, 6.4), (4.8, 8.2), 3),
    ("C strict  turns>=1", (4.8, 6.4), (5.4, 7.4), 1),
    ("D strict  turns>=3", (4.8, 6.4), (5.4, 7.4), 3),
    ("E strict  turns>=5", (4.8, 6.4), (5.4, 7.4), 5),
]


def parse_pdb_ca(path):
    ca = {}
    with open(path) as f:
        for line in f:
            if line[:6].strip() not in ("ATOM", "HETATM"):
                continue
            if line[12:16].strip().upper() != "CA":
                continue
            try:
                ca[int(line[22:26])] = np.array([float(line[30:38]), float(line[38:46]),
                                                 float(line[46:54])])
            except ValueError:
                continue
    return ca


def parse_cif_ca(path):
    cols, inA, coll = [], False, False
    ca = {}
    for raw in open(path):
        l = raw.strip()
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
                ca[int(r["_atom_site.label_seq_id"])] = np.array(
                    [float(r["_atom_site.Cartn_x"]), float(r["_atom_site.Cartn_y"]),
                     float(r["_atom_site.Cartn_z"])])
            except (KeyError, ValueError):
                continue
    return ca


def turns(ca, r, w3, w4):
    n = 0
    for j in range(r - 4, r + 1):
        if all(k in ca for k in (j, j + 3, j + 4)):
            d3 = float(np.linalg.norm(ca[j] - ca[j + 3]))
            d4 = float(np.linalg.norm(ca[j] - ca[j + 4]))
            if w3[0] <= d3 <= w3[1] and w4[0] <= d4 <= w4[1]:
                n += 1
    return n


def segments(ca, w3, w4, need, bridge=0):
    hel = sorted(r for r in ca if turns(ca, r, w3, w4) >= need)
    if not hel:
        return []
    out = [[hel[0], hel[0]]]
    for r in hel[1:]:
        if r - out[-1][1] - 1 <= bridge:
            out[-1][1] = r
        else:
            out.append([r, r])
    return [tuple(x) for x in out]


def seg_of(segs, r, slack=1):
    for rr in range(r - slack, r + slack + 1):
        for i, (s, e) in enumerate(segs):
            if s <= rr <= e:
                return i
    return None


targets = {"ref_PmTYR": parse_cif_ca(f"{BASE}/B2ZB02_taxID_1404_model.cif")}
print("PmTYR is scored on the B2ZB02 AF3 model (what the pipeline actually uses);")
print("the other four on their crystal PDBs.\n")
for name in CORE:
    if name not in targets:
        targets[name] = parse_pdb_ca(f"{FS}/{name}.pdb")

for label, w3, w4, need in DEFS:
    print(f"=== {label} ===")
    for name, ca in targets.items():
        segs = segments(ca, w3, w4, need)
        idx = [seg_of(segs, a) for a in ANCHORS[name]]
        ok = len({i for i in idx if i is not None}) == 4
        spans = []
        for a, i in zip(ANCHORS[name], idx):
            spans.append(f"{a}->{'none' if i is None else f'{segs[i][0]}-{segs[i][1]}'}")
        print(f"  {name:20} {'SEPARATE' if ok else 'MERGED  '}  n_seg={len(segs):3}  "
              + "  ".join(spans))
    print()
