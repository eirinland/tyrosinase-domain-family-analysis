#!/usr/bin/env python3
"""Check per-helix pLDDT for the 4 core PPO helices using Foldseek alignment mapping."""

import re
import os
import sys
import glob

REF_FIRST_RESI = {
    "ref_PmTYR": 4,
    "ref_8BBR_Vspinosum": 34,
    "ref_2Y9W_Abisporus": 17,
    "ref_5CE9_Jregia": 1,
    "ref_1BT3_Ibatatas": 1,
    "ref_5M8L_human": 81,
    "ref_1JS8_squid": 2503,
    "ref_I3D139_archaea": 6,
    "ref_A0A9N8ELP9_oomycota": 347,
}

_CORE_HELICES_PDB = {
    "ref_PmTYR": {"a1": (34, 46), "a2": (65, 83), "a3": (203, 211), "a4": (226, 244)},
    "ref_8BBR_Vspinosum": {"a1": (73, 83), "a2": (91, 110), "a3": (259, 265), "a4": (279, 300)},
    "ref_2Y9W_Abisporus": {"a1": (54, 61), "a2": (90, 113), "a3": (255, 267), "a4": (291, 309)},
    "ref_5CE9_Jregia": {"a1": (80, 91), "a2": (113, 131), "a3": (241, 246), "a4": (269, 283)},
    "ref_1BT3_Ibatatas": {"a1": (81, 92), "a2": (114, 133), "a3": (240, 247), "a4": (269, 287)},
    "ref_5M8L_human": {"a1": (184, 196), "a2": (220, 239), "a3": (376, 383), "a4": (399, 417)},
    "ref_1JS8_squid": {"a1": (2536, 2543), "a2": (2567, 2584), "a3": (2660, 2679), "a4": (2697, 2719)},
    "ref_I3D139_archaea": {"a1": (36, 47), "a2": (67, 85), "a3": (197, 205), "a4": (222, 240)},
    "ref_A0A9N8ELP9_oomycota": {"a1": (413, 424), "a2": (437, 456), "a3": (590, 601), "a4": (624, 642)},
}

CORE_HELICES = {}
for ref, helices in _CORE_HELICES_PDB.items():
    offset = REF_FIRST_RESI[ref]
    CORE_HELICES[ref] = {h: (s - offset, e - offset) for h, (s, e) in helices.items()}


def parse_cigar(cigar):
    return [(op, int(n)) for n, op in re.findall(r"(\d+)([MID])", cigar)]


def get_target_to_query_map(qstart, tstart, cigar):
    """Map target 0-based positions to query 0-based positions via cigar."""
    ops = parse_cigar(cigar)
    qpos = qstart
    tpos = tstart
    t2q = {}
    for op, length in ops:
        if op == "M":
            for _ in range(length):
                t2q[tpos] = qpos
                qpos += 1
                tpos += 1
        elif op == "D":
            tpos += length
        elif op == "I":
            qpos += length
    return t2q


def read_plddt_from_cif(cif_path):
    """Read per-residue pLDDT (CA atoms) from CIF. Returns list of pLDDT values
    in residue order (0-based index = Foldseek query position)."""
    plddts = []
    in_atom = False
    cols = []
    with open(cif_path) as f:
        for line in f:
            if line.startswith("_atom_site."):
                in_atom = True
                cols.append(line.strip().split(".")[1].split()[0])
                continue
            if in_atom and not line.startswith("_") and not line.startswith("#") and line.strip():
                tokens = line.split()
                bfac_idx = cols.index("B_iso_or_equiv")
                atom_idx = cols.index("label_atom_id")
                seq_idx = cols.index("label_seq_id")
                if tokens[seq_idx] == "." or tokens[atom_idx] != "CA":
                    continue
                plddts.append(float(tokens[bfac_idx]))
            elif in_atom and (line.startswith("#") or not line.strip()):
                break
    return plddts


def main():
    import argparse
    import csv as csvmod
    parser = argparse.ArgumentParser()
    parser.add_argument("tsv", help="Foldseek TSV with header and cigar column")
    parser.add_argument("--cif-dir", required=True, help="Directory with CIF files")
    parser.add_argument("-o", "--output", help="Output CSV path")
    args = parser.parse_args()

    all_hits = {}
    with open(args.tsv) as f:
        header = f.readline().strip().split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < len(header):
                continue
            query = parts[col["query"]]
            target = parts[col["target"]]
            tstart = int(parts[col["tstart"]])
            qstart = int(parts[col["qstart"]])
            qtmscore = float(parts[col["qtmscore"]])
            cigar = parts[col["cigar"]]
            acc = query.split("_taxID_")[0]
            all_hits.setdefault(acc, []).append({
                "query": query, "target": target, "tstart": tstart,
                "qstart": qstart, "qtmscore": qtmscore, "cigar": cigar,
            })

    cif_files = {}
    for fname in os.listdir(args.cif_dir):
        if fname.endswith(".cif"):
            acc = fname.split("_taxID_")[0]
            cif_files[acc] = os.path.join(args.cif_dir, fname)

    rows = []
    for acc in sorted(all_hits.keys()):
        hits = sorted(all_hits[acc], key=lambda x: -x["qtmscore"])
        best_info = hits[0]
        ref = best_info["target"]
        helices = CORE_HELICES.get(ref)
        if not helices or acc not in cif_files:
            continue

        plddt_list = read_plddt_from_cif(cif_files[acc])
        qstart = best_info["qstart"]
        t2q = get_target_to_query_map(qstart, best_info["tstart"], best_info["cigar"])

        helix_plddts = {}
        for hname in ["a1", "a2", "a3", "a4"]:
            h_start, h_end = helices[hname]
            vals = []
            for tpos in range(h_start, h_end + 1):
                qpos = t2q.get(tpos)
                if qpos is not None and 0 <= qpos < len(plddt_list):
                    vals.append(plddt_list[qpos])
            helix_plddts[hname] = sum(vals) / len(vals) if vals else 0.0

        rows.append({
            "accession": acc, "ref": ref,
            "a1_plddt": helix_plddts["a1"], "a2_plddt": helix_plddts["a2"],
            "a3_plddt": helix_plddts["a3"], "a4_plddt": helix_plddts["a4"],
            "min_helix_plddt": min(helix_plddts.values()),
        })

    if args.output:
        fields = ["accession", "ref", "a1_plddt", "a2_plddt", "a3_plddt", "a4_plddt", "min_helix_plddt"]
        with open(args.output, "w", newline="") as f:
            w = csvmod.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for r in rows:
                w.writerow({k: (f"{v:.1f}" if isinstance(v, float) else v) for k, v in r.items()})
        print("Wrote %s (%d structures)" % (args.output, len(rows)))

    lo = sum(1 for r in rows if r["min_helix_plddt"] < 70)
    mid = sum(1 for r in rows if 70 <= r["min_helix_plddt"] < 90)
    hi = sum(1 for r in rows if r["min_helix_plddt"] >= 90)
    print("Helix pLDDT summary: %d structures, %d min>=90, %d 70-90, %d <70" % (len(rows), hi, mid, lo))


if __name__ == "__main__":
    main()
