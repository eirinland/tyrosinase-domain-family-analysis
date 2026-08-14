#!/usr/bin/env python3
"""Investigate canonical structures that fail the 4/4 helix coverage check."""

import csv
import os
import sys

FILTER_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HELIX_CSV = os.path.join(FILTER_DIR, "helix_coverage.csv")
PLDDT_CSV = os.path.join(FILTER_DIR, "helix_plddt.csv")
CANON_CSV = os.path.join(FILTER_DIR, "canonical_trimmed.csv")
PASS2_TSV = os.path.join(FILTER_DIR, "foldseek_pass2.tsv")
CIF_DIR = os.path.join(FILTER_DIR, "cifs_trimmed")


def overall_plddt_from_cif(cif_path):
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
    return sum(plddts) / len(plddts) if plddts else 0.0, len(plddts)


def main():
    helix = {}
    with open(HELIX_CSV) as f:
        for r in csv.DictReader(f):
            helix[r["accession"]] = r

    canon = {}
    with open(CANON_CSV) as f:
        for r in csv.DictReader(f):
            canon[r["accession"]] = r

    targets = [acc for acc in helix if int(helix[acc]["n_helices"]) < 4 and canon.get(acc, {}).get("canonical") == "True"]
    targets.sort()
    print(f"Found {len(targets)} canonical structures with <4/4 helices\n")

    pass2 = {}
    with open(PASS2_TSV) as f:
        for line in f:
            parts = line.strip().split("\t")
            acc = parts[0].split("_taxID_")[0]
            if acc in targets:
                pass2.setdefault(acc, []).append({
                    "target": parts[1], "qstart": int(parts[2]), "qend": int(parts[3]),
                    "qlen": int(parts[4]), "tstart": int(parts[5]), "tend": int(parts[6]),
                    "tlen": int(parts[7]), "alntmscore": float(parts[8]),
                    "qtmscore": float(parts[9]), "alnlen": int(parts[12]),
                })

    out_rows = []
    for acc in targets:
        h = helix[acc]
        c = canon[acc]
        cif_path = None
        for fname in os.listdir(CIF_DIR):
            if fname.startswith(acc + "_"):
                cif_path = os.path.join(CIF_DIR, fname)
                break

        avg_plddt, n_res = 0.0, 0
        if cif_path:
            avg_plddt, n_res = overall_plddt_from_cif(cif_path)

        hits = pass2.get(acc, [])
        best_hit = max(hits, key=lambda x: x["qtmscore"]) if hits else {}

        row = {
            "accession": acc,
            "n_helices": h["n_helices"],
            "ref": h["ref"],
            "qtmscore": h["qtmscore"],
            "a1_cov": h["a1_cov"],
            "a2_cov": h["a2_cov"],
            "a3_cov": h["a3_cov"],
            "a4_cov": h["a4_cov"],
            "cu_dist": c.get("cu_dist", ""),
            "n_coord_his": c.get("n_coord_his", ""),
            "min_his_plddt": c.get("min_plddt", ""),
            "trimmed_res": n_res,
            "overall_plddt": f"{avg_plddt:.1f}",
            "alnlen": best_hit.get("alnlen", ""),
            "qlen": best_hit.get("qlen", ""),
        }
        out_rows.append(row)

        print(f"--- {acc} ---")
        print(f"  Helices: {h['n_helices']}/4  ref={h['ref']}  qtmscore={h['qtmscore']}")
        print(f"  Coverage: a1={h['a1_cov']} a2={h['a2_cov']} a3={h['a3_cov']} a4={h['a4_cov']}")
        print(f"  Canonical: Cu-Cu={c.get('cu_dist','?')} A, {c.get('n_coord_his','?')} His, min_plddt={c.get('min_plddt','?')}")
        print(f"  Trimmed: {n_res} res, overall pLDDT={avg_plddt:.1f}")
        if best_hit:
            print(f"  Best hit: {best_hit['target']} qtm={best_hit['qtmscore']:.3f} alnlen={best_hit['alnlen']} qlen={best_hit['qlen']}")
        print()

    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "canonical_low_helix.csv")
    fields = list(out_rows[0].keys())
    with open(outpath, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(out_rows)
    print(f"Wrote {outpath}")


if __name__ == "__main__":
    main()
