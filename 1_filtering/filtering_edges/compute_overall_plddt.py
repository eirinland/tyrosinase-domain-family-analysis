#!/usr/bin/env python3
"""Compute overall pLDDT for all trimmed PPO core structures."""

import os
import sys
import csv

CIF_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cifs_trimmed")


def plddt_from_cif(cif_path):
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
    files = sorted(f for f in os.listdir(CIF_DIR) if f.endswith(".cif"))
    print(f"Processing {len(files)} CIF files", file=sys.stderr)

    outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), "overall_plddt.csv")
    with open(outpath, "w", newline="") as out:
        w = csv.writer(out)
        w.writerow(["accession", "trimmed_res", "overall_plddt"])
        for i, fname in enumerate(files):
            acc = fname.split("_taxID_")[0]
            avg, nres = plddt_from_cif(os.path.join(CIF_DIR, fname))
            w.writerow([acc, nres, f"{avg:.1f}"])
            if (i + 1) % 5000 == 0:
                print(f"  {i+1}/{len(files)}", file=sys.stderr)

    print(f"Wrote {outpath}", file=sys.stderr)

    # Quick summary
    import statistics
    vals = []
    with open(outpath) as f:
        for r in csv.DictReader(f):
            vals.append(float(r["overall_plddt"]))
    print(f"\nOverall pLDDT summary ({len(vals)} structures):")
    print(f"  >=90: {sum(1 for v in vals if v >= 90)}")
    print(f"  70-90: {sum(1 for v in vals if 70 <= v < 90)}")
    print(f"  50-70: {sum(1 for v in vals if 50 <= v < 70)}")
    print(f"  <50: {sum(1 for v in vals if v < 50)}")
    print(f"  median: {statistics.median(vals):.1f}")


if __name__ == "__main__":
    main()
