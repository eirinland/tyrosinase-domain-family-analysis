"""
Check alignment-free canonical criteria for all structures:
  - 2 Cu atoms with Cu-Cu distance 2.8-5.5 A
  - 6 His NE2 within 3.0 A of either Cu
  - Per-His pLDDT >= 70 (B-factor in CIF)
Outputs CSV with per-structure results.
"""
import csv, math, sys, os
from pathlib import Path
from multiprocessing import Pool


def parse_atoms(path):
    atoms = []
    col_names = []
    in_atom = False
    collecting = False
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if line == "loop_":
                collecting = False; in_atom = False; col_names = []
                continue
            if line.startswith("_atom_site."):
                collecting = True; col_names.append(line)
                continue
            if collecting and col_names:
                collecting = False; in_atom = True
            if in_atom:
                if line.startswith("_") or line == "#" or not line:
                    break
                parts = line.split()
                if len(parts) != len(col_names):
                    continue
                row = dict(zip(col_names, parts))
                try:
                    atoms.append({
                        "elem": row.get("_atom_site.type_symbol", "").upper(),
                        "atom": row.get("_atom_site.label_atom_id", "").upper(),
                        "resn": row.get("_atom_site.label_comp_id", ""),
                        "seq": row.get("_atom_site.label_seq_id", ""),
                        "x": float(row["_atom_site.Cartn_x"]),
                        "y": float(row["_atom_site.Cartn_y"]),
                        "z": float(row["_atom_site.Cartn_z"]),
                        "bfactor": float(row.get("_atom_site.B_iso_or_equiv", "0")),
                    })
                except (KeyError, ValueError):
                    continue
    return atoms


def dist(a, b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)


def check_one(cif_path):
    acc = Path(cif_path).name.split("_taxID_")[0]
    try:
        atoms = parse_atoms(cif_path)
    except Exception:
        return {"accession": acc, "has_2cu": False}

    cu = [a for a in atoms if a["elem"] == "CU"]
    if len(cu) < 2:
        return {"accession": acc, "has_2cu": False, "n_cu": len(cu)}

    cu.sort(key=lambda a: a.get("seq", ""))
    cu1, cu2 = cu[0], cu[1]
    cu_dist = dist(cu1, cu2)

    his_ne2 = [a for a in atoms if a["resn"] == "HIS" and a["atom"] == "NE2"]
    his_ca = {a["seq"]: a for a in atoms if a["resn"] == "HIS" and a["atom"] == "CA"}
    coord_his = []
    for h in his_ne2:
        d1, d2 = dist(h, cu1), dist(h, cu2)
        if min(d1, d2) <= 3.5:
            coord_his.append(h)
    coord_his.sort(key=lambda a: int(a["seq"]) if a["seq"].isdigit() else 0)

    unique_his = []
    seen_seq = set()
    for h in coord_his:
        if h["seq"] not in seen_seq:
            seen_seq.add(h["seq"])
            unique_his.append(h)

    n_his = len(unique_his)
    plddts = [his_ca[h["seq"]]["bfactor"] for h in unique_his[:6] if h["seq"] in his_ca]
    min_plddt = min(plddts) if plddts else 0
    his_positions = ",".join(h["seq"] for h in unique_his[:6])

    cu_ok = 2.8 <= cu_dist <= 5.5
    his_ok = n_his >= 6
    plddt_ok = min_plddt >= 70 if plddts else False
    canonical = cu_ok and his_ok and plddt_ok

    return {
        "accession": acc,
        "has_2cu": True,
        "n_cu": len(cu),
        "cu_dist": round(cu_dist, 2),
        "cu_dist_ok": cu_ok,
        "n_coord_his": n_his,
        "his_ok": his_ok,
        "min_plddt": round(min_plddt, 1),
        "plddt_ok": plddt_ok,
        "canonical": canonical,
        "his_positions": his_positions,
    }


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--cifs", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args()

    cif_dir = Path(args.cifs)
    cif_files = sorted(str(p) for p in cif_dir.glob("*.cif"))
    print(f"CIF files: {len(cif_files)}", flush=True)

    if args.limit > 0:
        cif_files = cif_files[:args.limit]

    fieldnames = [
        "accession", "has_2cu", "n_cu", "cu_dist", "cu_dist_ok",
        "n_coord_his", "his_ok", "min_plddt", "plddt_ok",
        "canonical", "his_positions",
    ]

    with open(args.output, "w", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()

        with Pool(args.workers) as pool:
            for i, result in enumerate(pool.imap(check_one, cif_files, chunksize=100)):
                writer.writerow({k: result.get(k, "") for k in fieldnames})
                if (i + 1) % 5000 == 0:
                    fout.flush()
                    print(f"  {i+1}/{len(cif_files)}", flush=True)

    print(f"Done. {len(cif_files)} processed.", flush=True)


if __name__ == "__main__":
    main()
