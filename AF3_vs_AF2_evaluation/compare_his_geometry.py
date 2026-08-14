#!/usr/bin/env python3
"""Compare CuA/CuB His-cluster geometry between AF3 (with Cu) and AFDB (without Cu).

For each canonical AF3 structure:
1. Find 2 Cu, identify 6 coordinating His, assign 3 to CuA and 3 to CuB
2. Download matching AFDB CIF
3. In AFDB CIF: find same His residues by seq_id, measure NE2 centroid-centroid distance
4. Compare AF3 Cu-Cu distance vs AFDB His-centroid distance
"""

import csv
import math
import os
import sys
import urllib.request
import urllib.error
from multiprocessing import Pool
from pathlib import Path


def parse_atoms(path):
    atoms = []
    col_names = []
    in_atom = False
    collecting = False
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if line == "loop_":
                collecting = False
                in_atom = False
                col_names = []
                continue
            if line.startswith("_atom_site."):
                collecting = True
                col_names.append(line)
                continue
            if collecting and col_names:
                collecting = False
                in_atom = True
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
    return math.sqrt((a["x"] - b["x"])**2 + (a["y"] - b["y"])**2 + (a["z"] - b["z"])**2)


def centroid(atom_list):
    n = len(atom_list)
    if n == 0:
        return None
    return {
        "x": sum(a["x"] for a in atom_list) / n,
        "y": sum(a["y"] for a in atom_list) / n,
        "z": sum(a["z"] for a in atom_list) / n,
    }


def get_his_cu_assignment(af3_path):
    """From AF3 structure, get CuA/CuB His assignments.
    Returns: (cu_cu_dist, cuA_his_seqs, cuB_his_seqs, af3_centroid_dist) or None."""
    atoms = parse_atoms(af3_path)
    cu = [a for a in atoms if a["elem"] == "CU"]
    if len(cu) < 2:
        return None

    cu.sort(key=lambda a: a.get("seq", ""))
    cu1, cu2 = cu[0], cu[1]
    cu_dist = dist(cu1, cu2)
    if not (2.8 <= cu_dist <= 5.5):
        return None

    his_ne2 = [a for a in atoms if a["resn"] == "HIS" and a["atom"] == "NE2"]
    coord_his = []
    seen_seq = set()
    for h in his_ne2:
        if h["seq"] in seen_seq:
            continue
        d1, d2 = dist(h, cu1), dist(h, cu2)
        if min(d1, d2) <= 3.5:
            seen_seq.add(h["seq"])
            cu_label = "A" if d1 < d2 else "B"
            coord_his.append((h, cu_label))

    cuA_his = [h for h, l in coord_his if l == "A"]
    cuB_his = [h for h, l in coord_his if l == "B"]

    if len(cuA_his) < 2 or len(cuB_his) < 2:
        return None

    cA = centroid(cuA_his)
    cB = centroid(cuB_his)
    af3_centroid_dist = dist(cA, cB) if cA and cB else None

    return {
        "cu_cu_dist": cu_dist,
        "cuA_seqs": [h["seq"] for h in cuA_his],
        "cuB_seqs": [h["seq"] for h in cuB_his],
        "n_cuA_his": len(cuA_his),
        "n_cuB_his": len(cuB_his),
        "af3_centroid_dist": af3_centroid_dist,
    }


def download_afdb(acc, output_dir):
    out_path = os.path.join(output_dir, f"AF-{acc}-F1.cif")
    if os.path.exists(out_path):
        return out_path
    for ver in ("v6", "v4", "v3"):
        url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_{ver}.cif"
        try:
            urllib.request.urlretrieve(url, out_path)
            return out_path
        except (urllib.error.HTTPError, urllib.error.URLError):
            continue
    return None


def measure_afdb_his_distance(afdb_path, cuA_seqs, cuB_seqs):
    """Measure NE2 centroid-centroid distance in AFDB structure using known His assignments."""
    atoms = parse_atoms(afdb_path)
    his_ne2 = {a["seq"]: a for a in atoms if a["resn"] == "HIS" and a["atom"] == "NE2"}
    his_ca = {a["seq"]: a for a in atoms if a["resn"] == "HIS" and a["atom"] == "CA"}

    cuA_atoms = [his_ne2[s] for s in cuA_seqs if s in his_ne2]
    cuB_atoms = [his_ne2[s] for s in cuB_seqs if s in his_ne2]

    if len(cuA_atoms) < 2 or len(cuB_atoms) < 2:
        return None

    cA = centroid(cuA_atoms)
    cB = centroid(cuB_atoms)

    cuA_ca = [his_ca[s] for s in cuA_seqs if s in his_ca]
    cuB_ca = [his_ca[s] for s in cuB_seqs if s in his_ca]
    cuA_plddt = sum(a["bfactor"] for a in cuA_ca) / len(cuA_ca) if cuA_ca else 0
    cuB_plddt = sum(a["bfactor"] for a in cuB_ca) / len(cuB_ca) if cuB_ca else 0

    return {
        "afdb_centroid_dist": dist(cA, cB),
        "afdb_cuA_found": len(cuA_atoms),
        "afdb_cuB_found": len(cuB_atoms),
        "afdb_cuA_plddt": cuA_plddt,
        "afdb_cuB_plddt": cuB_plddt,
    }


def process_one(args):
    acc, af3_path, afdb_dir = args
    info = get_his_cu_assignment(af3_path)
    if info is None:
        return None

    afdb_path = download_afdb(acc, afdb_dir)
    if afdb_path is None:
        return {"accession": acc, "afdb_available": False, **info}

    afdb_result = measure_afdb_his_distance(afdb_path, info["cuA_seqs"], info["cuB_seqs"])
    if afdb_result is None:
        return {"accession": acc, "afdb_available": True, "afdb_his_missing": True, **info}

    delta = afdb_result["afdb_centroid_dist"] - info["af3_centroid_dist"]
    return {
        "accession": acc,
        "afdb_available": True,
        "afdb_his_missing": False,
        "af3_cu_cu": round(info["cu_cu_dist"], 2),
        "af3_centroid": round(info["af3_centroid_dist"], 2),
        "afdb_centroid": round(afdb_result["afdb_centroid_dist"], 2),
        "delta_centroid": round(delta, 2),
        "n_cuA_his": info["n_cuA_his"],
        "n_cuB_his": info["n_cuB_his"],
        "afdb_cuA_found": afdb_result["afdb_cuA_found"],
        "afdb_cuB_found": afdb_result["afdb_cuB_found"],
        "afdb_cuA_plddt": round(afdb_result["afdb_cuA_plddt"], 1),
        "afdb_cuB_plddt": round(afdb_result["afdb_cuB_plddt"], 1),
    }


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--canonical-csv", required=True,
                        help="canonical_criteria_all_ca.csv")
    parser.add_argument("--af3-dir", required=True,
                        help="Directory with AF3 CIF files (sqsh mount)")
    parser.add_argument("--afdb-dir", required=True,
                        help="Temporary directory for downloaded AFDB CIFs")
    parser.add_argument("--output", required=True,
                        help="Output CSV path")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()

    os.makedirs(args.afdb_dir, exist_ok=True)

    canonical_accs = []
    with open(args.canonical_csv) as f:
        for row in csv.DictReader(f):
            if row.get("canonical", "").lower() == "true":
                canonical_accs.append(row["accession"])
    print(f"Canonical structures: {len(canonical_accs)}", flush=True)

    af3_files = {}
    for fname in os.listdir(args.af3_dir):
        if fname.endswith(".cif"):
            acc = fname.split("_taxID_")[0]
            af3_files[acc] = os.path.join(args.af3_dir, fname)

    tasks = []
    for acc in canonical_accs:
        if acc in af3_files:
            tasks.append((acc, af3_files[acc], args.afdb_dir))

    print(f"Processing {len(tasks)} structures with AF3 CIFs...", flush=True)

    fieldnames = [
        "accession", "afdb_available", "afdb_his_missing",
        "af3_cu_cu", "af3_centroid", "afdb_centroid", "delta_centroid",
        "n_cuA_his", "n_cuB_his", "afdb_cuA_found", "afdb_cuB_found",
        "afdb_cuA_plddt", "afdb_cuB_plddt",
    ]

    results = []
    n_afdb_ok = 0
    with Pool(args.workers) as pool:
        for i, result in enumerate(pool.imap_unordered(process_one, tasks, chunksize=50)):
            if result is not None:
                results.append(result)
                if result.get("afdb_available"):
                    n_afdb_ok += 1
            if (i + 1) == 100 and n_afdb_ok == 0:
                print("WARNING: 0/100 AFDB downloads succeeded — check URL or network", flush=True)
                sys.exit(1)
            if (i + 1) % 1000 == 0:
                print(f"  {i+1}/{len(tasks)} done, {len(results)} results, "
                      f"{n_afdb_ok} AFDB downloaded", flush=True)

    results.sort(key=lambda r: -abs(r.get("delta_centroid", 0)))

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow({k: r.get(k, "") for k in fieldnames})

    n_avail = sum(1 for r in results if r.get("afdb_available"))
    n_with_data = sum(1 for r in results if r.get("afdb_centroid"))
    big_delta = [r for r in results if r.get("delta_centroid", 0) > 3]

    print(f"\nResults: {len(results)} total, {n_avail} AFDB available, {n_with_data} with His data")
    print(f"Structures with delta_centroid > 3 A: {len(big_delta)}")
    if big_delta:
        print("\nTop 20 largest His-centroid shifts (AF3 -> AFDB):")
        for r in big_delta[:20]:
            print(f"  {r['accession']}: AF3 centroid={r['af3_centroid']:.1f} A, "
                  f"AFDB centroid={r['afdb_centroid']:.1f} A, "
                  f"delta={r['delta_centroid']:.1f} A")


if __name__ == "__main__":
    main()
