#!/usr/bin/env python3
"""Run AllMetal3D on a batch of CIF structures and extract Cu predictions.

Usage: python run_metal3d.py --cif-dir DIR --output results.tsv [--start N --end M]
"""

import argparse
import csv
import glob
import os
import shutil
import subprocess
import sys
import tempfile

import gemmi


def cif_to_pdb(cif_path, pdb_path):
    """Convert CIF to PDB, stripping metals (Metal3D expects apo structure)."""
    doc = gemmi.cif.read(cif_path)
    block = doc.sole_block()
    st = gemmi.make_structure_from_block(block)

    metals = {"CU", "ZN", "FE", "MN", "CO", "NI", "MG", "CA"}
    for model in st:
        for chain in model:
            to_remove = []
            for i, res in enumerate(chain):
                if res.name in metals:
                    to_remove.append(i)
            for i in reversed(to_remove):
                del chain[i]

    st.write_pdb(pdb_path)
    return True


def parse_probes(probe_pdb):
    """Extract predicted metal sites from Metal3D probe PDB."""
    sites = []
    if not os.path.exists(probe_pdb):
        return sites
    with open(probe_pdb) as f:
        for line in f:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    occupancy = float(line[54:60])
                except ValueError:
                    continue
                element = line[76:80].strip() if len(line) > 78 else ""
                sites.append({"x": x, "y": y, "z": z,
                              "probability": occupancy, "element": element})
    return sites


def get_af3_cu_positions(cif_path):
    """Extract Cu atom positions from AF3 CIF."""
    doc = gemmi.cif.read(cif_path)
    block = doc.sole_block()
    st = gemmi.make_structure_from_block(block)

    cu_positions = []
    for model in st:
        for chain in model:
            for res in chain:
                if res.name == "CU":
                    for atom in res:
                        if atom.element.name == "Cu":
                            cu_positions.append({
                                "x": atom.pos.x,
                                "y": atom.pos.y,
                                "z": atom.pos.z,
                                "bfactor": atom.b_iso
                            })
    return cu_positions


def distance(a, b):
    return ((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)**0.5


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cif-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--end", type=int, default=None)
    args = parser.parse_args()

    cif_files = sorted(glob.glob(os.path.join(args.cif_dir, "*.cif")))
    if args.end is not None:
        cif_files = cif_files[args.start:args.end]
    else:
        cif_files = cif_files[args.start:]

    print(f"Processing {len(cif_files)} structures (index {args.start}-{args.start+len(cif_files)-1})")

    fieldnames = ["accession", "status", "af3_cu_index",
                  "af3_cu_x", "af3_cu_y", "af3_cu_z", "af3_cu_plddt",
                  "n_metal3d_sites", "n_metal3d_cu_sites",
                  "closest_metal3d_dist", "closest_metal3d_prob",
                  "closest_metal3d_element",
                  "closest_cu_dist", "closest_cu_prob", "error"]
    outfile = open(args.output, "w", newline="")
    writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t",
                            extrasaction="ignore")
    writer.writeheader()
    n_written = 0

    def write_row(row):
        writer.writerow(row)
        outfile.flush()
        nonlocal n_written
        n_written += 1

    for i, cif_path in enumerate(cif_files):
        basename = os.path.basename(cif_path).replace(".cif", "")
        accession = basename.split("_taxID_")[0] if "_taxID_" in basename else basename
        print(f"  [{i+1}/{len(cif_files)}] {accession}...", end=" ", flush=True)

        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path = os.path.join(tmpdir, "input.pdb")
            out_dir = os.path.join(tmpdir, "output")
            os.makedirs(out_dir)

            try:
                cif_to_pdb(cif_path, pdb_path)
            except Exception as e:
                print(f"CIF conversion failed: {e}")
                write_row({"accession": accession, "status": "cif_error",
                           "error": str(e)})
                continue

            try:
                result = subprocess.run(
                    ["allmetal3d", "-i", pdb_path, "-o", out_dir,
                     "--models", "allmetal3d", "-m", "fast", "-p", "0.1"],
                    capture_output=True, text=True, timeout=1200
                )
                if result.returncode != 0:
                    print(f"metal3d exit {result.returncode}: {result.stderr[-200:]}")
                    write_row({"accession": accession, "status": "metal3d_error",
                               "error": result.stderr[-200:]})
                    continue
            except subprocess.TimeoutExpired:
                print("timeout")
                write_row({"accession": accession, "status": "timeout"})
                continue
            except Exception as e:
                print(f"metal3d error: {e}")
                write_row({"accession": accession, "status": "metal3d_error",
                           "error": str(e)})
                continue

            # Find metal prediction output (*_metals.pdb)
            metal_files = glob.glob(os.path.join(out_dir, "*_metals.pdb"))
            if not metal_files:
                metal_files = glob.glob(os.path.join(out_dir, "*.pdb"))
            all_predicted = []
            for pf in metal_files:
                all_predicted.extend(parse_probes(pf))

            # Get AF3 Cu positions
            af3_cu = get_af3_cu_positions(cif_path)

            # Match: for each AF3 Cu, find closest Metal3D prediction
            cu_predicted = [p for p in all_predicted
                           if p.get("element", "").upper() == "CU"]
            n_metal3d_sites = len(all_predicted)
            n_cu_sites = len(cu_predicted)

            for j, cu in enumerate(af3_cu):
                row = {
                    "accession": accession,
                    "status": "ok",
                    "af3_cu_index": j + 1,
                    "af3_cu_x": f"{cu['x']:.2f}",
                    "af3_cu_y": f"{cu['y']:.2f}",
                    "af3_cu_z": f"{cu['z']:.2f}",
                    "af3_cu_plddt": f"{cu['bfactor']:.1f}",
                    "n_metal3d_sites": n_metal3d_sites,
                    "n_metal3d_cu_sites": n_cu_sites,
                }

                if all_predicted:
                    closest = min(all_predicted, key=lambda p: distance(cu, p))
                    row["closest_metal3d_dist"] = f"{distance(cu, closest):.2f}"
                    row["closest_metal3d_prob"] = f"{closest['probability']:.3f}"
                    row["closest_metal3d_element"] = closest.get("element", "")

                if cu_predicted:
                    closest_cu = min(cu_predicted, key=lambda p: distance(cu, p))
                    row["closest_cu_dist"] = f"{distance(cu, closest_cu):.2f}"
                    row["closest_cu_prob"] = f"{closest_cu['probability']:.3f}"

                write_row(row)

            if not af3_cu:
                write_row({"accession": accession, "status": "no_af3_cu"})

            print(f"ok ({n_metal3d_sites} metal sites, {n_cu_sites} Cu)")

    outfile.close()
    print(f"\nDone: {n_written} rows written to {args.output}")


if __name__ == "__main__":
    main()
