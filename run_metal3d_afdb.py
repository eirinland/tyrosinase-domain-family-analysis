"""
Run AllMetal3D on AFDB structures to predict Cu positions,
then compare predicted Cu-Cu distance with AF3 Cu-Cu distance.
"""
import argparse, csv, glob, os, subprocess, tempfile
import gemmi


def cif_to_pdb(cif_path, pdb_path):
    doc = gemmi.cif.read(cif_path)
    block = doc.sole_block()
    st = gemmi.make_structure_from_block(block)
    metals = {"CU", "ZN", "FE", "MN", "CO", "NI", "MG", "CA"}
    for model in st:
        for chain in model:
            to_remove = [i for i, res in enumerate(chain) if res.name in metals]
            for i in reversed(to_remove):
                del chain[i]
    st.write_pdb(pdb_path)


def parse_probes(probe_pdb):
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
                    prob = float(line[54:60])
                except ValueError:
                    continue
                elem = line[76:80].strip() if len(line) > 78 else ""
                sites.append({"x": x, "y": y, "z": z, "prob": prob, "elem": elem})
    return sites


def dist(a, b):
    return ((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)**0.5


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--afdb-dir", required=True)
    parser.add_argument("--af3-cu-csv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--end", type=int, default=None)
    args = parser.parse_args()

    af3_cu = {}
    with open(args.af3_cu_csv) as f:
        for r in csv.DictReader(f):
            if r["canonical"] == "True" and r["cu_dist"]:
                af3_cu[r["accession"]] = float(r["cu_dist"])

    cif_files = sorted(glob.glob(os.path.join(args.afdb_dir, "*.cif")))
    if args.end is not None:
        cif_files = cif_files[args.start:args.end]
    else:
        cif_files = cif_files[args.start:]
    print(f"Processing {len(cif_files)} structures (index {args.start}-{args.start+len(cif_files)-1})", flush=True)

    fieldnames = [
        "accession", "status",
        "af3_cu_dist",
        "n_predicted_cu", "predicted_cu_dist",
        "cu_dist_diff",
        "n_all_metals",
        "cu1_prob", "cu2_prob",
        "error",
    ]
    fout = open(args.output, "w", newline="")
    writer = csv.DictWriter(fout, fieldnames=fieldnames)
    writer.writeheader()

    for i, cif_path in enumerate(cif_files):
        basename = os.path.basename(cif_path)
        acc = basename.replace("AF-", "").replace("-F1-model_v6.cif", "").replace("-F1-model_v4.cif", "")

        row = {"accession": acc, "af3_cu_dist": af3_cu.get(acc, "")}

        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path = os.path.join(tmpdir, "input.pdb")
            out_dir = os.path.join(tmpdir, "output")
            os.makedirs(out_dir)

            try:
                cif_to_pdb(cif_path, pdb_path)
            except Exception as e:
                row["status"] = "cif_error"
                row["error"] = str(e)
                writer.writerow(row)
                fout.flush()
                print(f"  [{i+1}] {acc}: CIF error: {e}", flush=True)
                continue

            try:
                result = subprocess.run(
                    ["allmetal3d", "-i", pdb_path, "-o", out_dir,
                     "--models", "allmetal3d", "-m", "fast", "-p", "0.1"],
                    capture_output=True, text=True, timeout=1200
                )
                if result.returncode != 0:
                    row["status"] = "metal3d_error"
                    row["error"] = result.stderr[-300:].replace("\n", " ")
                    writer.writerow(row)
                    fout.flush()
                    print(f"  [{i+1}] {acc}: metal3d exit {result.returncode}: {result.stderr[-200:]}", flush=True)
                    continue
            except subprocess.TimeoutExpired:
                row["status"] = "timeout"
                writer.writerow(row)
                fout.flush()
                print(f"  [{i+1}] {acc}: timeout", flush=True)
                continue

            metal_files = glob.glob(os.path.join(out_dir, "*_metals.pdb"))
            if not metal_files:
                metal_files = glob.glob(os.path.join(out_dir, "*.pdb"))
            predicted = []
            for pf in metal_files:
                predicted.extend(parse_probes(pf))

            cu_sites = sorted(
                [p for p in predicted if p["elem"].upper() == "CU"],
                key=lambda p: -p["prob"]
            )

            row["n_all_metals"] = len(predicted)
            row["n_predicted_cu"] = len(cu_sites)

            if len(cu_sites) >= 2:
                row["cu1_prob"] = f"{cu_sites[0]['prob']:.3f}"
                row["cu2_prob"] = f"{cu_sites[1]['prob']:.3f}"
                cu_cu = dist(cu_sites[0], cu_sites[1])
                row["predicted_cu_dist"] = f"{cu_cu:.2f}"
                if acc in af3_cu:
                    row["cu_dist_diff"] = f"{cu_cu - af3_cu[acc]:.2f}"
                row["status"] = "ok"
            elif len(cu_sites) == 1:
                row["cu1_prob"] = f"{cu_sites[0]['prob']:.3f}"
                row["status"] = "only_1cu"
            else:
                row["status"] = "no_cu"

        writer.writerow(row)
        fout.flush()
        print(f"  [{i+1}/{len(cif_files)}] {acc}: {row['status']} "
              f"({row.get('n_predicted_cu', 0)} Cu)", flush=True)

    fout.close()
    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
