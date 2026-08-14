"""
Compare per-His pLDDT between AF3 (with Cu) and AFDB/AF2 (no Cu) structures.
Downloads AFDB CIFs and extracts pLDDT at the same His positions identified
from AF3 Cu coordination.
"""
import csv, math, sys, os, urllib.request
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


def get_his_plddt_af3(cif_path):
    """Get per-His pLDDT from AF3 CIF using Cu coordination."""
    atoms = parse_atoms(cif_path)
    cu = [a for a in atoms if a["elem"] == "CU"]
    if len(cu) < 2:
        return None
    cu.sort(key=lambda a: a.get("seq", ""))
    cu1, cu2 = cu[0], cu[1]

    his_ne2 = [a for a in atoms if a["resn"] == "HIS" and a["atom"] == "NE2"]
    coord_his = []
    seen = set()
    for h in his_ne2:
        if min(dist(h, cu1), dist(h, cu2)) <= 3.5 and h["seq"] not in seen:
            seen.add(h["seq"])
            coord_his.append(h)
    coord_his.sort(key=lambda a: int(a["seq"]) if a["seq"].isdigit() else 0)

    if len(coord_his) < 6:
        return None

    ca_bfac = {}
    for a in atoms:
        if a["atom"] == "CA" and a["seq"].isdigit():
            ca_bfac[int(a["seq"])] = a["bfactor"]

    result = []
    for h in coord_his[:6]:
        seq = int(h["seq"])
        plddt = ca_bfac.get(seq, h["bfactor"])
        result.append({"seq": seq, "plddt": plddt})
    return result


def get_plddt_at_positions(cif_path, positions):
    """Get pLDDT (B-factor of CA) at specific residue positions."""
    atoms = parse_atoms(cif_path)
    ca_bfac = {}
    for a in atoms:
        if a["atom"] == "CA" and a["seq"].isdigit():
            ca_bfac[int(a["seq"])] = a["bfactor"]

    result = []
    for pos in positions:
        plddt = ca_bfac.get(pos)
        result.append({"seq": pos, "plddt": plddt})
    return result


def download_afdb(acc, outdir):
    """Download AFDB CIF for a UniProt accession."""
    outpath = os.path.join(outdir, f"AF-{acc}-F1-model_v4.cif")
    if os.path.exists(outpath):
        return outpath
    url = f"https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.cif"
    try:
        urllib.request.urlretrieve(url, outpath)
        return outpath
    except Exception:
        return None


def process_one(args):
    acc, af3_cif, afdb_dir = args
    af3_his = get_his_plddt_af3(af3_cif)
    if af3_his is None:
        return None

    positions = [h["seq"] for h in af3_his]
    afdb_cif = download_afdb(acc, afdb_dir)
    if afdb_cif is None:
        return {"accession": acc, "status": "no_afdb"}

    afdb_his = get_plddt_at_positions(afdb_cif, positions)

    return {
        "accession": acc,
        "status": "ok",
        "his_positions": ",".join(str(p) for p in positions),
        "af3_plddts": ",".join(f"{h['plddt']:.1f}" for h in af3_his),
        "afdb_plddts": ",".join(f"{h['plddt']:.1f}" if h["plddt"] is not None else "NA" for h in afdb_his),
        "af3_mean": sum(h["plddt"] for h in af3_his) / 6,
        "af3_min": min(h["plddt"] for h in af3_his),
        "afdb_mean": sum(h["plddt"] for h in afdb_his if h["plddt"] is not None) / max(1, sum(1 for h in afdb_his if h["plddt"] is not None)),
        "afdb_min": min((h["plddt"] for h in afdb_his if h["plddt"] is not None), default=None),
        "afdb_missing": sum(1 for h in afdb_his if h["plddt"] is None),
    }


def main():
    import argparse, random
    parser = argparse.ArgumentParser()
    parser.add_argument("--cifs", required=True)
    parser.add_argument("--canonical-csv", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--afdb-dir", required=True)
    parser.add_argument("--sample", type=int, default=200)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    os.makedirs(args.afdb_dir, exist_ok=True)

    with open(args.canonical_csv) as f:
        canon = [r for r in csv.DictReader(f) if r["canonical"] == "True"]
    print(f"Canonical structures: {len(canon)}", flush=True)

    random.seed(args.seed)
    sample = random.sample(canon, min(args.sample, len(canon)))
    print(f"Sample size: {len(sample)}", flush=True)

    cif_dir = Path(args.cifs)
    acc_to_cif = {}
    for p in cif_dir.glob("*.cif"):
        acc = p.name.split("_taxID_")[0]
        acc_to_cif[acc] = str(p)

    work = [(r["accession"], acc_to_cif[r["accession"]], args.afdb_dir)
            for r in sample if r["accession"] in acc_to_cif]
    print(f"With CIF files: {len(work)}", flush=True)

    fieldnames = [
        "accession", "status", "his_positions",
        "af3_plddts", "afdb_plddts",
        "af3_mean", "af3_min", "afdb_mean", "afdb_min", "afdb_missing",
    ]

    results = []
    for i, w in enumerate(work):
        r = process_one(w)
        if r:
            results.append(r)
        if (i + 1) % 50 == 0:
            print(f"  {i+1}/{len(work)}", flush=True)

    with open(args.output, "w", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow({k: r.get(k, "") for k in fieldnames})

    ok = [r for r in results if r["status"] == "ok" and r["afdb_missing"] == 0]
    no_afdb = [r for r in results if r["status"] == "no_afdb"]
    missing = [r for r in results if r["status"] == "ok" and r["afdb_missing"] > 0]

    print(f"\n=== Summary ===")
    print(f"Compared: {len(ok)}")
    print(f"No AFDB entry: {len(no_afdb)}")
    print(f"AFDB missing His positions: {len(missing)}")

    if ok:
        af3_means = [r["af3_mean"] for r in ok]
        afdb_means = [r["afdb_mean"] for r in ok]
        af3_mins = [r["af3_min"] for r in ok]
        afdb_mins = [r["afdb_min"] for r in ok]
        print(f"\nAF3  mean His pLDDT: {sum(af3_means)/len(af3_means):.1f} (avg of per-structure means)")
        print(f"AFDB mean His pLDDT: {sum(afdb_means)/len(afdb_means):.1f}")
        print(f"AF3  min His pLDDT:  {sum(af3_mins)/len(af3_mins):.1f} (avg of per-structure mins)")
        print(f"AFDB min His pLDDT:  {sum(afdb_mins)/len(afdb_mins):.1f}")

        better = sum(1 for a3, a2 in zip(af3_means, afdb_means) if a3 > a2)
        worse = sum(1 for a3, a2 in zip(af3_means, afdb_means) if a3 < a2)
        equal = sum(1 for a3, a2 in zip(af3_means, afdb_means) if abs(a3 - a2) < 0.5)
        print(f"\nAF3 better: {better}, AF3 worse: {worse}, ~equal: {equal}")

        diffs = [a3 - a2 for a3, a2 in zip(af3_means, afdb_means)]
        diffs.sort()
        print(f"pLDDT difference (AF3-AFDB): min={diffs[0]:.1f}, median={diffs[len(diffs)//2]:.1f}, max={diffs[-1]:.1f}")

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
