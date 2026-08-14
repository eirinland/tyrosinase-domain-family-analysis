import csv, math, glob
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

def dist3d(a, b):
    return math.sqrt((a["x"]-b["x"])**2 + (a["y"]-b["y"])**2 + (a["z"]-b["z"])**2)

def check_one(args):
    cif_path, acc = args
    try:
        atoms = parse_atoms(cif_path)
    except Exception:
        return None
    cu = [a for a in atoms if a["elem"] == "CU"]
    cu.sort(key=lambda a: a.get("chain", "A"))
    if len(cu) < 2:
        return None
    cu1, cu2 = cu[0], cu[1]
    cu_cu = dist3d(cu1, cu2)
    if not (2.8 <= cu_cu <= 5.5):
        return None
    his_ne2 = [a for a in atoms if a["resn"] == "HIS" and a["atom"] == "NE2"]
    coord = []
    for a in his_ne2:
        if min(dist3d(a, cu1), dist3d(a, cu2)) <= 3.0:
            coord.append(a)
    coord.sort(key=lambda a: int(a["seq"]))
    if len(coord) < 6:
        return None
    plddts = [a["bfactor"] for a in coord[:6]]
    min_plddt = min(plddts)
    return {
        "accession": acc,
        "min_plddt": min_plddt,
        "plddts": plddts,
        "passes_plddt70": all(p >= 70 for p in plddts),
        "cu_cu": cu_cu,
    }

# Six-His-but-not-canonical candidates (the 278 group), from the v3 criteria table.
# (Was rebuilt from v2 bioinf_redo filter_results.csv + helix_coverage.tsv.)
HERE = Path(__file__).resolve().parent
has6 = set()
with open(HERE / "canonical_criteria_all_ca.csv") as f:
    for r in csv.DictReader(f):
        if r["canonical"] == "False" and int(r["n_coord_his"]) >= 6:
            has6.add(r["accession"])

cif_dir = Path("/mnt/models")
tasks = []
for p in sorted(cif_dir.glob("*.cif")):
    acc = p.name.split("_taxID_")[0]
    if acc in has6:
        tasks.append((str(p), acc))
print(f"Checking {len(tasks)} NC structures with 6+ His")

results = []
with Pool(processes=16) as pool:
    for r in pool.imap_unordered(check_one, tasks, chunksize=50):
        if r is not None:
            results.append(r)

pass70 = [r for r in results if r["passes_plddt70"]]
fail70 = [r for r in results if not r["passes_plddt70"]]

print(f"Valid results: {len(results)}")
print(f"Pass pLDDT >= 70 at all 6 His: {len(pass70)}")
print(f"Fail pLDDT >= 70: {len(fail70)}")

if fail70:
    print(f"\nFailing structures:")
    for r in sorted(fail70, key=lambda x: x["min_plddt"]):
        pstr = ", ".join(f"{p:.1f}" for p in r["plddts"])
        print(f"  {r[accession]}: min={r[min_plddt]:.1f}, plddts=[{pstr}]")

min_plddts = sorted([r["min_plddt"] for r in results])
print(f"\nMin-pLDDT distribution:")
print(f"  min={min_plddts[0]:.1f}, p5={min_plddts[len(min_plddts)//20]:.1f}, "
      f"p25={min_plddts[len(min_plddts)//4]:.1f}, median={min_plddts[len(min_plddts)//2]:.1f}, "
      f"max={min_plddts[-1]:.1f}")

print(f"\n=== SIX-HIS-BUT-NOT-CANONICAL ({len(results)}) ===")
print(f"pass per-His pLDDT>=70: {len(pass70)}  (would be canonical on geometry; FINAL rule still discards six-His fails)")
print(f"fail per-His pLDDT>=70: {len(fail70)}")
