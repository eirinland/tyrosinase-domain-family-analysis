import csv, math, glob
from pathlib import Path
from multiprocessing import Pool
from collections import Counter

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
            coord.append(int(a["seq"]))
    coord.sort()
    if len(coord) < 6:
        return None
    seqs = sorted(set(int(a["seq"]) for a in atoms if a["seq"].isdigit()))
    if not seqs:
        return None
    return {
        "accession": acc,
        "total_len": seqs[-1] - seqs[0] + 1,
        "first_his": coord[0],
        "first_res": seqs[0],
        "nterm_before_his1": coord[0] - seqs[0],
        "his_positions": coord[:6],
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
print(f"NC with 6 His to check: {len(has6)}, CIFs found: {len(tasks)}")

results = []
with Pool(processes=16) as pool:
    for r in pool.imap_unordered(check_one, tasks, chunksize=50):
        if r is not None:
            results.append(r)

nterms = sorted([r["nterm_before_his1"] for r in results])
print(f"Valid: {len(results)}")
if nterms:
    print(f"N-term before His1: min={nterms[0]}, median={nterms[len(nterms)//2]}, max={nterms[-1]}")
    bins = Counter()
    for n in nterms:
        if n < 10: bins["<10"] += 1
        elif n < 20: bins["10-19"] += 1
        elif n < 30: bins["20-29"] += 1
        elif n < 50: bins["30-49"] += 1
        else: bins["50+"] += 1
    for b in ["<10", "10-19", "20-29", "30-49", "50+"]:
        print(f"  {b}: {bins.get(b, 0)}")
    
    short = [r for r in results if r["nterm_before_his1"] < 10]
    print(f"\nShort N-term (<10 residues): {len(short)}")
    for s in sorted(short, key=lambda r: r["nterm_before_his1"]):
        print(f"  {s[accession]}: nterm={s[nterm_before_his1]}, len={s[total_len]}, his={s[his_positions]}")
