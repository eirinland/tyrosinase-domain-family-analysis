import csv, math, glob, os, sys
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

def check_one(cif_path):
    acc = Path(cif_path).name.split("_taxID_")[0]
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
    
    first_his = coord[0]
    last_his = coord[5]
    first_res = seqs[0]
    last_res = seqs[-1]
    
    return {
        "accession": acc,
        "total_len": last_res - first_res + 1,
        "first_res": first_res,
        "first_his": first_his,
        "last_his": last_his,
        "last_res": last_res,
        "nterm_before_his1": first_his - first_res,
        "cterm_after_his6": last_res - last_his,
        "his_positions": ",".join(str(x) for x in coord[:6]),
    }

def load_accessions(path):
    accs = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and line.lower() != "accession":
                accs.add(line.split(",")[0])
    return accs

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--cifs", required=True)
    parser.add_argument("--canonical", required=True,
                        help="final_pools/canonical_accessions.csv")
    parser.add_argument("--noncanonical", required=True,
                        help="final_pools/noncanonical_accessions.csv")
    parser.add_argument("--criteria", required=True,
                        help="canonical_criteria_all_ca.csv (six-His-but-not-canonical set)")
    parser.add_argument("--workers", type=int, default=16)
    args = parser.parse_args()

    canonical = load_accessions(args.canonical)
    promoted = load_accessions(args.noncanonical)   # non-canonical pool (valid core, <6 His)

    nc_6his = set()
    with open(args.criteria) as f:
        for r in csv.DictReader(f):
            if r["canonical"] == "False" and int(r["n_coord_his"]) >= 6:
                nc_6his.add(r["accession"])

    all_targets = canonical | promoted | nc_6his
    print(f"Canonical: {len(canonical)}, Non-canonical: {len(promoted)}, Six-His-fail: {len(nc_6his)}", flush=True)

    cif_dir = Path(args.cifs)
    acc_to_cif = {}
    for p in sorted(cif_dir.glob("*.cif")):
        acc = p.name.split("_taxID_")[0]
        if acc in all_targets:
            acc_to_cif[acc] = str(p)
    print(f"CIF files matched: {len(acc_to_cif)}", flush=True)

    cif_list = sorted(acc_to_cif.values())
    results = []
    done = 0
    with Pool(processes=args.workers) as pool:
        for result in pool.imap_unordered(check_one, cif_list, chunksize=200):
            done += 1
            if done % 5000 == 0:
                print(f"  {done}/{len(cif_list)}...", flush=True)
            if result is not None:
                results.append(result)
    print(f"Valid results: {len(results)}", flush=True)

    # Classify each result
    can_nterm = []
    promo_nterm = []
    nc_nterm = []
    for r in results:
        acc = r["accession"]
        nterm = r["nterm_before_his1"]
        if acc in canonical:
            can_nterm.append(r)
        elif acc in promoted:
            promo_nterm.append(r)
        elif acc in nc_6his:
            nc_nterm.append(r)

    print(f"\n=== EXISTING CANONICAL ({len(can_nterm)}) ===")
    nterms = sorted([r["nterm_before_his1"] for r in can_nterm])
    print(f"N-term before His1: min={nterms[0]}, p5={nterms[len(nterms)//20]}, "
          f"p25={nterms[len(nterms)//4]}, median={nterms[len(nterms)//2]}, "
          f"p75={nterms[3*len(nterms)//4]}, max={nterms[-1]}")
    bins = Counter()
    for n in nterms:
        if n < 10: bins["<10"] += 1
        elif n < 20: bins["10-19"] += 1
        elif n < 30: bins["20-29"] += 1
        elif n < 50: bins["30-49"] += 1
        else: bins["50+"] += 1
    for b in ["<10", "10-19", "20-29", "30-49", "50+"]:
        print(f"  {b}: {bins.get(b, 0)}")

    print(f"\n=== NON-CANONICAL POOL ({len(promo_nterm)}) ===")
    if promo_nterm:
        nterms = sorted([r["nterm_before_his1"] for r in promo_nterm])
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
        # Show the worst ones
        worst = sorted(promo_nterm, key=lambda r: r["nterm_before_his1"])[:10]
        print("Shortest N-term:")
        for w in worst:
            print(f"  {w[accession]}: nterm={w[nterm_before_his1]}, len={w[total_len]}, his={w[his_positions]}")

    print(f"\n=== SIX-HIS-BUT-NOT-CANONICAL ({len(nc_nterm)}) ===")
    if nc_nterm:
        nterms = sorted([r["nterm_before_his1"] for r in nc_nterm])
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
        worst = sorted(nc_nterm, key=lambda r: r["nterm_before_his1"])[:10]
        print("Shortest N-term:")
        for w in worst:
            print(f"  {w[accession]}: nterm={w[nterm_before_his1]}, len={w[total_len]}, his={w[his_positions]}")

    # Also check: what fraction of EXISTING canonical have <10 residues before His1?
    print(f"\n=== FRAGMENTATION SUMMARY ===")
    can_short = [r for r in can_nterm if r["nterm_before_his1"] < 10]
    promo_short = [r for r in promo_nterm if r["nterm_before_his1"] < 10]
    nc_short = [r for r in nc_nterm if r["nterm_before_his1"] < 10]
    print(f"Existing canonical <10 before His1: {len(can_short)}/{len(can_nterm)} ({100*len(can_short)/max(1,len(can_nterm)):.1f}%)")
    print(f"Non-canonical <10 before His1: {len(promo_short)}/{len(promo_nterm)} ({100*len(promo_short)/max(1,len(promo_nterm)):.1f}%)")
    print(f"Six-His-fail <10 before His1: {len(nc_short)}/{len(nc_nterm)} ({100*len(nc_short)/max(1,len(nc_nterm)):.1f}%)")

if __name__ == "__main__":
    main()
