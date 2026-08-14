"""Trim AF3 structures to PPO core domain(s) identified by Chainsaw+Foldseek.
Keeps all HETATM/ligand lines (Cu, etc.) regardless of trimming range."""
import csv, os, sys, argparse


def parse_ppo_range(ppo_range_str):
    if not ppo_range_str:
        return None
    parts = ppo_range_str.split("-")
    return int(parts[0]), int(parts[1])


def trim_cif(cif_path, out_path, res_start, res_end):
    with open(cif_path) as f:
        lines = f.readlines()

    in_atom_site = False
    col_names = []
    col_idx = {}
    header_lines = []
    atom_lines = []
    pre_lines = []
    post_lines = []
    past_atom_site = False

    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("_atom_site."):
            if not in_atom_site:
                in_atom_site = True
            col_name = line.strip().split(".")[1].split()[0]
            col_names.append(col_name)
            col_idx[col_name] = len(col_names) - 1
            header_lines.append(line)
            i += 1
            continue

        if in_atom_site and not line.startswith("_") and not line.startswith("#") and line.strip():
            tokens = line.split()
            seq_col = col_idx.get("label_seq_id")
            group_col = col_idx.get("group_PDB")
            if seq_col is not None and seq_col < len(tokens):
                seq_id = tokens[seq_col]
                if seq_id == ".":
                    # HETATM (ligands, metals) — always keep
                    atom_lines.append(line)
                elif res_start <= int(seq_id) <= res_end:
                    atom_lines.append(line)
            i += 1
            continue

        if in_atom_site and (line.startswith("#") or line.strip() == "" or line.startswith("_")):
            in_atom_site = False
            past_atom_site = True
            post_lines.append(line)
            i += 1
            continue

        if not past_atom_site and not in_atom_site:
            pre_lines.append(line)
        else:
            post_lines.append(line)
        i += 1

    with open(out_path, "w") as f:
        for l in pre_lines:
            f.write(l)
        for l in header_lines:
            f.write(l)
        for l in atom_lines:
            f.write(l)
        for l in post_lines:
            f.write(l)

    return len(atom_lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--domain-csv", required=True)
    parser.add_argument("--cif-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--accessions", default=None)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    domains = {}
    with open(args.domain_csv) as f:
        for r in csv.DictReader(f):
            if r["ppo_range"]:
                domains[r["accession"]] = r

    if args.accessions:
        acc_list = args.accessions.split(",")
    else:
        acc_list = sorted(domains.keys())

    cif_files = {}
    for fname in os.listdir(args.cif_dir):
        if fname.endswith(".cif"):
            acc = fname.split("_taxID_")[0]
            cif_files[acc] = os.path.join(args.cif_dir, fname)

    trimmed = 0
    skipped = 0
    for acc in acc_list:
        if acc not in domains:
            print("SKIP %s: no PPO domain" % acc)
            skipped += 1
            continue
        if acc not in cif_files:
            print("SKIP %s: no CIF file" % acc)
            skipped += 1
            continue

        rng = parse_ppo_range(domains[acc]["ppo_range"])
        if not rng:
            skipped += 1
            continue

        cif_name = os.path.basename(cif_files[acc])
        out_path = os.path.join(args.out_dir, cif_name.replace(".cif", "_trimmed.cif"))
        n_atoms = trim_cif(cif_files[acc], out_path, rng[0], rng[1])
        print("TRIM %s: %d-%d (%d atom lines)" % (acc, rng[0], rng[1], n_atoms))
        trimmed += 1

    print("\nTrimmed: %d, Skipped: %d" % (trimmed, skipped))


if __name__ == "__main__":
    main()
