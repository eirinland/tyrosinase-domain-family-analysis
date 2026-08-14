"""Extract chain A from PDB crystal structures and trim AF3 fill-ins.
Outputs clean single-chain PDB files for Foldseek database.
No external dependencies — parses mmCIF atom_site records directly.
"""
import os, sys

FOLDSEEK_DIR = os.path.dirname(os.path.abspath(__file__))
REFS_DIR = "/cluster/home/eirinlandsem/Super_reference_pipeline/refs"

refs = [
    (f"{FOLDSEEK_DIR}/8BBR.cif", "ref_8BBR_Vspinosum", "A", None, None),
    (f"{FOLDSEEK_DIR}/2Y9W.cif", "ref_2Y9W_Abisporus", "A", None, None),
    (f"{FOLDSEEK_DIR}/1BT3.cif", "ref_1BT3_Ibatatas", "A", None, None),
    (f"{FOLDSEEK_DIR}/5CE9.cif", "ref_5CE9_Jregia", "A", None, None),
    (f"{FOLDSEEK_DIR}/5M8L.cif", "ref_5M8L_human", "A", None, None),
    (f"{FOLDSEEK_DIR}/1JS8.cif", "ref_1JS8_squid", "A", None, None),
    (f"{REFS_DIR}/I3D139_taxID_859350_model.cif", "ref_I3D139_archaea", None, 6, 310),
    (f"{REFS_DIR}/A0A9N8ELP9_taxID_568900_model.cif", "ref_A0A9N8ELP9_oomycota", None, 347, 758),
]


def parse_cif_columns(path):
    """Find the column order in the _atom_site loop and return atom lines."""
    columns = []
    in_atom_site = False
    in_data = False
    atoms = []

    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("_atom_site."):
                in_atom_site = True
                col_name = line.split(".")[1].split()[0]
                columns.append(col_name)
                continue
            if in_atom_site and not line.startswith("_atom_site."):
                in_atom_site = False
                in_data = True
            if in_data:
                if line.startswith("#") or line.startswith("_") or line.startswith("loop_"):
                    break
                if line.strip():
                    atoms.append(line)

    return columns, atoms


def extract_chain(path, chain_id, start_res, end_res):
    """Extract atoms for a given chain (and optional residue range) as PDB lines."""
    columns, atoms = parse_cif_columns(path)

    col_idx = {name: i for i, name in enumerate(columns)}
    i_group = col_idx.get("group_PDB")
    i_type = col_idx.get("type_symbol")
    i_name = col_idx.get("label_atom_id")
    i_alt = col_idx.get("label_alt_id")
    i_resn = col_idx.get("label_comp_id")
    i_chain = col_idx.get("auth_asym_id", col_idx.get("label_asym_id"))
    i_seq = col_idx.get("auth_seq_id", col_idx.get("label_seq_id"))
    i_x = col_idx.get("Cartn_x")
    i_y = col_idx.get("Cartn_y")
    i_z = col_idx.get("Cartn_z")
    i_occ = col_idx.get("occupancy")
    i_bfac = col_idx.get("B_iso_or_equiv")

    pdb_lines = []
    serial = 0
    for atom_line in atoms:
        fields = atom_line.split()
        group = fields[i_group]
        if group not in ("ATOM", "HETATM"):
            continue

        ch = fields[i_chain]
        if chain_id and ch != chain_id:
            continue

        seq = int(fields[i_seq])
        if start_res and seq < start_res:
            continue
        if end_res and seq > end_res:
            continue

        alt = fields[i_alt] if fields[i_alt] != "." else " "
        if alt not in (" ", "A"):
            continue

        serial += 1
        atom_name = fields[i_name]
        resn = fields[i_resn]
        x = float(fields[i_x])
        y = float(fields[i_y])
        z = float(fields[i_z])
        occ = float(fields[i_occ]) if i_occ else 1.0
        bfac = float(fields[i_bfac]) if i_bfac else 0.0
        elem = fields[i_type] if i_type else atom_name[0]

        # PDB format atom name: 4 chars, left-justified if 4 chars, else starts at col 2
        if len(atom_name) < 4:
            atom_name_fmt = f" {atom_name:<3s}"
        else:
            atom_name_fmt = f"{atom_name:<4s}"

        pdb_lines.append(
            f"{group:<6s}{serial:5d} {atom_name_fmt}{alt if alt != ' ' else ' '}"
            f"{resn:>3s} {'A':1s}{seq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfac:6.2f}"
            f"          {elem:>2s}  "
        )

    return pdb_lines


for input_path, out_name, chain_id, start, end in refs:
    print(f"Processing {out_name}...")
    pdb_lines = extract_chain(input_path, chain_id, start, end)
    out_path = os.path.join(FOLDSEEK_DIR, f"{out_name}.pdb")
    with open(out_path, "w") as f:
        for line in pdb_lines:
            f.write(line + "\n")
        f.write("END\n")
    n_ca = sum(1 for l in pdb_lines if " CA " in l and l.startswith("ATOM"))
    print(f"  -> {out_path} ({n_ca} residues)")

print("\nDone. PmTYR must be added from user's crystal structure PDB.")
