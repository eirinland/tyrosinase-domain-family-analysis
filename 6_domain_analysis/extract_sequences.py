#!/usr/bin/env python3
"""Extract amino acid sequences from mmCIF files to FASTA.

Reads CIFs from the squashfuse-mounted directory, writes one FASTA file
per chunk for array job processing.
"""

import csv
import os
import re
import sys
from pathlib import Path

BASE = Path(__file__).parent.parent
CIF_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/af3_cifs")
OUTDIR = Path(__file__).parent / "sequences"
CHUNK_SIZE = int(sys.argv[2]) if len(sys.argv) > 2 else 500

THREE_TO_ONE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'MSE': 'M', 'UNK': 'X',
}


def extract_sequence(cif_path):
    """Extract protein sequence from mmCIF _atom_site records."""
    residues = {}
    in_atom_site = False
    columns = []
    
    with open(cif_path) as f:
        for line in f:
            if line.startswith("_atom_site."):
                in_atom_site = True
                columns.append(line.strip().split(".")[1])
                continue
            if in_atom_site and not line.startswith(("_", "#", "loop_")):
                if line.strip() == "":
                    continue
                parts = line.split()
                if len(parts) < len(columns):
                    continue
                row = dict(zip(columns, parts))
                if row.get("group_PDB") != "ATOM":
                    continue
                resname = row.get("label_comp_id", "")
                resseq = int(row.get("label_seq_id", "0"))
                if resseq > 0 and resname in THREE_TO_ONE:
                    residues[resseq] = THREE_TO_ONE[resname]
            elif in_atom_site and (line.startswith("#") or line.startswith("loop_")):
                in_atom_site = False
                columns = []

    if not residues:
        return ""
    seq = "".join(residues[k] for k in sorted(residues))
    return seq


def main():
    OUTDIR.mkdir(exist_ok=True)
    
    # Get all accessions
    canon_file = BASE / "canonical_criteria_all_ca.csv"
    accessions = []
    with open(canon_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            accessions.append(row["accession"])
    
    print(f"Total accessions: {len(accessions)}", file=sys.stderr)
    
    # Find CIF files
    cif_files = {}
    for fn in os.listdir(CIF_DIR):
        if fn.endswith(".cif"):
            acc = fn.split("_taxID_")[0]
            cif_files[acc] = CIF_DIR / fn
    
    print(f"Found {len(cif_files)} CIF files", file=sys.stderr)
    
    # Extract and chunk
    chunk_idx = 0
    seqs_in_chunk = 0
    outf = None
    n_extracted = 0
    
    for acc in accessions:
        cif = cif_files.get(acc)
        if not cif:
            continue
        
        if seqs_in_chunk == 0:
            if outf:
                outf.close()
            outf = open(OUTDIR / f"chunk_{chunk_idx:04d}.fasta", "w")
        
        seq = extract_sequence(str(cif))
        if seq:
            outf.write(f">{acc}\n{seq}\n")
            n_extracted += 1
            seqs_in_chunk += 1
        
        if seqs_in_chunk >= CHUNK_SIZE:
            seqs_in_chunk = 0
            chunk_idx += 1
        
        if n_extracted % 1000 == 0:
            print(f"  [{n_extracted}]", file=sys.stderr)
    
    if outf:
        outf.close()
    if seqs_in_chunk > 0:
        chunk_idx += 1
    
    print(f"\nExtracted {n_extracted} sequences into {chunk_idx} chunks", file=sys.stderr)


if __name__ == "__main__":
    main()
