"""Run AllMetal3D on AF2 structures and save metal prediction PDBs."""
import glob, os, subprocess, sys
import gemmi

EVAL_DIR = os.path.dirname(os.path.abspath(__file__))

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

for cif in sorted(glob.glob(os.path.join(EVAL_DIR, "*_AF2_*.cif"))):
    basename = os.path.basename(cif).replace(".cif", "")
    print(f"Processing {basename}...", flush=True)

    pdb_path = os.path.join(EVAL_DIR, f"{basename}.pdb")
    out_dir = os.path.join(EVAL_DIR, f"{basename}_metal3d")
    os.makedirs(out_dir, exist_ok=True)

    try:
        cif_to_pdb(cif, pdb_path)
    except Exception as e:
        print(f"  CIF error: {e}", flush=True)
        continue

    result = subprocess.run(
        ["allmetal3d", "-i", pdb_path, "-o", out_dir,
         "--models", "allmetal3d", "-m", "fast", "-p", "0.1"],
        capture_output=True, text=True, timeout=1200
    )
    if result.returncode != 0:
        print(f"  Metal3D error: {result.stderr[-200:]}", flush=True)
    else:
        files = os.listdir(out_dir)
        print(f"  Done: {files}", flush=True)

print("\nAll done.", flush=True)
