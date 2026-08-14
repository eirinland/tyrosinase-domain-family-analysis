#!/usr/bin/env bash
#SBATCH --job-name=m3d_5m8l
#SBATCH --account=nn1003k
#SBATCH --partition=accel
#SBATCH --gpus=1
#SBATCH --mem-per-gpu=120G
#SBATCH --time=02:00:00
#SBATCH --output=/cluster/work/projects/nn1003k/eirin/bioinf/bioinf_redo/2_canonical_analysis/supplementary/metal3d_dct_dhica/metal3d_5m8l_%j.log

set -euo pipefail

module load NRIS/GPU

SUBMITDIR=/cluster/work/projects/nn1003k/eirin/bioinf/bioinf_redo/2_canonical_analysis/supplementary/metal3d_dct_dhica
SIF=/cluster/projects/nn1003k/prog/allmetal3d/allmetal3d_gpu.sif
CIFDIR=$SUBMITDIR/pdb_5m8l

echo "Started: $(date)"

apptainer exec --nv --cleanenv \
    --bind "$CIFDIR:$CIFDIR:ro" \
    --bind "$SUBMITDIR:$SUBMITDIR" \
    "$SIF" \
    python3 -c "
import gemmi, subprocess, tempfile, os, glob

cif_path = '$CIFDIR/5m8l.cif'
doc = gemmi.cif.read(cif_path)
block = doc.sole_block()
st = gemmi.make_structure_from_block(block)

# Strip metals for Metal3D (expects apo)
metals = {'CU', 'ZN', 'FE', 'MN', 'CO', 'NI', 'MG', 'CA'}
cu_positions = []
for model in st:
    for chain in model:
        for res in chain:
            if res.name == 'CU':
                for atom in res:
                    if atom.element.name == 'Cu':
                        cu_positions.append((chain.name, atom.pos.x, atom.pos.y, atom.pos.z, atom.b_iso))
        to_remove = [i for i, res in enumerate(chain) if res.name in metals]
        for i in reversed(to_remove):
            del chain[i]

print(f'Found {len(cu_positions)} Cu atoms in PDB structure')
for ch, x, y, z, b in cu_positions:
    print(f'  Chain {ch}: ({x:.2f}, {y:.2f}, {z:.2f}) B={b:.1f}')

with tempfile.TemporaryDirectory() as tmpdir:
    pdb_path = os.path.join(tmpdir, '5m8l.pdb')
    out_dir = os.path.join(tmpdir, 'output')
    os.makedirs(out_dir)
    st.write_pdb(pdb_path)

    result = subprocess.run(
        ['allmetal3d', '-i', pdb_path, '-o', out_dir,
         '--models', 'allmetal3d', '-m', 'fast', '-p', '0.1'],
        capture_output=True, text=True, timeout=3600
    )
    print(f'Metal3D exit code: {result.returncode}')
    if result.stderr:
        print(f'stderr (last 500): {result.stderr[-500:]}')

    # Parse predictions
    metal_files = glob.glob(os.path.join(out_dir, '*_metals.pdb'))
    if not metal_files:
        metal_files = glob.glob(os.path.join(out_dir, '*.pdb'))
    print(f'Output files: {metal_files}')

    sites = []
    for pf in metal_files:
        with open(pf) as f:
            for line in f:
                if line.startswith('HETATM') or line.startswith('ATOM'):
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        prob = float(line[54:60])
                        elem = line[76:80].strip() if len(line) > 78 else ''
                        sites.append((x, y, z, prob, elem))
                    except ValueError:
                        pass
    print(f'Total predicted metal sites: {len(sites)}')

    # Report nearest prediction for each Cu
    from math import sqrt
    for ch, cx, cy, cz, b in cu_positions:
        print(f'\nCu in chain {ch} ({cx:.2f}, {cy:.2f}, {cz:.2f}):')
        if sites:
            dists = [(sqrt((cx-sx)**2+(cy-sy)**2+(cz-sz)**2), sp, se) for sx, sy, sz, sp, se in sites]
            dists.sort()
            print('  Nearest predictions:')
            for d, p, e in dists[:5]:
                print(f'    {e:4s}  dist={d:.2f} A  prob={p:.3f}')
"

echo "Done: $(date)"
