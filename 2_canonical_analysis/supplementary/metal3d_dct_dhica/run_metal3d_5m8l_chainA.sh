#!/usr/bin/env bash
#SBATCH --job-name=m3d_5m8l
#SBATCH --account=nn1003k
#SBATCH --partition=accel
#SBATCH --gpus=1
#SBATCH --mem-per-gpu=120G
#SBATCH --time=02:00:00
#SBATCH --output=/cluster/work/projects/nn1003k/eirin/bioinf/bioinf_redo/2_canonical_analysis/supplementary/metal3d_dct_dhica/metal3d_5m8l_chainA_%j.log

set -euo pipefail

module load NRIS/GPU

SUBMITDIR=/cluster/work/projects/nn1003k/eirin/bioinf/bioinf_redo/2_canonical_analysis/supplementary/metal3d_dct_dhica
SIF=/cluster/projects/nn1003k/prog/allmetal3d/allmetal3d_gpu.sif
CIFPATH=$SUBMITDIR/pdb_5m8l/5m8l.cif

echo "Started: $(date)"

apptainer exec --nv --cleanenv \
    --bind "$SUBMITDIR:$SUBMITDIR" \
    "$SIF" \
    python3 -c "
import gemmi, subprocess, tempfile, os, glob
from math import sqrt

st = gemmi.read_structure('$CIFPATH')
model = st[0]

# Keep only chain A
chains_to_remove = [ch.name for ch in model if ch.name != 'A']
for name in chains_to_remove:
    model.remove_chain(name)
print(f'Kept chain A only: {sum(1 for ch in model for _ in ch)} residues')

# Find Zn positions before stripping
zn_positions = []
metals = {'CU', 'ZN', 'FE', 'MN', 'CO', 'NI', 'MG', 'CA'}
for chain in model:
    for res in chain:
        if res.name == 'ZN':
            for atom in res:
                if atom.element.name == 'Zn':
                    zn_positions.append((chain.name, atom.pos.x, atom.pos.y, atom.pos.z, atom.b_iso))
                    print(f'  Zn: chain {chain.name} ({atom.pos.x:.2f}, {atom.pos.y:.2f}, {atom.pos.z:.2f}) B={atom.b_iso:.1f}')

# Strip metals
for chain in model:
    to_remove = [i for i, res in enumerate(chain) if res.name in metals]
    for i in reversed(to_remove):
        del chain[i]

# Also remove waters and other hetero
for chain in model:
    to_remove = [i for i, res in enumerate(chain) if res.name == 'HOH']
    for i in reversed(to_remove):
        del chain[i]

print(f'After stripping: {sum(1 for ch in model for _ in ch)} residues')

with tempfile.TemporaryDirectory() as tmpdir:
    pdb_path = os.path.join(tmpdir, '5m8l_chainA.pdb')
    out_dir = os.path.join(tmpdir, 'output')
    os.makedirs(out_dir)
    st.write_pdb(pdb_path)

    # Check file size
    sz = os.path.getsize(pdb_path)
    print(f'PDB file size: {sz} bytes')

    result = subprocess.run(
        ['allmetal3d', '-i', pdb_path, '-o', out_dir,
         '--models', 'allmetal3d', '-m', 'fast', '-p', '0.1'],
        capture_output=True, text=True, timeout=3600
    )
    print(f'Metal3D exit code: {result.returncode}')
    if result.stderr:
        lines = result.stderr.strip().split('\n')
        print(f'stderr (last 10 lines):')
        for l in lines[-10:]:
            print(f'  {l}')

    metal_files = glob.glob(os.path.join(out_dir, '*_metals.pdb'))
    if not metal_files:
        metal_files = glob.glob(os.path.join(out_dir, '*.pdb'))
    print(f'Output files: {[os.path.basename(f) for f in metal_files]}')

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

    print(f'\nTotal predicted metal sites: {len(sites)}')
    # Count by element
    from collections import Counter
    elem_counts = Counter(e for _, _, _, _, e in sites)
    print(f'By element: {dict(elem_counts)}')

    # Match to Zn positions
    for ch, zx, zy, zz, b in zn_positions:
        print(f'\nZn in chain {ch} ({zx:.2f}, {zy:.2f}, {zz:.2f}):')
        if sites:
            dists = [(sqrt((zx-sx)**2+(zy-sy)**2+(zz-sz)**2), sp, se) for sx, sy, sz, sp, se in sites]
            dists.sort()
            print('  5 nearest predictions:')
            for d, p, e in dists[:5]:
                print(f'    {e:4s}  dist={d:.2f} A  prob={p:.3f}')
"

echo "Done: $(date)"
