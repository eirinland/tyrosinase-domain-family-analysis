"""
PyMOL: colour PmTYR (B2ZB02) active-site residues with the SAME colours the
positions carry in crystal_annotation_heatmap.py (the column-header colours).

Each residue here uses the exact RGB from COL_COLORS in
crystal_annotation_heatmap.py, so a residue in PyMOL matches its column in the
figure 1:1. Numbering is the B2ZB02 AF3 model (PMTYR_POSITIONS in
extract_position_vectors.py).

Usage (PyMOL):
  load /path/to/B2ZB02_taxID_1404_model.cif, pmtyr
  run pymol_active_site.py
  # or just `run pymol_active_site.py` and it loads PMTYR_PATH below.
Then: ray 2400, 1800; png active_site.png, dpi=300
"""

from pymol import cmd

# Local AF3 model of PmTYR (adjust if needed).
import os
_HERE = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
PMTYR_PATH = os.path.join(_HERE, '..', '..', '1_filtering', 'B2ZB02_taxID_1404_model.cif')

# Heatmap position -> (B2ZB02 resi, COL_COLORS rgb from crystal_annotation_heatmap.py).
# Colours are copied verbatim so PyMOL == figure column headers.
HEATMAP_RESIDUES = {
    'Gly46':  (46,  (1.0,        0.95294118, 0.69019610)),
    'Phe65':  (65,  (0.80392158, 0.70588237, 0.85882354)),
    'Trp68':  (68,  (0.72156864, 0.94901961, 0.90196079)),
    'Glu195': (195, (1.0,        0.71764706, 0.75686275)),
    'Asn205': (205, (0.73725490, 0.88627451, 0.68235294)),
    'Arg209': (209, (0.52941176, 0.78039216, 0.64705882)),
    'Val218': (218, (1.0,        0.82352941, 0.65098039)),
    'Ala221': (221, (0.961,      0.678,      0.506)),
    'Phe227': (227, (0.80392158, 0.70588237, 0.85882354)),
    'His230': (230, (0.55686277, 0.77254903, 0.98823529)),
}
# 11th heatmap column ("Cys"/thioether). PmTYR has no thioether Cys, so this is
# only applied if a Cys SG sits within 3.5 A of the CuA His2 (resi 60) ring.
THIOETHER_COLOR = (1.0, 0.95294118, 0.69019610)

# Six coordinating His (superposition anchors) — neutral grey scaffold, not a heatmap column.
ANCHOR_HIS = [42, 60, 69, 204, 208, 231]


def setup_scene():
    cmd.bg_color('white')
    cmd.set('ray_shadows', 'off')
    cmd.set('specular', 0.2)
    cmd.set('ambient', 0.4)
    cmd.set('label_size', 14)
    cmd.set('label_color', 'black')
    cmd.set('label_font_id', 7)
    cmd.set('label_position', [0, 0, 3])


def visualise(obj='pmtyr'):
    setup_scene()

    cmd.show('cartoon', obj)
    cmd.color('grey80', obj)
    cmd.set('cartoon_transparency', 0.7, obj)

    # Cu atoms — orange spheres
    cmd.select('cu_atoms', f'{obj} and name CU')
    cmd.show('spheres', 'cu_atoms')
    cmd.color('orange', 'cu_atoms')
    cmd.set('sphere_scale', 0.2, 'cu_atoms')

    # Six coordinating His — grey scaffold sticks
    for resi in ANCHOR_HIS:
        sel = f'anchor_his_{resi}'
        cmd.select(sel, f'{obj} and resi {resi} and not name C+N+O')
        cmd.show('sticks', sel)
        cmd.color('grey50', sel)
        cmd.set('stick_radius', 0.15, sel)

    # Heatmap positions — coloured sticks with the figure's column colours
    for name, (resi, rgb) in HEATMAP_RESIDUES.items():
        cname = f'heat_{resi}'
        cmd.set_color(cname, list(rgb))
        sel = f'pos_{resi}'
        cmd.select(sel, f'{obj} and resi {resi}')
        cmd.show('sticks', sel)
        cmd.color(cname, sel)
        cmd.set('stick_radius', 0.2, sel)
        cmd.label(f'{obj} and resi {resi} and name CA', f'"{name}"')

    # Thioether Cys (heatmap "Cys" column) — only if one is geometrically present
    cmd.set_color('heat_thioether', list(THIOETHER_COLOR))
    cmd.select('thioether_cys', f'{obj} and resn CYS and name SG '
                                f'within 3.5 of ({obj} and resi 60 and name ND1+NE2+CG+CD2+CE1)')
    if cmd.count_atoms('thioether_cys') > 0:
        cmd.select('thioether_cys', 'byres thioether_cys')
        cmd.show('sticks', 'thioether_cys')
        cmd.color('heat_thioether', 'thioether_cys')
        cmd.set('stick_radius', 0.2, 'thioether_cys')
    else:
        print('No thioether Cys near CuA His2 (expected for PmTYR) — Cys column not drawn.')

    cmd.deselect()

    all_resi = ANCHOR_HIS + [r for r, _ in HEATMAP_RESIDUES.values()]
    cmd.select('active_site', f'{obj} and resi ' + '+'.join(map(str, all_resi)))
    cmd.zoom('active_site', buffer=5)
    cmd.delete('active_site')

    print('Active-site visualisation ready (colours match crystal_annotation_heatmap.py).')
    print('Render: ray 2400, 1800; png active_site.png, dpi=300')


if cmd.get_names('objects'):
    obj = cmd.get_names('objects')[0]
    print(f'Using loaded object: {obj}')
    visualise(obj)
else:
    cmd.load(PMTYR_PATH, 'pmtyr')
    visualise('pmtyr')
