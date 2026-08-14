"""
PyMOL: load PmTYR (B2ZB02) + one representative for EVERY residue seen at the
hypervariable Val218 position, drawn as a grey conservation-shaded scaffold with
ONLY position 218 picked out in colour. Companion to pymol_variants_active_site.py
but focused on a single position.

Val218 is the least-conserved active-site position (WT Val only 17%; Phe 19% is
actually more common). Each structure below is the cleanest-placed representative
of one Val218 identity (Cu-bound, well-superposed, minimal confounding
substitutions). Ordered by pool frequency, most common first.

Colouring scheme (identical to pymol_variants_active_site.py):
  * Six coordinating His ................... grey60, 50% transparent
  * Other active-site positions, by conservation tier .. grey70/80/90
  * Position 218 (the variable residue) .... PmTYR Val218 heat colour, labelled
    with its residue name (VAL on the reference, PHE/LEU/... on each variant)
  * Cu .................................... copper, 50% transparent
  * Side chains only (backbone N/O/C hidden)

Usage (PyMOL):
  run pymol_val218_active_site.py
  show_only val218_F_A0A2G2XM25     # flip to any structure
  show_only 0_B2ZB02                # PmTYR (Val, the wild type)
Render:
  ray 2400, 1800; png val218.png, dpi=300
"""

import os
from pymol import cmd, util

# Resolve the visualisation dir robustly (PyMOL `run` rewrites __file__).
_CANDIDATES = []
if '__file__' in globals():
    _CANDIDATES.append(os.path.dirname(os.path.abspath(__file__)))
_CANDIDATES.append(os.getcwd())
_CANDIDATES.append('/Users/eirinlandsem/Library/Mobile Documents/com~apple~CloudDocs/'
                   'Proteinkjemi_PhD/Skriving/Thesis/manuscripts:articles/Manuscript_TyrY/'
                   'New_bioinf/bioinf_redo/Super_reference_pipeline/2_canonical_analysis/'
                   'visualisation')
HERE = next((d for d in _CANDIDATES if os.path.isdir(os.path.join(d, 'Structures'))),
            _CANDIDATES[-1])
STRUCT_DIR = os.path.join(HERE, 'Structures')

# B2ZB02 position -> RGB (verbatim from pymol_active_site.py HEATMAP_RESIDUES).
HEATMAP_COLORS = {
    46:  (1.0,        0.95294118, 0.69019610),
    65:  (0.80392158, 0.70588237, 0.85882354),
    68:  (0.72156864, 0.94901961, 0.90196079),
    195: (1.0,        0.71764706, 0.75686275),
    205: (0.73725490, 0.88627451, 0.68235294),
    209: (0.52941176, 0.78039216, 0.64705882),
    218: (1.0,        0.82352941, 0.65098039),
    221: (0.961,      0.678,      0.506),
    227: (0.80392158, 0.70588237, 0.85882354),
    230: (0.55686277, 0.77254903, 0.98823529),
}
POS_LABELS = {46: 'Gly46', 65: 'Phe65', 68: 'Trp68', 195: 'Glu195', 205: 'Asn205',
              209: 'Arg209', 218: 'Val218', 221: 'Ala221', 227: 'Phe227', 230: 'His230'}
HEATMAP_POSITIONS = list(HEATMAP_COLORS.keys())

# Dataset-wide conservation tier per active-site position -> grey shade.
CONSERVATION = {
    227: 'high', 65: 'high', 68: 'high', 230: 'high', 195: 'high',
    221: 'intermediate', 205: 'intermediate',
    46: 'hypervariable', 218: 'hypervariable', 209: 'hypervariable',
}
TIER_GREY = {'high': 'grey70', 'intermediate': 'grey80', 'hypervariable': 'grey90'}

# Six coordinating His scaffold + CuA His2 (thioether Cys partner).
ANCHOR_HIS = [42, 60, 69, 204, 208, 231]
CUA_HIS2 = 60
FOCUS_POS = 218  # the only position coloured on each variant

# Representative per Val218 residue, ordered by pool frequency (most common first).
# (residue, accession) -- WT Val is the B2ZB02 reference, not repeated here.
VAL218_REPS = [
    ('F', 'A0A2G2XM25'), ('L', 'A0ABD1RDH2'), ('P', 'A0ABZ0KMA8'), ('T', 'A0A6J8BNB1'),
    ('I', 'A0ABR4MAT1'), ('G', 'A0A9D4JWL7'), ('A', 'A0A507FGD0'), ('N', 'A0ABW5LHB0'),
    ('S', 'A0A183C149'), ('M', 'A0A914WL47'), ('Q', 'A0A166CKA0'), ('Y', 'A0A9Q8PMF0'),
    ('K', 'W3XL66'),     ('R', 'A0A0W0FIC5'), ('H', 'A0AAV9H3L9'), ('E', 'A0A5N5QBY0'),
    ('D', 'A0A2J5HIG3'), ('C', 'A0A1I7XST7'), ('W', 'A0A8H7NBP5'),
]

REF_KEY = '0_B2ZB02'
# key -> positions to colour (only 218 on every variant).
DEFINING = {REF_KEY: []}
VARIANT_KEYS = []
for _res, _acc in VAL218_REPS:
    _k = f'val218_{_res}_{_acc}'
    DEFINING[_k] = [FOCUS_POS]
    VARIANT_KEYS.append(_k)

MAPPED = {}  # obj -> {position: variant resi}


def setup_scene():
    cmd.bg_color('white')
    cmd.set('valence', 0)
    cmd.set('ray_shadows', 'off')
    cmd.set('specular', 0.2)
    cmd.set('ambient', 0.4)
    cmd.set('label_size', 14)
    cmd.set('label_color', 'black')
    cmd.set('label_font_id', 7)
    cmd.set('label_position', [0, 0, 3])


def _ca_coord(sel):
    m = cmd.get_model(f'({sel}) and name CA')
    return m.atom[0].coord if m.atom else None


def _resn(sel):
    m = cmd.get_model(f'({sel}) and name CA')
    return m.atom[0].resn if m.atom else '?'


def map_position(obj, ref_resi, tol=3.0):
    """B2ZB02 position -> obj residue whose Ca is nearest (after cealign)."""
    ref_xyz = _ca_coord(f'pmtyr and resi {ref_resi}')
    if ref_xyz is None:
        return None
    near = f'{obj} and name CA within {tol} of (pmtyr and resi {ref_resi} and name CA)'
    best, best_d = None, 1e9
    for at in cmd.get_model(near).atom:
        d = sum((a - b) ** 2 for a, b in zip(at.coord, ref_xyz)) ** 0.5
        if d < best_d:
            best_d, best = d, at.resi
    return best


def _sticks(obj, resi, radius):
    cmd.show('sticks', f'{obj} and resi {resi} and not name C+N+O')
    cmd.set('stick_radius', radius, f'{obj} and resi {resi}')


def color_structure(obj, key, is_ref=False):
    color_pos = set(HEATMAP_POSITIONS) if is_ref else set(DEFINING[key])

    cmd.show('cartoon', obj)
    cmd.color('grey90', obj)
    cmd.set('cartoon_transparency', 0.75, obj)

    cu = f'{obj} and name CU'
    if cmd.count_atoms(cu) > 0:
        cmd.show('spheres', cu)
        util.cnc(cu)
        cmd.set('sphere_scale', 0.25, cu)
        cmd.set('sphere_transparency', 0.5, cu)

    his_resis = []
    for a in ANCHOR_HIS:
        rv = a if is_ref else map_position(obj, a)
        if rv is not None:
            _sticks(obj, rv, 0.18)
            cmd.color('grey60', f'{obj} and resi {rv} and not name C+N+O')
            cmd.set('stick_transparency', 0.5, f'{obj} and resi {rv}')
            his_resis.append(rv)

    if his_resis and cmd.count_atoms(cu) > 0:
        for rv in his_resis:
            ne2 = f'{obj} and resi {rv} and name NE2'
            cu_lig = f'({cu}) within 3.5 of ({ne2})'
            if cmd.count_atoms(cu_lig) > 0:
                cmd.bond(ne2, cu_lig)
                cmd.set_bond('stick_radius', 0.15, ne2, cu_lig)
                cmd.set_bond('stick_transparency', 0.5, ne2, cu_lig)
        cmd.show('sticks', cu)

    mapped = {}
    for pos in HEATMAP_POSITIONS:
        rv = pos if is_ref else map_position(obj, pos)
        if rv is None:
            continue
        mapped[pos] = rv
        if pos in color_pos:
            cname = f'heat_{pos}'
            cmd.set_color(cname, list(HEATMAP_COLORS[pos]))
            _sticks(obj, rv, 0.28)
            cmd.color(cname, f'{obj} and resi {rv} and not name C+N+O')
            # Reference = full position legend; variants = the 218 residue identity.
            lab = POS_LABELS[pos] if is_ref else _resn(f'{obj} and resi {rv}')
            cmd.label(f'{obj} and resi {rv} and name CA', f'"{lab}"')
        else:
            _sticks(obj, rv, 0.20)
            cmd.color(TIER_GREY[CONSERVATION[pos]], f'{obj} and resi {rv} and not name C+N+O')

    MAPPED[obj] = mapped
    print(f'{key}: mapped {len(mapped)}/10 positions, coloured {sorted(color_pos & mapped.keys())}.')


def show_only(key):
    """Disable every other object, enable only `key`, zoom its active site."""
    obj = 'pmtyr' if key == REF_KEY else key
    for k in cmd.get_names('objects'):
        cmd.disable(k)
    cmd.enable(obj)
    resis = list(MAPPED.get(obj, {}).values()) + \
            [v for v in (map_position(obj, a) for a in ANCHOR_HIS) if v]
    if resis:
        cmd.zoom(f'{obj} and resi ' + '+'.join(map(str, resis)), buffer=5)
    cmd.deselect()


def build():
    setup_scene()
    cmd.load(os.path.join(STRUCT_DIR, f'{REF_KEY}.cif'), 'pmtyr')
    color_structure('pmtyr', REF_KEY, is_ref=True)

    for key in VARIANT_KEYS:
        path = os.path.join(STRUCT_DIR, f'{key}.cif')
        if not os.path.exists(path):
            print(f'MISSING: {path}')
            continue
        cmd.load(path, key)
        try:
            cmd.cealign('pmtyr', key, object=f'aln_{key}')
            cmd.delete(f'aln_{key}')
        except Exception as e:
            print(f'{key}: cealign failed ({e}) -- positions may be mis-mapped.')
        color_structure(key, key)

    cmd.deselect()
    show_only(VARIANT_KEYS[0])
    print('\nReady. Val218 representatives (pool-frequency order):')
    print('  ' + '   '.join(VARIANT_KEYS))
    print('Flip with:  show_only <name>   e.g. show_only val218_F_A0A2G2XM25')
    print('Reference (WT Val):  show_only 0_B2ZB02')
    print('Render: ray 2400, 1800; png val218.png, dpi=300')


cmd.extend('show_only', show_only)
build()
