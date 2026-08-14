"""
PyMOL: load PmTYR (B2ZB02) + the 22 representative active-site variants and draw
each active site as a grey, conservation-shaded scaffold with ONLY the residue(s)
discussed in the manuscript text picked out in colour.

Colouring scheme (per structure):
  * Six coordinating His ................... grey60
  * Active-site positions, by dataset-wide conservation tier:
        highly conserved  (Phe227,Phe65,Trp68,His230,Glu195; 77-99%) -> grey70
        intermediate      (Ala221, Asn205) ............................ grey80
        hypervariable     (Gly46, Val218, Arg209) ..................... grey90
  * The text-mentioned residue(s) for that group .... coloured with the SAME
    colour the equivalent position carries in PmTYR (pymol_active_site.py /
    crystal_annotation_heatmap.py column colours).
  * Cu .................................... normal copper colour, 50% transparent
  * Backbone N/O/C atoms are hidden on every stick (side chains only).

B2ZB02 (0_B2ZB02) is the reference: every position is shown in its PmTYR colour
as the colour legend. Variants do NOT share B2ZB02 numbering, so each position is
transferred by cmd.cealign onto pmtyr + nearest-Ca mapping (tol 3.0 A).

Usage (PyMOL):
  run pymol_variants_active_site.py
  show_only 1_G2QC95        # flip to any structure
  show_only 0_B2ZB02        # the PmTYR colour legend
Render:
  ray 2400, 1800; png variant.png, dpi=300
"""

import os
from pymol import cmd, util

# Resolve the visualisation dir robustly: PyMOL's `run` sets __file__ to its own
# package path, so try __file__, then cwd, then the known absolute path, and pick
# whichever actually contains Structures/.
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
    227: 'high', 65: 'high', 68: 'high', 230: 'high', 195: 'high',  # 77-99%
    221: 'intermediate', 205: 'intermediate',
    46: 'hypervariable', 218: 'hypervariable', 209: 'hypervariable',
}
TIER_GREY = {'high': 'grey70', 'intermediate': 'grey80', 'hypervariable': 'grey90'}

# Six coordinating His scaffold + CuA His2 (thioether Cys partner).
ANCHOR_HIS = [42, 60, 69, 204, 208, 231]
CUA_HIS2 = 60

# Residue(s) the manuscript text discusses for each representative -> coloured.
# Numbering follows the order the substitutions appear in the manuscript:
# aromatics (Phe65, Trp68, Phe227) -> activity tuning (Glu195, His230)
# -> variable positions (Gly46, Asn205, Ala221).
DEFINING = {
    '0_B2ZB02':     [],                 # reference: coloured in full as the legend
    '1_G2QC95':     [65, 68], '2_Q2H7I7': [65, 68], '3_A0A8J2IPY8': [65],
    '4_A0ACC1SAF6': [68],  '5_A0AA39P123': [68],
    '6_A0A0D2BB96': [227],
    '7_A0A336U966': [195], '8_Q0CRX0':    [195], '9_A0A2B4SFS2': [195],
    '10_A0A2P4XQG5':[195], '11_R1EQM2':   [195],
    '12_W9YYG2':    [230], '13_U1TTP3':    [230], '14_A0A2R6XDH6':[230], '15_A0A6A3GT16':[230],
    '16_B1VTI5':    [46],  '17_D6RTB9':    [46],  '18_A0A397AZ25':[46],  '19_A0A8J5IJ88':[46],
    '20_A0A5N7APG9':[205], '21_T0JYN4':    [205],
    '22_A0A084B4D9':[221],
    '23_A0A9P9WW19':[195],   # Glu195->Tyr (fungal; Tyr Ca 0.29 A from PmTYR Glu195, pLDDT 99)
    '24_A0A8K0SIZ5':[209],   # Arg209->Glu (Stachybotrys, fungal; charge reversal at substrate-anchoring position)
    '25_A0A7H2XHW1':[209],   # Arg209->Cys (Acinetobacter, bacterial; thiol at substrate-anchoring position)
}
# Thioether state (from position_vectors.csv): 'C' genuine, 'C*' partial, '-' none.
THIO = {'12_W9YYG2': 'C', '13_U1TTP3': 'C', '14_A0A2R6XDH6': 'C*'}

REF_KEY = '0_B2ZB02'
VARIANT_KEYS = sorted([k for k in DEFINING if k != REF_KEY], key=lambda k: int(k.split('_')[0]))
MAPPED = {}  # obj -> {position: variant resi}


def setup_scene():
    cmd.bg_color('white')
    cmd.set('valence', 0)  # draw all bonds as single lines (no bond-order tick marks)
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
    """Show side-chain sticks (backbone N/O/C hidden)."""
    cmd.show('sticks', f'{obj} and resi {resi} and not name C+N+O')
    cmd.set('stick_radius', radius, f'{obj} and resi {resi}')


def color_structure(obj, key, is_ref=False):
    color_pos = set(HEATMAP_POSITIONS) if is_ref else set(DEFINING[key])

    cmd.show('cartoon', obj)
    cmd.color('grey90', obj)
    cmd.set('cartoon_transparency', 0.75, obj)

    # Cu: normal copper colour, 50% transparent spheres.
    cu = f'{obj} and name CU'
    if cmd.count_atoms(cu) > 0:
        cmd.show('spheres', cu)
        util.cnc(cu)  # restore element (copper) colour after the grey object colour
        cmd.set('sphere_scale', 0.25, cu)
        cmd.set('sphere_transparency', 0.5, cu)

    # Six coordinating His -> grey60 scaffold, 50% transparent sticks.
    his_resis = []
    for a in ANCHOR_HIS:
        rv = a if is_ref else map_position(obj, a)
        if rv is not None:
            _sticks(obj, rv, 0.18)
            cmd.color('grey60', f'{obj} and resi {rv} and not name C+N+O')
            cmd.set('stick_transparency', 0.5, f'{obj} and resi {rv}')
            his_resis.append(rv)

    # Coordination bonds: real His NE2 -> Cu bonds (canonical criterion <=3.5 A),
    # created with cmd.bond (the "Create Bond"/Shift-Cmd-T action). Radius and
    # transparency are applied per-bond with set_bond -- atom-level cmd.set does
    # NOT reliably reach these inter-residue bonds, so the Cu half was otherwise
    # rendering opaque and at the thick global default radius (0.25 > His 0.18).
    if his_resis and cmd.count_atoms(cu) > 0:
        for rv in his_resis:
            ne2 = f'{obj} and resi {rv} and name NE2'
            cu_lig = f'({cu}) within 3.5 of ({ne2})'
            if cmd.count_atoms(cu_lig) > 0:
                cmd.bond(ne2, cu_lig)
                cmd.set_bond('stick_radius', 0.15, ne2, cu_lig)       # < His 0.18
                cmd.set_bond('stick_transparency', 0.5, ne2, cu_lig)  # match His/Cu
        cmd.show('sticks', cu)

    # Ten active-site positions: text-mentioned -> colour; else grey by tier.
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
            cmd.label(f'{obj} and resi {rv} and name CA', f'"{POS_LABELS[pos]}"')
        else:
            _sticks(obj, rv, 0.20)
            cmd.color(TIER_GREY[CONSERVATION[pos]], f'{obj} and resi {rv} and not name C+N+O')

    # Genuine thioether Cys -> grey90 "other residue" (only where a real SG<=3.5 A
    # of the CuA His2 ring exists). C* partial / none are skipped.
    if THIO.get(key) == 'C':
        r60 = CUA_HIS2 if is_ref else map_position(obj, CUA_HIS2)
        if r60 is not None:
            ring = f'{obj} and resi {r60} and name ND1+NE2+CG+CD2+CE1'
            sg = f'{obj} and resn CYS and name SG within 3.5 of ({ring})'
            if cmd.count_atoms(sg) > 0:
                cys = f'byres ({sg})'
                cmd.show('sticks', f'({cys}) and not name C+N+O')
                cmd.color('grey90', f'({cys}) and not name C+N+O')
                cmd.set('stick_radius', 0.2, cys)

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
    print('\nReady. Flip structures with:  show_only <name>   e.g. show_only 1_G2QC95')
    print('Render: ray 2400, 1800; png variant.png, dpi=300')


cmd.extend('show_only', show_only)
build()
