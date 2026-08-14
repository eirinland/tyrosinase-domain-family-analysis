#!/usr/bin/env python3
"""Copper-anchored 4-helix-bundle core test (replaces M1 + cealign + pLDDT floor).

For each query: foldseek-align vs the 5 Cu-bearing PPO refs, pick the best-qtm
ref, helix-anchored Kabsch+ICP onto that ref's 4 core helices, then test -- WHERE
the copper is coordinated -- whether the query actually has each helix:

  for each of the 4 core helices, take its helix-borne His anchors (loop anchors
  dropped), find the nearest query Ca, and accept the helix iff >=1 anchor has
      nearest-Ca dist <= DMAX  AND  query residue locally helical  AND  Ca pLDDT >= PMIN
  core_ok = all 4 core helices accepted.

A half-bundle (a lobe's helices absent) fails because its anchors find no
close+helical+confident query backbone. Unlike the old M1 test, helicity is read
on the query's OWN backbone AT the Cu-anchored position, not at globally-aligned
positions that can thread through stray/bleedover residues.

Stdlib + numpy only. Helix ranges/anchors are auto-derived per ref from the Cu
site (3 nearest His NE2 per Cu) intersected with hard-coded core-helix ranges.
"""
import argparse, csv, math, os, re
from collections import defaultdict, Counter
import numpy as np

AA3 = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS',
       'MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
IMIDAZOLE = {'ND1','NE2','CD2','CE1'}
HELIX_D3 = (4.8, 6.4)   # Ca i->i+3 in an alpha helix
HELIX_D4 = (5.4, 7.4)   # Ca i->i+4

# Core-helix ranges in each ref's own residue numbering (from ppo_core_check.py).
CORE_HELICES = {
    "ref_PmTYR":          [(34, 46), (65, 83), (203, 211), (226, 244)],
    "ref_2Y9W_Abisporus": [(54, 61), (90, 113), (255, 267), (291, 309)],
    "ref_5CE9_Jregia":    [(80, 91), (113, 131), (241, 246), (269, 283)],
    "ref_1BT3_Ibatatas":  [(81, 92), (114, 133), (240, 247), (269, 287)],
    "ref_1JS8_squid":     [(2536, 2543), (2567, 2584), (2660, 2679), (2697, 2719)],
}


# ---------------- parsing ----------------
def parse_cif(path):
    atoms = []; cols = []; inA = False; coll = False
    with open(path) as f:
        for raw in f:
            l = raw.strip()
            if l == 'loop_':
                coll = False; inA = False; cols = []; continue
            if l.startswith('_atom_site.'):
                coll = True; cols.append(l); continue
            if coll and cols:
                coll = False; inA = True
            if inA:
                if l.startswith('_') or l == '#' or not l:
                    break
                p = l.split()
                if len(p) != len(cols):
                    continue
                r = dict(zip(cols, p))
                try:
                    atoms.append({'elem': r['_atom_site.type_symbol'].upper(),
                        'atom': r['_atom_site.label_atom_id'].upper(),
                        'resn': r['_atom_site.label_comp_id'].upper(),
                        'seq': r['_atom_site.label_seq_id'],
                        'bfac': float(r.get('_atom_site.B_iso_or_equiv', '0') or 0),
                        'x': float(r['_atom_site.Cartn_x']),
                        'y': float(r['_atom_site.Cartn_y']),
                        'z': float(r['_atom_site.Cartn_z'])})
                except (KeyError, ValueError):
                    continue
    return atoms


def parse_pdb(path):
    atoms = []
    with open(path) as f:
        for line in f:
            if line[:6].strip() not in ('ATOM', 'HETATM'):
                continue
            atom = line[12:16].strip().upper(); resn = line[17:20].strip().upper()
            seq = line[22:26].strip(); elem = line[76:78].strip().upper() or (atom[0] if atom else '')
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                bfac = float(line[60:66] or 0)
            except ValueError:
                continue
            atoms.append({'elem': elem, 'atom': atom, 'resn': resn, 'seq': seq,
                          'bfac': bfac, 'x': x, 'y': y, 'z': z})
    return atoms


# ---------------- geometry ----------------
def d3(a, b): return math.sqrt((a['x']-b['x'])**2 + (a['y']-b['y'])**2 + (a['z']-b['z'])**2)
def vec(a): return np.array([a['x'], a['y'], a['z']])
def in_helix(t, ranges): return any(s <= t <= e for s, e in ranges)


def ca_map(atoms):
    return {int(a['seq']): vec(a) for a in atoms
            if a['atom'] == 'CA' and a['seq'].lstrip('-').isdigit()}


def ca_bfac_map(atoms):
    return {int(a['seq']): a['bfac'] for a in atoms
            if a['atom'] == 'CA' and a['seq'].lstrip('-').isdigit()}


def ordered_seqs(atoms):
    return sorted({int(a['seq']) for a in atoms
                   if a['atom'] == 'CA' and a['seq'].lstrip('-').isdigit()})


def kabsch(P, Q):
    cP, cQ = P.mean(0), Q.mean(0)
    H = (P - cP).T @ (Q - cQ); U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T); R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cQ - R @ cP


def parse_cigar(c): return [(op, int(n)) for n, op in re.findall(r'(\d+)([MID])', c)]


def cigar_pairs(qstart, tstart, cigar, qseqs, tseqs):
    qi = qstart - 1; ti = tstart - 1; pairs = []
    for op, L in parse_cigar(cigar):
        if op == 'M':
            for _ in range(L):
                if 0 <= qi < len(qseqs) and 0 <= ti < len(tseqs):
                    pairs.append((qseqs[qi], tseqs[ti]))
                qi += 1; ti += 1
        elif op == 'I': qi += L
        elif op == 'D': ti += L
    return pairs


def helix_seed(qca, qseqs, qstart, tstart, cigar, ref_ca, ref_seqs, ranges):
    P = []; Q = []
    for q, t in cigar_pairs(qstart, tstart, cigar, qseqs, ref_seqs):
        if in_helix(t, ranges) and q in qca and t in ref_ca:
            P.append(qca[q]); Q.append(ref_ca[t])
    if len(P) < 8:
        return None
    R, t = kabsch(np.array(P), np.array(Q))
    return R, t


def icp_refine(qca, R, t, ref_ca, rcs, cutoff=4.0, iters=6):
    qseqs = list(qca.keys()); Qbase = np.array([qca[s] for s in qseqs])
    ref_core = np.array([ref_ca[s] for s in rcs]); last = 999.0
    for _ in range(iters):
        Qt = (R @ Qbase.T).T + t; P = []; Q = []
        for rc in ref_core:
            dd = np.linalg.norm(Qt - rc, axis=1); j = int(dd.argmin())
            if dd[j] <= cutoff:
                P.append(Qbase[j]); Q.append(rc)
        if len(P) < 8:
            break
        R, t = kabsch(np.array(P), np.array(Q))
        Pt = (R @ np.array(P).T).T + t
        last = float(np.sqrt(np.mean(np.sum((Pt - np.array(Q))**2, 1))))
    return R, t, last


# ---------------- reference prep ----------------
def get_cu(atoms):
    return [a for a in atoms if a['elem'] == 'CU' or a['atom'] == 'CU' or a['resn'] == 'CU']


def catalytic_pair(cus):
    if len(cus) < 2:
        return cus
    best = None; bestkey = None
    for i in range(len(cus)):
        for j in range(i + 1, len(cus)):
            dd = d3(cus[i], cus[j]); key = (0 if 2.5 <= dd <= 6.0 else 1, dd)
            if bestkey is None or key < bestkey:
                bestkey = key; best = (cus[i], cus[j])
    return list(best)


def detect_anchors(atoms, cus):
    ne2 = [a for a in atoms if a['resn'] == 'HIS' and a['atom'] == 'NE2'
           and a['seq'].lstrip('-').isdigit()]
    per = {0: [], 1: []}
    for a in ne2:
        d0 = d3(a, cus[0]); d1 = d3(a, cus[1]); k = 0 if d0 < d1 else 1
        per[k].append((min(d0, d1), int(a['seq'])))
    anchors = []
    for k in (0, 1):
        per[k].sort(); anchors += [s for _, s in per[k][:3]]
    return sorted(set(anchors))


def prep_ref(name, atoms, ranges):
    cus = catalytic_pair(get_cu(atoms))
    ref_ca = ca_map(atoms); ref_seqs = ordered_seqs(atoms)
    rcs = [s for s in ref_ca if in_helix(s, ranges)]
    anchors = detect_anchors(atoms, cus)
    # group anchors by which core helix they fall in (loop anchors excluded)
    helix_groups = defaultdict(list)
    for p in anchors:
        for hi, (s, e) in enumerate(ranges):
            if s <= p <= e:
                helix_groups[hi].append(p); break
    return dict(name=name, ref_ca=ref_ca, ref_seqs=ref_seqs, rcs=rcs, ranges=ranges,
                anchors=anchors, helix_groups=dict(helix_groups))


def load_refs(pmtyr, refdir):
    refs = {"ref_PmTYR": prep_ref("ref_PmTYR", parse_cif(pmtyr), CORE_HELICES["ref_PmTYR"])}
    for name, ranges in CORE_HELICES.items():
        if name == "ref_PmTYR":
            continue
        pdb = os.path.join(refdir, name + ".pdb")
        if os.path.exists(pdb):
            refs[name] = prep_ref(name, parse_pdb(pdb), ranges)
    return refs


def load_hits(fs_tsv, refs):
    hits = defaultdict(dict)
    with open(fs_tsv) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 9:
                continue
            acc = p[0].split('_taxID_')[0]; ref = p[1]
            try:
                qstart = int(p[2]); tstart = int(p[5]); qtm = float(p[7]); cigar = p[-1]
            except ValueError:
                continue
            for suf in ('.pdb', '.cif'):
                if ref.endswith(suf):
                    ref = ref[:-len(suf)]
            if ref not in refs and ref.endswith('_A'):
                ref = ref[:-2]
            if ref not in refs:
                continue
            d = hits[acc]
            if ref not in d or qtm > d[ref][3]:
                d[ref] = (qstart, tstart, cigar, qtm)
    return hits


# ---------------- core test ----------------
def _helical_turns(qca, r):
    """Count overlapping helical turns covering residue r (a turn at j covers
    j..j+4, so j in [r-4, r]); a turn is genuine when Ca i->i+3 and i->i+4 match
    alpha-helix spacing."""
    cnt = 0
    for j in range(r - 4, r + 1):
        if all(k in qca for k in (j, j + 3, j + 4)):
            dd3 = float(np.linalg.norm(qca[j] - qca[j + 3]))
            dd4 = float(np.linalg.norm(qca[j] - qca[j + 4]))
            if HELIX_D3[0] <= dd3 <= HELIX_D3[1] and HELIX_D4[0] <= dd4 <= HELIX_D4[1]:
                cnt += 1
    return cnt


def residue_helical(qca, r, need):
    """True iff residue r -- or an immediately adjacent residue r-1/r+1 -- sits in
    a genuine helix (>= `need` covering helical turns). The r+/-1 fallback recovers
    Cu-His anchors that map one residue into the loop at a helix cap; symmetric, so
    the off-by-one can land on either terminus. Distance/pLDDT gates are unchanged."""
    if r is None:
        return False
    return any(_helical_turns(qca, rr) >= need for rr in (r, r - 1, r + 1))


def nearest_q(anchor_seq, ref_ca, qca_t):
    rc = ref_ca[anchor_seq]; bd = 1e9; bs = None
    for s, c in qca_t.items():
        dd = float(np.linalg.norm(c - rc))
        if dd < bd:
            bd = dd; bs = s
    return bs, bd


def core_test(qca, qca_b, qseqs, fshit, ref, dmax, pmin, need):
    qstart, tstart, cigar, qtm = fshit
    seed = helix_seed(qca, qseqs, qstart, tstart, cigar, ref['ref_ca'], ref['ref_seqs'], ref['ranges'])
    if seed is None:
        return False, 0, [], None
    R, t = seed
    R, t, icp = icp_refine(qca, R, t, ref['ref_ca'], ref['rcs'])
    qca_t = {s: R @ c + t for s, c in qca.items()}
    detail = []; nok = 0
    for hi in range(4):
        anchs = ref['helix_groups'].get(hi, [])
        best = None; hok = False
        for p in anchs:
            bs, bd = nearest_q(p, ref['ref_ca'], qca_t)
            hel = residue_helical(qca, bs, need); pl = qca_b.get(bs, 0.0)
            sat = (bd <= dmax) and hel and (pl >= pmin)
            cand = dict(helix=hi + 1, anchor=p, qres=bs, dist=round(bd, 2),
                        helical=hel, plddt=round(pl, 1), sat=sat)
            if best is None or (sat, -bd) > (best['sat'], -best['dist']):
                best = cand
            if sat:
                hok = True
        if best is None:
            best = dict(helix=hi + 1, anchor=None, qres=None, dist=None,
                        helical=False, plddt=None, sat=False)
        detail.append(best)
        nok += int(hok)
    return nok == 4, nok, detail, icp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--cif-dir', required=True)
    ap.add_argument('--acc-list', required=True)
    ap.add_argument('--ref-dir', required=True)
    ap.add_argument('--pmtyr', required=True)
    ap.add_argument('--fs-tsv', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--dmax', type=float, default=4.0)
    ap.add_argument('--pmin', type=float, default=70.0)
    ap.add_argument('--need', type=int, default=2)
    ap.add_argument('--verbose', action='store_true')
    a = ap.parse_args()

    accs = [l.strip() for l in open(a.acc_list) if l.strip()]
    acc2cif = {fn.split('_taxID_')[0]: fn for fn in os.listdir(a.cif_dir) if fn.endswith('.cif')}
    refs = load_refs(a.pmtyr, a.ref_dir)
    hits = load_hits(a.fs_tsv, refs)

    rows = []; cnt = Counter()
    for acc in accs:
        fn = acc2cif.get(acc)
        if not fn:
            rows.append(dict(accession=acc, core_ok='', note='no_cif')); cnt['no_cif'] += 1; continue
        qa = parse_cif(os.path.join(a.cif_dir, fn))
        qca = ca_map(qa); qca_b = ca_bfac_map(qa); qseqs = ordered_seqs(qa)
        d = hits.get(acc, {})
        if not d:
            rows.append(dict(accession=acc, core_ok='', note='no_fs_hit')); cnt['no_fs_hit'] += 1; continue
        best_ref = max(d, key=lambda r: d[r][3])
        ok, nok, detail, icp = core_test(qca, qca_b, qseqs, d[best_ref], refs[best_ref],
                                         a.dmax, a.pmin, a.need)
        cnt['core_ok' if ok else 'core_fail'] += 1
        row = dict(accession=acc, best_ref=best_ref, best_qtm=round(d[best_ref][3], 3),
                   icp_rmsd=round(icp, 2) if icp is not None else '',
                   n_helix_ok=nok, core_ok=ok, note='')
        for hd in detail:
            h = hd['helix']
            row[f'a{h}_anchor'] = hd['anchor'] if hd['anchor'] is not None else ''
            row[f'a{h}_qres'] = hd['qres'] if hd['qres'] is not None else ''
            row[f'a{h}_dist'] = hd['dist'] if hd['dist'] is not None else ''
            row[f'a{h}_helical'] = int(hd['helical'])
            row[f'a{h}_plddt'] = hd['plddt'] if hd['plddt'] is not None else ''
            row[f'a{h}_sat'] = int(hd['sat'])
        rows.append(row)
        if a.verbose:
            ds = '  '.join(f"a{hd['helix']}:{'OK' if hd['sat'] else '--'}"
                           f"(d={hd['dist']},hel={int(hd['helical'])},pl={hd['plddt']})" for hd in detail)
            print(f"{acc:14} {best_ref:20} qtm={d[best_ref][3]:.3f} "
                  f"{'CORE_OK' if ok else f'FAIL({nok}/4)':9}  {ds}")

    fields = ['accession', 'best_ref', 'best_qtm', 'icp_rmsd', 'n_helix_ok', 'core_ok']
    for h in (1, 2, 3, 4):
        fields += [f'a{h}_anchor', f'a{h}_qres', f'a{h}_dist', f'a{h}_helical', f'a{h}_plddt', f'a{h}_sat']
    fields += ['note']
    with open(a.output, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader(); w.writerows(rows)

    print(f"\n=== core test (dmax={a.dmax} pmin={a.pmin} need={a.need}) on {len(accs)} ===")
    for k in ('core_ok', 'core_fail', 'no_fs_hit', 'no_cif'):
        if cnt[k]:
            print(f"  {k:10}: {cnt[k]}")
    print(f"  written {a.output}")


if __name__ == '__main__':
    main()
