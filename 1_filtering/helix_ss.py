"""Unified structural metric extraction for PPO filtering pipeline.

Aligns each structure onto PmTYR using PyMOL super, then extracts:
  1. Helix coverage (backbone presence + P-SEA helical character) for h1-h4
  2. Per-residue pLDDT at the 6 Cu-coordinating His positions (reference probe)
  3. Residue identity at the 6 His reference positions
  4. HIS-within-radius scan at each reference position
  5. Cu ref probe distances (AF3 Cu vs PmTYR Cu after superposition)
  6. Cu-Cu distance between the two AF3 Cu atoms
  7. Cu-imidazole coordination: HIS count near AF3 Cu and probe Cu positions

Usage: python3 helix_ss.py <cif_dir> <pmtyr.pdb> <structure_list.txt> <out.tsv> [threads]
"""
import sys, os, csv, warnings
warnings.filterwarnings("ignore")
import numpy as np
from multiprocessing import Pool
import pymol2
from Bio.PDB import PDBParser
import biotite.structure.io.pdbx as pdbx
import biotite.structure as struc

HELIX = {"helix1": range(37, 43), "helix2": range(67, 78),
         "helix3": range(204, 209), "helix4": range(228, 238)}
HIS6 = [42, 60, 69, 204, 208, 231]
PLDDT_SHELL = 5.0
HIS_SCAN_RADIUS = 3.5
CU_REF_RETRY = 10.0
CU_IMID_RADIUS = 3.0
IMID_ATOMS = {"CG", "ND1", "CD2", "CE1", "NE2"}

cifdir, pmtyr_path, listfile, out = sys.argv[1:5]
WORKERS = int(sys.argv[5]) if len(sys.argv) > 5 else 1

accs = [l.strip() for l in open(listfile) if l.strip()]
cif = {fn.split("_taxID_")[0]: os.path.join(cifdir, fn)
       for fn in os.listdir(cifdir) if fn.endswith(".cif")}

pp = PDBParser(QUIET=True)
pm = pp.get_structure("pm", pmtyr_path)[0]
REF_HELIX = {h: np.array([pm['A'][r]['CA'].coord for r in rs
                           if r in pm['A'] and 'CA' in pm['A'][r]])
             for h, rs in HELIX.items()}
REF_CU = np.array([a.coord.copy() for a in pm.get_atoms() if a.element == 'CU'])
REF_HIS6 = np.array([pm['A'][r]['CA'].coord for r in HIS6])

_w_pm = None
_w_ref_his6 = None
_w_ref_cu = None


def _init_worker(pmtyr):
    global _w_pm, _w_ref_his6, _w_ref_cu
    _w_pm = pymol2.PyMOL()
    _w_pm.start()
    _w_pm.cmd.load(pmtyr, "ref")
    _w_pm.cmd.remove("ref and not polymer and not elem CU")

    _w_ref_his6 = []
    for h in HIS6:
        m = _w_pm.cmd.get_model(f"ref and chain A and resi {h} and name CA")
        _w_ref_his6.append(np.array(m.atom[0].coord))
    _w_ref_his6 = np.array(_w_ref_his6)

    m = _w_pm.cmd.get_model("ref and elem CU")
    _w_ref_cu = np.array(sorted([a.coord for a in m.atom], key=lambda c: c[0]))


def _extract_atoms(model):
    """Extract CA, Cu, and HIS imidazole atoms from a PyMOL model."""
    cas, resnames, bfacs, rids = [], [], [], []
    cu_coords = []
    his_imid = []
    for atom in model.atom:
        if atom.name == "CA" and atom.hetatm == 0:
            cas.append(np.array(atom.coord))
            resnames.append(atom.resn)
            bfacs.append(atom.b)
            rids.append(int(atom.resi) if atom.resi.lstrip('-').isdigit() else 0)
        if atom.symbol == "CU":
            cu_coords.append(np.array(atom.coord))
        if atom.resn == "HIS" and atom.name in IMID_ATOMS and atom.hetatm == 0:
            resi = int(atom.resi) if atom.resi.lstrip('-').isdigit() else 0
            his_imid.append((resi, np.array(atom.coord)))
    if not cas:
        return None
    return {
        "cas": np.array(cas),
        "resnames": np.array(resnames),
        "bfacs": np.array(bfacs),
        "rids": np.array(rids),
        "cu": cu_coords,
        "his_imid": his_imid,
    }


def _extract_aligned(path, acc):
    """Load CIF, align with super, extract atom data from aligned structure."""
    _w_pm.cmd.load(path, "query")
    try:
        rms_info = _w_pm.cmd.super("query", "ref")
    except Exception:
        _w_pm.cmd.delete("query")
        return None, None

    model = _w_pm.cmd.get_model("query")
    _w_pm.cmd.delete("query")

    data = _extract_atoms(model)
    if data is None:
        return None, None
    data["rms"] = rms_info[0] if rms_info else 999.0
    return data, rms_info


def _count_his_near_cu(cu_coord, his_imid):
    """Count unique HIS residues with imidazole atoms within CU_IMID_RADIUS of cu."""
    resis = set()
    for resi, coord in his_imid:
        if np.linalg.norm(coord - cu_coord) <= CU_IMID_RADIUS:
            resis.add(resi)
    return len(resis)


def _max_cu_ref_dist(cu_coords):
    if len(cu_coords) < 2:
        return 999.0
    cu = sorted(cu_coords, key=lambda c: c[0])
    return max(float(min(np.linalg.norm(qcu - rcu) for rcu in _w_ref_cu)) for qcu in cu[:2])


def _align_best(acc, p):
    """Align with super. If Cu ends up far from reference, retry with just the
    Cu-bearing domain (residues within 30 A of Cu center)."""
    data, rms_info = _extract_aligned(p, acc)
    if data is None:
        return None

    best_dist = _max_cu_ref_dist(data["cu"])
    if best_dist <= CU_REF_RETRY:
        return data

    best_his = sum(1 for rc in _w_ref_his6
                   if np.linalg.norm(data["cas"] - rc, axis=1).min() <= 5.0)

    if len(data["cu"]) < 2:
        return data

    _w_pm.cmd.load(p, "query2")
    cu_model = _w_pm.cmd.get_model("query2 and elem CU")
    if len(cu_model.atom) < 2:
        _w_pm.cmd.delete("query2")
        return data

    _w_pm.cmd.select("cu_domain",
        "(query2 and polymer within 30 of (query2 and elem CU)) or (query2 and elem CU)")
    _w_pm.cmd.remove("query2 and not cu_domain")
    _w_pm.cmd.delete("cu_domain")

    try:
        _w_pm.cmd.super("query2", "ref")
    except Exception:
        _w_pm.cmd.delete("query2")
        return data

    model2 = _w_pm.cmd.get_model("query2")
    _w_pm.cmd.delete("query2")

    data2 = _extract_atoms(model2)
    if data2 is None:
        return data

    d = _max_cu_ref_dist(data2["cu"])
    h = sum(1 for rc in _w_ref_his6
            if np.linalg.norm(data2["cas"] - rc, axis=1).min() <= 5.0)

    if d < best_dist and h >= best_his:
        data2["rms"] = data["rms"]
        return data2
    return data


def _compute_metrics(data, cif_path):
    cas = data["cas"]
    resnames = data["resnames"]
    bfacs = data["bfacs"]
    rids = data["rids"]
    cu_coords = data["cu"]
    his_imid = data["his_imid"]

    # Active-site pLDDT (mean B-factor of CAs within 8 A of reference Cu)
    shell = np.zeros(len(cas), dtype=bool)
    for cu in REF_CU:
        shell |= np.linalg.norm(cas - cu, axis=1) <= 8.0
    plddt = float(np.mean(bfacs[shell])) if shell.any() else None

    # P-SEA secondary structure from biotite (uses original CIF, not aligned)
    f = pdbx.CIFFile.read(cif_path)
    arr = pdbx.get_structure(f, model=1)
    arr = arr[struc.filter_amino_acids(arr)]
    chains = np.unique(arr.chain_id)
    arr = arr[arr.chain_id == chains[0]]
    sse = struc.annotate_sse(arr)
    bt_rids = struc.get_residues(arr)[0]
    ss = dict(zip(bt_rids.tolist(), sse.tolist()))

    row = {"active_site_plddt": round(plddt, 2) if plddt is not None else ""}
    minc = 1.0
    minhc = 1.0
    for h, rc in REF_HELIX.items():
        pres, heli = [], []
        for point in rc:
            d = np.linalg.norm(cas - point, axis=1)
            j = d.argmin()
            near = d[j] <= 3.0
            pres.append(near)
            heli.append(near and ss.get(int(rids[j]), 'c') == 'a')
        cov = float(np.mean(pres))
        hcov = float(np.mean(heli))
        row[h + "_cov"] = round(cov, 2)
        row[h + "_hcov"] = round(hcov, 2)
        minc = min(minc, cov)
        minhc = min(minhc, hcov)
    row["min_cov"] = round(minc, 2)
    row["min_hcov"] = round(minhc, 2)

    # Cu ref probe distances
    if len(cu_coords) >= 2:
        qcu = sorted(cu_coords, key=lambda c: c[0])
        for i, q in enumerate(qcu[:2], 1):
            d = min(np.linalg.norm(q - rcu) for rcu in _w_ref_cu)
            row[f"cu{i}_ref_dist"] = round(float(d), 2)
        row["cu_cu_distance"] = round(float(np.linalg.norm(qcu[0] - qcu[1])), 3)
    else:
        row["cu_cu_distance"] = ""

    # Cu-imidazole coordination: HIS count near each Cu
    if len(cu_coords) >= 2:
        qcu = sorted(cu_coords, key=lambda c: c[0])
        row["af3_cu1_his"] = _count_his_near_cu(qcu[0], his_imid)
        row["af3_cu2_his"] = _count_his_near_cu(qcu[1], his_imid)
    else:
        row["af3_cu1_his"] = 0
        row["af3_cu2_his"] = 0

    if len(_w_ref_cu) >= 2:
        row["probe_cu1_his"] = _count_his_near_cu(_w_ref_cu[0], his_imid)
        row["probe_cu2_his"] = _count_his_near_cu(_w_ref_cu[1], his_imid)
    else:
        row["probe_cu1_his"] = 0
        row["probe_cu2_his"] = 0

    # Per-His pLDDT, residue identity, and HIS scan
    for k, ref_ca in enumerate(REF_HIS6):
        d = np.linalg.norm(cas - ref_ca, axis=1)
        j = d.argmin()
        if d[j] <= 5.0:
            row[f"his{HIS6[k]}_plddt"] = round(float(bfacs[j]), 2)
            row[f"his{HIS6[k]}_resname"] = resnames[j]
        else:
            row[f"his{HIS6[k]}_plddt"] = ""
            row[f"his{HIS6[k]}_resname"] = ""
        his_nearby = any(resnames[m] == "HIS" and d[m] <= HIS_SCAN_RADIUS
                         for m in range(len(resnames)))
        row[f"his{HIS6[k]}_has_his"] = "1" if his_nearby else "0"

    return row


def _process_one(acc):
    p = cif.get(acc)
    if not p:
        return {"accession": acc, "status": "nf"}
    try:
        data = _align_best(acc, p)
        if data is None:
            return {"accession": acc, "status": "err"}
        row = _compute_metrics(data, p)
        row["accession"] = acc
        row["status"] = "ok"
        return row
    except Exception:
        return {"accession": acc, "status": "err"}


fn = ["accession", "status", "active_site_plddt",
      "helix1_cov", "helix2_cov", "helix3_cov", "helix4_cov", "min_cov",
      "helix1_hcov", "helix2_hcov", "helix3_hcov", "helix4_hcov", "min_hcov",
      "cu1_ref_dist", "cu2_ref_dist", "cu_cu_distance",
      "af3_cu1_his", "af3_cu2_his", "probe_cu1_his", "probe_cu2_his",
      "his42_plddt", "his60_plddt", "his69_plddt",
      "his204_plddt", "his208_plddt", "his231_plddt",
      "his42_resname", "his60_resname", "his69_resname",
      "his204_resname", "his208_resname", "his231_resname",
      "his42_has_his", "his60_has_his", "his69_has_his",
      "his204_has_his", "his208_has_his", "his231_has_his"]

print(f"Processing {len(accs)} structures with {WORKERS} workers...", flush=True)

with open(out, "w", newline="") as fout:
    w = csv.DictWriter(fout, fieldnames=fn, delimiter="\t", extrasaction="ignore")
    w.writeheader()

    with Pool(WORKERS, initializer=_init_worker, initargs=(pmtyr_path,)) as pool:
        for i, row in enumerate(pool.imap_unordered(_process_one, accs), 1):
            w.writerow(row)
            fout.flush()
            if i % 5000 == 0:
                print(f"  {i}/{len(accs)} processed ...", flush=True)

print(f"done {out} {len(accs)}")
