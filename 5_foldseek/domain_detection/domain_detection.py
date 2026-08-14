#!/usr/bin/env python3
"""Unified domain detection for PPO structures.

Two modes:
  ntd  — N-terminal capping domain (sliding-window Kabsch)
  ctd  — C-terminal blocking domain (anchor-and-extend helix Kabsch)

Both extract the query domain on-the-fly from a CIF, run against all
target CIFs in a directory, and report RMSD + pLDDT + nearby aromatics.
"""
import argparse
import csv
import json
import os
import sys
from collections import Counter

import numpy as np
from Bio.PDB import MMCIFParser

AROMATICS = {"PHE", "TYR", "TRP"}


# ── CIF helpers ──────────────────────────────────────────────────────────

def get_ca_info(structure):
    """Return [(resid, resname, coord_array), ...] for CA atoms, first chain."""
    cas = []
    for model in structure:
        for chain in model:
            for res in chain:
                if "CA" in res:
                    cas.append((res.id[1], res.get_resname(),
                                res["CA"].get_vector().get_array()))
        break
    return cas


def get_plddt(cif_path):
    """Return {seq_id: pLDDT} from B-factor column, chain A CA atoms."""
    result = {}
    cols, in_loop = {}, False
    col_idx = 0
    with open(cif_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith("loop_"):
                in_loop, cols, col_idx = True, {}, 0
                continue
            if in_loop and line.startswith("_atom_site."):
                cols[line.split(".")[1]] = col_idx
                col_idx += 1
                continue
            if not in_loop or not cols:
                continue
            if line.startswith("ATOM") or line.startswith("HETATM"):
                parts = line.split()
                if len(parts) < col_idx:
                    continue
                if parts[cols.get("label_atom_id", -1)] != "CA":
                    continue
                if parts[cols.get("label_asym_id", -1)] != "A":
                    continue
                seq_id = int(parts[cols["label_seq_id"]])
                plddt = float(parts[cols["B_iso_or_equiv"]])
                result[seq_id] = plddt
            elif line.startswith("#") or (line.startswith("_") and
                                          not line.startswith("_atom_site")):
                in_loop, cols = False, {}
    return result


# ── Kabsch ───────────────────────────────────────────────────────────────

def kabsch_rmsd(q, t):
    q = q - q.mean(axis=0)
    t = t - t.mean(axis=0)
    H = q.T @ t
    U, S, Vt = np.linalg.svd(H)
    return np.sqrt(max(0, (np.sum(q**2) + np.sum(t**2) - 2 * np.sum(S)) / len(q)))


def kabsch_transform(q, t):
    q_mean, t_mean = q.mean(axis=0), t.mean(axis=0)
    qc, tc = q - q_mean, t - t_mean
    H = qc.T @ tc
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, q_mean, t_mean


def nearest_aromatic(t_cas, ref_pos, max_dist=20.0):
    best_d, best_resid, best_name = float("inf"), None, None
    for resid, resname, coords in t_cas:
        if resname not in AROMATICS:
            continue
        d = np.linalg.norm(np.array(coords) - ref_pos)
        if d < best_d:
            best_d, best_resid, best_name = d, resid, resname
    if best_d > max_dist:
        return None, None, None
    return best_resid, best_name, best_d


# ── NTD: sliding-window Kabsch ───────────────────────────────────────────

def sliding_window_detect(q_coords, t_cas, ref_resid_offset):
    """Slide query domain over target CAs, return best RMSD + position."""
    n_q = len(q_coords)
    t_coords = np.array([c for _, _, c in t_cas])
    n_t = len(t_coords)
    if n_t < n_q:
        return None, None, None

    best_rmsd, best_j = float("inf"), -1
    for j in range(n_t - n_q + 1):
        rmsd = kabsch_rmsd(q_coords, t_coords[j:j + n_q])
        if rmsd < best_rmsd:
            best_rmsd, best_j = rmsd, j

    start_resid = t_cas[best_j][0]
    end_j = min(best_j + n_q - 1, n_t - 1)

    R, q_mean, t_mean = kabsch_transform(q_coords, t_coords[best_j:best_j + n_q])
    return best_rmsd, start_resid, (R, q_mean, t_mean)


# ── CTD: anchor-and-extend ───────────────────────────────────────────────

def anchor_extend_detect(q_helix_coords, anchor_indices, flank_indices,
                         helices, anchor_start, anchor_span,
                         anchor_helix_offsets, t_cas, search_radius):
    t_coords = np.array([c for _, _, c in t_cas])
    n_t = len(t_coords)
    if n_t < anchor_span:
        return None, None, None, None, 0, {}

    q_anchor = np.array([q_helix_coords[anchor_indices[0]][0].mean(axis=0)])
    idx_arr = np.array(anchor_helix_offsets)
    q_anchor_coords = np.concatenate([q_helix_coords[i] for i in anchor_indices])

    best_rmsd, best_j = float("inf"), -1
    for j in range(n_t - anchor_span + 1):
        t_helix = t_coords[j + idx_arr]
        rmsd = kabsch_rmsd(q_anchor_coords, t_helix)
        if rmsd < best_rmsd:
            best_rmsd, best_j = rmsd, j

    helix_positions = {}
    for i in anchor_indices:
        offset = helices[i]["start"] - anchor_start
        helix_positions[i] = best_j + offset

    for i in flank_indices:
        expected_offset = helices[i]["start"] - anchor_start
        expected_j = best_j + expected_offset
        hlen = len(q_helix_coords[i])
        lo = max(0, expected_j - search_radius)
        hi = min(n_t - hlen, expected_j + search_radius)
        if lo > hi:
            continue
        best_f_rmsd, best_f_j = float("inf"), -1
        for j in range(lo, hi + 1):
            rmsd = kabsch_rmsd(q_helix_coords[i], t_coords[j:j + hlen])
            if rmsd < best_f_rmsd:
                best_f_rmsd, best_f_j = rmsd, j
        if best_f_j >= 0:
            helix_positions[i] = best_f_j

    n_found = len(helix_positions)
    q_found = np.concatenate([q_helix_coords[i] for i in sorted(helix_positions)])
    t_found = np.concatenate([t_coords[helix_positions[i]:helix_positions[i] + len(q_helix_coords[i])]
                               for i in sorted(helix_positions)])
    rmsd = kabsch_rmsd(q_found, t_found) if len(q_found) >= 3 else None

    q_anchor_only = np.concatenate([q_helix_coords[i] for i in anchor_indices if i in helix_positions])
    t_anchor_only = np.concatenate([t_coords[helix_positions[i]:helix_positions[i] + len(q_helix_coords[i])]
                                     for i in anchor_indices if i in helix_positions])
    anchor_rmsd = kabsch_rmsd(q_anchor_only, t_anchor_only) if len(q_anchor_only) >= 3 else None

    R, q_mean, t_mean = kabsch_transform(q_found, t_found)
    start_resid = t_cas[min(helix_positions.values())][0]
    return rmsd, anchor_rmsd, start_resid, (R, q_mean, t_mean), n_found, helix_positions


# ── NTD mode ─────────────────────────────────────────────────────────────

def run_ntd(args):
    cif_parser = MMCIFParser(QUIET=True)
    ref_struct = cif_parser.get_structure("ref", args.ref_cif)
    ref_cas = get_ca_info(ref_struct)
    ref_map = {r: (n, np.array(c)) for r, n, c in ref_cas}

    q_coords = np.array([ref_map[r][1] for r in range(args.domain_start, args.domain_end + 1)
                          if r in ref_map])
    n_q = len(q_coords)
    print(f"NTD query: {args.ref_cif} residues {args.domain_start}-{args.domain_end} ({n_q} CAs)")

    ref_trp_offset = args.trp_offset
    print(f"Trp search: domain_end + 1 to domain_end + {args.trp_search_range}")

    targets = sorted(f for f in os.listdir(args.targets) if f.endswith(".cif"))
    print(f"Targets: {len(targets)}\n")

    rows = []
    for ti, fname in enumerate(targets):
        name = fname.replace(".cif", "")
        path = os.path.join(args.targets, fname)
        try:
            struct = cif_parser.get_structure("t", path)
            t_cas = get_ca_info(struct)
            n_t = len(t_cas)

            rmsd, start_resid, transform = sliding_window_detect(q_coords, t_cas, 0)
            if rmsd is None:
                rows.append({"name": name, "domain_rmsd": None, "domain_start": None,
                             "domain_plddt": None, "trp_resid": None, "trp_dist": None,
                             "n_ca": n_t})
                continue

            plddt_dict = get_plddt(path)
            t_resids = [r for r, _, _ in t_cas]
            best_j = t_resids.index(start_resid) if start_resid in t_resids else 0
            domain_resids = [t_cas[best_j + k][0] for k in range(n_q) if best_j + k < n_t]
            dom_plddt_vals = [plddt_dict[r] for r in domain_resids if r in plddt_dict]
            dom_plddt = np.mean(dom_plddt_vals) if dom_plddt_vals else None

            # Trp search after domain end
            domain_end_j = best_j + n_q
            trp_resid, trp_dist = None, None
            for k in range(domain_end_j, min(domain_end_j + args.trp_search_range, n_t)):
                resid, resname, coords = t_cas[k]
                if resname == "TRP":
                    R, q_mean, t_mean = transform
                    ref_trp_pos = R @ (ref_map.get(args.domain_end + ref_trp_offset, (None, np.zeros(3)))[1] - q_mean) + t_mean
                    trp_dist = np.linalg.norm(np.array(coords) - np.array(t_cas[k][2]))
                    trp_resid = resid
                    break

            rows.append({"name": name, "domain_rmsd": rmsd, "domain_start": start_resid,
                         "domain_plddt": dom_plddt, "trp_resid": trp_resid,
                         "trp_dist": trp_dist, "n_ca": n_t})
        except Exception as e:
            print(f"  ERROR {name}: {e}", file=sys.stderr, flush=True)
            rows.append({"name": name, "domain_rmsd": None, "domain_start": None,
                         "domain_plddt": None, "trp_resid": None, "trp_dist": None,
                         "n_ca": -1})

        if (ti + 1) % 100 == 0:
            print(f"  Processed {ti + 1}/{len(targets)}...", flush=True)

    with open(args.output, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sequence_id", "domain_rmsd", "domain_start", "domain_plddt",
                     "trp_resid", "trp_dist", "n_ca"])
        for r in sorted(rows, key=lambda x: x["name"]):
            def fmt(v):
                return f"{v:.3f}" if isinstance(v, float) else ("NA" if v is None else v)
            w.writerow([r["name"], fmt(r["domain_rmsd"]), r.get("domain_start", "NA"),
                        fmt(r["domain_plddt"]), r.get("trp_resid", "NA"),
                        fmt(r.get("trp_dist")), r["n_ca"]])

    valid = [r for r in rows if r["domain_rmsd"] is not None]
    rmsds = sorted(r["domain_rmsd"] for r in valid)
    print(f"\n{len(valid)}/{len(rows)} targets with results")
    if rmsds:
        for c in [1.0, 2.0, 3.0, 4.0, 5.0]:
            n = sum(1 for r in rmsds if r <= c)
            print(f"  RMSD <= {c:.1f}: {n}/{len(valid)} ({100 * n / len(valid):.1f}%)")
        print(f"  Median: {rmsds[len(rmsds) // 2]:.3f}")

    has_dom = [r for r in valid if r["domain_rmsd"] <= 2.0
               and r["domain_start"] is not None
               and r["domain_start"] <= 80]
    n_with_trp = sum(1 for r in has_dom if r["trp_resid"] is not None)
    print(f"\nWith domain (RMSD <= 2.0, start <= 80): {len(has_dom)}")
    print(f"  With Trp post-domain: {n_with_trp}/{len(has_dom)}")

    if has_dom and has_dom[0].get("domain_plddt") is not None:
        plddts = sorted(r["domain_plddt"] for r in has_dom if r["domain_plddt"])
        if plddts:
            print(f"  Domain pLDDT: median {plddts[len(plddts)//2]:.1f}, "
                  f"min {plddts[0]:.1f}, max {plddts[-1]:.1f}")


# ── CTD mode ─────────────────────────────────────────────────────────────

def run_ctd(args):
    helices = json.loads(args.helices)
    anchor_indices = [int(x) for x in args.anchor_idx.split(",")]
    flank_indices = [i for i in range(len(helices)) if i not in anchor_indices]

    for i, h in enumerate(helices):
        h["length"] = h["end"] - h["start"] + 1
        tag = "ANCHOR" if i in anchor_indices else "flank"
        print(f"  Helix {i}: res {h['start']}-{h['end']} ({h['length']} CAs) [{tag}]")

    anchor_start = min(helices[i]["start"] for i in anchor_indices)
    anchor_end = max(helices[i]["end"] for i in anchor_indices)
    anchor_span = anchor_end - anchor_start + 1

    anchor_helix_offsets = []
    for i in anchor_indices:
        for r in range(helices[i]["start"], helices[i]["end"] + 1):
            anchor_helix_offsets.append(r - anchor_start)

    cif_parser = MMCIFParser(QUIET=True)
    query_struct = cif_parser.get_structure("q", args.ref_cif)
    all_q_cas = get_ca_info(query_struct)
    q_map = {r: (n, np.array(c)) for r, n, c in all_q_cas}

    q_helix_coords = []
    for h in helices:
        coords = np.array([q_map[r][1] for r in range(h["start"], h["end"] + 1)])
        q_helix_coords.append(coords)

    ref_ca = q_map[args.ref_resid][1]
    ref_name = q_map[args.ref_resid][0]
    ref_in_anchor_offset = args.ref_resid - anchor_start
    print(f"\nReference residue: {ref_name}{args.ref_resid}")
    print(f"Anchor span: {anchor_start}-{anchor_end} ({anchor_span} res, "
          f"{len(anchor_helix_offsets)} helix CAs)")

    targets = sorted(f for f in os.listdir(args.targets) if f.endswith(".cif"))
    print(f"Targets: {len(targets)}\n")

    rows = []
    for ti, fname in enumerate(targets):
        name = fname.replace(".cif", "")
        path = os.path.join(args.targets, fname)
        try:
            struct = cif_parser.get_structure("t", path)
            t_cas = get_ca_info(struct)
            n_t = len(t_cas)
            t_coords = np.array([c for _, _, c in t_cas])

            rmsd, anchor_rmsd, start_resid, transform, n_found, hpos = anchor_extend_detect(
                q_helix_coords, anchor_indices, flank_indices,
                helices, anchor_start, anchor_span,
                anchor_helix_offsets, t_cas, args.search_radius)

            if rmsd is None:
                rows.append({"name": name, "helix_rmsd": None, "anchor_rmsd": None, "domain_start": None,
                             "aligned_resid": None, "aligned_resname": None,
                             "nearest_arom_resid": None, "nearest_arom_name": None,
                             "nearest_arom_dist": None, "domain_plddt": None,
                             "n_ca": n_t, "n_helices_found": 0})
                continue

            R, q_mean, t_mean = transform
            ref_transformed = R @ (ref_ca - q_mean) + t_mean

            # Aligned residue at reference position
            best_j_anchor = hpos[anchor_indices[0]]
            ref_j = best_j_anchor + ref_in_anchor_offset
            aligned_resid = t_cas[ref_j][0] if 0 <= ref_j < n_t else None
            aligned_resname = t_cas[ref_j][1] if 0 <= ref_j < n_t else None

            arom_resid, arom_name, arom_dist = nearest_aromatic(t_cas, ref_transformed)

            # pLDDT over matched helix positions
            plddt_dict = get_plddt(path)
            dom_plddt_vals = []
            for i in sorted(hpos):
                j = hpos[i]
                for k in range(helices[i]["length"]):
                    if j + k < n_t:
                        r = t_cas[j + k][0]
                        if r in plddt_dict:
                            dom_plddt_vals.append(plddt_dict[r])
            dom_plddt = np.mean(dom_plddt_vals) if dom_plddt_vals else None

            rows.append({"name": name, "helix_rmsd": rmsd, "anchor_rmsd": anchor_rmsd, "domain_start": start_resid,
                         "aligned_resid": aligned_resid, "aligned_resname": aligned_resname,
                         "nearest_arom_resid": arom_resid, "nearest_arom_name": arom_name,
                         "nearest_arom_dist": arom_dist, "domain_plddt": dom_plddt,
                         "n_ca": n_t, "n_helices_found": n_found})
        except Exception as e:
            print(f"  ERROR {name}: {e}", file=sys.stderr, flush=True)
            rows.append({"name": name, "helix_rmsd": None, "anchor_rmsd": None, "domain_start": None,
                         "aligned_resid": None, "aligned_resname": None,
                         "nearest_arom_resid": None, "nearest_arom_name": None,
                         "nearest_arom_dist": None, "domain_plddt": None,
                         "n_ca": -1, "n_helices_found": 0})

        if (ti + 1) % 100 == 0:
            print(f"  Processed {ti + 1}/{len(targets)}...", flush=True)

    with open(args.output, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sequence_id", "helix_rmsd", "anchor_rmsd", "domain_start", "aligned_resid",
                     "aligned_resname", "nearest_arom_resid", "nearest_arom_name",
                     "nearest_arom_dist", "domain_plddt", "n_ca", "n_helices_found"])
        for r in sorted(rows, key=lambda x: x["name"]):
            def fmt(v):
                return f"{v:.3f}" if isinstance(v, float) else ("NA" if v is None else v)
            w.writerow([r["name"], fmt(r["helix_rmsd"]), fmt(r.get("anchor_rmsd")), r.get("domain_start", "NA"),
                        r.get("aligned_resid", "NA"), r.get("aligned_resname", "NA"),
                        r.get("nearest_arom_resid", "NA"), r.get("nearest_arom_name", "NA"),
                        fmt(r.get("nearest_arom_dist")), fmt(r.get("domain_plddt")),
                        r["n_ca"], r["n_helices_found"]])

    valid = [r for r in rows if r["helix_rmsd"] is not None]
    rmsds = sorted(r["helix_rmsd"] for r in valid)
    print(f"\n{len(valid)}/{len(rows)} targets with results")
    if rmsds:
        for c in [1.0, 2.0, 3.0, 4.0, 5.0]:
            n = sum(1 for r in rmsds if r <= c)
            print(f"  RMSD <= {c:.1f}: {n}/{len(valid)} ({100 * n / len(valid):.1f}%)")
        print(f"  Median: {rmsds[len(rmsds) // 2]:.3f}")

    a_rmsds = sorted(r["anchor_rmsd"] for r in valid if r["anchor_rmsd"] is not None)
    if a_rmsds:
        print(f"\nAnchor-only RMSD (2-helix, {len(anchor_helix_offsets)} CAs):")
        for c in [1.0, 2.0, 3.0, 4.0, 5.0]:
            n = sum(1 for r in a_rmsds if r <= c)
            print(f"  RMSD <= {c:.1f}: {n}/{len(a_rmsds)} ({100 * n / len(a_rmsds):.1f}%)")
        print(f"  Median: {a_rmsds[len(a_rmsds) // 2]:.3f}")

    has_ctd = [r for r in valid if r["helix_rmsd"] <= 4.0
               and r["domain_start"] is not None
               and r["domain_start"] >= r["n_ca"] // 2 - 50]
    print(f"\nWith CTD (RMSD <= 4.0, C-terminal half): {len(has_ctd)}")

    if has_ctd:
        n_full = sum(1 for r in has_ctd if r["n_helices_found"] == len(helices))
        print(f"  All {len(helices)} helices found: {n_full}/{len(has_ctd)}")

        aligned_counts = Counter(r["aligned_resname"] for r in has_ctd if r["aligned_resname"])
        print(f"\n  Gatekeeper residue at ref position {args.ref_resid}:")
        for aa, count in aligned_counts.most_common():
            print(f"    {aa:3s}: {count:4d} ({100 * count / len(has_ctd):.1f}%)")

        for cutoff in [3, 5, 8]:
            n = sum(1 for r in has_ctd
                    if r["nearest_arom_dist"] is not None and r["nearest_arom_dist"] <= cutoff)
            print(f"  Aromatic within {cutoff} A: {n}/{len(has_ctd)}")

        if has_ctd[0].get("domain_plddt") is not None:
            plddts = sorted(r["domain_plddt"] for r in has_ctd if r["domain_plddt"])
            if plddts:
                print(f"  Domain pLDDT: median {plddts[len(plddts)//2]:.1f}, "
                      f"min {plddts[0]:.1f}, max {plddts[-1]:.1f}")

    species = set()
    for r in (has_ctd if has_ctd else []):
        parts = r["name"].split("_taxID_")
        if len(parts) > 1:
            species.add(parts[1].split("_")[0])
    if species:
        print(f"  Distinct taxIDs: {len(species)}")


# ── CLI ──────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="mode", required=True)

    # NTD subcommand
    p_ntd = sub.add_parser("ntd", help="N-terminal capping domain detection")
    p_ntd.add_argument("--ref-cif", required=True,
                       help="Reference CIF with the capping domain")
    p_ntd.add_argument("--domain-start", type=int, required=True,
                       help="First residue of capping domain in reference")
    p_ntd.add_argument("--domain-end", type=int, required=True,
                       help="Last residue of capping domain in reference")
    p_ntd.add_argument("--trp-offset", type=int, default=6,
                       help="Expected Trp position relative to domain end (default 6)")
    p_ntd.add_argument("--trp-search-range", type=int, default=30,
                       help="Search window past domain end for Trp (default 30)")
    p_ntd.add_argument("--targets", required=True)
    p_ntd.add_argument("--output", required=True)

    # CTD subcommand
    p_ctd = sub.add_parser("ctd", help="C-terminal blocking domain detection")
    p_ctd.add_argument("--ref-cif", required=True,
                       help="Reference CIF with the C-terminal domain")
    p_ctd.add_argument("--helices", required=True,
                       help='JSON: [{"start":N,"end":N}, ...] — helix ranges in reference')
    p_ctd.add_argument("--anchor-idx", required=True,
                       help="Comma-separated 0-based indices of anchor helices")
    p_ctd.add_argument("--ref-resid", type=int, required=True,
                       help="Reference residue to track (gatekeeper position)")
    p_ctd.add_argument("--search-radius", type=int, default=15)
    p_ctd.add_argument("--targets", required=True)
    p_ctd.add_argument("--output", required=True)

    args = parser.parse_args()
    if args.mode == "ntd":
        run_ntd(args)
    elif args.mode == "ctd":
        run_ctd(args)


if __name__ == "__main__":
    main()
