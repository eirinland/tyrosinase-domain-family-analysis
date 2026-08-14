#!/usr/bin/env python3
"""Check coverage of 4 core PPO helices from Foldseek pass 2 alignments."""

import re
import sys

# Core helices in PDB residue numbers, with first_resi for 0-based conversion
# Foldseek uses 0-based sequential positions, so subtract first_resi
REF_FIRST_RESI = {
    "ref_PmTYR": 4,
    "ref_8BBR_Vspinosum": 34,
    "ref_2Y9W_Abisporus": 17,
    "ref_5CE9_Jregia": 1,
    "ref_1BT3_Ibatatas": 1,
    "ref_5M8L_human": 81,
    "ref_1JS8_squid": 2503,
    "ref_I3D139_archaea": 6,
    "ref_A0A9N8ELP9_oomycota": 347,
}

_CORE_HELICES_PDB = {
    "ref_PmTYR": {"a1": (34, 46), "a2": (65, 83), "a3": (203, 211), "a4": (226, 244)},
    "ref_8BBR_Vspinosum": {"a1": (73, 83), "a2": (91, 110), "a3": (259, 265), "a4": (279, 300)},
    "ref_2Y9W_Abisporus": {"a1": (54, 61), "a2": (90, 113), "a3": (255, 267), "a4": (291, 309)},
    "ref_5CE9_Jregia": {"a1": (80, 91), "a2": (113, 131), "a3": (241, 246), "a4": (269, 283)},
    "ref_1BT3_Ibatatas": {"a1": (81, 92), "a2": (114, 133), "a3": (240, 247), "a4": (269, 287)},
    "ref_5M8L_human": {"a1": (184, 196), "a2": (220, 239), "a3": (376, 383), "a4": (399, 417)},
    "ref_1JS8_squid": {"a1": (2536, 2543), "a2": (2567, 2584), "a3": (2660, 2679), "a4": (2697, 2719)},
    "ref_I3D139_archaea": {"a1": (36, 47), "a2": (67, 85), "a3": (197, 205), "a4": (222, 240)},
    "ref_A0A9N8ELP9_oomycota": {"a1": (413, 424), "a2": (437, 456), "a3": (590, 601), "a4": (624, 642)},
}

CORE_HELICES = {}
for ref, helices in _CORE_HELICES_PDB.items():
    offset = REF_FIRST_RESI[ref]
    CORE_HELICES[ref] = {
        h: (s - offset, e - offset) for h, (s, e) in helices.items()
    }


def parse_cigar(cigar):
    return [(op, int(n)) for n, op in re.findall(r"(\d+)([MID])", cigar)]


def get_aligned_target_positions(tstart, cigar):
    ops = parse_cigar(cigar)
    tpos = tstart
    aligned = set()
    for op, length in ops:
        if op == "M":
            for i in range(length):
                aligned.add(tpos)
                tpos += 1
        elif op == "D":
            tpos += length
    return aligned


def helix_coverage(aligned_positions, helix_start, helix_end):
    helix_len = helix_end - helix_start + 1
    covered = sum(1 for p in range(helix_start, helix_end + 1) if p in aligned_positions)
    return covered / helix_len


def check_all(tsv_path, output_csv=None, verbose=False):
    all_hits = {}
    with open(tsv_path) as f:
        header = f.readline().strip().split("\t")
        col = {name: i for i, name in enumerate(header)}
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < len(header):
                continue
            query = parts[col["query"]]
            target = parts[col["target"]]
            tstart = int(parts[col["tstart"]])
            qtmscore = float(parts[col["qtmscore"]])
            cigar = parts[col["cigar"]]
            acc = query.split("_taxID_")[0]
            all_hits.setdefault(acc, []).append({
                "target": target, "tstart": tstart,
                "qtmscore": qtmscore, "cigar": cigar,
            })

    results = []
    for acc in sorted(all_hits.keys()):
        hits = sorted(all_hits[acc], key=lambda x: -x["qtmscore"])
        best_result = None
        for info in hits:
            ref = info["target"]
            helices = CORE_HELICES.get(ref)
            if not helices:
                continue
            aligned = get_aligned_target_positions(info["tstart"], info["cigar"])
            covs = {h: helix_coverage(aligned, s, e) for h, (s, e) in helices.items()}
            n_ok = sum(1 for c in covs.values() if c >= 0.5)
            candidate = {
                "acc": acc, "ref": ref, "qtmscore": info["qtmscore"],
                "a1": covs["a1"], "a2": covs["a2"],
                "a3": covs["a3"], "a4": covs["a4"],
                "n_helices": n_ok,
            }
            if best_result is None or n_ok > best_result["n_helices"] or (
                n_ok == best_result["n_helices"] and info["qtmscore"] > best_result["qtmscore"]
            ):
                best_result = candidate
            if n_ok == 4:
                break
        if best_result:
            results.append(best_result)

    if verbose:
        hdr = "{:<16s} {:<24s} {:>6s}  {:>5s} {:>5s} {:>5s} {:>5s}  {:>6s}"
        print(hdr.format("accession", "ref", "qtm", "a1", "a2", "a3", "a4", "result"))
        print("-" * 90)
        for r in results:
            status = f"{r['n_helices']}/4"
            print(
                f"{r['acc']:<16s} {r['ref']:<24s} {r['qtmscore']:6.3f}"
                f"  {r['a1']:5.2f} {r['a2']:5.2f}"
                f" {r['a3']:5.2f} {r['a4']:5.2f}"
                f"  {status}"
            )

    if output_csv:
        import csv
        with open(output_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=[
                "accession", "ref", "qtmscore",
                "a1_cov", "a2_cov", "a3_cov", "a4_cov", "n_helices",
            ])
            w.writeheader()
            for r in results:
                w.writerow({
                    "accession": r["acc"], "ref": r["ref"],
                    "qtmscore": f"{r['qtmscore']:.4f}",
                    "a1_cov": f"{r['a1']:.3f}",
                    "a2_cov": f"{r['a2']:.3f}",
                    "a3_cov": f"{r['a3']:.3f}",
                    "a4_cov": f"{r['a4']:.3f}",
                    "n_helices": r["n_helices"],
                })

    n4 = sum(1 for r in results if r["n_helices"] == 4)
    n3 = sum(1 for r in results if r["n_helices"] == 3)
    n2 = sum(1 for r in results if r["n_helices"] <= 2)
    print(f"\nSummary: {len(results)} structures, {n4} with 4/4, {n3} with 3/4, {n2} with <=2/4")
    return results


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("tsv", help="Foldseek pass 2 TSV with backtrace cigar")
    parser.add_argument("--output", "-o", help="Output CSV path")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()
    check_all(args.tsv, output_csv=args.output, verbose=args.verbose)


if __name__ == "__main__":
    main()
