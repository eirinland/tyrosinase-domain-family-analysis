#!/usr/bin/env python3
"""Summarize the multi-method benchmark: per-stratum method calls, control checks,
pairwise method agreement, and the disagreement queue for human inspection."""
import csv, sys
from collections import Counter, defaultdict

# A method "calls core present" (structure belongs in a pool, not discarded):
CALL_THRESHOLDS = {
    "m1": "M1 foldseek-helix core_ok",
    "m2": "M2 ref-norm TM >=0.50",
    "m3": "M3 super cov>=0.55 & rmsd<=5.0",
    "m4": "M4 cealign cov>=0.55 & rmsd<=6.0",
    "m5": "M5 intrinsic 4-helix bundle",
    "m6": "M6 bundle AND PPO-identity",
}


def fnum(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def calls(r):
    return {
        "m1": r["m1_core_ok"] == "True",
        "m2": r["m2_best_ttm"] != "" and fnum(r["m2_best_ttm"]) >= 0.50,
        "m3": fnum(r["m3_super_cov"]) >= 0.55 and fnum(r["m3_super_rmsd"], 99) <= 5.0,
        "m4": fnum(r["m4_ce_cov"]) >= 0.55 and fnum(r["m4_ce_rmsd"], 99) <= 6.0,
        "m5": r["m5_bundle"] == "True",
        "m6": r["m6_combined"] == "True",
    }


def main():
    res_path, out_txt, out_dis = sys.argv[1], sys.argv[2], sys.argv[3]
    rows = list(csv.DictReader(open(res_path), delimiter="\t"))
    methods = list(CALL_THRESHOLDS)
    L = []

    L.append("PPO QUALITY-FILTERING METHOD BENCHMARK")
    L.append("=" * 70)
    L.append(f"structures scored: {len(rows)}")
    L.append("")
    L.append("Method = core present (belongs in canonical or non-canonical pool):")
    for m in methods:
        L.append(f"  {m.upper()}  {CALL_THRESHOLDS[m]}")
    L.append("")

    # per-stratum call rates
    L.append("PER-STRATUM: # called core-present / n   (expected in brackets)")
    L.append("-" * 70)
    by_str = defaultdict(list)
    for r in rows:
        by_str[r["stratum"]].append(r)
    order = ["ctrl_ref", "ctrl_microbispora", "ctrl_char", "canon_disagree",
             "canon_rand", "noncanon_rand", "noncanon_lowqtm",
             "discard_fulllen", "discard_partial", "discard_nohit"]
    header = "%-20s %4s  " % ("stratum", "n") + "  ".join("%5s" % m.upper() for m in methods) + "   expected"
    L.append(header)
    for st in order + [s for s in by_str if s not in order]:
        rs = by_str.get(st)
        if not rs:
            continue
        n = len(rs)
        cc = {m: sum(calls(r)[m] for r in rs) for m in methods}
        exp = Counter(r["expected"] for r in rs).most_common(1)[0][0]
        L.append("%-20s %4d  " % (st, n) + "  ".join("%5d" % cc[m] for m in methods) + f"   [{exp}]")
    L.append("")

    # control checks
    L.append("CONTROL CHECKS")
    L.append("-" * 70)
    refs = [r for r in rows if r["stratum"] == "ctrl_ref"]
    L.append(f"  9 reference cores: all 6 methods call core-present? "
             + ("YES" if all(all(calls(r).values()) for r in refs) else
                "NO -> " + ", ".join(r["accession"] for r in refs if not all(calls(r).values()))))
    mb = [r for r in rows if r["stratum"] == "ctrl_microbispora"]
    if mb:
        r = mb[0]; c = calls(r)
        L.append(f"  Microbispora A0A8H9LF69 (must be NON-CANONICAL = core present, geom fail):")
        L.append(f"     core-present calls: " + " ".join(f"{m.upper()}={'Y' if c[m] else 'N'}" for m in methods)
                 + f"   | Cu-Cu={r['cu_dist']} His={r['n_his']} pLDDT={r['min_plddt']} pool={r['pool']}")
    for r in [r for r in rows if r["stratum"] == "ctrl_char"]:
        c = calls(r)
        L.append(f"  {r['accession']} (char, pool={r['pool']}): "
                 + " ".join(f"{m.upper()}={'Y' if c[m] else 'N'}" for m in methods))
    L.append("")

    # pairwise agreement
    L.append("PAIRWISE METHOD AGREEMENT (% of structures with same call)")
    L.append("-" * 70)
    L.append("       " + "  ".join("%5s" % m.upper() for m in methods))
    for m1 in methods:
        cells = []
        for m2 in methods:
            agree = sum(calls(r)[m1] == calls(r)[m2] for r in rows) / len(rows) * 100
            cells.append("%5.0f" % agree)
        L.append("%-6s " % m1.upper() + "  ".join(cells))
    L.append("")

    # total core-present per method
    L.append("TOTAL called core-present (out of %d):" % len(rows))
    for m in methods:
        L.append(f"  {m.upper()}: {sum(calls(r)[m] for r in rows)}")
    L.append("")

    # disagreement queue
    dis = [r for r in rows if len(set(calls(r).values())) > 1]
    L.append(f"DISAGREEMENT QUEUE: {len(dis)} structures where methods disagree "
             f"(-> inspect these)")
    L.append(f"   written to: {out_dis}")
    L.append("")
    L.append("=" * 70)
    L.append("Next: open the disagreement queue in PyMOL, label each as "
             "canonical / non-canonical / junk, then pick the method (or threshold) "
             "that best reproduces your labels.")

    open(out_txt, "w").write("\n".join(L) + "\n")
    print("\n".join(L))

    cols = ["accession", "stratum", "expected", "pool", "qlen", "cu_dist", "n_his",
            "min_plddt", "m1_core_ok", "m1_nhel", "m1_qtm", "m2_best_ttm",
            "m3_super_cov", "m3_super_rmsd", "m4_ce_cov", "m4_ce_rmsd",
            "m5_nlonghelix", "m5_bundle", "m6_combined"]
    with open(out_dis, "w", newline="\n") as f:
        wr = csv.writer(f, delimiter="\t")
        wr.writerow(cols + ["calls", "n_call"])
        for r in sorted(dis, key=lambda r: (r["stratum"], r["accession"])):
            c = calls(r)
            wr.writerow([r.get(k, "") for k in cols]
                        + ["".join(m[1] if c[m] else "." for m in methods),
                           sum(c.values())])


if __name__ == "__main__":
    main()
