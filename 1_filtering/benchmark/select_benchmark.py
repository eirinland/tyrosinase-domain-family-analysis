#!/usr/bin/env python3
"""Assemble a stratified benchmark set spanning the canonical / non-canonical /
discarded decision boundary, for evaluating quality-filtering methods.

Strata (each accession tagged with an `expected` prior label; the user verifies
by eye):
  ctrl_microbispora  A0A8H9LF69  Microbispora non-canonical anchor (must stay non-canon)
  ctrl_char          Q92396 / P00440 / B2ZB02  characterized PPOs, current pool
  canon_disagree     the 64 canonical pool members the helix-check scored <4/4
  canon_rand         random canonical (positive control: core present)
  noncanon_lowqtm    non-canonical scored 4/4 but global qTM<0.3  (possible FALSE POSITIVE / junk)
  noncanon_rand      random non-canonical
  discard_fulllen    discarded but full length (>=260 res)  (possible FALSE NEGATIVE / real)
  discard_partial    discarded, short (<200 res)  (expected genuine junk: half-bundle)
  discard_nohit      discarded with no Foldseek hit to any ref (expected junk)

The 9 reference PDBs are added as positive controls by run_benchmark.sh, not here.
"""
import argparse, csv, random, sys

SEED = 13
CONTROLS = {
    "A0A8H9LF69": ("ctrl_microbispora", "noncanonical"),
    "Q92396":     ("ctrl_char",         "noncanonical"),
    "P00440":     ("ctrl_char",         "canonical"),
    "B2ZB02":     ("ctrl_char",         "canonical"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", required=True)        # three_pool_assignment.csv
    ap.add_argument("--core", required=True)        # ppo_core_results.csv (candidates)
    ap.add_argument("--disagree", required=True)    # canonical_helix_disagreements.csv
    ap.add_argument("--bestbt", required=True)      # foldseek_multiref_best_bt.tsv
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    rng = random.Random(SEED)

    pool = {r["accession"]: r["pool"] for r in csv.DictReader(open(a.pool))}
    core = {r["accession"]: r for r in csv.DictReader(open(a.core))}

    # accession -> qlen, best_qtm from the multi-ref best-hit table
    qlen, hit_qtm = {}, {}
    for r in csv.DictReader(open(a.bestbt), delimiter="\t"):
        acc = r["query"].split("_taxID_")[0]
        try:
            qlen[acc] = int(r["qlen"]); hit_qtm[acc] = float(r["qtmscore"])
        except (ValueError, KeyError):
            pass

    disagree = []
    for r in csv.DictReader(open(a.disagree)):
        acc = r["accession"]
        if r.get("best_ref", "") != "NO_HIT":
            disagree.append(acc)

    rows, seen = [], set()

    def add(acc, stratum, expected):
        if acc in seen or acc not in pool:
            return
        seen.add(acc)
        rows.append({"accession": acc, "stratum": stratum, "expected": expected,
                     "source_pool": pool.get(acc, "NA"),
                     "qlen": qlen.get(acc, ""), "best_qtm": f"{hit_qtm.get(acc,'')}"})

    # controls first
    for acc, (st, exp) in CONTROLS.items():
        add(acc, st, exp)
    # canonical helix disagreements (all)
    for acc in disagree:
        add(acc, "canon_disagree", "canonical")
    # noncanon low-qTM 4/4 (all) — possible false positives
    nc_low = [acc for acc, r in core.items()
              if pool.get(acc) == "non_canonical" and r["core_ok"] == "True"
              and float(r["best_qtm"]) < 0.3]
    for acc in sorted(nc_low):
        add(acc, "noncanon_lowqtm", "UNKNOWN")
    # discarded full length (>=260) with a hit — possible false negatives
    disc_full = [acc for acc in core
                 if pool.get(acc) == "discarded" and qlen.get(acc, 0) >= 260]
    rng.shuffle(disc_full)
    for acc in disc_full[:30]:
        add(acc, "discard_fulllen", "UNKNOWN")
    # discarded short (<200) with a hit — expected junk
    disc_part = [acc for acc in core
                 if pool.get(acc) == "discarded" and 0 < qlen.get(acc, 0) < 200]
    rng.shuffle(disc_part)
    for acc in disc_part[:20]:
        add(acc, "discard_partial", "discard")
    # discarded with NO foldseek hit — expected junk
    disc_nohit = [acc for acc, p in pool.items()
                  if p == "discarded" and acc not in core]
    rng.shuffle(disc_nohit)
    for acc in disc_nohit[:10]:
        add(acc, "discard_nohit", "discard")
    # random canonical / non-canonical positive controls
    canon = [acc for acc, p in pool.items() if p == "canonical" and acc not in seen]
    rng.shuffle(canon)
    for acc in canon[:20]:
        add(acc, "canon_rand", "canonical")
    ncan = [acc for acc, p in pool.items() if p == "non_canonical" and acc not in seen]
    rng.shuffle(ncan)
    for acc in ncan[:20]:
        add(acc, "noncanon_rand", "noncanonical")

    with open(a.output, "w", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=["accession", "stratum", "expected",
                                          "source_pool", "qlen", "best_qtm"])
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    print("benchmark set:", len(rows), file=sys.stderr)
    for st, n in sorted(Counter(r["stratum"] for r in rows).items()):
        print(f"  {st}: {n}", file=sys.stderr)


if __name__ == "__main__":
    main()
