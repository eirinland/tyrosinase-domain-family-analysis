#!/usr/bin/env python3
"""Final 3-pool PPO structural classification (M7 copper-anchored core-helix rule).

Pools (32,069 structures):
  canonical        iff canonical_criteria == True
                   (2 Cu, Cu-Cu 2.8-5.5 A, 6 His NE2 <=3.5 A, per-His Ca pLDDT >=70)
  else >=6 His     -> discarded   (has the full six-His canonical site but fails the
                                   geometry/pLDDT criteria; no analytical home -- not
                                   canonical, and the non-canonical analysis only studies
                                   His-site divergence, so six-His structures do not belong)
  else (<6 His)    -> non-canonical iff the copper-anchored core-helix test passes
                       (core_helix_check.py: foldseek-align the query vs the 5 Cu-bearing
                        refs, pick best-qtm ref, helix-anchored Kabsch + ICP onto that
                        ref's 4 core helices, then require the query's OWN backbone to be
                        helical AND close AND confident at all 4 core helices at the
                        copper-anchored position -- core_ok = all 4 accepted)
                   else discarded.

This M7 rule REPLACES the 2026-06-14 (M1 + cealign + pLDDT-floor) core gate for the
<6-His failed-canonical scope only. Canonical and the six-His fails are untouched, so the
canonical pool (21,893) is unchanged; only the non-canonical/discarded split moves.

Inputs (Olivia):
  ../canonical_criteria_all_ca.csv                 canonical flag + n_coord_his
  ./core_helix_filter/core_helix_results.tsv       M7 copper-anchored core test (core_ok);
                                                   one row per <6-His failed-canonical query
                                                   (9,898 rows). Missing accession = out of
                                                   scope or no result = core fail = discard.
Outputs -> ./final_pools/
"""
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
CRIT = os.path.join(HERE, "..", "canonical_criteria_all_ca.csv")
CORE = os.path.join(HERE, "core_helix_filter", "core_helix_results.tsv")
OUT  = os.path.join(HERE, "final_pools")

def fnum(x, d=None):
    try: return float(x)
    except (TypeError, ValueError): return d

os.makedirs(OUT, exist_ok=True)

# M7 copper-anchored core verdict (missing accession = out of scope / no result = fail)
core = {}
with open(CORE) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        core[r["accession"]] = r

rows = list(csv.DictReader(open(CRIT)))
labels = ["canonical", "noncanonical", "discarded"]
counts = {l: 0 for l in labels}
reasons = {}
out_rows = []

for r in rows:
    acc = r["accession"]
    is_canon = r["canonical"] == "True"
    nhis = int(fnum(r["n_coord_his"], 0))
    cr = core.get(acc)
    ok = bool(cr) and cr.get("core_ok") == "True"

    if is_canon:
        pool, why = "canonical", "canonical_geometry"
    elif nhis >= 6:
        pool, why = "discarded", "sixHis_fail_geometry"
    elif ok:
        pool, why = "noncanonical", "corehelix_core_ok"
    else:
        pool, why = "discarded", "corehelix_core_fail"

    counts[pool] += 1
    reasons[why] = reasons.get(why, 0) + 1
    out_rows.append((acc, pool, r["canonical"], nhis, ok,
                     (cr or {}).get("best_ref", ""),
                     (cr or {}).get("best_qtm", ""),
                     (cr or {}).get("n_helix_ok", ""),
                     why))

with open(os.path.join(OUT, "three_pool_assignment_final.csv"), "w", newline="") as f:
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["accession", "pool", "canonical", "n_his", "core_ok",
                "best_ref", "best_qtm", "n_helix_ok", "reason"])
    w.writerows(out_rows)

# accession lists: plain LF (downstream shell loops break on CRLF)
for l in labels:
    with open(os.path.join(OUT, f"{l}_accessions.csv"), "w", newline="\n") as f:
        f.write("accession\n")
        for a, p, *_ in out_rows:
            if p == l:
                f.write(a + "\n")

print(f"total {len(out_rows)}")
for l in labels:
    print(f"  {l:13} {counts[l]:6}  ({100*counts[l]/len(out_rows):.1f}%)")
print("\nreason breakdown:")
for why, n in sorted(reasons.items(), key=lambda kv: -kv[1]):
    print(f"  {why:24} {n}")
