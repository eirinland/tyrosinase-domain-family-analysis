#!/usr/bin/env python3
"""Summarise hmmscan domain hits per accession and per pool/group.

Reads domtblout files from hmmscan_results/, merges with pool assignments
and Chainsaw results. Outputs:
  domain_hits_per_accession.tsv  — all significant domain hits per protein
  domain_architecture.tsv        — one row per accession, domain string
  domain_summary_by_pool.tsv     — domain frequencies per pool
"""

import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path(__file__).parent.parent
OUTDIR = Path(__file__).parent
HMMSCAN_DIR = OUTDIR / "hmmscan_results"

POOLS = BASE / "1_filtering" / "final_pools" / "three_pool_assignment_final.csv"
CHAINSAW = OUTDIR / "chainsaw_results_all.csv"

PPO_PFAMS = {"PF00264", "PF18132"}  # Tyrosinase + Tyrosinase_C (PPO-associated)


def parse_domtblout(path):
    """Parse HMMER domtblout format, yield per-domain hits."""
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 23:
                continue
            if float(parts[12]) > 1e-3:
                continue
            yield {
                "target": parts[0],       # Pfam name
                "accession": parts[1],     # Pfam accession (PFxxxxx.xx)
                "query": parts[3],         # protein accession
                "evalue": float(parts[6]), # full sequence E-value
                "score": float(parts[7]),
                "dom_evalue": float(parts[12]),
                "dom_score": float(parts[13]),
                "env_from": int(parts[19]),
                "env_to": int(parts[20]),
                "desc": " ".join(parts[22:]),
            }


def main():
    # Load pool assignments
    pools = {}
    with open(POOLS) as f:
        for row in csv.DictReader(f):
            pools[row["accession"]] = row["pool"]

    # Load Chainsaw
    chainsaw = {}
    if CHAINSAW.exists():
        with open(CHAINSAW) as f:
            for row in csv.DictReader(f):
                chainsaw[row["accession"]] = {
                    "nres": int(row["nres"]),
                    "ndom": int(row["ndom"]),
                }

    # Parse all domtblout files
    print("Parsing hmmscan results...", file=sys.stderr)
    hits = defaultdict(list)  # accession -> [hit, ...]
    for fn in sorted(os.listdir(HMMSCAN_DIR)):
        if not fn.startswith("domtbl_"):
            continue
        for h in parse_domtblout(HMMSCAN_DIR / fn):
            pfam_acc = h["accession"].split(".")[0]  # PF00264.xx -> PF00264
            hits[h["query"]].append({
                "pfam_id": pfam_acc,
                "pfam_name": h["target"],
                "evalue": h["dom_evalue"],
                "score": h["dom_score"],
                "env_from": h["env_from"],
                "env_to": h["env_to"],
                "desc": h["desc"],
            })

    print(f"Hits for {len(hits)} accessions", file=sys.stderr)

    # Deduplicate overlapping hits (keep best score per Pfam per region)
    for acc in hits:
        by_pfam = defaultdict(list)
        for h in hits[acc]:
            by_pfam[h["pfam_id"]].append(h)
        deduped = []
        for pfam_id, domain_hits in by_pfam.items():
            domain_hits.sort(key=lambda x: x["evalue"])
            kept = []
            for h in domain_hits:
                overlap = False
                for k in kept:
                    ov = min(h["env_to"], k["env_to"]) - max(h["env_from"], k["env_from"])
                    span = max(h["env_to"] - h["env_from"], k["env_to"] - k["env_from"])
                    if ov > 0.5 * span:
                        overlap = True
                        break
                if not overlap:
                    kept.append(h)
            deduped.extend(kept)
        hits[acc] = sorted(deduped, key=lambda x: x["env_from"])

    # Write per-accession hits
    out1 = OUTDIR / "domain_hits_per_accession.tsv"
    with open(out1, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["accession", "pool", "length", "ndom_chainsaw",
                     "pfam_id", "pfam_name", "env_from", "env_to", "dom_evalue", "description"])
        for acc in sorted(hits):
            pool = pools.get(acc, "unknown")
            cs = chainsaw.get(acc, {})
            for h in hits[acc]:
                w.writerow([acc, pool, cs.get("nres", ""), cs.get("ndom", ""),
                           h["pfam_id"], h["pfam_name"],
                           h["env_from"], h["env_to"], f"{h['evalue']:.1e}", h["desc"]])
    print(f"Wrote {out1}", file=sys.stderr)

    # Write domain architecture (one row per accession)
    out2 = OUTDIR / "domain_architecture.tsv"
    with open(out2, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["accession", "pool", "length", "ndom_chainsaw", "n_pfam_domains",
                     "architecture", "extra_domains", "extra_pfam_ids"])
        for acc in sorted(pools):
            pool = pools[acc]
            cs = chainsaw.get(acc, {})
            acc_hits = hits.get(acc, [])
            arch = " | ".join(f"{h['pfam_name']}({h['env_from']}-{h['env_to']})" for h in acc_hits)
            extra = [h for h in acc_hits if h["pfam_id"] not in PPO_PFAMS]
            extra_str = " | ".join(f"{h['pfam_name']}({h['env_from']}-{h['env_to']})" for h in extra)
            extra_ids = ";".join(sorted(set(h["pfam_id"] for h in extra)))
            w.writerow([acc, pool, cs.get("nres", ""), cs.get("ndom", ""),
                       len(acc_hits), arch or "no_pfam_hit",
                       extra_str or "none", extra_ids or "none"])
    print(f"Wrote {out2}", file=sys.stderr)

    # Write summary by pool
    out3 = OUTDIR / "domain_summary_by_pool.tsv"
    pool_extra = defaultdict(Counter)  # pool -> Counter of pfam_ids
    pool_counts = Counter()
    pool_with_extra = Counter()
    for acc, pool in pools.items():
        pool_counts[pool] += 1
        acc_hits = hits.get(acc, [])
        extra = [h for h in acc_hits if h["pfam_id"] not in PPO_PFAMS]
        if extra:
            pool_with_extra[pool] += 1
        for h in extra:
            pool_extra[pool][f"{h['pfam_id']} ({h['pfam_name']})"] += 1

    with open(out3, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["pool", "n_total", "n_with_extra_domains", "pct_with_extra",
                     "top_extra_domains"])
        for pool in ["canonical", "noncanonical", "discarded"]:
            n = pool_counts[pool]
            ne = pool_with_extra[pool]
            top = "; ".join(f"{name}: {ct} ({100*ct/n:.1f}%)"
                           for name, ct in pool_extra[pool].most_common(10))
            w.writerow([pool, n, ne, f"{100*ne/n:.1f}" if n else "0", top])
            print(f"\n{pool} ({n}): {ne} with extra domains ({100*ne/n:.1f}%)", file=sys.stderr)
            for name, ct in pool_extra[pool].most_common(10):
                print(f"  {name}: {ct} ({100*ct/n:.1f}%)", file=sys.stderr)

    print(f"\nWrote {out3}", file=sys.stderr)


if __name__ == "__main__":
    main()
