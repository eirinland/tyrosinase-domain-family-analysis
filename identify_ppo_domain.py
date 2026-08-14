"""
Identify PPO core domain by cross-referencing Foldseek alignment
regions (qstart-qend) with Chainsaw domain boundaries.
"""
import csv, sys
from collections import defaultdict


def parse_chopping(chopping_str):
    if not chopping_str or chopping_str == "NULL":
        return []
    domains = []
    for i, dom_str in enumerate(chopping_str.split(",")):
        segs = []
        for seg in dom_str.split("_"):
            a, b = seg.split("-")
            segs.append((int(a), int(b)))
        domains.append((i, segs))
    return domains


def domain_residue_set(domain):
    _, segs = domain
    res = set()
    for a, b in segs:
        res.update(range(a, b + 1))
    return res


def overlap_fraction(domain, qstart, qend):
    aligned = set(range(qstart, qend + 1))
    dom_res = domain_residue_set(domain)
    if not dom_res:
        return 0.0
    return len(dom_res & aligned) / len(dom_res)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--foldseek", required=True)
    parser.add_argument("--chainsaw", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-overlap", type=float, default=0.5,
                        help="Min fraction of domain overlapping aligned region to count as PPO")
    args = parser.parse_args()

    chainsaw = {}
    with open(args.chainsaw) as f:
        for r in csv.DictReader(f):
            chainsaw[r["accession"]] = r

    foldseek = {}
    with open(args.foldseek) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            acc = r["query"].split("_taxID_")[0]
            foldseek[acc] = r

    print(f"Chainsaw: {len(chainsaw)}, Foldseek: {len(foldseek)}", flush=True)

    fieldnames = [
        "accession", "ndom", "chopping", "ppo_domains", "ppo_range",
        "non_ppo_domains", "non_ppo_range",
        "alntmscore", "qtmscore", "qstart", "qend", "qlen",
    ]

    with open(args.output, "w", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=fieldnames)
        writer.writeheader()

        for acc in sorted(chainsaw):
            fs = foldseek.get(acc)
            cw = chainsaw[acc]
            ndom = int(cw["ndom"])
            domains = parse_chopping(cw["chopping"])

            row = {
                "accession": acc,
                "ndom": ndom,
                "chopping": cw["chopping"],
            }

            if not fs:
                row["ppo_domains"] = "no_foldseek"
                writer.writerow(row)
                continue

            qstart = int(fs["qstart"])
            qend = int(fs["qend"])
            row["alntmscore"] = fs["alntmscore"]
            row["qtmscore"] = fs["qtmscore"]
            row["qstart"] = qstart
            row["qend"] = qend
            row["qlen"] = fs["qlen"]

            if ndom == 0:
                row["ppo_domains"] = "no_chainsaw"
                writer.writerow(row)
                continue

            ppo_doms = []
            non_ppo_doms = []
            for d in domains:
                frac = overlap_fraction(d, qstart, qend)
                if frac >= args.min_overlap:
                    ppo_doms.append(d)
                else:
                    non_ppo_doms.append(d)

            if ppo_doms:
                ppo_res = set()
                for d in ppo_doms:
                    ppo_res.update(domain_residue_set(d))
                row["ppo_domains"] = "+".join("d%d" % d[0] for d in ppo_doms)
                row["ppo_range"] = "%d-%d" % (min(ppo_res), max(ppo_res))
            else:
                row["ppo_domains"] = "none"

            if non_ppo_doms:
                non_ppo_res = set()
                for d in non_ppo_doms:
                    non_ppo_res.update(domain_residue_set(d))
                row["non_ppo_domains"] = "+".join("d%d" % d[0] for d in non_ppo_doms)
                row["non_ppo_range"] = "%d-%d" % (min(non_ppo_res), max(non_ppo_res))

            writer.writerow(row)

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
