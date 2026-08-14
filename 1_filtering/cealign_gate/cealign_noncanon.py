#!/usr/bin/env python3
"""PyMOL cealign of candidate structures vs the 9 PPO references.

Replicates exactly the M4 cealign block of the benchmark run_methods.py:
  - ref_nca[ref] = number of CA residues in the reference PDB (chain A/blank)
  - for each query: cmd.cealign(ref+" and name CA", "mob and name CA")
    cov = alignment_length / ref_nca[ref]; keep the best-cov reference hit.
Benchmark gate (compare_gates.py): cov >= 0.55 AND rmsd <= 5.0.
"""
import argparse, csv, os, sys


def ref_nca_count(pdb):
    n = 0
    for line in open(pdb):
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        if line[21] not in (" ", "A"):
            continue
        n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif-dir", required=True)
    ap.add_argument("--ref-dir", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()

    refs = [f[:-4] for f in sorted(os.listdir(a.ref_dir)) if f.endswith(".pdb")]
    from pymol import cmd
    cmd.feedback("disable", "all", "everything")
    ref_nca = {}
    for ref in refs:
        p = os.path.join(a.ref_dir, ref + ".pdb")
        ref_nca[ref] = ref_nca_count(p)
        cmd.load(p, ref)
    print(f"loaded {len(refs)} refs", file=sys.stderr)

    queries = []
    for fn in sorted(os.listdir(a.cif_dir)):
        if fn.endswith(".cif"):
            queries.append((fn.split("_taxID_")[0], os.path.join(a.cif_dir, fn)))

    w = csv.writer(open(a.output, "w", newline="\n"))
    w.writerow(["accession", "best_ce_cov", "best_ce_rmsd", "best_ce_naln", "best_ce_ref"])
    for i, (acc, path) in enumerate(queries, 1):
        cmd.load(path, "mob")
        best = None
        for ref in refs:
            try:
                rc = cmd.cealign(ref + " and name CA", "mob and name CA")
                cov = rc["alignment_length"] / ref_nca[ref]
                if best is None or cov > best[0]:
                    best = (cov, rc["RMSD"], rc["alignment_length"], ref)
            except Exception:
                pass
        cmd.delete("mob")
        if best:
            w.writerow([acc, f"{best[0]:.3f}", f"{best[1]:.2f}", best[2], best[3]])
        else:
            w.writerow([acc, "", "", "", ""])
        if i % 50 == 0:
            print(f"  ...{i}/{len(queries)}", file=sys.stderr)
    print(f"cealign done: {len(queries)} structures", file=sys.stderr)


if __name__ == "__main__":
    main()
