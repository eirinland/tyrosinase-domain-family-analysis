#!/usr/bin/env python3
"""Summarise seed variability across AF3 diffusion samples."""

import csv
import sys
from collections import defaultdict
from pathlib import Path


def mean(vals):
    return sum(vals) / len(vals) if vals else None

def sd(vals):
    if len(vals) < 2:
        return None
    m = mean(vals)
    return (sum((v - m) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5

def fmt(v):
    return f'{v:.4f}' if v is not None else ''


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <results_dir> <pools.csv>", file=sys.stderr)
        sys.exit(1)

    results_dir = Path(sys.argv[1])
    pools_path = Path(sys.argv[2])

    pools = {}
    with open(pools_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            pools[row['accession']] = row['pool']

    data = defaultdict(list)
    for tsv in sorted(results_dir.glob('seed_geometry_*.tsv')):
        with open(tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row.get('error') and row['error'] not in ('', 'None'):
                    continue
                data[row['accession']].append(row)

    output = results_dir / 'seed_variability_summary.tsv'
    fields = [
        'accession', 'pool', 'n_samples',
        'cu_cu_mean', 'cu_cu_sd', 'cu_cu_min', 'cu_cu_max',
        'cu1_plddt_mean', 'cu1_plddt_sd',
        'cu2_plddt_mean', 'cu2_plddt_sd',
        'n_coordinating_his_values', 'his_agreement',
        'min_coord_his_plddt_mean', 'min_coord_his_plddt_sd',
        'n_canonical', 'all_canonical', 'any_canonical',
    ]

    with open(output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        writer.writeheader()

        for acc in sorted(data):
            rows = data[acc]
            samples = [r for r in rows if r['sample'] != 'top']
            if not samples:
                samples = rows

            n = len(samples)
            cu_cu = [float(r['cu_cu_dist']) for r in samples
                     if r.get('cu_cu_dist') not in ('', 'None', None)]
            cu1_p = [float(r['cu1_plddt']) for r in samples
                     if r.get('cu1_plddt') not in ('', 'None', None)]
            cu2_p = [float(r['cu2_plddt']) for r in samples
                     if r.get('cu2_plddt') not in ('', 'None', None)]
            n_his = [int(r['n_coordinating_his']) for r in samples
                     if r.get('n_coordinating_his') not in ('', 'None', None)]
            min_hp = [float(r['min_coord_his_ca_plddt']) for r in samples
                      if r.get('min_coord_his_ca_plddt') not in ('', 'None', None)]
            canonical = [r['canonical'] == 'True' for r in samples]

            his_vals = sorted(set(n_his))

            writer.writerow({
                'accession': acc,
                'pool': pools.get(acc, 'unknown'),
                'n_samples': n,
                'cu_cu_mean': fmt(mean(cu_cu)),
                'cu_cu_sd': fmt(sd(cu_cu)),
                'cu_cu_min': fmt(min(cu_cu)) if cu_cu else '',
                'cu_cu_max': fmt(max(cu_cu)) if cu_cu else '',
                'cu1_plddt_mean': fmt(mean(cu1_p)),
                'cu1_plddt_sd': fmt(sd(cu1_p)),
                'cu2_plddt_mean': fmt(mean(cu2_p)),
                'cu2_plddt_sd': fmt(sd(cu2_p)),
                'n_coordinating_his_values': ','.join(str(v) for v in his_vals),
                'his_agreement': 'yes' if len(his_vals) <= 1 else 'no',
                'min_coord_his_plddt_mean': fmt(mean(min_hp)),
                'min_coord_his_plddt_sd': fmt(sd(min_hp)),
                'n_canonical': sum(canonical),
                'all_canonical': 'yes' if all(canonical) else 'no',
                'any_canonical': 'yes' if any(canonical) else 'no',
            })

    # Print summary by pool
    pool_stats = defaultdict(lambda: {'n': 0, 'his_agree': 0, 'canon_agree': 0,
                                       'cu_cu_sds': [], 'cu1_plddt_sds': []})
    for acc, rows in data.items():
        pool = pools.get(acc, 'unknown')
        samples = [r for r in rows if r['sample'] != 'top']
        if not samples:
            samples = rows

        ps = pool_stats[pool]
        ps['n'] += 1

        n_his = [int(r['n_coordinating_his']) for r in samples
                 if r.get('n_coordinating_his') not in ('', 'None', None)]
        if len(set(n_his)) <= 1:
            ps['his_agree'] += 1

        canonical = [r['canonical'] == 'True' for r in samples]
        if len(set(canonical)) <= 1:
            ps['canon_agree'] += 1

        cu_cu = [float(r['cu_cu_dist']) for r in samples
                 if r.get('cu_cu_dist') not in ('', 'None', None)]
        s = sd(cu_cu)
        if s is not None:
            ps['cu_cu_sds'].append(s)

        cu1_p = [float(r['cu1_plddt']) for r in samples
                 if r.get('cu1_plddt') not in ('', 'None', None)]
        s = sd(cu1_p)
        if s is not None:
            ps['cu1_plddt_sds'].append(s)

    print(f"\nSeed variability summary ({len(data)} accessions):")
    for pool in ['canonical', 'noncanonical', 'discarded']:
        ps = pool_stats.get(pool)
        if not ps or ps['n'] == 0:
            continue
        n = ps['n']
        print(f"\n  {pool} ({n}):")
        print(f"    His count agreement:      {ps['his_agree']}/{n} ({100*ps['his_agree']/n:.1f}%)")
        print(f"    Canonical agreement:       {ps['canon_agree']}/{n} ({100*ps['canon_agree']/n:.1f}%)")
        if ps['cu_cu_sds']:
            med_sd = sorted(ps['cu_cu_sds'])[len(ps['cu_cu_sds'])//2]
            print(f"    Cu-Cu dist SD (median):    {med_sd:.4f} A")
        if ps['cu1_plddt_sds']:
            med_sd = sorted(ps['cu1_plddt_sds'])[len(ps['cu1_plddt_sds'])//2]
            print(f"    Cu1 pLDDT SD (median):     {med_sd:.2f}")

    print(f"\nWrote {output}")


if __name__ == '__main__':
    main()
