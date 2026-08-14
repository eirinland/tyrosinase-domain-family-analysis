#!/usr/bin/env python3
"""Submit AFDB structures to AlphaFill web API and compare Cu placement to AF3."""

import csv
import json
import math
import os
import sys
import time
import urllib.request
import urllib.error

ALPHAFILL_URL = "https://alphafill.eu/v1/aff"


def submit_structure(cif_path):
    with open(cif_path) as f:
        cif_content = f.read()
    data = json.dumps({"structure": cif_content}).encode()
    req = urllib.request.Request(
        ALPHAFILL_URL, data=data,
        headers={"Content-Type": "application/json"}, method="POST",
    )
    r = urllib.request.urlopen(req, timeout=60)
    return json.loads(r.read())


def poll_status(job_id, max_wait=300):
    url = "%s/%s/status" % (ALPHAFILL_URL, job_id)
    for _ in range(max_wait // 5):
        try:
            r = urllib.request.urlopen(url, timeout=10)
            status = json.loads(r.read())
            if status.get("status") in ("finished", "ok"):
                return True
            if status.get("status") == "error":
                return False
        except Exception:
            pass
        time.sleep(5)
    return False


def get_result_cu_atoms(job_id):
    url = "%s/%s" % (ALPHAFILL_URL, job_id)
    r = urllib.request.urlopen(url, timeout=60)
    data = r.read().decode()
    cu_atoms = []
    for line in data.split("\n"):
        parts = line.split()
        if len(parts) > 12 and parts[0] == "HETATM" and parts[2] == "CU":
            try:
                cu_atoms.append({
                    "x": float(parts[10]),
                    "y": float(parts[11]),
                    "z": float(parts[12]),
                    "chain": parts[6] if len(parts) > 6 else "?",
                })
            except (ValueError, IndexError):
                continue
    return cu_atoms


def dist(a, b):
    return math.sqrt((a["x"] - b["x"])**2 + (a["y"] - b["y"])**2 + (a["z"] - b["z"])**2)


def best_cu_pair(cu_atoms):
    """Find the Cu pair closest to canonical range (2.8-5.5 A)."""
    if len(cu_atoms) < 2:
        return None, None
    best_dist = None
    best_pair = None
    for i in range(len(cu_atoms)):
        for j in range(i + 1, len(cu_atoms)):
            d = dist(cu_atoms[i], cu_atoms[j])
            if best_dist is None or abs(d - 4.0) < abs(best_dist - 4.0):
                best_dist = d
                best_pair = (i, j)
    return best_dist, best_pair


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--metal3d-csv", required=True)
    parser.add_argument("--afdb-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--delay", type=float, default=2.0,
                        help="Seconds between submissions")
    args = parser.parse_args()

    m3d = {}
    with open(args.metal3d_csv) as f:
        for r in csv.DictReader(f):
            if r["status"] == "ok" and r["predicted_cu_dist"]:
                m3d[r["accession"]] = {
                    "af3_cu_dist": float(r["af3_cu_dist"]),
                    "m3d_cu_dist": float(r["predicted_cu_dist"]),
                }

    print("Metal3D structures to process: %d" % len(m3d))

    afdb_files = {}
    for fname in os.listdir(args.afdb_dir):
        if fname.endswith(".cif"):
            acc = fname.split("-")[1] if fname.startswith("AF-") else None
            if acc:
                afdb_files[acc] = os.path.join(args.afdb_dir, fname)

    tasks = [(acc, afdb_files[acc]) for acc in m3d if acc in afdb_files]
    print("AFDB CIFs found: %d" % len(tasks))

    results = []
    for i, (acc, cif_path) in enumerate(sorted(tasks)):
        sys.stdout.write("  [%d/%d] %s ... " % (i + 1, len(tasks), acc))
        sys.stdout.flush()
        try:
            resp = submit_structure(cif_path)
            job_id = resp["id"]
            if resp.get("status") == "too-many-requests":
                print("rate limited, waiting 60s")
                time.sleep(60)
                resp = submit_structure(cif_path)
                job_id = resp["id"]

            ok = poll_status(job_id)
            if not ok:
                print("timeout/error")
                results.append({"accession": acc, "alphafill_status": "error"})
                continue

            cu_atoms = get_result_cu_atoms(job_id)
            cu_dist, pair = best_cu_pair(cu_atoms)

            info = m3d[acc]
            row = {
                "accession": acc,
                "alphafill_status": "ok",
                "af3_cu_dist": round(info["af3_cu_dist"], 2),
                "m3d_cu_dist": round(info["m3d_cu_dist"], 2),
                "alphafill_n_cu": len(cu_atoms),
                "alphafill_cu_dist": round(cu_dist, 2) if cu_dist else "",
                "alphafill_diff": round(cu_dist - info["af3_cu_dist"], 2) if cu_dist else "",
                "m3d_diff": round(info["m3d_cu_dist"] - info["af3_cu_dist"], 2),
            }
            results.append(row)
            if cu_dist:
                print("ok, %d Cu, best pair=%.2f A (AF3=%.2f)" % (
                    len(cu_atoms), cu_dist, info["af3_cu_dist"]))
            else:
                print("ok, %d Cu (no pair)" % len(cu_atoms))

        except Exception as e:
            print("error: %s" % e)
            results.append({"accession": acc, "alphafill_status": "error_%s" % type(e).__name__})

        time.sleep(args.delay)

    fieldnames = [
        "accession", "alphafill_status", "af3_cu_dist", "m3d_cu_dist",
        "alphafill_n_cu", "alphafill_cu_dist", "alphafill_diff", "m3d_diff",
    ]
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in results:
            w.writerow({k: r.get(k, "") for k in fieldnames})

    n_ok = sum(1 for r in results if r.get("alphafill_cu_dist"))
    print("\nDone: %d submitted, %d with Cu pairs" % (len(results), n_ok))
    if n_ok:
        diffs = [abs(float(r["alphafill_diff"])) for r in results if r.get("alphafill_diff")]
        m3d_diffs = [abs(float(r["m3d_diff"])) for r in results if r.get("m3d_diff")]
        print("AlphaFill |diff| from AF3: median=%.2f, mean=%.2f" % (
            sorted(diffs)[len(diffs)//2], sum(diffs)/len(diffs)))
        print("Metal3D   |diff| from AF3: median=%.2f, mean=%.2f" % (
            sorted(m3d_diffs)[len(m3d_diffs)//2], sum(m3d_diffs)/len(m3d_diffs)))


if __name__ == "__main__":
    main()
