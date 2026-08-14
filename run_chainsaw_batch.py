"""
Batch Chainsaw domain segmentation for PPO structures.
Processes structures sequentially (model loaded once) with per-structure output.
"""
import csv, os, sys, shutil, time, logging
from pathlib import Path

CHAINSAW_DIR = "/cluster/home/eirinlandsem/Super_reference_pipeline/chainsaw"
sys.path.insert(0, CHAINSAW_DIR)

os.environ["LOGLEVEL"] = "WARNING"
logging.basicConfig(level=logging.WARNING, stream=sys.stderr,
                    format="%(asctime)s | %(levelname)s | %(message)s")

from get_predictions import predict, load_model

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--cifs", required=True)
    parser.add_argument("--accessions", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--scratch", default="/tmp/chainsaw_scratch")
    args = parser.parse_args()

    with open(args.accessions) as f:
        target_accs = set(line.strip() for line in f if line.strip())
    print(f"Target accessions: {len(target_accs)}", flush=True)

    cif_dir = Path(args.cifs)
    acc_to_cif = {}
    for p in sorted(cif_dir.glob("*.cif")):
        acc = p.name.split("_taxID_")[0]
        if acc in target_accs:
            acc_to_cif[acc] = str(p)
    print(f"CIF files matched: {len(acc_to_cif)}", flush=True)

    scratch = Path(args.scratch)
    scratch.mkdir(parents=True, exist_ok=True)

    model_dir = os.path.join(CHAINSAW_DIR, "saved_models", "model_v3")
    print("Loading model...", flush=True)
    model = load_model(model_dir=model_dir)
    print("Model loaded.", flush=True)

    outpath = Path(args.output)
    fieldnames = ["accession", "nres", "ndom", "chopping", "confidence"]
    
    # Resume from existing output
    done_accs = set()
    if outpath.exists():
        with open(outpath) as f:
            for r in csv.DictReader(f):
                done_accs.add(r["accession"])
        print(f"Resuming: {len(done_accs)} already done", flush=True)
    
    write_header = not outpath.exists() or len(done_accs) == 0
    outfile = open(outpath, "a", newline="")
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    if write_header:
        writer.writeheader()

    todo = sorted(set(acc_to_cif.keys()) - done_accs)
    print(f"To process: {len(todo)}", flush=True)

    t0 = time.time()
    for i, acc in enumerate(todo):
        cif_path = acc_to_cif[acc]
        local_cif = scratch / Path(cif_path).name
        try:
            shutil.copy2(cif_path, local_cif)
            result = predict(model, str(local_cif))
            writer.writerow({
                "accession": acc,
                "nres": result.nres,
                "ndom": result.ndom,
                "chopping": result.chopping if result.chopping else "NULL",
                "confidence": f"{result.confidence:.3f}" if result.confidence else "",
            })
        except Exception as e:
            writer.writerow({
                "accession": acc,
                "nres": 0,
                "ndom": -1,
                "chopping": f"ERROR:{str(e)[:80]}",
                "confidence": "",
            })
        finally:
            for f in scratch.glob("*"):
                try:
                    f.unlink()
                except Exception:
                    pass
        
        if (i + 1) % 100 == 0:
            outfile.flush()
            elapsed = time.time() - t0
            rate = (i + 1) / elapsed
            eta = (len(todo) - i - 1) / rate
            print(f"  {i+1}/{len(todo)} ({rate:.1f}/s, ETA {eta/60:.0f}min)", flush=True)

    outfile.close()
    elapsed = time.time() - t0
    print(f"\nDone. {len(todo)} processed in {elapsed/60:.1f}min", flush=True)

if __name__ == "__main__":
    main()
