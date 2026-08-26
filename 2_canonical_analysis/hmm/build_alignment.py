#!/usr/bin/env python3
"""Regenerate the PF00264 profile alignment of the query sequences, and write it in
the compact form the analyses actually consume.

Why this exists
---------------
`run_hmmalign.sh` produces `all_hmmalign.afa`, which is 782 MB and therefore not
deposited; `query.fasta` was also excluded from the repository as bulk data. Between
them that made novelty_pipeline.py stages G and I (Table S4's `hmm_agreement.tsv`)
and the non-canonical six-position vector comparison impossible for a reader to
re-run. `query.fasta.gz` is now deposited, and this script turns it back into the
alignment with no external tools and no 782 MB intermediate.

What it writes
--------------
`hmm_match_columns.tsv.gz`: one row per sequence, `accession <TAB> 214 characters`,
one character per PF00264.26 match state -- upper-case residue, or `-` where the
sequence has a deletion. Insert columns are dropped, which is the only part of the
alignment no downstream analysis here uses. About 7 MB uncompressed, 2 MB on disk,
and it is the complete input for:
  * novelty_pipeline.py stages G and I  (`--hmm-cols`)
  * 3_noncanonical_analysis/hmm_vector_check/  (the 504-of-1,060 claim)

hmmalign aligns each sequence to the profile independently, so the match-state
assignment for a sequence does not depend on which other sequences are in the run.
That is what lets this stream in chunks (bounded memory) and still give exactly the
columns a single 32,753-sequence `hmmalign` call gives. Pass `--chunk 0` to align in
one pass and `--afa-out FILE` to also write the full .afa in hmmalign's format.

Requires pyhmmer (`pip install pyhmmer`), which embeds the same HMMER3 code as the
hmmalign binary; if pyhmmer is missing and an `hmmalign` binary is on PATH, use
run_hmmalign.sh instead and pass its output here with `--from-afa`.

Usage
-----
    python3 build_alignment.py                       # -> hmm_match_columns.tsv.gz
    python3 build_alignment.py --chunk 0 --afa-out all_hmmalign.afa
    python3 build_alignment.py --from-afa all_hmmalign.afa
"""
import argparse
import gzip
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def find_fasta(explicit):
    if explicit:
        return explicit
    for cand in ("query.fasta.gz", "query.fasta"):
        p = os.path.join(HERE, cand)
        if os.path.exists(p):
            return p
    sys.exit("no query.fasta(.gz) next to this script; pass --fasta")


def read_fasta(path):
    name, buf = None, []
    with open_maybe_gz(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(buf)
                name, buf = line[1:].split("|")[0].strip(), []
            else:
                buf.append(line.strip())
    if name:
        yield name, "".join(buf)


def match_columns_from_afa(path, n_match):
    """Extract the match-state string from an hmmalign .afa: match columns are the
    upper-case / '-' positions, insert columns are lower-case / '.'."""
    def collapse(al):
        return "".join(c for c in al if c == "-" or c.isupper())
    name, buf = None, []
    with open_maybe_gz(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    yield name, collapse("".join(buf))
                name, buf = line[1:].split()[0].split("|")[0], []
            else:
                buf.append(line)
    if name:
        yield name, collapse("".join(buf))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hmm", default=os.path.join(HERE, "PF00264.hmm"))
    ap.add_argument("--fasta", default=None, help="query.fasta or query.fasta.gz")
    ap.add_argument("--out", default=os.path.join(HERE, "hmm_match_columns.tsv.gz"))
    ap.add_argument("--afa-out", default=None,
                    help="also write the full .afa (requires --chunk 0)")
    ap.add_argument("--from-afa", default=None,
                    help="skip alignment; derive the columns from an existing .afa")
    ap.add_argument("--chunk", type=int, default=2000,
                    help="sequences per alignment batch; 0 = single pass (needs ~1 GB)")
    ap.add_argument("--ref-acc", default="B2ZB02",
                    help="reference sequence whose residue numbering is mapped to "
                         "match states (PmTYR)")
    ap.add_argument("--ref-map-out", default=os.path.join(HERE, "reference_map.tsv"),
                    help="reference residue number -> match state table")
    a = ap.parse_args()

    if a.from_afa:
        n_written = 0
        with gzip.open(a.out, "wt") as out:
            out.write("accession\tmatch_columns\n")
            for name, cols in match_columns_from_afa(a.from_afa, None):
                out.write(f"{name}\t{cols}\n")
                n_written += 1
        print(f"wrote {a.out} from {a.from_afa} ({n_written:,} sequences)")
        return

    try:
        import pyhmmer
    except ImportError:
        sys.exit("pyhmmer not installed (pip install pyhmmer). Alternatively run "
                 "run_hmmalign.sh and re-run this script with --from-afa.")

    fasta = find_fasta(a.fasta)
    alphabet = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.plan7.HMMFile(a.hmm) as fh:
        hmm = fh.read()
    name = hmm.name.decode() if isinstance(hmm.name, bytes) else hmm.name
    acc = hmm.accession.decode() if isinstance(hmm.accession, bytes) else hmm.accession
    print(f"profile {name} {acc}  M={hmm.M} match states")
    print(f"sequences from {os.path.basename(fasta)}")

    if a.afa_out and a.chunk:
        sys.exit("--afa-out needs --chunk 0 (a chunked run has per-chunk insert columns)")

    def align(batch):
        seqs = [pyhmmer.easel.TextSequence(name=n.encode(), sequence=s).digitize(alphabet)
                for n, s in batch]
        msa = pyhmmer.hmmalign(hmm, seqs, trim=False)
        names = [n.decode() if isinstance(n, bytes) else n for n in msa.names]
        rf = msa.reference
        rf = rf.decode() if isinstance(rf, bytes) else (rf or "")
        if rf:
            cols = [i for i, c in enumerate(rf) if c != "."]
        else:
            cols = [i for i, c in enumerate(msa.alignment[0])
                    if c != "." and not c.islower()]
        assert len(cols) == hmm.M, (len(cols), hmm.M)
        return names, msa, cols

    n_written = 0
    ref_map = []
    with gzip.open(a.out, "wt") as out:
        out.write("accession\tmatch_columns\n")
        batch = []
        pending = []

        def flush():
            nonlocal n_written, batch
            if not batch:
                return
            names, msa, cols = align(batch)
            colset = set(cols)
            for nm, al in zip(names, msa.alignment):
                out.write(f"{nm}\t{''.join(al[i] for i in cols).upper()}\n")
                n_written += 1
                if nm == a.ref_acc and not ref_map:
                    # walk the reference's own aligned row once: residue number ->
                    # match state (None when the residue sits in an insert column)
                    num, ms = 0, 0
                    for i, c in enumerate(al):
                        is_match = i in colset
                        if is_match:
                            ms += 1
                        if c in "-.":
                            continue
                        num += 1
                        ref_map.append((num, c.upper(), ms if is_match else ""))
            if a.afa_out:
                pending.append((names, msa))
            batch = []
            print(f"  {n_written:,} aligned", end="\r", flush=True)

        for rec in read_fasta(fasta):
            batch.append(rec)
            if a.chunk and len(batch) >= a.chunk:
                flush()
        flush()

    print(f"\nwrote {a.out} ({n_written:,} sequences x {hmm.M} match states)")
    if ref_map:
        with open(a.ref_map_out, "w") as fh:
            fh.write("resnum\tresidue\tmatch_state\n")
            for num, res, ms in ref_map:
                fh.write(f"{num}\t{res}\t{ms}\n")
        n_ins = sum(1 for _, _, ms in ref_map if ms == "")
        print(f"wrote {a.ref_map_out} ({a.ref_acc}: {len(ref_map)} residues, "
              f"{n_ins} in insert columns)")
    elif a.ref_acc:
        print(f"WARNING: reference {a.ref_acc} not found in {os.path.basename(fasta)}; "
              f"no reference map written")
    if a.afa_out and pending:
        names, msa = pending[0]
        with open(a.afa_out, "w") as fh:
            for nm, al in zip(names, msa.alignment):
                fh.write(f">{nm}\n{al}\n")
        print(f"wrote {a.afa_out}")


if __name__ == "__main__":
    main()
