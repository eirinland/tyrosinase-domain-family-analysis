#!/bin/bash
# Regenerate the Gly46_Y per-locus neighbourhoods in THIS directory.
# The fetch scripts read/write in their own parent dir, so they are copied in
# here first; running them in 4_genome_neighbourhood/ would overwrite the
# non-canonical GNA outputs.
set -eu
HERE="$(cd "$(dirname "$0")" && pwd)"
cp "$HERE/../2_fetch_genome_context.py" "$HERE/../3_fetch_neighbourhoods.py" "$HERE/"
cd "$HERE"
python3 2_fetch_genome_context.py
python3 3_fetch_neighbourhoods.py --flank 10 --workers 4
rm -f 2_fetch_genome_context.py 3_fetch_neighbourhoods.py
echo "done: $(( $(wc -l < neighbourhoods.tsv) - 1 )) flanking genes, $(cut -f1 neighbourhoods.tsv | tail -n +2 | sort -u | wc -l) loci"
