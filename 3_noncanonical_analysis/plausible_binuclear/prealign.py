#!/usr/bin/env python3
"""Pre-align the 66 plausible-binuclear CIFs onto PmTYR and bake the transformed
coordinates back into each CIF, so inspect_plausible_binuclear.pml can load them
already overlaid (no per-structure `super` at load -> fast PyMOL). PmTYR itself is
the reference frame and is left untouched. Run once after (re)fetching the CIFs:

    pymol -ckq prealign.py        # headless

Raw (un-aligned) coordinates remain recoverable from the AF3 model CIFs; these
labeled copies are inspection artifacts whose only purpose is the PmTYR overlay.
"""
import os, glob
from pymol import cmd

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)

cmd.set("retain_order", 1)
cmd.load("PmTYR_B2ZB02.cif", "PmTYR")

cifs = sorted(f for f in glob.glob("*.cif") if f != "PmTYR_B2ZB02.cif")
print(f"{len(cifs)} structures to pre-align", flush=True)

nbad = 0
for fn in cifs:
    obj = fn[:-4]
    cmd.load(fn, obj)
    try:
        rmsd = cmd.super(obj, "PmTYR")[0]
    except Exception as e:
        print("super FAILED:", obj, e, flush=True); nbad += 1; rmsd = -1.0
    cmd.save(fn, obj)          # overwrite with transformed coords
    cmd.delete(obj)
    print(f"{obj:34} rmsd={rmsd:.2f}", flush=True)

print(f"done, {nbad} failures", flush=True)
