#!/usr/bin/env python3
"""Print the convex hulls already found, one per line, for `webb-toroids --new-hulls FILE`.

A hull is identified by what results.tsv records about it: symmetry group order, circumradius,
inradius and the list of hull face sizes. Every (R)(A)(Q)(T) row of every runs/*/results.tsv
(and runs/*/combined/results.tsv) counts. Usage: known_hulls.py [RUNS_DIR] > hulls.txt
"""
import csv
import glob
import os
import sys


def main():
    runs = sys.argv[1] if len(sys.argv) > 1 else "runs"
    hulls = set()
    for f in glob.glob(os.path.join(runs, "*", "results.tsv")) + glob.glob(os.path.join(runs, "*", "combined", "results.tsv")):
        try:
            for r in csv.DictReader(open(f), delimiter="\t"):
                if r.get("class") == "RAQT":
                    hulls.add("\t".join((r["group"], r["rcirc"], r["rin"], r["hull_faces"])))
        except (OSError, KeyError):
            continue
    for h in sorted(hulls):
        print(h)


if __name__ == "__main__":
    main()
