#!/usr/bin/env python3
"""Write a headerless isolate<TAB>population TSV (pixy's expected
--populations format), grouping isolates by visit label parsed from their
own names via pipeline_helpers.group_isolates_by_visit.

Usage:
    populations_file.py --patient NAME --isolates ISO1 ISO2 ... --out PATH
"""

import argparse
from pathlib import Path

from pipeline_helpers import group_isolates_by_visit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient", required=True)
    parser.add_argument("--isolates", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    groups = group_isolates_by_visit(args.patient, args.isolates)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for visit in sorted(groups):
            for iso in groups[visit]:
                f.write(f"{iso}\t{visit}\n")


if __name__ == "__main__":
    main()
