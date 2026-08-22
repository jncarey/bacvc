#!/usr/bin/env python3
"""Write a headerless isolate<TAB>population TSV (pixy's expected
--populations format), grouping isolates by visit label parsed from their
own names via pipeline_helpers.group_isolates_by_visit.

--collapse-visits/--collapse-as let the caller merge specific visit labels
into one canonical population label -- e.g. a patient whose "V1A"/"V1B"
labels are actually the same specimen arrayed/sequenced twice (not two
distinct visits), where without collapsing, whichever replicate a caller's
own dedup step happened to keep would produce a meaningless artificial
split. This project policy knowledge (which patients/labels need collapsing)
belongs to the caller, not this generic script.

Usage:
    populations_file.py --patient NAME --isolates ISO1 ISO2 ... --out PATH
        [--collapse-visits V1A V1B --collapse-as V1]
"""

import argparse
from pathlib import Path

from pipeline_helpers import group_isolates_by_visit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patient", required=True)
    parser.add_argument("--isolates", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--collapse-visits", nargs="+", default=[],
                         help="Visit labels to merge into one canonical population label")
    parser.add_argument("--collapse-as", default=None,
                         help="Canonical population label for --collapse-visits")
    args = parser.parse_args()

    groups = group_isolates_by_visit(args.patient, args.isolates)
    collapse_set = set(args.collapse_visits)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        for visit in sorted(groups):
            label = args.collapse_as if visit in collapse_set else visit
            for iso in groups[visit]:
                f.write(f"{iso}\t{label}\n")


if __name__ == "__main__":
    main()
