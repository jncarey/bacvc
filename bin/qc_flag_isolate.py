#!/usr/bin/env python3
"""Write one isolate's QC-flag row from fastp/flagstat/mosdepth output.

Thin CLI over pipeline_helpers.qc_flag_isolate -- see that function's
docstring for the two flags computed (low_coverage, mapping_failure).

Usage:
    qc_flag_isolate.py --isolate NAME --fastp-json PATH --flagstat PATH \\
        --mosdepth-summary PATH --min-mean-depth N --min-mapped-pct N \\
        --out PATH
"""

import argparse
from pathlib import Path

from pipeline_helpers import qc_flag_isolate

FIELDS = ["isolate", "bases_after_filter", "mean_depth", "mapped_pct",
          "low_coverage", "mapping_failure"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--isolate", required=True)
    parser.add_argument("--fastp-json", required=True)
    parser.add_argument("--flagstat", required=True)
    parser.add_argument("--mosdepth-summary", required=True)
    parser.add_argument("--min-mean-depth", type=float, required=True)
    parser.add_argument("--min-mapped-pct", type=float, required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    result = qc_flag_isolate(
        args.isolate, args.fastp_json, args.flagstat, args.mosdepth_summary,
        args.min_mean_depth, args.min_mapped_pct,
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as f:
        f.write("\t".join(FIELDS) + "\n")
        f.write("\t".join(str(result[k]) for k in FIELDS) + "\n")

    print(f"  QC: {args.isolate} mean_depth={result['mean_depth']:.2f} "
          f"mapped_pct={result['mapped_pct']:.2f} "
          f"low_coverage={result['low_coverage']} "
          f"mapping_failure={result['mapping_failure']}")


if __name__ == "__main__":
    main()
