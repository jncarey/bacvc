#!/usr/bin/env python3
"""
Generate a per-population SNP alignment FASTA from pop_tables/snp.gt.tsv.

Each row of snp.gt.tsv is a SNP/MNP site with columns:
    CHR  POS  TYPE  REF  ALT  <isolate_1>  <isolate_2>  ...
where each isolate cell is "0" (REF), "1" (ALT), or "-" (missing).

For each isolate, this script concatenates a per-row string:
    "0" -> REF string for that row
    "1" -> ALT string for that row
    "-" -> "N" * len(REF)
Rows where len(REF) != len(ALT) are skipped with a warning (snp/mnp from
freebayes are always aligned length).

Usage:
    python 05_generate_snp_alignment.py --pop <POP> --gt-tab <PATH> --out <PATH>
"""

import argparse
import csv
import sys
from pathlib import Path


def write_snp_alignment(gt_tab_path, out_path, pop_name):
    gt_tab_path = Path(gt_tab_path)
    out_path = Path(out_path)

    with open(gt_tab_path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if header[:5] != ["CHR", "POS", "TYPE", "REF", "ALT"]:
            print(f"ERROR: unexpected header in {gt_tab_path}: {header[:5]}",
                  file=sys.stderr)
            sys.exit(1)
        isolates = header[5:]

        sequences = ["" for _ in isolates]
        n_sites = 0
        n_skipped = 0
        seq_len = 0
        for row in reader:
            ref = row[3]
            alt = row[4]
            if len(ref) != len(alt):
                print(f"  WARNING: skipping row CHR={row[0]} POS={row[1]} "
                      f"REF={ref} ALT={alt} (length mismatch)", file=sys.stderr)
                n_skipped += 1
                continue
            cells = row[5:]
            gap = "N" * len(ref)
            for i, cell in enumerate(cells):
                if cell == "0":
                    sequences[i] += ref
                elif cell == "1":
                    sequences[i] += alt
                else:
                    sequences[i] += gap
            n_sites += 1
            seq_len += len(ref)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        for iso, seq in zip(isolates, sequences):
            fh.write(f">{iso}\n{seq}\n")

    print(f"[{pop_name}] Wrote {out_path} "
          f"({len(isolates)} isolates × {seq_len} bp from {n_sites} sites"
          + (f"; skipped {n_skipped} length-mismatch rows" if n_skipped else "")
          + ")")


def main():
    parser = argparse.ArgumentParser(
        description="Generate a per-population SNP alignment FASTA from snp.gt.tsv."
    )
    parser.add_argument("--pop", required=True, help="Population name (e.g. IA04)")
    parser.add_argument("--gt-tab", required=True,
                        help="Path to pop_tables/snp.gt.tsv")
    parser.add_argument("--out", required=True,
                        help="Output FASTA path (e.g. aln/IA04.snps.aln.fa)")
    args = parser.parse_args()
    write_snp_alignment(args.gt_tab, args.out, args.pop)


if __name__ == "__main__":
    main()
