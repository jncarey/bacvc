#!/usr/bin/env python3
"""
04_generate_variant_vcf.py

Repeat-exclusion + mosdepth-mask pass that builds the canonical multi-sample
VCF. Its output is the single source of truth for both the TreeTime VCF and the
population variant tables.

Input: a merged VCF produced with `bcftools merge --missing-to-ref` over ALL
variant types (snp/mnp/ins/del/complex), so samples without a row at a given
site default to 0/0. Other possible GT values:
    1/1, 2/2, ...   confident ALT calls (kept from per-sample filtered VCFs)
    ./.             upstream setGT-masked cells (low QUAL/DP/AO frac/het)

This pass walks every record and:
  - drops the whole record if its REF span overlaps a self-alignment repeat
    interval (--repeat-bed); applied once here so every downstream artifact
    shares an identical excluded set;
  - for each 0/0 cell, looks up the sample's mosdepth coverage at that position
    and demotes to ./. when coverage was below mask_min_dp. Cells that are not
    0/0 are never touched;
  - drops the whole record as a high-missingness population-level site when, after
    masking, at least --max-missing-frac of isolates are missing ("-") AND the
    record's INFO/TYPE is purely snp/mnp. Indel/complex and multiallelic mixed
    records (any allele not snp/mnp) are never eligible, so indels are never
    affected; this propagates to both the TreeTime VCF and the snp tables, which
    derive from this canonical VCF. A high-missing snp sharing a multiallelic
    record with an indel/complex allele therefore survives (uncommon).

If --missingness-tsv is given, every record reaching the missingness check (i.e.
after repeat exclusion, after masking) emits one QC row, so the file is the full
per-site "-" distribution and flags exactly what the threshold removed.

Per-sample coverage is held as a bit-packed callable mask (1 bit per genome
position), so all sample masks together fit in ~1 GB even at N=1000 isolates
on a 6.8 Mb bacterial genome.

Output is plain VCF on stdout; the shell wrapper bgzips + tabix-indexes.

Usage:
    python 04_generate_variant_vcf.py \\
        --merged <merged.raw.vcf.gz> \\
        --var-dir <results/populations/POP/variants> \\
        --mask-min-dp <int> \\
        > <merged.masked.vcf>
"""

import argparse
import gzip
import re
import sys
from datetime import datetime
from pathlib import Path

from pipeline_helpers import (
    build_callable_mask,
    classify_isolates,
    is_callable,
    is_missing_gt,
    load_repeat_bed,
    overlaps_repeat,
)


CONTIG_RE = re.compile(r'^##contig=<.*?ID=([^,>]+).*?length=(\d+)')


def log(msg):
    print(f"[{datetime.now():%H:%M:%S}] {msg}", file=sys.stderr, flush=True)


def opener(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def info_type_tokens(info):
    """Return the comma-separated INFO/TYPE tokens (per-allele types) as a list,
    or [] when absent. A merged multiallelic record carries one token per ALT."""
    for kv in info.split(";"):
        if kv.startswith("TYPE="):
            return kv[len("TYPE="):].split(",")
    return []


def parse_header(merged_path):
    """Return (chrom_lengths, isolates) from the merged VCF header."""
    chrom_lengths = {}
    isolates = []
    with opener(merged_path) as fh:
        for line in fh:
            if line.startswith("##contig"):
                m = CONTIG_RE.match(line.rstrip("\n"))
                if m:
                    chrom_lengths[m.group(1)] = int(m.group(2))
            elif line.startswith("#CHROM"):
                isolates = line.rstrip("\n").split("\t")[9:]
                break
            elif not line.startswith("#"):
                break
    return chrom_lengths, isolates


def main():
    parser = argparse.ArgumentParser(
        description="Demote 0/0 -> ./. in a merged VCF where mosdepth coverage "
                    "is below threshold."
    )
    parser.add_argument("--merged", required=True,
                        help="Merged multi-sample VCF (.vcf or .vcf.gz)")
    parser.add_argument("--var-dir", required=True,
                        help="Directory holding <iso>/<iso>.per-base.bed.gz per sample")
    parser.add_argument("--mask-min-dp", type=int, required=True,
                        help="Minimum mosdepth coverage to keep 0/0; below -> ./.")
    parser.add_argument("--repeat-bed", default=None,
                        help="Self-alignment repeat BED; any record overlapping a "
                             "repeat interval is dropped. Omit to disable.")
    parser.add_argument("--max-missing-frac", type=float, default=None,
                        help="Drop a snp/mnp record when (after masking) at least "
                             "this fraction of isolates are missing ('-'). Omit to "
                             "disable the drop (QC is still emitted).")
    parser.add_argument("--pop", default=None,
                        help="Population name; only used to label PR/PO columns in "
                             "the missingness QC TSV.")
    parser.add_argument("--missingness-tsv", default=None,
                        help="Write the per-site '-' distribution to this path "
                             "(one row per record reaching the missingness check).")
    args = parser.parse_args()

    merged_path = Path(args.merged)
    var_dir = Path(args.var_dir)
    mask_min_dp = args.mask_min_dp
    max_missing_frac = args.max_missing_frac

    repeat_idx = None
    if args.repeat_bed:
        repeat_idx = load_repeat_bed(args.repeat_bed)
        n_iv = sum(len(s) for s, _ in repeat_idx.values())
        log(f"Loaded repeat BED {args.repeat_bed}: "
            f"{n_iv} intervals across {len(repeat_idx)} contig(s)")

    chrom_lengths, isolates = parse_header(merged_path)
    if not chrom_lengths:
        print("ERROR: no ##contig lines found in merged VCF header",
              file=sys.stderr)
        sys.exit(1)
    if not isolates:
        print("ERROR: no #CHROM line in merged VCF", file=sys.stderr)
        sys.exit(1)

    total_bp = sum(chrom_lengths.values())
    per_sample_kb = (total_bp + 7) // 8 / 1024
    log(f"Reference: {len(chrom_lengths)} contigs, {total_bp:,} bp")
    log(f"Building callable masks for {len(isolates)} samples "
        f"(~{per_sample_kb:.0f} KB each, "
        f"~{per_sample_kb * len(isolates) / 1024:.1f} MB total)")

    callable_masks = {}
    # Progress every ~5% of samples (and always at least every 50).
    progress_step = max(1, min(50, len(isolates) // 20))
    for i_iso, iso in enumerate(isolates, start=1):
        bed = var_dir / iso / f"{iso}.per-base.bed.gz"
        if not bed.exists():
            log(f"WARNING: no depth file for {iso} ({bed}); "
                f"its 0/0 cells will not be masked")
            callable_masks[iso] = None
            continue
        callable_masks[iso] = build_callable_mask(bed, chrom_lengths, mask_min_dp)
        if i_iso % progress_step == 0 or i_iso == len(isolates):
            log(f"  built {i_iso}/{len(isolates)} masks")

    # Population-level missingness filter + QC setup. PR/PO column indices are
    # precomputed so each record's per-timepoint missing count is a cheap lookup.
    n_total = len(isolates)
    if args.pop:
        pr_isolates, po_isolates = classify_isolates(args.pop, isolates, log=log)
    else:
        pr_isolates, po_isolates = [], []
    pr_idx = [i for i, iso in enumerate(isolates) if iso in set(pr_isolates)]
    po_idx = [i for i, iso in enumerate(isolates) if iso in set(po_isolates)]
    if max_missing_frac is not None:
        log(f"Missingness filter: drop snp/mnp sites with >= "
            f"{max_missing_frac:.2f} of {n_total} isolates missing")
    else:
        log("Missingness filter: disabled (QC only)")

    qc_fh = None
    if args.missingness_tsv:
        qc_fh = open(args.missingness_tsv, "w")
        qc_fh.write("\t".join([
            "CHR", "POS", "REF", "ALT", "TYPE", "N_TOTAL", "N_MISSING",
            "FRAC_MISSING", "PR_TOTAL", "PR_MISSING", "PO_TOTAL", "PO_MISSING",
            "EXCLUDED",
        ]) + "\n")
        log(f"Writing per-site '-' distribution to {args.missingness_tsv}")

    log("Walking merged VCF and applying mask")
    n_masked = 0
    n_kept_ref = 0
    n_alt = 0
    n_already_missing = 0
    n_missingness_excluded = 0
    n_rows = 0
    n_repeat_excluded = 0
    out = sys.stdout

    with opener(merged_path) as fh:
        for line in fh:
            if line.startswith("#"):
                out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            # Repeat exclusion: drop the whole record if its REF span overlaps a
            # self-alignment repeat interval. Applied once here so every artifact
            # derived from this canonical VCF shares an identical excluded set.
            if repeat_idx is not None and overlaps_repeat(repeat_idx, chrom, pos, len(ref)):
                n_repeat_excluded += 1
                continue
            fmt_keys = fields[8].split(":")
            try:
                gt_idx = fmt_keys.index("GT")
            except ValueError:
                out.write(line)
                continue

            for i, iso in enumerate(isolates):
                cell = fields[9 + i]
                subfields = cell.split(":")
                gt = subfields[gt_idx] if gt_idx < len(subfields) else "./."

                if gt in ("./.", "."):
                    n_already_missing += 1
                    continue
                if gt not in ("0/0", "0|0", "0"):
                    # 1/1, 2/2, etc. -- trusted ALT calls; never demoted.
                    n_alt += 1
                    continue

                masks = callable_masks.get(iso)
                if masks is None:
                    n_kept_ref += 1
                    continue
                if not is_callable(masks, chrom, pos):
                    subfields[gt_idx] = "./."
                    fields[9 + i] = ":".join(subfields)
                    n_masked += 1
                else:
                    n_kept_ref += 1

            # Population-level missingness: count "-" cells from the *final*
            # (post-mask) genotypes, decide whether to drop the site, and emit QC.
            def _cell_missing(col):
                sub = fields[9 + col].split(":")
                return is_missing_gt(sub[gt_idx] if gt_idx < len(sub) else "./.")

            n_missing = sum(1 for i in range(n_total) if _cell_missing(i))
            frac_missing = n_missing / n_total if n_total else 0.0
            type_tokens = info_type_tokens(fields[7])
            is_snp_mnp = bool(type_tokens) and all(
                t in ("snp", "mnp") for t in type_tokens)
            excluded = (
                is_snp_mnp
                and max_missing_frac is not None
                and frac_missing >= max_missing_frac
            )

            if qc_fh is not None:
                pr_missing = sum(1 for i in pr_idx if _cell_missing(i))
                po_missing = sum(1 for i in po_idx if _cell_missing(i))
                qc_fh.write("\t".join([
                    chrom, str(pos), ref, fields[4],
                    ",".join(type_tokens), str(n_total), str(n_missing),
                    f"{frac_missing:.4f}", str(len(pr_idx)), str(pr_missing),
                    str(len(po_idx)), str(po_missing),
                    "1" if excluded else "0",
                ]) + "\n")

            if excluded:
                n_missingness_excluded += 1
                continue

            out.write("\t".join(fields) + "\n")
            n_rows += 1
            if n_rows % 10000 == 0:
                log(f"  processed {n_rows:,} sites")

    if qc_fh is not None:
        qc_fh.close()

    log(f"Done. Sites written:              {n_rows}")
    log(f"  records dropped (repeat region):{n_repeat_excluded}")
    log(f"  records dropped (high-missing snp site): {n_missingness_excluded}")
    log(f"  0/0 -> ./. masked (low cov):    {n_masked}")
    log(f"  0/0 kept (covered REF):         {n_kept_ref}")
    log(f"  ALT calls left untouched:       {n_alt}")
    log(f"  ./. preserved (upstream-masked):{n_already_missing}")


if __name__ == "__main__":
    main()
