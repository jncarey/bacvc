#!/usr/bin/env bash
#
# normalize_annotate_vcf.sh
#
# Normalize + snpEff-annotate a raw freebayes VCF into a per-isolate
# snps.norm.annot.vcf. Standalone, parameterized twin of config.sh's
# normalize_annotate_vcf() bash function, for callers (e.g. Snakemake rules)
# that don't source config.sh's globals.
#
# IMPORTANT: no FMT/DP filter is applied here. Low-coverage variants are
# retained so the DP decision is made once, downstream (the +setGT
# FMT/DP<VT_MIN_DP mask plus the mosdepth MASK_MIN_DP mask). Only
# QUAL>=filt_min_qual is applied. FORMAT/DP is kept in the output so that
# downstream step can still gate on it.
#
# Must run in an environment with bcftools, vt, and snpEff already on PATH
# (e.g. via `conda run -n vc` or a Snakemake conda: directive). Writes
# atomically: a mid-pipe failure under `set -o pipefail` leaves the existing
# output intact rather than a truncated VCF that a downstream merge would
# silently pick up.
#
# Usage:
#   normalize_annotate_vcf.sh <raw_vcf> <ref_fa> <ref_name> <out_vcf> \
#       <snpeff_config> <snpeff_datadir> <filt_min_qual>

set -euo pipefail

if [[ $# -ne 7 ]]; then
    echo "Usage: $0 <raw_vcf> <ref_fa> <ref_name> <out_vcf> <snpeff_config> <snpeff_datadir> <filt_min_qual>" >&2
    exit 1
fi

raw_vcf="$1"
ref_fa="$2"
ref_name="$3"
out_vcf="$4"
snpeff_config="$5"
snpeff_datadir="$6"
filt_min_qual="$7"

tmp="${out_vcf}.tmp.$$"

bcftools view --include "QUAL>=${filt_min_qual}" "${raw_vcf}" \
    | vt normalize -r "${ref_fa}" - \
    | bcftools annotate \
        --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' \
    | snpEff ann \
        -noLog -noStats \
        -no-downstream -no-upstream -no-utr \
        -c "${snpeff_config}" \
        -dataDir "${snpeff_datadir}" \
        "${ref_name}" \
    > "${tmp}"

mv "${tmp}" "${out_vcf}"
