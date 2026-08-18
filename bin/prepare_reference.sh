#!/usr/bin/env bash
#
# prepare_reference.sh
#
# Library entry point for indexing one reference genome for variant calling,
# with no config.sh dependency: every path is an explicit argument. Leaner
# than 02_prepare_reference.sh's per-population loop body in one respect:
# it doesn't copy/symlink a snpEff config or GFF into the output directory,
# since callers already have those paths from build_snpeff_db.sh and pass
# them explicitly to normalize_annotate_vcf.sh / 04_generate_variant_tables.py
# -- no script actually reads a config or GFF back out of this directory.
#
# Usage:
#   prepare_reference.sh <ref_fna> <out_dir> [minimap2_selfalign_opts]
#
# Output (under <out_dir>/):
#   ref.fa, ref.fa.fai, ref.fa.bwt (+ other bwa index files)
#   ref.repeats.bed - self-alignment repeat mask (see "Population-level site
#                     exclusions" in README.md)

set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: $0 <ref_fna> <out_dir> [minimap2_selfalign_opts]" >&2
    exit 1
fi

REF_FNA="$1"
OUT_DIR="$2"
MINIMAP2_SELFALIGN_OPTS="${3:--D -P -w19 -m200}"

if [[ ! -f "${REF_FNA}" ]]; then
    echo "ERROR: Reference FASTA not found: ${REF_FNA}" >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"

cp "${REF_FNA}" "${OUT_DIR}/ref.fa"
echo "Copied ${REF_FNA} -> ${OUT_DIR}/ref.fa"

echo "Running bwa index..."
bwa index "${OUT_DIR}/ref.fa"

echo "Running samtools faidx..."
samtools faidx "${OUT_DIR}/ref.fa"

# Self-alignment repeat BED: minimap2-align the reference to itself to find
# repetitive / self-similar regions, emit as a merged BED. Variants in these
# regions are excluded from the canonical merged VCF downstream (see
# "Population-level site exclusions" in README.md). Built from ref.fa so
# contig names match the VCFs/GFF by construction.
echo "Building self-alignment repeat BED (minimap2 ${MINIMAP2_SELFALIGN_OPTS})..."
minimap2 ${MINIMAP2_SELFALIGN_OPTS} "${OUT_DIR}/ref.fa" "${OUT_DIR}/ref.fa" 2>/dev/null \
    | cut -f 1,3,4 \
    | sort -k1,1 -k2,2n \
    | bedtools slop -i stdin -g "${OUT_DIR}/ref.fa.fai" -r 1 -l 0 \
    | bedtools merge -i stdin \
    > "${OUT_DIR}/ref.repeats.bed"

REPEAT_N=$(wc -l < "${OUT_DIR}/ref.repeats.bed")
REPEAT_BP=$(awk '{sum += $3 - $2} END {print sum + 0}' "${OUT_DIR}/ref.repeats.bed")
REF_BP=$(awk '{sum += $2} END {print sum + 0}' "${OUT_DIR}/ref.fa.fai")
REPEAT_PCT=$(awk -v m="${REPEAT_BP}" -v r="${REF_BP}" 'BEGIN {printf "%.2f", (r ? 100 * m / r : 0)}')
echo "Repeat BED: ${OUT_DIR}/ref.repeats.bed (${REPEAT_N} intervals, ${REPEAT_BP}/${REF_BP} bp = ${REPEAT_PCT}%)"
