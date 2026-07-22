#!/usr/bin/env bash
#
# 02_prepare_reference.sh
#
# Prepare a population-specific reference for variant calling:
#   - Copy bakta .fna to reference directory
#   - bwa index
#   - samtools faidx
#   - Symlink snpEff database
#
# Usage:
#   conda activate snippy
#   bash 02_prepare_reference.sh <POP> [POP ...]
#   # e.g.: bash 02_prepare_reference.sh IA01 IA02 UW03

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

if [[ $# -eq 0 ]]; then
    echo "Usage: $0 <POP> [POP ...]" >&2
    echo "Populations: ${!POP_TO_REF[*]}" >&2
    exit 1
fi

for POP in "$@"; do
    echo "=== Preparing reference for population ${POP} ==="

    REF_NAME="$(get_ref_name "${POP}")"
    REF_FNA="$(get_ref_fasta "${POP}")"
    REF_GFF="$(get_ref_gff "${POP}")"
    POP_DIR="$(get_pop_dir "${POP}")"
    REF_DIR="${POP_DIR}/reference"

    # Validate inputs
    if [[ ! -f "${REF_FNA}" ]]; then
        echo "ERROR: Reference FASTA not found: ${REF_FNA}" >&2
        exit 1
    fi
    if [[ ! -f "${REF_GFF}" ]]; then
        echo "ERROR: Reference GFF not found: ${REF_GFF}" >&2
        echo "  Run 01_build_snpeff_db.sh first" >&2
        exit 1
    fi

    mkdir -p "${REF_DIR}"

    # Copy reference FASTA (bakta .fna with contig_* names matching GFF3)
    if [[ ! -f "${REF_DIR}/ref.fa" ]]; then
        cp "${REF_FNA}" "${REF_DIR}/ref.fa"
        echo "  Copied ${REF_FNA} -> ref.fa"
    else
        echo "  ref.fa already present."
    fi

    # Copy snpEff genes.gff (flat CDS-only GFF built from bakta .gbff)
    if [[ ! -f "${REF_DIR}/ref.gff" ]]; then
        cp "${REF_GFF}" "${REF_DIR}/ref.gff"
        echo "  Copied ${REF_GFF} -> ref.gff"
    else
        echo "  ref.gff already present."
    fi

    # bwa index
    if [[ ! -f "${REF_DIR}/ref.fa.bwt" ]]; then
        echo "  Running bwa index..."
        bwa index "${REF_DIR}/ref.fa"
    else
        echo "  bwa index already present."
    fi

    # samtools faidx
    if [[ ! -f "${REF_DIR}/ref.fa.fai" ]]; then
        echo "  Running samtools faidx..."
        samtools faidx "${REF_DIR}/ref.fa"
    else
        echo "  samtools faidx already present."
    fi

    # Self-alignment repeat BED: minimap2-align the reference to itself to find
    # repetitive / self-similar regions, emit as a merged BED. Variants in these
    # regions are excluded from the canonical merged VCF in step 04. Built from
    # ref.fa so contig names match the VCFs/GFF by construction.
    REPEAT_BED="$(get_repeat_bed "${POP}")"
    if [[ -f "${REPEAT_BED}" ]]; then
        echo "  Repeat BED already present: ${REPEAT_BED}"
    fi
    echo "  Building self-alignment repeat BED (minimap2 ${MINIMAP2_SELFALIGN_OPTS})..."
    minimap2 ${MINIMAP2_SELFALIGN_OPTS} "${REF_DIR}/ref.fa" "${REF_DIR}/ref.fa" 2>/dev/null \
        | cut -f 1,3,4 \
        | sort -k1,1 -k2,2n \
        | bedtools slop -i stdin -g "${REF_DIR}/ref.fa.fai" -r 1 -l 0 \
        | bedtools merge -i stdin \
        > "${REPEAT_BED}"
    REPEAT_N=$(wc -l < "${REPEAT_BED}")
    REPEAT_BP=$(awk '{sum += $3 - $2} END {print sum + 0}' "${REPEAT_BED}")
    REF_BP=$(awk '{sum += $2} END {print sum + 0}' "${REF_DIR}/ref.fa.fai")
    REPEAT_PCT=$(awk -v m="${REPEAT_BP}" -v r="${REF_BP}" 'BEGIN {printf "%.2f", (r ? 100 * m / r : 0)}')
    echo "  Repeat BED: ${REPEAT_BED} (${REPEAT_N} intervals, ${REPEAT_BP}/${REF_BP} bp = ${REPEAT_PCT}%)"

    # Symlink snpEff config and database
    if [[ -f "${SNPEFF_CONFIG}" ]]; then
        ln -sf "${SNPEFF_CONFIG}" "${REF_DIR}/snpeff.config"
        echo "  Linked snpEff config"
    else
        echo "  WARNING: snpEff config not found at ${SNPEFF_CONFIG}" >&2
        echo "  Run 01_build_snpeff_db.sh first" >&2
    fi

    echo "  DONE: ${POP} -> ${REF_NAME}"
    echo ""
done
