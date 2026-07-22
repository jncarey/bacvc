#!/usr/bin/env bash
#
# 04_run_generate_variant_tables.sh
#
# Generate population variant tables from the canonical merged VCF built by
# 04_run_generate_variant_vcf.sh (${POP}.merged.vcf.gz). Both artifacts derive
# from that single VCF, so the tables and the TreeTime VCF agree on every shared
# (CHR,POS,REF,ALT) key and the self-alignment repeat exclusion carries through.
#
# Pipeline per population:
#   1. bcftools norm -m-   split the canonical (multi-allelic) VCF into biallelic
#      rows -> one (CHR,POS,REF,ALT) key per row, matching the table model.
#   2. snpEff ann          re-annotate the biallelic rows (same DB/config as
#      03_worker.sh) so GT, INFO/TYPE and INFO/ANN all come from the same records
#      -- annotation keys agree with genotype keys by construction.
#   3. 04_generate_variant_tables.py  emit var/snp/indel/disruptive tables.
#
# Can be submitted as an SGE job (held until the VCF builder completes) or run
# directly.
#
# Usage:
#   bash 04_run_generate_variant_tables.sh              # all populations
#   bash 04_run_generate_variant_tables.sh <POP> [POP]  # specific population(s)
#   # SGE submission (single pop) -- VC_CONFIG must be passed explicitly:
#   # SGE spools the submitted script before executing it, so `dirname "$0"`
#   # inside the job would not point at this repo.
#   qsub -v POP=IA01,VC_CONFIG=/abs/path/to/vc/bin/config.sh \
#       04_run_generate_variant_tables.sh

#$ -S /bin/bash
#$ -cwd
#$ -V

set -euo pipefail

if [[ -n "${VC_CONFIG:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${VC_CONFIG}")" && pwd)"
else
    # Fallback for direct (non-qsub) invocation, where $0 still resolves.
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
fi
source "${SCRIPT_DIR}/config.sh"

# Accept POP(s) as arguments, environment variable, or default to all 9
if [[ $# -gt 0 ]]; then
    POPS=("$@")
elif [[ -n "${POP:-}" ]]; then
    POPS=("${POP}")
else
    POPS=("${!POP_TO_REF[@]}")
fi

for POP in "${POPS[@]}"; do
    REF_NAME="$(get_ref_name "${POP}")"
    POP_DIR="$(get_pop_dir "${POP}")"
    REF_GFF="${POP_DIR}/reference/ref.gff"
    MERGED_GZ="${POP_DIR}/aln/${POP}.merged.vcf.gz"
    OUT_DIR="${POP_DIR}/pop_tables"

    echo "=== Generating variant tables for ${POP} ==="
    echo "  Reference: ${REF_NAME}"
    echo "  Canonical: ${MERGED_GZ}"
    echo "  Output:    ${OUT_DIR}"

    if [[ ! -f "${MERGED_GZ}" ]]; then
        echo "WARNING: canonical VCF not found: ${MERGED_GZ}" >&2
        echo "  Run 04_run_generate_variant_vcf.sh ${POP} first. Skipping." >&2
        continue
    fi

    # This step re-runs snpEff, so the genome must be in the snpEff config.
    # Fail loud (rather than a cryptic Java stack trace) if it is missing.
    if ! grep -q "^${REF_NAME}\.genome" "${SNPEFF_CONFIG}" 2>/dev/null; then
        echo "ERROR: '${REF_NAME}.genome' not found in ${SNPEFF_CONFIG}" >&2
        echo "       Run 01_build_snpeff_db.sh first. Skipping ${POP}." >&2
        continue
    fi

    WORK=$(mktemp -d -t "vt-${POP}-XXXXXX")
    # shellcheck disable=SC2064
    trap "rm -rf '${WORK}'" EXIT
    ANNOT_VCF="${WORK}/${POP}.biallelic.annot.vcf"

    conda run -n "${CONDA_ENV}" bash -c "
        set -euo pipefail

        # --- Step 1+2: split multi-allelics, re-annotate with snpEff ---
        echo '  [1/2] bcftools norm -m- | snpEff ann...'
        bcftools norm -m- '${MERGED_GZ}' 2>/dev/null \
          | snpEff ann \
                -noLog -noStats \
                -no-downstream -no-upstream -no-utr \
                -c '${SNPEFF_CONFIG}' \
                -dataDir '${SNPEFF_DIR}' \
                '${REF_NAME}' \
          > '${ANNOT_VCF}'

        # --- Step 3: build tables ---
        echo '  [2/2] generating tables...'
        python3 '${SCRIPT_DIR}/04_generate_variant_tables.py' \
            --pop '${POP}' \
            --merged-vcf '${ANNOT_VCF}' \
            --gff '${REF_GFF}' \
            --out-dir '${OUT_DIR}'
    "

    rm -rf "${WORK}"
    trap - EXIT

    echo "=== Done: ${POP} ==="
    echo ""
done
