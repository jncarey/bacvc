#!/usr/bin/env bash
#
# 04_run_generate_variant_vcf.sh
#
# Build the canonical multi-sample VCF per population (the single source of
# truth for both the TreeTime VCF and the variant tables in
# 04_run_generate_variant_tables.sh), directly from per-isolate
# snps.norm.annot.vcf files.
#
# Pipeline per population:
#   1. Per-sample setGT on snps.norm.annot.vcf (ALL variant types kept):
#        - bcftools +setGT cells to ./. on:
#              QUAL < VT_MIN_QUAL
#            | FMT/DP < VT_MIN_DP
#            | het GT
#            | hom-alt with FMT/AO/FMT/DP < VT_MIN_AO_FRAC
#            | hom-ref with FMT/RO/FMT/DP < VT_MIN_RO_FRAC
#        - strip INFO down to INFO/TYPE, keep only FORMAT/GT.
#   2. bcftools merge --missing-to-ref  -> samples absent at a site default
#      to 0/0; samples that had a row keep whatever GT setGT produced
#      (1/1, 2/2, ..., or ./.). Multi-allelic rows carry per-allele TYPE.
#   3. 04_generate_variant_vcf.py: drop records overlapping the self-alignment
#      repeat BED (ref.repeats.bed), then mosdepth-mask -- for each 0/0 cell,
#      demote to ./. when per-base coverage < MASK_MIN_DP. Never touches cells
#      that aren't 0/0 (./. and ALT calls preserved).
#   4. bgzip + tabix -> canonical merged VCF (all types, INFO/TYPE retained).
#   5. Derive the TreeTime VCF: split multi-allelics, keep INFO/TYPE in
#      {snp,mnp}, rejoin multi-allelics (preserving 2/2 etc.), strip INFO.
#
# Filtering matches config.sh: VT_MIN_DP, VT_MIN_QUAL, VT_MIN_AO_FRAC,
# VT_MIN_RO_FRAC, MASK_MIN_DP, EXCLUDE_REPEATS.
#
# Outputs: ${POP_DIR}/aln/${POP}.merged.vcf.gz   (canonical, all types) + .tbi
#          ${POP_DIR}/aln/${POP}.treetime.vcf.gz (snp/mnp subset)        + .tbi
#
# Usage:
#   bash 04_run_generate_variant_vcf.sh              # all populations
#   bash 04_run_generate_variant_vcf.sh <POP> [POP]  # specific population(s)
#   # SGE submission (single pop) -- VC_CONFIG must be passed explicitly:
#   # SGE spools the submitted script before executing it, so `dirname "$0"`
#   # inside the job would not point at this repo.
#   qsub -v POP=IA01,VC_CONFIG=/abs/path/to/vc/bin/config.sh \
#       04_run_generate_variant_vcf.sh

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

if [[ $# -gt 0 ]]; then
    POPS=("$@")
elif [[ -n "${POP:-}" ]]; then
    POPS=("${POP}")
else
    POPS=("${!POP_TO_REF[@]}")
fi

# Genotype-level mask expression for +setGT.
SETGT_EXPR="QUAL<${VT_MIN_QUAL} | FMT/DP<${VT_MIN_DP} | GT=\"het\""
SETGT_EXPR+=" | (GT=\"AA\" & FMT/DP>0 & FMT/AO/FMT/DP<${VT_MIN_AO_FRAC})"
SETGT_EXPR+=" | (GT=\"RR\" & FMT/DP>0 & FMT/RO/FMT/DP<${VT_MIN_RO_FRAC})"

for POP in "${POPS[@]}"; do
    POP_DIR="$(get_pop_dir "${POP}")"
    VAR_DIR="${POP_DIR}/variants"
    ALN_DIR="${POP_DIR}/aln"
    MERGED_GZ="${ALN_DIR}/${POP}.merged.vcf.gz"
    OUT_GZ="${ALN_DIR}/${POP}.treetime.vcf.gz"
    REPEAT_BED="$(get_repeat_bed "${POP}")"

    echo "=== Building canonical merged VCF + TreeTime VCF for ${POP} ==="
    echo "  Variants: ${VAR_DIR}"
    echo "  Canonical: ${MERGED_GZ}"
    echo "  TreeTime:  ${OUT_GZ}"

    if [[ ! -d "${VAR_DIR}" ]]; then
        echo "WARNING: ${VAR_DIR} not found. Skipping ${POP}." >&2
        continue
    fi

    # Resolve repeat-exclusion: fail loud if enabled but the BED is missing,
    # since silently shipping unfiltered calls defeats the purpose.
    REPEAT_ARG=""
    if [[ "${EXCLUDE_REPEATS}" == "true" ]]; then
        if [[ ! -f "${REPEAT_BED}" ]]; then
            echo "ERROR: EXCLUDE_REPEATS=true but repeat BED not found: ${REPEAT_BED}" >&2
            echo "       Run 02_prepare_reference.sh ${POP} first." >&2
            exit 1
        fi
        REPEAT_ARG="--repeat-bed ${REPEAT_BED}"
        echo "  Repeats:   ${REPEAT_BED} (excluded)"
    else
        echo "  Repeats:   disabled (EXCLUDE_REPEATS=${EXCLUDE_REPEATS})"
    fi

    # Population SNP-site missingness filter. The per-site '-' distribution is
    # always written; the actual drop is gated on EXCLUDE_HIGH_MISSING_SNP_SITES.
    MISSING_TSV="${ALN_DIR}/${POP}.snp_site_missingness.tsv"
    MISSING_ARG=""
    if [[ "${EXCLUDE_HIGH_MISSING_SNP_SITES}" == "true" ]]; then
        MISSING_ARG="--max-missing-frac ${SNP_SITE_MAX_MISSING_FRAC}"
        echo "  Missing:   drop snp/mnp sites with frac '-' >= ${SNP_SITE_MAX_MISSING_FRAC}"
    else
        echo "  Missing:   filter disabled (EXCLUDE_HIGH_MISSING_SNP_SITES=${EXCLUDE_HIGH_MISSING_SNP_SITES}); QC only"
    fi

    # Discover completed isolates (.done marker + snps.norm.annot.vcf present)
    ISOLATES=()
    while IFS= read -r d; do
        iso=$(basename "$d")
        if [[ -f "${VAR_DIR}/${iso}/snps.norm.annot.vcf" && -f "${VAR_DIR}/${iso}/.done" ]]; then
            ISOLATES+=("${iso}")
        fi
    done < <(find "${VAR_DIR}" -mindepth 1 -maxdepth 1 -type d | sort)

    if [[ ${#ISOLATES[@]} -eq 0 ]]; then
        echo "WARNING: no completed isolates in ${VAR_DIR}. Skipping ${POP}." >&2
        continue
    fi
    echo "  Isolates: ${#ISOLATES[@]}"

    mkdir -p "${ALN_DIR}"

    WORK=$(mktemp -d -t "treetime-${POP}-XXXXXX")
    # shellcheck disable=SC2064
    trap "rm -rf '${WORK}'" EXIT

    conda run -n "${CONDA_ENV}" bash -c "
        set -euo pipefail

        # --- Step 1: per-sample setGT (all variant types; keep INFO/TYPE) ---
        echo '  [1/5] Per-sample setGT (${#ISOLATES[@]} samples)...'
        for ISO in ${ISOLATES[*]}; do
            bcftools +setGT \"${VAR_DIR}/\${ISO}/snps.norm.annot.vcf\" -Ou \
                    -- -t q -n . -i '${SETGT_EXPR}' \
              | bcftools annotate -x '^INFO/TYPE,^FORMAT/GT' -Oz \
                    -o \"${WORK}/\${ISO}.vcf.gz\"
            bcftools index -t \"${WORK}/\${ISO}.vcf.gz\"
        done

        # --- Step 2: merge (absent samples default to 0/0) ---
        echo '  [2/5] bcftools merge --missing-to-ref...'
        VCF_LIST=''
        for ISO in ${ISOLATES[*]}; do
            VCF_LIST=\"\${VCF_LIST} ${WORK}/\${ISO}.vcf.gz\"
        done
        # shellcheck disable=SC2086
        bcftools merge --missing-to-ref \
            -Oz -o \"${WORK}/merged.raw.vcf.gz\" \
            \${VCF_LIST}
        bcftools index -t \"${WORK}/merged.raw.vcf.gz\"

        # --- Step 3: repeat exclusion + mosdepth mask + missingness filter ---
        echo '  [3/5] repeat exclusion + mosdepth mask + missingness filter...'
        python3 '${SCRIPT_DIR}/04_generate_variant_vcf.py' \
            --merged \"${WORK}/merged.raw.vcf.gz\" \
            --var-dir \"${VAR_DIR}\" \
            --mask-min-dp '${MASK_MIN_DP}' \
            --pop '${POP}' \
            --missingness-tsv \"${MISSING_TSV}\" \
            ${REPEAT_ARG} ${MISSING_ARG} \
            > \"${WORK}/merged.masked.vcf\"

        # --- Step 4: persist canonical merged VCF (all types, INFO/TYPE kept) ---
        echo '  [4/5] bgzip + tabix canonical merged VCF...'
        bgzip -c \"${WORK}/merged.masked.vcf\" > \"${MERGED_GZ}\"
        tabix -f -p vcf \"${MERGED_GZ}\"

        # --- Step 5: derive TreeTime VCF (snp/mnp only, multi-allelic preserved) ---
        # split -> keep snp/mnp alleles -> rejoin multi-allelic -> strip INFO.
        echo '  [5/5] deriving TreeTime VCF (snp/mnp subset)...'
        bcftools norm -m- \"${MERGED_GZ}\" 2>/dev/null \
          | bcftools view -i 'INFO/TYPE=\"snp\" || INFO/TYPE=\"mnp\"' \
          | bcftools norm -m+ 2>/dev/null \
          | bcftools annotate -x INFO -Oz -o \"${OUT_GZ}\"
        tabix -f -p vcf \"${OUT_GZ}\"
    "

    rm -rf "${WORK}"
    trap - EXIT

    echo "=== Done: ${POP} ==="
    echo ""
done
