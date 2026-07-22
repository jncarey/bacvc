#!/usr/bin/env bash
#
# run_population.sh
#
# Main entry point for the variant calling pipeline.
# Runs all stages for one or more populations:
#   01: Build snpEff database
#   02: Prepare reference (bwa index, samtools faidx, snpEff link, repeat BED)
#   03: Per-isolate variant calling (SGE array: fastp -> bwa -> freebayes -> snpEff)
#   04: Build the canonical merged VCF + derive outputs (held until 03 completes):
#         04a: canonical merged VCF + TreeTime VCF  (04_run_generate_variant_vcf.sh)
#         04b: variant tables, held on 04a          (04_run_generate_variant_tables.sh)
#   05: Generate per-population SNP alignment (held until 04b completes)
#
# Both step-04 artifacts derive from the single canonical merged VCF, so the
# TreeTime VCF and the tables stay consistent, and the self-alignment repeat
# exclusion (BED from step 02) is applied once when building it.
#
# Steps 01 and 02 are skipped if outputs already exist.
#
# Usage:
#   bash run_population.sh [--dry-run] [--from-step N] [POP ...]
#   # No args: runs every population in POP_TO_REF (config.sh)
#   # --from-step N (1-5): skip stages before N. Submitted stages have their
#   #   hold dependency dropped when the prior stage is skipped.
#   # e.g.: bash run_population.sh IA01 IA02 UW03
#   #       bash run_population.sh --dry-run IA01
#   #       bash run_population.sh --from-step 4 IA01   # rebuild VCF + tables + alignment
#   #       bash run_population.sh --from-step 5 IA01   # just rebuild the SNP alignment
#
# Every step 03-05 job is submitted with `-v POP=...,VC_CONFIG=<abs path to
# this repo's bin/config.sh>`. VC_CONFIG is required, not cosmetic: SGE copies
# the submitted script into a spool directory before executing it
# (/var/spool/.../job_scripts/<job-id>), so inside the job `$0`/`dirname "$0"`
# point at that spool copy, not this repo -- a self-locating `source
# "$(dirname "$0")/config.sh"` would silently source nothing. Passing the
# resolved path explicitly is what makes the job portable across clones.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
source "${SCRIPT_DIR}/config.sh"

DRY_RUN=false
FROM_STEP=1
POPS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true; shift ;;
        --from-step) FROM_STEP="$2"; shift 2 ;;
        --from-step=*) FROM_STEP="${1#*=}"; shift ;;
        *) POPS+=("$1"); shift ;;
    esac
done

if ! [[ "${FROM_STEP}" =~ ^[1-5]$ ]]; then
    echo "ERROR: --from-step must be 1, 2, 3, 4, or 5 (got '${FROM_STEP}')" >&2
    exit 1
fi

if [[ "${FROM_STEP}" -gt 1 ]]; then
    echo "Starting from step ${FROM_STEP} (skipping earlier stages)"
fi

if [[ ${#POPS[@]} -eq 0 ]]; then
    POPS=("${!POP_TO_REF[@]}")
    echo "No populations specified, running all: ${POPS[*]}"
fi

# ===================================================================
# Step 01: Build snpEff databases (shared across all populations)
# ===================================================================
if [[ "${FROM_STEP}" -le 1 ]]; then
    # Check if all needed databases exist
    SNPEFF_NEEDED=()
    for POP in "${POPS[@]}"; do
        REF_NAME="$(get_ref_name "${POP}")"
        if [[ ! -f "${SNPEFF_DIR}/${REF_NAME}/snpEffectPredictor.bin" ]]; then
            SNPEFF_NEEDED+=("${REF_NAME}")
        fi
    done

    if [[ ${#SNPEFF_NEEDED[@]} -gt 0 ]]; then
        echo "============================================"
        echo "Step 01: Building snpEff databases"
        echo "  Missing: ${SNPEFF_NEEDED[*]}"
        echo "============================================"
        if [[ "${DRY_RUN}" == "true" ]]; then
            echo "  [DRY RUN] Would run: bash ${SCRIPT_DIR}/01_build_snpeff_db.sh ${SNPEFF_NEEDED[*]}"
        else
            conda run -n "${CONDA_ENV}" \
                bash "${SCRIPT_DIR}/01_build_snpeff_db.sh" "${SNPEFF_NEEDED[@]}"
        fi
        echo ""
    else
        echo "snpEff databases already built for all requested populations."
    fi
fi

# ===================================================================
# Per-population stages
# ===================================================================
for POP in "${POPS[@]}"; do
    echo "============================================"
    echo "Population: ${POP}"
    echo "============================================"

    REF_NAME="$(get_ref_name "${POP}")"
    POP_DIR="$(get_pop_dir "${POP}")"
    REF_DIR="${POP_DIR}/reference"
    LOG_DIR="${POP_DIR}/logs"

    mkdir -p "${POP_DIR}" "${LOG_DIR}"

    # --- Prerequisite check when resuming from a later stage ---
    # Verifies that outputs from the skipped upstream stage are present, so we
    # don't submit a doomed job (e.g. step 04 against an empty variants/ dir).
    MISSING=()
    case "${FROM_STEP}" in
        2)
            [[ -f "${SNPEFF_DIR}/${REF_NAME}/snpEffectPredictor.bin" ]] \
                || MISSING+=("${SNPEFF_DIR}/${REF_NAME}/snpEffectPredictor.bin (run step 01)")
            ;;
        3)
            [[ -f "${REF_DIR}/ref.fa.bwt" ]] \
                || MISSING+=("${REF_DIR}/ref.fa.bwt (run step 02)")
            ;;
        4)
            [[ -f "${POP_DIR}/input.tab" ]] \
                || MISSING+=("${POP_DIR}/input.tab (run step 03)")
            DONE_COUNT=$(find "${POP_DIR}/variants" -name .done 2>/dev/null | wc -l)
            [[ "${DONE_COUNT}" -gt 0 ]] \
                || MISSING+=("any ${POP_DIR}/variants/*/.done marker (run step 03)")
            if [[ "${EXCLUDE_REPEATS}" == "true" ]]; then
                [[ -f "$(get_repeat_bed "${POP}")" ]] \
                    || MISSING+=("$(get_repeat_bed "${POP}") (run step 02)")
            fi
            ;;
        5)
            [[ -f "${POP_DIR}/pop_tables/snp.gt.tab" ]] \
                || MISSING+=("${POP_DIR}/pop_tables/snp.gt.tab (run step 04)")
            ;;
    esac
    if [[ ${#MISSING[@]} -gt 0 ]]; then
        echo "ERROR: --from-step ${FROM_STEP} requires outputs from earlier stages." >&2
        echo "  Missing for ${POP}:" >&2
        printf '    - %s\n' "${MISSING[@]}" >&2
        exit 1
    fi

    # --- Step 02: Prepare reference if not already done ---
    # Re-run when any step-02 output is missing -- including ref.repeats.bed,
    # which step 04 requires when EXCLUDE_REPEATS=true. 02_prepare_reference.sh is
    # idempotent (it skips the bwa index / faidx when those already exist), so a
    # re-run only fills in what's absent (e.g. the repeat BED on an old reference).
    if [[ "${FROM_STEP}" -le 2 ]]; then
        NEED_PREP=false
        [[ -f "${REF_DIR}/ref.fa.bwt" ]] || NEED_PREP=true
        if [[ "${EXCLUDE_REPEATS}" == "true" && ! -f "$(get_repeat_bed "${POP}")" ]]; then
            NEED_PREP=true
        fi
        if [[ "${NEED_PREP}" == "true" ]]; then
            echo "Step 02: Preparing reference..."
            if [[ "${DRY_RUN}" == "true" ]]; then
                echo "  [DRY RUN] Would run: bash ${SCRIPT_DIR}/02_prepare_reference.sh ${POP}"
            else
                conda run -n "${CONDA_ENV}" \
                    bash "${SCRIPT_DIR}/02_prepare_reference.sh" "${POP}"
            fi
        else
            echo "Step 02: Reference already prepared."
        fi
    fi

    # Track hold dependencies; populated only when an upstream stage is submitted
    # in this iteration.
    HOLD_FOR_04=""
    HOLD_FOR_05=""

    # --- Step 03: Discover FASTQ files and submit worker array job ---
    if [[ "${FROM_STEP}" -le 3 ]]; then
        echo ""
        echo "Discovering FASTQ files..."
        generate_input_tab "${POP}" >&2
        INPUT_TAB="${POP_DIR}/input.tab"
        N_SAMPLES=$(wc -l < "${INPUT_TAB}" 2>/dev/null || echo 0)

        if [[ "${N_SAMPLES}" -eq 0 ]]; then
            echo "WARNING: No samples found for ${POP}. Skipping." >&2
            continue
        fi

        echo "  ${N_SAMPLES} samples to process"

        echo ""
        echo "Step 03: Submitting SGE array job (${N_SAMPLES} tasks, max ${MAX_CONCURRENT} concurrent)..."

        WORKER_CMD="qsub \
            -N vc_${POP} \
            -t 1-${N_SAMPLES} \
            -tc ${MAX_CONCURRENT} \
            -l mfree=${SGE_MEMORY} \
            -l h_rt=${SGE_WALLTIME} \
            -pe serial ${BWA_THREADS} \
            -o ${LOG_DIR} \
            -e ${LOG_DIR} \
            -v POP=${POP},VC_CONFIG=${SCRIPT_DIR}/config.sh \
            ${SCRIPT_DIR}/03_worker.sh"

        if [[ "${DRY_RUN}" == "true" ]]; then
            echo "  [DRY RUN] ${WORKER_CMD}"
        else
            WORKER_JOB_ID=$(${WORKER_CMD} | grep -oP 'Your job(-array)? \K[0-9]+' | head -1)
            echo "  Submitted worker array job: ${WORKER_JOB_ID}"
        fi
        HOLD_FOR_04="-hold_jid vc_${POP}"
    fi

    # --- Step 04: Build canonical merged VCF, then derive TreeTime VCF + tables ---
    # 04a (tv_) builds the canonical merged VCF and the TreeTime VCF, held on the
    # worker. 04b (vt_) builds the variant tables from that canonical VCF, held on
    # 04a so the two artifacts stay consistent.
    if [[ "${FROM_STEP}" -le 4 ]]; then
        echo ""
        echo "Step 04a: Submitting canonical/TreeTime VCF job${HOLD_FOR_04:+ (held on vc_${POP})}..."

        VCF_CMD="qsub \
            -N tv_${POP} \
            ${HOLD_FOR_04} \
            -l mfree=${SGE_MEMORY} \
            -l h_rt=2:00:00 \
            -j y \
            -o ${LOG_DIR} \
            -v POP=${POP},VC_CONFIG=${SCRIPT_DIR}/config.sh \
            ${SCRIPT_DIR}/04_run_generate_variant_vcf.sh"

        if [[ "${DRY_RUN}" == "true" ]]; then
            echo "  [DRY RUN] ${VCF_CMD}"
        else
            VCF_JOB_ID=$(${VCF_CMD} | grep -oP 'Your job(-array)? \K[0-9]+' | head -1)
            echo "  Submitted canonical/TreeTime VCF job: ${VCF_JOB_ID}"
        fi

        echo "Step 04b: Submitting variant table job (held on tv_${POP})..."

        TABLE_CMD="qsub \
            -N vt_${POP} \
            -hold_jid tv_${POP} \
            -l mfree=${SGE_MEMORY} \
            -l h_rt=4:00:00 \
            -j y \
            -o ${LOG_DIR} \
            -v POP=${POP},VC_CONFIG=${SCRIPT_DIR}/config.sh \
            ${SCRIPT_DIR}/04_run_generate_variant_tables.sh"

        if [[ "${DRY_RUN}" == "true" ]]; then
            echo "  [DRY RUN] ${TABLE_CMD}"
        else
            TABLE_JOB_ID=$(${TABLE_CMD} | grep -oP 'Your job(-array)? \K[0-9]+' | head -1)
            echo "  Submitted variant table job: ${TABLE_JOB_ID}"
        fi
        HOLD_FOR_05="-hold_jid vt_${POP}"
    fi

    # --- Step 05: Submit SNP alignment generation ---
    if [[ "${FROM_STEP}" -le 5 ]]; then
        echo ""
        echo "Step 05: Submitting SNP alignment job${HOLD_FOR_05:+ (held on vt_${POP})}..."

        ALN_CMD="qsub \
            -N aln_${POP} \
            ${HOLD_FOR_05} \
            -l mfree=${SGE_MEMORY} \
            -l h_rt=2:00:00 \
            -j y \
            -o ${LOG_DIR} \
            -v POP=${POP},VC_CONFIG=${SCRIPT_DIR}/config.sh \
            ${SCRIPT_DIR}/05_run_generate_snp_alignment.sh"

        if [[ "${DRY_RUN}" == "true" ]]; then
            echo "  [DRY RUN] ${ALN_CMD}"
        else
            ALN_JOB_ID=$(${ALN_CMD} | grep -oP 'Your job(-array)? \K[0-9]+' | head -1)
            echo "  Submitted SNP alignment job: ${ALN_JOB_ID}"
        fi
    fi

    echo ""
done

echo "============================================"
echo "All populations submitted."
if [[ "${DRY_RUN}" == "true" ]]; then
    echo "(DRY RUN - no jobs were submitted)"
fi
echo "Monitor with: qstat -u \$USER"
echo "============================================"
