#!/usr/bin/env bash
#
# 03_worker.sh
#
# Per-isolate variant calling worker (SGE array job task).
# Performs the full pipeline for a single sample:
#   1. fastp QC
#   2. bwa mem alignment with samclip + duplicate removal
#   3. freebayes variant calling
#   4. vt normalize + snpEff annotation
#   5. mosdepth per-base depth (for low-coverage masking in step 04)
#   6. samtools flagstat (mapping rate, for the QC gate below)
#   7. per-isolate QC flag (low_coverage / mapping_failure -- see
#      pipeline_helpers.qc_flag_isolate; does not stop the pipeline, just
#      records the verdict for the caller to act on)
#
# Usage (SGE array):
#   Called by run_population.sh; not meant to be run directly.
#   Reads task info from input.tab using SGE_TASK_ID.
#
# Environment variables (set by run_population.sh via qsub -v):
#   POP          - population name (e.g. IA01)
#   VC_CONFIG    - absolute path to this repo's bin/config.sh. SGE spools the
#                  submitted script before executing it, so `dirname "$0"`
#                  inside the job does NOT point at this repo -- VC_CONFIG is
#                  how the job finds config.sh regardless of clone location.
#   SGE_TASK_ID  - array task index (1-based, set by SGE)

#$ -S /bin/bash
#$ -cwd
#$ -V

set -euo pipefail

if [[ -n "${VC_CONFIG:-}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${VC_CONFIG}")" && pwd)"
else
    # Fallback for direct (non-qsub) invocation, e.g. `bash 03_worker.sh` for
    # local debugging, where $0 still resolves correctly.
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
fi
source "${SCRIPT_DIR}/config.sh"

# Activate the pipeline env for bwa, samtools, samclip, freebayes, bcftools, vt, snpEff
eval "$(conda shell.bash hook 2>/dev/null)"
conda activate "${CONDA_ENV}"

# ---------------------------------------------------------------------------
# Resolve sample from input.tab
# ---------------------------------------------------------------------------
POP_DIR="$(get_pop_dir "${POP}")"
INPUT_TAB="${POP_DIR}/input.tab"

if [[ ! -f "${INPUT_TAB}" ]]; then
    echo "ERROR: input.tab not found: ${INPUT_TAB}" >&2
    exit 1
fi

LINE=$(sed -n "${SGE_TASK_ID}p" "${INPUT_TAB}")
if [[ -z "${LINE}" ]]; then
    echo "ERROR: No entry for SGE_TASK_ID=${SGE_TASK_ID} in ${INPUT_TAB}" >&2
    exit 1
fi

SAMPLE_ID=$(echo "${LINE}" | cut -f1)
R1=$(echo "${LINE}" | cut -f2)
R2=$(echo "${LINE}" | cut -f3)

echo "=== Worker: ${SAMPLE_ID} (task ${SGE_TASK_ID}) ==="
echo "  Population: ${POP}"
echo "  R1: ${R1}"
echo "  R2: ${R2}"

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REF_NAME="$(get_ref_name "${POP}")"
REF_DIR="${POP_DIR}/reference"
REF_FA="${REF_DIR}/ref.fa"

TRIM_DIR="${POP_DIR}/trimmed_reads/${SAMPLE_ID}"
VAR_DIR="${POP_DIR}/variants/${SAMPLE_ID}"

mkdir -p "${TRIM_DIR}" "${VAR_DIR}"

R1_TRIMMED="${TRIM_DIR}/${SAMPLE_ID}_R1.trimmed.fastq.gz"
R2_TRIMMED="${TRIM_DIR}/${SAMPLE_ID}_R2.trimmed.fastq.gz"
BAM="${VAR_DIR}/${SAMPLE_ID}.bam"

# Skip if already done
if [[ -f "${VAR_DIR}/.done" ]]; then
    echo "  Already completed (${VAR_DIR}/.done exists). Skipping."
    exit 0
fi

# Validate reference
if [[ ! -f "${REF_FA}" ]]; then
    echo "ERROR: Reference not prepared. Run 02_prepare_reference.sh ${POP} first." >&2
    exit 1
fi

# ===================================================================
# Step 1: fastp QC
# ===================================================================
echo "  [1/7] fastp QC..."

conda run -n "${CONDA_ENV}" \
    fastp \
        --in1 "${R1}" \
        --in2 "${R2}" \
        --out1 "${R1_TRIMMED}" \
        --out2 "${R2_TRIMMED}" \
        --length_required "${FASTP_LENGTH_REQUIRED}" \
        --json "${TRIM_DIR}/${SAMPLE_ID}.fastp.json" \
        --html "${TRIM_DIR}/${SAMPLE_ID}.fastp.html" \
        --thread "${BWA_THREADS}" \
        ${FASTP_EXTRA_OPTS}

echo "  [1/7] fastp done."

# ===================================================================
# Step 2: bwa mem -> samclip -> sort -> fixmate -> sort -> markdup
# ===================================================================
echo "  [2/7] Alignment (bwa mem + samclip + markdup)..."

RG="@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}"

bwa mem \
    -t "${BWA_THREADS}" \
    -R "${RG}" \
    "${REF_FA}" \
    "${R1_TRIMMED}" "${R2_TRIMMED}" \
    | samclip --max "${SAMCLIP_MAX_SOFT}" --ref "${REF_FA}.fai" \
    | samtools sort -n -@ "${BWA_THREADS}" - \
    | samtools fixmate -m - - \
    | samtools sort -@ "${BWA_THREADS}" - \
    | samtools markdup -r -s - - \
    > "${BAM}"

samtools index "${BAM}"

echo "  [2/7] Alignment done. BAM: ${BAM}"

# Clean up trimmed reads
if [[ "${REMOVE_TRIMMED}" == "true" ]]; then
    rm -f "${R1_TRIMMED}" "${R2_TRIMMED}"
    rm -f "${TRIM_DIR}/${SAMPLE_ID}.fastp.json" "${TRIM_DIR}/${SAMPLE_ID}.fastp.html"
    rmdir "${TRIM_DIR}" 2>/dev/null || true
    echo "  Cleaned up trimmed reads."
fi

# ===================================================================
# Step 3: freebayes variant calling
# ===================================================================
echo "  [3/7] Freebayes variant calling..."

freebayes \
    -p "${FB_PLOIDY}" \
    -P 0 \
    -C 2 \
    -F "${FB_MIN_ALT_FRAC}" \
    --min-coverage "${FB_MIN_COV}" \
    --min-repeat-entropy "${FB_MIN_REPEAT_ENTROPY}" \
    -q "${FB_BASE_QUAL}" \
    -m "${FB_MAP_QUAL}" \
    --strict-vcf \
    -f "${REF_FA}" \
    "${BAM}" \
    > "${VAR_DIR}/snps.raw.vcf"

echo "  [3/7] Freebayes done. Raw variants: ${VAR_DIR}/snps.raw.vcf"

# ===================================================================
# Step 4: Normalize + snpEff annotation
# ===================================================================
echo "  [4/7] Normalizing and annotating with snpEff..."
# No DP filter here: low-coverage variants are retained so the DP decision is
# made once, at stage 04 (+setGT FMT/DP<VT_MIN_DP mask + mosdepth MASK_MIN_DP
# mask). The normalize+annotate pipe lives in normalize_annotate_vcf()
# (config.sh) so this stage's pipe is defined in exactly one place.
normalize_annotate_vcf \
    "${VAR_DIR}/snps.raw.vcf" "${REF_FA}" "${REF_NAME}" \
    "${VAR_DIR}/snps.norm.annot.vcf"

echo "  [4/7] Annotated VCF: ${VAR_DIR}/snps.norm.annot.vcf"

# ===================================================================
# Step 5: Per-base depth (consumed by step 04 for low-coverage masking, and
# by step 7's QC gate for mean depth)
# ===================================================================
echo "  [5/7] mosdepth per-base depth..."

# mosdepth ships in the pipeline conda env (environment.yml), which is already
# active here via the snippy conda env.
mosdepth \
    --fast-mode \
    --threads "${BWA_THREADS}" \
    "${VAR_DIR}/${SAMPLE_ID}" \
    "${BAM}"

# Drop the per-base distribution histogram; keep the per-base bedGraph +
# index (step 04's coverage mask) and summary.txt (step 7's mean depth).
rm -f "${VAR_DIR}/${SAMPLE_ID}.mosdepth.global.dist.txt"

echo "  [5/7] Depth: ${VAR_DIR}/${SAMPLE_ID}.per-base.bed.gz"

# ===================================================================
# Step 6: samtools flagstat (mapping rate, consumed by step 7's QC gate).
# Must run before the BAM cleanup below.
# ===================================================================
echo "  [6/7] samtools flagstat..."

samtools flagstat "${BAM}" > "${VAR_DIR}/${SAMPLE_ID}.flagstat.txt"

echo "  [6/7] Flagstat: ${VAR_DIR}/${SAMPLE_ID}.flagstat.txt"

# ===================================================================
# Step 7: Per-isolate QC flag (low_coverage / mapping_failure)
# ===================================================================
echo "  [7/7] QC flag..."

python3 "${SCRIPT_DIR}/qc_flag_isolate.py" \
    --isolate "${SAMPLE_ID}" \
    --fastp-json "${TRIM_DIR}/${SAMPLE_ID}.fastp.json" \
    --flagstat "${VAR_DIR}/${SAMPLE_ID}.flagstat.txt" \
    --mosdepth-summary "${VAR_DIR}/${SAMPLE_ID}.mosdepth.summary.txt" \
    --min-mean-depth "${MIN_MEAN_DEPTH}" \
    --min-mapped-pct "${MIN_MAPPED_PCT}" \
    --out "${VAR_DIR}/${SAMPLE_ID}.qc_flag.tsv"

echo "  [7/7] QC flag: ${VAR_DIR}/${SAMPLE_ID}.qc_flag.tsv"

# ===================================================================
# Cleanup
# ===================================================================
if [[ "${REMOVE_BAMS}" == "true" ]]; then
    rm -f "${BAM}" "${BAM}.bai"
    echo "  Cleaned up BAM files."
fi

# Mark completion
touch "${VAR_DIR}/.done"
echo "  DONE: ${SAMPLE_ID}"
