#!/usr/bin/env bash
#
# config.sh - Central configuration for the variant calling pipeline
#
# Source this file from other pipeline scripts:
#   source "$(dirname "$0")/config.sh"

# === Paths ===
VC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Inputs (symlinks under vc/data/ point to the canonical source locations)
DATA_DIR="${VC_DIR}/data"
FASTQ_BASE_DIR="${DATA_DIR}/allreads"
BAKTA_GFF_DIR="${DATA_DIR}/bakta_gff"

# Outputs (everything the pipeline produces lives under vc/results/)
RESULTS_DIR="${VC_DIR}/results"
SNPEFF_DIR="${RESULTS_DIR}/snpeff"
SNPEFF_CONFIG="${SNPEFF_DIR}/snpeff.config"

# === Population -> Reference genome mapping ===
# Each population uses its own reference genome.
declare -A POP_TO_REF=(
    [IA01]="IA01_autocycler"
    [IA02]="IA02_autocycler"
    [IA03]="IA03_autocycler"
    [IA04]="IA04_autocycler"
    [UW01]="UW01_autocycler"
    [UW03]="UW03_autocycler"
    [UW04]="UW04_autocycler"
    [UW05]="UW05_autocycler"
    [UW06]="UW06_autocycler"
    [POP]="POP_autocycler"   # bundled demo dataset, see example/
)

# === Conda environment ===
# Single env built from environment.yml (repo root):
# `mamba env create -f environment.yml`. Covers every stage's tools (fastp in
# 03_worker.sh; bwa/samtools/bcftools/freebayes/vt/snpEff/samclip elsewhere).
CONDA_ENV="vc"

# === fastp QC ===
FASTP_LENGTH_REQUIRED=50
FASTP_EXTRA_OPTS=""

# === bwa mem ===
BWA_THREADS=1

# === samclip ===
SAMCLIP_MAX_SOFT=10   # max soft-clip length before clipping

# === freebayes ===
FB_PLOIDY=2
FB_MIN_COV=10
FB_MIN_ALT_FRAC=0.05
FB_MIN_REPEAT_ENTROPY=1.0
FB_BASE_QUAL=13
FB_MAP_QUAL=60

# === VCF pre-filter (03_worker.sh, before snpEff) ===
FILT_MIN_QUAL=100

# === Variant table generation ===
# VT_MIN_DP is the single DP knob for the whole pipeline. It now feeds three
# places (all at stage 04, never at stage 03): the mosdepth low-coverage mask
# (MASK_MIN_DP below), the +setGT genotype mask in 04_run_generate_variant_vcf.sh,
# and the variant-table genotype DP gate in 04_generate_variant_tables.py.
VT_MIN_DP=15
VT_MIN_QUAL=300
VT_MIN_AO_FRAC=0.8
VT_MIN_RO_FRAC=0.8

# === Per-base depth masking (step 04) ===
# Mask any sample/segregating-site cell with mosdepth depth < this threshold.
# Catches the asymmetric case where a sample has no VCF record at a site
# called by another isolate but is itself low-coverage there.
# Aligned with VT_MIN_DP so a "called 0" and a "defaulted 0" use the same bar.
MASK_MIN_DP="${VT_MIN_DP}"

# === Self-alignment repeat masking (step 02 + step 04) ===
# Each reference is self-aligned with minimap2 (step 02) to find repetitive /
# self-similar regions, emitted as ${REF_DIR}/ref.repeats.bed. Variants falling
# in those regions are excluded from the canonical merged VCF in step 04, so the
# exclusion is applied once and shared by every downstream artifact.
EXCLUDE_REPEATS=true
MINIMAP2_SELFALIGN_OPTS="-D -P -w19 -m200"

# === Population SNP-site missingness filter (step 04) ===
# Drop a snp/mnp site from the canonical VCF (and thus the TreeTime VCF and the
# SNP tables / alignment) when at least this fraction of isolates are missing
# ('-') after low-coverage masking. Applies to snp/mnp records only; indels are
# never affected. A per-site '-' distribution is always written to
# ${ALN_DIR}/${POP}.snp_site_missingness.tsv regardless of this toggle.
EXCLUDE_HIGH_MISSING_SNP_SITES=true
SNP_SITE_MAX_MISSING_FRAC=0.80   # exclude when frac_missing >= this

# === SGE cluster ===
MAX_CONCURRENT=60
SGE_MEMORY="8G"
SGE_WALLTIME="48:00:00"

# === Cleanup ===
REMOVE_TRIMMED=true
REMOVE_BAMS=true

# === Output directory structure ===
# Each population gets:
#   ${RESULTS_DIR}/populations/${POP}/
#     reference/              - bwa-indexed reference + snpEff link
#     trimmed_reads/<sample>/ - fastp output (cleaned up if REMOVE_TRIMMED)
#     variants/<sample>/      - VCFs, BAMs
#     pop_tables/             - population variant tables
#     logs/                   - per-sample log files

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# Get reference genome name for a population
get_ref_name() {
    local pop="$1"
    local ref="${POP_TO_REF[$pop]:-}"
    if [[ -z "$ref" ]]; then
        echo "ERROR: No reference mapping for population '${pop}'" >&2
        return 1
    fi
    echo "$ref"
}

# Get path to reference FASTA (.fna from bakta, contig names match GFF3)
get_ref_fasta() {
    local ref_name
    ref_name="$(get_ref_name "$1")" || return 1
    echo "${BAKTA_GFF_DIR}/${ref_name}/${ref_name}.fna"
}

# Get path to reference GFF3 (snpEff genes.gff built from .gbff)
get_ref_gff() {
    local ref_name
    ref_name="$(get_ref_name "$1")" || return 1
    echo "${SNPEFF_DIR}/${ref_name}/genes.gff"
}

# Get population output directory
get_pop_dir() {
    echo "${RESULTS_DIR}/populations/$1"
}

# Get path to the self-alignment repeat BED for a population's reference.
# Generated by 02_prepare_reference.sh, consumed by step 04.
get_repeat_bed() {
    echo "$(get_pop_dir "$1")/reference/ref.repeats.bed"
}

# Normalize + snpEff-annotate a raw freebayes VCF into a per-isolate
# snps.norm.annot.vcf. Used by stage 03 (03_worker.sh).
#
# IMPORTANT: no FMT/DP filter is applied here. Low-coverage variants are
# retained so the DP decision is made once, at stage 04 (the +setGT
# FMT/DP<VT_MIN_DP mask plus the mosdepth MASK_MIN_DP mask). Only
# QUAL>=FILT_MIN_QUAL is applied. FORMAT/DP is kept in the output so stage 04
# can still gate on it.
#
# Must run inside CONDA_ENV (bcftools, vt, snpEff). Writes atomically:
# a mid-pipe failure under `set -o pipefail` leaves the existing output intact
# rather than a truncated VCF that stage 04 would silently merge.
#
# Usage: normalize_annotate_vcf <raw_vcf> <ref_fa> <ref_name> <out_vcf>
normalize_annotate_vcf() {
    local raw_vcf="$1" ref_fa="$2" ref_name="$3" out_vcf="$4"
    local tmp="${out_vcf}.tmp.$$"
    bcftools view --include "QUAL>=${FILT_MIN_QUAL}" "${raw_vcf}" \
        | vt normalize -r "${ref_fa}" - \
        | bcftools annotate \
            --remove '^INFO/TYPE,^INFO/DP,^INFO/RO,^INFO/AO,^INFO/AB,^FORMAT/GT,^FORMAT/DP,^FORMAT/RO,^FORMAT/AO,^FORMAT/QR,^FORMAT/QA,^FORMAT/GL' \
        | snpEff ann \
            -noLog -noStats \
            -no-downstream -no-upstream -no-utr \
            -c "${SNPEFF_CONFIG}" \
            -dataDir "${SNPEFF_DIR}" \
            "${ref_name}" \
        > "${tmp}"
    mv "${tmp}" "${out_vcf}"
}

# Discover FASTQ input tab for a population
# Writes input.tab: SAMPLE_ID <tab> R1_path <tab> R2_path
generate_input_tab() {
    local pop="$1"
    local pop_dir
    pop_dir="$(get_pop_dir "$pop")"
    local input_tab="${pop_dir}/input.tab"

    mkdir -p "${pop_dir}"
    > "${input_tab}"

    local count=0
    for timepoint in "pre" "post"; do
        local fq_dir="${FASTQ_BASE_DIR}/${pop}${timepoint}"
        if [[ ! -d "${fq_dir}" ]]; then
            echo "  Note: ${fq_dir} does not exist (no ${timepoint} samples for ${pop})"
            continue
        fi
        for r1 in "${fq_dir}"/*_1.fastq.gz; do
            [[ -f "$r1" ]] || continue
            local r2="${r1/_1.fastq.gz/_2.fastq.gz}"
            if [[ ! -f "$r2" ]]; then
                echo "  WARNING: Missing R2 for $(basename "$r1")" >&2
                continue
            fi
            local sample_id
            sample_id="$(basename "$r1" _1.fastq.gz)"
            printf '%s\t%s\t%s\n' "${sample_id}" "${r1}" "${r2}" >> "${input_tab}"
            count=$((count + 1))
        done
    done

    echo "  Generated ${input_tab} with ${count} samples"
    echo "${count}"
}
