#!/usr/bin/env bash
#
# 05_run_generate_snp_alignment.sh
#
# Wrapper to build per-population SNP alignment FASTAs from pop_tables/snp.gt.tsv.
# Output: ${POP_DIR}/aln/${POP}.snps.aln.fa
#
# Usage:
#   bash 05_run_generate_snp_alignment.sh              # all populations
#   bash 05_run_generate_snp_alignment.sh <POP> [POP]  # specific population(s)
#   # SGE submission (single pop) -- VC_CONFIG must be passed explicitly:
#   # SGE spools the submitted script before executing it, so `dirname "$0"`
#   # inside the job would not point at this repo.
#   qsub -v POP=IA01,VC_CONFIG=/abs/path/to/vc/bin/config.sh \
#       05_run_generate_snp_alignment.sh

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

for POP in "${POPS[@]}"; do
    POP_DIR="$(get_pop_dir "${POP}")"
    GT_TAB="${POP_DIR}/pop_tables/snp.gt.tsv"
    ALN_DIR="${POP_DIR}/aln"
    OUT_FA="${ALN_DIR}/${POP}.snps.aln.fa"

    echo "=== Generating SNP alignment for ${POP} ==="
    echo "  Input:  ${GT_TAB}"
    echo "  Output: ${OUT_FA}"

    if [[ ! -f "${GT_TAB}" ]]; then
        echo "WARNING: ${GT_TAB} not found. Skipping ${POP}." >&2
        continue
    fi

    mkdir -p "${ALN_DIR}"

    conda run -n "${CONDA_ENV}" \
        python3 "${SCRIPT_DIR}/05_generate_snp_alignment.py" \
            --pop "${POP}" \
            --gt-tab "${GT_TAB}" \
            --out "${OUT_FA}"

    echo "=== Done: ${POP} ==="
    echo ""
done
