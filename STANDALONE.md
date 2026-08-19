# Standalone mode (`run_population.sh`)

This is the pipeline's original, self-contained way to run: it discovers its
own FASTQ inputs from a fixed directory convention, submits an SGE array job
for the per-isolate stage plus `-hold_jid`-chained downstream stages, and
writes output to a fixed directory tree inside this repo. It requires an SGE
cluster and nothing else about your project's own layout or scheduling.

If you're instead calling this repo's scripts directly from your own
orchestrator (e.g. a project-specific Snakemake pipeline with its own
isolate discovery and scheduling), see [README.md](README.md) — that's the
library interface, documented separately from this file. Both modes share
the same underlying variant-calling/table-generation logic and semantics
(see README.md's "Variant table semantics" and "Population-level site
exclusions" sections, referenced below rather than duplicated here).

## Requirements

- **SGE cluster** with `qsub` — stages 03–05 are submitted as SGE (array)
  jobs with `-hold_jid` dependency chaining. There is no non-SGE / local
  execution mode for this mode specifically.
- **conda or mamba**, to build the pipeline's tool environment (see Setup).
- **Python 3** (standard library only) to run `bin/*.py` — no `pip install`
  needed.

## Setup

```bash
git clone <this-repo-url> vc
cd vc

# One conda env for every tool the pipeline calls (fastp, bwa, samtools,
# bcftools, freebayes, vt, snpEff, samclip, minimap2, bedtools, seqtk):
mamba env create -f environment.yml     # creates env "vc"
conda activate vc

# Point the pipeline at your own input data (see "Input data contract"
# below). data/allreads and data/bakta_gff are git-ignored. If your reads
# already sit in a single directory laid out as <POP>{pre,post}/, symlink it
# wholesale:
ln -s /path/to/your/reads      data/allreads
ln -s /path/to/your/bakta_gff  data/bakta_gff

# Otherwise, make data/allreads/data/bakta_gff plain directories and symlink
# (or copy) in one <POP>{pre,post}/ or <ref>/ entry at a time -- this is what
# the Quickstart below does for the bundled example.
```

`bin/config.sh`'s `CONDA_ENV` points at the `vc` env created above by default;
every stage runs its tools via `conda run`/`conda activate` against it.

## Input data contract

- **Paired-end FASTQ reads** under
  `${FASTQ_BASE_DIR}/<POP>{pre,post}/<sample>_{1,2}.fastq.gz`
  (`FASTQ_BASE_DIR` is set in `bin/config.sh`, defaults to `data/allreads`).
  `<POP>` is a short population code (e.g. `POP1`); `pre`/`post` are literal
  suffixes selecting the *discovery* subdirectory only — `run_population.sh`
  uses them purely to find FASTQs, not to classify isolates. `<sample>`
  becomes the isolate ID used everywhere downstream (BAM read groups, VCF
  sample names, table columns).
- **Isolate naming for frequency grouping** works identically to library
  mode — see README.md's "Variant table semantics" for the `<POP>_V<visit>_`
  naming requirement. It's independent of which `pre`/`post` directory a
  sample's reads were found in — the bundled `example/` dataset's isolates
  are named `POP_V1_*`/`POP_V2_*` accordingly, even though its reads still
  live under the `POPpre`/`POPpost` discovery directories (see Quickstart
  below).
- **Reference genomes**: one [bakta](https://github.com/oschwengers/bakta)
  annotation per population reference isolate, providing both `<ref>.fna`
  (assembly FASTA) and `<ref>.gbff` (GenBank flatfile with CDS features) under
  `${BAKTA_GFF_DIR}/<ref>/` (`BAKTA_GFF_DIR` defaults to `data/bakta_gff`).
  Contig names in `<ref>.fna` must match the `LOCUS` names in `<ref>.gbff`
  (bakta output satisfies this by construction).
- **Population → reference mapping** declared in `POP_TO_REF` in
  `bin/config.sh` (e.g. `POP1 → POP1_autocycler`).

## Quickstart with the bundled example

A small demo dataset ships in `example/` — 4 real isolates (2 at visit `V1`,
2 at visit `V2`), reads subset to those overlapping a single ~400 kb window of
the reference (`contig_1:1-400000`) so per-base depth in that window matches
the isolates' real sequencing depth (~65,000-80,000 read pairs/isolate,
tens of MB total), paired with the real (full-size) bakta reference. The rest
of the genome is uncovered by design (masked as missing downstream), but the
covered window has real depth, so the demo exercises genuine variant calling
and the QUAL/DP thresholds, not just plumbing. It runs the entire pipeline
end-to-end in a few minutes, so you can validate your environment/cluster
setup before pointing the pipeline at your own data.

```bash
conda activate vc

# Wire the bundled example into the input locations the pipeline expects:
ln -s "$(pwd)/example/reads/POPpre"      data/allreads/POPpre
ln -s "$(pwd)/example/reads/POPpost"     data/allreads/POPpost
ln -s "$(pwd)/example/bakta_gff/POP_autocycler" data/bakta_gff/POP_autocycler

bash bin/run_population.sh POP
```

`bin/config.sh`'s `POP_TO_REF` already maps `POP → POP_autocycler`.
Monitor with `qstat -u $USER`; once the `aln_POP` job finishes, check
`results/populations/POP/{aln,pop_tables}/` for output (never committed —
see `.gitignore`). `example/results/POP/` holds a pre-generated copy of that
output for reference, so you can compare without running anything.

## Stages

| # | Script | Purpose |
|---|--------|---------|
| 01 | `bin/01_build_snpeff_db.sh` | Build a snpEff database per reference genome from its `.gbff`/`.fna`. Writes `snpeff/<ref>/` and a shared `snpeff/snpeff.config`. Uses `bin/gbff_to_snpeff_gff.py` to make a flat CDS-only GFF that snpEff and the table generator both consume; that script optionally applies `bin/gene_name_overrides.tsv` if present (cohort-specific locus-tag reconciliation — harmless/no-op for a fresh dataset). |
| 02 | `bin/02_prepare_reference.sh` | For each population, copy the reference `.fna` to `populations/<POP>/reference/ref.fa`, `bwa index`, `samtools faidx`, copy `genes.gff` → `ref.gff`, symlink the snpEff config, and minimap2-self-align the reference to emit `reference/ref.repeats.bed` (repetitive regions excluded in stage 04). |
| 03 | `bin/03_worker.sh` | Per-isolate SGE array task: fastp QC → `bwa mem` → `samclip` → `samtools sort/fixmate/markdup` → `freebayes` → normalize + snpEff-annotate (via `config.sh`'s `normalize_annotate_vcf()` bash function — `bcftools view` QUAL filter → `vt normalize` → `bcftools annotate` field cleanup → `snpEff ann`) → `mosdepth` → `samtools flagstat` → per-isolate QC flag (`qc_flag_isolate.py`, README.md). Output: `populations/<POP>/variants/<sample>/snps.norm.annot.vcf` (+ QC files below). `bin/normalize_annotate_vcf.sh`, documented in README.md, is a standalone, parameterized twin of that same function's logic. |
| 04a | `bin/04_run_generate_variant_vcf.sh` + `04_generate_variant_vcf.py` | After all workers finish, build the **canonical merged multi-sample VCF** (`aln/<POP>.merged.vcf.gz`) — see README.md's "Population-level site exclusions" for what this step drops and why. Derive the TreeTime VCF (`aln/<POP>.treetime.vcf.gz`, snp/mnp subset) from it. |
| 04b | `bin/04_run_generate_variant_tables.sh` + `04_generate_variant_tables.py` | Held on 04a. Derive population tables from the same canonical VCF — see README.md's "Variant table semantics". Re-runs snpEff, so the reference's `snpeff.config` entry from stage 01 must be present (the wrapper fails loud otherwise). |
| 05 | `bin/05_run_generate_snp_alignment.sh` + `05_generate_snp_alignment.py` | Held on 04b. Build the per-population SNP alignment FASTA from `pop_tables/snp.gt.tsv`. |

Both stage-04 artifacts derive from the single canonical merged VCF, so the
TreeTime VCF and the tables agree on every shared `(CHR,POS,REF,ALT)` snp/mnp
key and share one repeat-exclusion step.

**Note:** if a population has zero confident variant calls across every
isolate (e.g. far too little sequencing depth), `bcftools merge` in stage 04a
can crash with a segfault on the resulting empty per-sample VCFs, rather than
producing an empty (but valid) merged VCF. This is a `bcftools`-side edge
case on degenerate input, not something this pipeline can paper over --
if you hit it, check per-isolate `variants/<sample>/snps.norm.annot.vcf`
record counts first (`grep -vc '^#'`) to confirm whether any isolate has
confident calls at all.

`bin/run_population.sh` is the entry point: it discovers FASTQ inputs,
generates `populations/<POP>/input.tab`, submits the worker array, then submits
stage 04a (`-hold_jid` on the worker array), 04b (held on 04a) and 05 (held on
04b).

## Configuration

All knobs live in `bin/config.sh`:

- Input paths (`FASTQ_BASE_DIR`, `BAKTA_GFF_DIR`) — both default to
  symlinks under `vc/data/` (`vc/data/allreads`, `vc/data/bakta_gff`) —
  and the output root `RESULTS_DIR` (defaults to `vc/results/`, with
  `SNPEFF_DIR` and the populations tree underneath it).
- `POP_TO_REF` mapping.
- Conda env: `CONDA_ENV` (defaults to the single `vc` env from
  `environment.yml`; see Setup).
- Tool parameters: `BWA_THREADS`, `SAMCLIP_MAX_SOFT`, freebayes thresholds
  (`FB_*`), VCF pre-filter (`FILT_MIN_QUAL`, `MASK_MIN_DP`), table filters
  (`VT_*`), repeat exclusion (`EXCLUDE_REPEATS`, `MINIMAP2_SELFALIGN_OPTS`),
  SNP-site missingness exclusion (`EXCLUDE_HIGH_MISSING_SNP_SITES`,
  `SNP_SITE_MAX_MISSING_FRAC`), per-isolate QC gate (`MIN_MEAN_DEPTH`,
  `MIN_MAPPED_PCT` — see README.md's `qc_flag_isolate.py` entry).
- SGE submission (`MAX_CONCURRENT`, `SGE_MEMORY`, `SGE_WALLTIME`).
- Cleanup toggles (`REMOVE_TRIMMED`, `REMOVE_BAMS`).

## Usage

```bash
conda activate vc

# Run all populations (submits SGE array per population)
bash bin/run_population.sh

# Run specific populations
bash bin/run_population.sh POP1 POP2

# Dry run (prints what would be submitted)
bash bin/run_population.sh --dry-run POP1
```

Stages 01 and 02 are skipped automatically when their outputs already exist.

## Output layout

```
vc/results/                            # ${RESULTS_DIR}
├── snpeff/
│   ├── snpeff.config
│   └── <ref>/                         # snpEff DB per reference
└── populations/
    └── <POP>/
        ├── input.tab                  # sample_id <tab> R1 <tab> R2
        ├── reference/                 # ref.fa (+bwa index +fai), ref.gff, snpeff.config, ref.repeats.bed
        ├── trimmed_reads/<sample>/    # removed if REMOVE_TRIMMED=true
        ├── variants/<sample>/
        │   ├── snps.raw.vcf
        │   ├── snps.norm.annot.vcf
        │   ├── <sample>.per-base.bed.gz   # mosdepth per-base depth
        │   ├── <sample>.per-base.bed.gz.csi
        │   ├── <sample>.mosdepth.summary.txt  # mosdepth mean depth (kept, not ephemeral)
        │   ├── <sample>.flagstat.txt      # samtools flagstat, before BAM cleanup
        │   ├── <sample>.qc_flag.tsv       # qc_flag_isolate.py verdict
        │   └── .done
        ├── aln/
        │   ├── <POP>.merged.vcf.gz        # canonical merged VCF (all types) +tbi
        │   ├── <POP>.treetime.vcf.gz      # snp/mnp subset for TreeTime +tbi
        │   ├── <POP>.snp_site_missingness.tsv  # per-site '-' distribution + EXCLUDED flag
        │   └── <POP>.snps.aln.fa          # SNP alignment FASTA (stage 05)
        ├── pop_tables/
        │   ├── var.annot.tsv
        │   ├── snp.{gt,annot}.tsv
        │   ├── indel.{gt,annot}.tsv
        │   └── disruptive.annot.tsv
        └── logs/
```

`results/` is git-ignored in its entirety — nothing under it is ever committed.
