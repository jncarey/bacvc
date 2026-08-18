# Variant Calling Pipeline (vc)

A variant calling pipeline of within-host bacterial isolates with 
population-specific reference genomes, adapted from [snippy](https://github.com/tseemann/snippy). 
Population specific reference genomes were assembled de novo and annotated with [bakta](https://github.com/oschwengers/bakta).

## Requirements

- **SGE cluster** with `qsub` — stages 03–05 are submitted as SGE (array)
  jobs with `-hold_jid` dependency chaining. There is no non-SGE / local
  execution mode.
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
  uses them purely to find FASTQs, not to classify isolates (see below).
  `<sample>` becomes the isolate ID used everywhere downstream (BAM read
  groups, VCF sample names, table columns).
- **Isolate naming for frequency grouping**: independent of which `pre`/`post`
  directory a sample's reads were found in, `<sample>` must literally start
  with `<POP>_V<visit>_` (e.g. `009-007_V1A_10` for population/patient
  `009-007`, visit `V1A`) for `bin/pipeline_helpers.py`'s
  `group_isolates_by_visit()` to place it in a visit group. Isolates that
  don't match (wrong population prefix, or no parseable `_V<visit>_` token)
  are excluded from every per-visit frequency column — reported as a
  `WARNING` on stderr, not a hard failure (see "Variant table semantics"
  below). **The bundled `example/` dataset predates this convention** (its
  isolates are named `POPpr*`/`POPpo*`, not `POP_V*_*`) — see the Quickstart
  note below.
- **Reference genomes**: one [bakta](https://github.com/oschwengers/bakta)
  annotation per population reference isolate, providing both `<ref>.fna`
  (assembly FASTA) and `<ref>.gbff` (GenBank flatfile with CDS features) under
  `${BAKTA_GFF_DIR}/<ref>/` (`BAKTA_GFF_DIR` defaults to `data/bakta_gff`).
  Contig names in `<ref>.fna` must match the `LOCUS` names in `<ref>.gbff`
  (bakta output satisfies this by construction).
- **Population → reference mapping** declared in `POP_TO_REF` in
  `bin/config.sh` (e.g. `POP1 → POP1_autocycler`).

## Quickstart with the bundled example

A small demo dataset ships in `example/` — 4 real isolates (2 pre-, 2
post-treatment), reads subset to those overlapping a single ~400 kb window of
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
see `.gitignore`).

> **Stale as of the `group_isolates_by_visit` rewrite:** the demo's isolates
> are named `POPpr*`/`POPpo*`, which no longer matches the `POP_V<visit>_*`
> naming the table/VCF generators now group by (see "Input data contract").
> The Quickstart still *runs* end-to-end, but every isolate will land in
> `pipeline_helpers.py`'s `unmatched` bucket — you'll see a `WARNING: ...
> isolates didn't match population 'POP' + a parseable visit label` for all
> 4, and the output tables will carry only the overall `FRQ` column, no
> per-visit `V*_CT`/`V*_FRQ` breakdown. Renaming the demo's isolates (and
> re-deriving `example/results/`) to match is a follow-up, not done yet.

## Stages

| # | Script | Purpose |
|---|--------|---------|
| 01 | `bin/01_build_snpeff_db.sh` | Build a snpEff database per reference genome from its `.gbff`/`.fna`. Writes `snpeff/<ref>/` and a shared `snpeff/snpeff.config`. Uses `bin/gbff_to_snpeff_gff.py` to make a flat CDS-only GFF that snpEff and the table generator both consume; that script optionally applies `bin/gene_name_overrides.tsv` if present (cohort-specific locus-tag reconciliation — harmless/no-op for a fresh dataset). |
| 02 | `bin/02_prepare_reference.sh` | For each population, copy the reference `.fna` to `populations/<POP>/reference/ref.fa`, `bwa index`, `samtools faidx`, copy `genes.gff` → `ref.gff`, symlink the snpEff config, and minimap2-self-align the reference to emit `reference/ref.repeats.bed` (repetitive regions excluded in stage 04). |
| 03 | `bin/03_worker.sh` | Per-isolate SGE array task: fastp QC → `bwa mem` → `samclip` → `samtools sort/fixmate/markdup` → `freebayes` → normalize + snpEff-annotate (via `config.sh`'s `normalize_annotate_vcf()` bash function — QUAL filter → `vt normalize` → `snpEff ann`). Output: `populations/<POP>/variants/<sample>/snps.norm.annot.vcf`. `bin/normalize_annotate_vcf.sh` is a standalone, parameterized twin of that same function's logic, for callers that don't source `config.sh` (e.g. an external Snakemake pipeline calling this repo's scripts directly — see "Reuse as a library" below). |
| 04a | `bin/04_run_generate_variant_vcf.sh` + `04_generate_variant_vcf.py` | After all workers finish, build the **canonical merged multi-sample VCF** (`aln/<POP>.merged.vcf.gz`): per-sample `+setGT` masking → `bcftools merge` → drop records overlapping `ref.repeats.bed` → mosdepth low-coverage mask → drop snp/mnp sites with ≥80% isolates missing (writes a per-site `-` distribution to `aln/<POP>.snp_site_missingness.tsv`, with a column pair per visit group — see "Variant table semantics"). Derive the TreeTime VCF (`aln/<POP>.treetime.vcf.gz`, snp/mnp subset) from it. |
| 04b | `bin/04_run_generate_variant_tables.sh` + `04_generate_variant_tables.py` | Held on 04a. Derive population tables from the same canonical VCF (`bcftools norm -m-` → `snpEff ann` → genotype/annotation tables for all/SNP/INDEL/disruptive variants, with per-visit allele counts and frequencies — see "Variant table semantics"). Re-runs snpEff, so the reference's `snpeff.config` entry from stage 01 must be present (the wrapper fails loud otherwise). |
| 05 | `bin/05_run_generate_snp_alignment.sh` + `05_generate_snp_alignment.py` | Held on 04b. Build the per-population SNP alignment FASTA from `pop_tables/snp.gt.tab`. |

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

### Population-level site exclusions (canonical VCF, stage 04a)

Two population-level exclusions are applied once when building the canonical
VCF, so every downstream artifact (TreeTime VCF, tables, SNP alignment) inherits
an identical excluded set:

1. **Repeat exclusion** — any record overlapping the minimap2 self-alignment
   mask (`reference/ref.repeats.bed`) is dropped (all variant types).
2. **SNP-site missingness exclusion** — a record is dropped when its
   `INFO/TYPE` is purely `snp`/`mnp` **and**, after low-coverage masking, at
   least `SNP_SITE_MAX_MISSING_FRAC` (default 0.80) of isolates are missing
   (`-`). Scoped to snp/mnp only, so **indels are never affected**. A per-site
   `-` distribution (overall + one `{VISIT}_TOTAL`/`{VISIT}_MISSING` column
   pair per visit group found in the isolate names — see "Variant table
   semantics" — with an `EXCLUDED` flag) is always written to
   `aln/<POP>.snp_site_missingness.tsv`, regardless of the toggle.

`bin/run_population.sh` is the entry point: it discovers FASTQ inputs,
generates `populations/<POP>/input.tab`, submits the worker array, then submits
stage 04a (`-hold_jid` on the worker array), 04b (held on 04a) and 05 (held on
04b).

## Reuse as a library

`run_population.sh`'s SGE array + `-hold_jid` chaining is one way to drive
this pipeline, not the only one. A project with its own orchestrator (e.g. a
Snakemake pipeline with its own isolate-discovery, its own reference-per-patient
mapping, its own scheduling) can call the underlying building blocks directly
instead, without sourcing `config.sh` or going through `run_population.sh` at
all:

- `bin/gbff_to_snpeff_gff.py <gbff>` — flat CDS-only GFF from a bakta `.gbff`.
- `bin/normalize_annotate_vcf.sh <raw_vcf> <ref_fa> <ref_name> <out_vcf>
  <snpeff_config> <snpeff_datadir> <filt_min_qual>` — the standalone twin of
  stage 03's `normalize_annotate_vcf()` bash function; takes every path/param
  as an explicit argument instead of reading them off `config.sh` globals.
- `bin/04_generate_variant_vcf.py` and `bin/04_generate_variant_tables.py` —
  plain argparse CLIs (see `--help`); `--pop` only affects visit-group
  labeling (see "Variant table semantics"), so any string works as long as
  isolate names are prefixed with it.
- `bin/05_generate_snp_alignment.py` — likewise a standalone CLI.
- `bin/pipeline_helpers.py` — shared logic (callable-mask building, repeat-BED
  overlap, missing-GT detection, `group_isolates_by_visit`) importable
  directly if a caller is doing its own thing in Python rather than shelling
  out to the scripts above.

The external orchestrator is responsible for its own scheduling (SGE, local,
whatever) and its own input discovery; this repo's contract in that mode is
just "given these files/paths, produce this pipeline stage's output."

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
  `SNP_SITE_MAX_MISSING_FRAC`).
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
        │   └── .done
        ├── aln/
        │   ├── <POP>.merged.vcf.gz        # canonical merged VCF (all types) +tbi
        │   ├── <POP>.treetime.vcf.gz      # snp/mnp subset for TreeTime +tbi
        │   ├── <POP>.snp_site_missingness.tsv  # per-site '-' distribution + EXCLUDED flag
        │   └── <POP>.snps.aln.fa          # SNP alignment FASTA (stage 05)
        ├── pop_tables/
        │   ├── var.annot.tab
        │   ├── snp.{gt,annot}.tab
        │   ├── indel.{gt,annot}.tab
        │   └── disruptive.annot.tab
        └── logs/
```

`results/` is git-ignored in its entirety — nothing under it is ever committed.

## Variant table semantics

`04_generate_variant_tables.py` emits one row per `(CHR, POS, REF, ALT)`
tuple. Per-isolate genotypes are encoded `1` (hom-ALT, passes
DP/QUAL/AO-frac thresholds), `0` (hom-REF, passes RO-frac threshold), or
`-` (het / missing / failed filter). Annotation columns come from
snpEff's `ANN=` field reconciled against the bakta GFF for stable
`LOCUS_TAG` / `GENE` / `PRODUCT`.

Frequency columns are **dynamic per-visit**, not the fixed pre/post pair
older versions of this pipeline used. `bin/pipeline_helpers.py`'s
`group_isolates_by_visit(pop_name, isolates)` parses each isolate's real
name for a `^<pop_name>_V<visit>_` prefix (e.g. `009-007_V1A_10` → visit
`V1A`) and buckets it there; isolates that don't match are excluded from
every group (see "Input data contract"). The annot tables (`var.annot.tab`,
`snp.annot.tab`, `indel.annot.tab`, `disruptive.annot.tab`) then get one
overall `FRQ` column followed by a `{VISIT}_CT`/`{VISIT}_FRQ` pair **per
visit label actually present in the data**, sorted alphabetically — so the
column set differs population to population depending on how many distinct
visits its isolates span (a population with only `V1`/`V1A` isolates gets
two CT/FRQ pairs; one spanning `V1`, `V1A`, `V1B`, `V2` gets four). This
replaces the old fixed `PR_CT`/`PO_CT`/`PR_FRQ`/`PO_FRQ` columns entirely —
there is no longer a hardcoded two-timepoint assumption anywhere in the
table/VCF generators.
