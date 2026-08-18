# Variant Calling Pipeline (vc)

A variant calling pipeline of within-host bacterial isolates with
population-specific reference genomes, adapted from [snippy](https://github.com/tseemann/snippy).
Population specific reference genomes were assembled de novo and annotated with [bakta](https://github.com/oschwengers/bakta).

This repo can be used two ways:

- **As a library** — an external orchestrator (e.g. a project's own
  Snakemake pipeline, with its own isolate discovery, reference mapping, and
  scheduling) calls the building-block scripts below directly. **This
  document describes that mode.**
- **Standalone** — the pipeline's original, self-contained SGE-orchestrated
  mode (`run_population.sh`), which discovers its own FASTQ inputs from a
  fixed directory convention and writes to a fixed output tree inside this
  repo. Still fully packaged in this repo — see
  [STANDALONE.md](STANDALONE.md).

Both modes share the same underlying variant-calling and table-generation
logic, so "Variant table semantics" and "Population-level site exclusions"
below apply regardless of which one you use.

## Requirements

- **conda or mamba**, to build the pipeline's tool environment (see Setup).
- **Python 3** (standard library only) to run `bin/*.py` — no `pip install`
  needed.
- Your own orchestrator and scheduling. This repo's library scripts don't
  require or assume SGE — each one just takes explicit file paths as CLI
  arguments and produces the named output file(s); how you invoke them
  (Snakemake, a plain shell loop, SGE, whatever) is entirely up to you.
  (Contrast with standalone mode, which does require SGE — see
  [STANDALONE.md](STANDALONE.md).)

## Setup

```bash
git clone <this-repo-url> vc
cd vc

# One conda env for every tool the pipeline calls (fastp, bwa, samtools,
# bcftools, freebayes, vt, snpEff, samclip, minimap2, bedtools, seqtk):
mamba env create -f environment.yml     # creates env "vc"
conda activate vc
```

There's no input-data symlink convention to set up in this mode (contrast
with standalone mode's `data/allreads`/`data/bakta_gff`) — every script below
takes explicit paths as CLI arguments, so your orchestrator passes whatever
paths your own project actually uses.

## Building blocks

Every stage (01–05) has a config.sh-independent entry point, so an external
orchestrator never needs to reimplement pipeline logic itself.

- **`bin/gbff_to_snpeff_gff.py <gbff>`** — flat CDS-only GFF from a bakta
  `.gbff`, printed to stdout. Consumed by both a snpEff database build and
  `04_generate_variant_tables.py`'s annotation lookup. Optionally applies
  `bin/gene_name_overrides.tsv` if present (cohort-specific locus-tag
  reconciliation — harmless/no-op for a fresh dataset).

- **`bin/build_snpeff_db.sh <ref_name> <gbff> <fna> <out_datadir>`** — builds
  one reference's snpEff database under `<out_datadir>/<ref_name>/`
  (`genes.gff`, `sequences.fa`, snpEff's build artifacts) plus a
  **self-contained** `snpeff.config` covering just that one genome. This is
  the config.sh-free twin of `01_build_snpeff_db.sh`'s per-sample loop body:
  that script appends every sample's genome entry to one shared
  `results/snpeff/snpeff.config` file, which only works safely because it
  runs samples sequentially in a single process. Here, each reference's
  config is independent, so building several references as separate
  parallel jobs (e.g. one Snakemake rule instance per reference) never races
  on a shared file — and no aggregation step is needed afterward, since the
  config this produces already has everything `snpEff build`/`snpEff ann`
  need for that one reference. Pass the resulting
  `<out_datadir>/<ref_name>/snpeff.config` as the `<snpeff_config>` argument
  below and to `04_generate_variant_tables.py`.

- **`bin/prepare_reference.sh <ref_fna> <out_dir>
  [minimap2_selfalign_opts]`** — indexes one reference for alignment:
  `bwa index` + `samtools faidx` on a copy of `<ref_fna>` at
  `<out_dir>/ref.fa`, plus the self-alignment repeat BED (see "Population-
  level site exclusions" below). Config.sh-free twin of
  `02_prepare_reference.sh`'s per-population loop body, minus the
  snpEff-config symlink step — not needed here, since no script reads a
  config or GFF back out of this directory; callers already have those
  paths from `build_snpeff_db.sh` above and pass them explicitly wherever
  they're needed.

- **`bin/normalize_annotate_vcf.sh <raw_vcf> <ref_fa> <ref_name> <out_vcf>
  <snpeff_config> <snpeff_datadir> <filt_min_qual>`** — normalize +
  snpEff-annotate a raw freebayes VCF into a per-isolate
  `snps.norm.annot.vcf`: `bcftools view` QUAL filter → `vt normalize` →
  `bcftools annotate` field cleanup → `snpEff ann`. This is the standalone,
  parameterized twin of standalone mode's `config.sh`-sourcing
  `normalize_annotate_vcf()` bash function — same logic, every path/param
  passed explicitly instead of read off global config. No FMT/DP filter is
  applied here; low-coverage variants are retained so that decision is made
  once, downstream, via the mosdepth-based masking `04_generate_variant_vcf.py`
  applies.

- **`bin/04_generate_variant_vcf.py`** (`--help` for full args) — builds the
  **canonical merged multi-sample VCF** from a `bcftools merge
  --missing-to-ref` of every isolate's `+setGT`-masked VCF: drops records
  overlapping a repeat BED, applies a mosdepth-based low-coverage mask, and
  optionally drops snp/mnp sites with excessive missingness. See
  "Population-level site exclusions" below for what it drops and why.
  `--pop` only affects visit-group labeling in the missingness QC TSV (see
  "Variant table semantics") — any string works as long as isolate names are
  prefixed with it.

  > **Note:** if every input VCF has zero records (e.g. far too little
  > sequencing depth across a whole population), `bcftools merge
  > --missing-to-ref` — which your orchestrator runs to build this script's
  > `--merged` input — can crash with a segfault, rather than producing an
  > empty (but valid) merged VCF. This is a `bcftools`-side edge case on
  > degenerate input; if you hit it, check per-isolate
  > `snps.norm.annot.vcf` record counts first (`grep -vc '^#'`) to confirm
  > whether any isolate has confident calls at all, and consider
  > substituting `bcftools view --print-header` for the merge in that case
  > (a real merge, just of an all-empty input).

- **`bin/04_generate_variant_tables.py`** (`--help` for full args) —
  genotype/annotation tables (all variants, SNPs, INDELs, disruptive
  variants) from that canonical VCF, with per-visit allele counts and
  frequencies — see "Variant table semantics" below. Re-runs `snpEff ann`
  internally, so it needs the reference's snpEff config/database available
  (whatever `build_snpeff_db.sh` above produced).

- **`bin/05_generate_snp_alignment.py`** (`--help` for full args) — SNP
  alignment FASTA from `04_generate_variant_tables.py`'s `snp.gt.tab` output.

- **`bin/pipeline_helpers.py`** — shared logic (callable-mask building,
  repeat-BED overlap, missing-GT detection, `group_isolates_by_visit`)
  importable directly if you're doing your own thing in Python rather than
  shelling out to the scripts above.

## Output layout

Library mode has no fixed output convention — every script takes its output
paths as explicit arguments, so directory names are entirely your
orchestrator's choice. The tree below illustrates one real layout (pa_promise's
`Snakefile.variantcalling`), to make the building blocks' outputs concrete:

```
<snpeff_datadir>/                        # arg to build_snpeff_db.sh
└── <ref_name>/
    ├── genes.gff
    ├── sequences.fa
    ├── snpeff.config                    # self-contained, one genome only
    └── snpEffectPredictor.bin           # (+ other snpEff build artifacts)

<ref_out_dir>/                           # one per reference; arg to prepare_reference.sh
├── ref.fa (+ .fai, .bwt, .pac, .sa, .amb, .ann)
└── ref.repeats.bed

<isolates_dir>/<isolate>/                # one per isolate
├── qc/                                  # fastp
│   └── <isolate>.fastp.{json,html}      # trimmed fastq.gz is typically temp(), removed after alignment
├── <isolate>.bam (+ .bai)               # typically temp() -- removed once freebayes/mosdepth finish
├── <isolate>.snps.raw.vcf               # freebayes
├── <isolate>.snps.norm.annot.vcf        # normalize_annotate_vcf.sh
├── <isolate>.setgt.vcf.gz (+ .tbi)      # typically temp() -- removed once the population merge finishes
└── <isolate>.per-base.bed.gz (+ .csi)   # mosdepth per-base depth

<population_out_dir>/                    # one per population/patient
├── aln/
│   ├── <pop>.merged.vcf.gz (+ .tbi)     # canonical merged VCF, 04_generate_variant_vcf.py
│   ├── <pop>.treetime.vcf.gz (+ .tbi)   # snp/mnp subset for TreeTime
│   ├── <pop>.snp_site_missingness.tsv   # per-site '-' distribution + EXCLUDED flag
│   └── <pop>.snps.aln.fa                # 05_generate_snp_alignment.py
└── pop_tables/                          # 04_generate_variant_tables.py
    ├── var.annot.tab
    ├── snp.{gt,annot}.tab
    ├── indel.{gt,annot}.tab
    └── disruptive.annot.tab
```

`.bam`/`.setgt.vcf.gz` are marked "typically temp()" because that's this
orchestrator's choice, not something the building-block scripts impose —
nothing stops a different orchestrator from keeping them.

## Population-level site exclusions

Two exclusions are applied once, in `04_generate_variant_vcf.py`, when
building the canonical merged VCF, so every downstream artifact (TreeTime
VCF, tables, SNP alignment) inherits an identical excluded set:

1. **Repeat exclusion** — any record overlapping the `--repeat-bed` mask
   (typically a minimap2 self-alignment of the reference) is dropped (all
   variant types).
2. **SNP-site missingness exclusion** — a record is dropped when its
   `INFO/TYPE` is purely `snp`/`mnp` **and**, after low-coverage masking, at
   least `--max-missing-frac` (a common choice: 0.80) of isolates are
   missing (`-`). Scoped to snp/mnp only, so **indels are never affected**.
   A per-site `-` distribution (overall + one `{VISIT}_TOTAL`/
   `{VISIT}_MISSING` column pair per visit group found in the isolate names
   — see "Variant table semantics" — with an `EXCLUDED` flag) is always
   written to `--missingness-tsv`, regardless of whether the drop itself is
   enabled.

Both stage-04-equivalent artifacts (the TreeTime VCF and the tables) derive
from the single canonical merged VCF, so they agree on every shared
`(CHR,POS,REF,ALT)` snp/mnp key and share this one repeat-exclusion step.

## Variant table semantics

`04_generate_variant_tables.py` emits one row per `(CHR, POS, REF, ALT)`
tuple. Per-isolate genotypes are encoded `1` (hom-ALT, passes
DP/QUAL/AO-frac thresholds), `0` (hom-REF, passes RO-frac threshold), or
`-` (het / missing / failed filter). Annotation columns come from
snpEff's `ANN=` field reconciled against the bakta GFF for stable
`LOCUS_TAG` / `GENE` / `PRODUCT`.

Frequency columns are **dynamic per-visit**, not a fixed pre/post pair.
`bin/pipeline_helpers.py`'s `group_isolates_by_visit(pop_name, isolates)`
parses each isolate's real name for a `^<pop_name>_V<visit>_` prefix (e.g.
`009-007_V1A_10` → visit `V1A`) and buckets it there; isolates whose names
don't match — wrong population prefix, or no parseable `_V<visit>_` token —
are excluded from every per-visit frequency column, reported as a
`WARNING` on stderr (not a hard failure). The annot tables (`var.annot.tab`,
`snp.annot.tab`, `indel.annot.tab`, `disruptive.annot.tab`) then get one
overall `FRQ` column followed by a `{VISIT}_CT`/`{VISIT}_FRQ` pair **per
visit label actually present in the data**, sorted alphabetically — so the
column set differs population to population depending on how many distinct
visits its isolates span (a population with only `V1`/`V1A` isolates gets
two CT/FRQ pairs; one spanning `V1`, `V1A`, `V1B`, `V2` gets four). There is
no hardcoded two-timepoint (pre/post) assumption anywhere in the table/VCF
generators.
