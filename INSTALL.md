# TARVa — installation and use

TARVa analyses RNA variants in matched RNA-seq and WGS data: it calls variants in
both, keeps the RNA-specific ones, builds a database from them, and runs the
global and local analyses.

It works on **any matched RNA + WGS dataset with two conditions**. Nothing is
specific to a cohort: condition labels, covariate columns and individual
identifiers are read from your own metadata file and carried through under their
own names.

---

## Which install do you need?

| You want to | Install |
|---|---|
| Build a database and run the analyses, from variant calls you already have | **Analysis only** — pip, works anywhere |
| Run the whole thing from FASTQs | **Full pipeline** — conda, adds the aligners and callers |

The analysis half is pure Python. The upstream half shells out to `STAR`, `bwa`,
`samtools`, `gatk`, `stringtie`, `gffcompare`, `gffread`, `fastqc` and
`trimmomatic`, which pip cannot install.

### Analysis only

Needs Python ≥ 3.7 and nothing else. (3.7 is end-of-life — a newer interpreter
is better where you have the choice, but nothing here requires one.)

```bash
pip install git+https://github.com/kkingwoo/TARVa
```

To work on the code, clone first and install in place:

```bash
git clone https://github.com/kkingwoo/TARVa
cd TARVa
pip install -e .
```

### Full pipeline

```bash
git clone https://github.com/kkingwoo/TARVa
cd TARVa
conda env create -f environment.yml
conda activate tarva
```

That installs TARVa itself as well. Check it:

```bash
tarva --help
```

Every command is also reachable as `python -m tarva.cli ...`, which is useful on
a cluster where the console script is not on `PATH`.

---

## What you supply

**A metadata CSV**, one row per individual. TARVa reads it verbatim — every
column is kept under its own name — and you tell it which two columns matter:

- `--id-column` the individual identifier. Must be unique and non-blank.
- `--condition-column` the grouping variable. Must have **exactly two** distinct
  values. Which one is the case group is `--case-condition`; without it the
  alphabetically first is used, and either way the choice is recorded.

Every other column is available to your own scripts straight from the
`individual` table. Covariate and bias testing is deliberately **not** part of
the workflow — which covariates exist and which test suits each is specific to
your cohort, so it belongs in a script you control.

**A reference annotation** (GTF) and, for the full pipeline, a reference genome
plus the FASTQs.

---

## Running it

### 1. Build the database

```bash
tarva dbpipe all \
    --gtf merged.gtf \
    --metadata metadata.csv \
    --rna-unique rna_unique/ \
    --stringtie-dir assembly/ \
    --db TARVa.db
```

Four steps in one command — annotation, then the cohort, then the variant calls,
then expression. The order is required: calls are resolved to genes the
annotation wrote, and to individuals the cohort step wrote.

Each step also runs on its own (`tarva dbpipe annotate --help`, and so on), which
is what you want when re-running one part. `--overwrite` recomputes instead of
refusing to write over existing rows.

### 2. Run the analyses

```bash
tarva analysispipe all \
    --db TARVa.db \
    --output results/ \
    --min-depth 20
```

Six steps: the per-gene and per-call-type tests with the gene subsets, the
threshold sensitivity check, expression against variant load, unique versus
common positions, the local per-position narrowing, and the filtering cascade.
Every CSV the figures read lands in `--output`.

### 3. Draw the figures

```bash
python figures/Fig1_filter_cascade.py --input-dir results/
python figures/Fig5_uniqueness_forest.py --input-dir results/
```

The other two are R Markdown:

```r
rmarkdown::render("figures/Fig3_calltype_boxplot.Rmd",
                  params = list(input_dir = "results/"))
rmarkdown::render("figures/Fig4_tpm_variant_correlation.Rmd",
                  params = list(input_dir = "results/"))
```

R figures need `tidyverse` and `patchwork`.

No figure computes anything. Every number plotted comes from a step's CSV, so a
figure and the database cannot disagree.

### The full pipeline, from FASTQs

Upstream of the database, in order, each with its own `--help`:

```
tarva sort-input       sort FASTQs into the expected tree
tarva prep-reference   decompress, index, build STAR/BWA indexes
tarva qc               read lengths, FastQC, Trimmomatic
tarva stringtiepipe    first-pass alignment, 3-round StringTie assembly
tarva wgseqpipe        WGS alignment, GATK preprocessing, PON, variants
tarva rnaseqpipe       two-pass RNA alignment, preprocessing, variants
tarva analysispipe rna-unique   RNA variants with no WGS counterpart
```

`rna-unique` produces the input `dbpipe` consumes. It is not part of
`analysispipe all`, which only runs steps that read the database.

---

## Reading the results

Everything is in the SQLite database, queryable directly:

```bash
sqlite3 TARVa.db "SELECT * FROM run_meta ORDER BY written_at"
```

`run_meta` records every threshold, cohort size, input path and decision each
step made — including which p-value basis a step used and whether it ran on a
filtered gene set or fell back to all genes. If a number in a figure needs
explaining, it is in there.

### When a step reports nothing

This is a normal outcome, not a failure. A step adjusts its p-values, and adopts
a result only when enough of the tests reach the threshold; when neither the
adjusted nor the nominal values do, it says so and selects nothing rather than
returning a handful of borderline tests as a result. The message names the usual
causes: too few individuals per condition, a dataset without the signal, or a
`--min-depth` strict enough to have excluded the supporting calls.


---

## License

MIT. See `LICENSE`.
