# `tarva dbpipe` — building the TARVa database

**Last updated 2026-10-15.** Operational reference for the database build: what to
run, what each step writes, and what is still undecided. The *reasoning* behind the
schema lives in [README_database_plan.md](../README_database_plan.md); this file is
the how.

---

## 1. Steps

Five steps, run in order. Each is the **single writer** of its own tables, which is
what makes the build restartable one step at a time and keeps two steps from
disagreeing about a column.

| Step | Writes | Indexes | Reads | Status |
|---|---|---|---|---|
| `annotate` | `gene`, `transcript`, `exon` | `gene_locus_idx`, `gene_name_idx`, `trans_gene_idx`, `trans_locus_idx`, `trans_span_idx` | `merged.gtf`, attribute TSV, assembly report | **written; verified on mouse 2026-10-15** |
| `load-samples` | `sample` | `sample_rid_idx` | `metadata.csv` | **written** |
| `load-calls` | `site`, `call` | `site_gene_idx`, `call_rid_idx`, `call_type_idx` | `<id>_rna_unique.vcf`, `sample`, `gene`, `run_meta` | **written** |
| `load-expression` | `expression` | `expr_tid_idx` | `t_data.ctab` per sample, `sample`, `transcript` | **written** |

`gene_test` and `gene_class` are declared in [schema.py](schema.py) but written by
**`tarva analysispipe gene-tests`**, not by a `dbpipe` step. That step loads nothing —
it reads the tables above and computes — so it belongs with the analyses. `schema.py`
is the database's schema rather than dbpipe's private one, and the
one-writer-per-table rule is about tables, not packages.

Each step builds the indexes for its own tables as it finishes them, so no
table is left unindexed because a later step has not run.

"Verified" means `sanity_checks/verify_annotate.py` passes against the built
database — exon arithmetic over every transcript, `exon_num` contiguity and
`trans_offset` accumulation, strand ordering, and referential integrity that
`foreign_keys=OFF` does not enforce at insert time. The mouse run loaded 43,037
genes, 145,447 transcripts and 1,634,917 exons over 20 contigs, with an 84.5%
symbol-join biotype rate.

`run_meta` is written by every step, each owning keys under its own prefix
(`annotate.*`, `load_calls.*`, …).

## 2. Schema

Eleven tables. Full DDL in [schema.py](schema.py).

| Table | Key | Holds |
|---|---|---|
| `individual` | `individualID` | the metadata CSV, verbatim — one row per individual, every column under its own name |
| `gene` | `ref_id` | reference gene or novel assembly locus; `gene_span`, `biotype`, `novel`, `n_transcripts` |
| `transcript` | `tid` | `ref_id`, `locus_id`, span, spliced length, biotype, TSL, CCDS |
| `exon` | `(tid, exon_num)` | coordinates plus `trans_offset` |
| `site` | `(chrom, pos)` | one row per distinct called position: `ref_base` (the single genomic base), `strand`, `ref_id` |
| `call` | `(chrom, pos, individualID)` | one row per individual per position: `ref`, `alt1`, `a1_ad`, `dp`, `a1type`, `filter` |
| `expression` | `(individualID, tid)` | `FPKM`, `TPM` |
| `gene_test` | `(ref_id, a1type)` | test statistics, `complete_in`, observed/padded counts, and the `sum_ad`/`sum_dp` pair step 4 aggregates |
| `gene_class` | `ref_id` | `edit_dist`, `gene_group`, per-condition alternate-allele counts |
| `calltype_test` | `(gene_set, a1type)` | step 1b's across-gene test per call type, with arm sizes |
| `run_meta` | `key` | thresholds, cohort sizes, versions, input paths, per-step counts |

`a1prop`, `single_iso` (= `n_transcripts = 1`) and the legacy
`All_Info_TopGenes_tab` column names are **views** over these, not stored columns,
so they cannot drift from `a1_ad`, `dp`, `gene_span` and `n_transcripts`.

`n_transcripts` *is* stored, despite being a count over `transcript`, because it
is step 3 of the resolution rule below and `load_calls` consults it once per
contested position — a column read beats a `GROUP BY` over 147,908 transcripts.

### Sizing, from the published human database

| | rows |
|---|---|
| `call` | 37,797,345 |
| `site` | 9,613,214 |
| calls per site | 3.93 |
| calls per sample | 640,633 |

Roughly 3 GB, against the published database's ~80 GB — the difference being the
441,809,342-row per-transcriptome-base table this schema replaces with
`exon.trans_offset`.

## 3. Conventions the code depends on

- **`gene_span = gene_end - gene_start + 1`**, the contiguous genomic span over all of
  a gene's transcripts, introns included. This is the `a1prop` denominator and it
  reproduces the published values. Coordinates are the GTF's own 1-based inclusive
  ones throughout, so an exon's length is `exon_end - exon_start + 1`.
- **`trans_offset`** is the number of transcript bases preceding an exon in transcript
  order, so a transcript coordinate is one addition:
  `+ trans_pos = trans_offset + (pos - exon_start) + 1`,
  `- trans_pos = trans_offset + (exon_end - pos) + 1`.
- **Attributes are parsed by key, never by position.** A reference-matched transcript
  line ends in `ref_gene_id` and a novel one ends in `transcript_id`, so positional
  indexing gives novel transcripts a `ref_id` equal to their own `tid`.
- **The attribute TSV joins on the unversioned accession.** `merged.gtf` and the
  attribute table routinely disagree on version, and a versioned join fails silently
  for every accession whose version moved. `--min-attrs-match` turns that into a hard
  failure instead of a database of NULL biotypes.
- **Version stripping is conservative**: a trailing `.N` comes off only when the stem
  has no dot and ends in a digit, so `ENSMUSG00000000001.5` and `XM_006495550.4` are
  stripped while `MSTRG.203` and `MSTRG.203.1` are left alone. `--keep-versions`
  disables it.
- **The GFF3 ID namespace comes off first, and unconditionally.** StringTie lifts
  the reference GFF3's `ID=` attributes verbatim, so `merged.gtf` carries
  `ref_gene_id "gene-Xkr4"` and `transcript_id "rna-XM_006495550.5"`. The prefix is
  a namespace, not part of the accession, and it is also what stops version
  stripping from firing at all — `-` is not a word character, so without this
  `tid_version` is NULL for every transcript in the assembly. The strip is *not*
  keyed to feature type, because the prefix does not encode it: an NCBI gene with
  no separate RNA record arrives as `transcript_id "gene-Gm57489"`, the gene
  feature serving as its own transcript. Verified on this assembly: the prefix came
  off all 139,383 `ref_gene_id` values, which is every transcript that has one.

  A stripped identifier is checked for a *remaining* GFF3 namespace only
  (`cds-`, `exon-`, a doubled `gene-`), never for a general `letters-` pattern —
  mouse symbols legitimately look like that. MGI names nuclear tRNA genes
  `n-TKctt7`, and 510 of this assembly's genes are such symbols.
- **Records this build cannot key are dropped and counted**, not fatal. No
  `transcript_id` means the record cannot be joined or coordinate-mapped; no
  `gene_id` means it was never placed at a locus, and `locus_id` is what the
  resolution rule groups on. An empty attribute value counts as absent — NCBI
  writes `transcript_id ""` on gene lines, which `str.extract` returns as an empty
  string rather than NaN. The `annotate.n_dropped_*` keys in `run_meta` are the
  only record of how many.
- **The assembly report's columns are read from its own header line**, not assumed
  positionally. NCBI closes the comment block with
  `# Sequence-Name<TAB>Sequence-Role<TAB>…`, so the names are in the file; they are
  normalized to lowercase-with-underscores and matched by name. Same rule as GTF
  attributes, same reason: a positional layout breaks silently the moment a column is
  added or reordered, and a missing name should say which rather than misalign
  everything after it.
- **The contig set is decided at `annotate` and recorded**, in `run_meta` as
  `annotate.contigs_kept`, because the choice is database-wide. `load_calls` reads
  that key rather than re-deriving the rule.
- **`call.alt1` can be `*`**, VCF's spanning-deletion allele, which HaplotypeCaller
  emits where an upstream deletion overlaps the position. 192 of this cohort's 235,635
  calls. They are kept rather than skipped — the row is a true statement that a call
  existed there — and `a1type` is NULL for them, which excludes them from every
  `GROUP BY a1type` analysis. Do not compute `a1_ad / dp` over them: that ratio is the
  fraction of reads supporting "overlapped by a deletion", not an editing proportion.
- **`TPM` is normalized over the whole quantified transcriptome**, before any filtering.
  TPM is a within-sample composition, so the denominator is every row of that sample's
  ctab; dropping unannotated transcripts first would renormalize over a subset and
  change every value silently.
- **`site.ref_base` and `call.ref` are different things**, as they are in
  `All_Info_TopGenes_tab`. `ref_base` is the single genomic base *at* `pos` — the first
  character of the VCF REF, and the only part that is a function of the coordinate.
  `ref` is the call's own REF allele in full. For an indel those differ, because REF
  spans the affected bases: one sample's deletion is `REF=TGG ALT=T` while another's
  insertion at the same coordinate is `REF=T ALT=TAA`. They share the anchor base and
  nothing else. Keeping `ref` on `call` is what makes an indel row interpretable —
  `alt1='T'` beside `ref_base='T'` reads as no variant at all.
- **The whole metadata file becomes the `individual` table, verbatim.** Every column
  keeps its own name — ROSMAP's clinical columns, the mouse cohort's `age`,
  `WGS`/`RNA`/`*_barcode` so a row traces back to its FASTQs. Nothing is renamed,
  reshaped or dropped: loading is not the place to narrow the data, and which
  covariates to *analyse* is named at the analysis step, where that knowledge lives.
- **It is `individual`, not `sample`.** One row per individual, who contributes a WGS
  sample and an RNA sample; legacy's `sample_tab` named the wrong thing. The
  identifier stays `individualID` in `individual`, `call` and `expression` alike, so
  the `rid` / `rids` / `individualID` three-way muddle in §4 of the plan does not
  survive into the rebuild.
- **`--id-column` and `--condition-column` identify, they do not rename.** Both names
  go into `run_meta`, and consumers read them back through
  `schema.individual_columns` rather than assuming. They are quoted when interpolated
  into SQL, since a CSV may legitimately have a column called `group`.
- **`individual` is the one table not declared in `schema.py`** — its columns belong to
  the dataset, so declaring them would mean generating DDL per dataset. A `UNIQUE`
  index on the identifier column supplies what a declared primary key would: one row
  per individual, and `call.individualID`'s foreign key resolvable. Because the write
  is `if_exists="replace"`, this is also the one step that cannot leave a stale table
  definition behind.
- **Rows are inserted with `pandas.to_sql`** against the already-declared table
  (`if_exists="append"`). The one exception is `run_meta`, which is an upsert from a
  dict — `to_sql` can only append or replace a whole table, and replacing would drop
  the declared schema and wipe other steps' keys.
- **Build pragmas:** `journal_mode=OFF`, `synchronous=OFF` during the load, since a
  crashed build is re-run from the VCFs. Nothing needs restoring afterwards — those
  modes are properties of the build *connection*, not of the file, since only WAL is
  recorded in the database header. The finished database is already an ordinary file
  that any reader opens under its own default, including on `/projects`.
- **Reads are concurrent, writes are not.** Many threads may read; each opens its own
  connection through `schema.connect_read_only`, which uses a `mode=ro` URI and sets
  `busy_timeout`. Never share a connection across threads.
- **Each step indexes its own tables**, as soon as its rows are in and never before
  — an index maintained during a multi-million-row insert costs far more than one
  built once over the finished table. `ANALYZE` is per-table for a sharper reason:
  `sqlite_stat1` is a snapshot that nothing maintains, so a bare `ANALYZE` after
  `annotate` would record *"`call` has 0 rows"*, `load_calls` would then insert
  37.8M, and nothing would correct it. The planner trusts that record
  unconditionally, which makes stale statistics worse than none — with no `stat1`
  row it falls back to a 1M-row guess, at least the right order of magnitude. This
  matters because §7.2 of the plan chose a *view*: every analysis query expands into
  a 5-way join the planner has to order correctly, every time.

## 4. Resolving a locus to one gene

A position can fall inside more than one gene's span. `site.ref_id` records the
winner; the alternatives stay discoverable with a range query, so nothing is stored
twice:

```sql
SELECT ref_id, gene_name, biotype, strand, gene_start, gene_end
FROM gene WHERE chrom = ? AND gene_start <= ? AND gene_end >= ?;
```

Resolution order:

1. **not novel** — a reference gene (has `ref_gene_id`) beats an assembly-only locus
2. **biotype** — `protein_coding` → `lncRNA` → other functional RNA → other / pseudogene
3. **most transcripts** — counted in the assembly, not the reference
4. **largest span** — deliberate and conservative
5. **lowest accession** — determinism only, so re-runs agree

A transcript's gene comes from `ref_gene_id` when it has one. Without one:

| its locus has | gene |
|---|---|
| exactly one reference gene | inherits it — the novel-isoform case |
| no reference gene | the locus itself, `novel = 1` |
| two or more | `NULL` — chimeric locus, no unambiguous parent |

All four counts are reported and recorded in `run_meta`.

## 5. Running it

```bash
module load anaconda3
source activate tarva
cd $ROOT/TARVaCleanupTest/NewTARVaCode
```

The smoke test, which loads every contig and leaves `biotype` NULL:

```bash
python -m tarva.cli dbpipe annotate --gtf $ROOT/FirstPassAssembly/merged.gtf \
  --db $ROOT/TARVa.db
```

The mouse run proper. **The database is `TARVaMouse.db`, not `TARVa.db`** — the
latter is the published human database behind the manuscript and must not be written
to:

```bash
python -m tarva.cli dbpipe annotate \
  --gtf $ROOT/FirstPassAssembly/merged.gtf --db $ROOT/TARVaMouse.db \
  --attrs $ROOT/TestData/GRCm39/MM_KnownAttrs_UCSCBrowser_261011.tsv \
  --attrs-gene-name-column transcriptName \
  --assembly-report $ROOT/TestData/GRCm39/GCF_000001635.27_GRCm39_assembly_report.txt \
  --overwrite
```

Then verify what it loaded. Read-only, four hard checks, exit status 0 only when all
of them pass:

```bash
python sanity_checks/verify_annotate.py --db $ROOT/TARVaMouse.db
```

or `sbatch sanity_checks/verify_annotate.slurm`.

Then the samples. The column defaults match what `tarva.metadata` requires, so the
pipeline's own metadata file needs nothing but the paths:

```bash
python -m tarva.cli dbpipe load-samples \
  --metadata $ROOT/metadata.csv --db $ROOT/TARVaMouse.db
```

Then the calls. The cohort comes from the `individual` table and the contig restriction
from `run_meta`, so neither is named again:

```bash
python -m tarva.cli dbpipe load-calls \
  --rna-unique $ROOT/RNA_Unique --db $ROOT/TARVaMouse.db
```

Then expression, from StringTie round 3's Ballgown tables:

```bash
python -m tarva.cli dbpipe load-expression \
  --stringtie-dir $ROOT/FirstPassAssembly --db $ROOT/TARVaMouse.db
```

Finally the tests, which read only the tables above — and so live in
`analysispipe`, not here:

```bash
python -m tarva.cli analysispipe gene-tests --db $ROOT/TARVaMouse.db
```

### The entry criterion

A `(gene, a1type)` pair is tested only if **at least one condition has every one of
its samples carrying an eligible call** (`dp >= --min-depth`) for that gene and that
type. Pairs where neither condition is complete are dropped before anything else.

That is what makes the zero-padding coherent: the complete condition's arm has no
padding by construction, so the padding lives entirely in the other arm and the test
asks whether a gene is modified across one whole group and largely absent from the
other. `gene_test.complete_in` records which condition qualified.

The criterion is strict — it is what reduces tens of thousands of genes to a few
hundred tested pairs — and it is **asymmetric when the groups differ in size**, since
a smaller group is easier to complete. For the mouse cohort (3 male, 2 female) expect
female-complete pairs to dominate.

`--attrs`, `--transcriptome-fasta` and `--assembly-report` are all optional. Without
them `biotype`, `tsl` and `ccds` stay NULL, the transcript-length check is skipped,
and every contig loads; the three tables still load completely either way.

`--overwrite` **drops and recreates** `gene`, `transcript` and `exon` rather than
deleting their rows, so it survives a schema change — `create_schema` uses
`CREATE TABLE IF NOT EXISTS`, and a `DELETE` would leave an existing database with a
table missing whatever column was just added. It takes the indexes with it and
rebuilds them. To start over completely, delete the `.db` file; there is nothing
else to clean up.

Other flags: `--keep-versions`, `--chunk-size`, `--min-attrs-match`,
`--drop-molecules`, `--overwrite`, and
`--attrs-{transcript,gene,gene-name,gene-type,transcript-type,tsl,ccds}-column` for
an attribute table whose headers differ from the UCSC `knownAttrs` defaults.

### Mouse versus human identifier spaces

The mouse test pipeline is **NCBI RefSeq** end to end
(`GCF_000001635.27_GRCm39_genomic.*`), so `merged.gtf` carries RefSeq identifiers —
`NC_000067.7` contigs, `ref_gene_id "gene-Xkr4"`, `transcript_id
"rna-XM_006495550.5"`. Note what the reference GFF3 keys its genes on: after the
namespace strip, `gene.ref_id` holds a bare **gene symbol** (`Xkr4`) for the mouse
run where the human run holds an ENSG accession. Same column, different kind of
identifier per annotation; `run_meta` records which.

`MM_KnownAttrs_UCSCBrowser_261011.tsv` is Ensembl-keyed (`ENSMUST`/`ENSMUSG`), so it
cannot be joined on accession. It can be joined on **symbol**, because the crosswalk
is already inside the file: `transcriptName` is `<symbol>-<number>` (`Gnai3-201`), and
stripping the suffix gives the same MGI symbols `merged.gtf` carries in `gene_name`.
That is what `--attrs-gene-name-column transcriptName` does. No external crosswalk is
needed, and none would help — there is no string rule relating an `ENSMUST` to its
`ENSMUSG`, the matching low numbers in the oldest entries being a historical artefact
of parallel id minting rather than a rule.

Two consequences. The join is **gene-level only**: `transcript_type`, `tsl` and `ccds`
stay NULL for mouse, because `kgID` is an `ENSMUST` and there is no `ENSMUST`
counterpart to a Gnomon-predicted `XM_`. And the miss mode is **symbol synonymy**,
silent by nature — the same disagreement as human's `DDX11L2` versus `DDX11L16`. The
reported match rate is the thing to read.

This is enough for the resolution rule, which needs `gene.biotype` and nothing else,
and it means the mouse biotype labels come from the same vocabulary as the published
human ones. **So the biotype step of the resolution rule is exercised by the mouse
run** — an earlier version of this section said it would not be.

`gencode.vM39.transcripts.fa.gz` still cannot serve as `--transcriptome-fasta`: its
records are `ENSMUST`-keyed and its exon arithmetic is GENCODE's, so a mismatch would
mean the two annotations differ, not that this build is wrong. Build a RefSeq-named
FASTA from the GFF already in place if the length check is wanted:

```bash
module load anaconda3
source activate stringtie
gffread -w GRCm39_refseq_transcripts.fa -g GCF_000001635.27_GRCm39_genomic.fna \
  GCF_000001635.27_GRCm39_genomic.gff
```

Check its first header before trusting the result. `check_fasta` runs record names
through the same normalization, so `>rna-XM_006495550.5` and `>XM_006495550.5` both
work — but a `0 of 147,908 present` count is what a failed normalization looks like,
and it is reported rather than silent.

The human data needs none of this: its Gencode-based assembly and
`knownAttrs_UCSCBrowser_260906.tsv` are the same identifier space, so the default
accession join applies.

### Contig selection

`--assembly-report` is the only way to express "primary assembly only". The
accession prefix cannot: `NC_` means *assembled molecule*, which covers the
mitochondrion, so keeping `NC_` drops the scaffolds and keeps chrM. The report's
`Sequence-Role` separates `assembled-molecule` from `unplaced-scaffold` and
`unlocalized-scaffold`, and `Assigned-Molecule` names `1`–`19`, `X`, `Y`, `MT`.

`--drop-molecules` defaults to `MT,Y`. **chrY was dropped in the published human
analysis**, so this reproduces that rather than deviating from it. chrX is kept and
carries no sex confound: X-inactivation equalizes X-linked expression, and `a1prop`
is a depth-normalized proportion measured on RNA *after* `rna-unique` has subtracted
the WGS calls — so germline heterozygosity, the only thing male hemizygosity would
alter, is already gone. Whether chrM and the scaffolds were also excluded from the
published set is **unconfirmed**, and both change the gene set the reproduction gate
runs over.

The report is indexed under every name it gives a sequence — RefSeq accession,
GenBank accession, sequence name and UCSC-style name — so one file classifies this
RefSeq-keyed assembly and a GENCODE-keyed `chr1` one alike.

## 6. Checks the build runs

| Check | Behaviour on failure |
|---|---|
| `Σ(exon lengths) == trans_len` per transcript | counted, written to `annotate_fasta_mismatches.csv`, warned — not fatal |
| attribute-table match rate | fatal below `--min-attrs-match`, measured on whichever join was attempted |
| transcript or exon record with no `transcript_id` | counted, warned, dropped |
| transcript record with no `gene_id` (so no locus) | counted, warned, dropped; reference genes lost entirely are named |
| contig not an assembled molecule, or in `--drop-molecules` | counted per reason, dropped |
| no contig surviving contig selection at all | fatal — names both identifier spaces |
| no transcript record carrying a `transcript_id` at all | fatal — the file yielded nothing keyable |
| identifier still carrying a GFF3 namespace (`cds-`, `exon-`, a doubled `gene-`) after stripping | warned, left as-is |
| transcripts with no exon records | warned, dropped |
| exons belonging to no transcript | warned, dropped |
| genes with transcripts on several contigs or both strands | warned, first kept |
| duplicate transcript or exon records | warned, collapsed |

The two fatal cases are both "this file yielded nothing usable", not "a record was
malformed". Individual unkeyable records are dropped by design.

## 7. Open items

1. ~~**Is `exon_number` transcript order or coordinate order?**~~ **Settled —
   coordinate order.** `sanity_checks/check_gtf_conventions.py` over
   `merged.gtf`: of 47,535 multi-exon transcripts, **23,699 minus-strand ascending,
   0 descending, 0 mixed**, and 23,836 plus-strand ascending. `exon_number` is
   coordinate order on both strands, so the strand-aware recomputation in
   `build_exons` is correct and load-bearing — trusting the attribute would have
   mirrored every one of those 23,699.
2. ~~**Do `strand`, `ref_id` and `exonic` belong on `site` or on `call`?**~~ **Void.**
   The question only existed under the rule that strand follows the VCF substitution,
   which made strand a function of the per-call `alt1`. That rule is **wrong**: it
   assumes a call is an edit in order to decide its orientation, and the published
   `a1type` counts disprove it, carrying **5,717 `T-C` and 6,611 `G-A`** rows as
   categories of their own, which substitution-inferred strand could not have
   produced. Strand is the covering gene's, so all three are functions of position
   and stay on `site`.
3. ~~**Strand for substitutions that do not imply one.**~~ **Void** — nothing implies
   strand, and nothing needs to. Every variant class takes the covering gene's strand:
   substitutions that look like edits, substitutions that do not, and indels alike. No
   fallback, no special case.
4. ~~**Is the `exonic` flag wanted?**~~ **Dropped.** Not needed to run the pipeline as
   legacy intended, so it is not declared rather than declared and left NULL. If the
   reproduction gate misses, this is the first thing to re-examine: legacy's effective
   set may have been exon-only, while `site.ref_id` is assigned by gene **span**,
   introns included — which is what the published `a1prop` denominator uses.
5. ~~**Does `expression` stay in core?**~~ **Yes, settled** — `load-expression` is
   written. Both consumers, TPM-correlation and deconvolution, are in `Other` and now
   read the table instead of re-walking the ctab files, and the `build_crosswalk` path
   mismatch is fixed once rather than twice. Keyed on `tid`, not `ref_id`: the ctab is
   per transcript, and gene-level expression is a SUM over `transcript` rather than a
   second stored copy. `varcount` stays derivable from `call`.
6. **Reproduction gate.** The published database should reproduce
   **582 / 24 / 72 / 17 / 33**. Whether the rebuild can is open: the published results
   rest on the 860,394 of 37,797,345 `sample_tab` rows that carried annotation, and
   how much of that 97.7% shortfall is `genomic='yes'` plus sub-depth plus non-variant
   rows — rather than anything else — has not been measured.

## 8. Not in this pipeline

HGVS and consequence coordinates, the external A-I database checks (REDIportal, AD
editome), the TPM-correlation analysis and the DSA deconvolution all live in `Other`.
They read the database; none of them writes to it.
