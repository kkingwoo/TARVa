"""
dbpipe annotate: load gene, transcript and exon from the merged GTF

    `tarva dbpipe annotate --gtf <merged.gtf> --db <TARVa.db>
        [--attrs <knownAttrs.tsv>] [--transcriptome-fasta <transcripts.fa>]
        [--assembly-report <GCF_..._assembly_report.txt>]
        [--drop-molecules MT,Y] [--keep-versions] [--chunk-size 500000]
        [--attrs-transcript-column kgID] [--attrs-gene-column geneId]
        [--attrs-gene-name-column transcriptName]
        [--attrs-gene-type-column geneType]
        [--attrs-transcript-type-column transcriptType]
        [--attrs-tsl-column supportLevel] [--attrs-ccds-column ccdsId]
        [--min-attrs-match 0.5] [--overwrite]`

The annotation layer, and the first step of the build. Everything else in
dbpipe depends on these three tables.

A locus with no reference gene at all becomes ONE gene row keyed on the locus
id with `novel = 1`. A novel transcript inside a locus carrying exactly one
reference gene inherits it -- that is the novel-isoform-of-a-known-gene case,
and the common one. 

Records this build cannot key are dropped and counted, not fatal: a record
with no transcript_id cannot be joined or coordinate-mapped, and one with no
gene_id has not been placed at a locus. An annotation is allowed to contain
records this build has no use for. The `annotate.n_dropped_*` keys are the
only record of how many, so they are worth reading after a run.

Writes: gene, transcript, exon. Records annotate.* keys in run_meta.
"""

import csv
import gzip
import re

from pathlib import Path

import numpy as np
import pandas as pd

from tarva.dbpipe import schema

GTF_COLUMNS = ["chrom", "source", "feature", "start", "end",
               "score", "strand", "frame", "attrs"]

VERSIONED = re.compile(r"^([A-Za-z][A-Za-z0-9_]*[0-9])\.([0-9]+)$")

GFF3_ID_PREFIX = re.compile(r"^(?:gene|rna|id)-")

GFF3_NAMESPACES = re.compile(r"^(?:gene|rna|id|cds|exon)-")

TRANSCRIPT_NAME_SUFFIX = re.compile(r"-[0-9]+$")

# The assembly report's column names are IN the file: NCBI ends its comment
# block with a tab-separated header line, "# Sequence-Name<TAB>Sequence-Role
# <TAB>...". So they are read from there and matched BY NAME.
ASSEMBLY_ROLE_COLUMN = "sequence_role"
ASSEMBLY_MOLECULE_COLUMN = "assigned_molecule"
ASSEMBLY_NAME_COLUMNS = ("refseq_accn", "genbank_accn",
                         "sequence_name", "ucsc_style_name")

ASSEMBLED_MOLECULE = "assembled-molecule"

# "na" is how the report marks a name a sequence does not have.
ASSEMBLY_ABSENT = "na"

DEFAULT_CHUNK_SIZE = 500000
DEFAULT_MIN_ATTRS_MATCH = 0.5
DEFAULT_DROP_MOLECULES = "MT,Y"

FASTA_MISMATCH_CSV = "annotate_fasta_mismatches.csv"


def add_arguments(parser):
    parser.add_argument("--gtf", required=True,
                        help="Merged StringTie GTF (merged.gtf from `stringtiepipe assemble`). "
                             "May be gzipped")
    parser.add_argument("--db", required=True,
                        help="SQLite database to create or add to")
    parser.add_argument("--attrs", metavar="TSV",
                        help="Annotation attribute table (UCSC knownAttrs or equivalent), "
                             "supplying biotype, TSL and CCDS. Joined on the UNVERSIONED "
                             "accession, or on the gene symbol with "
                             "--attrs-gene-name-column. ")
    parser.add_argument("--transcriptome-fasta", metavar="FASTA",
                        help="Transcript FASTA for transcript length and site coordinate check. ")
    parser.add_argument("--assembly-report", metavar="TXT",
                        help="NCBI *_assembly_report.txt for this reference. Drops unplaced and unlocalized "
                             "scaffolds. ")
    parser.add_argument("--drop-molecules", default=DEFAULT_DROP_MOLECULES,
                        help=f"Assembled molecules to exclude as well, by name, "
                             f"case-insensitive (default: {DEFAULT_DROP_MOLECULES}, for example MT or Y genes). Requires "
                             f"--assembly-report")
    parser.add_argument("--keep-versions", action="store_true",
                        help="Keep accession versions in ref_id and tid. Check the reported match rate. ")
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE,
                        help=f"GTF lines per chunk (default: {DEFAULT_CHUNK_SIZE})")
    parser.add_argument("--min-attrs-match", type=float, default=DEFAULT_MIN_ATTRS_MATCH,
                        help=f"Fail if fewer than this fraction of transcripts match --attrs "
                             f"(default: {DEFAULT_MIN_ATTRS_MATCH}). ")
    parser.add_argument("--attrs-transcript-column", default="kgID",
                        help="--attrs column holding the transcript accession (default: kgID)")
    parser.add_argument("--attrs-gene-column", default="geneId",
                        help="--attrs column holding the gene accession (default: geneId)")
    parser.add_argument("--attrs-gene-name-column", default=None,
                        help="--attrs column holding the gene SYMBOL, which joins on symbol "
                             "instead of on the gene accession. For a UCSC knownAttrs table "
                             "that is transcriptName, whose values are <symbol>-<number>; the "
                             "suffix is stripped. Use it when the annotation and the attribute "
                             "table are in different accession spaces but share symbols -- an "
                             "Ensembl-keyed table against a RefSeq-keyed assembly. Gene "
                             "biotype only: transcript_type, tsl and ccds stay NULL. ")
    parser.add_argument("--attrs-gene-type-column", default="geneType",
                        help="--attrs column holding the gene biotype (default: geneType)")
    parser.add_argument("--attrs-transcript-type-column", default="transcriptType",
                        help="--attrs column holding the transcript biotype "
                             "(default: transcriptType)")
    parser.add_argument("--attrs-tsl-column", default="supportLevel",
                        help="--attrs column holding transcript support level "
                             "(default: supportLevel)")
    parser.add_argument("--attrs-ccds-column", default="ccdsId",
                        help="--attrs column holding CCDS membership (default: ccdsId)")
    parser.add_argument("--overwrite", action="store_true",
                        help="Replace existing gene/transcript/exon rows instead of failing")


#----------------------------------------------------
# identifiers
#----------------------------------------------------

def strip_id_prefix(series, label=None):
    """
    Drop the GFF3 ID namespace from an accession series.

    """
    stripped = series.str.replace(GFF3_ID_PREFIX, "", regex=True)
    if label is not None:
        # `series.notna()` guard: in pandas NaN != NaN is True, so a
        # bare comparison counts every absent identifier as one that changed.
        changed = int((stripped.ne(series) & series.notna()).sum())
        if changed:
            schema.report(f"  {label}: GFF3 ID prefix stripped from {changed:,} of "
                          f"{len(series):,} identifier(s)")
        odd = stripped[stripped.str.contains(GFF3_NAMESPACES, na=False)].unique()
        if len(odd):
            schema.warn(f"{label}: {len(odd):,} distinct identifier(s) still carry a "
                        f"GFF3 namespace after stripping "
                        f"(e.g. {', '.join(map(str, odd[:3]))}) -- left as-is, so any "
                        f"join on them will miss")
    return stripped


def split_versions(series, keep=False, label=None):
    """
    (stem, version) for an accession series. See VERSIONED for the rule.

    """
    series = strip_id_prefix(series, label=label)
    if keep:
        return series, pd.Series([None] * len(series), index=series.index, dtype=object)
    extracted = series.str.extract(VERSIONED)
    stem = extracted[0].where(extracted[0].notna(), series)
    return stem, extracted[1]


#----------------------------------------------------
# the GTF
#----------------------------------------------------

def _open_text(path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def parse_gtf(path, chunk_size):
    """
    Read transcript and exon records, extracting attributes by key.

    """
    path = Path(path)
    if not path.is_file():
        raise SystemExit(f"tarva dbpipe annotate: --gtf not found: {path}")

    keys = ("gene_id", "transcript_id", "ref_gene_id", "gene_name")
    transcript_frames, exon_frames = [], []
    lines_read = 0

    with _open_text(path) as handle:
        reader = pd.read_csv(handle, sep="\t", header=None, names=GTF_COLUMNS,
                             dtype=str, chunksize=chunk_size, na_filter=False,
                             engine="c")
        for number, chunk in enumerate(reader, start=1):
            lines_read += len(chunk)
            chunk = chunk[chunk["feature"].isin(("transcript", "exon"))].copy()
            if chunk.empty:
                continue
            for key in keys:
                chunk[key] = chunk["attrs"].str.extract(
                        rf'\b{key} "([^"]*)"', expand=False)
                chunk[key] = chunk[key].replace("", np.nan)
            chunk["start"] = pd.to_numeric(chunk["start"], errors="coerce")
            chunk["end"] = pd.to_numeric(chunk["end"], errors="coerce")

            is_transcript = chunk["feature"].eq("transcript")
            transcript_frames.append(chunk.loc[is_transcript, [
                "chrom", "strand", "start", "end",
                "gene_id", "transcript_id", "ref_gene_id", "gene_name"]])
            exon_frames.append(chunk.loc[~is_transcript, [
                "chrom", "strand", "start", "end", "transcript_id"]])
            schema.report(f"  GTF: chunk {number} -- {lines_read:,} lines read")

    if not transcript_frames:
        raise SystemExit(
                f"tarva dbpipe annotate: no transcript records in {path}. Is this a GTF?")

    transcripts = pd.concat(transcript_frames, ignore_index=True)
    exons = pd.concat(exon_frames, ignore_index=True) if exon_frames else pd.DataFrame(
            columns=["chrom", "strand", "start", "end", "transcript_id"])

    # A record with no transcript_id cannot be keyed, joined or coordinate-
    # mapped, so it is dropped.
    no_tid = transcripts["transcript_id"].isna()
    n_transcripts_no_tid = int(no_tid.sum())
    if n_transcripts_no_tid:
        schema.warn(f"{n_transcripts_no_tid:,} transcript record(s) have no "
                    f"transcript_id -- dropped")
        transcripts = transcripts.loc[~no_tid].copy()

    no_tid_exon = exons["transcript_id"].isna()
    n_exons_no_tid = int(no_tid_exon.sum())
    if n_exons_no_tid:
        schema.warn(f"{n_exons_no_tid:,} exon record(s) have no transcript_id -- dropped")
        exons = exons.loc[~no_tid_exon].copy()

    if transcripts.empty:
        raise SystemExit(
                f"tarva dbpipe annotate: no transcript record in {path} carries a "
                f"transcript_id. Attributes are read by key, so either this is not a GTF "
                f"or its transcript records are keyed differently.")

    schema.report(f"GTF: {len(transcripts):,} transcript and {len(exons):,} exon record(s) "
                  f"from {lines_read:,} lines")
    return transcripts, exons, {
        "transcripts_no_tid": n_transcripts_no_tid,
        "exons_no_tid": n_exons_no_tid,
        }


def build_exons(exons, keep_versions):
    """
    Exon table with transcript-order numbering and cumulative offsets.

    exon_num is derived from coordinate rank and strand, NOT from the GTF's
    exon_number attribute, which is coordinate order on both strands.
    Coordinates stay 1-based inclusive as the GTF gives them, so an exon's
    length is (exon_end - exon_start + 1).
    """
    exons = exons.dropna(subset=["start", "end", "transcript_id"]).copy()
    exons["tid"], _ = split_versions(exons["transcript_id"], keep=keep_versions)
    exons = exons.rename(columns={"start": "exon_start", "end": "exon_end"})
    exons["exon_start"] = exons["exon_start"].astype("int64")
    exons["exon_end"] = exons["exon_end"].astype("int64")

    before = len(exons)
    exons = exons.drop_duplicates(subset=["tid", "exon_start", "exon_end"])
    if len(exons) != before:
        schema.warn(f"{before - len(exons):,} duplicate exon record(s) collapsed")

    # mergesort: stable, so equal coordinates keep file order rather than
    # shuffling between runs.
    exons = exons.sort_values(["tid", "exon_start"], kind="mergesort")
    exons["length"] = exons["exon_end"] - exons["exon_start"] + 1

    grouped = exons.groupby("tid", sort=False)
    coord_rank = grouped.cumcount()
    n_exons = grouped["exon_start"].transform("size")
    # '+': lowest coordinate is exon 1. '-': highest coordinate is exon 1.
    exons["exon_num"] = np.where(
            exons["strand"].values == "-",
            n_exons.values - coord_rank.values,
            coord_rank.values + 1,
            )

    exons = exons.sort_values(["tid", "exon_num"], kind="mergesort")
    running = exons.groupby("tid", sort=False)["length"].cumsum()
    exons["trans_offset"] = running - exons["length"]

    trans_len = exons.groupby("tid", sort=False)["length"].sum().rename("trans_len")
    return exons, trans_len


def resolve_genes(transcripts):
    """
    Assign each transcript its gene, and report the four cases.

    A transcript with ref_gene_id keeps it. Without one:
      - its locus has exactly one reference gene  -> inherit it
      - its locus has no reference gene           -> the locus IS the gene
      - its locus has two or more                 -> NULL, no unambiguous parent
    """
    has_ref = transcripts["ref_id"].notna()
    with_ref = transcripts.loc[has_ref]
    refs_per_locus = with_ref.groupby("locus_id")["ref_id"].nunique()
    first_ref = (with_ref.drop_duplicates("locus_id")
                 .set_index("locus_id")["ref_id"])

    locus = transcripts["locus_id"]
    n_refs = locus.map(refs_per_locus).fillna(0).astype("int64")

    inherit = (~has_ref) & n_refs.eq(1)
    own_gene = (~has_ref) & n_refs.eq(0)
    ambiguous = (~has_ref) & n_refs.ge(2)

    resolved = transcripts["ref_id"].copy()
    resolved[inherit] = locus[inherit].map(first_ref)
    resolved[own_gene] = locus[own_gene]
    transcripts["ref_id"] = resolved

    counts = {
        "matched": int(has_ref.sum()),
        "inherited": int(inherit.sum()),
        "novel_loci_transcripts": int(own_gene.sum()),
        "ambiguous": int(ambiguous.sum()),
        }
    novel_gene_ids = set(locus[own_gene].unique())
    return transcripts, counts, novel_gene_ids


def build_gene(transcripts, novel_gene_ids):
    """
    One row per gene. gene_span is the contiguous span, introns included.

    """
    assigned = transcripts.dropna(subset=["ref_id"])
    grouped = assigned.groupby("ref_id", sort=True)

    gene = pd.DataFrame({
        "chrom": grouped["chrom"].first(),
        "strand": grouped["strand"].first(),
        "gene_start": grouped["trans_start"].min().astype("int64"),
        "gene_end": grouped["trans_end"].max().astype("int64"),
        "n_transcripts": grouped.size().astype("int64"),
        "ref_id_version": grouped["ref_id_version"].first(),
        "gene_name": grouped["gene_name"].first(),
        }).reset_index()
    gene["gene_span"] = gene["gene_end"] - gene["gene_start"] + 1
    gene["novel"] = gene["ref_id"].isin(novel_gene_ids).astype("int64")

    mixed_chrom = grouped["chrom"].nunique()
    scattered = mixed_chrom[mixed_chrom > 1]
    if len(scattered):
        schema.warn(f"{len(scattered):,} gene(s) have transcripts on more than one contig "
                    f"(e.g. {', '.join(map(str, scattered.index[:3]))}) -- keeping the first")
    mixed_strand = grouped["strand"].nunique()
    flipped = mixed_strand[mixed_strand > 1]
    if len(flipped):
        schema.warn(f"{len(flipped):,} gene(s) have transcripts on both strands "
                    f"(e.g. {', '.join(map(str, flipped.index[:3]))}) -- keeping the first")
    return gene


#----------------------------------------------------
# contig selection
#----------------------------------------------------

def read_report_columns(path):
    """
    The column names from the last comment line of an assembly report.

    """
    header = None
    with _open_text(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            header = line
    if header is None:
        raise SystemExit(
                f"tarva dbpipe annotate: {path} has no comment block, so there is no "
                f"header line to read column names from. Is this an NCBI "
                f"*_assembly_report.txt?")
    names = header.lstrip("#").rstrip("\n").split("\t")
    return [re.sub(r"[^a-z0-9]+", "_", name.strip().lower()).strip("_")
            for name in names]


def load_assembly_report(path):
    """
    {contig name: (role, molecule)} for every sequence in an assembly.

    Columns come from the file's own header line, matched by name.
    """
    path = Path(path)
    if not path.is_file():
        raise SystemExit(f"tarva dbpipe annotate: --assembly-report not found: {path}")

    columns = read_report_columns(path)
    wanted = (ASSEMBLY_ROLE_COLUMN, ASSEMBLY_MOLECULE_COLUMN) + ASSEMBLY_NAME_COLUMNS
    missing = [name for name in wanted if name not in columns]
    if missing:
        raise SystemExit(
                f"tarva dbpipe annotate: {path} has no column(s) {', '.join(missing)}. "
                f"The header line gives: {', '.join(columns)}. Is this an NCBI "
                f"*_assembly_report.txt?")

    table = pd.read_csv(path, sep="\t", comment="#", header=None,
                        names=columns, dtype=str, na_filter=False)
    roles = table[ASSEMBLY_ROLE_COLUMN]
    molecules = table[ASSEMBLY_MOLECULE_COLUMN]
    classified = {}
    for column in ASSEMBLY_NAME_COLUMNS:
        for name, role, molecule in zip(table[column], roles, molecules):
            if name and name != ASSEMBLY_ABSENT:
                classified[name] = (role, molecule)
    schema.report(f"--assembly-report: {len(table):,} sequence(s), "
                  f"{len(classified):,} name(s) from {path.name}")
    return classified


def filter_contigs(transcripts, exons, args):
    """
    Restrict the annotation to the primary assembly, minus --drop-molecules.

    """
    report_map = load_assembly_report(args.assembly_report)
    drop = {m.strip().upper() for m in args.drop_molecules.split(",") if m.strip()}

    present = sorted(set(transcripts["chrom"].dropna()))
    kept, reasons = [], {}
    for contig in present:
        entry = report_map.get(contig)
        if entry is None:
            reasons[contig] = "not named in the assembly report"
        elif entry[0] != ASSEMBLED_MOLECULE:
            reasons[contig] = entry[0]
        elif entry[1].upper() in drop:
            reasons[contig] = f"molecule {entry[1]}"
        else:
            kept.append(contig)

    if not kept:
        raise SystemExit(
                f"tarva dbpipe annotate: no contig in the annotation survived "
                f"--assembly-report selection. The GTF names contigs "
                f"{', '.join(present[:3])}; the report names them "
                f"{', '.join(list(report_map)[:3])}. Those are different identifier "
                f"spaces -- check that the report matches the reference the assembly was "
                f"built against.")

    keep = set(kept)
    on_dropped = ~transcripts["chrom"].isin(keep)
    n_dropped = int(on_dropped.sum())

    by_reason = {}
    for contig, reason in reasons.items():
        by_reason[reason] = by_reason.get(reason, 0) + int(
                (transcripts["chrom"] == contig).sum())
    for reason, n in sorted(by_reason.items(), key=lambda pair: -pair[1]):
        schema.report(f"  contigs: {n:,} transcript(s) dropped -- {reason}")

    if n_dropped:
        transcripts = transcripts.loc[~on_dropped].copy()
        exons = exons.loc[exons["chrom"].isin(keep)].copy()
    schema.report(f"contigs: kept {len(kept):,} of {len(present):,} present "
                  f"({', '.join(kept[:5])}{', ...' if len(kept) > 5 else ''})")

    return transcripts, exons, kept, n_dropped


#----------------------------------------------------
# --attrs
#----------------------------------------------------

def load_attrs(path, args):
    """
    Read the attribute table and key it on unversioned accessions.
    """
    path = Path(path)
    if not path.is_file():
        raise SystemExit(f"tarva dbpipe annotate: --attrs not found: {path}")

    wanted = {
        "tid_raw": args.attrs_transcript_column,
        "gene_raw": args.attrs_gene_column,
        "gene_type": args.attrs_gene_type_column,
        "transcript_type": args.attrs_transcript_type_column,
        "tsl": args.attrs_tsl_column,
        "ccds": args.attrs_ccds_column,
        }

    table = pd.read_csv(path, sep="\t", dtype=str, na_filter=False)
    table.columns = [str(c).lstrip("#").strip() for c in table.columns]

    absent = {name: column for name, column in wanted.items()
              if column not in table.columns}
    required = ("tid_raw", "gene_raw")
    fatal = [wanted[name] for name in required if name in absent]
    if fatal:
        raise SystemExit(
                f"tarva dbpipe annotate: {path} has no column(s) {', '.join(fatal)}. "
                f"Columns found: {', '.join(table.columns)}. Set the --attrs-*-column flags "
                f"if this table names them differently.")
    for name, column in absent.items():
        if name not in required:
            schema.warn(f"{path.name} has no column {column!r} -- {name} will be NULL")

    frame = pd.DataFrame(index=table.index)
    for name, column in wanted.items():
        frame[name] = table[column] if column in table.columns else None

    frame["tid"], _ = split_versions(frame["tid_raw"].astype(str), keep=False)
    frame["ref_id"], _ = split_versions(frame["gene_raw"].astype(str), keep=False)
    frame["tsl"] = pd.to_numeric(frame["tsl"], errors="coerce")

    if args.attrs_gene_name_column:
        if args.attrs_gene_name_column not in table.columns:
            raise SystemExit(
                    f"tarva dbpipe annotate: {path} has no column "
                    f"{args.attrs_gene_name_column}. Columns found: "
                    f"{', '.join(table.columns)}")
        frame["symbol"] = (table[args.attrs_gene_name_column]
                           .str.replace(TRANSCRIPT_NAME_SUFFIX, "", regex=True))
        schema.report(f"--attrs: {frame['symbol'].nunique():,} distinct gene symbol(s) "
                      f"from {args.attrs_gene_name_column}")

    frame = frame.replace({"": None})
    schema.report(f"--attrs: {len(frame):,} row(s), "
                  f"{frame['ref_id'].nunique():,} gene(s) from {path.name}")
    return frame


def apply_attrs(gene, transcript, attrs, min_match, by_symbol=False):
    """
    Merge biotype/TSL/CCDS on, and refuse a join that does not allow matching.

    """
    if by_symbol:
        per_gene = (attrs.dropna(subset=["symbol"])
                    .drop_duplicates("symbol").set_index("symbol"))
        gene["biotype"] = gene["gene_name"].map(per_gene["gene_type"])
        matched, total, unit = int(gene["biotype"].notna().sum()), len(gene), "gene"
        schema.warn("--attrs joined on gene symbol: transcript_type, tsl and ccds stay "
                    "NULL, as the two transcript model sets do not correspond")
    else:
        per_transcript = attrs.drop_duplicates("tid").set_index("tid")
        transcript["transcript_type"] = transcript["tid"].map(per_transcript["transcript_type"])
        transcript["tsl"] = transcript["tid"].map(per_transcript["tsl"])
        transcript["ccds"] = transcript["tid"].map(per_transcript["ccds"])

        per_gene = attrs.dropna(subset=["ref_id"]).drop_duplicates("ref_id").set_index("ref_id")
        gene["biotype"] = gene["ref_id"].map(per_gene["gene_type"])
        matched, total, unit = (int(transcript["transcript_type"].notna().sum()),
                                len(transcript), "transcript")

    rate = matched / total if total else 0.0
    schema.report(f"--attrs: matched {matched:,} of {total:,} {unit}(s) ({rate:.1%})")
    if rate < min_match:
        raise SystemExit(
                f"tarva dbpipe annotate: only {rate:.1%} of {unit}s matched --attrs, "
                f"below --min-attrs-match {min_match:.0%}. The two files are probably in "
                f"different identifier spaces, sucha as an Ensembl attribute table against an NCBI "
                f"RefSeq annotation. Pass --attrs-gene-name-column to join on "
                f"the gene symbol instead. ")

    known = gene["biotype"].notna()
    schema.report(f"--attrs: {int(known.sum()):,} of {len(gene):,} gene(s) have a biotype")
    return gene, transcript


#----------------------------------------------------
# the integrity check
#----------------------------------------------------

def fasta_lengths(path):
    """
    {accession: sequence length} from a transcript FASTA, streamed.

    """
    lengths = {}
    name, total = None, 0
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = total
                name = line[1:].strip().split("|")[0].split()[0] if len(line) > 1 else ""
                total = 0
            else:
                total += len(line.strip())
    if name is not None:
        lengths[name] = total
    return lengths


def check_fasta(transcript, path, keep_versions, out_dir):
    """
    Assert sum(exon lengths) == sequence length, per transcript.

    """
    path = Path(path)
    if not path.is_file():
        raise SystemExit(f"tarva dbpipe annotate: --transcriptome-fasta not found: {path}")

    lengths = fasta_lengths(path)
    stems, _ = split_versions(pd.Series(list(lengths.keys()), dtype=object),
                              keep=keep_versions)
    by_stem = dict(zip(stems, lengths.values()))

    observed = transcript["tid"].map(by_stem)
    present = observed.notna()
    disagree = present & (observed != transcript["trans_len"])

    schema.report(f"FASTA check: {int(present.sum()):,} of {len(transcript):,} transcript(s) "
                  f"present in {path.name}; {int(disagree.sum()):,} length mismatch(es)")
    if int(disagree.sum()):
        report_path = Path(out_dir) / FASTA_MISMATCH_CSV
        rows = transcript.loc[disagree, ["tid", "trans_len"]].copy()
        rows["fasta_len"] = observed[disagree].astype("int64")
        rows.to_csv(report_path, index=False)
        schema.warn(f"exon arithmetic disagrees with the FASTA for "
                    f"{int(disagree.sum()):,} transcript(s) -- see {report_path}. Every "
                    f"transcript coordinate for those is suspect.")
    return int(present.sum()), int(disagree.sum())


#----------------------------------------------------
# writing
#----------------------------------------------------

WRITE_CHUNK = 100000

ANNOTATE_TABLES = ("gene", "transcript", "exon")


def write_tables(con, gene, transcript, exon, overwrite):
    gene_columns = ["ref_id", "ref_id_version", "gene_name", "biotype", "novel",
                    "chrom", "strand", "gene_start", "gene_end", "gene_span",
                    "n_transcripts"]
    transcript_columns = ["tid", "tid_version", "ref_id", "locus_id", "chrom", "strand",
                          "trans_start", "trans_end", "trans_len", "transcript_type",
                          "tsl", "ccds"]
    exon_columns = ["tid", "exon_num", "exon_start", "exon_end", "trans_offset"]

    existing = con.execute("SELECT COUNT(*) FROM gene").fetchone()[0]
    if existing and not overwrite:
        raise SystemExit(
                f"tarva dbpipe annotate: gene already holds {existing:,} row(s). Pass "
                f"--overwrite to replace the annotation, which also invalidates any site, "
                f"call or gene_test rows keyed to it.")
    if existing:
        # DROP and recreate rather than DELETE, so --overwrite survives a
        # schema change.
        for table in reversed(ANNOTATE_TABLES):
            con.execute(f"DROP TABLE IF EXISTS {table}")
        con.commit()
        schema.create_schema(con, ANNOTATE_TABLES)
        schema.report(f"--overwrite: dropped and recreated "
                      f"{', '.join(ANNOTATE_TABLES)}")

    # to_sql against the already-declared table: if_exists="append" inserts
    # into the schema create_schema() made rather than inventing one, and
    # pandas handles the NaN -> NULL and numpy -> sqlite conversions. Default
    # method (executemany under the hood) rather than "multi", which builds one
    # statement per chunk and would exceed SQLite's parameter limit at this
    # chunk size and column count.
    for table, frame, columns in (("gene", gene, gene_columns),
                                  ("transcript", transcript, transcript_columns),
                                  ("exon", exon, exon_columns)):
        frame.reindex(columns=columns).to_sql(
                table, con, if_exists="append", index=False, chunksize=WRITE_CHUNK)
        con.commit()
        schema.report(f"{table}: {len(frame):,} row(s) written")

    # Index table immediately after build
    schema.create_indexes(con, ANNOTATE_TABLES)


#----------------------------------------------------
# entry point
#----------------------------------------------------

def run(args):
    db_path = Path(args.db)
    out_dir = db_path.parent

    transcripts, exons, parse_counts = parse_gtf(args.gtf, args.chunk_size)

    transcripts = transcripts.rename(columns={"start": "trans_start", "end": "trans_end"})
    transcripts = transcripts.dropna(subset=["trans_start", "trans_end"]).copy()
    transcripts["trans_start"] = transcripts["trans_start"].astype("int64")
    transcripts["trans_end"] = transcripts["trans_end"].astype("int64")

    # Contig selection before build_exons, so exon numbering and offsets are
    # never computed for rows about to be dropped.
    n_contigs_present = int(transcripts["chrom"].nunique())
    if args.assembly_report:
        transcripts, exons, contigs_kept, n_contig_dropped = filter_contigs(
                transcripts, exons, args)
    else:
        contigs_kept, n_contig_dropped = None, 0
        schema.warn(f"no --assembly-report: all {n_contigs_present:,} contig(s) in the "
                    f"annotation are loaded, including any unplaced or unlocalized "
                    f"scaffolds, and both chrM and chrY")

    transcripts["tid"], transcripts["tid_version"] = split_versions(
            transcripts["transcript_id"], keep=args.keep_versions, label="transcript_id")
    transcripts["ref_id"], transcripts["ref_id_version"] = split_versions(
            transcripts["ref_gene_id"], keep=args.keep_versions, label="ref_gene_id")

    no_locus = transcripts["gene_id"].isna()
    n_no_locus = int(no_locus.sum())
    lost_genes = set()
    if n_no_locus:
        gone = transcripts.loc[no_locus]
        with_ref = int(gone["ref_id"].notna().sum())
        lost_genes = (set(gone["ref_id"].dropna())
                      - set(transcripts.loc[~no_locus, "ref_id"].dropna()))
        schema.warn(f"{n_no_locus:,} transcript(s) have no gene_id and so no locus "
                    f"-- dropped ({with_ref:,} of them carried a ref_gene_id)")
        if lost_genes:
            schema.warn(f"{len(lost_genes):,} reference gene(s) appear ONLY on those "
                        f"dropped transcripts and will be absent from `gene` "
                        f"(e.g. {', '.join(sorted(lost_genes)[:3])})")
        transcripts = transcripts.loc[~no_locus].copy()
    transcripts["locus_id"] = transcripts["gene_id"]

    before = len(transcripts)
    transcripts = transcripts.drop_duplicates(subset=["tid"])
    if len(transcripts) != before:
        schema.warn(f"{before - len(transcripts):,} duplicate transcript record(s) collapsed")

    exon, trans_len = build_exons(exons, args.keep_versions)

    orphan_exons = ~exon["tid"].isin(set(transcripts["tid"]))
    if int(orphan_exons.sum()):
        schema.warn(f"{int(orphan_exons.sum()):,} exon(s) belong to no transcript record "
                    f"-- dropped")
        exon = exon.loc[~orphan_exons]
        trans_len = trans_len.loc[trans_len.index.isin(set(transcripts["tid"]))]

    transcripts["trans_len"] = transcripts["tid"].map(trans_len)
    exonless = transcripts["trans_len"].isna()
    if int(exonless.sum()):
        schema.warn(f"{int(exonless.sum()):,} transcript(s) have no exon records -- dropped")
        transcripts = transcripts.loc[~exonless].copy()
    transcripts["trans_len"] = transcripts["trans_len"].astype("int64")

    transcripts, case_counts, novel_gene_ids = resolve_genes(transcripts)
    gene = build_gene(transcripts, novel_gene_ids)

    gene["biotype"] = None
    transcripts["transcript_type"] = None
    transcripts["tsl"] = None
    transcripts["ccds"] = None
    attrs_rows = 0
    by_symbol = bool(args.attrs_gene_name_column)
    if args.attrs:
        attrs = load_attrs(args.attrs, args)
        attrs_rows = len(attrs)
        gene, transcripts = apply_attrs(gene, transcripts, attrs,
                                        args.min_attrs_match, by_symbol=by_symbol)

    fasta_present, fasta_bad = 0, 0
    if args.transcriptome_fasta:
        fasta_present, fasta_bad = check_fasta(
                transcripts, args.transcriptome_fasta, args.keep_versions, out_dir)

    con = schema.connect_for_build(db_path)
    try:
        schema.create_schema(con)
        write_tables(con, gene, transcripts, exon, args.overwrite)
        schema.set_meta(con, "annotate", {
            "schema.version": schema.SCHEMA_VERSION,
            "tarva.version": schema.tarva_version(),
            "annotate.gtf": Path(args.gtf).resolve(),
            "annotate.attrs": Path(args.attrs).resolve() if args.attrs else None,
            "annotate.attrs_rows": attrs_rows,
            "annotate.attrs_joined_on": ("symbol" if by_symbol else "accession") if args.attrs
                                        else None,
            "annotate.transcriptome_fasta": (Path(args.transcriptome_fasta).resolve()
                                             if args.transcriptome_fasta else None),
            "annotate.assembly_report": (Path(args.assembly_report).resolve()
                                         if args.assembly_report else None),
            "annotate.drop_molecules": args.drop_molecules if args.assembly_report else None,
            "annotate.contigs_kept": ",".join(contigs_kept) if contigs_kept else None,
            "annotate.n_contigs_present": n_contigs_present,
            "annotate.n_contigs_kept": len(contigs_kept) if contigs_kept else None,
            "annotate.n_transcripts_dropped_contig": n_contig_dropped,
            "annotate.keep_versions": int(bool(args.keep_versions)),
            "annotate.n_genes": len(gene),
            "annotate.n_genes_novel": int(gene["novel"].sum()),
            "annotate.n_transcripts": len(transcripts),
            "annotate.n_exons": len(exon),
            "annotate.n_transcripts_matched": case_counts["matched"],
            "annotate.n_transcripts_inherited": case_counts["inherited"],
            "annotate.n_transcripts_novel_locus": case_counts["novel_loci_transcripts"],
            "annotate.n_transcripts_ambiguous": case_counts["ambiguous"],
            "annotate.n_dropped_no_tid": parse_counts["transcripts_no_tid"],
            "annotate.n_dropped_exons_no_tid": parse_counts["exons_no_tid"],
            "annotate.n_dropped_no_locus": n_no_locus,
            "annotate.n_genes_absent_no_locus": len(lost_genes),
            "annotate.fasta_checked": fasta_present,
            "annotate.fasta_mismatches": fasta_bad,
            })
    finally:
        con.close()

    schema.report("")
    schema.report("annotate complete")
    schema.report(f"  genes                      {len(gene):,} "
                  f"({int(gene['novel'].sum()):,} novel loci)")
    schema.report(f"  transcripts                {len(transcripts):,}")
    schema.report(f"  exons                      {len(exon):,}")
    schema.report(f"  reference-matched          {case_counts['matched']:,}")
    schema.report(f"  inherited (novel isoform)  {case_counts['inherited']:,}")
    schema.report(f"  in a purely novel locus    {case_counts['novel_loci_transcripts']:,}")
    schema.report(f"  ambiguous (chimeric locus) {case_counts['ambiguous']:,}")
    if contigs_kept:
        schema.report(f"  contigs                    {len(contigs_kept):,} kept of "
                      f"{n_contigs_present:,}")
    else:
        schema.report(f"  contigs                    {n_contigs_present:,} "
                      f"(no selection applied)")
    schema.report(f"  dropped, no transcript_id  {parse_counts['transcripts_no_tid']:,}")
    schema.report(f"  dropped, no gene_id        {n_no_locus:,}")
    schema.report(f"  dropped, contig excluded   {n_contig_dropped:,}")
    schema.report(f"Wrote {db_path}")
    return 0
