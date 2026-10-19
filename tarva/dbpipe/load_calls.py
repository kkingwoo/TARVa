"""
dbpipe load-calls: load the RNA-specific variant calls into the site and call
tables.

Each distinct position is resolved to one gene and strand, and each call is
given its transcript-oriented call type, which is what the analysis steps read.

    `tarva dbpipe load-calls --rna-unique <dir> --db <TARVa.db>
        [--min-depth 20] [--overwrite]`

Input:  one VCF per individual at
        <rna-unique>/<condition>/<rid>/<rid>_rna_unique.vcf, as `analysispipe
        rna-unique` writes it, and a database that `annotate` and
        `load-samples` have written.
Output: the site and call tables, and load_calls.* keys in run_meta. No depth
        filter is applied here; --min-depth is recorded for the analysis steps.
"""

import gzip
import heapq
import re

from pathlib import Path

import pandas as pd

from tarva.dbpipe import schema

LOAD_CALLS_TABLES = ("site", "call")

RNA_UNIQUE_SUFFIX = "_rna_unique.vcf"
DEFAULT_MIN_DEPTH = 20
WRITE_CHUNK = 100000

COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}

# Step 2 of the locus-resolution rule. Lower wins. IG_*/TR_*_gene segments are
# coding, so they rank with protein_coding; their pseudogene counterparts fall
# to the bottom tier with everything else named *pseudogene*.
BIOTYPE_TIERS = (
    (0, ("protein_coding",)),
    (1, ("lncRNA",)),
    (2, ("miRNA", "snoRNA", "snRNA", "scaRNA", "rRNA", "misc_RNA",
         "ribozyme", "sRNA", "vault_RNA", "Mt_tRNA", "Mt_rRNA")),
    )
BIOTYPE_TIER = {name: tier for tier, names in BIOTYPE_TIERS for name in names}
OTHER_TIER = 3
UNKNOWN_TIER = 4

# Provenance the rna-unique output carries in its own header. Copied into
# run_meta so what was masked travels with the data rather than living only in
# rna_unique_summary.csv beside the VCFs.
PROVENANCE_HEADER = re.compile(r"^##(tarva_rna_unique_[A-Za-z0-9_]+)=(.*)$")


def add_arguments(parser):
    parser.add_argument("--rna-unique", required=True, metavar="DIR",
                        help="The --output directory of `analysispipe rna-unique`. Each "
                             "individual's VCF is read from "
                             "<dir>/<condition>/<rid>/<rid>_rna_unique.vcf, with a "
                             "recursive search as a fallback")
    parser.add_argument("--db", required=True,
                        help="Database that `annotate` and `load-samples` have written")
    parser.add_argument("--min-depth", type=int, default=DEFAULT_MIN_DEPTH,
                        help=f"Depth threshold RECORDED in run_meta, not applied here "
                             f"(default: {DEFAULT_MIN_DEPTH}). ")
    parser.add_argument("--overwrite", action="store_true",
                        help="Drop and reload site and call instead of failing")


#----------------------------------------------------
# inputs
#----------------------------------------------------

def _open_vcf(path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def load_cohort(con):
    """
    [(individualID, condition)] from the individual table.

    Read from the database rather than from a metadata file, so
    call.individualID cannot name someone `individual` does not have. 
    """
    table = schema.INDIVIDUAL_TABLE
    schema.require_table(con, table, "tarva dbpipe load-samples")
    id_column, condition_column = schema.individual_columns(con)
    query = (f'SELECT "{id_column}", "{condition_column}" FROM {table} '
             f'ORDER BY "{id_column}"')
    return [(str(individual), str(condition))
            for individual, condition in con.execute(query)]


def resolve_vcf(root, rid, condition):
    """The documented path first, then a recursive search for the filename."""
    expected = Path(root) / condition / rid / f"{rid}{RNA_UNIQUE_SUFFIX}"
    if expected.is_file():
        return expected
    found = sorted(Path(root).rglob(f"{rid}{RNA_UNIQUE_SUFFIX}"))
    if len(found) == 1:
        schema.warn(f"{rid}: not at {expected}, using {found[0]}")
        return found[0]
    if len(found) > 1:
        raise SystemExit(
                f"tarva dbpipe load-calls: {len(found)} files named "
                f"{rid}{RNA_UNIQUE_SUFFIX} under {root} "
                f"({', '.join(str(p) for p in found[:3])}). Point --rna-unique at one "
                f"run's output directory.")
    raise SystemExit(
            f"tarva dbpipe load-calls: no {rid}{RNA_UNIQUE_SUFFIX} under {root}. "
            f"Expected {expected}. Has `analysispipe rna-unique` run for this cohort?")


def load_kept_contigs(con):
    """annotate.contigs_kept, or None when annotate applied no selection."""
    value = schema.get_meta(con, "annotate.contigs_kept")
    if not value:
        schema.warn("run_meta has no annotate.contigs_kept -- annotate ran without "
                    "--assembly-report, so no contig restriction is applied here either")
        return None
    return set(value.split(","))


def read_provenance(path):
    """The ##tarva_rna_unique_* headers one output VCF carries."""
    provenance = {}
    with _open_vcf(path) as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            found = PROVENANCE_HEADER.match(line.rstrip("\n"))
            if found:
                provenance[found.group(1)] = found.group(2)
    return provenance


#----------------------------------------------------
# VCF records
#----------------------------------------------------

def sample_column(header_fields, rid, path):
    """The index of this individual's sample column in the #CHROM line."""
    names = header_fields[9:]
    if not names:
        raise SystemExit(
                f"tarva dbpipe load-calls: {path} has no sample column in its #CHROM "
                f"line. It is not a genotyped VCF.")
    if rid in names:
        return 9 + names.index(rid)
    if len(names) == 1:
        schema.warn(f"{rid}: the VCF's sample column is named {names[0]!r}, not "
                    f"{rid!r} -- using it, and storing individualID as {rid!r}. "
                    f"Confirm these "
                    f"are the same individual.")
        return 9
    raise SystemExit(
            f"tarva dbpipe load-calls: {path} has sample columns "
            f"{', '.join(names)}, none of which is {rid!r}. Cannot tell which belongs "
            f"to this individual.")


def parse_records(path, rid, keep_contigs, counts):
    """
    Yield (chrom, pos, ref, alt1, a1_ad, dp, filter) per usable record.

    `ref` is the VCF REF allele in full, not the single base at pos. For an
    indel those differ: REF spans the affected bases, so a deletion is
    REF=TGG ALT=T while an insertion at the same coordinate is REF=T ALT=TAA.
    site.ref_base takes ref[0].
    """
    column = None
    with _open_vcf(path) as handle:
        for line in handle:
            if line.startswith("#"):
                if line.startswith("#CHROM"):
                    column = sample_column(line.rstrip("\n").split("\t"), rid, path)
                continue
            if column is None:
                raise SystemExit(
                        f"tarva dbpipe load-calls: {path} has records before its "
                        f"#CHROM line. It is not a valid VCF.")
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= column:
                continue
            counts["records"] += 1

            chrom = fields[0]
            if keep_contigs is not None and chrom not in keep_contigs:
                counts["off_contig"] += 1
                continue
            try:
                pos = int(fields[1])
            except ValueError:
                continue

            keys = fields[8].split(":")
            values = fields[column].split(":")
            if len(values) < len(keys):
                counts["no_format"] += 1
                continue
            try:
                depths = values[keys.index("AD")].split(",")
                dp = int(values[keys.index("DP")])
                a1_ad = int(depths[1])
            except (ValueError, IndexError):
                counts["no_format"] += 1
                continue

            counts["written"] += 1
            yield (chrom, pos, fields[3].upper(), fields[4].split(",")[0].upper(),
                   a1_ad, dp, fields[6] or ".")


#----------------------------------------------------
# resolving a position to one gene
#----------------------------------------------------

def biotype_tier(biotype):
    if biotype is None:
        return UNKNOWN_TIER
    return BIOTYPE_TIER.get(biotype, OTHER_TIER)


def load_genes(con):
    """
    {chrom: [(gene_start, gene_end, sort_key, ref_id, strand)]}, start-sorted.

    sort_key is the resolution rule's tiers 1-5 as a tuple, so choosing a
    winner is min() rather than five comparisons:
      1. not novel        2. biotype        3. most transcripts
      4. largest span     5. lowest accession (determinism only)
    Spans are negated so that larger wins under min().
    """
    schema.require_table(con, "gene", "tarva dbpipe annotate")
    query = ("SELECT chrom, gene_start, gene_end, ref_id, strand, biotype, "
             "novel, n_transcripts FROM gene ORDER BY chrom, gene_start")
    by_chrom = {}
    for (chrom, start, end, ref_id, strand, biotype,
         novel, n_transcripts) in con.execute(query):
        key = (int(novel), biotype_tier(biotype), -int(n_transcripts),
               -(int(end) - int(start) + 1), str(ref_id))
        by_chrom.setdefault(chrom, []).append(
                (int(start), int(end), key, str(ref_id), strand))
    return by_chrom


def resolve_sites(sites, genes):
    """
    {(chrom, pos): (ref_id, strand)} by a coordinate sweep.

    A position can fall inside several genes' spans. Sites and genes are both
    coordinate-ordered, so one pass with a heap of active genes keyed on
    gene_end finds every covering gene without an interval tree. The
    alternatives stay discoverable with a range query against `gene`, so
    nothing is stored twice.
    """
    resolved = {}
    by_chrom = {}
    for chrom, pos in sites:
        by_chrom.setdefault(chrom, []).append(pos)

    for chrom, positions in by_chrom.items():
        candidates = genes.get(chrom, [])
        positions.sort()
        active, index = [], 0
        for pos in positions:
            while index < len(candidates) and candidates[index][0] <= pos:
                start, end, key, ref_id, strand = candidates[index]
                heapq.heappush(active, (end, key, ref_id, strand))
                index += 1
            while active and active[0][0] < pos:
                heapq.heappop(active)
            if active:
                _, _, ref_id, strand = min(active, key=lambda item: item[1])
                resolved[(chrom, pos)] = (ref_id, strand)
            else:
                resolved[(chrom, pos)] = (None, None)
    return resolved


#----------------------------------------------------
# a1type
#----------------------------------------------------

def classify(ref, alt1, strand):
    """
    The transcript-oriented substitution label for one call.

    `ref` is the call's own VCF REF allele, not site.ref_base, because the
    ins/del decision is a length comparison between the two alleles as the
    record reported them.

    Transcript-level and strand-resolved.

    """
    if len(ref) != len(alt1):
        return "ins" if len(alt1) > len(ref) else "del"
    if len(ref) != 1:
        # Equal-length, multi-base: an MNV, which has no call-type label. NULL
        # rather than mislabelled as an indel.
        return None
    if strand not in ("+", "-"):
        return None
    if ref not in COMPLEMENT or alt1 not in COMPLEMENT:
        return None
    if strand == "-":
        ref, alt1 = COMPLEMENT[ref], COMPLEMENT[alt1]
    return f"{ref}-{alt1}"


#----------------------------------------------------
# writing
#----------------------------------------------------

def write_frame(con, table, rows, columns):
    if not rows:
        return 0
    frame = pd.DataFrame(rows, columns=columns)
    frame.to_sql(table, con, if_exists="append", index=False, chunksize=WRITE_CHUNK)
    con.commit()
    return len(frame)


def check_reloadable(con, overwrite):
    """
    Refuse early, before two passes over the VCFs, not after them.

    """
    existing = con.execute("SELECT COUNT(*) FROM call").fetchone()[0]
    site_rows = con.execute("SELECT COUNT(*) FROM site").fetchone()[0]
    if (existing or site_rows) and not overwrite:
        raise SystemExit(
                f"tarva dbpipe load-calls: call already holds {existing:,} row(s) and "
                f"site {site_rows:,}. Pass --overwrite to reload them, which also "
                f"invalidates any gene_test rows computed from them.")


def prepare_tables(con):
    existing = con.execute("SELECT COUNT(*) FROM call").fetchone()[0]
    site_rows = con.execute("SELECT COUNT(*) FROM site").fetchone()[0]
    if existing or site_rows:
        for table in reversed(LOAD_CALLS_TABLES):
            con.execute(f"DROP TABLE IF EXISTS {table}")
        con.commit()
        schema.create_schema(con, LOAD_CALLS_TABLES)
        schema.report(f"--overwrite: dropped and recreated "
                      f"{', '.join(LOAD_CALLS_TABLES)}")


#----------------------------------------------------
# entry point
#----------------------------------------------------

def run(args):
    db_path = Path(args.db)
    if not db_path.is_file():
        raise SystemExit(
                f"tarva dbpipe load-calls: database not found: {db_path}. Run "
                f"`dbpipe annotate` and `dbpipe load-samples` first.")

    con = schema.connect_for_build(db_path)
    try:
        schema.create_schema(con)
        check_reloadable(con, args.overwrite)
        cohort = load_cohort(con)
        keep_contigs = load_kept_contigs(con)
        genes = load_genes(con)
        schema.report(f"cohort: {len(cohort)} individual(s) from `individual`")
        schema.report(f"genes: {sum(len(v) for v in genes.values()):,} spans over "
                      f"{len(genes):,} contig(s)")

        paths = [(rid, resolve_vcf(args.rna_unique, rid, condition))
                 for rid, condition in cohort]

        # Pass 1: distinct positions, so every site's gene and strand is known
        # before any a1type is computed.
        sites = {}
        totals = dict(records=0, off_contig=0, no_format=0, usable=0)
        for rid, path in paths:
            counts = dict(records=0, off_contig=0, no_format=0, written=0)
            for chrom, pos, ref, _, _, _, _ in parse_records(
                    path, rid, keep_contigs, counts):
                # ref[0] only: an indel's full REF differs between samples at
                # one coordinate by design, but its anchor base cannot. A
                # disagreement here means two samples report different genomic
                # bases at the same position, which is a real problem rather
                # than a VCF representation detail.
                anchor = ref[0]
                previous = sites.setdefault((chrom, pos), anchor)
                if previous != anchor:
                    schema.warn(f"{chrom}:{pos} has reference base {previous} in one "
                                f"sample and {anchor} in another -- keeping {previous}. "
                                f"The same coordinate cannot have two genomic bases, so "
                                f"check these VCFs are against the same reference.")
            for key in ("records", "off_contig", "no_format"):
                totals[key] += counts[key]
            totals["usable"] += counts["written"]
            schema.report(f"  {rid}: {counts['records']:,} record(s), "
                          f"{counts['written']:,} usable, "
                          f"{counts['off_contig']:,} off-contig, "
                          f"{counts['no_format']:,} without AD/DP")

        if not sites:
            raise SystemExit(
                    f"tarva dbpipe load-calls: no usable records in any VCF under "
                    f"{args.rna_unique}. Check the off-contig and AD/DP counts above.")

        schema.report(f"sites: {len(sites):,} distinct position(s) from "
                      f"{totals['usable']:,} call(s)")
        resolved = resolve_sites(sites.keys(), genes)
        with_gene = sum(1 for value in resolved.values() if value[0] is not None)
        schema.report(f"sites: {with_gene:,} of {len(sites):,} fall inside a gene span "
                      f"({with_gene / len(sites):.1%})")

        prepare_tables(con)

        site_rows = [(chrom, pos, resolved[(chrom, pos)][1], anchor,
                      resolved[(chrom, pos)][0])
                     for (chrom, pos), anchor in sites.items()]
        written_sites = write_frame(con, "site", site_rows,
                                    ["chrom", "pos", "strand", "ref_base", "ref_id"])
        schema.report(f"site: {written_sites:,} row(s) written")

        # Pass 2: the calls, one sample at a time.
        written_calls = 0
        type_counts = {}
        for rid, path in paths:
            rows = []
            pass2 = dict(records=0, off_contig=0, no_format=0, written=0)
            for chrom, pos, ref, alt1, a1_ad, dp, vcf_filter in parse_records(
                    path, rid, keep_contigs, pass2):
                strand = resolved[(chrom, pos)][1]
                # This call's own REF, not site.ref_base: the ins/del decision
                # is a length comparison between the two alleles as the record
                # reported them.
                a1type = classify(ref, alt1, strand)
                type_counts[a1type] = type_counts.get(a1type, 0) + 1
                rows.append((chrom, pos, rid, ref, alt1, a1_ad, dp, a1type, vcf_filter))
            written_calls += write_frame(
                    con, "call", rows,
                    ["chrom", "pos", "individualID", "ref", "alt1", "a1_ad", "dp",
                     "a1type", "filter"])
            schema.report(f"  {rid}: {len(rows):,} call(s) written")
        schema.report(f"call: {written_calls:,} row(s) written")

        schema.create_indexes(con, LOAD_CALLS_TABLES)

        provenance = read_provenance(paths[0][1]) if paths else {}
        meta = {
            "load_calls.rna_unique_dir": Path(args.rna_unique).resolve(),
            "load_calls.min_depth": args.min_depth,
            "load_calls.n_samples": len(cohort),
            "load_calls.n_sites": written_sites,
            "load_calls.n_calls": written_calls,
            "load_calls.n_sites_in_gene": with_gene,
            "load_calls.n_records_read": totals["records"],
            "load_calls.n_dropped_off_contig": totals["off_contig"],
            "load_calls.n_dropped_no_format": totals["no_format"],
            "load_calls.contigs_restricted": int(keep_contigs is not None),
            }
        # The masking rna-unique applied, lifted out of the VCF header so it
        # travels with the data instead of only beside it.
        for key, value in provenance.items():
            meta[f"load_calls.{key}"] = value
        schema.set_meta(con, "load_calls", meta)
    finally:
        con.close()

    schema.report("")
    schema.report("load-calls complete")
    schema.report(f"  sites                     {written_sites:,}")
    schema.report(f"  calls                     {written_calls:,}")
    if written_sites:
        schema.report(f"  calls per site            "
                      f"{written_calls / written_sites:.2f}")
    for label, count in sorted(type_counts.items(),
                               key=lambda pair: -pair[1])[:16]:
        schema.report(f"  {str(label):<26}{count:,}")
    schema.report(f"Wrote {db_path}")
    return 0
