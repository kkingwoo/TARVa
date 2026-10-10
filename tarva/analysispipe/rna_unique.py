"""
analysispipe rna-unique: RNA variants with no WGS counterpart, panel-masked

    `tarva analysispipe rna-unique --metadata <meta.csv>
        --wgs-variants <wgs_vcf_dir> --rna-variants <rna_vcf_dir>
        --output <out_dir> [--pon <pon.vcf.gz>]
        [--germline-resource <af-only.vcf.gz>] [--jobs N] [--overwrite]`

Each individual's WGS and RNA VCF records are matched, and every RNA variant's
coordinates (chrom, pos) that have no WGS record are written to
<output>/<group>/<id>/<id>_rna_unique.vcf. This is the RNA-specific set stored
into a database for downstream analysis.

Panel membership is read straight from pon.vcf.gz -- every record in a
CreateSomaticPanelOfNormals output is a panel site, the minimum sample count
having been applied when it was built, so there is no threshold to re-apply.
Germline membership is read the same way from whatever --germline-resource
`wgseqpipe variants` was given; pass the same file, or the flagging will not
correspond to what was actually suppressed.

**--pon and --germline-resource are independently optional.** Pass neither,
either, or both, matching whatever `wgseqpipe variants` was actually given:

    neither   No masking. Correct only if the WGS calls used no panel and no
              germline resource; otherwise the sites Mutect2 suppressed are
              absent from the WGS VCFs and every RNA call at one of them is
              counted RNA-specific with no genomic baseline behind it.
    --pon     Panel positions are REMOVED from the output, every sample alike.
              A panel site is a recurrent artefact, so removing it
              symmetrically means it can neither appear as RNA-specific nor
              create a between-group asymmetry.
    --germline-resource
              Known germline positions are KEPT and FILTER-tagged
              germline_known. They are suppressed upstream of the panel, so
              --pon cannot cover them, and they are real DNA variants rather
              than artefacts -- so they are marked and the decision is left to
              whatever consumes this file.
    both      As above; a position in both is REMOVED, the panel taking
              precedence, so removed_by_pon and germline_flagged never count
              the same record twice.

Panel positions are deleted rather than flagged, so an output VCF built with
--pon cannot be told from one built without it by inspection. Both resources
are therefore recorded twice: as ##tarva_rna_unique_* headers in every output
VCF, and as the pon_vcf / pon_positions / germline_vcf / germline_positions
columns of rna_unique_summary.csv. A zero removed_by_pon beside a named file
with a nonzero index then means the panel genuinely had no overlap, not that no
panel was applied.
"""

import csv
import gzip
import sys

from array import array
from bisect import bisect_left
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from tarva.metadata import load_metadata

WGS_VCF_SUFFIX = "_contaminates_filtered.vcf.gz"
RNA_VCF_SUFFIX = "_varcall_flt.vcf"
RNA_UNIQUE_SUFFIX = "_rna_unique.vcf"

GERMLINE_FILTER_ID = "germline_known"
GERMLINE_FILTER_HEADER = (
    f'##FILTER=<ID={GERMLINE_FILTER_ID},Description="Position present in the '
    f'--germline-resource VCF. Flagged and retained, not removed.">'
    )

SUMMARY_FIELDS = [
    "individualID", "group", "status", "rna_unique_vcf", "rna_records",
    "wgs_positions", "removed_by_wgs", "pon_vcf", "pon_positions",
    "removed_by_pon", "germline_vcf", "germline_positions",
    "germline_flagged", "written", "error",
    ]


def add_arguments(parser):
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV (see metadata.py)")
    parser.add_argument("--wgs-variants", required=True,
                        help="The --output directory `wgseqpipe variants` wrote; VCFs are read "
                             f"from <dir>/<group>/<id>/<id>{WGS_VCF_SUFFIX}")
    parser.add_argument("--rna-variants", required=True,
                        help="The --output directory `rnaseqpipe variants` wrote; VCFs are read "
                             f"from <dir>/<group>/<id>/<id>{RNA_VCF_SUFFIX}")
    parser.add_argument("--output", required=True,
                        help="Directory for RNA-unique VCFs, as <output>/<group>/<id>/")
    parser.add_argument("--pon", metavar="VCF",
                        help="pon.vcf.gz from `tarva wgseqpipe pon`. Positions in it are "
                             "excluded from every sample's output. Optional, and independent "
                             "of --germline-resource, but strongly recommended whenever the "
                             "WGS calls were made with a panel -- without it, panel sites "
                             "reappear as RNA-specific")
    parser.add_argument("--germline-resource", metavar="VCF",
                        help="The af-only VCF `wgseqpipe variants` used. RNA records at "
                             "these positions are FLAGGED and kept, not removed: the "
                             "germline suppression happens upstream of the panel, so these "
                             "sites are absent from the WGS VCF but are not in the panel "
                             "either, and a position-only comparison would otherwise count "
                             "them RNA-specific. Optional, and independent of --pon")
    parser.add_argument("--jobs", type=int, default=1,
                        help="Samples to process concurrently (default: 1). This is I/O "
                             "bound, so a few jobs help and many do not")
    parser.add_argument("--overwrite", action="store_true",
                        help="Reprocess samples that already have an RNA-unique VCF")


def _open_vcf(path):
    """VCF text handle, gzipped or not. Never decompresses to disk."""
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path)


def load_position_index(vcf_path):
    """
    {chrom: sorted array of positions}, for membership tests on a large VCF.

    """
    index = {}
    with _open_vcf(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 2)
            if len(fields) < 2:
                continue
            try:
                chrom, pos = fields[0], int(fields[1])
            except ValueError:
                continue
            positions = index.get(chrom)
            if positions is None:
                positions = index[chrom] = array("l")
            elif pos < positions[-1]:
                raise ValueError(
                        f"{vcf_path} is not coordinate-sorted -- {chrom}:{pos} follows "
                        f"{chrom}:{positions[-1]}. Sort it with "
                        f"`gatk SortVcf --SEQUENCE_DICTIONARY <ref>.dict`."
                        )
            positions.append(pos)
    return index


def position_in(index, chrom, pos):
    """Membership test against a load_position_index() result."""
    positions = index.get(chrom)
    if not positions:
        return False
    i = bisect_left(positions, pos)
    return i < len(positions) and positions[i] == pos


def index_size(index):
    return sum(len(v) for v in index.values())


def load_resource(vcf_path, label, required=False):
    """
    One shared position index, loaded once for each individual.
    """
    if not vcf_path:
        if required:
            raise SystemExit(f"tarva analysispipe rna-unique: {label} is required")
        print(f"{label}: none given")
        return {}
    vcf_path = Path(vcf_path)
    if not vcf_path.is_file():
        raise SystemExit(f"tarva analysispipe rna-unique: {label} not found: {vcf_path}")
    try:
        index = load_position_index(vcf_path)
    except ValueError as exc:
        raise SystemExit(f"tarva analysispipe rna-unique: {label}: {exc}")
    print(f"{label}: {index_size(index):,} positions across {len(index)} contig(s) "
          f"from {vcf_path.name}")
    return index


def rna_contig_names(args, meta):
    """
    Contig names declared by the first readable RNA VCF, for the mismatch check.

    """
    for sample in meta.samples:
        rna_vcf = _sample_vcf(args.rna_variants, sample, RNA_VCF_SUFFIX)
        if not rna_vcf.is_file():
            continue
        contigs = set()
        try:
            with _open_vcf(rna_vcf) as fh:
                for line in fh:
                    if not line.startswith("#"):
                        break
                    if line.startswith("##contig=<ID="):
                        name = line.split("##contig=<ID=", 1)[1].split(",", 1)[0]
                        contigs.add(name.rstrip(">\n"))
        except OSError:
            continue
        if contigs:
            return contigs
    return set()


def warn_contig_mismatch(label, index, rna_contigs):
    """
    Catch an Ensembl/RefSeq naming mismatch before it masks nothing silently.

    A resource whose contigs are named differently from the RNA calls misses on
    every lookup: removed_by_pon comes out 0, the run reports success, and no
    masking was applied at all -- indistinguishable from not passing the flag.
    This project has already hit the same mismatch with the MGP known-sites VCF
    (Ensembl `1` versus RefSeq `NC_000067.7`).
    """
    if not index or not rna_contigs:
        return
    if set(index) & rna_contigs:
        return
    print(f"WARNING: {label} shares NO contig names with the RNA VCFs "
          f"(e.g. {sorted(index)[:3]} vs {sorted(rna_contigs)[:3]}). Every lookup will "
          f"miss, so nothing will be masked or flagged. Check for an Ensembl/RefSeq "
          f"naming mismatch.", file=sys.stderr)


def provenance_headers(args, wgs_vcf, pon_index, germline_index):
    """
    What was applied to this file, as VCF header lines.

    """
    pon = Path(args.pon).resolve() if args.pon else None
    germline = Path(args.germline_resource).resolve() if args.germline_resource else None
    return [
        f"##tarva_rna_unique_wgs={wgs_vcf}",
        (f"##tarva_rna_unique_pon={pon} ({index_size(pon_index):,} positions, removed)"
         if pon else "##tarva_rna_unique_pon=none"),
        (f"##tarva_rna_unique_germline={germline} "
         f"({index_size(germline_index):,} positions, flagged)"
         if germline else "##tarva_rna_unique_germline=none"),
        ]


def _tag_filter(fields, filter_id):
    """
    Add `filter_id` to a VCF record's FILTER column, preserving what is there.

    """
    existing = fields[6]
    if existing in ("", ".", "PASS"):
        fields[6] = filter_id
    elif filter_id not in existing.split(";"):
        fields[6] = f"{existing};{filter_id}"
    return fields


def _sample_vcf(root, sample, suffix):
    return (Path(root) / sample.group / sample.individual_id
            / f"{sample.individual_id}{suffix}")


def _rna_unique_one(sample, args, pon_index, germline_index):
    """
    One sample's comparison. Failures become FAILED rows, never exceptions.
    """
    out_dir = Path(args.output) / sample.group / sample.individual_id
    out_vcf = out_dir / f"{sample.individual_id}{RNA_UNIQUE_SUFFIX}"
    # Which resources were configured is a fact about the run, not about this
    # sample's success, so it goes on skipped and FAILED rows too.
    row = {"individualID": sample.individual_id, "group": sample.group,
           "status": "FAILED", "rna_unique_vcf": "", "rna_records": 0,
           "wgs_positions": 0, "removed_by_wgs": 0,
           "pon_vcf": str(Path(args.pon).resolve()) if args.pon else "",
           "pon_positions": index_size(pon_index), "removed_by_pon": 0,
           "germline_vcf": (str(Path(args.germline_resource).resolve())
                            if args.germline_resource else ""),
           "germline_positions": index_size(germline_index),
           "germline_flagged": 0, "written": 0, "error": ""}
    try:
        if out_vcf.is_file() and not args.overwrite:
            print(f"{sample.individual_id}: {out_vcf.name} already present, skipping")
            return {**row, "status": "skipped", "rna_unique_vcf": str(out_vcf)}

        wgs_vcf = _sample_vcf(args.wgs_variants, sample, WGS_VCF_SUFFIX)
        rna_vcf = _sample_vcf(args.rna_variants, sample, RNA_VCF_SUFFIX)
        for label, path in (("WGS", wgs_vcf), ("RNA", rna_vcf)):
            if not path.is_file():
                raise FileNotFoundError(f"{label} VCF not found: {path}")

        wgs_index = load_position_index(wgs_vcf)
        out_dir.mkdir(parents=True, exist_ok=True)
        partial = out_vcf.with_name(out_vcf.name + ".partial")
        provenance = provenance_headers(args, wgs_vcf, pon_index, germline_index)

        counts = dict(rna_records=0, removed_by_wgs=0, removed_by_pon=0,
                      germline_flagged=0, written=0)
        with _open_vcf(rna_vcf) as src, open(partial, "w") as dst:
            for line in src:
                # Header copied verbatim so the output stays a valid VCF with
                # the RNA call set's own contigs, INFO and FORMAT definitions.
                # Provenance and the germline FILTER are declared just before
                # #CHROM, which must be the last header line.
                if line.startswith("#"):
                    if line.startswith("#CHROM"):
                        for header in provenance:
                            dst.write(header + "\n")
                        if germline_index:
                            dst.write(GERMLINE_FILTER_HEADER + "\n")
                    dst.write(line)
                    continue
                fields = line.split("\t", 2)
                if len(fields) < 2:
                    continue
                try:
                    chrom, pos = fields[0], int(fields[1])
                except ValueError:
                    continue
                counts["rna_records"] += 1
                if position_in(wgs_index, chrom, pos):
                    counts["removed_by_wgs"] += 1
                    continue
                # The panel takes precedence over the germline resource: a
                # position in both is removed, not flagged-and-kept, so the two
                # counts never cover the same record.
                if position_in(pon_index, chrom, pos):
                    counts["removed_by_pon"] += 1
                    continue
                # Germline sites are flagged, not dropped.
                if position_in(germline_index, chrom, pos):
                    # Only the matching records pay for a full field split.
                    all_fields = line.rstrip("\n").split("\t")
                    if len(all_fields) >= 7:
                        dst.write("\t".join(_tag_filter(all_fields, GERMLINE_FILTER_ID))
                                  + "\n")
                    else:
                        dst.write(line)
                    counts["germline_flagged"] += 1
                    counts["written"] += 1
                    continue
                dst.write(line)
                counts["written"] += 1
        partial.replace(out_vcf)

        return {**row, "status": "done", "rna_unique_vcf": str(out_vcf),
                "wgs_positions": index_size(wgs_index), **counts}
    except Exception as exc:
        return {**row, "error": f"{type(exc).__name__}: {exc}"}


def run(args):
    meta = load_metadata(args.metadata)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_csv = output_dir / "rna_unique_summary.csv"

    if not args.pon:
        print("WARNING: no --pon given. If the WGS calls were made with a panel of "
              "normals, panel sites are absent from the WGS VCFs and every RNA call at "
              "one of them will be counted RNA-specific with no genomic baseline behind "
              "it. Pass the pon.vcf.gz that `wgseqpipe variants` used.", file=sys.stderr)
    if not args.germline_resource:
        print("WARNING: no --germline-resource given. If the WGS calls were made with "
              "one, known germline sites are absent from the WGS VCFs and are NOT in the "
              "panel either -- they are suppressed before the panel is built -- so --pon "
              "does not cover them and they will be counted RNA-specific. Pass the "
              "af-only VCF that `wgseqpipe variants` used.", file=sys.stderr)
    pon_index = load_resource(args.pon, "Panel of normals")
    germline_index = load_resource(args.germline_resource, "Germline resource")

    # Only worth reading an RNA header when there is something to compare to.
    if pon_index or germline_index:
        rna_contigs = rna_contig_names(args, meta)
        warn_contig_mismatch("Panel of normals", pon_index, rna_contigs)
        warn_contig_mismatch("Germline resource", germline_index, rna_contigs)

    print(f"RNA-unique: {len(meta.samples)} sample(s), {args.jobs} at a time")
    rows = []
    with open(summary_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        fh.flush()
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(_rna_unique_one, s, args, pon_index, germline_index)
                       for s in meta.samples]
            for future in as_completed(futures):
                row = future.result()
                writer.writerow(row)
                fh.flush()
                rows.append(row)

    failed = [r for r in rows if r["status"] == "FAILED"]
    done = [r for r in rows if r["status"] == "done"]
    print(f"\nRNA-unique: {len(done)} done, {len(failed)} failed")
    print(f"  --pon:                {args.pon or 'none'}")
    print(f"  --germline-resource:  {args.germline_resource or 'none'}")
    if done:
        print(f"{'sample':<12}{'RNA':>12}{'-WGS':>12}{'-PON':>12}"
              f"{'unique':>12}{'germline':>12}")
        for row in sorted(done, key=lambda r: r["individualID"]):
            print(f"{row['individualID']:<12}{row['rna_records']:>12,}"
                  f"{row['removed_by_wgs']:>12,}{row['removed_by_pon']:>12,}"
                  f"{row['written']:>12,}{row['germline_flagged']:>12,}")
        print("germline = of the unique records, how many carry "
              f"FILTER={GERMLINE_FILTER_ID} (flagged, still present)")
    print(f"Wrote summary: {summary_csv}")
    for row in failed:
        print(f"  FAILED {row['individualID']}: {row['error']}", file=sys.stderr)
    if failed:
        raise SystemExit(
                f"tarva analysispipe rna-unique: {len(failed)} sample(s) failed "
                f"-- see {summary_csv}"
                )
    return rows
