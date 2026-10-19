"""
Shared GATK preprocessing chain for the WGS and RNA variant-calling paths.

    FastqToSam -> SortSam(queryname, aligned) -> MergeBamAlignment
    -> MarkDuplicates -> SortSam(coordinate, +index)
    -> [RNA: SplitNCigarReads -> AddOrReplaceReadGroups]
    -> BaseRecalibrator -> ApplyBQSR

One chain for both paths, so nothing structural is duplicated between them.
`wgseqpipe/preprocess.py` and `rnaseqpipe/preprocess.py` are thin per-sample
callers over run_preprocessing_chain(), differing only in `is_rna`.

Four things worth knowing before changing any of it:

- **The ALIGNED BAM is what gets queryname-sorted.** MergeBamAlignment needs
  both inputs in queryname order, and the uBAM is already there, since
  FastqToSam's SORT_ORDER defaults to queryname.
- **MarkDuplicates runs bare**, with default parameters, on both paths.
- **Every BAM that reaches BaseRecalibrator is indexed**, so the final SortSam
  and AddOrReplaceReadGroups both pass --CREATE_INDEX true.
- **Read groups are parsed from the Illumina header at runtime** rather than
  from filenames, which needs no extra metadata columns, and the complete
  group -- ID, SM, LB, PL, PU -- is set once at FastqToSam.
  AddOrReplaceReadGroups then runs on the RNA path only, after
  SplitNCigarReads; `--read-group-repair never` skips it, since with a
  complete group already on the uBAM it only rewrites values the BAM carries.

Every BAM-producing step is resumable: it runs `samtools quickcheck` on its
own output, skips if the file is already complete, and deletes and redoes it
if not -- a chain this long will get killed at the wall clock sooner or later.

Input:  trimmed FASTQs, the aligned BAM for a sample, the reference, and
        optionally one or more --known-sites VCFs. Without them the chain
        stops before BQSR and reports that, so it is testable before a
        known-sites file exists.
Output: <id>_bqsr.bam and its index, or the last BAM the chain reached when
        BQSR was skipped, plus the per-sample metrics each tool writes.
"""

import csv
import gzip
import re
import subprocess
import sys

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

from tarva.metadata import DEFAULT_PAIRED_TOKENS, load_metadata, parse_paired_tokens
from tarva.reference import ensure_unzipped_fasta
from tarva.slurm_utils import (parse_directives_dict, split_directives, submit_job,
                               write_sbatch_script)

DEFAULT_JAVA_HEAP = "32g"
ILLUMINA_PLATFORM = "illumina"
READ_GROUP_REPAIR_CHOICES = ("auto", "always", "never")

# GATK4 bundles Picard, and the two halves spell the temp-directory option
# differently: Picard-derived tools take --TMP_DIR, engine tools take --tmp-dir.
# Passing the wrong spelling is a hard argument-parsing failure, not a warning.
PICARD_TOOLS = frozenset({
    "FastqToSam", "SortSam", "MergeBamAlignment", "MarkDuplicates",
    "AddOrReplaceReadGroups", "CreateSequenceDictionary",
    })


@dataclass
class ReadGroup:
    """One @RG line's worth of fields, plus where they came from."""
    id: str
    sample: str
    library: str
    platform_unit: str
    platform: str = ILLUMINA_PLATFORM
    description: str = ""
    source: str = ""

    def fastq_to_sam_args(self):
        args = ["-SM", self.sample, "-RG", self.id, "-LB", self.library,
                "-PL", self.platform, "-PU", self.platform_unit]
        if self.description:
            args += ["-DS", self.description]
        return args


@dataclass
class PreprocessResult:
    sample_id: str
    status: str
    final_bam: str = ""
    recal_table: str = ""
    duplicate_metrics: str = ""
    read_group: ReadGroup = None
    read_group_source: str = ""
    steps_run: list = field(default_factory=list)
    steps_skipped: list = field(default_factory=list)
    error: str = ""


class GatkStepError(RuntimeError):
    """A GATK tool exited non-zero. Carries the tool name for the summary row."""

    def __init__(self, tool, returncode):
        super().__init__(f"gatk {tool} exited {returncode}")
        self.tool = tool
        self.returncode = returncode


def resolve_read_group_repair(choice, is_rna):
    """
    Turn --read-group-repair into a boolean.

    "auto" runs AddOrReplaceReadGroups on the RNA path only, after
    SplitNCigarReads, and never on the WGS path.
    """
    if choice == "always":
        return True
    if choice == "never":
        return False
    return bool(is_rna)


def run_gatk(tool, tool_args, temp_dir=None, java_heap=DEFAULT_JAVA_HEAP):
    """
    Run one GATK tool, raising GatkStepError on failure.

    --java-options -Xmx is set explicitly: GATK's default heap is a fraction of
    physical memory, which on a shared node is both unpredictable and usually
    far less than SLURM granted. MarkDuplicates and SortSam on a whole-genome
    BAM will OOM or thrash without it.
    """
    cmd = ["gatk"]
    if java_heap:
        cmd += ["--java-options", f"-Xmx{java_heap}"]
    cmd += [tool] + [str(a) for a in tool_args]
    if temp_dir:
        cmd += ["--TMP_DIR" if tool in PICARD_TOOLS else "--tmp-dir", str(temp_dir)]
    print("Running: " + " ".join(cmd), flush=True)
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise GatkStepError(tool, result.returncode)


def bam_is_complete(bam):
    """
    True if `bam` exists and passes samtools quickcheck.

    quickcheck reads the header and the EOF block only, so it is instant even
    on a 40gb BAM, and it is specifically good at catching the failure mode
    this chain produces: a BAM truncated by a wall-clock kill, which looks
    perfectly fine to an os.path.exists test.
    """
    bam = Path(bam)
    if not bam.is_file():
        return False
    result = subprocess.run(["samtools", "quickcheck", "-v", str(bam)],
                            capture_output=True, text=True)
    if result.returncode == 0:
        return True
    print(f"Incomplete BAM, will redo: {bam} ({result.stdout.strip()})")
    return False


def _reuse_or_clear(output, label, result):
    """
    Resume support: True means "already done, skip".

    A present-but-truncated output is removed rather than left in place, so the
    retry starts clean instead of appending to a corrupt file.
    """
    output = Path(output)
    if bam_is_complete(output):
        print(f"{label}: {output.name} already complete, skipping")
        result.steps_skipped.append(label)
        return True
    if output.is_file():
        output.unlink()
    result.steps_run.append(label)
    return False


def _first_fastq_header(fastq_path):
    """
    The first line of a FASTQ, gzipped or not.

    Reads line 1 specifically rather than scanning for a line starting with
    '@': a quality string may legitimately begin with '@', so a scan can match
    the wrong line on a malformed file.
    """
    fastq_path = Path(fastq_path)
    opener = gzip.open if fastq_path.suffix == ".gz" else open
    with opener(fastq_path, "rt") as fh:
        header = fh.readline().strip()
    if not header.startswith("@"):
        raise ValueError(f"{fastq_path} does not start with a FASTQ header line")
    return header


def parse_read_group(fastq_path, sample_id, description=None):
    """
    Derive an @RG from the first read's header.

    An Illumina CASAVA 1.8+ header is
        @<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y> <read>:<filtered>:<control>:<index>
    giving ID=<flowcell>.<lane> and PU=<flowcell>.<lane>.<index>, the GATK
    convention.

    LB has no header field at all, so it is derived from sample_id.

    ENA and SRA frequently rewrite headers to `@ERR2814677.1 1/1`, which has
    none of these fields. That is not an error -- BQSR only needs the groups to
    be consistent and PL to be set -- so we fall back to one synthetic group per
    sample and record that in `source`, rather than silently emitting a PU built
    from a read serial number.
    """
    header = _first_fastq_header(fastq_path)
    name, _, comment = header.lstrip("@").partition(" ")
    fields = name.split(":")

    if len(fields) >= 7:
        flowcell, lane = fields[2], fields[3]
        comment_fields = comment.split(":") if comment else []
        barcode = comment_fields[3] if len(comment_fields) >= 4 else ""
        return ReadGroup(
                id=f"{flowcell}.{lane}",
                sample=sample_id,
                library=f"{sample_id}-lib",
                platform_unit=f"{flowcell}.{lane}" + (f".{barcode}" if barcode else ""),
                description=description or "",
                source="illumina-header",
                )

    accession = name.split(".")[0] or sample_id
    return ReadGroup(
            id=f"{sample_id}.1",
            sample=sample_id,
            library=f"{sample_id}-lib",
            platform_unit=f"{accession}.1",
            description=description or "",
            source=f"fallback (header {header[:60]!r} is not Illumina-format)",
            )


def fastq_to_ubam(fastq1, fastq2, out_bam, read_group, result, temp_dir=None,
                  java_heap=DEFAULT_JAVA_HEAP):
    """
    FastqToSam. Output is queryname-sorted by Picard default, which is exactly
    what MergeBamAlignment wants, so no SortSam follows it.
    """
    if _reuse_or_clear(out_bam, "FastqToSam", result):
        return out_bam
    run_gatk("FastqToSam",
          ["-F1", fastq1, "-F2", fastq2, "-O", out_bam] + read_group.fastq_to_sam_args(),
          temp_dir=temp_dir, java_heap=java_heap)
    return out_bam


def sort_queryname(in_bam, out_bam, result, temp_dir=None, java_heap=DEFAULT_JAVA_HEAP):
    """
    SortSam -SO queryname on the ALIGNED BAM.

    Both wgseqpipe/align.py and the two-pass STAR step leave coordinate order,
    so this is unavoidable at present. Emitting queryname order from the
    aligners instead would save a full sort of a ~40gb BAM per sample; that is
    noted as open in the status doc and is this chain's call to make.
    """
    if _reuse_or_clear(out_bam, "SortSam(queryname)", result):
        return out_bam
    run_gatk("SortSam", ["-I", in_bam, "-O", out_bam, "-SO", "queryname"],
          temp_dir=temp_dir, java_heap=java_heap)
    return out_bam


def merge_bam_alignment(aligned_bam, unmapped_bam, out_bam, ref_fasta, result,
                        temp_dir=None, java_heap=DEFAULT_JAVA_HEAP):
    """
    MergeBamAlignment: carries the uBAM's read groups onto the alignments.

    This is why no -R read group is passed to bwa mem or STAR -- setting it
    there would be redundant and could conflict with what arrives here.
    """
    if _reuse_or_clear(out_bam, "MergeBamAlignment", result):
        return out_bam
    run_gatk("MergeBamAlignment",
          ["-ALIGNED", aligned_bam, "-UNMAPPED", unmapped_bam, "-O", out_bam,
           "-R", ref_fasta],
          temp_dir=temp_dir, java_heap=java_heap)
    return out_bam


def mark_duplicates(in_bam, out_bam, metrics_file, result, temp_dir=None,
                    java_heap=DEFAULT_JAVA_HEAP, extra_args=None):
    """
    MarkDuplicates with default parameters, the same on both paths.

    extra_args exists for a deliberate per-path override, should one ever be
    wanted.
    """
    if _reuse_or_clear(out_bam, "MarkDuplicates", result):
        return out_bam
    run_gatk("MarkDuplicates",
          ["-I", in_bam, "-O", out_bam, "-M", metrics_file] + list(extra_args or []),
          temp_dir=temp_dir, java_heap=java_heap)
    return out_bam


def sort_coordinate_and_index(in_bam, out_bam, result, temp_dir=None,
                              java_heap=DEFAULT_JAVA_HEAP):
    """
    SortSam -SO coordinate --CREATE_INDEX true.

    The index is not optional: BaseRecalibrator and every downstream caller
    need it.
    """
    if _reuse_or_clear(out_bam, "SortSam(coordinate)", result):
        return out_bam
    run_gatk("SortSam",
          ["-I", in_bam, "-O", out_bam, "-SO", "coordinate", "--CREATE_INDEX", "true"],
          temp_dir=temp_dir, java_heap=java_heap)
    return out_bam


def split_ncigar_reads(in_bam, out_bam, ref_fasta, result, temp_dir=None,
                       java_heap=DEFAULT_JAVA_HEAP):
    """
    SplitNCigarReads: RNA only. Splits reads spanning splice junctions at their
    N operators, which the DNA-oriented callers downstream cannot interpret.

    A GATK engine tool, so it writes its own output index -- no CREATE_INDEX.
    """
    if _reuse_or_clear(out_bam, "SplitNCigarReads", result):
        return out_bam
    run_gatk("SplitNCigarReads", ["-R", ref_fasta, "-I", in_bam, "-O", out_bam],
          temp_dir=temp_dir, java_heap=java_heap)
    return out_bam


def add_or_replace_read_groups(in_bam, out_bam, read_group, result, temp_dir=None,
                               java_heap=DEFAULT_JAVA_HEAP):
    """
    AddOrReplaceReadGroups, run on the RNA path after SplitNCigarReads.

    The complete group is already set at FastqToSam, so this rewrites values
    the BAM carries; it is kept as the RNA path's behaviour and skippable with
    `--read-group-repair never`.

    --CREATE_INDEX true matters here: this is the last BAM before
    BaseRecalibrator on the RNA path.
    """
    if _reuse_or_clear(out_bam, "AddOrReplaceReadGroups", result):
        return out_bam
    args = ["-I", in_bam, "-O", out_bam,
            "--RGSM", read_group.sample, "--RGID", read_group.id,
            "--RGPL", read_group.platform, "--RGLB", read_group.library,
            "--RGPU", read_group.platform_unit,
            "--CREATE_INDEX", "true"]
    if read_group.description:
        args += ["--RGDS", read_group.description]
    run_gatk("AddOrReplaceReadGroups", args, temp_dir=temp_dir, java_heap=java_heap)
    return out_bam


def normalize_known_sites(known_sites):
    """
    Coerce --known-sites into a list of paths.

    A bare string is wrapped rather than iterated: `for vcf in "a.vcf"` would
    otherwise emit one --known-sites per character, which fails in a way that
    takes a while to read.
    """
    if known_sites is None:
        return []
    if isinstance(known_sites, (str, Path)):
        return [known_sites]
    return list(known_sites)


def base_recalibrate(in_bam, out_table, ref_fasta, known_sites, temp_dir=None,
                     java_heap=DEFAULT_JAVA_HEAP):
    """
    BaseRecalibrator over one or more known-sites VCFs.

    GATK accepts --known-sites repeatedly and takes the union of the sites,
    which is how the human best practice pairs dbSNP with a known-indel set,
    and how a SNP and an indel release are paired for any other organism.
    Taking a sequence also makes "no known sites at all" an empty list rather
    than a special case.

    The table is small and cheap to redo, so unlike the BAM steps this always
    runs rather than resuming.
    """
    args = ["-I", in_bam, "-R", ref_fasta, "-O", out_table]
    for vcf in normalize_known_sites(known_sites):
        args += ["--known-sites", vcf]
    run_gatk("BaseRecalibrator", args, temp_dir=temp_dir, java_heap=java_heap)
    return out_table


def apply_bqsr(in_bam, recal_table, out_bam, ref_fasta, result, temp_dir=None,
               java_heap=DEFAULT_JAVA_HEAP):
    if _reuse_or_clear(out_bam, "ApplyBQSR", result):
        return out_bam
    run_gatk("ApplyBQSR",
          ["-R", ref_fasta, "-I", in_bam, "--bqsr-recal-file", recal_table,
           "-O", out_bam],
          temp_dir=temp_dir, java_heap=java_heap)
    return out_bam


def _discard(paths, keep):
    """
    Delete intermediates once the chain has passed them.

    Not optional housekeeping at this scale: the chain produces six or seven
    BAMs per sample and a WGS BAM is ~40gb, so keeping them all is well over
    200gb per sample. --keep-intermediates turns it off for debugging.

    Both index spellings are swept: Picard's CREATE_INDEX writes <stem>.bai,
    GATK engine tools and samtools write <name>.bam.bai.
    """
    if keep:
        return
    for path in paths:
        path = Path(path)
        for candidate in (path, path.with_suffix(".bai"),
                          path.with_name(path.name + ".bai")):
            if candidate.is_file():
                print(f"Removing intermediate: {candidate.name}")
                candidate.unlink()


def run_preprocessing_chain(fastq1, fastq2, aligned_bam, sample_id, ref_fasta,
                            known_sites, work_dir, is_rna, description=None,
                            mark_duplicates_extra_args=None, temp_dir=None,
                            java_heap=DEFAULT_JAVA_HEAP, keep_intermediates=False,
                            read_group_repair="auto"):
    """
    The whole per-sample chain. Returns a PreprocessResult; raises nothing.

    `known_sites` is a sequence of VCF paths, or empty. When empty the chain
    runs everything up to and including the RNA steps and then stops, leaving
    that BAM as final: BaseRecalibrator cannot run without known sites, and the
    mouse substitute for hg38 dbsnp138 -- the Mouse Genomes Project release --
    needs its contigs renamed from Ensembl to RefSeq style first. This keeps the
    rest of the chain testable before that is done, and the result records BQSR
    as skipped so no row silently looks complete.
    """
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(temp_dir) if temp_dir else work_dir / "tmp"
    temp_dir.mkdir(parents=True, exist_ok=True)
    repair = resolve_read_group_repair(read_group_repair, is_rna)
    known_sites = normalize_known_sites(known_sites)

    result = PreprocessResult(sample_id=sample_id, status="FAILED")
    p = work_dir / sample_id
    ubam = Path(f"{p}_unaln.bam")
    qnsort = Path(f"{p}_aln_qnsort.bam")
    merged = Path(f"{p}_merged.bam")
    markdup = Path(f"{p}_markdup.bam")
    metrics = Path(f"{p}_markdup_metrics.txt")
    coordsort = Path(f"{p}_coordsort.bam")
    splitn = Path(f"{p}_splitncigar.bam")
    regrouped = Path(f"{p}_rg.bam")
    recal_table = Path(f"{p}_recal.table")
    bqsr_bam = Path(f"{p}_bqsr.bam")

    try:
        read_group = parse_read_group(fastq1, sample_id, description=description)
        result.read_group = read_group
        result.read_group_source = read_group.source
        print(f"{sample_id}: @RG ID={read_group.id} PU={read_group.platform_unit} "
              f"LB={read_group.library} [{read_group.source}]")

        common = {"temp_dir": temp_dir, "java_heap": java_heap}
        fastq_to_ubam(fastq1, fastq2, ubam, read_group, result, **common)
        sort_queryname(aligned_bam, qnsort, result, **common)
        merge_bam_alignment(qnsort, ubam, merged, ref_fasta, result, **common)
        _discard([ubam, qnsort], keep_intermediates)

        mark_duplicates(merged, markdup, metrics, result,
                        extra_args=mark_duplicates_extra_args, **common)
        _discard([merged], keep_intermediates)
        result.duplicate_metrics = str(metrics)

        current = sort_coordinate_and_index(markdup, coordsort, result, **common)
        _discard([markdup], keep_intermediates)

        if is_rna:
            current = split_ncigar_reads(current, splitn, ref_fasta, result, **common)
            _discard([coordsort], keep_intermediates)
        if repair:
            previous = current
            current = add_or_replace_read_groups(current, regrouped, read_group,
                                                 result, **common)
            _discard([previous], keep_intermediates)
        else:
            result.steps_skipped.append("AddOrReplaceReadGroups")

        if not known_sites:
            print(f"{sample_id}: no known-sites VCF given, stopping before "
                  f"BaseRecalibrator/ApplyBQSR. {Path(current).name} is final.")
            result.steps_skipped += ["BaseRecalibrator", "ApplyBQSR"]
            result.final_bam = str(current)
            result.status = "preprocessed (no BQSR)"
            return result

        result.steps_run.append("BaseRecalibrator")
        base_recalibrate(current, recal_table, ref_fasta, known_sites, **common)
        result.recal_table = str(recal_table)
        apply_bqsr(current, recal_table, bqsr_bam, ref_fasta, result, **common)
        _discard([current], keep_intermediates)

        result.final_bam = str(bqsr_bam)
        result.status = "preprocessed"
        return result
    except GatkStepError as exc:
        result.error = str(exc)
        return result
    except Exception as exc:
        result.error = f"{type(exc).__name__}: {exc}"
        return result


# ----------------------------------------------------------------------------
# Shared command driver.
#
# The GATK chain above is the consolidation that matters, but the scaffolding
# around it -- metadata, sample selection, the summary CSV, the job array --
# would otherwise be duplicated verbatim between wgseqpipe/preprocess.py and
# rnaseqpipe/preprocess.py. It lives here instead, parameterised by PipelineSpec,
# so each caller is a docstring, a spec, and two one-line functions.
# ----------------------------------------------------------------------------

PREPROCESS_SUMMARY_FIELDS = [
    "individualID", "group", "status", "final_bam", "recal_table",
    "duplicate_metrics", "read_group_id", "read_group_source",
    "steps_run", "steps_skipped", "error",
    ]


@dataclass(frozen=True)
class PipelineSpec:
    """Everything that differs between the WGS and RNA preprocessing commands."""
    name: str
    subcommand: tuple
    is_rna: bool
    aligned_suffix: str
    summary_name: str
    job_name: str

    def trimmed_pair(self, sample, qc_output, paired_tokens):
        if self.is_rna:
            return sample.trimmed_rna_pair(qc_output, paired_tokens)
        return sample.trimmed_wgs_pair(qc_output, paired_tokens)


def preprocess_summary_row(sample, result):
    return {
        "individualID": sample.individual_id,
        "group": sample.group,
        "status": result.status,
        "final_bam": result.final_bam,
        "recal_table": result.recal_table,
        "duplicate_metrics": result.duplicate_metrics,
        "read_group_id": result.read_group.id if result.read_group else "",
        "read_group_source": result.read_group_source,
        "steps_run": ";".join(result.steps_run),
        "steps_skipped": ";".join(result.steps_skipped),
        "error": result.error,
        }


def add_preprocess_arguments(parser, spec, ref_help, aligned_help, known_sites_help):
    """
    The full flag set for a preprocess subcommand.

    `ref_help`, `aligned_help` and `known_sites_help` are passed in because they
    are the only text that differs between the two paths; everything else is
    identical by construction rather than by careful copying.
    """
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV (see metadata.py)")
    parser.add_argument("--qc-output", required=True,
                        help="The --output directory `tarva qc` wrote; trimmed reads are "
                             f"read from its trimmed/{'RNA' if spec.is_rna else 'WGS'}/ tree")
    parser.add_argument("--aligned", required=True, help=aligned_help)
    parser.add_argument("--output", required=True,
                        help="Directory for preprocessed BAMs, as "
                             "<output>/<group>/<individualID>/")
    parser.add_argument("--ref", required=True, help=ref_help)
    parser.add_argument("--paired", type=parse_paired_tokens, default=DEFAULT_PAIRED_TOKENS,
                        metavar="FWD,REV",
                        help="Substrings identifying the forward and reverse mate in trimmed "
                             "FASTQ filenames (default: _1,_2)")
    parser.add_argument("--known-sites", action="append", metavar="VCF",
                        help=known_sites_help + " Repeatable: pass it once per VCF "
                             "and GATK takes the union, as the human best practice "
                             "pairs dbSNP with a known-indel set")
    parser.add_argument("--work-dir",
                        help="Where to write the decompressed FASTA if --ref is gzipped "
                             "(default: alongside --ref)")
    parser.add_argument("--temp-dir", metavar="DIR",
                        help="GATK/Picard temp directory (default: <output>/<group>/<id>/tmp). "
                             "Several steps spill as much as the BAM itself, and the JVM "
                             "default of /tmp is often small on a compute node -- point this "
                             "at scratch")
    parser.add_argument("--java-heap", default=DEFAULT_JAVA_HEAP, metavar="SIZE",
                        help=f"GATK --java-options -Xmx (default: {DEFAULT_JAVA_HEAP}). "
                             "Keep this below the job's --mem: the JVM needs headroom "
                             "beyond the heap")
    parser.add_argument("--keep-intermediates", action="store_true",
                        help="Keep every intermediate BAM. Several times the final BAM "
                             "per sample -- debugging only")
    parser.add_argument("--read-group-repair", choices=READ_GROUP_REPAIR_CHOICES,
                        default="auto",
                        help="Whether to run AddOrReplaceReadGroups late in the chain. "
                             "'auto' (default) is on for RNA, off for WGS. 'never' "
                             "skips it -- the complete @RG is already set at "
                             "FastqToSam, so it only rewrites a large BAM. 'always' "
                             "runs it on both paths")
    parser.add_argument("--description", help="@RG DS field, applied to every sample")
    parser.add_argument("--jobs", type=int, default=1,
                        help="Samples to preprocess concurrently within one process "
                             "(default: 1). Each runs its own JVM at --java-heap, so this "
                             "multiplies memory. Prefer --array on a cluster")
    parser.add_argument("--sample-index", type=int, metavar="N",
                        help="Preprocess only sample N (0-based, metadata row order). Set "
                             "automatically from $SLURM_ARRAY_TASK_ID by --array")
    parser.add_argument("--overwrite", action="store_true",
                        help="Reprocess samples that already have a final BAM")
    parser.add_argument("--hpc", action="store_true",
                        help="Submit as a SLURM job instead of running locally")
    parser.add_argument("--array", action="store_true",
                        help="With --hpc, submit a job array of one task per sample. The "
                             "chain is many hours per sample, so this is the intended mode")
    parser.add_argument("--array-concurrency", type=int, metavar="N",
                        help="Cap simultaneously running array tasks (SLURM's %%N suffix)")
    parser.add_argument("--conda-env", default="bio_qc",
                        help="Conda environment activated in the generated SLURM script "
                             "(default: bio_qc -- provides gatk and samtools)")
    parser.add_argument("--modules", default="",
                        help="Comma-separated HPC modules to load. Empty by default, "
                             "since gatk and samtools both come from the conda env and "
                             "a module's version may be older than the env's")
    parser.add_argument("--directives-dict", dest="directives_dict",
                        type=parse_directives_dict, default={}, metavar="DICT",
                        help="SLURM sbatch directives as a Python dict literal, only used "
                             "with --hpc. Budget --mem above --java-heap. Example: "
                             "\"{'partition': 'Orion', 'mem': '48gb', 'time': '24:00:00', "
                             "'cpus_per_task': 4}\"")


def _require_reference_companions(fasta):
    """
    GATK needs a .fai and a .dict beside the FASTA; fail fast if either is absent.

    Both are quick to build, but building them from inside a job array would have
    every task racing to write the same two files. prep-reference already makes
    them, so this points there rather than doing it here -- the same reasoning as
    wgseqpipe align declining to build a missing BWA index.
    """
    fasta = Path(fasta)
    fai = fasta.with_name(fasta.name + ".fai")
    stem = fasta.name
    for suffix in (".fa", ".fasta", ".fna"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    sequence_dict = fasta.with_name(f"{stem}.dict")
    missing = [str(p) for p in (fai, sequence_dict) if not p.is_file()]
    if missing:
        raise SystemExit(
                f"tarva preprocess: {fasta} is missing {missing}. GATK needs both a "
                f"FASTA index and a sequence dictionary. Run "
                f"`tarva prep-reference --ref {fasta}` first -- it builds both and is "
                f"quick, but doing it here would have concurrent array tasks racing to "
                f"write the same files."
                )


def _require_known_sites(known_sites):
    """
    Every known-sites VCF must exist and be indexed before the chain starts.

    BaseRecalibrator is the second-to-last step, so a typo'd path or a missing
    index would otherwise surface only after every BAM step had already run --
    hours per sample, thrown away at the end.

    A .csi index is deliberately not accepted. htsjdk, which GATK reads features
    through, wants .tbi or Tribble's .idx; .csi is what bcftools writes and what
    the Mouse Genomes Project ships, so this is the likely case rather than an
    exotic one.
    """
    problems = []
    for path in known_sites:
        path = Path(path)
        if not path.is_file():
            problems.append(f"{path} does not exist")
            continue
        if any(Path(f"{path}{ext}").is_file() for ext in (".tbi", ".idx")):
            continue
        extra = (" (a .csi is present, but GATK does not read .csi)"
                 if Path(f"{path}.csi").is_file() else "")
        problems.append(f"{path} has no .tbi or .idx beside it{extra}")
    if problems:
        raise SystemExit(
                "tarva preprocess: --known-sites problem(s):\n  "
                + "\n  ".join(problems)
                + "\nIndex a VCF with `gatk IndexFeatureFile -I <vcf>`. Contig names "
                  "must also match --ref: the Mouse Genomes Project release is "
                  "Ensembl-style and needs renaming to RefSeq accessions first."
                )


def _select_samples(meta, sample_index):
    if sample_index is None:
        return meta.samples
    if not 0 <= sample_index < len(meta.samples):
        raise SystemExit(
                f"tarva preprocess: --sample-index {sample_index} is out of range for the "
                f"{len(meta.samples)} sample(s) in --metadata "
                f"(valid: 0..{len(meta.samples) - 1})."
                )
    return [meta.samples[sample_index]]


def _summary_path(output_dir, spec, samples, sample_index):
    """One writer, one file -- see wgseqpipe/align.py for why array tasks can't share."""
    if sample_index is None:
        return Path(output_dir) / spec.summary_name
    return (Path(output_dir) / "summary_parts"
            / f"{sample_index:03d}_{samples[0].individual_id}.csv")


def _merge_summary_parts(output_dir, spec, total_samples):
    """Last array task to finish assembles the canonical summary from the parts."""
    parts = sorted((Path(output_dir) / "summary_parts").glob("*.csv"))
    combined = Path(output_dir) / spec.summary_name
    if len(parts) < total_samples:
        print(f"Summary parts: {len(parts)}/{total_samples}. {combined.name} will be "
              f"written by the last array task to finish.")
        return None
    partial = combined.with_name(combined.name + ".partial")
    with open(partial, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=PREPROCESS_SUMMARY_FIELDS)
        writer.writeheader()
        for part in parts:
            with open(part, newline="") as part_fh:
                for row in csv.DictReader(part_fh):
                    writer.writerow(row)
    partial.replace(combined)
    print(f"Wrote combined summary: {combined}")
    return combined


def _gb(value):
    """Leading number of a SLURM size string ('48gb', '32G', 150) in gb, or None."""
    match = re.match(r"^\s*(\d+)\s*([a-zA-Z]*)", str(value))
    if not match:
        return None
    amount, unit = int(match.group(1)), match.group(2).lower()
    if unit.startswith("m"):
        return amount / 1024
    if unit.startswith("t"):
        return amount * 1024
    return amount


def _warn_heap_vs_mem(args):
    """
    The JVM needs real memory beyond -Xmx -- metaspace, thread stacks, Picard's
    off-heap buffers. A --mem equal to or below --java-heap is an OOM kill
    waiting to happen, and with --jobs > 1 the heap is multiplied.
    """
    mem_gb = _gb(args.directives_dict.get("mem")) if args.directives_dict else None
    heap_gb = _gb(args.java_heap)
    if mem_gb is None or heap_gb is None:
        return
    needed = heap_gb * args.jobs
    if mem_gb <= needed:
        print(f"WARNING: --mem {args.directives_dict.get('mem')} is not above "
              f"--java-heap {args.java_heap}"
              + (f" x --jobs {args.jobs}" if args.jobs > 1 else "")
              + f" (= {needed:g}gb of heap). The JVM needs headroom beyond -Xmx, so "
                f"SLURM will likely OOM-kill this. Raise --mem or lower --java-heap.",
              file=sys.stderr)


def _preprocess_one(sample, spec, args, ref):
    """One sample. Exceptions become FAILED rows, never aborting the pool."""
    sample_dir = Path(args.output) / sample.group / sample.individual_id
    try:
        fastq1, fastq2 = spec.trimmed_pair(sample, args.qc_output, args.paired)
        aligned_bam = (Path(args.aligned) / sample.group / sample.individual_id
                       / f"{sample.individual_id}{spec.aligned_suffix}")
        if not aligned_bam.is_file():
            raise FileNotFoundError(
                    f"aligned BAM not found: {aligned_bam}. Run `tarva "
                    f"{' '.join(spec.subcommand[:-1])} align` first, and check --aligned "
                    f"points at that command's --output."
                    )
        if args.overwrite:
            for stale in sample_dir.glob(f"{sample.individual_id}_*.bam"):
                stale.unlink()

        result = run_preprocessing_chain(
                fastq1=fastq1,
                fastq2=fastq2,
                aligned_bam=aligned_bam,
                sample_id=sample.individual_id,
                ref_fasta=ref,
                known_sites=args.known_sites or [],
                work_dir=sample_dir,
                is_rna=spec.is_rna,
                description=args.description,
                temp_dir=Path(args.temp_dir) / sample.individual_id if args.temp_dir else None,
                java_heap=args.java_heap,
                keep_intermediates=args.keep_intermediates,
                read_group_repair=args.read_group_repair,
                )
        return preprocess_summary_row(sample, result)
    except Exception as exc:
        return preprocess_summary_row(
                sample,
                PreprocessResult(sample_id=sample.individual_id, status="FAILED",
                                 error=f"{type(exc).__name__}: {exc}"),
                )


def run_preprocess_command(args, spec):
    """Local (non---hpc) execution: fan out over samples, flush the summary per sample."""
    meta = load_metadata(args.metadata)
    samples = _select_samples(meta, args.sample_index)
    ref = ensure_unzipped_fasta(args.ref, work_dir=args.work_dir)
    _require_reference_companions(ref)
    _require_known_sites(args.known_sites or [])

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_csv = _summary_path(output_dir, spec, samples, args.sample_index)
    summary_csv.parent.mkdir(parents=True, exist_ok=True)

    print(f"Reference: {ref}")
    print(f"{spec.name.upper()} preprocessing: {len(samples)} of {len(meta.samples)} "
          f"sample(s) ({', '.join(s.individual_id for s in samples)}), "
          f"{args.jobs} at a time, -Xmx{args.java_heap} each")
    if args.known_sites:
        print(f"Known sites ({len(args.known_sites)}): "
              + ", ".join(str(v) for v in args.known_sites))
    else:
        print("No --known-sites: BaseRecalibrator and ApplyBQSR will be skipped and the "
              "coordinate-sorted BAM left as final.")

    rows = []
    with open(summary_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=PREPROCESS_SUMMARY_FIELDS)
        writer.writeheader()
        fh.flush()
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(_preprocess_one, sample, spec, args, ref)
                       for sample in samples]
            for future in as_completed(futures):
                row = future.result()
                writer.writerow(row)
                fh.flush()
                rows.append(row)

    failed = [r for r in rows if r["status"] == "FAILED"]
    done = [r for r in rows if r["status"] != "FAILED"]
    print(f"\n{spec.name.upper()} preprocessing: {len(done)} completed, {len(failed)} failed")
    print(f"Wrote summary: {summary_csv}")
    if args.sample_index is not None:
        _merge_summary_parts(output_dir, spec, len(meta.samples))
    for row in rows:
        if row["read_group_source"].startswith("fallback"):
            print(f"WARNING: {row['individualID']} read group came from the fallback path "
                  f"({row['read_group_source']}) -- one synthetic group for the whole "
                  f"sample, so BQSR cannot model per-lane effects.", file=sys.stderr)
    for row in failed:
        print(f"  FAILED {row['individualID']}: {row['error']}", file=sys.stderr)
    if failed:
        raise SystemExit(
                f"tarva {' '.join(spec.subcommand)}: {len(failed)} sample(s) failed "
                f"-- see {summary_csv}"
                )
    return rows


def submit_preprocess_hpc(args, spec):
    """Generate and submit the sbatch script, as a job array when --array is given."""
    command_parts = ["python", "-u", "-m", "tarva.cli", *spec.subcommand,
                     "--metadata", str(args.metadata),
                     "--qc-output", str(args.qc_output),
                     "--aligned", str(args.aligned),
                     "--output", str(args.output),
                     "--ref", str(args.ref),
                     "--paired", ",".join(args.paired),
                     "--java-heap", str(args.java_heap),
                     "--read-group-repair", str(args.read_group_repair),
                     "--jobs", str(args.jobs)]
    # Repeatable, so forwarded one at a time rather than through the pairs below.
    for vcf in args.known_sites or []:
        command_parts += ["--known-sites", str(vcf)]
    for flag, value in (("--work-dir", args.work_dir),
                        ("--temp-dir", args.temp_dir),
                        ("--description", args.description)):
        if value:
            command_parts += [flag, str(value)]
    if args.keep_intermediates:
        command_parts.append("--keep-intermediates")
    if args.overwrite:
        command_parts.append("--overwrite")

    _warn_heap_vs_mem(args)
    sbatch_kwargs, extra_directives = split_directives(args.directives_dict)
    # GATK here is effectively single-threaded; extra cores mostly buy the JVM
    # room for concurrent GC rather than faster tools.
    sbatch_kwargs.setdefault("cpus_per_task", 4 * args.jobs)

    log_suffix = ""
    if args.array:
        if args.jobs > 1:
            print(f"WARNING: --array already gives each sample its own task, so --jobs "
                  f"{args.jobs} multiplies every task's memory demand. Use --jobs 1 "
                  f"with --array.", file=sys.stderr)
        n_samples = len(load_metadata(args.metadata).samples)
        array_spec = f"0-{n_samples - 1}"
        if args.array_concurrency:
            array_spec += f"%{args.array_concurrency}"
        extra_directives = dict(extra_directives)
        extra_directives.setdefault("array", array_spec)
        command_parts += ["--sample-index", "$SLURM_ARRAY_TASK_ID"]
        log_suffix = "_%a"
        print(f"Submitting a job array of {n_samples} task(s): --array={array_spec}")

    log_dir = Path(args.output) / "logs"
    script_path = write_sbatch_script(
            job_name=spec.job_name,
            commands=" ".join(command_parts),
            log_dir=log_dir,
            script_path=log_dir / f"{spec.job_name}.slurm",
            modules=["anaconda3"] + [m.strip() for m in args.modules.split(",") if m.strip()],
            conda_env=args.conda_env,
            extra_directives=extra_directives,
            log_suffix=log_suffix,
            **sbatch_kwargs,
            )
    print(f"Wrote SLURM script: {script_path}")
    print(submit_job(script_path))


def run_preprocess(args, spec):
    """Entry point both preprocess modules delegate to."""
    if args.hpc:
        submit_preprocess_hpc(args, spec)
        return
    run_preprocess_command(args, spec)
