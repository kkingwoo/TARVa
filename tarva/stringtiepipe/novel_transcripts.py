"""
stringtiepipe novel-transcripts: novel-transcript FASTA, merged with the genome

    `tarva stringtiepipe novel-transcripts --assembly <assembly_dir>
        --gtf <reference_annotation.gff[.gz]> --ref <genome.fa[.gz]>
        --output <dir> [--gffcompare-bin PATH] [--gffread-bin PATH]
        [--known-class-codes CODES] [--overwrite] [--work-dir DIR]
        [--hpc] [--conda-env NAME] [--modules M1,M2] [--directives-dict "{...}"]`

Four stages:

    1.) gffcompare the StringTie merged.gtf against the reference annotation,
        which is what assigns each transcript a class_code
    2.) select transcripts whose class_code is not in --known-class-codes
        (default "="), writing novel.gtf
    3.) gffread -w to extract those transcripts' spliced sequences
    4.) concatenate the genome FASTA and the novel sequences into
        combined_reference.fa

Note what stage 4 implies: novel transcript sequences become additional
contigs alongside the genome they were derived from, so reads from those
regions can map to either copy. It is called out here because it affects
multi-mapping rates downstream.

Input:  the assembly directory holding merged.gtf, the reference annotation
        and the genome FASTA. gffcompare uses the same gclib parser as
        StringTie, so --gtf wants the annotation `stringtiepipe assemble` was
        given -- for NCBI RefSeq the .gff, not the .gtf, which gclib rejects.
Output: gffcompare/, novel.gtf, novel_transcripts.fa and
        combined_reference.fa with its .fai and .dict. That combined FASTA is
        a reference in its own right -- the two-pass RNA alignment maps
        against it -- but its STAR index is NOT built here, so
        `tarva prep-reference --ref <combined> --star-index DIR` is still
        required before alignment.
"""

import csv
import shutil
import subprocess

from collections import Counter
from pathlib import Path

from tarva.reference import (ensure_fasta_index, ensure_sequence_dictionary,
                             ensure_unzipped, ensure_unzipped_fasta)
from tarva.stringtiepipe.assemble import MERGED_GTF_NAME
from tarva.slurm_utils import parse_directives_dict, split_directives, submit_job, write_sbatch_script

DEFAULT_GFFCOMPARE_BIN = "gffcompare"
DEFAULT_GFFREAD_BIN = "gffread"
DEFAULT_KNOWN_CLASS_CODES = ("=",)
GFFCOMPARE_PREFIX = "gffcmp"
NOVEL_GTF_NAME = "novel.gtf"
NOVEL_FASTA_NAME = "novel_transcripts.fa"
COMBINED_FASTA_NAME = "combined_reference.fa"
CLASS_CODE_COUNTS_NAME = "class_code_counts.csv"


def add_arguments(parser):
    parser.add_argument("--assembly", required=True,
                        help="The --output directory `stringtiepipe assemble` wrote; "
                             f"{MERGED_GTF_NAME} is read from its root")
    parser.add_argument("--gtf", required=True,
                        help="Reference annotation gffcompare compares against; decompressed "
                             "first if it ends in .gz. Use the same annotation given to assemble")
    parser.add_argument("--ref", required=True,
                        help="Genome FASTA the novel sequences are extracted from and merged "
                             "with; decompressed first if it ends in .gz")
    parser.add_argument("--output", required=True, help="Directory for all output of this step")
    parser.add_argument("--gffcompare-bin", default=DEFAULT_GFFCOMPARE_BIN,
                        help="gffcompare executable (default: from PATH). Pass an absolute "
                             "path when it lives in a different conda env than --conda-env, "
                             "e.g. ~/.conda/envs/stringtie/bin/gffcompare")
    parser.add_argument("--gffread-bin", default=DEFAULT_GFFREAD_BIN,
                        help="gffread executable (default: from PATH); same note as "
                             "--gffcompare-bin")
    parser.add_argument("--known-class-codes", default=",".join(DEFAULT_KNOWN_CLASS_CODES),
                        metavar="CODES",
                        help="Comma-separated gffcompare class codes treated as KNOWN; every "
                             "other code is novel (default: '='). Widen this to e.g. '=,c,j' "
                             "to also treat contained and junction-sharing transcripts as known")
    parser.add_argument("--work-dir",
                        help="Where to write decompressed copies of --gtf/--ref if they are "
                             "gzipped (default: alongside each source)")
    parser.add_argument("--overwrite", action="store_true",
                        help="Redo stages whose output is already present")
    parser.add_argument("--hpc", action="store_true",
                        help="Submit as a SLURM job instead of running locally")
    parser.add_argument("--conda-env", default="bio_qc",
                        help="Conda environment activated in the generated SLURM script "
                             "(default: bio_qc -- runs the CLI itself)")
    parser.add_argument("--modules", default="",
                        help="Comma-separated HPC modules to load (default: none -- gffcompare "
                             "and gffread come from their --*-bin paths)")
    parser.add_argument("--directives-dict", dest="directives_dict",
                        type=parse_directives_dict, default={}, metavar="DICT",
                        help="SLURM sbatch directives as a Python dict literal, only used with "
                             "--hpc. This step is light: \"{'partition': 'Orion', "
                             "'mem': '16gb', 'time': '02:00:00', 'cpus_per_task': 2}\"")


def _attribute(attributes, key):
    """Pull one GTF attribute value out of the 9th column."""
    marker = f'{key} "'
    start = attributes.find(marker)
    if start == -1:
        return None
    start += len(marker)
    end = attributes.find('"', start)
    return attributes[start:end] if end != -1 else None


def run_gffcompare(gffcompare_bin, reference_gtf, merged_gtf, output_dir, overwrite=False):
    """
    Stage 1: annotate merged.gtf against the reference, yielding class codes.

    gffcompare is what assigns class_code; `stringtie --merge` does not emit it,
    which is why this stage exists rather than reading merged.gtf directly.
    """
    compare_dir = Path(output_dir) / "gffcompare"
    annotated = compare_dir / f"{GFFCOMPARE_PREFIX}.annotated.gtf"
    if annotated.is_file() and not overwrite:
        print(f"gffcompare output already present: {annotated}")
        return annotated
    compare_dir.mkdir(parents=True, exist_ok=True)
    cmd = [str(gffcompare_bin), "-r", str(reference_gtf),
           "-o", str(compare_dir / GFFCOMPARE_PREFIX), str(merged_gtf)]
    print("Running gffcompare: " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    if not annotated.is_file():
        raise SystemExit(
                f"novel-transcripts: gffcompare exited 0 but {annotated} is missing"
                )
    return annotated


def select_novel_transcripts(annotated_gtf, novel_gtf, known_class_codes):
    """
    Stage 2: write the subset of `annotated_gtf` whose class_code is novel.

    Two passes on purpose. gffcompare carries class_code on transcript lines,
    but gffread needs the exon lines too, so pass one maps transcript_id ->
    class_code and pass two emits every line belonging to a selected
    transcript. Filtering line-by-line on class_code alone would drop the exons
    and produce an unusable GTF.
    """
    annotated_gtf, novel_gtf = Path(annotated_gtf), Path(novel_gtf)
    class_codes, counts = {}, Counter()
    with open(annotated_gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            transcript_id = _attribute(fields[8], "transcript_id")
            code = _attribute(fields[8], "class_code")
            if transcript_id and code and transcript_id not in class_codes:
                class_codes[transcript_id] = code
                counts[code] += 1

    novel_ids = {tid for tid, code in class_codes.items() if code not in known_class_codes}
    if not novel_ids:
        raise SystemExit(
                f"novel-transcripts: no transcript in {annotated_gtf} has a class_code outside "
                f"{sorted(known_class_codes)}, so there is nothing novel to extract. Class "
                f"codes seen: {dict(counts)}"
                )

    written = 0
    with open(annotated_gtf) as src, open(novel_gtf, "w") as dst:
        for line in src:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                continue
            if _attribute(fields[8], "transcript_id") in novel_ids:
                dst.write(line)
                written += 1

    print(f"Selected {len(novel_ids)} novel of {len(class_codes)} transcripts "
          f"({written} GTF lines) -> {novel_gtf}")
    return novel_ids, counts


def write_class_code_counts(counts, known_class_codes, path):
    """Per-class-code tally, so the novelty call is auditable after the fact."""
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["class_code", "transcripts", "treated_as"])
        for code, count in sorted(counts.items(), key=lambda kv: -kv[1]):
            writer.writerow([code, count, "known" if code in known_class_codes else "novel"])
    print(f"Wrote class-code tally: {path}")


def run_gffread(gffread_bin, genome_fasta, novel_gtf, novel_fasta, overwrite=False):
    """Stage 3: extract spliced sequences for the selected transcripts."""
    novel_fasta = Path(novel_fasta)
    if novel_fasta.is_file() and novel_fasta.stat().st_size > 0 and not overwrite:
        print(f"Novel transcript FASTA already present: {novel_fasta}")
        return novel_fasta
    cmd = [str(gffread_bin), "-w", str(novel_fasta), "-g", str(genome_fasta), str(novel_gtf)]
    print("Running gffread: " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    if not novel_fasta.is_file() or novel_fasta.stat().st_size == 0:
        raise SystemExit(
                f"novel-transcripts: gffread exited 0 but {novel_fasta} is missing or empty"
                )
    return novel_fasta


def _last_byte(path):
    with open(path, "rb") as fh:
        if fh.seek(0, 2) == 0:
            return b""
        fh.seek(-1, 2)
        return fh.read(1)


def combine_fastas(genome_fasta, novel_fasta, combined, overwrite=False):
    """
    Stage 4: genome + novel sequences in one FASTA.

    A missing trailing newline on either source would concatenate the next
    file's '>' header onto the previous sequence line and silently corrupt both
    records, so each part is newline-terminated explicitly. Streams to .partial
    and renames, for the same reason ensure_unzipped does.
    """
    combined = Path(combined)
    if combined.is_file() and combined.stat().st_size > 0 and not overwrite:
        print(f"Combined reference already present: {combined}")
        return combined
    partial = combined.with_name(combined.name + ".partial")
    print(f"Combining {genome_fasta} + {novel_fasta} -> {combined}")
    with open(partial, "wb") as dst:
        for source in (genome_fasta, novel_fasta):
            with open(source, "rb") as src:
                shutil.copyfileobj(src, dst, length=16 * 1024 * 1024)
            if _last_byte(source) != b"\n":
                dst.write(b"\n")
    partial.replace(combined)
    return combined


def run(args):
    if args.hpc:
        _submit_hpc(args)
        return

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    merged_gtf = Path(args.assembly) / MERGED_GTF_NAME
    if not merged_gtf.is_file():
        raise SystemExit(
                f"novel-transcripts: {merged_gtf} not found. Run `stringtiepipe assemble` "
                f"(at least --steps round1,merge) first."
                )

    known_class_codes = tuple(
            c.strip() for c in str(args.known_class_codes).split(",") if c.strip()
            )
    if not known_class_codes:
        raise SystemExit("novel-transcripts: --known-class-codes was empty")

    reference_gtf = ensure_unzipped(args.gtf, work_dir=args.work_dir,
                                    description="annotation GTF")
    genome_fasta = ensure_unzipped_fasta(args.ref, work_dir=args.work_dir)

    annotated = run_gffcompare(args.gffcompare_bin, reference_gtf, merged_gtf,
                               output_dir, overwrite=args.overwrite)
    novel_gtf = output_dir / NOVEL_GTF_NAME
    _, counts = select_novel_transcripts(annotated, novel_gtf, known_class_codes)
    write_class_code_counts(counts, known_class_codes, output_dir / CLASS_CODE_COUNTS_NAME)

    novel_fasta = run_gffread(args.gffread_bin, genome_fasta, novel_gtf,
                              output_dir / NOVEL_FASTA_NAME, overwrite=args.overwrite)
    combined = combine_fastas(genome_fasta, novel_fasta,
                              output_dir / COMBINED_FASTA_NAME, overwrite=args.overwrite)

    # The .fai and .dict take seconds and every GATK step downstream refuses to
    # start without them, so build them here rather than leaving a manual step
    # that only surfaces once `rnaseqpipe preprocess` has already been queued.
    # Safe inline: this step is a single process, not a job array, so nothing
    # races to write the same file.
    ensure_fasta_index(combined)
    ensure_sequence_dictionary(combined)

    print(f"\nCombined reference: {combined}")
    print("  .fai and .dict written alongside it.")
    print("This is a new reference and still needs its own STAR index before alignment:")
    print(f"  tarva prep-reference --ref {combined} --star-index <dir> --threads N --hpc ...")


def _submit_hpc(args):
    command_parts = [
        # -u: unbuffered, so a job killed at the wall clock still leaves its
        # progress in the .o rather than discarding the buffer.
        "python", "-u", "-m", "tarva.cli", "stringtiepipe", "novel-transcripts",
        "--assembly", str(args.assembly),
        "--gtf", str(args.gtf),
        "--ref", str(args.ref),
        "--output", str(args.output),
        "--gffcompare-bin", str(args.gffcompare_bin),
        "--gffread-bin", str(args.gffread_bin),
        "--known-class-codes", str(args.known_class_codes),
        ]
    if args.work_dir:
        command_parts += ["--work-dir", str(args.work_dir)]
    if args.overwrite:
        command_parts.append("--overwrite")

    sbatch_kwargs, extra_directives = split_directives(args.directives_dict)
    log_dir = Path(args.output) / "logs"
    script_path = write_sbatch_script(
            job_name="tarva_novel_transcripts",
            commands=" ".join(command_parts),
            log_dir=log_dir,
            script_path=log_dir / "tarva_novel_transcripts.slurm",
            modules=["anaconda3"] + [m.strip() for m in args.modules.split(",") if m.strip()],
            conda_env=args.conda_env,
            extra_directives=extra_directives,
            **sbatch_kwargs,
            )
    print(f"Wrote SLURM script: {script_path}")
    print(submit_job(script_path))
