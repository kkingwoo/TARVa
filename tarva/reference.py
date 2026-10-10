"""
prep-reference: decompress, index, and build aligner indexes for a reference

    `tarva prep-reference --ref <genome.fa[.gz]> [--star-index DIR] [--gtf FILE]
        [--bwa-index] [--threads N] [--work-dir DIR] [--force-rebuild]
        [--hpc] [--conda-env NAME] [--modules M1,M2] [--directives-dict "{...}"]`

STAR indexes are built WITHOUT the annotation GTF by default. The annotation
is inserted at mapping time instead (--sjdbGTFfile/--sjdbOverhang), which
keeps one index usable across annotations -- including the StringTie-merged
GTF produced later in the pipeline, which would otherwise force a rebuild.
"""

import gzip
import math
import os
import shutil
import subprocess
import sys

from pathlib import Path

from tarva.slurm_utils import parse_directives_dict, split_directives, submit_job, write_sbatch_script

DEFAULT_STAR_MODULE, DEFAULT_BWA_MODULE = "star/2.7.11b", "bwa/0.7.17"
FASTA_SUFFIXES = (".fa", ".fasta", ".fna")
STAR_INDEX_SENTINELS = ("SAindex", "Genome", "genomeParameters.txt")


def _run(cmd, description):
    print(f"{description}: " + " ".join(str(c) for c in cmd))
    subprocess.run([str(c) for c in cmd], check=True)


def ensure_unzipped(path, work_dir=None, description="file"):
    """
    Return a path to an uncompressed copy of `path`, decompressing if needed.

    """
    path = Path(path)
    if not path.is_file():
        raise SystemExit(f"tarva: {description} not found: {path}")
    if path.suffix != ".gz":
        return path

    dest_dir = Path(work_dir) if work_dir else path.parent
    dest_dir.mkdir(parents=True, exist_ok=True)
    target = dest_dir / path.stem
    if target.is_file() and target.stat().st_size > 0:
        print(f"Already decompressed: {target}")
        return target

    
    partial = target.with_name(f"{target.name}.{os.getpid()}.partial")
    print(f"Decompressing {path} -> {target}")
    with gzip.open(path, "rb") as src, open(partial, "wb") as dst:
        shutil.copyfileobj(src, dst, length=16 * 1024 * 1024)
    partial.replace(target)
    return target


def ensure_unzipped_fasta(ref, work_dir=None):
    """."""
    return ensure_unzipped(ref, work_dir=work_dir, description="reference FASTA")


def ensure_fasta_index(fasta):
    """samtools faidx"""
    fasta = Path(fasta)
    fai = fasta.with_name(fasta.name + ".fai")
    if fai.is_file():
        return fai
    _run(["samtools", "faidx", fasta], "Indexing FASTA")
    return fai


def ensure_sequence_dictionary(fasta):
    """
    gatk CreateSequenceDictionary.

    """
    fasta = Path(fasta)
    stem = fasta.name
    for suffix in FASTA_SUFFIXES:
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    sequence_dict = fasta.with_name(f"{stem}.dict")
    if sequence_dict.is_file():
        return sequence_dict
    _run(["gatk", "CreateSequenceDictionary", "-R", fasta, "-O", sequence_dict],
         "Creating sequence dictionary")
    return sequence_dict


def ensure_bwa_index(fasta):
    """bwa index -- returns the prefix bwa mem expects (the FASTA path itself)."""
    fasta = Path(fasta)
    if all(fasta.with_name(fasta.name + ext).is_file()
           for ext in (".bwt", ".pac", ".ann", ".amb", ".sa")):
        print(f"BWA index already present: {fasta}")
        return fasta
    _run(["bwa", "index", fasta], "Building BWA index")
    return fasta


def genome_length(fasta_index):
    """Total reference length, summed from column 2 of the .fai."""
    total = 0
    with open(fasta_index) as fh:
        for line in fh:
            fields = line.split("\t")
            if len(fields) > 1:
                total += int(fields[1])
    if total <= 0:
        raise SystemExit(
                f"tarva prep-reference: {fasta_index} gives a total sequence length of "
                f"{total} -- the reference is empty or its index is malformed. Delete the "
                f".fai and re-run so samtools faidx rebuilds it."
                )
    return total


def sa_index_nbases(total_length):
    """
    STAR's --genomeSAindexNbases, per the manual's min(14, log2(L)/2 - 1).

    STAR fails in confusing ways on small genomes.
    """
    return max(4, min(14, round(math.log2(total_length) / 2 - 1)))


def ensure_star_index(fasta, index_dir, threads=1, gtf=None, sjdb_overhang=None,
                      sa_index_nbases_override=None, force_rebuild=False):
    """
    STAR --runMode genomeGenerate.

    """
    fasta, index_dir = Path(fasta), Path(index_dir)
    existing = [s for s in STAR_INDEX_SENTINELS if (index_dir / s).is_file()]
    if len(existing) == len(STAR_INDEX_SENTINELS) and not force_rebuild:
        version = _star_index_version(index_dir)
        print(f"STAR index already present: {index_dir}"
              + (f" (versionGenome {version})" if version else ""))
        print("  If STAR later refuses this index as incompatible, re-run with --force-rebuild.")
        return index_dir
    if existing and not force_rebuild:
        raise SystemExit(
                f"tarva prep-reference: {index_dir} looks like a partially built STAR index "
                f"(found {existing}, expected all of {list(STAR_INDEX_SENTINELS)}) -- a "
                f"previous genomeGenerate was probably killed. Re-run with --force-rebuild "
                f"to discard it and start over."
                )
    if force_rebuild and index_dir.is_dir():
        print(f"--force-rebuild: removing {index_dir}")
        shutil.rmtree(index_dir)

    fai = ensure_fasta_index(fasta)
    nbases = sa_index_nbases_override or sa_index_nbases(genome_length(fai))
    index_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "STAR", "--runMode", "genomeGenerate",
        "--genomeDir", index_dir,
        "--genomeFastaFiles", fasta,
        "--runThreadN", threads,
        "--genomeSAindexNbases", nbases,
        "--genomeLoad", "NoSharedMemory",
        ]
    if gtf:
        cmd += ["--sjdbGTFfile", gtf]
        if sjdb_overhang:
            cmd += ["--sjdbOverhang", sjdb_overhang]
    _run(cmd, "Building STAR index")
    return index_dir


def _star_index_version(index_dir):
    params = Path(index_dir) / "genomeParameters.txt"
    if not params.is_file():
        return None
    try:
        with open(params) as fh:
            for line in fh:
                if line.startswith("versionGenome"):
                    return line.split(None, 1)[1].strip()
    except OSError:
        return None
    return None


def add_arguments(parser):
    parser.add_argument("--ref", required=True,
                        help="Reference genome FASTA; decompressed first if it ends in .gz")
    parser.add_argument("--work-dir",
                        help="Where to write the decompressed FASTA if --ref is gzipped "
                             "(default: alongside --ref; set this if ref directory is "
                             "read-only or short on space")
    parser.add_argument("--star-index", metavar="DIR",
                        help="Build (or verify) a STAR genome index in DIR")
    parser.add_argument("--bwa-index", action="store_true",
                        help="Build (or verify) the BWA index alongside the FASTA")
    parser.add_argument("--gtf",
                        help="Normally omitted -- "
                             "the annotation is inserted at mapping time so one index serves "
                             "multiple annotations, including the StringTie-merged GTF")
    parser.add_argument("--sjdb-overhang", type=int,
                        help="Only meaningful with --gtf; max read length minus 1")
    parser.add_argument("--genome-sa-index-nbases", type=int,
                        help="Override STAR --genomeSAindexNbases (default: derived from "
                             "genome length per STAR's min(14, log2(L)/2 - 1))")
    parser.add_argument("--threads", type=int, default=1,
                        help="Threads for STAR genomeGenerate (default: 1)")
    parser.add_argument("--force-rebuild", action="store_true",
                        help="Discard and rebuild the STAR index")
    parser.add_argument("--hpc", action="store_true",
                        help="Submit as a SLURM job instead of running locally")
    parser.add_argument("--conda-env", default="bio_qc",
                        help="Conda environment to activate in the generated SLURM script "
                             "(default: bio_qc -- provides samtools and gatk)")
    parser.add_argument("--modules", default=f"{DEFAULT_STAR_MODULE},{DEFAULT_BWA_MODULE}",
                        help="Comma-separated HPC modules to load before running "
                             f"(default: {DEFAULT_STAR_MODULE},{DEFAULT_BWA_MODULE})")
    parser.add_argument("--directives-dict", dest="directives_dict",
                        type=parse_directives_dict, default={}, metavar="DICT",
                        help="SLURM sbatch directives as a Python dict literal, only used "
                             "with --hpc. Example: \"{'partition': "
                             "'Mamba', 'mem': '64gb', 'time': '04:00:00', 'cpus_per_task': 8}\"")


def run(args):
    if args.hpc:
        _submit_hpc(args)
        return
    if not args.star_index and not args.bwa_index:
        print("Neither --star-index nor --bwa-index given: preparing FASTA, .fai and .dict only",
              file=sys.stderr)

    fasta = ensure_unzipped_fasta(args.ref, work_dir=args.work_dir)
    gtf = ensure_unzipped(args.gtf, work_dir=args.work_dir,
                          description="annotation GTF") if args.gtf else None
    ensure_fasta_index(fasta)
    ensure_sequence_dictionary(fasta)
    if args.bwa_index:
        ensure_bwa_index(fasta)
    if args.star_index:
        ensure_star_index(
                fasta, args.star_index,
                threads=args.threads,
                gtf=gtf,
                sjdb_overhang=args.sjdb_overhang,
                sa_index_nbases_override=args.genome_sa_index_nbases,
                force_rebuild=args.force_rebuild,
                )
    print(f"Reference ready: {fasta}")


def _submit_hpc(args):
    # -u: unbuffered, so a job killed at the wall clock still leaves its
    # progress in the .o rather than discarding the buffer.
    command_parts = ["python", "-u", "-m", "tarva.cli", "prep-reference",
                     "--ref", str(args.ref), "--threads", str(args.threads)]
    for flag, value in (("--work-dir", args.work_dir), ("--star-index", args.star_index),
                        ("--gtf", args.gtf), ("--sjdb-overhang", args.sjdb_overhang),
                        ("--genome-sa-index-nbases", args.genome_sa_index_nbases)):
        if value is not None:
            command_parts += [flag, str(value)]
    if args.bwa_index:
        command_parts.append("--bwa-index")
    if args.force_rebuild:
        command_parts.append("--force-rebuild")

    sbatch_kwargs, extra_directives = split_directives(args.directives_dict)
    sbatch_kwargs.setdefault("cpus_per_task", args.threads)
    log_dir = Path(args.star_index or Path(args.ref).parent) / "logs"
    script_path = write_sbatch_script(
            job_name="tarva_prep_reference",
            commands=" ".join(command_parts),
            log_dir=log_dir,
            script_path=log_dir / "tarva_prep_reference.slurm",
            modules=["anaconda3"] + [m.strip() for m in args.modules.split(",") if m.strip()],
            conda_env=args.conda_env,
            extra_directives=extra_directives,
            **sbatch_kwargs,
            )
    print(f"Wrote SLURM script: {script_path}")
    print(submit_job(script_path))
