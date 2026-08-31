"""
sort-input: copy raw FASTQ files into the <WGS|RNA>/<group>/<adapter>/ tree
    that `tarva qc` expects.

    `tarva sort-input --metadata <meta.csv> --input <input_dir> --output <sorted_dir>`

    For each sample and sequencing type (WGS, RNA), files whose name starts
    with that sample's WGS/RNA prefix are located anywhere under --input
    Original files will be left in their path and copied-only
    into:
        <output>/WGS/<group>/<adapter_stem>/
        <output>/RNA/<group>/<adapter_stem>/
    where <adapter_stem> is the filename (without .fa) of the adapter FASTA
    resolved from that sample's wgs_barcode/rna_barcode metadata cell.

    Files are copied concurrently across --threads workers
    
    This command's --output is --input to `tarva qc`.
"""

import shutil
import sys

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from tarva.metadata import load_metadata
from tarva.adapters import resolve_adapter

FASTQ_EXTENSIONS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def add_arguments(parser):
    parser.add_argument("--metadata", required=True, help="Path to metadata CSV (see metadata.py)")
    parser.add_argument(
            "--input", required=True,
            help="Directory containing raw FASTQ files in whatever layout they exist today (searched recursively)",
            )
    parser.add_argument(
            "--output", required=True,
            help="Destination root for the sorted WGS|RNA/group/adapter tree -- pass this as --input to `tarva qc`",
            )
    parser.add_argument(
            "--threads", type=int, default=1,
            help="Number of files to copy concurrently (default: 1)",
            )


def find_raw_files(raw_root, pattern):
    pattern = str(pattern).strip()
    if not pattern:
        raise ValueError("Empty WGS/RNA prefix -- cannot match any files")
    matches = Path(raw_root).rglob(f"{pattern}*")
    return sorted(p for p in matches if p.is_file() and any(p.name.endswith(ext) for ext in FASTQ_EXTENSIONS))

def _copy_one(job):
    src,dest = job
    shutil.copy2(src, dest)
    return src, dest


def run(args):
    meta = load_metadata(args.metadata)
    raw_root = Path(args.input)
    output_root = Path(args.output)

    copy_jobs = []
    for sample in meta.samples:
        for seq_type, pattern, barcode in (
                ("WGS", sample.wgs_pattern, sample.wgs_barcode),
                ("RNA", sample.rna_pattern, sample.rna_barcode),
            ):
            files = find_raw_files(raw_root, pattern)
            if not files:
                print(
                        f"WARNING: no {seq_type} files found for sample {sample.individual_id} "
                        f"(pattern: {pattern!r}, searched under {raw_root})",
                        file=sys.stderr,
                        )
                continue

            adapter = resolve_adapter(barcode)
            dest_dir = output_root / seq_type / sample.group / adapter.stem
            dest_dir.mkdir(parents=True, exist_ok=True)

            for f in files:
                copy_jobs.append((f, dest_dir / f.name))
    copied = 0
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = [pool.submit(_copy_one, job) for job in copy_jobs]
        for future in as_completed(futures):
            src, dest = future.result()
            copied +=1
            print(f"Copied {src} -> {dest}")

    print(f"Done -- copied {copied} file(s) into {output_root}")
