"""
Parses and validates the input metadata CSV, and resolves each sample's FASTQ
paths before and after trimming.

Input:  the metadata CSV described below, plus the directory trees
        `tarva sort-input` and `tarva qc` write.
Output: a Metadata object holding one Sample per row, with the per-sample file
        lookups the pipeline steps call.

The CSV must have the following required columns:
    group        - contains exactly TWO distinct condition or tissue values
                   Example: tumor, peritumoral
    individualID - unique identifier for each individual
                   Example: LungF2
    WGS          - FASTQ filename prefix in all WGS files that belong to individualID
                   Example: files = [ERR000032_1.fastq.gz, ERR000032_2.fastq.gz]
                            WGS value = ERR000032
    RNA          - FASTQ filename prefix in all RNA-seq files that belong to individual
                   Example: file = SQXR24385.fq
                            RNA value = SQXR24385
    wgs_barcode  - adapter for this sample's WGS reads: a bundled adapter
                   number/name (see `tarva qc --list-adapters`) or a path to
                   a custom adapter FASTA. Also determines the folder this
                   sample's WGS files are expected to live in (see below).
    rna_barcode  - same, for RNA-seq reads.

Input FASTQ files (produced by `tarva sort-input`) are expected to live in a
directory tree of the form:
    <input>/WGS/<group>/<adapter_stem>/*.fastq.gz
    <input>/RNA/<group>/<adapter_stem>/*.fastq.gz
where <adapter_stem> is the filename (without .fa) of the adapter FASTA
resolved from wgs_barcode/rna_barcode.

After `tarva qc`, trimmed reads live in a parallel tree rooted at that
command's --output:
    <qc_output>/trimmed/WGS/<group>/<adapter_stem>/*.paired.fastq.gz
    <qc_output>/trimmed/RNA/<group>/<adapter_stem>/*.paired.fastq.gz
Resolve these with Sample.trimmed_*_files()/trimmed_*_pair(), not the
pre-trim wgs_files()/rna_files() -- the roots differ (`qc`'s --output vs
`sort-input`'s) and the trimmed glob must exclude the unpaired mates.
"""

import argparse

from dataclasses import dataclass, field
from pathlib import Path
import pandas as pd

from tarva.adapters import resolve_adapter

REQUIRED_COLUMNS = ["group", "individualID", "WGS", "RNA", "wgs_barcode", "rna_barcode"]
TRIMMED_SUFFIX = ".paired.fastq.gz"
DEFAULT_PAIRED_TOKENS = ("_1", "_2")


def parse_paired_tokens(value):
    """
    argparse `type=` callable for --paired: "FWD,REV" -> ("FWD", "REV").

    These are the substrings that distinguish the forward from the reverse mate
    in a FASTQ filename -- "_1,_2" for ERR2814681_1.fastq.gz, "_R1_001,_R2_001",
    "END1,END2". A flag rather than a constant, because the naming varies
    between datasets.
    """
    parts = [p.strip() for p in str(value).split(",")]
    if len(parts) != 2 or not all(parts):
        raise argparse.ArgumentTypeError(
                f"--paired must be two non-empty comma-separated tokens, "
                f"e.g. '_1,_2' (got: {value!r})"
                )
    if parts[0] == parts[1]:
        raise argparse.ArgumentTypeError(
                f"--paired forward and reverse tokens must differ (both were {parts[0]!r})"
                )
    return (parts[0], parts[1])


def _match_files(pattern, directory):
    directory = Path(directory)
    pattern = str(pattern).strip()
    if not pattern:
        raise ValueError("Empty WGS/RNA prefix -- cannot match any files")
    return sorted(directory.glob(f"{pattern}*"))


def _seq_dir(input_root, seq_type, group, adapter_stem):
    return Path(input_root) / seq_type / group / adapter_stem


def _trimmed_dir(qc_output, seq_type, group, adapter_stem):
    return Path(qc_output) / "trimmed" / seq_type / group / adapter_stem


def _match_trimmed(pattern, directory):
    """
    Trimmomatic's paired output for one sample/seq-type.

    The `*.paired.fastq.gz` glob both excludes the unpaired mates and stays
    agnostic about the `.fastq.fastq` doubling in the current output names, so
    fixing that doesn't invalidate an existing trimmed tree.
    """
    directory = Path(directory)
    pattern = str(pattern).strip()
    if not pattern:
        raise ValueError("Empty WGS/RNA prefix -- cannot match any files")
    return sorted(directory.glob(f"{pattern}*{TRIMMED_SUFFIX}"))


def _classify_pair(files, pattern, paired_tokens, description):
    """
    Split `files` into (forward, reverse) using `paired_tokens`.

    Tokens are matched against the part of the filename following the metadata
    prefix, so a prefix that itself ends in a token-like string can't be
    mistaken for a mate marker. Fails loudly rather than falling back on sort
    order, which silently mispairs any scheme where the reverse mate sorts first.
    """
    fwd_token, rev_token = paired_tokens
    forward, reverse, unmatched = [], [], []
    for f in files:
        remainder = f.name[len(pattern):] if f.name.startswith(pattern) else f.name
        is_fwd, is_rev = fwd_token in remainder, rev_token in remainder
        if is_fwd and is_rev:
            raise ValueError(
                    f"{description}: {f.name} contains both --paired tokens "
                    f"({fwd_token!r} and {rev_token!r}) after the prefix {pattern!r} -- "
                    f"cannot tell which mate it is. Choose tokens that appear in only "
                    f"one mate's filename."
                    )
        if is_fwd:
            forward.append(f)
        elif is_rev:
            reverse.append(f)
        else:
            unmatched.append(f)

    if len(forward) != 1 or len(reverse) != 1:
        raise ValueError(
                f"{description}: expected exactly one forward and one reverse FASTQ, found "
                f"{len(forward)} matching {fwd_token!r} ({[f.name for f in forward]}) and "
                f"{len(reverse)} matching {rev_token!r} ({[f.name for f in reverse]})"
                + (f"; unmatched: {[f.name for f in unmatched]}" if unmatched else "")
                + ". Adjust --paired to match this dataset's mate naming."
                )
    return forward[0], reverse[0]


@dataclass
class Sample:
    individual_id: str
    group: str
    wgs_pattern: str
    rna_pattern: str
    wgs_barcode: str
    rna_barcode: str
    extra: dict = field(default_factory=dict)

    def wgs_dir(self, input_root):
        return _seq_dir(input_root, "WGS", self.group, resolve_adapter(self.wgs_barcode).stem)

    def rna_dir(self, input_root):
        return _seq_dir(input_root, "RNA", self.group, resolve_adapter(self.rna_barcode).stem)

    def wgs_files(self, input_root):
        return _match_files(self.wgs_pattern, self.wgs_dir(input_root))

    def rna_files(self, input_root):
        return _match_files(self.rna_pattern, self.rna_dir(input_root))

    def trimmed_wgs_dir(self, qc_output):
        return _trimmed_dir(qc_output, "WGS", self.group, resolve_adapter(self.wgs_barcode).stem)

    def trimmed_rna_dir(self, qc_output):
        return _trimmed_dir(qc_output, "RNA", self.group, resolve_adapter(self.rna_barcode).stem)

    def trimmed_wgs_files(self, qc_output):
        return _match_trimmed(self.wgs_pattern, self.trimmed_wgs_dir(qc_output))

    def trimmed_rna_files(self, qc_output):
        return _match_trimmed(self.rna_pattern, self.trimmed_rna_dir(qc_output))

    def trimmed_wgs_pair(self, qc_output, paired_tokens=DEFAULT_PAIRED_TOKENS):
        return _classify_pair(
                self.trimmed_wgs_files(qc_output), self.wgs_pattern, paired_tokens,
                f"sample {self.individual_id} (WGS, {self.trimmed_wgs_dir(qc_output)})",
                )

    def trimmed_rna_pair(self, qc_output, paired_tokens=DEFAULT_PAIRED_TOKENS):
        return _classify_pair(
                self.trimmed_rna_files(qc_output), self.rna_pattern, paired_tokens,
                f"sample {self.individual_id} (RNA, {self.trimmed_rna_dir(qc_output)})",
                )


@dataclass
class Metadata:
    samples: list
    groups: list
    df: pd.DataFrame = None

    def samples_in_group(self, group):
        return [s for s in self.samples if s.group == group]


def load_metadata(path):
    """
    Reads and validates the metadata CSV at `path`.

    Raises ValueError if:
        1.) required columns are missing
        2.) the 'group' column does not contain exactly 2 distinct values
        3.) any row is missing a required field
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Metadata file not found: {path}")

    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip() for c in df.columns]

    missing_cols = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_cols:
        raise ValueError(
                f"Metadata file {path} is missing required column(s): {missing_cols}. "
                f"Required columns are: {REQUIRED_COLUMNS}"
                )

    for col in df.columns:
        df[col] = df[col].fillna("").str.strip()

    missing_mask = df[REQUIRED_COLUMNS].isna() | (df[REQUIRED_COLUMNS] == "")
    if missing_mask.any().any():
        bad_rows = df.index[missing_mask.any(axis=1)] + 2
        bad_cols = [c for c in REQUIRED_COLUMNS if missing_mask[c].any()]
        raise ValueError(
                f"Metadata file {path}: missing required value(s) in column(s) {bad_cols} "
                f"at line(s) {list(bad_rows)}"
                )

    groups_seen = list(dict.fromkeys(df["group"]))
    if len(groups_seen) != 2:
        raise ValueError(
                f"Metadata file {path}: the 'group' column must contain exactly 2 "
                f"distinct values, found {len(groups_seen)}: {groups_seen}"
                )

    extra_cols = [c for c in df.columns if c not in REQUIRED_COLUMNS]
    samples = [
        Sample(
            individual_id=row.individualID,
            group=row.group,
            wgs_pattern=row.WGS,
            rna_pattern=row.RNA,
            wgs_barcode=row.wgs_barcode,
            rna_barcode=row.rna_barcode,
            extra={c: getattr(row, c) for c in extra_cols},
            )
        for row in df.itertuples(index=False)
        ]

    return Metadata(samples=samples, groups=groups_seen, df=df)
