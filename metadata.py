"""
Parses and validates input metadata CSV.

Input file must have the following required columns:
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
"""

from dataclasses import dataclass, field
from pathlib import Path
import pandas as pd

from tarva.adapters import resolve_adapter

REQUIRED_COLUMNS = ["group", "individualID", "WGS", "RNA", "wgs_barcode", "rna_barcode"]


def _match_files(pattern, directory):
    directory = Path(directory)
    pattern = str(pattern).strip()
    if not pattern:
        raise ValueError("Empty WGS/RNA prefix -- cannot match any files")
    return sorted(directory.glob(f"{pattern}*"))


def _seq_dir(input_root, seq_type, group, adapter_stem):
    return Path(input_root) / seq_type / group / adapter_stem


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
