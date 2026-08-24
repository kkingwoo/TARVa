"""
Parses and validates input metadata CSV.

Input file must have the following columns, in order:
    [0] group - contains exactly TWO distinct condition or tissue values
        Example: tumor, peritumoral

    [1] individualID - unique identifier for each individual
        Example: LungF2

    [2] WGS - FASTQ filename prefix in all WGS files that belong to individualID 
        Example:files = [ERR000032_1.fastq.gz, ERR000032_2.fastq.gz] 
                WGS value = ERR000032

    [3] RNA - FASTQ filename prefix in all RNA-seq files that belong to individual
        Example: file = SQXR24385.fq
                 RNA value = SQXR24385
"""

from dataclasses import dataclass, field
from pathlib import Path
import pandas as pd

REQUIRED_COLUMNS = ["group","individualID","WGS","RNA"]

def _match_files(pattern, directory):
    directory = Path(directory)
    prefix = str(prefix.strip())
    if not prefix:
        raise ValueError("Empty WGS/RNA prefix -- cannot match any files")
    return sorted(director.glob(f"{prefix}*"))

@dataclass
class Sample:
    individual_id: str
    group: str
    wgs_pattern: str
    rna_pattern: str
    extra: dict = field(default_factory=dict)

    def wgs_files(self, directory):
        return _match_files(self.wgs_pattern, directory)

    def rna_files(self, directory):
        return _match_files(self.rna_pattern, directory)



@dataclass
class Metadata:
    samples: list
    groups: list

    def samples_in_group(self,group):
        return [s for s in self.samples if s.group == group]

def load_metadata(path):
    """ 
    Reads can validates metadata CSV at 'path'.

    Raise ValueError if:
        1.) required columns are missing
        2.) column does not contain exactly 2 distinct values
        3.) any rows are missing a required field
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Metadata file not found: {path}")
    
    df = pd.read_csv(path, dtype=str)
    df.columns = [c.strip() for c in df.columns]

    missing_cols = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing_cols:
        raise ValueError(
                f"Metadata file {path} is missing required column(s): {missing_cols}"
                f"Required columns are: {REQUIRED_COLUMNS}"
                )
    
    for col in df.columns:
        df[col] = df[col.str.strip()]

    missing_mask = df[REQUIRED_COLUMNS].isna() | (df[REQUIRED_COLUMNS] == "")
    if missing_mask.any().any():
        bad_rows = df.index[missing_mask.any(axis=1)] + 2
        bad_cols = [c for c in REQUIRED_COLUMNS if missing_mask[c].any()]
        raise ValueError(
                f"Metadata file {path}: missing required value(s) in column(s) {bad_cols}"
                f"at line(s) {list(bad_rows)}"
                )
    groups_seen = list(dict.fromkeys(df"group"[])

    if len(groups_seen) !=2:
        raise ValueError(
            f"Metadata file {path}: the 'group' column must contain exactly 2 distinct groups"
            f"distinct values, found {len(groups_seen)}: {groups_seen}"
            )
    extra_cols = [c for c in df.columns if c not in REQUIRED_COLUMNS]

    samples = [
        Sample(
            individual_id=row.individualID,
            group=row.group,
            wgs_pattern=row.WGS,
            rna_pattern=row.RNA,
            extra={c: getattr(row,c) for c in extra_cols},
            )
        for row in df.itertuples(index=False)
        ]

    return Metadata(samples=samples, groups=groups_seen, df=df)


