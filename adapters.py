"""
Shared adapter FASTA resolution, used by `tarva qc` and `tarva sort-input`.
"""

import sys
from pathlib import Path

ADAPTERS_DIR = Path(__file__).parent / "adapters"
ILLUMINA_ADAPTER_DOCS_URL = "https://support-docs.illumina.com/SHARE/AdapterSequences/Content/AdapterSeq.htm"


def list_available_adapters():
    return sorted(ADAPTERS_DIR.glob("*.fa"))


def print_available_adapters():
    adapters = list_available_adapters()
    if not adapters:
        print(f"No adapter FASTA files found in {ADAPTERS_DIR}.", file=sys.stderr)
        return
    print("Available adapter FASTA files:")
    for path in adapters:
        print(f"  {path.name}")
    print(
            "Put one of these filenames in each sample's "
            "wgs_barcode/rna_barcode metadata column, or a path to your own adapter FASTA."
            )
    print(
            "If none of these fit your library prep, Illumina's official adapter "
            f"sequence reference is here: {ILLUMINA_ADAPTER_DOCS_URL}. You can also put a "
            "custom adapter file in the path if needed."
            )
    print(
            "To use a custom adapter -- whether downloaded from the Illumina link above, "
            "sourced from somewhere else, or created yourself -- either:\n"
            f"  1.) place the .fa file directly in {ADAPTERS_DIR} so it shows up in this "
            "list for any metadata.csv, or\n"
            "  2.) leave it wherever it already is and put that file's exact path (not "
            "just the filename) in the wgs_barcode/rna_barcode cell.\n"
            "Either way, the filename must end in '.fa'."
            )


def resolve_adapter(adapter_arg):
    if adapter_arg is None or not str(adapter_arg).strip():
        raise ValueError(
                "No adapter specified. Every sample's wgs_barcode/rna_barcode metadata "
                "cell must contain the filename of an adapter FASTA, ending in '.fa' -- "
                "either one already bundled in tarva/adapters/ (run `tarva qc "
                "--list-adapters` to see what's there) or a path to your own adapter FASTA."
                )
    adapter_arg = str(adapter_arg).strip()
    if not adapter_arg.endswith(".fa"):
        raise ValueError(
                f"Adapter '{adapter_arg}' is invalid -- wgs_barcode/rna_barcode metadata "
                f"values must end in '.fa' (e.g. 'TruSeq3-PE-2.fa'), whether it names a "
                f"bundled adapter or a path to your own adapter FASTA."
                )

    for candidate in list_available_adapters():
        if candidate.name == adapter_arg:
            return candidate

    path = Path(adapter_arg)
    if path.is_file():
        return path

    bundled = list_available_adapters()
    bundled_names = ", ".join(p.name for p in bundled) if bundled else "(none found)"
    raise FileNotFoundError(
            f"Adapter '{adapter_arg}' was not found -- checked the bundled adapters in "
            f"{ADAPTERS_DIR} ({bundled_names}) and '{adapter_arg}' as a file path.\n"
            f"To fix this, either:\n"
            f"  1.) change wgs_barcode/rna_barcode to one of the bundled names above, or\n"
            f"  2.) add your own adapter FASTA. Wherever you got it from -- downloaded "
            f"from Illumina's reference ({ILLUMINA_ADAPTER_DOCS_URL}), pulled from "
            f"another source, or built by hand -- put the file at either:\n"
            f"       - {ADAPTERS_DIR / adapter_arg}  (bundles it alongside the others, "
            f"usable by any metadata.csv), or\n"
            f"       - any path of your choosing, then set the metadata cell to that "
            f"exact path (e.g. '/home/you/adapters/{adapter_arg}').\n"
            f"The file itself must end in '.fa' and be a valid Trimmomatic adapter FASTA."
            )
