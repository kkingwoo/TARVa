"""
Shared plumbing for the R-backed steps of analysispipe.

Each step resolves its parameters in Python, merges them into
<output>/params.json, then runs `Rscript <NN_step>.R <output>/params.json`.
The R scripts take no other arguments, so every path and setting is decided in
exactly one place, and params.json doubles as the provenance record for the
run. Parameters accumulate rather than overwrite, so a later step can read a
setting an earlier step resolved -- and any step can be re-run by hand with
`Rscript R/<NN_step>.R <output>/params.json`.

Note the generated SLURM script loads only the R module -- deliberately no
conda activation, since conda's LD_LIBRARY_PATH can shadow the system
libraries a site-built R was compiled against.

Input:  nothing -- this module is imported by the R-backed steps.
Output: <output>/params.json, and the Rscript invocation or submitted job.
"""

import json
import subprocess

from datetime import datetime
from pathlib import Path

from tarva.slurm_utils import parse_directives_dict, split_directives, submit_job, write_sbatch_script

R_DIR = Path(__file__).parent / "R"


def add_common_arguments(parser):
    """
    Arguments every R-backed step takes. Steps that read or write the TARVa
    database add their own --db; it is deliberately not here, so a step that
    does not touch sample data does not advertise a flag it ignores.
    """
    parser.add_argument("--output", required=True, help="Directory to write results into")
    parser.add_argument(
            "--r-module", default=None, metavar="NAME",
            help=(
                "R module to `module load` in the generated SLURM script, e.g. 'R/4.5.1'. "
                "Pin the exact version -- personal R libraries are version-scoped, so a "
                "different module sees none of the installed packages. Only used with --hpc."
                ),
            )
    parser.add_argument(
            "--r-libs", default=None, metavar="PATH",
            help=(
                "Value to export as R_LIBS_USER before running Rscript, if the personal "
                "library is not already on the path non-interactively. Only used with --hpc."
                ),
            )
    parser.add_argument(
            "--hpc", action="store_true",
            help="Submit this step as a SLURM job instead of running Rscript locally",
            )
    parser.add_argument(
            "--directives-dict", dest="directives_dict", type=parse_directives_dict, default={},
            metavar="DICT",
            help=(
                "SLURM sbatch directives as a Python dict literal, only used with --hpc. "
                "See `tarva qc --help` for recognized keys. Example: "
                "\"{'partition': 'Orion', 'mem': '32gb', 'time': '02:00:00'}\""
                ),
            )


def r_script_path(script_name):
    path = R_DIR / script_name
    if not path.is_file():
        raise SystemExit(
                f"tarva analysispipe: bundled R script not found: {path}. Is the package "
                f"directory complete?"
                )
    return path


def write_params(args, step, extra):
    """
    Merge `extra` into <output>/params.json and return its path. Existing keys
    from earlier steps are preserved so params.json accumulates the full run.
    """
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    params_path = output_dir / "params.json"

    params = json.loads(params_path.read_text()) if params_path.is_file() else {}
    params["output_dir"] = str(output_dir.resolve())
    params[f"{step}_run_at"] = datetime.now().isoformat(timespec="seconds")
    params["r_module"] = args.r_module

    db = getattr(args, "db", None)
    if db:
        params["db"] = str(Path(db).resolve())

    params.update(extra)

    params_path.write_text(json.dumps(params, indent=2, sort_keys=True) + "\n")
    print(f"Wrote parameters: {params_path}")
    return params_path


def run_r_step(args, step, script_name, extra, job_name=None):
    params_path = write_params(args, step, extra)
    script = r_script_path(script_name)
    if args.hpc:
        _submit_hpc(args, script, params_path, job_name or f"tarva_{step}")
    else:
        cmd = ["Rscript", str(script), str(params_path)]
        print("Running: " + " ".join(cmd))
        subprocess.run(cmd, check=True)


def _submit_hpc(args, script, params_path, job_name):
    output_dir = Path(args.output)
    log_dir = output_dir / "logs"

    commands = []
    if args.r_libs:
        commands.append(f"export R_LIBS_USER={args.r_libs}")
    commands.append(f"Rscript {script} {params_path}")

    sbatch_kwargs, extra_directives = split_directives(args.directives_dict)

    script_path = write_sbatch_script(
            job_name=job_name,
            commands=commands,
            log_dir=log_dir,
            script_path=output_dir / f"{job_name}.slurm",
            modules=[args.r_module] if args.r_module else None,
            extra_directives=extra_directives,
            **sbatch_kwargs,
            )
    print(f"Wrote SLURM script: {script_path}")
    print(submit_job(script_path))
