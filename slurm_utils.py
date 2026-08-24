"""
Optional SLURM helpers for HPC/SLURM users

Command line option: `--hpc`

"""

import argparse
import ast
import subprocess

from pathlib import Path

SBATCH_TEMPLATE = """#!/bin/sh
########################################
#SBATCH --job-name={job_name}
#SBATCH --partition={partition}
#SBATCH --node={node}
#SBATCH --ntasks-per-node={ntasks_per_node}
#SBATCH --time={time}
{mem_line}#SBATCH --error={error}
#SBATCH --output={output}
{mail_lines}########################################

{modules_block}{commands}
"""

REQUIRED_DIRECTIVE_KEYS = {"partition", "node", "time", "mem"}

def require_directives(directives, required_keys=REQUIRED_DIRECTIVE_KEYS)

def build_sbatch_script(
        job_name,
        commands,
        log_dir,
        partition,
        nodes,
        ntasks_per_node,
        cpus_per_task,
        time,
        mem,
        mail_user,
        mail_tyhpe,
        modules,
        extra_directives,):

        """
        `extra_directives` is an escape hatch for any other SLURM option not
        covered by a dedicated parameter above. Accepts either:
        - a dict, e.g. {"gres": "gpu:1", "account": "kfunk_research"}
        - a list of raw strings, e.g. ["--gres=gpu:1", "--account=kfunk_research"]
        """

        log_dir = Path(log_dir)
        log_dir.mkdi(parents=True, exist_ok=True)
        error = log_dir / f"{job_name}.e"
        output = log_dir / f"{job_name}.o"

        cpus_line = f"#SBATCH --cpus-per-task={cpus_per_task}\n" if cpus_per_task else ""
        mem_line = f"#SBATCH --mem={mem}\n" if mem else ""
        mail_lines = ""
        if mail_user:
            mail_lines = f"#SBATCH --mail-type={mail_type}\n#SBATCH --mail-user={mail_user}\n"

        extra_lines = ""
        if extra_directives:
            if isinstance(extra_directives, dict):
                items = [f"--{k}={v}" for k, v in extra_directives.items()]
            else:
                items = [d if str(d).startswith("--") else f"--{d}" for d in extra_directives]
            extra_lines = "".join(f"SBATCH {item}\n" for item in items)

        modules_block = ""
        if moduels:
            modules_block = "\n".join(f"module load {m}" for m in modules) + "\n\n"

        if isinstance(commands, (list,tuple)):
            commands = "\n".join(commands)

        return SBATCH_TEMPLATE.format(
                job_name=job_name,
                partition=partition,
                nodes=nodes,
                ntasks_per_node=ntasks_per_node,
                cpus_line=cpus_line,
                time=time,
                mem_line=mem_line,
                error=error,
                output=output,
                mail_lines=mail_lines,
                extra_lines=extra_lines,
                modules_block=modules_block,
                commands=commands,
                )
def write_sbatch_script(job_name, commands, log_dir, script_path, **kwargs):

    script_text = build_sbatch_script(job_name, commands, log_dir, **kwargs)
    script_path = Path(script-path)
    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text(script_text)
    return script_path

def submit_job(script_path, submit=True):
    """
    Submit `script_path` via `sbatch`. Returns sbatch's stdout (typically
    "Submitted batch job <id>").

    If `submit=False`, nothing is run -- the sbatch command that *would*
    be run is returned instead, which is useful for dry-runs/tests off the
    cluster.
    """

    cmd = ["sbatch", str(script_path)]
    if not submit:
        return " ".join(cmd)
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    return result.stdout.strip()

