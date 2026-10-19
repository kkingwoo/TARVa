"""
Running a pipeline's steps in order, from one command -- what `dbpipe all` and
`analysispipe all` are built on.

Each step module owns its own `add_arguments` and `run`. A runner declares the
UNION of those arguments and hands each step a namespace built from that
step's own defaults, overridden by whatever the runner was given. Nothing is
restated, so a step that gains a flag gets it in the combined command for free
and a runner cannot drift from the steps it runs.

Arguments shared by several steps are declared ONCE and reach all of them:
one --db, one --min-depth, one --output, one --overwrite across the whole run.
Where two steps declare the same option the first step's spelling wins.

Input:  nothing -- this module is imported by the `all` runners.
Output: the combined parser, and the exit status of the steps it ran.
"""

import argparse

from tarva.dbpipe import schema


class _Recorder:
    """
    Stands in for an ArgumentParser to capture what a step would declare.

    Only `add_argument` is used by any step's add_arguments, so recording the
    calls is enough to learn a step's argument names and defaults without
    reaching into argparse's internals.
    """

    def __init__(self):
        self.calls = []

    def add_argument(self, *args, **kwargs):
        self.calls.append((args, kwargs))
        return None


def _dest_of(args, kwargs):
    """
    The attribute name argparse would give this argument.

    argparse uses an explicit dest when given one, otherwise the FIRST long
    option with dashes turned into underscores. Positionals are returned as
    their own name; no step here declares any.
    """
    if "dest" in kwargs:
        return kwargs["dest"]
    for option in args:
        if option.startswith("--"):
            return option[2:].replace("-", "_")
    for option in args:
        if not option.startswith("-"):
            return option.replace("-", "_")
    return None


def _default_of(kwargs):
    """The value argparse would leave in place when the flag is absent."""
    if "default" in kwargs:
        return kwargs["default"]
    if kwargs.get("action") == "store_true":
        return False
    if kwargs.get("action") == "store_false":
        return True
    return None


def step_defaults(module):
    """{dest: default} for one step, read from the step itself."""
    recorder = _Recorder()
    module.add_arguments(recorder)
    defaults = {}
    for args, kwargs in recorder.calls:
        dest = _dest_of(args, kwargs)
        if dest is not None:
            defaults[dest] = _default_of(kwargs)
    return defaults


def add_union(parser, steps):
    """
    Declare every argument the steps between them accept, once each.

    Returns {option: step name} for the options that were shared, which the
    runner reports so it is visible that one flag drove several steps.
    """
    seen, shared = set(), {}
    for name, module in steps:
        recorder = _Recorder()
        module.add_arguments(recorder)
        for args, kwargs in recorder.calls:
            options = [option for option in args if option.startswith("-")]
            if any(option in seen for option in options):
                shared.setdefault(options[0], []).append(name)
                continue
            seen.update(options)
            parser.add_argument(*args, **kwargs)
    return shared


def run_steps(args, steps, label):
    """
    Each step in turn, with its own namespace. Stops at the first failure.

    A step raising SystemExit stops the run with its own message, which is
    what should happen: the steps are ordered because each one's inputs are an
    earlier one's outputs, so continuing past a failure would compound it.
    """
    total = len(steps)
    for index, (name, module) in enumerate(steps, start=1):
        defaults = step_defaults(module)
        values = {dest: getattr(args, dest, default)
                  for dest, default in defaults.items()}

        schema.report("")
        schema.report("=" * 72)
        schema.report(f"{label} [{index}/{total}] {name}")
        schema.report("=" * 72)
        module.run(argparse.Namespace(**values))

    schema.report("")
    schema.report(f"{label}: all {total} step(s) complete")
    return 0
