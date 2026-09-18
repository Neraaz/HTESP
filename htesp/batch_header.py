#!/usr/bin/env python
"""Build a starting ``batch.header`` by asking the scheduler what exists.

``mainprogram jobscript`` turns ``batch.header`` into the per-stage submission
scripts, but the header itself has always been the user's to write -- and the
one shipped with the examples says ``--partition=dense``, which exists on the
machine it was written for and nowhere else.  Copying an example therefore
produces a job the scheduler rejects before anything runs.

This module fills in what SLURM and Lmod actually report and leaves everything
else as a marked ``TODO``.  It deliberately does **not** produce a header that
claims to be ready to submit: the launcher, the queue a campaign belongs in and
the executable's own flags are site and study decisions, and a plausible-looking
guess at them is worse than a blank marked "fill this in".

What is asked, and of what:

``sinfo``       partitions, their CPU counts, time limits and GRES
``sacctmgr``    the accounts this user may charge to
``scontrol``    whether GRES is configured at all (``GresTypes``)
``$LMOD_CMD``   the available ``qe``/``vasp`` modules, preferring Lmod's default

Every probe is optional: on a machine with no SLURM the header still comes out,
with those fields as TODO.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

#: module names each code is known by, most specific first
MODULE_NAMES = {
    "qe": ("quantum-espresso", "quantum_espresso", "espresso", "qe"),
    "vasp": ("vasp",),
}

#: the executable each code is driven by, used to check a module really provides it
CODE_EXECUTABLE = {"qe": "pw.x", "vasp": "vasp_std"}

#: MPI launchers, in the order they are preferred when several exist
LAUNCHERS = ("ibrun", "srun", "mpirun", "mpiexec")

TIMEOUT = 25


def _run(command: list[str]) -> str:
    """Run ``command``, returning stdout (plus stderr, which Lmod writes to).

    Never raises: a missing scheduler is a fact to report, not an error.
    """
    try:
        proc = subprocess.run(command, capture_output=True, text=True,
                              timeout=TIMEOUT, check=False)
    except (OSError, subprocess.SubprocessError):
        return ""
    return (proc.stdout or "") + (proc.stderr or "")


def partitions() -> list[dict]:
    """Partitions from ``sinfo``; the default one carries ``default: True``."""
    if shutil.which("sinfo") is None:
        return []
    out = _run(["sinfo", "-h", "-o", "%P|%c|%l|%G"])
    found, seen = [], set()
    for line in out.splitlines():
        parts = line.split("|")
        if len(parts) < 4:
            continue
        name, cpus, timelimit, gres = (p.strip() for p in parts[:4])
        default = name.endswith("*")
        name = name.rstrip("*")
        if not name or name in seen:
            continue
        seen.add(name)
        found.append({"name": name, "cpus": cpus, "time": timelimit,
                      "gres": "" if gres in ("(null)", "") else gres,
                      "default": default})
    return found


def accounts() -> list[str]:
    """Accounts this user may charge, from ``sacctmgr``."""
    if shutil.which("sacctmgr") is None:
        return []
    out = _run(["sacctmgr", "-nP", "show", "assoc",
                "user=" + os.environ.get("USER", ""), "format=Account"])
    names, seen = [], set()
    for line in out.splitlines():
        name = line.split("|")[0].strip()
        if name and name not in seen and not name.startswith("sacctmgr"):
            seen.add(name)
            names.append(name)
    return names


def gres_configured() -> bool:
    """True when the cluster configures GRES at all.

    Vista is the case that makes this necessary: it is a Grace-Hopper GPU
    machine whose ``GresTypes`` is ``(null)`` and whose nodes report
    ``Gres=(null)``.  GPUs come from choosing the ``gh`` partition, and an
    emitted ``--gres=gpu:1`` would be a flag the scheduler rejects.  So GRES is
    written only where the scheduler says GRES exists.
    """
    if shutil.which("scontrol") is None:
        return False
    for line in _run(["scontrol", "show", "config"]).splitlines():
        if line.strip().startswith("GresTypes"):
            value = line.split("=", 1)[-1].strip()
            return bool(value) and value != "(null)"
    return False


def modules(code: str) -> list[str]:
    """Available modules for ``code``; Lmod's default (``(D)``) comes first.

    ``module`` is a shell function, not a program, so this goes through
    ``$LMOD_CMD`` -- and Lmod prints its listing on **stderr** while stdout
    carries shell code to eval, which is why :func:`_run` keeps both.

    The ``(D)`` marker is preferred over the highest version number: it is the
    site's deliberate choice, and a site that defaults to an older build has a
    reason.
    """
    lmod = os.environ.get("LMOD_CMD")
    if not lmod or not Path(lmod).exists():
        return []
    names = MODULE_NAMES.get(code, ())
    listing = _run([lmod, "bash", "avail"])
    found, default, previous = [], None, None
    # Lmod pads the default marker away from the name --
    # "vasp/6.4.3                (D)" -- so a whitespace split makes "(D)" its
    # own token, marking the module *before* it.  Handle both spellings.
    for token in re.split(r"\s{2,}|\n", listing):
        token = token.strip()
        if not token:
            continue
        if token == "(D)":
            if previous:
                default = previous
            continue
        is_default = token.endswith("(D)")
        name = token[:-3].strip() if is_default else token
        stem = name.split("/", 1)[0].lower()
        if stem not in names:
            previous = None
            continue
        if name not in found:
            found.append(name)
        previous = name
        if is_default:
            default = name
    if default:
        found.remove(default)
        found.insert(0, default)
    return found


def launcher() -> str | None:
    """The MPI launcher this machine appears to use."""
    for name in LAUNCHERS:
        if shutil.which(name):
            return name
    return None


def build(code: str, partition: str | None = None, account: str | None = None,
          time: str = "1-0", module: str | None = None) -> str:
    """Return the text of a starting ``batch.header`` for ``code``."""
    code = code.lower()
    if code not in MODULE_NAMES:
        raise ValueError("code must be 'qe' or 'vasp', not {!r}".format(code))

    parts = partitions()
    chosen = None
    if partition:
        chosen = next((p for p in parts if p["name"] == partition),
                      {"name": partition, "cpus": "", "time": "", "gres": ""})
    elif parts:
        chosen = next((p for p in parts if p["default"]), parts[0])

    accts = accounts()
    account = account or (accts[0] if accts else None)
    mods = modules(code)
    module = module or (mods[0] if mods else None)
    run = launcher()
    executable = CODE_EXECUTABLE[code]

    sbatch = ["#SBATCH"]
    if chosen:
        sbatch.append("--partition=" + chosen["name"])
    else:
        sbatch.append("--partition=<TODO>")
    if account:
        sbatch.append("--account=" + account)
    cpus = (chosen or {}).get("cpus", "").rstrip("+")
    sbatch.append("--nodes=1")
    sbatch.append("--ntasks-per-node={}".format(cpus or "<TODO>"))
    sbatch.append("--time=" + time)
    if gres_configured() and (chosen or {}).get("gres"):
        sbatch.append("--gres=" + chosen["gres"].split("(")[0])

    lines = ["#!/bin/bash", "",
             "# Generated by 'mainprogram jobscript --init-header {}'."
             .format(code),
             "# Values below were read from this machine; anything marked TODO",
             "# could not be, and is yours to fill in.", "",
             " ".join(sbatch),
             "# TODO: --nodes and --time are placeholders -- nothing on this",
             "#       machine says how big or how long your study is.", ""]

    if parts and len(parts) > 1:
        others = ", ".join("{} ({} cpus, {})".format(p["name"], p["cpus"], p["time"])
                           for p in parts[:6])
        lines += ["# partitions on this cluster: " + others,
                  "#   the default is not always the right one -- a *-dev queue",
                  "#   is for testing, not for a campaign.", ""]
    if len(accts) > 1:
        lines += ["# accounts you may charge: " + ", ".join(accts), ""]

    if module:
        lines.append("module load {}".format(module))
        if len(mods) > 1:
            lines.append("#   also available: " + ", ".join(mods[1:]))
    else:
        lines.append("# TODO: module load <your {} module>  "
                     "(none detected on this machine)".format(code.upper()))
    lines.append("")

    # No run command goes in the file: 'mainprogram jobscript' copies this
    # header verbatim and appends its own line, built from config.json.  A
    # command written here would be run as well, before that one.
    lines += ["# 'mainprogram jobscript' appends the {} command here, built"
              .format(executable),
              "# from config.json:"]
    if run:
        lines += ['#     "parallel_command": "{}",   (found on $PATH)'.format(run),
                  '#     "nproc": "{}"'.format(cpus or "<TODO>")]
        if run == "ibrun":
            lines += ["#   ibrun is handed no process count -- it runs the whole",
                      "#   allocation, so nproc only sizes the #SBATCH line above."]
    else:
        lines += ["# TODO: no MPI launcher (srun / ibrun / mpirun) found on",
                  "#       $PATH -- set parallel_command to yours."]
    lines.append("")
    return "\n".join(lines)


def write(path, code: str, force: bool = False, **kwargs) -> int:
    """Write :func:`build` to ``path``.  Returns a process exit status."""
    target = Path(path)
    if target.exists() and not force:
        print("{} already exists; pass --force to overwrite it".format(target))
        print("  (it is usually hand-tuned -- look before replacing it)")
        return 1
    text = build(code, **kwargs)
    target.write_text(text)
    print("wrote {}".format(target))
    for line in text.splitlines():
        if line.startswith("#SBATCH") or line.startswith("module load"):
            print("  " + line)
    todos = sum(1 for line in text.splitlines() if "TODO" in line)
    if todos:
        print("  {} TODO line(s) to fill in before submitting".format(todos))
    return 0
