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

#: module names each code is known by.  Sites spell Quantum ESPRESSO every way
#: there is -- "qe/7.3" here, "QuantumESPRESSO/7.1" at the next site -- and a
#: name this list misses means the header comes out with a TODO where the
#: module load belongs, on a machine that has the module.  Matching is
#: case-insensitive on the part before the "/".
MODULE_NAMES = {
    "qe": ("quantum-espresso", "quantum_espresso", "quantumespresso",
           "quantum-espresso-gpu", "qe-gpu", "espresso", "qe"),
    "vasp": ("vasp", "vasp-gpu", "vasp_gpu"),
}

#: spellings accepted for the code itself, so `--init-header QuantumEspresso`
#: works as well as `--init-header qe`
CODE_ALIASES = {
    "qe": "qe",
    "quantumespresso": "qe",
    "quantum-espresso": "qe",
    "quantum_espresso": "qe",
    "espresso": "qe",
    "pw": "qe",
    "vasp": "vasp",
}


def normalise_code(code: str) -> str:
    """Map any accepted spelling of a code onto ``"qe"`` or ``"vasp"``.

    Raises
    ------
    ValueError
        When the spelling is not one HTESP writes headers for.
    """
    key = str(code).strip().lower().replace(" ", "")
    try:
        return CODE_ALIASES[key]
    except KeyError:
        raise ValueError(
            "code must be one of {}, not {!r}"
            .format(", ".join(sorted(CODE_ALIASES)), code)) from None

#: the executable each code is driven by, used to check a module really provides it
CODE_EXECUTABLE = {"qe": "pw.x", "vasp": "vasp_std"}

#: MPI launchers, in the order they are preferred when several exist
LAUNCHERS = ("ibrun", "srun", "mpirun", "mpiexec")

#: ``name`` or ``name/version`` and nothing else.  Used to tell a prerequisite
#: line apart from the prose that follows it in ``module spider`` output.
_MODULE_SPEC = re.compile(r"^[A-Za-z0-9_.+-]+(/[A-Za-z0-9_.+-]+)*$")

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
    try:
        names = MODULE_NAMES[normalise_code(code)]
    except ValueError:
        return []
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


def loaded_modules() -> list[str]:
    """What this shell already has loaded, as ``name/version`` strings."""
    lmod = os.environ.get("LMOD_CMD")
    if not lmod or not Path(lmod).exists():
        return []
    out = []
    for line in _run([lmod, "bash", "--terse", "list"]).splitlines():
        line = line.strip()
        if not line or line.endswith(":") or "=" in line or line.startswith("export"):
            continue
        out.append(line)
    return out


def _version_key(name: str) -> tuple:
    """Sort key for ``qe/7.3`` / ``vasp/5.4.4.pl2``: numbers numerically."""
    _, _, version = name.partition("/")
    parts = []
    for chunk in re.split(r"[._-]", version):
        parts.append((0, int(chunk)) if chunk.isdigit() else (1, chunk))
    return tuple(parts)


def latest(names) -> str | None:
    """The highest-versioned of *names*.

    ``modules()`` puts Lmod's ``(D)`` default first, which is the right thing
    to *load*.  This is for probing: the newest build is the one whose help
    text describes the toolchain the site currently expects, and on a site
    that defaults to an older build the newer one is still what a new study
    should be told about.
    """
    names = list(names)
    return max(names, key=_version_key) if names else None


def prerequisites(module: str) -> tuple[list[str], str]:
    """Modules that must be loaded before ``module`` can be.

    On a **hierarchical** Lmod site -- which most HPC centres now are -- an
    application module is not visible until the compiler and MPI it was built
    against are loaded.  ``qe/7.3`` on this machine lives under
    ``/opt/apps/nvidia24/openmpi5/modulefiles``, so ``module load qe`` in a job
    script fails with "these module(s) exist but cannot be loaded as
    requested" unless ``nvidia`` and ``openmpi`` came first.  It works
    interactively only because the login shell already has them.

    ``module help`` is asked **first**, and its answer is taken whenever it
    gives one.  That text is written by whoever built the module: on
    Bridges-2 it says ``> module load intel-oneapi
    QuantumEspresso/7.5-intel``, naming a runtime dependency that Lmod's
    hierarchy does not model at all.  It is a statement of intent, where
    ``spider`` is a derivation from the module tree.

    ``module spider`` is consulted only when help names no prerequisite --
    which is the common case, since most help text says no more than
    ``module load qe/7.3``.  Vista is exactly that: help adds nothing, and
    spider is what reveals ``nvidia cuda openmpi``.  So "help first" is a
    preference, not an exclusion; dropping spider when help merely *exists*
    would break every hierarchical site.

    Returns
    -------
    (names, line)
        ``names`` are unversioned, to be loaded in order; ``line`` is the exact
        versioned combination spider gave, for the comment above them, and is
        empty when the answer came from help.
    """
    lmod = os.environ.get("LMOD_CMD")
    if not lmod or not Path(lmod).exists() or not module:
        return [], ""
    # `module help qe/7.3` and `module help qe` both work -- the second
    # resolves to the site default -- but a site that names its builds
    # `QuantumEspresso/7.5-intel` may only carry the help on one of them, so
    # try the exact build first and fall back to the bare name.
    helptext = _run([lmod, "bash", "help", module])
    stem = module.split("/", 1)[0]
    if stem != module and "module load" not in helptext:
        helptext = _run([lmod, "bash", "help", stem])
    from_help = _help_prerequisites(module, helptext)
    if from_help:
        return from_help, ""

    text = _run([lmod, "bash", "spider", module])
    marker = "You will need to load all module(s) on any one of the lines below"
    if marker not in text:
        # A flat site says "This module can be loaded directly" and has no
        # hierarchy block -- but spider's own Help copy may still spell out a
        # toolchain that `module help` did not.
        return _help_prerequisites(module, text), ""
    block = text.split(marker, 1)[1]
    options = []
    for line in block.splitlines()[1:]:
        line = line.strip()
        if not line:
            if options:                      # the list ends at the blank line
                break
            continue
        if line.startswith("-") or "=" in line or ":" in line:
            break
        parts = line.split()
        # The "Help:" section that follows carries prose and its own
        # "module load qe/7.3" line -- which is the *incomplete* advice this
        # function exists to correct, and which would parse as the three
        # "modules" module, load and qe.  Only accept a line whose every token
        # is a bare module spec.
        if not all(_MODULE_SPEC.match(part) for part in parts):
            break
        options.append(parts)
    if not options:
        return [], ""
    have = set(loaded_modules())
    best = max(options, key=lambda parts: len(have.intersection(parts)))
    names = [part.split("/", 1)[0] for part in best]
    for extra in _help_prerequisites(module, text):
        if extra not in names:
            names.append(extra)
    return names, "  ".join(best)


def _help_prerequisites(module: str, text: str) -> list[str]:
    """Modules the module's own Help text says to load alongside it.

    Lmod's hierarchy does not express every dependency.  Bridges-2 reports
    ``This module can be loaded directly: module load
    QuantumEspresso/7.5-intel`` -- no hierarchy at all -- while the Help
    underneath says::

        To load the module type
        > module load intel-oneapi QuantumEspresso/7.5-intel

    ``intel-oneapi`` carries the Intel MPI and MKL runtimes that build is
    linked against.  Loading QE alone puts ``pw.x`` on ``PATH`` and then fails
    at run time on a missing shared library, which is a far more confusing
    failure than a module that refuses to load.

    Only the tokens *before* the module itself are taken, and only from a line
    that actually names it; the Help of a hierarchical module says plainly
    ``module load qe/7.3``, which yields nothing, as it should.
    """
    stem = module.split("/", 1)[0].lower()
    for line in text.splitlines():
        line = line.strip().lstrip(">$ ").strip()
        if not line.startswith("module load "):
            continue
        parts = line[len("module load "):].split()
        if not all(_MODULE_SPEC.match(part) for part in parts):
            continue
        before = []
        for part in parts:
            if part.split("/", 1)[0].lower() == stem:
                return [name.split("/", 1)[0] for name in before]
            before.append(part)
    return []


def launcher() -> str | None:
    """The MPI launcher this machine appears to use."""
    for name in LAUNCHERS:
        if shutil.which(name):
            return name
    return None


def build(code: str, partition: str | None = None, account: str | None = None,
          time: str = "1-0", module: str | None = None) -> str:
    """Return the text of a starting ``batch.header`` for ``code``."""
    code = normalise_code(code)

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
        # A hierarchical site hides the application module until its compiler
        # and MPI are loaded, so those have to come first or the job dies on
        # "these module(s) exist but cannot be loaded as requested".
        # Probe the newest build: its help text describes the toolchain the
        # site currently expects.  What gets *loaded* is still the bare name,
        # so Lmod resolves it to the site default.
        needs, exact = prerequisites(latest(mods) or module)
        if needs:
            if exact:
                lines.append("# {} is built against a specific compiler/MPI; "
                             "'module spider {}'".format(module, module))
                lines.append("# reports this combination, so load it first:")
                lines.append("#   " + exact)
            else:
                lines.append("# {}'s own 'module spider' help says to load "
                             "these alongside it".format(module))
                lines.append("# (the runtimes it is linked against):")
            lines.append("module load " + " ".join(needs))
        # Load the bare name, not "qe/7.3": Lmod then resolves it to whatever
        # the site has marked default, so the header keeps working when 7.3 is
        # retired -- which it will be, long before anyone edits this file
        # again.  The versions found are listed below so the choice can still
        # be pinned by hand when a study needs one exact build.
        stem = module.split("/", 1)[0]
        lines.append("module load {}".format(stem))
        if mods:
            lines.append("#   versions available now: " + ", ".join(mods))
            lines.append("#   (pin one by writing it out: module load {})"
                         .format(mods[0]))
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


#: printed after every generated header.  The probes answer what this machine
#: reports; they cannot know what a particular build needs at run time, and a
#: module chain that is one module short fails inside the job, minutes after
#: the queue accepted it, with an error that names a shared library rather
#: than a module.
WRITE_WARNING = """
CHECK THIS FILE BEFORE SUBMITTING WITH IT.
  It was assembled from what SLURM and Lmod report here, which is a starting
  point and not a working job script.  In particular:
    * make sure every module the build needs is loaded, including its
      dependencies -- the toolchain written above comes from 'module spider'
      and a site can have requirements Lmod does not model;
    * check it in a login shell first:
          source {path} && which {code}
      If that prints nothing, the module chain is incomplete.
    * the node count and wall time are placeholders, and the partition is
      whichever this machine offered -- neither knows the size of your study.
"""


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
    print(WRITE_WARNING.format(path=target, code=CODE_EXECUTABLE[normalise_code(code)]))
    return 0
