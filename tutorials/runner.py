#!/usr/bin/env python
"""Execution engine for the HTESP tutorial runner.

Responsibilities, in the order they happen:

1. :func:`preflight` -- check *everything* that can be checked before the first
   job is submitted, and report all of it at once.  A campaign that dies two
   hours in because ``MP_API_KEY`` was never exported has wasted two hours.
2. :func:`seed_workdir` -- build an isolated work directory per tutorial.
   ``examples/`` is read-only input: nothing is ever written there.
3. :class:`TutorialRunner` -- run each step as a ``mainprogram`` *subprocess*,
   wait for the cluster jobs it submitted, verify the artefacts it claimed to
   produce, and checkpoint after every single step.

``htesp`` is never imported here.  ``mainprogram`` is a subprocess so that one
tutorial blowing up cannot take the driver with it, and so that the captured
log file is a real file the stop report can point at.
"""
from __future__ import annotations

import functools
import logging
import threading
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from tutorials.catalog import (CATALOG, EXAMPLES, PACKAGE_ROOT, Step,
                               Tutorial, iters, readme_for,
                               searched_for_examples)
from tutorials.state import (BLOCKED, DONE, FAILED, PENDING, RUNNING, SKIPPED,
                             RunState, StepState, TutorialState, step_key)
from tutorials.workdirs import (changed_since, missing_artifacts,
                                patch_input_in, relax_output_present,
                                seed_workdir, snapshot)

LOG = logging.getLogger("htesp.tutorials")

#: There are no modes left.  The runner prepares inputs and never runs a DFT
#: calculation: every step is invoked as ``mainprogram ... --dry-run``, so
#: nothing is submitted and nothing is deleted.  ``--no-dft`` went first (one
#: mode more than the runner needed explaining), then real mode, which existed
#: to drive a campaign -- that is ``mainprogram``'s job, not this driver's.
#: The name survives only as the label written into ``state.json``.
DRY_RUN = "dry-run"

#: default seconds a single ``mainprogram`` call may take
#: how long a whole tutorial may take before the runner gives up on it.
#:
#: Measured over a healthy sweep, the slowest tutorial (QE/4, two live OQMD
#: calls) totalled 58.6 seconds and no tutorial exceeded a minute, so this is
#: tight but real work fits inside it.
DEFAULT_TUTORIAL_TIMEOUT = 60

#: how long one `mainprogram` call may take before the runner kills it.
#:
#: This was six hours, which is how a hung step went unnoticed: OQMD stalled
#: an `oqmd-download` for nineteen minutes and counting, alive with two open
#: sockets and no output, because `qmpy_rester` builds a bare
#: `requests.Session()` with no timeout and a stalled connection blocks for
#: ever.  With no DFT run here, no step does real work for minutes, so a
#: short limit turns a hang into a reported failure instead of a wait.
DEFAULT_STEP_TIMEOUT = 60


def _as_text(blob) -> str:
    """Decode whatever a subprocess handed back; None becomes ``""``."""
    if blob is None:
        return ""
    if isinstance(blob, bytes):
        return blob.decode("utf-8", "replace")
    return str(blob)


def mp_api_key() -> str | None:
    """The Materials Project key, resolved the way HTESP itself resolves it.

    FIX: this used to be ``os.environ.get("MP_API_KEY")`` and nothing else,
    which contradicted the rest of the package.  ``htesp-check --set_mp_api``
    exists precisely so nobody has to export the variable -- an exported key
    is lost by a batch job, a nohup-ed sweep or a new terminal -- and it
    writes ``~/.config/htesp/credentials``.  A correctly configured machine
    therefore had eight tutorials skip with "MP_API_KEY is not set" while
    ``mainprogram search`` run by hand in the same directory worked.

    ``htesp.config.api_key()`` is the one resolver: environment, then the
    credentials file, then ``config.json`` -- returning None for the shipped
    ``use_your_API_KEY`` placeholder.  It is imported lazily because
    ``tutorials/`` is meant to stay importable on a machine with no
    scientific stack, and falls back to the environment if it cannot be.
    """
    try:
        from htesp.config import api_key
    except Exception:                          # noqa: BLE001 - no htesp here
        return os.environ.get("MP_API_KEY", "").strip() or None
    try:
        return api_key()
    except Exception:                          # noqa: BLE001 - unreadable config
        return os.environ.get("MP_API_KEY", "").strip() or None


def child_env(root: Path = PACKAGE_ROOT) -> dict[str, str]:
    """Environment for a ``mainprogram`` subprocess.

    The repository root is prepended to ``PYTHONPATH`` so that the driver works
    in a checkout that has not been ``pip install``-ed -- and so that, when it
    has been, the tutorials are run against *this* tree rather than whatever
    else is on the path.
    """
    env = dict(os.environ)
    existing = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = (f"{root}{os.pathsep}{existing}" if existing else str(root))
    return env


@dataclass
class RunOptions:
    """Everything the runner needs that is not the catalogue."""

    workdir: Path
    resume: bool = True
    poll_interval: float = 60.0
    job_timeout: float = 24 * 3600
    step_timeout: float = DEFAULT_STEP_TIMEOUT
    #: seconds a whole tutorial may take.  A step is never given more than
    #: what is left of it, so an unbounded network call is cut short by the
    #: tutorial's own budget rather than running until the step limit.
    tutorial_timeout: float = DEFAULT_TUTORIAL_TIMEOUT
    workers: int | None = None
    from_step: str | None = None
    #: how many tutorials to run at once.  1 is sequential, as it has always
    #: been.  See TutorialRunner._run_group for why this is safe and where the
    #: time actually goes.
    jobs: int = 1
    #: which dependency iteration to run: None runs the whole selection in one pass
    #: (what every earlier version did), an int runs exactly that iteration, and
    #: "next" runs the lowest-numbered iteration the checkpoint has not finished.
    iteration: int | str | None = None
    verbose: bool = False
    python: str = sys.executable
    examples: Path = EXAMPLES
    #: what to do with ``tutorial_runs/`` when the run ends.  "all" keeps every
    #: work directory (the default, and what every earlier version did) and
    #: "none" keeps none; the command line spells these ``--keep_output yes``
    #: and ``--keep_output no``.  "failed" -- keep only the directories worth
    #: looking at -- is honoured here for callers driving the runner directly,
    #: but is deliberately not a command-line choice: one flag, two answers.
    #: Logs, the report and the checkpoint are never removed; they are the
    #: record of the run, and they are small.
    keep: str = "all"

    @property
    def runs_root(self) -> Path:
        return self.workdir / "tutorial_runs"

    @property
    def logs_root(self) -> Path:
        return self.workdir / "logs"


@dataclass
class Problem:
    """One preflight finding.  ``fatal`` problems stop the run before it starts."""

    level: str          # "error" | "warning"
    message: str
    fix: str = ""

    @property
    def fatal(self) -> bool:
        return self.level == "error"

    def __str__(self) -> str:
        return f"[{self.level}] {self.message}" + (f"  ({self.fix})" if self.fix else "")


# --------------------------------------------------------------------------- #
#  preflight
# --------------------------------------------------------------------------- #
#: the enumlib executables pymatgen's EnumlibAdaptor looks for.  It needs one
#: of the enumerators *and* one of the structure makers, but any of these being
#: present means enumlib was built, which is what preflight is really asking.
ENUMLIB_TOOLS = ("enum.x", "multienum.x", "makestr.x", "makeStr.py")



@functools.lru_cache(maxsize=1)
def potcars_available() -> bool:
    """Can pymatgen produce a POTCAR on this machine?

    True when ``PMG_VASP_PSP_DIR`` points at a directory holding a potential
    family.  Cached: the answer cannot change during a run, and it is asked
    once per VASP step.
    """
    try:
        from pymatgen.core import SETTINGS
    except ImportError:
        return False
    root = SETTINGS.get("PMG_VASP_PSP_DIR")
    if not root:
        return False
    base = Path(root)
    return any((base / name).is_dir() for name in (
        "POT_GGA_PAW_PBE", "POT_GGA_PAW_PBE_52", "POT_GGA_PAW_PBE_54",
        "POT_GGA_PAW_PW91", "POT_LDA_PAW", "POT_LDA_PAW_52", "POT_LDA_PAW_54"))


def preflight(codes: Sequence[str], options: RunOptions,
              catalog: dict[str, Tutorial] | None = None) -> list[Problem]:
    """Check every prerequisite of the selected tutorials, all at once.

    Returns the findings rather than raising, so the caller can print the whole
    list.  Nothing here touches the network or the scheduler.
    """
    catalog = CATALOG if catalog is None else catalog
    out: list[Problem] = []
    tutorials = [catalog[c] for c in codes if c in catalog]
    unknown = [c for c in codes if c not in catalog]
    if unknown:
        out.append(Problem("error", f"unknown tutorial code(s): {', '.join(unknown)}",
                           "run with --list to see the catalogue"))

    if not options.examples.is_dir():
        looked = "; ".join(searched_for_examples())
        out.append(Problem(
            "error", f"the example tree {options.examples} is missing",
            "examples/ is 185 MB and is not shipped inside the wheel, so it is "
            "not beside an installed package. Run from the repository, or pass "
            "--examples /path/to/HTESP/examples, or set $HTESP_EXAMPLES. "
            f"Looked in: {looked}"))
        return out

    # examples/ is read-only input; writing runs into it would mix generated
    # output into the reference tree and is almost always `--workdir examples/`
    # typed when `--examples` was meant
    try:
        inside = options.workdir.resolve().is_relative_to(options.examples.resolve())
    except (OSError, ValueError):                    # pragma: no cover
        inside = False
    if inside:
        out.append(Problem(
            "error", f"--workdir {options.workdir} is inside the example tree",
            "examples/ is read-only input; --workdir is where runs are written. "
            "Point --workdir somewhere else, and use --examples to say where "
            "the example tree is."))

    for tut in tutorials:
        if not tut.directory.is_dir():
            out.append(Problem("error", f"{tut.code}: {tut.directory} does not exist"))
        if tut.stub:
            out.append(Problem("warning", f"{tut.code}: {tut.note or 'stub tutorial'}",
                               "it will be attempted but is expected to fail"))

    codes_used = {t.dft for t in tutorials}
    for dft in sorted(codes_used):
        header = options.examples / dft / "batch.header"
        if not header.is_file():
            out.append(Problem("error", f"{header} is missing",
                               "every submission script is built from it"))
    if "QE" in codes_used:
        pp = options.examples / "QE" / "pp"
        if not pp.is_dir() or not any(pp.glob("*.upf")):
            out.append(Problem("error", f"no pseudopotentials in {pp}",
                               "QE tutorials cannot build an input without them"))

    needs_key = sorted({t.code for t in tutorials
                        if any(s.needs_api_key for s in t.steps)})
    if needs_key and not mp_api_key():
        # A warning, not an error: the run is still useful without a key --
        # the four tutorials that need one are skipped and the rest proceed.
        out.append(Problem(
            "warning",
            "no Materials Project API key is configured, and these "
            "tutorials query the Materials "
            f"Project: {', '.join(needs_key)}",
            "htesp-check --set_mp_api <your key>  (it is verified and stored "
            "in ~/.config/htesp/credentials, which survives batch jobs; "
            "$MP_API_KEY still wins when set)"))

    # enumlib is a separate C/Fortran package that pymatgen's EnumlibAdaptor
    # shells out to; nothing pip-installs it.  This is a warning in every mode,
    # never an error: only the magnetic-ordering tutorials need it, and failing
    # preflight would block the other forty for the sake of two.
    needs_enumlib = sorted({t.code for t in tutorials
                            if any(s.needs_enumlib for s in t.steps)})
    if needs_enumlib and not any(shutil.which(tool) for tool in ENUMLIB_TOOLS):
        out.append(Problem(
            "warning",
            "enumlib is not on PATH, and these tutorials enumerate magnetic "
            f"orderings with it: {', '.join(needs_enumlib)}",
            "they will fail with \"EnumlibAdaptor requires the executables "
            "'enum.x' or 'multienum.x' and 'makestr.x' or 'makeStr.py'\". "
            "Build them from https://github.com/msg-byu/enumlib and put them "
            "on PATH, or leave these tutorials out with "
            f"--skip {','.join(needs_enumlib)}"))

    probe = _probe_mainprogram(options)
    if probe:
        out.append(probe)
    return out


def _probe_mainprogram(options: RunOptions) -> Problem | None:
    """Confirm ``python -m htesp`` can at least start."""
    try:
        proc = subprocess.run([options.python, "-m", "htesp", "--version"],
                              capture_output=True, text=True, timeout=120,
                              env=child_env(),
                              cwd=os.fspath(options.workdir if options.workdir.is_dir()
                                            else Path.cwd()))
    except (OSError, subprocess.SubprocessError) as exc:
        return Problem("error", f"could not run 'python -m htesp': {exc}")
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout).strip().splitlines()[-1:] or [""]
        return Problem("error", f"'python -m htesp --version' exited "
                                f"{proc.returncode}: {detail[0]}",
                       "install the package: pip install -e .")
    return None


def econv_converged(workdir: Path) -> bool | None:
    """Has the relaxation loop converged?  ``None`` when it cannot be told.

    ``mainprogram e0`` writes ``econv.csv`` with an ``niteration`` column; the
    relaxation tutorial says to repeat ``2`` -> ``3`` -> ``e0`` until that
    number drops below 3.
    """
    import csv

    path = workdir / "econv.csv"
    try:
        rows = list(csv.DictReader(path.read_text().splitlines()))
    except (OSError, ValueError):
        return None
    values = []
    for row in rows:
        try:
            values.append(int(str(row.get("niteration", "")).strip()))
        except (TypeError, ValueError):
            return None
    if not values:
        return None
    return max(values) < 3


#: convergence probes a :class:`~tutorials.catalog.Loop` can name
PROBES = {"econv": econv_converged}


# --------------------------------------------------------------------------- #
#  the runner
# --------------------------------------------------------------------------- #
class TutorialRunner:
    """Run a selection of tutorials, checkpointing after every step."""

    def __init__(self, codes: Sequence[str], options: RunOptions,
                 catalog: dict[str, Tutorial] | None = None,
                 state: RunState | None = None):
        self.codes = list(codes)
        self.options = options
        self.catalog = CATALOG if catalog is None else catalog
        self.state = state if state is not None else RunState(
            path=options.workdir / "state.json")
        self.state.mode = DRY_RUN
        # One checkpoint file, one writer at a time.  save() is already atomic
        # (temp file + replace), so this only stops two threads serialising the
        # same dict at once and racing to replace.
        self._lock = threading.Lock()
        # Per-tutorial deadline.  Thread-local because several tutorials can
        # be in flight at once under --jobs.
        self._deadline = threading.local()

    # -- entry point -------------------------------------------------------- #
    def run(self) -> RunState:
        """Run every selected tutorial; always returns a saved checkpoint."""
        self.options.runs_root.mkdir(parents=True, exist_ok=True)
        self.options.logs_root.mkdir(parents=True, exist_ok=True)
        codes = list(self.codes)
        try:
            self._run_group(codes)
        except KeyboardInterrupt:
            self.state.interrupted = True
            self._mark_running_as_interrupted()
            LOG.error("interrupted -- the checkpoint is at %s", self.state.path)
        self._save()
        self._clean_work_dirs()
        return self.state

    def _how_to_run_it(self, tutorial: Tutorial) -> str:
        """``"; to run it for real, follow <path>"`` -- or nothing.

        Only for steps that genuinely need a DFT run.  A step skipped because
        the machine lacks a POTCAR, an API key or enumlib already says what to
        install, and replacing that with "read the README" would be a
        downgrade; a step skipped because the one before it was skipped names
        that step, which is the actual cause.

        VASP/21 has no instructions anywhere -- neither its own directory nor
        a QE counterpart, since the QE tree has no IFermi tutorial -- so it
        gets no pointer rather than a path that does not exist.
        """
        readme = readme_for(tutorial, self.catalog)
        if readme is None:
            return ""
        try:
            where = readme.relative_to(PACKAGE_ROOT)
        except ValueError:                           # pragma: no cover
            where = readme
        return f"; to run it for real, follow {where}"

    def _time_left(self) -> float | None:
        """Seconds left in this tutorial's budget, or None when unbounded."""
        at = getattr(self._deadline, "at", None)
        return None if at is None else at - time.time()

    def _tag(self, code: str) -> str:
        """``"QE/9  "`` when tutorials run in parallel, ``""`` when they do not.

        The step lines carry no tutorial code, because sequentially the
        heading above them says which tutorial it is.  Interleaved, that
        heading is meaningless and the output becomes unreadable.
        """
        if max(1, int(getattr(self.options, "jobs", 1) or 1)) == 1:
            return ""
        return f"{code:<9s} "

    def _save(self) -> None:
        """Write the checkpoint under the lock (threads share one file)."""
        with self._lock:
            self.state.save()

    def _run_group(self, codes: Sequence[str]) -> None:
        """Run *codes*, in parallel when asked, respecting dependencies.

        Nearly all of a sweep's wall time is spent waiting on other people's
        servers: the database front ends query Materials Project, OQMD and
        AFLOW over the network, and in one measured run seven `search` steps
        accounted for over six of the eight minutes elapsed -- the slowest a
        single 163-second call -- while every local step finished in under two
        seconds.  Those tutorials are independent of each other, so the waiting
        can overlap.

        Threads, not processes: every step is a `subprocess.run` of
        `mainprogram`, and the interpreter releases the GIL for the duration,
        so threads give full overlap while keeping one shared checkpoint.

        Dependencies are honoured by running one dependency level at a time --
        a tutorial that seeds from another's *work directory* cannot start
        before it finishes.  The graph is two levels deep, and the slow
        database tutorials are all in the first, so a level barrier costs
        almost nothing next to a full dependency-aware scheduler.
        """
        jobs = max(1, int(getattr(self.options, "jobs", 1) or 1))
        if jobs == 1 or len(codes) < 2:
            for code in codes:
                self._run_tutorial(code)
            return

        from concurrent.futures import ThreadPoolExecutor

        for level in iters(codes, self.catalog):
            if len(level) == 1:
                self._run_tutorial(level[0])
                continue
            width = min(jobs, len(level))
            LOG.info("running %d tutorial(s) %d at a time", len(level), width)
            with ThreadPoolExecutor(max_workers=width) as pool:
                list(pool.map(self._run_tutorial, level))

    def _how_to_run_it(self, tutorial: Tutorial) -> str:
        """``"; to run it for real, follow <path>"`` -- or nothing.

        Only for steps that genuinely need a DFT run.  A step skipped because
        the machine lacks a POTCAR, an API key or enumlib already says what to
        install, and replacing that with "read the README" would be a
        downgrade; a step skipped because the one before it was skipped names
        that step, which is the actual cause.

        VASP/21 has no instructions anywhere -- neither its own directory nor
        a QE counterpart, since the QE tree has no IFermi tutorial -- so it
        gets no pointer rather than a path that does not exist.
        """
        readme = readme_for(tutorial, self.catalog)
        if readme is None:
            return ""
        try:
            where = readme.relative_to(PACKAGE_ROOT)
        except ValueError:                           # pragma: no cover
            where = readme
        return f"; to run it for real, follow {where}"

    def _time_left(self) -> float | None:
        """Seconds left in this tutorial's budget, or None when unbounded."""
        at = getattr(self._deadline, "at", None)
        return None if at is None else at - time.time()

    def _tag(self, code: str) -> str:
        """``"QE/9  "`` when tutorials run in parallel, ``""`` when they do not.

        The step lines carry no tutorial code, because sequentially the
        heading above them says which tutorial it is.  Interleaved, that
        heading is meaningless and the output becomes unreadable.
        """
        if max(1, int(getattr(self.options, "jobs", 1) or 1)) == 1:
            return ""
        return f"{code:<9s} "

    def _save(self) -> None:
        """Write the checkpoint under the lock (threads share one file)."""
        with self._lock:
            self.state.save()

    def _run_group(self, codes: Sequence[str]) -> None:
        """Run *codes*, in parallel when asked, respecting dependencies.

        Nearly all of a sweep's wall time is spent waiting on other people's
        servers: the database front ends query Materials Project, OQMD and
        AFLOW over the network, and in one measured run seven `search` steps
        accounted for over six of the eight minutes elapsed -- the slowest a
        single 163-second call -- while every local step finished in under two
        seconds.  Those tutorials are independent of each other, so the waiting
        can overlap.

        Threads, not processes: every step is a `subprocess.run` of
        `mainprogram`, and the interpreter releases the GIL for the duration,
        so threads give full overlap while keeping one shared checkpoint.

        Dependencies are honoured by running one dependency level at a time --
        a tutorial that seeds from another's *work directory* cannot start
        before it finishes.  The graph is two levels deep, and the slow
        database tutorials are all in the first, so a level barrier costs
        almost nothing next to a full dependency-aware scheduler.
        """
        jobs = max(1, int(getattr(self.options, "jobs", 1) or 1))
        if jobs == 1 or len(codes) < 2:
            for code in codes:
                self._run_tutorial(code)
            return

        from concurrent.futures import ThreadPoolExecutor

        for level in iters(codes, self.catalog):
            if len(level) == 1:
                self._run_tutorial(level[0])
                continue
            width = min(jobs, len(level))
            LOG.info("running %d tutorial(s) %d at a time", len(level), width)
            with ThreadPoolExecutor(max_workers=width) as pool:
                list(pool.map(self._run_tutorial, level))

    # -- iters -------------------------------------------------------------- #
    def _iter_to_run(self, plan: list[list[str]]) -> tuple[int | None, list[str] | None]:
        """Which iteration this invocation covers, and the codes in it.

        ``(None, every code)`` is the historical single-pass run; that is still
        the default: nothing is ever queued, so there is nothing to wait
        between iterations *for*.
        """
        wanted = self.options.iteration
        if wanted is None:
            return None, list(self.codes)
        if wanted == "next":
            number = self.state.next_iter(plan)
            if number is None:
                LOG.info("every iteration is done (%d of %d); nothing to run",
                         len(plan), len(plan))
                return None, None
        else:
            number = int(wanted)
            if not 0 <= number < len(plan):
                LOG.error("iteration %d does not exist: this selection has iters 0-%d",
                          number, len(plan) - 1)
                return None, None
        return number, list(plan[number])

    def _close_iter(self, record, plan: list[list[str]]) -> None:
        """Record how the iteration ended, and say what to run next.

        A iteration is DONE only when every tutorial in it is; anything else leaves
        it FAILED, so ``--iter next`` offers it again rather than stepping over
        it onto structures that were never produced.
        """
        record.finished = time.time()
        bad = [code for code in record.codes
               if self.state.tutorials.get(code, TutorialState(code=code)).status
               not in (DONE, SKIPPED)]
        if self.state.interrupted:
            record.status, record.reason = FAILED, "interrupted with Ctrl-C"
        elif bad:
            record.status = FAILED
            record.reason = "did not finish: " + ", ".join(bad)
        else:
            record.status = DONE
        following = record.number + 1
        if record.status != DONE:
            LOG.error("iteration %d did not finish (%s); fix those, then re-run "
                      "--iteration %d", record.number, record.reason, record.number)
        elif following < len(plan):
            LOG.info("iteration %d done.  Wait for its jobs to finish, then run "
                     "'htesp-tutorials --resume --iter next' for iteration %d (%s)",
                     record.number, following, ", ".join(plan[following]))
        else:
            LOG.info("iteration %d done -- that was the last one.", record.number)

    def _clean_work_dirs(self) -> None:
        """Drop work directories the user asked not to keep.

        A sweep leaves 42 of these, and after a green run they are just bulk --
        but they are also the only place a failure can be examined, so
        ``--keep failed`` (and the ``all`` default) err towards keeping them.
        Cleanup happens once, at the very end: a tutorial's work directory is
        the seed for everything that depends on it, so removing one mid-run
        would starve its dependents.  ``--resume`` after a cleanup re-runs the
        tutorials whose directories went away, which is why this never touches
        ``state.json``, the logs or the report.
        """
        keep = getattr(self.options, "keep", "all")
        if keep == "all" or self.state.interrupted:
            return
        kept_states = {FAILED, BLOCKED} if keep == "failed" else set()
        removed = 0
        for code, record in self.state.tutorials.items():
            if record.status in kept_states:
                continue
            target = self.options.runs_root / code.replace("/", "-")
            if not target.is_dir():
                continue
            try:
                shutil.rmtree(target)
                removed += 1
            except OSError as exc:                 # pragma: no cover
                LOG.warning("could not remove %s: %s", target, exc)
        if removed:
            spelling = {"none": "--keep_output no", "failed": "keep=failed"}
            LOG.info("removed %d work director%s (%s); logs, report and "
                     "checkpoint are kept", removed,
                     "y" if removed == 1 else "ies",
                     spelling.get(keep, keep))

    def _mark_running_as_interrupted(self) -> None:
        for tut in self.state.tutorials.values():
            for step in tut.steps.values():
                if step.status == RUNNING:
                    step.status = FAILED
                    step.reason = "interrupted with Ctrl-C while this step was running"
                    step.finished = time.time()
                    step.duration = step.finished - (step.started or step.finished)
            if tut.status == RUNNING:
                tut.status = FAILED
                tut.reason = tut.reason or "interrupted with Ctrl-C"

    # -- one tutorial ------------------------------------------------------- #
    def _run_tutorial(self, code: str) -> None:
        tutorial = self.catalog[code]
        record = self.state.tutorial(code, tutorial.title)
        blocker = self._blocked_by(tutorial)
        if blocker:
            record.status, record.blocked_by = BLOCKED, blocker
            record.reason = (f"{blocker} did not finish, and {code} starts from "
                             "what it produces")
            LOG.warning("%-9s BLOCKED by %s", code, blocker)
            self._save()
            return

        record.status, record.started = RUNNING, time.time()
        # A tutorial may ask for longer than the run-wide budget; OQMD does,
        # because its own searches take most of a minute on a good day.
        budget = float(tutorial.timeout
                       or getattr(self.options, "tutorial_timeout", 0) or 0)
        self._deadline.at = (time.time() + budget) if budget > 0 else None
        self._deadline.budget = budget
        try:
            workdir = seed_workdir(tutorial, self.options)
        except (OSError, FileNotFoundError) as exc:
            record.status, record.reason = FAILED, f"could not seed the work directory: {exc}"
            record.finished = time.time()
            LOG.error("%-9s FAILED while seeding: %s", code, exc)
            self._save()
            return
        record.workdir = str(workdir)
        if tutorial.stub:
            record.reason = "stub: cannot run as shipped -- " + tutorial.note
        LOG.info("%-9s %s", code, tutorial.title)

        failed = self._execute_steps(tutorial, record, workdir)
        # A tutorial whose failure is most likely someone else's server gets
        # another run at it, with a fresh budget.  Steps that already passed
        # are not repeated: the checkpoint says they are done, so a retry
        # picks up at the step that ran out of time.
        for attempt in range(2, max(1, tutorial.attempts) + 1):
            if failed is None or not self._out_of_time(failed):
                break
            LOG.warning("%-9s attempt %d of %d: %s", code, attempt,
                        tutorial.attempts, failed.reason)
            record.steps.pop(failed.key, None)
            self._deadline.at = (time.time() + budget) if budget > 0 else None
            self._deadline.budget = budget
            failed = self._execute_steps(tutorial, record, workdir, resume=True)
        record.finished = time.time()
        # A service that did not answer is not a defect anyone reading this
        # report can act on.  After every attempt has been spent waiting on
        # it, record the tutorial as skipped rather than failed -- the same
        # treatment an absent POTCAR or API key already gets -- so a red run
        # still means something is actually wrong.  Only a timeout is
        # forgiven: a wrong answer from OQMD is still a failure.
        gave_up_on_service = (failed is not None and tutorial.flaky_service
                              and self._out_of_time(failed))
        if gave_up_on_service:
            record.status = SKIPPED
            record.reason = (
                "{} did not respond within {:.0f}s, on {} attempt(s) -- "
                "skipped, not failed: the service is outside HTESP's control "
                "and its client has no request timeout".format(
                    tutorial.topic.upper(), budget, tutorial.attempts))
            LOG.warning("%-9s SKIPPED -- %s", code, record.reason)
        else:
            record.status = FAILED if failed else DONE
        if not failed and not tutorial.stub:
            record.reason = ""
        if failed and not gave_up_on_service:
            record.reason = f"stopped at step {failed.index}/{len(tutorial.steps)} " \
                            f"({failed.step_id}): {failed.reason}"
        self._save()

    @staticmethod
    def _out_of_time(step: StepState) -> bool:
        """Did this step fail because the clock ran out, not the command?"""
        reason = (step.reason or "").lower()
        return "budget" in reason or "timed out" in reason

    def _blocked_by(self, tutorial: Tutorial) -> str:
        """The dependency that stops this tutorial running, or ``""``.

        A dependency SKIPPED because its service did not answer is not a
        blocker.  `data-combine` merges what the three database tutorials
        downloaded, and OQMD being unreachable should not stop it: seeding
        falls back to that tutorial's reference, which holds the same files.
        Blocking there would fail a tutorial for a reason that is neither its
        own nor HTESP's.
        """
        for dep in tutorial.depends_on:
            if dep not in self.codes:
                continue
            record = self.state.tutorials.get(dep)
            if record is not None and record.status == DONE:
                continue
            if (record is not None and record.status == SKIPPED
                    and self.catalog[dep].flaky_service):
                LOG.info("%-9s %s was skipped (%s did not answer); its "
                         "reference will be used instead",
                         tutorial.code, dep, self.catalog[dep].topic.upper())
                continue
            return dep
        return ""

    def _execute_steps(self, tutorial: Tutorial, record: TutorialState,
                       workdir: Path, resume: bool = False) -> StepState | None:
        """Run the step list, honouring the convergence loop.  Returns the failure."""
        steps = tutorial.steps
        loop = tutorial.loop
        loop_first = min((i for i, s in enumerate(steps)
                          if loop and s.id in loop.steps), default=-1)
        loop_last = max((i for i, s in enumerate(steps)
                         if loop and s.id in loop.steps), default=-1)
        started = self.options.from_step is None
        index, cycle = 0, 1
        while index < len(steps):
            step = steps[index]
            if not started:
                if step.id == self.options.from_step:
                    started = True
                else:
                    self._record_skip(record, step, index + 1, cycle, workdir,
                                      f"before --from {self.options.from_step}")
                    index += 1
                    continue
            state = self._run_step(tutorial, record, step, index + 1, cycle,
                                   workdir, force_resume=resume)
            if state.status == FAILED:
                return state
            if loop and index == loop_last and cycle < loop.max_cycles:
                probe = PROBES.get(loop.probe, lambda _wd: None)
                if probe(workdir) is False:
                    cycle += 1
                    index = loop_first
                    LOG.info("%-9s relaxation not converged; loop cycle %d",
                             tutorial.code, cycle)
                    continue
            index += 1
        return None

    # -- one step ----------------------------------------------------------- #
    def _record_skip(self, record: TutorialState, step: Step, index: int, cycle: int,
                     workdir: Path, reason: str) -> StepState:
        key = step_key(step.id, cycle)
        state = StepState(step_id=step.id, key=key, index=index, cycle=cycle,
                          status=SKIPPED, reason=reason, workdir=str(workdir),
                          expected=list(step.artifacts))
        record.steps[key] = state
        self._save()
        LOG.info("    %s%2d. %-18s SKIPPED (%s)",
                 self._tag(record.code), index, step.id, reason)
        return state

    def _log_path(self, tutorial: Tutorial, step: Step, index: int, cycle: int) -> Path:
        folder = self.options.logs_root / tutorial.workdir_name
        folder.mkdir(parents=True, exist_ok=True)
        suffix = "" if cycle <= 1 else f"-cycle{cycle}"
        return folder / f"{index:02d}-{step.id}{suffix}.log"

    def _command(self, step: Step) -> list[str]:
        cmd = [self.options.python, "-m", "htesp", str(step.command), *step.args]
        if self.options.workers:
            cmd += ["--workers", str(self.options.workers)]
        # Always: this runner prepares inputs and never submits.
        cmd.append("--dry-run")
        if self.options.verbose:
            cmd.append("-v")
        return cmd

    def _run_step(self, tutorial: Tutorial, record: TutorialState, step: Step,
                  index: int, cycle: int, workdir: Path,
                  force_resume: bool = False) -> StepState:
        key = step_key(step.id, cycle)
        # force_resume: a retry must not repeat the steps that already passed,
        # even under --restart -- re-running a 39-second search to retry the
        # download after it would spend the new budget on work already done.
        if (self.options.resume or force_resume) and self.state.is_done(
                tutorial.code, key):
            existing = record.steps[key]
            LOG.info("    %s%2d. %-18s %s (resumed)",
                     self._tag(tutorial.code), index, step.id, existing.status)
            return existing
        left = self._time_left()
        if left is not None and left <= 0:
            state = self._record_skip(
                record, step, index, cycle, workdir,
                "the tutorial ran out of its {:.0f}s budget before this step"
                .format(getattr(self._deadline, "budget", 0) or 0))
            state.status = FAILED
            return state
        if step.needs_relax_output and not relax_output_present(workdir,
                                                                tutorial.dft):
            return self._record_skip(
                record, step, index, cycle, workdir,
                "needs a finished relaxation, and none is present")
        if step.needs_dft_output:
            return self._record_skip(
                record, step, index, cycle, workdir,
                "reads the output of a real DFT run, which this runner never "
                "performs" + self._how_to_run_it(tutorial))
        if step.needs_api_key and not mp_api_key():
            return self._record_skip(
                record, step, index, cycle, workdir,
                "no Materials Project API key is configured -- set one with "
                "'htesp-check --set_mp_api <your key>', or export MP_API_KEY")
        if step.needs_potcar and not potcars_available():
            return self._record_skip(
                record, step, index, cycle, workdir,
                "VASP POTCARs are not configured, so the inputs this step "
                "writes would be unusable -- run 'htesp-check "
                "--config_vasp_pot /path/to/POT_GGA_PAW_PBE'")
        upstream = self._skipped_upstream(record, step)
        if upstream:
            return self._record_skip(
                record, step, index, cycle, workdir,
                f"step '{upstream}', whose output it reads, was skipped")

        patch_input_in(workdir, step.input_patch, tutorial.dft)
        before = snapshot(workdir)
        log_path = self._log_path(tutorial, step, index, cycle)
        state = StepState(step_id=step.id, key=key, index=index, cycle=cycle,
                          status=RUNNING, workdir=str(workdir),
                          expected=list(step.artifacts), log=str(log_path),
                          started=time.time())
        record.steps[key] = state
        self._save()

        if step.is_callable:
            self._run_callable(step, state, workdir, log_path)
        else:
            state.command = self._command(step)
            self._run_subprocess(step, state, workdir, log_path)

        if state.status != FAILED:
            self._after_step(step, state, workdir, before)
        state.finished = time.time()
        state.duration = state.finished - (state.started or state.finished)
        self._save()
        mark = " "
        LOG.info("    %s%2d. %-18s %s%s %5.1fs", self._tag(tutorial.code),
                 index, step.id, state.status.upper(), mark, state.duration)
        return state

    @staticmethod
    def _skipped_upstream(record: TutorialState, step: Step) -> str:
        """The first step in ``step.after`` that did not actually run."""
        for other in step.after:
            state = record.steps.get(step_key(other))
            if state is not None and state.status == SKIPPED:
                return other
        return ""

    def _run_subprocess(self, step: Step, state: StepState, workdir: Path,
                        log_path: Path) -> None:
        timeout = step.timeout or self.options.step_timeout
        remaining = self._time_left()
        if remaining is not None:
            timeout = min(timeout, max(remaining, 1.0))
        header = [f"# command : {' '.join(state.command)}",
                  f"# cwd     : {workdir}", f"# started : {time.ctime()}", ""]
        # FIX: write the header *before* running, not after.  Every branch
        # below wrote the log only once the child had returned, so a step that
        # hung or was interrupted left an empty log directory and the stop
        # report had nothing to show -- which is the one thing this runner
        # exists to provide.  `htesp-tutorials` sat on QE/4 for 33 minutes and
        # logs/QE-4/ was empty afterwards.
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text("\n".join(header))
        try:
            # FIX: stdin=DEVNULL.  The child inherited the runner's stdin, so
            # anything that prompts blocks for ever with nobody watching.
            # qmpy_rester.get_oqmd_phases() defaults to verbose=True and calls
            # input('Proceed? [Y/n]:'); that is exactly what stalled QE/4.
            # With no stdin such a prompt fails immediately and is reported.
            proc = subprocess.run(state.command, cwd=os.fspath(workdir),
                                  capture_output=True, text=True, timeout=timeout,
                                  stdin=subprocess.DEVNULL, env=child_env())
        except KeyboardInterrupt:
            state.status, state.exit_code = FAILED, None
            state.reason = "interrupted with Ctrl-C while this step was running"
            log_path.write_text("\n".join(header) + f"\n# {state.reason}\n")
            raise
        except subprocess.TimeoutExpired as exc:
            state.status, state.exit_code = FAILED, None
            # the tutorial's own budget, not the run-wide default: OQMD asks
            # for 100s and quoting 1440s here was simply wrong
            budget = float(getattr(self._deadline, "budget", 0) or 0)
            if budget and timeout < (step.timeout or self.options.step_timeout):
                state.reason = ("killed after {:.0f}s: the tutorial's {:.0f}s "
                                "budget ran out while this step was running"
                                .format(timeout, budget))
            else:
                state.reason = f"timed out after {timeout:.0f}s"
            # FIX: `subprocess.run(text=True)` still hands TimeoutExpired
            # *bytes*, so this concatenation raised
            # "TypeError: can only concatenate str (not "bytes") to str"
            # and the timeout escaped as a traceback, leaving the tutorial
            # RUNNING and skipping the retry.  It never showed up while the
            # step limit was six hours, because nothing ever timed out.
            body = _as_text(exc.stdout) + _as_text(exc.stderr)
            log_path.write_text("\n".join(header) + body + f"\n# {state.reason}\n")
            return
        except OSError as exc:
            state.status, state.reason = FAILED, f"could not start the command: {exc}"
            log_path.write_text("\n".join(header) + f"\n# {state.reason}\n")
            return
        state.exit_code = proc.returncode
        log_path.write_text("\n".join(header) + proc.stdout +
                            ("\n--- stderr ---\n" + proc.stderr if proc.stderr else "") +
                            f"\n# exit code: {proc.returncode}\n")
        if proc.returncode != 0:
            state.status = FAILED
            state.reason = f"mainprogram {step.command} exited {proc.returncode}"
            return
        state.status = DONE

    def _run_callable(self, step: Step, state: StepState, workdir: Path,
                      log_path: Path) -> None:
        state.command = [step.describe()]
        try:
            code = int(step.command(workdir, step) or 0)   # type: ignore[operator]
        except Exception as exc:                            # noqa: BLE001 - reported
            state.status, state.exit_code = FAILED, None
            state.reason = f"{type(exc).__name__}: {exc}"
            log_path.write_text(f"# {step.describe()}\n# {state.reason}\n")
            return
        state.exit_code = code
        log_path.write_text(f"# {step.describe()}\n# exit code: {code}\n")
        state.status = DONE if code == 0 else FAILED
        if code:
            state.reason = f"{step.describe()} returned {code}"

    def _after_step(self, step: Step, state: StepState, workdir: Path,
                    before: dict | None = None) -> None:
        """Verify what the step left behind.

        This used to wait on the cluster first -- poll ``squeue`` until the
        step's jobs left the queue, ask ``sacct`` how each one ended, then
        read the output for convergence markers.  None of that has anything to
        do with a runner that never submits, and it went with real mode.
        """
        if before is not None:
            state.produced = changed_since(before, snapshot(workdir))
        state.missing = missing_artifacts(workdir, step.artifacts)
        if state.missing:
            state.status = FAILED
            state.reason = ("the command exited 0 but produced nothing: no match for "
                            + ", ".join(state.missing))
