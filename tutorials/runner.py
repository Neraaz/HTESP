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
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from tutorials.catalog import (CATALOG, EXAMPLES, PACKAGE_ROOT, Step,
                               Tutorial, searched_for_examples, waves)
from tutorials.state import (BLOCKED, DONE, FAILED, PENDING, RUNNING, SKIPPED,
                             RunState, StepState, TutorialState, step_key)
from tutorials.workdirs import (collect_job_ids, missing_artifacts,
                                patch_input_in, relaxation_converged,
                                seed_workdir, wait_for_jobs)

LOG = logging.getLogger("htesp.tutorials")

#: the three ways to run a tutorial
DRY_RUN, NO_DFT, REAL = "dry-run", "no-dft", "real"

#: default seconds a single ``mainprogram`` call may take
DEFAULT_STEP_TIMEOUT = 6 * 3600


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
    mode: str = DRY_RUN
    resume: bool = True
    poll_interval: float = 60.0
    job_timeout: float = 24 * 3600
    step_timeout: float = DEFAULT_STEP_TIMEOUT
    workers: int | None = None
    from_step: str | None = None
    #: which dependency wave to run: None runs the whole selection in one pass
    #: (what every earlier version did), an int runs exactly that wave, and
    #: "next" runs the lowest-numbered wave the checkpoint has not finished.
    wave: int | str | None = None
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
    if needs_key and not os.environ.get("MP_API_KEY"):
        level = "warning" if options.mode == DRY_RUN else "error"
        out.append(Problem(
            level,
            "MP_API_KEY is not set, and these tutorials query the Materials "
            f"Project: {', '.join(needs_key)}",
            "export MP_API_KEY=... (the key was removed from every config.json)"))

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

    if options.mode == REAL:
        for tool, why in (("sbatch", "submitting jobs"), ("squeue", "waiting for jobs")):
            if shutil.which(tool) is None:
                out.append(Problem("error", f"{tool} is not on PATH ({why})",
                                   "run with --dry-run or --no-dft on a laptop"))
        wanted = {"QE": "pw.x", "VASP": "vasp_std"}
        for dft in sorted(codes_used):
            if shutil.which(wanted[dft]) is None:
                out.append(Problem("error",
                                   f"{wanted[dft]} is not on PATH but {dft} "
                                   "tutorials were selected"))
    elif options.mode == NO_DFT and shutil.which("squeue") is None:
        out.append(Problem("warning", "squeue is not on PATH",
                           "--no-dft never submits anything, so this is harmless"))

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
        self.state.mode = options.mode

    # -- entry point -------------------------------------------------------- #
    def run(self) -> RunState:
        """Run every selected tutorial; always returns a saved checkpoint."""
        self.options.runs_root.mkdir(parents=True, exist_ok=True)
        self.options.logs_root.mkdir(parents=True, exist_ok=True)
        plan = waves(self.codes, self.catalog)
        number, codes = self._wave_to_run(plan)
        if codes is None:                      # nothing left to do
            self.state.save()
            return self.state
        record = self.state.wave(number, codes) if number is not None else None
        if record is not None:
            record.status, record.started = RUNNING, time.time()
            LOG.info("wave %d of %d: %d tutorial(s) -- %s",
                     number, len(plan) - 1, len(codes), ", ".join(codes))
            self.state.save()
        try:
            for code in codes:
                self._run_tutorial(code)
        except KeyboardInterrupt:
            self.state.interrupted = True
            self._mark_running_as_interrupted()
            LOG.error("interrupted -- the checkpoint is at %s", self.state.path)
        if record is not None:
            self._close_wave(record, plan)
        self.state.save()
        self._clean_work_dirs()
        return self.state

    # -- waves -------------------------------------------------------------- #
    def _wave_to_run(self, plan: list[list[str]]) -> tuple[int | None, list[str] | None]:
        """Which wave this invocation covers, and the codes in it.

        ``(None, every code)`` is the historical single-pass run; that is still
        the default, because in --dry-run and --no-dft nothing is ever queued
        and there is nothing to wait between waves *for*.
        """
        wanted = self.options.wave
        if wanted is None:
            return None, list(self.codes)
        if wanted == "next":
            number = self.state.next_wave(plan)
            if number is None:
                LOG.info("every wave is done (%d of %d); nothing to run",
                         len(plan), len(plan))
                return None, None
        else:
            number = int(wanted)
            if not 0 <= number < len(plan):
                LOG.error("wave %d does not exist: this selection has waves 0-%d",
                          number, len(plan) - 1)
                return None, None
        return number, list(plan[number])

    def _close_wave(self, record, plan: list[list[str]]) -> None:
        """Record how the wave ended, and say what to run next.

        A wave is DONE only when every tutorial in it is; anything else leaves
        it FAILED, so ``--wave next`` offers it again rather than stepping over
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
            LOG.error("wave %d did not finish (%s); fix those, then re-run "
                      "--wave %d", record.number, record.reason, record.number)
        elif following < len(plan):
            LOG.info("wave %d done.  Wait for its jobs to finish, then run "
                     "'htesp-tutorials --resume --wave next' for wave %d (%s)",
                     record.number, following, ", ".join(plan[following]))
        else:
            LOG.info("wave %d done -- that was the last one.", record.number)

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
            self.state.save()
            return

        record.status, record.started = RUNNING, time.time()
        try:
            workdir = seed_workdir(tutorial, self.options)
        except (OSError, FileNotFoundError) as exc:
            record.status, record.reason = FAILED, f"could not seed the work directory: {exc}"
            record.finished = time.time()
            LOG.error("%-9s FAILED while seeding: %s", code, exc)
            self.state.save()
            return
        record.workdir = str(workdir)
        if tutorial.stub:
            record.reason = "stub: cannot run as shipped -- " + tutorial.note
        LOG.info("%-9s %s", code, tutorial.title)

        failed = self._execute_steps(tutorial, record, workdir)
        record.finished = time.time()
        record.status = FAILED if failed else DONE
        if not failed and not tutorial.stub:
            record.reason = ""
        if failed:
            record.reason = f"stopped at step {failed.index}/{len(tutorial.steps)} " \
                            f"({failed.step_id}): {failed.reason}"
        self.state.save()

    def _blocked_by(self, tutorial: Tutorial) -> str:
        for dep in tutorial.depends_on:
            if dep not in self.codes:
                continue
            record = self.state.tutorials.get(dep)
            if record is None or record.status not in (DONE,):
                return dep
        return ""

    def _execute_steps(self, tutorial: Tutorial, record: TutorialState,
                       workdir: Path) -> StepState | None:
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
            state = self._run_step(tutorial, record, step, index + 1, cycle, workdir)
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
        self.state.save()
        LOG.info("    %2d. %-18s SKIPPED (%s)", index, step.id, reason)
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
        if self.options.mode == DRY_RUN or (self.options.mode == NO_DFT and step.submits):
            cmd.append("--dry-run")
        if self.options.verbose:
            cmd.append("-v")
        return cmd

    def _run_step(self, tutorial: Tutorial, record: TutorialState, step: Step,
                  index: int, cycle: int, workdir: Path) -> StepState:
        key = step_key(step.id, cycle)
        if self.options.resume and self.state.is_done(tutorial.code, key):
            existing = record.steps[key]
            LOG.info("    %2d. %-18s %s (resumed)", index, step.id, existing.status)
            return existing
        if step.needs_dft_output and self.options.mode == DRY_RUN:
            return self._record_skip(record, step, index, cycle, workdir,
                                     "reads the output of a real DFT run, which "
                                     "--dry-run never produces")
        if step.needs_api_key and not os.environ.get("MP_API_KEY"):
            return self._record_skip(record, step, index, cycle, workdir,
                                     "MP_API_KEY is not set")
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
        log_path = self._log_path(tutorial, step, index, cycle)
        state = StepState(step_id=step.id, key=key, index=index, cycle=cycle,
                          status=RUNNING, workdir=str(workdir),
                          expected=list(step.artifacts), log=str(log_path),
                          started=time.time())
        record.steps[key] = state
        self.state.save()

        if step.is_callable:
            self._run_callable(step, state, workdir, log_path)
        else:
            state.command = self._command(step)
            self._run_subprocess(step, state, workdir, log_path)

        if state.status != FAILED:
            self._after_step(step, state, workdir)
        state.finished = time.time()
        state.duration = state.finished - (state.started or state.finished)
        self.state.save()
        mark = "!" if state.unverifiable else " "
        LOG.info("    %2d. %-18s %s%s %5.1fs", index, step.id,
                 state.status.upper(), mark, state.duration)
        return state

    @staticmethod
    def _skipped_upstream(record: TutorialState, step: Step) -> str:
        """The first step in ``step.after`` that did not actually run."""
        for other in step.after:
            state = record.steps.get(step_key(other))
            if state is not None and state.status == SKIPPED:
                return other
        return ""

    def _unconverged(self, step: Step, workdir: Path) -> str:
        """``""`` when *step*'s relaxations converged, else what is wrong.

        A DFT code can exit 0, and its job can be COMPLETED, without the
        structure having relaxed: QE stops at ``nstep`` and VASP at ``NSW``,
        both leaving a full set of output files that satisfy an artefact glob.
        Everything downstream then runs on a structure that is not relaxed, and
        produces numbers that look entirely reasonable.  The strings checked
        here are the same ones htesp/workflow.py uses to decide whether a
        relaxation is finished.
        """
        if not step.check_converged or self.options.mode != REAL:
            return ""
        problems = []
        for pattern in step.job_dirs or ():
            for directory in sorted(workdir.glob(pattern)):
                verdict = relaxation_converged(directory)
                if verdict is False:
                    problems.append(directory.relative_to(workdir).as_posix())
        if not problems:
            return ""
        return ("the run finished but did not converge in "
                + ", ".join(problems[:5])
                + (f" (and {len(problems) - 5} more)" if len(problems) > 5 else "")
                + " -- everything downstream would start from an unrelaxed "
                  "structure")

    def _run_subprocess(self, step: Step, state: StepState, workdir: Path,
                        log_path: Path) -> None:
        timeout = step.timeout or self.options.step_timeout
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
            state.reason = f"timed out after {timeout:.0f}s"
            body = (exc.stdout or "") + (exc.stderr or "")
            log_path.write_text("\n".join(header) + str(body) + f"\n# {state.reason}\n")
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
        if self.options.mode == NO_DFT and step.submits:
            # mainprogram logs through the logging module, i.e. to stderr
            merged = (proc.stdout + "\n" + proc.stderr).splitlines()
            state.would_submit = [line.strip() for line in merged
                                  if "[dry-run]" in line]
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

    def _after_step(self, step: Step, state: StepState, workdir: Path) -> None:
        """Wait for the cluster, then verify the declared artefacts."""
        if step.submits and self.options.mode == REAL:
            state.jobs = collect_job_ids(workdir, step.job_dirs, state.started or 0.0)
            if not state.jobs:
                state.unverifiable = True
                state.reason = ("no job ids were recorded in "
                                f"{', '.join(step.job_dirs) or '<no stage dirs>'}"
                                "/.htesp_job.json -- nothing appears to have been "
                                "submitted")
            else:
                ok, problem = wait_for_jobs(state.jobs, self.options.poll_interval,
                                            self.options.job_timeout)
                if not ok and "timed out" in problem:
                    state.status, state.reason = FAILED, problem
                    return
                if not ok:
                    state.unverifiable, state.reason = True, problem
                elif problem:
                    # the jobs finished, and sacct says at least one did not
                    # finish *well*.  That is a failure, not an unverifiable:
                    # the scheduler said so plainly.
                    state.status, state.reason = FAILED, problem
                    return
            unconverged = self._unconverged(step, workdir)
            if unconverged:
                state.status, state.reason = FAILED, unconverged
                return
        state.missing = missing_artifacts(workdir, step.artifacts)
        if state.missing:
            state.status = FAILED
            state.reason = ("the command exited 0 but produced nothing: no match for "
                            + ", ".join(state.missing))
