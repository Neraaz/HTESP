#!/usr/bin/env python
"""Work directories, ``input.in`` edits, cluster-job waiting and verification.

``examples/`` is read-only input to the tutorial runner, so every tutorial gets
its own directory under ``tutorial_runs/`` seeded with the files it needs: the
shared example configuration, whatever the tutorials it depends on produced,
and its own shipped files.  This module does that seeding, keeps ``input.in``
in step with the catalogue, waits on the jobs the workflow layer recorded in
``<stage dir>/.htesp_job.json``, and answers the question a step is judged by:
did the artefacts it promised actually appear?

Split out of :mod:`tutorials.runner` so that seeding and verification can be
tested without constructing a runner.
"""
from __future__ import annotations

import fnmatch
import logging
import os
import shutil
import subprocess
import tarfile
import time
from pathlib import Path
from typing import TYPE_CHECKING, Iterable, Sequence

from tutorials.catalog import SEED_EXCLUDE, Seed, Tutorial

if TYPE_CHECKING:                       # pragma: no cover - typing only
    from tutorials.runner import RunOptions

LOG = logging.getLogger("htesp.tutorials")


# --------------------------------------------------------------------------- #
#  seeding
# --------------------------------------------------------------------------- #
def _excluded(name: str) -> bool:
    return any(fnmatch.fnmatch(name, pattern) for pattern in SEED_EXCLUDE)


def _place(src: Path, dst: Path, overwrite: bool, link: bool) -> bool:
    """Copy (or symlink) *src* to *dst*; returns True when something was placed."""
    if dst.exists() and not overwrite:
        return False
    dst.parent.mkdir(parents=True, exist_ok=True)
    if link and src.is_dir():
        if dst.is_symlink() or dst.exists():
            return False
        dst.symlink_to(src.resolve(), target_is_directory=True)
        return True
    if src.is_dir():
        shutil.copytree(src, dst, dirs_exist_ok=True, symlinks=True)
        return True
    if dst.is_dir():
        shutil.rmtree(dst) if overwrite else None
    shutil.copy2(src, dst)
    return True


def _seed_source(seed: Seed, tutorial: Tutorial, options: "RunOptions") -> Path | None:
    if seed.source == "code":
        return options.examples / tutorial.dft
    if seed.source == "self":
        return tutorial.directory
    if seed.source == "archive":
        return None
    return options.runs_root / seed.source.replace("/", "-")


def _seed_from_archive(seed: Seed, tutorial: Tutorial, workdir: Path,
                       overwrite: bool | None = None) -> list[str]:
    """Flatten matching members of ``reference*.tar.gz`` into the work directory."""
    placed: list[str] = []
    if overwrite is None:
        overwrite = seed.overwrite
    for archive in sorted(tutorial.directory.glob("reference*.tar.gz")):
        try:
            with tarfile.open(archive) as tar:
                for member in tar.getmembers():
                    if not member.isfile():
                        continue
                    name = Path(member.name).name
                    if not any(fnmatch.fnmatch(name, p) for p in seed.patterns):
                        continue
                    target = workdir / name
                    if target.exists() and not overwrite:
                        continue
                    handle = tar.extractfile(member)
                    if handle is None:
                        continue
                    target.write_bytes(handle.read())
                    placed.append(name)
        except (OSError, tarfile.TarError) as exc:
            LOG.warning("%s: could not read %s (%s)", tutorial.code, archive.name, exc)

    # FIX: some tutorials ship their reference material *unpacked* -- QE/8 has
    # reference.tar.gz but VASP/8 has a plain reference/ directory, and the
    # .cif files this seed exists to provide live inside it.  Globbing only
    # for the tarball left the work directory with no .cif at all, so
    # `mainprogram download` in fromcif mode exited 0 having written nothing
    # and the runner reported "produced no artifacts".
    for folder in sorted(tutorial.directory.glob("reference*")):
        if not folder.is_dir():
            continue
        for member in sorted(folder.rglob("*")):
            if not member.is_file():
                continue
            name = member.name
            if not any(fnmatch.fnmatch(name, p) for p in seed.patterns):
                continue
            target = workdir / name
            if target.exists() and not overwrite:
                continue
            try:
                shutil.copy2(member, target)
            except OSError as exc:
                LOG.warning("%s: could not copy %s (%s)", tutorial.code, member, exc)
                continue
            placed.append(name)
    return placed


def seed_workdir(tutorial: Tutorial, options: "RunOptions") -> Path:
    """Create and populate ``tutorial_runs/<code>/`` and return it.

    Seeds are applied in catalogue order: the shared example configuration
    first, then whatever the tutorials this one depends on produced, then the
    tutorial's own shipped files (which therefore win for ``config.json`` and
    friends).  Documentation and ``reference*`` archives are never copied --
    those are the expected answer, not input.
    """
    workdir = options.runs_root / tutorial.workdir_name
    # FIX: seeding must not clobber what an earlier run produced.  Seeds
    # normally overwrite, so that a fresh work directory gets the tutorial's
    # shipped files -- but on --resume (and --from, and --only against an
    # existing tree) the steps that would have regenerated those files are
    # skipped, and the shipped copy wins instead.  examples/QE/tutorial6 ships
    # an mpid-list.in full of *OQMD* ids left over from the author's own run,
    # so resuming QE/6 fed OQMD ids to a Materials Project lookup and the
    # tutorial "failed" for a reason that had nothing to do with the code.
    # When the directory already exists and we are resuming, seed only what is
    # missing.  --restart removes it first, so a fresh run is unaffected.
    resuming = workdir.is_dir() and getattr(options, "resume", False)
    workdir.mkdir(parents=True, exist_ok=True)
    for seed in tutorial.seeds:
        if seed.source == "archive":
            _seed_from_archive(seed, tutorial, workdir, overwrite=not resuming)
            continue
        source = _seed_source(seed, tutorial, options)
        if source is None or not source.is_dir():
            if seed.required:
                raise FileNotFoundError(
                    f"{tutorial.code}: required seed source {source} is missing")
            continue
        for pattern in seed.patterns:
            for path in sorted(source.glob(pattern)):
                if seed.source == "self" and _excluded(path.name):
                    continue
                name = dict(seed.rename).get(path.name, path.name)
                overwrite = seed.overwrite and not resuming
                _place(path, workdir / name, overwrite, seed.link)
    _ensure_config(tutorial, workdir, options)
    _ensure_batch_header(tutorial, workdir, options)
    return workdir


def _ensure_config(tutorial: Tutorial, workdir: Path, options: "RunOptions") -> None:
    """Some tutorials ship ``config.json_search`` &c. but no plain ``config.json``."""
    if (workdir / "config.json").is_file():
        return
    for candidate in sorted(tutorial.directory.glob("config.json_*")):
        shutil.copy2(candidate, workdir / "config.json")
        LOG.info("%s: used %s as config.json", tutorial.code, candidate.name)
        return


def _ensure_batch_header(tutorial: Tutorial, workdir: Path,
                         options: "RunOptions") -> None:
    """Write a ``batch.header`` this cluster will actually accept.

    The header seeded from ``examples/<code>/`` says ``--partition=dense`` and
    loads no module: it was written for one machine and exists nowhere else,
    so every submitting step of a real run would be rejected by ``sbatch``
    before any calculation started.  Where SLURM is present, replace it with
    one built from what this machine reports -- partition, account, cores per
    node, and the ``qe``/``vasp`` module Lmod offers.

    Nothing is generated on a machine without ``sinfo``: there the probes have
    nothing to say, the generated header would be all ``TODO``, and a laptop
    ``--dry-run`` submits nothing anyway, so the shipped header is the more
    useful thing to leave in place.

    The launcher is written into ``config.json`` rather than the header,
    because ``mainprogram jobscript`` builds the run line from
    ``job_script.parallel_command`` -- see :mod:`htesp.generate_submission`.
    """
    if shutil.which("sinfo") is None:
        return
    from htesp import batch_header

    target = workdir / "batch.header"
    try:
        text = batch_header.build(tutorial.dft)
    except (ValueError, OSError) as exc:            # pragma: no cover
        LOG.warning("%s: could not build batch.header (%s); keeping the "
                    "shipped one", tutorial.code, exc)
        return
    if target.is_file() and target.read_text() == text:
        return
    target.write_text(text)
    LOG.debug("%s: wrote batch.header for this cluster", tutorial.code)
    _match_launcher(workdir)


def _match_launcher(workdir: Path) -> None:
    """Point ``job_script.parallel_command`` at this machine's launcher.

    A header that loads the right module is still useless if the run line says
    ``mpirun`` on a machine whose MPI is driven by ``ibrun`` or ``srun``.
    ``nproc`` is left alone: it is the study's choice, not the machine's.
    """
    import json

    from htesp import batch_header

    found = batch_header.launcher()
    path = workdir / "config.json"
    if not found or not path.is_file():
        return
    try:
        data = json.loads(path.read_text())
    except (OSError, ValueError):                   # pragma: no cover
        return
    job = data.get("job_script")
    if not isinstance(job, dict) or job.get("parallel_command") == found:
        return
    job["parallel_command"] = found
    try:
        path.write_text(json.dumps(data, indent=4))
    except OSError:                                 # pragma: no cover
        pass


# --------------------------------------------------------------------------- #
#  input.in
# --------------------------------------------------------------------------- #
#: ``input.in`` is six lines: start, end, "nkpt kcut", track file, plot types, DFT
INPUT_IN_DEFAULT = ("1", "2", "200 0", "mpid.in", "phband", "DFT = QE")

#: track files to fall back on, in order, when the one named does not exist
TRACK_FALLBACKS = ("mpid.in", "mpid-list.in")


def patch_input_in(workdir: Path, patch, dft: str = "QE") -> list[str]:
    """Apply an :class:`~tutorials.catalog.InputPatch` to ``<workdir>/input.in``."""
    path = workdir / "input.in"
    try:
        lines = path.read_text().splitlines()
    except OSError:
        lines = []
    lines = list(lines) + list(INPUT_IN_DEFAULT[len(lines):])
    lines[5] = f"DFT = {dft}"
    if patch is not None:
        if patch.start is not None:
            lines[0] = str(patch.start)
        if patch.end is not None:
            lines[1] = str(patch.end)
        if patch.track is not None:
            lines[3] = patch.track
        if patch.plot is not None:
            lines[4] = patch.plot
    # The shared examples/<code>/input.in names mpid-list.in, but a tutorial
    # seeded from the relaxation hub only has mpid.in.  Naming a track file that
    # is not there makes mainprogram exit 2 before it does anything, which is a
    # seeding artefact rather than a tutorial failure -- so fall back.
    if not (workdir / lines[3]).is_file():
        for candidate in TRACK_FALLBACKS:
            if (workdir / candidate).is_file():
                LOG.debug("input.in: %s is missing, using %s", lines[3], candidate)
                lines[3] = candidate
                break
    path.write_text("\n".join(lines) + "\n")
    return lines


# --------------------------------------------------------------------------- #
#  cluster jobs
# --------------------------------------------------------------------------- #
def collect_job_ids(workdir: Path, patterns: Iterable[str],
                    since: float = 0.0) -> list[str]:
    """Job ids the workflow layer recorded in ``<stage dir>/.htesp_job.json``.

    The format is ``{tag: [{"job": "<id>", "time": <epoch>}, ...]}``; only ids
    recorded at or after *since* count, so a rerun does not wait for the jobs of
    the previous attempt.
    """
    import json

    found: list[str] = []
    seen: set[str] = set()
    for pattern in patterns:
        for stage in sorted(workdir.glob(pattern)):
            store = stage / ".htesp_job.json"
            try:
                data = json.loads(store.read_text())
            except (OSError, ValueError):
                continue
            for runs in data.values():
                for entry in runs:
                    job = str(entry.get("job", "")).strip()
                    if job and job not in seen and float(entry.get("time", 0)) >= since:
                        seen.add(job)
                        found.append(job)
    return found


#: `sacct` states that mean the job ran to the end.  COMPLETED is the only
#: unambiguous success; the rest are ways of not finishing, and each has been
#: seen to leave partial output that satisfies an artefact glob.
JOB_STATE_OK = ("COMPLETED",)

#: what each failure state means, in the words the report should use
JOB_STATE_MEANING = {
    "TIMEOUT": "hit the wall time",
    "OUT_OF_MEMORY": "ran out of memory",
    "CANCELLED": "was cancelled",
    "NODE_FAIL": "lost its node",
    "PREEMPTED": "was preempted",
    "FAILED": "exited non-zero",
}


def job_states(jobs: Sequence[str]) -> dict[str, str]:
    """``sacct`` state per job id, or ``{}`` when sacct cannot answer.

    Leaving the queue is not succeeding: COMPLETED, FAILED, TIMEOUT, CANCELLED
    and OUT_OF_MEMORY all make a job id disappear from ``squeue`` alike.  A
    relaxation killed at the wall time still leaves an OUTCAR and a CONTCAR
    behind, so the artefact check passes and everything downstream then runs on
    a structure that was never relaxed.
    """
    if not jobs or shutil.which("sacct") is None:
        return {}
    proc = subprocess.run(
        ["sacct", "-n", "-X", "-P", "-o", "JobID,State", "-j", ",".join(jobs)],
        capture_output=True, text=True, check=False)
    out = {}
    for line in proc.stdout.splitlines():
        parts = line.split("|")
        if len(parts) < 2:
            continue
        # "CANCELLED by 12345" -- the reason is not part of the state
        out[parts[0].strip()] = parts[1].strip().split()[0] if parts[1].strip() else ""
    return out


def check_job_states(jobs: Sequence[str]) -> str:
    """``""`` when every job completed, else what went wrong.

    An empty result also covers "sacct could not tell us": accounting is not
    enabled everywhere, and a machine without it is not a machine where every
    job failed.  :func:`wait_for_jobs` reports that case as unverifiable.
    """
    states = job_states(jobs)
    if not states:
        return ""
    bad = []
    for job in jobs:
        state = states.get(str(job))
        if state is None or state in JOB_STATE_OK:
            continue
        meaning = JOB_STATE_MEANING.get(state)
        bad.append(f"{job} {state}" + (f" ({meaning})" if meaning else ""))
    if not bad:
        return ""
    return "the scheduler reports " + ", ".join(bad)


def wait_for_jobs(jobs: Sequence[str], poll: float, timeout: float,
                  sleep=time.sleep) -> tuple[bool, str]:
    """Poll ``squeue`` until none of *jobs* is queued or running, then ask
    ``sacct`` how each one ended.

    Returns ``(finished, problem)``.  ``problem`` is non-empty when the wait
    could not be carried out -- no ``squeue`` on PATH, or the timeout expired --
    in which case the caller must treat the step as *unverifiable* rather than
    quietly successful.  A job that finished badly is a failure, not an
    unverifiable: the scheduler said so plainly.
    """
    if not jobs:
        return True, ""
    if shutil.which("squeue") is None:
        return False, ("squeue is not on PATH, so the driver cannot tell whether "
                       f"jobs {', '.join(jobs)} finished")
    deadline = time.time() + timeout
    while True:
        proc = subprocess.run(["squeue", "-h", "-o", "%i", "-j", ",".join(jobs)],
                              capture_output=True, text=True, check=False)
        if not proc.stdout.strip():
            return True, check_job_states(jobs)
        if time.time() > deadline:
            still = " ".join(proc.stdout.split())
            return False, (f"timed out after {timeout:.0f}s waiting for jobs "
                           f"{still or ', '.join(jobs)}")
        sleep(poll)


# --------------------------------------------------------------------------- #
#  verification
# --------------------------------------------------------------------------- #
#: what "this relaxation finished" looks like in each code's output.  These are
#: the same markers htesp/workflow.py greps for (see `_further_relax_input_one`
#: and `_phonopy_energies`), kept here as one table rather than a second copy of
#: the logic that could drift from it.
CONVERGED_MARKERS = {
    # VASP prints this once the ionic loop has met EDIFFG
    "OUTCAR": "reached required accuracy - stopping structural energy minimisation",
    # QE prints "JOB DONE." on any clean exit and the BFGS line only on success
    "relax.out": "End of BFGS Geometry Optimization",
    "scf.out": "JOB DONE.",
}

#: a relaxation that stopped because it ran out of steps, per code.  Finding one
#: of these is proof of non-convergence even when a success marker is absent for
#: an unrelated reason.
EXHAUSTED_MARKERS = (
    "The maximum number of steps has been reached",     # QE
    "convergence NOT achieved",                         # QE scf
)


def relaxation_converged(directory: Path) -> bool | None:
    """Did the relaxation in *directory* converge?

    Returns ``True``, ``False``, or ``None`` when there is nothing to judge by
    -- no recognised output file, which is not the same as a failure and is
    left to the artefact check to report.

    A job that COMPLETED is not a relaxation that converged: both codes stop
    at their ionic-step limit and exit cleanly, leaving CONTCAR/OUTCAR or a
    full ``relax.out`` behind.  The artefact globs match either way, so without
    this every dependent tutorial would start from an unrelaxed cell.
    """
    verdict = None
    for name, marker in CONVERGED_MARKERS.items():
        for path in sorted(directory.glob(name)):
            try:
                text = path.read_text(errors="replace")
            except OSError:                             # pragma: no cover
                continue
            if marker in text:
                return True
            verdict = False
            if any(bad in text for bad in EXHAUSTED_MARKERS):
                return False
    return verdict



def missing_artifacts(workdir: Path, patterns: Sequence[str]) -> list[str]:
    """Which of *patterns* match nothing under *workdir*."""
    return [pattern for pattern in patterns if not list(workdir.glob(pattern))]
