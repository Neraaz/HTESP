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
    if seed.source in ("archive", "reference-output"):
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
        if seed.source == "reference-output":
            _seed_reference_output(seed, tutorial, workdir)
            continue
        source = _seed_source(seed, tutorial, options)
        if source is None or not source.is_dir():
            if seed.required:
                raise FileNotFoundError(
                    f"{tutorial.code}: required seed source {source} is missing")
            continue
        placed = 0
        for pattern in seed.patterns:
            for path in sorted(source.glob(pattern)):
                if seed.source == "self" and _excluded(path.name):
                    continue
                name = dict(seed.rename).get(path.name, path.name)
                overwrite = seed.overwrite and not resuming
                _place(path, workdir / name, overwrite, seed.link)
                placed += 1
        if not placed and "/" in seed.source:
            # The tutorial this one seeds from produced nothing -- OQMD was
            # skipped because the service did not answer, say.  Its reference
            # holds the same files, so use those rather than blocking a
            # tutorial whose own subject is combining the three databases.
            _seed_from_dependency_reference(seed, workdir,
                                            overwrite=not resuming)
    _ensure_config(tutorial, workdir, options)
    _ensure_batch_header(tutorial, workdir, options)
    return workdir


def _seed_from_dependency_reference(seed: Seed, workdir: Path,
                                    overwrite: bool = True) -> list[str]:
    """Seed from another tutorial's *reference* when its work directory is bare.

    ``data-combine`` merges what the Materials Project, OQMD and AFLOW
    tutorials each downloaded, so it seeds from all three work directories.
    When one of them did not run -- OQMD is skipped when the service does not
    answer -- that directory is empty and the combine tutorial would be
    blocked by a failure that is not its own and not HTESP's.

    The missing tutorial's ``reference*.tar.gz`` contains the same files it
    would have written (``mpid-list.in``, ``mpid.in``, ``scf_dir/``, the
    ``R<id>-<compound>/`` directories), so they are used instead.  This runs
    *only* when the work directory yielded nothing: a tutorial that really
    ran always wins over its reference.
    """
    from tutorials.catalog import CATALOG

    tutorial = CATALOG.get(seed.source)
    if tutorial is None:
        return []
    placed: list[str] = []

    def wanted(parts: tuple) -> Path | None:
        """Map ``reference/scf_dir/x.in`` to ``<workdir>/scf_dir/x.in``."""
        if len(parts) < 2 or not parts[0].startswith("reference"):
            return None
        rest = parts[1:]
        if not any(fnmatch.fnmatch(rest[0], pattern)
                   for pattern in seed.patterns):
            return None
        return workdir.joinpath(*rest)

    for archive in sorted(tutorial.directory.glob("reference*.tar.gz")):
        try:
            with tarfile.open(archive) as tar:
                for member in tar.getmembers():
                    if not member.isfile():
                        continue
                    target = wanted(Path(member.name).parts)
                    if target is None or (target.exists() and not overwrite):
                        continue
                    handle = tar.extractfile(member)
                    if handle is None:
                        continue
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(handle.read())
                    placed.append(target.relative_to(workdir).as_posix())
        except (OSError, tarfile.TarError) as exc:      # pragma: no cover
            LOG.warning("could not read %s (%s)", archive.name, exc)

    for folder in sorted(tutorial.directory.glob("reference*")):
        if not folder.is_dir():
            continue
        for path in sorted(folder.rglob("*")):
            if not path.is_file():
                continue
            target = wanted(path.relative_to(folder.parent).parts)
            if target is None or (target.exists() and not overwrite):
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(path.read_bytes())
            placed.append(target.relative_to(workdir).as_posix())

    if placed:
        LOG.info("%s produced nothing; seeded %d file(s) from its reference "
                 "instead", seed.source, len(placed))
    return placed


def _seed_reference_output(seed: Seed, tutorial: Tutorial,
                           workdir: Path) -> list[str]:
    """Place the reference relaxation output into ``R<mpid>-<compound>/relax/``.

    Seventeen steps read a finished relaxation and nothing else -- the total
    energy, the relaxed structure, the cell to build the next input from.
    They have always been skipped because no DFT runs here, yet the answer is
    sitting in the tutorial's own reference.  Putting just that output back
    lets them run for real.

    Only the files in :data:`~tutorials.catalog.REFERENCE_OUTPUT` are taken,
    and only from a ``relax/`` directory.  ``econv.csv`` and
    ``scf_dir/scf-relax-*.in`` -- the things those steps must *produce* --
    stay out, or a step would pass by finding an artefact it never wrote.

    Never overwrites: a resumed run, or one where a step has already written
    its own output, keeps what is there.
    """
    wanted = set(seed.patterns)
    placed: list[str] = []

    def take(parts: tuple, read) -> None:
        """*parts* is the reference path; place it under the same R*/relax/.

        The destination directory cannot be found by globbing the work
        directory: VASP/9 ships an ``R<mpid>-<compound>/`` but QE/9 does not
        -- there the directory is created by ``mainprogram 1``, which runs
        after seeding.  The reference path carries the name, so use it.
        """
        name = parts[-1]
        if name not in wanted or "relax" not in parts:
            return
        index = len(parts) - 1 - parts[::-1].index("relax")
        if index == 0:
            return                                 # no R* component to anchor to
        target = workdir / parts[index - 1] / "relax" / name
        if target.exists():
            return
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(read())
        placed.append(target.relative_to(workdir).as_posix())

    for archive in sorted(tutorial.directory.glob("reference*.tar.gz")):
        try:
            with tarfile.open(archive) as tar:
                for member in tar.getmembers():
                    if not member.isfile():
                        continue
                    handle = tar.extractfile(member)
                    if handle is None:
                        continue
                    data = handle.read()
                    take(Path(member.name).parts, lambda d=data: d)
        except (OSError, tarfile.TarError) as exc:      # pragma: no cover
            LOG.warning("%s: could not read %s (%s)",
                        tutorial.code, archive.name, exc)

    for folder in sorted(tutorial.directory.glob("reference*")):
        if not folder.is_dir():
            continue
        for path in sorted(folder.rglob("*")):
            if path.is_file():
                take(path.relative_to(folder.parent).parts,
                     lambda p=path: p.read_bytes())

    if placed:
        LOG.info("%s: seeded reference relaxation output (%s)",
                 tutorial.code, ", ".join(sorted({Path(p).name for p in placed})))
    return placed


def snapshot(workdir: Path) -> dict:
    """``{relative path: (size, mtime)}`` for every file under *workdir*.

    Symlinked directories are not followed: QE work directories link
    ``pp/`` at the shared pseudopotential tree, and walking it would add
    hundreds of files that no step ever writes.
    """
    found: dict = {}
    stack = [workdir]
    while stack:
        current = stack.pop()
        try:
            entries = list(os.scandir(current))
        except OSError:                              # pragma: no cover
            continue
        for entry in entries:
            if entry.is_symlink():
                continue
            if entry.is_dir():
                stack.append(Path(entry.path))
                continue
            try:
                stat = entry.stat()
            except OSError:                          # pragma: no cover
                continue
            key = Path(entry.path).relative_to(workdir).as_posix()
            found[key] = (stat.st_size, stat.st_mtime_ns)
    return found


def changed_since(before: dict, after: dict) -> list[str]:
    """Paths that appeared, or whose size or mtime moved."""
    return sorted(name for name, stamp in after.items()
                  if before.get(name) != stamp)


def relax_output_present(workdir: Path, dft: str) -> bool:
    """Is there a finished relaxation under ``workdir`` to read?

    What the seeded reference provides, and what a real run would leave
    behind, are the same files -- so a step that needs only a relaxation can
    ask this instead of asking which mode it is running in.
    """
    from tutorials.catalog import REFERENCE_OUTPUT

    names = REFERENCE_OUTPUT.get(dft.upper(), ())
    return any(list(workdir.glob(f"R*-*/relax/{name}")) for name in names)


def _ensure_config(tutorial: Tutorial, workdir: Path, options: "RunOptions") -> None:
    """Some tutorials ship ``config.json_search`` &c. but no plain ``config.json``."""
    if (workdir / "config.json").is_file():
        return
    for candidate in sorted(tutorial.directory.glob("config.json_*")):
        shutil.copy2(candidate, workdir / "config.json")
        LOG.info("%s: used %s as config.json", tutorial.code, candidate.name)
        return


#: built headers, keyed by code.  The header depends on the *cluster* and the
#: code, never on the tutorial, but seeding calls this once per work directory
#: -- 42 times in a full sweep, at ~0.75 s each, because every call shells out
#: to sinfo, sacctmgr, scontrol and Lmod avail/help/spider.  Building once per
#: code turns 30 seconds of probing into 1.5.
_HEADER_CACHE: dict = {}


def _header_for(dft: str, code: str = "") -> str | None:
    """``batch.header`` text for this code, built at most once per run."""
    if dft in _HEADER_CACHE:
        return _HEADER_CACHE[dft]
    from htesp import batch_header

    try:
        text = batch_header.build(dft)
    except (ValueError, OSError) as exc:            # pragma: no cover
        LOG.warning("%s: could not build batch.header (%s); keeping the "
                    "shipped one", code or dft, exc)
        text = None
    _HEADER_CACHE[dft] = text
    return text


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
    this runner submits nothing anyway, so the shipped header is the more
    useful thing to leave in place.

    The launcher is written into ``config.json`` rather than the header,
    because ``mainprogram jobscript`` builds the run line from
    ``job_script.parallel_command`` -- see :mod:`htesp.generate_submission`.
    """
    if shutil.which("sinfo") is None:
        return
    target = workdir / "batch.header"
    text = _header_for(tutorial.dft, tutorial.code)
    if text is None:
        return
    if target.is_file() and target.read_text() == text:
        return
    target.write_text(text)
    LOG.debug("%s: wrote batch.header for this cluster", tutorial.code)
    _warn_about_headers()
    _match_launcher(workdir)


#: the header warning belongs in front of a real campaign, but once -- a
#: 42-tutorial sweep writes 42 of these.
_HEADER_WARNING_SHOWN = False


def _warn_about_headers() -> None:
    global _HEADER_WARNING_SHOWN
    if _HEADER_WARNING_SHOWN:
        return
    _HEADER_WARNING_SHOWN = True
    LOG.warning(
        "batch.header is generated per work directory from what SLURM and "
        "Lmod report here.  Read one before a real run and make sure every "
        "module the build needs is loaded, dependencies included -- a chain "
        "that is one module short fails inside the job, not at submission.")


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
#  verification
# --------------------------------------------------------------------------- #
def missing_artifacts(workdir: Path, patterns: Sequence[str]) -> list[str]:
    """Which of *patterns* match nothing under *workdir*."""
    return [pattern for pattern in patterns if not list(workdir.glob(pattern))]
