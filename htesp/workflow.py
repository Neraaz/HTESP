#!/usr/bin/env python
"""
HTESP workflow layer -- a single-file Python replacement for ``src/bash/*``.

Every bash scan script of the package is a method of :class:`HTESPWorkflow`.
The method keeps the *name* and the *argument contract* of the script it
replaces (``start``, ``end`` exclusive, track file, optional extra argument),
so ``mainprogram.py`` can keep dispatching to it unchanged.

Design notes
------------

*   **Every top-level "for material" loop is a multiprocessing loop.**
    :meth:`HTESPWorkflow.map` fans the per-material body out over
    ``min(len(materials), workers)`` processes and falls back to a plain serial
    loop when only one core (or one material) is available, or when
    ``HTESP_WORKERS=1`` / ``--workers 1`` is given.  Results come back ordered
    by the material index, so every list/CSV written after the loop
    (``result.csv``, ``econv.csv``, ``mpid-finished.in`` ...) is deterministic
    and no longer depends on scheduling order -- something the bash layer's
    in-loop counters could never guarantee.

*   **No shared mutable state between iterations.**  The bash layer wrote
    fixed-name scratch files (``mass.dat``, ``qpoint.dat``, ``kpoint.dat``,
    ``scf_dir/kpathlines.dat``, ``BZ.pdf``, ``temp*.in`` ...) into the project
    root, which is what made a naive ``&``/GNU-parallel port unsafe.  Here each
    material body runs inside its own private scratch directory
    (:meth:`HTESPWorkflow.scratch`) that carries a *local* ``scf_dir``; the
    helper modules keep their existing "write next to me" behaviour and the
    produced artefacts are moved to their final destination afterwards.

*   **No ``os.chdir`` leaking across iterations.**  :func:`pushd` always
    restores the previous directory, and a failed ``cd`` raises instead of
    letting the body (``rm -r ...``, ``sbatch ...``) run in the project root --
    the single most destructive behaviour of the bash layer.

*   **Repeated bash fragments became one method each**, e.g.
    :meth:`qe_section`, :meth:`final_coordinates`, :meth:`prefix_of`,
    :meth:`nbnd_heuristic`, :meth:`read_mesh_file`, :meth:`stage_and_submit`,
    :meth:`Scheduler.submit`.  The ~60% of the bash layer that was
    copy-pasted boilerplate exists exactly once here.

*   **Python helpers are imported, not exec'd.**  ``qe_input.py``,
    ``kpath.py``, ``elph.py``, ``band.py``, ``dos.py``, ``q2r.py``,
    ``matdyn.py``, ``matdyn_dos.py``, ``vasp_process.py``, ``plot.py`` &c. are
    imported once per process and their ``main()`` (or the function their
    ``if __name__ == "__main__"`` block calls) is invoked directly through
    :func:`run_helper`, which temporarily installs the ``sys.argv`` those
    modules expect.  No ``os.system`` round-trip, no interpreter start-up per
    material.

*   **Job ids are captured.**  ``sbatch --parsable`` ids are recorded in
    ``<stage dir>/.htesp_job.json`` so status checks no longer have to grep
    ``squeue`` output for a compound name (``B``, ``C``, ``Si`` ... used to
    match half the queue).

Usage
-----

    from htesp.workflow import HTESPWorkflow
    wf = HTESPWorkflow()                       # reads ./input.in
    wf.relax_scan(1, 5, "mpid.in")
    wf.create_inputs(1, 5, "mpid.in", nkpt=50)

or from the command line::

    python workflow.py relax-scan 1 5 mpid.in
    python workflow.py create-inputs 1 5 mpid.in 50 --workers 8
    python workflow.py phonopy-scan 1 5 mpid.in 1

Written as the Python port of the HTESP bash layer
(original bash by Niraj K. Nepal, Ph.D.).
"""

from __future__ import annotations

import contextlib
import csv
import glob
import importlib
import io
import json
import logging
import multiprocessing as mp
import os
import re
import functools
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

#: directory this module lives in (used to locate packaged data files)
_HERE = Path(__file__).resolve().parent

LOG = logging.getLogger("htesp")

RY_TO_EV = 13.605698
BOHR_TO_ANG = 0.52917725
AU_TO_ANG3 = 0.14818453          # bohr^3 -> angstrom^3 (as used by phonopy-scan)
PV_UNIT_CONV = 0.00062415        # kbar * A^3 -> eV
EV_TO_THZ = 241.79905043

#: directories the workflow creates below the project root
PROJECT_DIRS = ("scf_dir", "elph_dir", "q2r_dir", "matdyn_dir", "kpath", "plots")

#: name of the private scratch tree (one sub-directory per material per stage)
SCRATCH_ROOT = ".htesp_scratch"

#: optional marker line in ``batch.header``.  When it is absent the command
#: is appended instead, which is what ``mainprogram jobscript`` has always
#: done and what every shipped ``batch.header`` expects.
PLACEHOLDER = "submission here"


# --------------------------------------------------------------------------- #
#  small generic helpers
# --------------------------------------------------------------------------- #
@contextlib.contextmanager
def pushd(path: os.PathLike | str):
    """``cd path`` for the duration of the block, always restoring the old cwd.

    Replaces the 84 unguarded ``cd X ... cd ../../`` pairs of the bash layer.
    A missing directory raises here instead of silently leaving the body to run
    in the project root.
    """
    previous = Path.cwd()
    os.chdir(os.fspath(path))
    try:
        yield Path.cwd()
    finally:
        os.chdir(previous)


@contextlib.contextmanager
def argv(*args: Any):
    """Temporarily install ``sys.argv`` for a helper module that parses it."""
    saved = sys.argv
    sys.argv = [str(a) for a in args]
    try:
        yield
    finally:
        sys.argv = saved


def run_helper(module: str, entry: str, *args: Any, capture: bool = False):
    """Import ``module`` and call ``entry`` with ``sys.argv = [module, *args]``.

    This is the replacement for every ``os.system("<helper>.py a b c")`` call in
    the bash layer.  ``capture=True`` returns whatever the helper printed on
    stdout (needed for ``qe_axsf2cellpos.py``, which writes its result there).
    """
    mod = importlib.import_module(module if "." in module else "htesp." + module)
    func = getattr(mod, entry)
    buffer = io.StringIO()
    with argv(f"{module}.py", *args):
        if capture:
            with contextlib.redirect_stdout(buffer):
                func()
            return buffer.getvalue()
        func()
    return None


def read_text(path: os.PathLike | str, default: str = "") -> str:
    try:
        return Path(path).read_text(errors="replace")
    except OSError:
        return default


def read_lines(path: os.PathLike | str) -> list[str]:
    return read_text(path).splitlines()


def write_lines(path: os.PathLike | str, lines: Iterable[str]) -> None:
    Path(path).write_text("\n".join(lines) + ("\n" if lines else ""))


def grep(pattern: str, path: os.PathLike | str, regex: bool = False) -> list[str]:
    """``grep pattern path`` -> matching lines (empty list if the file is gone)."""
    lines = read_lines(path)
    if regex:
        rx = re.compile(pattern)
        return [ln for ln in lines if rx.search(ln)]
    return [ln for ln in lines if pattern in ln]


def grep_count(pattern: str, path: os.PathLike | str, regex: bool = False) -> int:
    """``grep pattern path | wc -l``."""
    return len(grep(pattern, path, regex=regex))


def field_of(line: str, index: int, default: str = "") -> str:
    """awk ``{print $index}`` (1-based), tolerant of short lines."""
    parts = line.split()
    return parts[index - 1] if 0 < index <= len(parts) else default


def last_match_field(pattern: str, path, index: int, default: str = "") -> str:
    """``grep pattern file | tail -n 1 | awk '{print $index}'``."""
    hits = grep(pattern, path)
    return field_of(hits[-1], index, default) if hits else default


def to_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(str(value).strip().rstrip(",").strip("'\"")))
    except (TypeError, ValueError):
        return default


def to_float(value: Any, default: float = float("nan")) -> float:
    try:
        return float(str(value).strip().rstrip(",").strip("'\""))
    except (TypeError, ValueError):
        return default


def is_number(value: Any) -> bool:
    return bool(re.fullmatch(r"[+-]?[0-9]+\.?[0-9]*([eEdD][+-]?[0-9]+)?", str(value).strip()))


def ensure_dir(path: os.PathLike | str) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def copy_file(src, dst, missing_ok: bool = True) -> bool:
    """``cp src dst`` -- returns False instead of raising when src is absent."""
    src, dst = Path(src), Path(dst)
    if not src.exists():
        if missing_ok:
            return False
        raise FileNotFoundError(src)
    if dst.is_dir():
        dst = dst / src.name
    ensure_dir(dst.parent)
    if src.is_dir():
        shutil.copytree(src, dst, dirs_exist_ok=True)
    else:
        shutil.copy2(src, dst)
    return True


def move_file(src, dst, missing_ok: bool = True) -> bool:
    src, dst = Path(src), Path(dst)
    if not src.exists():
        if missing_ok:
            return False
        raise FileNotFoundError(src)
    if dst.is_dir():
        dst = dst / src.name
    ensure_dir(dst.parent)
    if dst.exists():
        if dst.is_dir():
            shutil.rmtree(dst)
        else:
            dst.unlink()
    shutil.move(os.fspath(src), os.fspath(dst))
    return True


def remove(path, recursive: bool = False) -> None:
    """``rm`` / ``rm -r`` that never escapes to a wildcard in the wrong cwd."""
    p = Path(path)
    if p.is_symlink() or p.is_file():
        p.unlink()
    elif p.is_dir() and recursive:
        shutil.rmtree(p, ignore_errors=True)


def remove_glob(pattern: str, recursive: bool = False) -> None:
    for hit in glob.glob(pattern):
        remove(hit, recursive=recursive)


def concat(destination, *sources) -> None:
    """``cat a b c > destination`` (missing pieces are skipped, not fatal)."""
    chunks = []
    for src in sources:
        if src is None:
            continue
        if isinstance(src, (list, tuple)):
            chunks.append("\n".join(str(s) for s in src) + "\n")
        elif Path(src).is_file():
            chunks.append(read_text(src))
    text = "".join(chunks)
    if text and not text.endswith("\n"):
        text += "\n"
    ensure_dir(Path(destination).parent)
    Path(destination).write_text(text)


# --------------------------------------------------------------------------- #
#  material / track file
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Material:
    """One ``v<N> <mpid> <compound>`` entry of a track file."""

    index: int
    mpid: str
    compound: str
    root: Path

    @property
    def name(self) -> str:
        return f"R{self.mpid}-{self.compound}"

    @property
    def dir(self) -> Path:
        return self.root / self.name

    def sub(self, stage: str) -> Path:
        return self.dir / stage

    @property
    def scf_template(self) -> Path:
        return self.root / "scf_dir" / f"scf-{self.mpid}.in"

    @property
    def scf_relaxed(self) -> Path:
        return self.root / "scf_dir" / f"scf-relax-{self.mpid}-{self.compound}.in"

    def phonon_folder(self, prefer: str = "phonon") -> str:
        """``phonon/`` when it exists, otherwise ``calc/`` (bash convention)."""
        other = "calc" if prefer == "phonon" else "phonon"
        return prefer if self.sub(prefer).is_dir() else other


@dataclass
class Result:
    """What one material's body reports back to the (serial) finish phase."""

    index: int
    mpid: str = ""
    compound: str = ""
    status: str = "ok"
    message: str = ""
    job: str | None = None
    rows: dict = field(default_factory=dict)
    log: list[str] = field(default_factory=list)

    def say(self, *parts: Any) -> None:
        """Buffer a line; the parent process replays it in material order."""
        line = " ".join(str(p) for p in parts)
        self.log.append(line)
        LOG.debug(line)


def read_track_file(path, start: int, end: int, root: Path) -> list[Material]:
    """Parse ``v<N> <mpid> <compound>`` lines for ``start <= N < end``.

    ``end`` is *exclusive*, matching ``for ((ii=$1; ii<$2; ++ii))``.  Entries
    that the track file does not contain are reported once and skipped instead
    of producing an ``R-`` path (the bash layer happily built ``R-/relax`` and
    then ran ``rm -r`` in the project root).
    """
    entries: dict[int, tuple[str, str]] = {}
    for raw in read_lines(path):
        parts = raw.split()
        if len(parts) < 3 or not parts[0].startswith("v"):
            continue
        try:
            idx = int(parts[0][1:])
        except ValueError:
            continue
        entries[idx] = (parts[1], parts[2])

    materials: list[Material] = []
    for idx in range(int(start), int(end)):
        if idx not in entries:
            LOG.warning("v%s not found in %s -- skipped", idx, path)
            continue
        mpid, compound = entries[idx]
        materials.append(Material(idx, mpid, compound, root))
    return materials


def write_track_file(path, rows: Sequence[tuple[str, str]]) -> None:
    """Write ``v1 ... vN`` with deterministic numbering, once, after the loop."""
    if not rows:
        remove(path)
        return
    write_lines(path, [f"v{i} {mpid} {comp}" for i, (mpid, comp) in enumerate(rows, 1)])


# --------------------------------------------------------------------------- #
#  scheduler
# --------------------------------------------------------------------------- #
class Scheduler:
    """Job submission, job-id capture and queue queries.

    Replaces ``jobscript.sh`` (``qerun`` / ``vasprun``) and the five places
    where that ladder was re-inlined at a different relative depth.  The
    ``CALC_VISIBLE_WITH_*`` marker files are read **once** from the project root
    rather than looked up at ``../../``, ``../../../`` and ``../../../../``.
    """

    def __init__(self, root: Path, command: str = "sbatch", dry_run: bool = False,
                 throttle: float = 0.25):
        self.root = Path(root)
        self.command = command
        self.dry_run = dry_run
        self.throttle = throttle
        self.name_style = self._detect_name_style()

    def _detect_name_style(self) -> str:
        if (self.root / "CALC_VISIBLE_WITH_ID").is_file():
            return "id"
        if (self.root / "CALC_VISIBLE_WITH_NAME").is_file():
            return "name"
        if (self.root / "CALC_VISIBLE_WITH_ID-NAME").is_file():
            return "id-name"
        return "plain"

    def job_name(self, mpid: str, compound: str, tag: str | None = None) -> str:
        stem = {"id": mpid, "name": compound, "id-name": f"{mpid}-{compound}"}.get(
            self.name_style, "run" if tag is None else f"run-{tag}")
        if self.name_style == "plain":
            return f"{stem}.sh" if tag is None else f"run-{tag}.sh"
        return f"{stem}.sh" if tag is None else f"{stem}-{tag}.sh"

    def submit(self, script: str, cwd, mpid: str = "", compound: str = "",
               tag: str | None = None, rename: bool = True) -> str | None:
        """Submit ``script`` from ``cwd``; returns the job id (or None).

        ``rename`` reproduces the ``mv run-<tag>.sh <name>.sh`` behaviour that
        makes the job visible under a useful name in ``squeue``.
        """
        cwd = Path(cwd)
        source = cwd / script
        if not source.is_file():
            LOG.warning("submission script %s not found in %s", script, cwd)
            return None

        target = source
        if rename:
            wanted = cwd / self.job_name(mpid, compound, tag)
            if wanted != source:
                shutil.move(os.fspath(source), os.fspath(wanted))
            target = wanted

        if self.dry_run:
            LOG.info("[dry-run] %s %s (in %s)", self.command, target.name, cwd)
            return None

        try:
            proc = subprocess.run(
                [self.command, "--parsable", target.name],
                cwd=os.fspath(cwd), capture_output=True, text=True, check=False)
        except FileNotFoundError:
            LOG.error("%s not found -- is the scheduler available on this host?",
                      self.command)
            return None

        if proc.returncode != 0:
            # --parsable is Slurm-specific; retry plain, then give up loudly.
            proc = subprocess.run([self.command, target.name], cwd=os.fspath(cwd),
                                  capture_output=True, text=True, check=False)
        if proc.returncode != 0:
            LOG.error("submission failed in %s: %s", cwd, proc.stderr.strip())
            return None

        job_id = ""
        for token in proc.stdout.split():
            if token.strip().split(";")[0].isdigit():
                job_id = token.strip().split(";")[0]
                break
        self.record(cwd, job_id, tag or "run")
        if self.throttle:
            time.sleep(self.throttle)
        return job_id or None

    @staticmethod
    def record(cwd, job_id: str, tag: str) -> None:
        if not job_id:
            return
        store = Path(cwd) / ".htesp_job.json"
        try:
            data = json.loads(store.read_text())
        except (OSError, ValueError):
            data = {}
        data.setdefault(tag, []).append({"job": job_id, "time": time.time()})
        store.write_text(json.dumps(data, indent=1))

    @staticmethod
    def job_ids(cwd) -> list[str]:
        try:
            data = json.loads((Path(cwd) / ".htesp_job.json").read_text())
        except (OSError, ValueError):
            return []
        return [entry["job"] for runs in data.values() for entry in runs]

    def is_running(self, cwd) -> bool:
        """True when any job recorded for ``cwd`` is still queued or running.

        This replaces ``squeue -o '%.100j' | grep "$B"``, which matched any job
        whose name merely *contained* the compound name.
        """
        ids = self.job_ids(cwd)
        if not ids:
            return False
        try:
            proc = subprocess.run(["squeue", "-h", "-o", "%i", "-j", ",".join(ids)],
                                  capture_output=True, text=True, check=False)
        except FileNotFoundError:
            return False
        return bool(proc.stdout.strip())


# --------------------------------------------------------------------------- #
#  input.in / qpoint.in / kpoint.in / pressure.in / charge.in
# --------------------------------------------------------------------------- #
from htesp.inputin import InputIn  # one implementation, shared with the CLI


def read_mesh_file(path, base: Sequence[int], default_divisor: int = 2
                   ) -> tuple[list[int], list[int]]:
    """Parse ``qpoint.in`` / ``kpoint.in``; returns ``(mesh, shift)``.

    One implementation for the five copies the bash layer carried.  Accepted
    contents: ``f`` (divide ``base`` by f), ``n1 n2 n3``, ``f s1 s2 s3`` and
    ``n1 n2 n3 s1 s2 s3``.  A missing file means ``base / default_divisor``.
    """
    base = [to_int(b, 1) for b in base]
    shift = [0, 0, 0]
    tokens = read_text(path).split()
    if not tokens:
        return [max(1, b // default_divisor) for b in base], shift

    if len(tokens) == 1:
        frac = to_float(tokens[0], float(default_divisor)) or default_divisor
        return [max(1, int(b / frac)) for b in base], shift
    if len(tokens) == 3:
        return [to_int(t, 1) for t in tokens], shift
    if len(tokens) == 4:
        frac = to_float(tokens[0], float(default_divisor)) or default_divisor
        return ([max(1, int(b / frac)) for b in base],
                [to_int(t, 0) for t in tokens[1:4]])
    if len(tokens) >= 6:
        return [to_int(t, 1) for t in tokens[:3]], [to_int(t, 0) for t in tokens[3:6]]

    LOG.warning("%s: give either a divisor or an explicit mesh", path)
    return [max(1, b // default_divisor) for b in base], shift


def read_indexed_file(path) -> list[tuple[str, str]]:
    """``v<N> <value>`` files: ``pressure.in``, ``charge.in`` -> ``[(key, value)]``."""
    rows = []
    for line in read_lines(path):
        parts = line.split()
        if len(parts) >= 2 and parts[0].startswith("v"):
            rows.append((parts[0], parts[1]))
    return rows


# --------------------------------------------------------------------------- #
#  Quantum ESPRESSO text handling (the only real "logic" of the bash layer)
# --------------------------------------------------------------------------- #
class QEText:
    """Section extraction / editing for ``pw.x`` inputs and outputs.

    The bash layer re-derived section boundaries with
    ``sed -n '/A/,/B/p' | sed '$d'`` 26 times and hard-coded the
    ``Begin final coordinates`` block layout in 8 places.  All of that lives
    here once, and the relax-vs-vc-relax difference (no ``CELL_PARAMETERS``
    block for ``calculation='relax'``) is handled instead of silently deleting
    the first atoms.
    """

    CARDS = ("ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS",
             "CELL_PARAMETERS", "OCCUPATIONS", "CONSTRAINTS", "ATOMIC_FORCES")

    # ---- reading ---------------------------------------------------------- #
    @staticmethod
    def section(text: str, start_pat: str, end_pat: str | None = None,
                include_end: bool = False) -> list[str]:
        """Lines from the first ``start_pat`` up to (not including) ``end_pat``."""
        lines = text.splitlines()
        out: list[str] = []
        started = False
        for line in lines:
            if not started:
                if start_pat in line:
                    started = True
                    out.append(line)
                continue
            if end_pat is not None and end_pat in line:
                if include_end:
                    out.append(line)
                return out
            out.append(line)
        return out

    @classmethod
    def card(cls, text: str, name: str) -> list[str]:
        """The named card, up to the next card keyword."""
        lines = text.splitlines()
        out: list[str] = []
        started = False
        for line in lines:
            head = line.strip().split()[0] if line.strip() else ""
            if not started:
                if head == name:
                    started = True
                    out.append(line)
                continue
            if head in cls.CARDS and head != name:
                break
            out.append(line)
        while out and not out[-1].strip():
            out.pop()
        return out

    @classmethod
    def header(cls, text: str) -> list[str]:
        """Everything before the first card keyword (the namelists)."""
        out = []
        for line in text.splitlines():
            head = line.strip().split()[0] if line.strip() else ""
            if head in cls.CARDS:
                break
            out.append(line)
        return out

    @staticmethod
    def value(text: str, key: str, default: str = "") -> str:
        """``key = value,`` from a namelist (first occurrence)."""
        match = re.search(rf"^\s*{re.escape(key)}\s*=\s*([^,\n]+)", text, re.M)
        return match.group(1).strip().rstrip(",").strip() if match else default

    @classmethod
    def prefix(cls, text: str) -> str:
        """``prefix = 'x',`` -> ``'x'`` (quotes kept: helpers expect them)."""
        raw = cls.value(text, "prefix", "''")
        return raw if raw.startswith("'") else f"'{raw}'"

    @staticmethod
    def kmesh(text: str) -> tuple[list[int], list[int]]:
        """``K_POINTS automatic`` mesh and shift."""
        lines = text.splitlines()
        for i, line in enumerate(lines):
            if line.strip().startswith("K_POINTS") and i + 1 < len(lines):
                nums = [to_int(t, 0) for t in lines[i + 1].split()]
                nums += [0] * (6 - len(nums))
                return nums[:3], nums[3:6]
        return [1, 1, 1], [0, 0, 0]

    # ---- editing ---------------------------------------------------------- #
    @staticmethod
    def drop(lines: Sequence[str], *patterns: str) -> list[str]:
        """``sed '/pat/d'`` for each pattern."""
        return [ln for ln in lines if not any(p in ln for p in patterns)]

    @staticmethod
    def drop_namelist(lines: Sequence[str], name: str) -> list[str]:
        """Remove ``&NAME ... /`` (used to strip ``&IONS``/``&CELL`` for scf)."""
        out, skipping = [], False
        for line in lines:
            stripped = line.strip()
            if not skipping and stripped.upper().startswith(f"&{name.upper()}"):
                skipping = True
                continue
            if skipping:
                if stripped in ("/", "&END") or stripped.startswith("/"):
                    skipping = False
                continue
            out.append(line)
        return out

    @staticmethod
    def insert_after(lines: Sequence[str], anchor: str, *new: str) -> list[str]:
        """``sed '/anchor/a text'`` (after the first match only)."""
        out, done = [], False
        for line in lines:
            out.append(line)
            if not done and anchor in line:
                out.extend(new)
                done = True
        if not done:
            LOG.debug("anchor %r not found -- %r not inserted", anchor, new)
        return out

    @staticmethod
    def replace(lines: Sequence[str], old: str, new: str) -> list[str]:
        return [ln.replace(old, new) for ln in lines]

    @classmethod
    def set_key(cls, lines: Sequence[str], namelist: str, key: str, value: Any
                ) -> list[str]:
        """Set ``key = value`` inside ``&NAMELIST``, replacing any old entry."""
        cleaned = [ln for ln in lines
                   if not re.match(rf"^\s*{re.escape(key)}\s*=", ln)]
        return cls.insert_after(cleaned, f"&{namelist.upper()}",
                                f"  {key} = {value},")

    # ---- pw.x output ------------------------------------------------------ #
    @staticmethod
    def final_coordinates(scf_out) -> list[str]:
        """``Begin/End final coordinates`` -> ``CELL_PARAMETERS`` + positions.

        Bash did ``sed '$d' | sed '1,4d' | sed '5d'``, which assumes the
        vc-relax layout (volume/density lines, then the cell block).  For
        ``calculation = 'relax'`` that deleted the header *and the first two
        atoms*.  Here the block is parsed by keyword, so both work and an
        unconverged run yields an empty list instead of a corrupt file.
        """
        text = read_text(scf_out)
        if "Begin final coordinates" not in text:
            return []
        block = QEText.section(text, "Begin final coordinates",
                               "End final coordinates")[1:]
        out: list[str] = []
        keep = False
        for line in block:
            head = line.strip().split()[0] if line.strip() else ""
            if head in ("CELL_PARAMETERS", "ATOMIC_POSITIONS"):
                keep = True
            if keep:
                out.append(line.replace("(angstrom)", "angstrom")
                               .replace("(crystal)", "crystal")
                               .replace("(alat=", "alat=") if head in
                           ("CELL_PARAMETERS", "ATOMIC_POSITIONS") else line)
        while out and not out[-1].strip():
            out.pop()
        return out

    @staticmethod
    def last_cell_block(scf_out, natoms: int) -> list[str]:
        """Last ``CELL_PARAMETERS (angstrom)`` + positions of an unfinished run."""
        lines = read_lines(scf_out)
        starts = [i for i, ln in enumerate(lines) if "CELL_PARAMETERS" in ln]
        if not starts:
            return []
        begin = starts[-1]
        block = lines[begin:begin + natoms + 6]
        return [ln.replace("(angstrom)", "angstrom").replace("(crystal)", "crystal")
                for ln in block if ln.strip()]

    @staticmethod
    def nelec(scf_out) -> int:
        return to_int(last_match_field("number of electrons", scf_out, 5), 0)

    @staticmethod
    def natoms(scf_out) -> int:
        return to_int(last_match_field("number of atoms/cell", scf_out, 5), 0)

    @staticmethod
    def total_energy(scf_out) -> float:
        hits = [ln for ln in read_lines(scf_out) if ln.strip().startswith("!")]
        return to_float(field_of(hits[-1], 5), float("nan")) if hits else float("nan")

    @staticmethod
    def iterations(scf_out) -> int:
        return len([ln for ln in read_lines(scf_out) if ln.strip().startswith("!")])

    @staticmethod
    def volume(scf_out) -> float:
        hits = grep("volume", scf_out)
        return to_float(field_of(hits[-1], 4), float("nan")) if hits else float("nan")

    @staticmethod
    def pressure(scf_out) -> float:
        hits = grep("(kbar)     P=", scf_out)
        return to_float(field_of(hits[-1], 6), float("nan")) if hits else float("nan")

    @staticmethod
    def is_relaxed(scf_out) -> bool:
        return grep_count("Final scf calculation at the relaxed structure.", scf_out) > 0

    @staticmethod
    def max_phonon_freq(elph_out) -> float:
        hits = grep("THz", elph_out)
        return to_float(field_of(hits[-1], 5), 0.0) if hits else 0.0


# --------------------------------------------------------------------------- #
#  the workflow object
# --------------------------------------------------------------------------- #
class HTESPWorkflow:
    """All former ``src/bash`` scripts, as methods, parallel by default.

    Parameters
    ----------
    root
        Project directory (the one holding ``input.in``); defaults to the cwd.
    workers
        Size of the process pool used for every per-material loop.  ``None``
        means ``min(os.cpu_count(), 8)``; ``1`` disables multiprocessing.  The
        environment variable ``HTESP_WORKERS`` overrides the default.
    dry_run
        Build every input file but never call the scheduler.
    """

    # ------------------------------------------------------------------ setup
    def __init__(self, root: os.PathLike | str = ".", workers: int | None = None,
                 dry_run: bool = False, submit_command: str = "sbatch",
                 keep_scratch: bool = False, log_level: int = logging.INFO):
        self.root = Path(root).resolve()
        self.dry_run = bool(dry_run)
        self.keep_scratch = bool(keep_scratch)
        self.submit_command = submit_command
        self.workers = self._resolve_workers(workers)

        if not LOG.handlers:
            logging.basicConfig(level=log_level, format="%(message)s")
        LOG.setLevel(log_level)

        #: results of the most recent :meth:`map` call (see :attr:`failed_count`)
        self.failures: list[Result] = []
        self.skipped: list[Result] = []
        self.last_total = 0

        self.config = self._load_config()
        self.input = InputIn.load(self.root / "input.in", self.config)
        self.scheduler = Scheduler(self.root, command=submit_command,
                                   dry_run=self.dry_run)

    @staticmethod
    def _resolve_workers(workers: int | None) -> int:
        if workers is None:
            workers = to_int(os.environ.get("HTESP_WORKERS", 0), 0)
        if not workers:
            workers = min(os.cpu_count() or 1, 8)
        return max(1, int(workers))

    def _load_config(self) -> dict:
        """The effective configuration for this project (see :mod:`htesp.config`)."""
        from htesp.config import config as _config
        return _config(self.root)

    # ---------------------------------------------------------------- picking
    def __getstate__(self) -> dict:
        """Everything on the instance is plain data, so workers can unpickle it."""
        state = self.__dict__.copy()
        state["scheduler"] = None
        return state

    def __setstate__(self, state: dict) -> None:
        self.__dict__.update(state)
        self.scheduler = Scheduler(self.root, command=self.submit_command,
                                   dry_run=self.dry_run)

    # ------------------------------------------------------------- bookkeeping
    @property
    def dft(self) -> str:
        return self.input.dft

    @property
    def is_vasp(self) -> bool:
        return self.input.is_vasp

    def materials(self, start: int, end: int, track: str | None = None
                  ) -> list[Material]:
        """Resolve ``[start, end)`` of ``track`` (default: ``input.in`` line 4).

        The bash layer located the track file with
        ``find * -name $3 | tail -n 1``: a full-tree walk on every invocation
        that could pick a nested copy under ``completed/`` or
        ``R*/pressure/``.  Here the name is resolved against the project root
        and only searched for if that fails.
        """
        name = track or self.input.track
        path = self.root / name
        if not path.is_file():
            hits = sorted(self.root.glob(f"*/{name}")) + sorted(self.root.glob(f"*/*/{name}"))
            if not hits:
                raise FileNotFoundError(f"track file {name!r} not found under {self.root}")
            path = hits[0]
            LOG.warning("track file %s taken from %s", name, path)
        return read_track_file(path, start, end, self.root)

    def banner(self, *message: str) -> None:
        rule = "-" * 111
        LOG.info(rule)
        for line in message:
            LOG.info(line)
        LOG.info(rule)

    # ------------------------------------------------------------ parallelism
    def map(self, body: Callable[[Material], Result], materials: Sequence[Material],
            parallel: bool = True) -> list[Result]:
        """Run ``body`` for every material, in parallel, ordered results back.

        This is *the* top-level loop of the package.  ``body`` must be a bound
        method of this instance (so it pickles), must not rely on the process
        cwd persisting between calls, and must not write fixed-name files into
        the project root -- use :meth:`scratch` for that.
        """
        if not materials:
            return []

        workers = min(self.workers, len(materials)) if parallel else 1
        if workers <= 1:
            results = [self._guard(body, m) for m in materials]
        else:
            ctx = mp.get_context("fork" if hasattr(os, "fork") else "spawn")
            LOG.info("running %d materials on %d workers", len(materials), workers)
            with ctx.Pool(processes=workers) as pool:
                results = pool.map(_pool_entry, [(self, body.__func__.__name__, m)
                                                 for m in materials])

        results.sort(key=lambda r: r.index)
        for res in results:
            for line in res.log:
                LOG.info(line)
            if res.status == "failed":
                LOG.error("%s %s: %s", res.mpid, res.compound, res.message)

        # Failure accounting.  The bash layer discarded every exit code, so a
        # stage in which every single material blew up still "succeeded" and
        # the next stage ran on nothing.  The counts are kept on the instance
        # so the CLI can turn them into a non-zero exit status.
        self.failures = [res for res in results if res.status == "failed"]
        self.skipped = [res for res in results if res.status == "skipped"]
        self.last_total = len(results)
        if self.failures:
            LOG.error("%d of %d materials failed", len(self.failures), len(results))
        if self.skipped:
            LOG.info("%d of %d materials skipped", len(self.skipped), len(results))
        return results

    @property
    def failed_count(self) -> int:
        """How many materials failed in the most recent :meth:`map` call."""
        return len(getattr(self, "failures", []))

    def failure_summary(self) -> str:
        """One line per failed material, for the caller's error report."""
        return "\n".join(f"  {res.mpid} {res.compound}: {res.message}"
                          for res in getattr(self, "failures", []))

    def _guard(self, body: Callable[[Material], Result], material: Material) -> Result:
        """Run one material's body; a crash fails that material only."""
        result = Result(material.index, material.mpid, material.compound)
        try:
            produced = body(material)
            if isinstance(produced, Result):
                produced.index = material.index
                produced.mpid = produced.mpid or material.mpid
                produced.compound = produced.compound or material.compound
                return produced
        except Exception as exc:                      # noqa: BLE001 - per-material isolation
            result.status = "failed"
            result.message = f"{type(exc).__name__}: {exc}"
            LOG.exception("material %s-%s failed", material.mpid, material.compound)
        return result

    # --------------------------------------------------------------- scratch
    @contextlib.contextmanager
    def scratch(self, material: Material, stage: str, copy: Sequence[str] = (),
                local_scf_dir: bool = True):
        """A private working directory for one material.

        Everything the helper modules write with a fixed name (``mass.dat``,
        ``qpoint.dat``, ``kpoint.dat``, ``BZ.pdf``, ``scf_dir/kpathlines.dat``,
        ``elph-<id>-<comp>.in`` ...) lands here instead of in the shared project
        root, which is what makes the loop safe to parallelise.

        ``copy`` names files, relative to the project root, that must be visible
        inside the scratch directory (``config.json`` and ``qpoint.in`` are
        always copied when present).
        """
        work = self.root / SCRATCH_ROOT / f"{stage}-{material.mpid}-{material.compound}"
        if work.exists():
            shutil.rmtree(work, ignore_errors=True)
        ensure_dir(work)
        if local_scf_dir:
            ensure_dir(work / "scf_dir")
        for name in ("config.json", "qpoint.in", "kpoint.in", *copy):
            copy_file(self.root / name, work / name)
        try:
            with pushd(work):
                yield work
        finally:
            if not self.keep_scratch:
                shutil.rmtree(work, ignore_errors=True)

    def collect_scratch(self, work: Path, patterns: dict[str, str]) -> list[str]:
        """Move artefacts out of a scratch dir: ``{glob: destination dir}``."""
        moved = []
        for pattern, destination in patterns.items():
            target = ensure_dir(self.root / destination)
            for hit in sorted(work.glob(pattern)):
                move_file(hit, target / hit.name)
                moved.append(hit.name)
        return moved

    # ------------------------------------------------------ destructive steps
    def remove(self, path, recursive: bool = False) -> None:
        """``remove()``, but a no-op under ``--dry-run``.

        ``mainprogram 20`` (clean-scan) and ``mainprogram 28`` (pressure-reset)
        are the two stages that delete real results.  They used to delete
        unconditionally, so ``--dry-run`` destroyed exactly what the user was
        checking.  Every deletion in those stages goes through here.
        """
        if self.dry_run:
            LOG.info("[dry-run] would remove %s%s", path, " (recursive)" if recursive else "")
            return
        remove(path, recursive=recursive)

    def remove_glob(self, pattern: str, recursive: bool = False) -> None:
        """:meth:`remove` for a glob; also a no-op under ``--dry-run``."""
        hits = glob.glob(pattern)
        if self.dry_run:
            for hit in hits:
                LOG.info("[dry-run] would remove %s", hit)
            return
        for hit in hits:
            remove(hit, recursive=recursive)

    # -------------------------------------------------------- shared fragments
    def ensure_project_dirs(self, *names: str) -> None:
        for name in (names or PROJECT_DIRS):
            ensure_dir(self.root / name)

    def prefix_of(self, material: Material) -> str:
        """``prefix = 'x',`` from the relax input, falling back to the template.

        16 copies of ``grep "prefix = " ... | sed 's/.$//'`` in bash.
        """
        for candidate in (material.sub("relax") / "scf.in", material.scf_template,
                          material.scf_relaxed):
            if Path(candidate).is_file():
                prefix = QEText.prefix(read_text(candidate))
                if prefix not in ("''", ""):
                    return prefix
        return f"'{material.compound}'"

    @staticmethod
    def nbnd_heuristic(nelec: int, soc: bool = False, spin: bool = False,
                       style: str = "qe") -> int:
        """One place for what bash spelled out three (mutually inconsistent) ways.

        ``style='qe'``    : +20 with SOC/spin, +10 below 30 electrons, else nelec
        ``style='vasp'``  : +20 with SOC, +10 below 30 electrons, else 1.5*n/2
        ``style='epw'``   : 2*n with SOC, else n+10
        """
        if style == "epw":
            return nelec * 2 if soc else nelec + 10
        if soc or spin:
            return nelec + 20
        if nelec < 30:
            return nelec + 10
        return nelec if style == "qe" else int(nelec * 1.50 / 2)

    def qmesh_for(self, kmesh: Sequence[int]) -> list[int]:
        """q-mesh from ``qpoint.in`` (or ``k-mesh / 2`` when it is absent)."""
        mesh, _ = read_mesh_file(self.root / "qpoint.in", kmesh, default_divisor=2)
        return mesh

    def stage_and_submit(self, material: Material, stage: str, script: str,
                         tag: str, files: dict[str, str] | None = None,
                         result: Result | None = None) -> Result:
        """``mkdir -p <stage>; cp <files>; cp <script>; cd; sbatch; cd ../../``

        The submit block that the bash layer repeated in 28 places.
        ``files`` maps *destination name inside the stage dir* to a path
        relative to the project root.
        """
        result = result or Result(material.index, material.mpid, material.compound)
        target = ensure_dir(material.sub(stage))
        if self.is_vasp:                       # see _submit_vasp for why
            from htesp.write_potcar import stage_potcar

            stage_potcar(target)
        for destination, source in (files or {}).items():
            if not copy_file(self.root / source, target / destination):
                result.say(f"  missing input {source} -- {destination} not staged")
        if not copy_file(self.root / script, target / script):
            result.status = "skipped"
            result.message = f"{script} not found in project root"
            return result
        result.say(f" submitting jobs in {material.name}/{stage}")
        result.job = self.scheduler.submit(script, target, material.mpid,
                                           material.compound, tag)
        return result



    # ------------------------------------------------------------------------ #
    #  helper-module shims (the former ``os.system("<helper>.py ...")`` calls)
    # ------------------------------------------------------------------------ #
    def run_qe_input(self, mpid: str) -> None:
        """``qe_input.py <mpid>`` -- download/convert an MP structure to scf.in."""
        importlib.import_module("htesp.qe_input").qe_input(mpid)

    def run_kpath(self, mode: str, *args: Any) -> None:
        """``kpath.py point|line ...`` (plus the ``KPT_OPT`` marker its main writes)."""
        if self.config.get("kpt_opt", False):
            Path("KPT_OPT").touch()
        run_helper("kpath", "main", mode, *args)

    def run_elph(self, material: Material, prefix: str) -> None:
        """``elph.py <mpid> <comp> <prefix>`` -- imported, mode from config."""
        elph = importlib.import_module("htesp.elph")
        mode = self.config.get("elph_mode", "serial")
        elph.elph_in(material.mpid, material.compound, prefix, elphmode=mode)
        elph.irr_q(material.mpid, material.compound)

    def run_vasp_process(self, what: str) -> None:
        """``vasp_process.py POSCAR|eigen|conventional|symmetrize|<scf.in>``."""
        run_helper("vasp_process", "main", what)

    def run_plot(self, *args: Any) -> None:
        run_helper("plot", "main", *args)

    def run_scftocif(self, *args: Any) -> None:
        run_helper("scftocif", "main", *args)

    def run_create_epw(self, *args: Any) -> None:
        run_helper("create_epw_inputs", "main", *args)

    def run_create_wt(self, *args: Any) -> None:
        run_helper("create_wt_inputs", "main", *args)

    # ======================================================================== #
    #  1  relax-scan
    # ======================================================================== #
    def relax_scan(self, start: int, end: int, track: str | None = None, *_) -> list[Result]:
        """``mainprogram 1`` -- submit the first structural relaxation."""
        self.banner("Submitting crystal structure relaxation")
        results = self.map(self._relax_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _relax_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)

        if self.is_vasp:
            if not material.dir.is_dir():
                res.status = "skipped"
                res.message = f"{material.name} not present"
                return res
            relax = ensure_dir(material.sub("relax"))
            copy_file(self.root / "vdw_kernel.bindat", relax / "vdw_kernel.bindat")
            return self._submit_vasp(material, "relax", res)

        # --- QE -------------------------------------------------------------
        if not material.scf_template.is_file():
            res.say(f" scf_dir/scf-{material.mpid}.in missing -- fetching from MP")
            with pushd(self.root):
                self.run_qe_input(material.mpid)
        if not material.scf_template.is_file():
            res.status = "skipped"
            res.message = f"scf-{material.mpid}.in could not be created"
            return res

        relax = ensure_dir(material.sub("relax"))
        lines = read_lines(material.scf_template)
        lines = QEText.drop(lines, "pseudo_dir")
        # bash anchored on the exact line "calculation = 'vc-relax'," and lost
        # pseudo_dir for every other calculation type; anchor on &CONTROL.
        lines = QEText.insert_after(lines, "&CONTROL", "pseudo_dir = '../../pp/',")
        write_lines(relax / "scf.in", lines)
        return self.stage_and_submit(material, "relax", "run-scf.sh", "scf", result=res)

    def _submit_vasp(self, material: Material, stage: str, res: Result,
                     script: str = "run-vasp.sh") -> Result:
        """``cp run-vasp.sh <stage>/run.sh; cd; vasprun; cd ../../``."""
        target = ensure_dir(material.sub(stage))
        # FIX: build the POTCAR here if it is missing.  Only the *download*
        # path (htesp/vasp_input.py) wrote one, so a stage directory that was
        # seeded rather than downloaded -- which is every tutorial that starts
        # from a prepared R<mpid>-<compound>/relax/, and any directory a user
        # assembled by hand -- went to the scheduler with INCAR, KPOINTS and
        # POSCAR but no POTCAR, and VASP stopped immediately.  stage_potcar
        # never raises: with no POTCARs configured it explains how and the
        # inputs are still written.
        from htesp.write_potcar import stage_potcar

        stage_potcar(target)
        if not copy_file(self.root / script, target / "run.sh"):
            res.status = "skipped"
            res.message = f"{script} not found in project root"
            return res
        res.say(f" submitting jobs in {material.name}/{stage}")
        res.job = self.scheduler.submit("run.sh", target, material.mpid,
                                        material.compound, None)
        return res

    # ======================================================================== #
    #  2  further-relax-input
    # ======================================================================== #
    def further_relax_input(self, start: int, end: int, track: str | None = None,
                            *_) -> list[Result]:
        """``mainprogram 2`` -- harvest the relaxed structure into a new input."""
        self.banner(f"Extracting relaxed structure, and updating {self.dft} input file")
        materials = self.materials(start, end, track)
        results = self.map(self._further_relax_input_one, materials)
        # the tracking list is written once, after the loop, with deterministic
        # numbering -- the bash counter made it submission-order dependent.
        write_track_file(self.root / "mpid-list-not-relaxed.in",
                         [(r.mpid, r.compound) for r in results
                          if r.rows.get("not_relaxed")])
        LOG.info("complete")
        return results

    def _further_relax_input_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        res.rows = {"mpid": material.mpid, "compound": material.compound}
        relax = material.sub("relax")

        if self.is_vasp:
            outcar, incar = relax / "OUTCAR", relax / "INCAR"
            if not outcar.is_file():
                res.status = "skipped"
                res.message = "OUTCAR not found"
                return res
            converged = grep_count(
                "reached required accuracy - stopping structural energy minimisation",
                outcar) > 0
            niter = grep_count("y  w", outcar)
            nsw_zero = grep_count("NSW = 0", incar) > 0

            if converged:
                nposcar = len(glob.glob(os.fspath(relax / "POSCAR*")))
                noutcar = len(glob.glob(os.fspath(relax / "OUTCAR*")))
                res.say(f"{nposcar} POSCARs found")
                copy_file(relax / "POSCAR", relax / f"POSCAR{nposcar}")
                copy_file(outcar, relax / f"OUTCAR{noutcar}")
                copy_file(relax / "CONTCAR", relax / "POSCAR")
                if (relax / "CHGCAR").is_file():
                    remove_glob(os.fspath(relax / "CHG*"))
                    remove(relax / "WAVECAR")
                if niter < 2:
                    write_lines(incar, [re.sub(r"^\s*NSW\s*=.*", "NSW = 0", ln)
                                        for ln in read_lines(incar)])
            else:
                copy_file(relax / "CONTCAR", relax / "POSCAR")
                if (relax / "CHGCAR").is_file():
                    remove_glob(os.fspath(relax / "CHG*"))
                    remove(relax / "WAVECAR")
                if nsw_zero:
                    (relax / "NSW_0_DETECTED").touch()
            return res

        # --- QE -------------------------------------------------------------
        scf_out = relax / "scf.out"
        if not material.dir.is_dir():
            res.status = "skipped"
            res.message = f"{material.name} not present"
            return res
        if not scf_out.is_file():
            res.status = "skipped"
            res.message = "relax/scf.out not found -- has the job run?"
            return res

        niter = QEText.iterations(scf_out)
        if niter and niter < 3:
            res.say("Structure already relaxed")
            return res

        if QEText.is_relaxed(scf_out):
            res.say(f"Relaxed structure found for {material.mpid} {material.compound}")
            nscf = len(glob.glob(os.fspath(relax / "scf.out*")))
            copy_file(scf_out, relax / f"scf.out{nscf}")
            copy_file(relax / "scf.in", relax / f"scf.in{nscf}")

            block = QEText.final_coordinates(scf_out)
            if not block:
                res.status = "failed"
                res.message = "final coordinates block could not be parsed"
                return res
            template = read_text(material.scf_template)
            header = QEText.header(template) + QEText.card(template, "ATOMIC_SPECIES")
            kpoints = QEText.card(template, "K_POINTS")
            positions = QEText.section("\n".join(block), "ATOMIC_POSITIONS")
            cell = QEText.section("\n".join(block), "CELL_PARAMETERS",
                                  "ATOMIC_POSITIONS")
            write_lines(material.scf_relaxed, header + positions + kpoints + cell)
            res.say(f"  wrote {material.scf_relaxed.name}")
        else:
            res.say(f"Not relaxed: {material.mpid} {material.compound} (time out)")
            scf_in = relax / "scf.in"
            header = QEText.header(read_text(scf_in)) + \
                QEText.card(read_text(scf_in), "ATOMIC_SPECIES")
            natoms = QEText.natoms(scf_out)
            cell_block = QEText.last_cell_block(scf_out, natoms)
            if not cell_block:
                res.status = "failed"
                res.message = "no CELL_PARAMETERS block in scf.out -- input left untouched"
                res.rows["not_relaxed"] = True
                return res
            copy_file(scf_in, relax / "scf-initial.in")
            copy_file(scf_out, relax / "scf-initial.out")
            kpoints = QEText.card(read_text(scf_in), "K_POINTS")
            write_lines(scf_in, header + cell_block + kpoints)
            res.rows["not_relaxed"] = True
        return res

    # ======================================================================== #
    #  3 / 5 / 6  resubmission of scf-type runs
    # ======================================================================== #
    def further_relax_scan(self, start: int, end: int, track: str | None = None,
                           *_) -> list[Result]:
        """``mainprogram 3`` -- resubmit the relaxation with the updated input."""
        self.banner("Resubmitting crystal relaxation calculations.")
        results = self.map(self._further_relax_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _further_relax_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        relax = material.sub("relax")

        if self.is_vasp:
            if (relax / "NSW_0_DETECTED").is_file():
                res.say("Already NSW = 0 found")
                res.status = "skipped"
                return res
            return self._submit_vasp(material, "relax", res)

        if not material.dir.is_dir():
            res.status = "skipped"
            return res
        niter = QEText.iterations(relax / "scf.out")
        if niter and niter < 3:
            res.say("Structure already relaxed")
            res.status = "skipped"
            return res
        if not copy_file(material.scf_relaxed, relax / "scf.in"):
            res.status = "skipped"
            res.message = f"{material.scf_relaxed.name} missing -- run process 2 first"
            return res
        return self.stage_and_submit(material, "relax", "run-scf.sh", "scf", result=res)

    def fine_scan(self, start: int, end: int, track: str | None = None, *_) -> list[Result]:
        """``mainprogram 5`` -- scf on the fine (doubled) k-mesh, in ``calc/``."""
        self.banner("Submitting QE scf calculations with fine k-mesh, required for "
                    "interpolating el-ph coupling (EPC) matrices")
        results = self.map(self._fine_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _fine_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        if not material.dir.is_dir():
            res.status = "skipped"
            return res
        source = f"scf_dir/scf-{material.mpid}-{material.compound}-fit.in"
        return self.stage_and_submit(material, "calc", "run-scf.sh", "scf",
                                     files={"scf.in": source}, result=res)

    def coarse_scan(self, start: int, end: int, track: str | None = None, *_) -> list[Result]:
        """``mainprogram 6`` -- scf on the coarse k-mesh (keeps the fine results)."""
        self.banner("Submitting QE scf calculation with coarse grid")
        results = self.map(self._coarse_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _coarse_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        if not material.dir.is_dir():
            res.status = "skipped"
            return res
        calc = ensure_dir(material.sub("calc"))
        copy_file(calc / "scf.in", calc / "scf-fit.in")
        copy_file(calc / "scf.out", calc / "scf-fit.out")
        source = f"scf_dir/scf-{material.mpid}-{material.compound}.in"
        return self.stage_and_submit(material, "calc", "run-scf.sh", "scf",
                                     files={"scf.in": source}, result=res)

    # ======================================================================== #
    #  4  create-inputs
    # ======================================================================== #
    def create_inputs(self, start: int, end: int, track: str | None = None,
                      nkpt: int | None = None, *_) -> list[Result]:
        """``mainprogram 4`` -- every downstream QE input from the relaxed scf.out.

        This was the worst shared-state offender of the bash layer: it wrote
        ``mass.dat``, ``qpoint.dat`` (appended!), ``kpoint.dat``, ``BZ.pdf`` and
        ``scf_dir/kpathlines.dat`` into the project root and then let four
        helper scripts read them back by name.  Each material now builds inside
        its own scratch directory, so the loop parallelises safely.
        """
        self.ensure_project_dirs("kpath", "elph_dir", "scf_dir", "matdyn_dir", "q2r_dir")
        self.args = {"nkpt": int(nkpt or self.input.nkpt), "kcut": self.input.kcut}
        self.banner("Creating input files for 'mainprogram process' with process = 5 - 18")
        results = self.map(self._create_inputs_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _create_inputs_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(f"material_id = {material.mpid}, compound = {material.compound}")

        if not material.scf_template.is_file():
            with pushd(self.root):
                self.run_qe_input(material.mpid)
        if not material.dir.is_dir():
            res.status = "skipped"
            res.message = f"{material.name} folder not found"
            return res

        scf_out = material.sub("relax") / "scf.out"
        block = QEText.final_coordinates(scf_out)
        if not block:
            res.status = "skipped"
            res.message = "no relaxed structure in relax/scf.out"
            return res

        template = read_text(material.scf_template)
        relax_in = read_text(material.sub("relax") / "scf.in")
        prefix = QEText.prefix(relax_in) if relax_in else QEText.prefix(template)
        nelec = QEText.nelec(scf_out)
        soc = "lspinorb = .true." in relax_in
        spin = "nspin = 2" in relax_in
        nbnd = self.nbnd_heuristic(nelec, soc=soc, spin=spin, style="qe")
        res.say(f"Number of bands for the band structure run "
                f"[electrons = {nelec}]: {nbnd}")

        kmesh, kshift = QEText.kmesh(template)
        fine = [2 * k for k in kmesh]
        qmesh = self.qmesh_for(kmesh)
        ecut = QEText.value(template, "ecutwfc", "?")
        res.say(f"Ecut: {ecut} Ry, K-mesh: {kmesh}, q-mesh: {qmesh}")

        # ---- namelist headers (bash: three sed pipelines per material) ------
        base = QEText.drop_namelist(QEText.drop_namelist(QEText.header(template),
                                                         "IONS"), "CELL")
        base = QEText.replace(base, "'vc-relax'", "'scf'")

        header_scf = list(base)
        header_band = QEText.replace(base, "'scf'", "'bands'")
        header_band = QEText.drop(header_band, "conv_thr")
        header_band = QEText.insert_after(header_band, "&SYSTEM", f"  nbnd={nbnd},")
        header_band = QEText.insert_after(header_band, "&ELECTRONS", "  conv_thr = 1d-10,")
        header_fit = QEText.insert_after(header_scf, "occupations", "  la2F = .true.,")

        species = QEText.card(template, "ATOMIC_SPECIES")
        kpoint_card = ["K_POINTS automatic",
                       " ".join(str(v) for v in kmesh + kshift)]
        kpoint_fine = ["K_POINTS automatic", " ".join(str(v) for v in fine) + " 0 0 0"]
        ntype = max(0, len(species) - 1)
        masses = [field_of(ln, 2) for ln in species[1:] if ln.strip()]

        mpid, comp = material.mpid, material.compound
        with self.scratch(material, "create-inputs",
                          copy=[f"scf_dir/scf-{mpid}.in"]) as work:
            # helpers such as elph.py(parallel_irr) read R<id>-<comp>/calc/...
            link = work / material.name
            if not link.exists():
                with contextlib.suppress(OSError):
                    link.symlink_to(material.dir, target_is_directory=True)
            ensure_dir(work / "elph_dir")

            sdir = work / "scf_dir"
            write_lines(sdir / f"scf-{mpid}-{comp}.in",
                        header_scf + species + kpoint_card + block)

            # k-path (writes scf_dir/kpathlines.dat, kspecial-points.dat, BZ.pdf)
            kcut = self.args.get("kcut", 0)
            if kcut:
                res.say("Partial k-path mesh will be used for the phonon bandstructure")
            self.run_kpath("point", f"scf_dir/scf-{mpid}-{comp}.in",
                           self.args["nkpt"], kcut, 0)
            kpathlines = read_lines(sdir / "kpathlines.dat")

            write_lines(sdir / f"scf-{mpid}-{comp}-band.in",
                        header_band + species + kpathlines + block)
            fit_lines = header_fit + species + kpoint_fine + block
            write_lines(sdir / f"scf-{mpid}-{comp}-fit.in", fit_lines)

            dos_lines = QEText.replace(fit_lines, "'scf'", "'nscf'")
            dos_lines = QEText.replace(dos_lines, "'smearing'", "'tetrahedra'")
            dos_lines = QEText.drop(dos_lines, "la2F", "degaus", "smearing", "conv_thr")
            dos_lines = QEText.insert_after(dos_lines, "&ELECTRONS", "  conv_thr = 1d-10,")
            write_lines(sdir / f"scf-{mpid}-{comp}-dos.in", dos_lines)

            # the three fixed-name files the helpers read back, now private
            write_lines(work / "mass.dat", masses)
            write_lines(work / "qpoint.dat", [" ".join(str(q) for q in qmesh)])
            write_lines(work / "kpoint.dat",
                        [f"v{material.index} {mpid} {comp} {prefix} "
                         f"{' '.join(str(q) for q in qmesh)} {ntype}"])

            self.run_elph(material, prefix)
            run_helper("band", "main", mpid, comp, prefix)
            run_helper("phonband", "phonband_in", mpid, comp, prefix)
            run_helper("dos", "dos_in", mpid, comp, prefix)
            run_helper("q2r", "q2r_in", mpid, comp, prefix)
            run_helper("matdyn", "matdyn_in", mpid, comp, prefix)
            run_helper("matdyn_dos", "matdyn_dos", mpid, comp, prefix)

            # ---- publish the artefacts --------------------------------------
            for name in (f"scf-{mpid}-{comp}.in", f"scf-{mpid}-{comp}-band.in",
                         f"scf-{mpid}-{comp}-fit.in", f"scf-{mpid}-{comp}-dos.in",
                         f"band-{mpid}-{comp}.in", f"bandproj-{mpid}-{comp}.in",
                         f"phonband-{mpid}-{comp}.in", f"dos-{mpid}-{comp}.in",
                         f"pdos-{mpid}-{comp}.in"):
                copy_file(sdir / name, self.root / "scf_dir" / name)
            self.collect_scratch(work, {
                f"elph-{mpid}-{comp}.in": "elph_dir",
                "elph_dir/*": "elph_dir",
                f"q2r-{mpid}-{comp}.in": "q2r_dir",
                f"dynmat-{mpid}-{comp}.in": "q2r_dir",
                f"matdyn-{mpid}-{comp}.in": "matdyn_dir",
                f"matdyn-{mpid}-{comp}-dos.in": "matdyn_dir",
            })
            move_file(work / "BZ.pdf", self.root / "kpath" / f"BZ-{mpid}-{comp}.pdf")
            copy_file(sdir / "kpathlines.dat",
                      self.root / "kpath" / f"kpath-{mpid}-{comp}.dat")
            copy_file(sdir / "kspecial-points.dat",
                      self.root / "kpath" / f"kspecial-{mpid}-{comp}.dat")
        return res

    # ======================================================================== #
    #  7  ph-scan  (+ the shared el-ph state machine used by checkph)
    # ======================================================================== #
    ELPH_DONE = "done"
    ELPH_WALLTIME = "walltime"
    ELPH_UNCONVERGED = "unconverged"
    ELPH_FRESH = "fresh"
    ELPH_NOT_STARTED = "not-started"

    def elph_status(self, material: Material, folder: str = "calc") -> dict:
        """Classify one el-ph run.  One implementation for ph-scan + checkph.

        The bash versions of this classification disagreed with each other and
        mis-filed the most common failure: ``dync == 0`` (only ``.dyn0``
        written, i.e. the job died at start-up) fell into the "unconverged"
        branch and got ``alpha_mix`` instead of ``recover=.true.``.
        """
        calc = material.sub(folder)
        elph_out = calc / "elph.out"
        info = {"state": self.ELPH_NOT_STARTED, "dyn_done": 0, "dyn_total": 0,
                "job_done": False, "unconverged": False, "qbreaksym": False,
                "fft": False, "scale_sym_ops": False, "lambda": False}
        if not elph_out.is_file():
            return info

        info["unconverged"] = grep_count("No convergence has been achieved", elph_out) > 0
        info["job_done"] = grep_count("JOB DONE", elph_out) > 0
        info["qbreaksym"] = grep_count("q-mesh breaks symmetry", elph_out) > 0
        info["fft"] = (grep_count("FFT grid incompatible with symmetry", elph_out) > 0
                       or grep_count("incompatible FFT grid", elph_out) > 0)
        info["scale_sym_ops"] = grep_count("Error in routine scale_sym_ops (8):",
                                           elph_out) > 0
        info["lambda"] = (calc / "lambda.out").is_file()

        dyn0 = calc / f"{material.compound}.dyn0"
        # bash: `[[ 0 -eq "" ]]` is TRUE, so a missing dyn0 read as "completed".
        info["dyn_total"] = to_int(read_lines(dyn0)[1], -1) if dyn0.is_file() and \
            len(read_lines(dyn0)) > 1 else -1
        info["dyn_done"] = len(glob.glob(os.fspath(calc / "*.dyn*"))) - 1

        if info["dyn_total"] >= 0 and info["dyn_done"] == info["dyn_total"]:
            info["state"] = self.ELPH_DONE
        elif info["unconverged"]:
            info["state"] = self.ELPH_UNCONVERGED
        elif info["dyn_done"] > 0 and info["dyn_done"] < info["dyn_total"]:
            info["state"] = self.ELPH_WALLTIME
        elif info["dyn_done"] <= 0:
            info["state"] = self.ELPH_FRESH
        else:
            info["state"] = self.ELPH_UNCONVERGED
        return info

    def ph_scan(self, start: int, end: int, track: str | None = None, *_) -> list[Result]:
        """``mainprogram 7`` -- submit/restart the ``ph.x`` el-ph calculation."""
        self.banner("Submitting/resubmitting el-ph coupling calculations")
        results = self.map(self._ph_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _ph_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        elph_dir = self.root / "elph_dir"
        tag = f"{material.mpid}-{material.compound}"

        if (elph_dir / f"PARALLEL_q-{tag}").is_file():
            return self._ph_parallel_q(material, res)
        if (elph_dir / f"PARALLEL_irr-{tag}").is_file():
            return self._ph_parallel_irr(material, res)
        return self._ph_serial(material, res)

    def _ph_serial(self, material: Material, res: Result) -> Result:
        calc = material.sub("calc")
        if not material.dir.is_dir():
            res.status = "skipped"
            return res
        master = self.root / "elph_dir" / f"elph-{material.mpid}-{material.compound}.in"
        if (self.root / "elph_dir" /
                f"{material.mpid}-{material.compound}-freq.dat").is_file():
            res.say("imaginary frequencies recorded -- staging a larger smearing input")
            copy_file(master, calc / "elph.in")
            copy_file(self.root / "run-elph.sh", calc / "run-elph.sh")
            res.status = "skipped"
            return res

        info = self.elph_status(material)
        state = info["state"]
        if state == self.ELPH_DONE:
            res.say("el-ph calculation completed, do nothing")
            res.status = "skipped"
            return res

        if state == self.ELPH_WALLTIME:
            res.say("calculation not completed due to walltime -- recover=.true.")
            lines = QEText.drop(read_lines(calc / "elph.in"), "recover")
            write_lines(calc / "elph.in",
                        QEText.insert_after(lines, "&inputph", "  recover=.true.,"))
            nelph = len(glob.glob(os.fspath(calc / "elph.out*")))
            copy_file(calc / "elph.out", calc / f"elph.out{nelph}")
        elif state in (self.ELPH_FRESH, self.ELPH_NOT_STARTED):
            res.say(f" submitting fresh jobs in {material.name}")
            if not copy_file(master, calc / "elph.in"):
                res.status = "skipped"
                res.message = f"{master.name} missing -- run process 4 first"
                return res
        else:
            res.say("restarting old unconverged job")
            lines = read_lines(master)
            if any("alpha_mix" in ln for ln in lines):
                res.say("Change nmix_ph in <QE>/PHonon/PH/phq_readin.f90 "
                        "(comment the nmix_ph range check) to use nmix_ph=8")
                lines = QEText.drop(lines, "alpha_mix")
                write_lines(master, lines)
                lines = QEText.insert_after(lines, "&inputph",
                                            "  alpha_mix=0.3, nmix_ph=8,")
            else:
                lines = QEText.insert_after(lines, "&inputph", "  alpha_mix=0.3,")
            # the master input in elph_dir stays pristine; only the run copy changes
            write_lines(calc / "elph.in", lines)

        return self.stage_and_submit(material, "calc", "run-elph.sh", "elph", result=res)

    def _ph_parallel_q(self, material: Material, res: Result) -> Result:
        """``elph_dir/PARALLEL_q-*`` -- one job per q point."""
        calc = material.sub("calc")
        tag = f"{material.mpid}-{material.compound}"
        marker = self.root / "elph_dir" / f"PARALLEL_q-{tag}"
        nq = to_int(field_of(read_lines(marker)[0] if read_lines(marker) else "", 1), 0)

        dyn0 = calc / f"{material.compound}.dyn0"
        incomplete = True
        if dyn0.is_file():
            sizes = [len(read_lines(p))
                     for p in sorted(glob.glob(os.fspath(calc / f"{material.compound}.dyn*")))]
            incomplete = any(size == 0 for size in sizes[:-1]) if sizes else True

        if not incomplete:
            res.say("Calculations completed. Resubmitting without start_q/last_q")
            lines = QEText.drop(read_lines(calc / "elph_1.in"), "start_q", "last_q")
            write_lines(calc / "elph.in",
                        QEText.insert_after(lines, "&inputph", "  recover=.true.,"))
            return self.stage_and_submit(material, "calc", "run-elph.sh", "elph",
                                         result=res)

        for iq in range(1, nq + 1):
            out = calc / f"elph_{iq}.out"
            if out.is_file() and grep_count("JOB DONE", out) > 0:
                res.say(f"Convergence has been achieved for q{iq}")
                continue
            if out.is_file():
                nelph = len(glob.glob(os.fspath(calc / f"elph_{iq}.out*")))
                copy_file(calc / f"elph_{iq}.in", calc / f"elph_{iq}-{nelph}.in")
                copy_file(out, calc / f"elph_{iq}-{nelph}.out")
                lines = read_lines(calc / f"elph_{iq}.in")
                write_lines(calc / f"elph_{iq}.in",
                            QEText.insert_after(lines, "&inputph", "  recover=.true.,"))
            else:
                copy_file(self.root / "elph_dir" / f"elph-{tag}-{iq}.in",
                          calc / f"elph_{iq}.in")
            self._write_run_script(calc, "run-elph.sh", "elph", f"elph_{iq}")
            res.job = self.scheduler.submit("run-elph.sh", calc, material.mpid,
                                            material.compound, "elph")
        return res

    def _write_run_script(self, target: Path, script: str, old: str, new: str) -> None:
        """``sed "s/<old>/<new>/g" run-X.sh`` -- retarget a submission script."""
        source = target / "temp.sh"
        if not source.is_file():
            copy_file(self.root / script, source)
        write_lines(target / script,
                    [ln.replace(old, new) for ln in read_lines(source)])

    def _ph_parallel_irr(self, material: Material, res: Result) -> Result:
        """``elph_dir/PARALLEL_irr-*`` -- one job per (q, irrep), then reassembly."""
        calc = material.sub("calc")
        comp = material.compound
        tag = f"{material.mpid}-{comp}"
        marker = self.root / "elph_dir" / f"PARALLEL_irr-{tag}"
        irr_of_q = {to_int(field_of(ln, 1).lstrip("v")): to_int(field_of(ln, 2))
                    for ln in read_lines(marker) if ln.strip()}
        res.say("parallel over q and irr")

        if glob.glob(os.fspath(calc / f"{comp}.dyn*")):
            copy_file(calc / f"{comp}.dyn0", calc / "only_init_dynmat")
            copy_file(calc / "elph.out", calc / "only_init_out")
            copy_file(calc / "elph.in", calc / "elph_init_out.in")

        submitted = calc / "ELPH_Q_IR_SUBMITTED"
        if not submitted.is_file():
            submitted.touch()
            for iq, nirr in sorted(irr_of_q.items()):
                for irr in range(1, nirr + 1):
                    res.say(f"q-point: {iq}, irreducible representation: {irr}")
                    copy_file(self.root / "elph_dir" / f"elph-{tag}-{iq}-{irr}.in",
                              calc / f"elph-{iq}-{irr}.in")
                    self._write_run_script(calc, "run-elph.sh", "elph", f"elph-{iq}-{irr}")
                    sub = calc / f"{iq}-{irr}"
                    remove(sub, recursive=True)
                    ensure_dir(sub / "_ph0" / f"{comp}.phsave")
                    for wfc in sorted(glob.glob(os.fspath(calc / f"{comp}.wfc*"))):
                        with contextlib.suppress(OSError):
                            (sub / Path(wfc).name).symlink_to(Path(wfc).resolve())
                    for shared in (f"{comp}.save", f"{comp}.xml"):
                        with contextlib.suppress(OSError):
                            (sub / shared).symlink_to((calc / shared).resolve())
                    copy_file(calc / "_ph0" / f"{comp}.phsave",
                              sub / "_ph0" / f"{comp}.phsave")
                    if iq > 1 and irr > 2:
                        with open(calc / "run-elph.sh", "a") as handle:
                            handle.write(f"rm -r {iq}-{irr}/_ph0/{comp}.q_{iq}\n")
                    res.job = self.scheduler.submit("run-elph.sh", calc, material.mpid,
                                                    comp, "elph")
            return res

        # ---- reassembly pass -------------------------------------------------
        phsave = calc / "_ph0" / f"{comp}.phsave"
        ensure_dir(phsave)
        for iq, nirr in sorted(irr_of_q.items()):
            if iq > 1:
                res.say(f"Replacing _ph0/{comp}.q_{iq} with {iq}-1/_ph0/{comp}.q_{iq}")
                remove(calc / "_ph0" / f"{comp}.q_{iq}", recursive=True)
                if not copy_file(calc / f"{iq}-1" / "_ph0" / f"{comp}.q_{iq}",
                                 calc / "_ph0" / f"{comp}.q_{iq}"):
                    copy_file(calc / f"{iq}-2" / "_ph0" / f"{comp}.q_{iq}",
                              calc / "_ph0" / f"{comp}.q_{iq}")
            else:
                copy_file(calc / "1-1" / "_ph0" / f"{comp}.aldv1",
                          calc / "_ph0" / f"{comp}.aldv1")
            for irr in range(1, nirr + 1):
                src = calc / f"{iq}-{irr}" / "_ph0" / f"{comp}.phsave"
                for name in (f"dynmat.{iq}.{irr}.xml", f"elph.{iq}.{irr}.xml"):
                    if not copy_file(src / name, phsave / name):
                        res.say(f"{name} not found in {src}")
            copy_file(calc / f"{iq}-1" / "_ph0" / f"{comp}.phsave" / f"dynmat.{iq}.0.xml",
                      phsave / f"dynmat.{iq}.0.xml")
        if not (calc / "_ph0" / f"{comp}.aldv1").is_file():
            copy_file(calc / "1-2" / "_ph0" / f"{comp}.aldv1",
                      calc / "_ph0" / f"{comp}.aldv1")
        copy_file(calc / "1-1" / "_ph0" / f"{comp}.phsave" / "tensors.xml",
                  phsave / "tensors.xml")

        lines = QEText.drop(read_lines(calc / "elph-1-1.in"),
                            "start_q", "last_q", "start_irr", "last_irr", "outdir")
        write_lines(calc / "elph.in",
                    QEText.insert_after(lines, "&inputph", "  outdir='./',"))
        return self.stage_and_submit(material, "calc", "run-elph.sh", "elph", result=res)

    # ======================================================================== #
    #  checkph / checkfreq
    # ======================================================================== #
    def phcheck_scan(self, start: int, end: int, track: str | None = None,
                     *_) -> list[Result]:
        """``mainprogram checkph`` -- classify every el-ph run, write the 4 lists."""
        self.banner("Checking status of el-ph calculations of process = 7")
        materials = self.materials(start, end, track)
        results = self.map(self._phcheck_one, materials)

        buckets = {"finished": [], "finished-no-lambda": [],
                   "not-converged": [], "not-completed": []}
        counters = {"total": 0, "done": 0, "unconverged": 0, "qbreaksym": 0,
                    "scale_sym_ops": 0, "fft": 0}
        for res in results:
            rows = res.rows
            if not rows:
                continue
            counters["total"] += 1
            for key in ("done", "unconverged", "qbreaksym", "scale_sym_ops", "fft"):
                counters[key] += int(bool(rows.get(key)))
            for bucket in buckets:
                if rows.get(bucket):
                    buckets[bucket].append((res.mpid, res.compound))

        LOG.info("######################### Check these files ######################")
        names = {"finished": "mpid-list-elph-finished.in",
                 "finished-no-lambda": "mpid-list-elph-finished-no-lambda.in",
                 "not-converged": "mpid-list-elph-not-converged.in",
                 "not-completed": "mpid-list-elph-not-completed.in"}
        for bucket, rows in buckets.items():
            # only write a list that has content -- bash `touch`ed the temp files
            # unconditionally and therefore overwrote yesterday's lists with
            # empty ones on every run.
            if rows:
                write_track_file(self.root / names[bucket], rows)
                LOG.info("  %s  (%d entries)", names[bucket], len(rows))

        not_completed = counters["total"] - counters["done"] - counters["unconverged"] \
            - counters["qbreaksym"] - counters["scale_sym_ops"] - counters["fft"]
        LOG.info("*" * 107)
        LOG.info("Total calculations: %d", counters["total"])
        LOG.info("Total completed: %d", counters["done"])
        LOG.info("Not converged: %d", counters["unconverged"])
        LOG.info("qbreaksym error: %d", counters["qbreaksym"])
        LOG.info("scale_sym_ops error: %d", counters["scale_sym_ops"])
        LOG.info("FFT grid incompatible error: %d", counters["fft"])
        LOG.info("Calculations not completed: %d", max(0, not_completed))
        LOG.info("all done")
        LOG.info("*" * 107)
        return results

    def _phcheck_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        tag = f"{material.mpid}-{material.compound}"
        calc = material.sub("calc")

        if (self.root / "elph_dir" / f"PARALLEL_irr-{tag}").is_file():
            marker = self.root / "elph_dir" / f"PARALLEL_irr-{tag}"
            for line in read_lines(marker):
                iq, nirr = to_int(field_of(line, 1).lstrip("v")), to_int(field_of(line, 2))
                done = sum(1 for irr in range(1, nirr + 1)
                           if grep_count("Convergence", calc / f"elph-{iq}-{irr}.out") > 0)
                xml = len(glob.glob(os.fspath(calc / f"{iq}-*" / "_ph0" / "*" /
                                              f"elph.{iq}.*.xml")))
                res.say(f"{done}/{nirr} irreps finished for q-point {iq}, "
                        f"{xml} elph*.xml files found")
            return res

        if not material.dir.is_dir():
            res.status = "skipped"
            return res
        res.say(f"Checking status for: v{material.index} {material.mpid} "
                f"{material.compound}")
        info = self.elph_status(material)
        if info["state"] == self.ELPH_NOT_STARTED:
            res.say(f"STATUS: el-ph coupling calculation has not started for {tag}")
            res.status = "skipped"
            return res

        rows: dict[str, Any] = {}
        nbnd = last_match_field("number of Kohn-Sham states=", calc / "scf.out", 5)
        nelm = last_match_field("number of electrons", material.sub("relax") / "scf.out", 5)
        res.say(f"nband: {nbnd}, nelm: {nelm}")

        if info["qbreaksym"]:
            res.say("STATUS: q break symmetry error occured")
            write_lines(self.root / "scf_dir" / f"{material.mpid}-qbreaksym", ["qbreaksym"])
            rows["qbreaksym"] = True
        if info["scale_sym_ops"]:
            ensure_dir(self.root / "scale_sym_ops")
            write_lines(self.root / "scale_sym_ops" / f"{material.mpid}-scale_sym_ops",
                        ["scale_sym_ops error"])
            rows["scale_sym_ops"] = True
        if info["fft"]:
            res.say("STATUS: FFT grid incompatible with symmetry -- raise ecutrho/ecutwfc")
            write_lines(self.root / "scf_dir" / f"{material.mpid}-fft-grid", ["FFT grid"])
            rows["fft"] = True

        if info["state"] == self.ELPH_DONE:
            rows["done"] = True
            if info["job_done"]:
                rows["finished"] = True
                if not info["lambda"]:
                    res.say("STATUS: el-ph calculation completed but lambda.out not found")
                    rows["finished-no-lambda"] = True
                else:
                    res.say("STATUS: el-ph calculation completed and lambda.out found")
            else:
                res.say("STATUS: el-ph calculation is completing soon")
        elif info["state"] == self.ELPH_WALLTIME:
            if self.scheduler.is_running(calc):
                res.say("STATUS: Calculation is in progress")
            else:
                res.say("STATUS: calculation not completed due to walltime")
                rows["not-completed"] = True
        elif info["state"] == self.ELPH_FRESH:
            res.say("STATUS: submit a fresh job")
            # bash tested `-f _ph0` for a *directory*, so the stale tree survived
            remove(calc / "_ph0", recursive=True)
        if info["unconverged"]:
            res.say("STATUS: not converged")
            rows["unconverged"] = True
            rows["not-converged"] = True

        res.say(f"{info['dyn_done']} out of {info['dyn_total']} .dyn files present")
        res.rows = rows
        return res

    def checkfreq_scan(self, start: int, end: int, track: str | None = None,
                       *_) -> list[Result]:
        """``mainprogram checkfreq`` -- flag imaginary phonon frequencies."""
        self.banner("Checking imaginary frequencies")
        results = self.map(self._checkfreq_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _checkfreq_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        calc = material.sub("calc")
        freq_gp = calc / f"{material.compound}.freq.gp"
        target = self.root / "elph_dir" / \
            f"{material.mpid}-{material.compound}-freq.dat"

        with self.scratch(material, "checkfreq", local_scf_dir=False) as work:
            if freq_gp.is_file():
                run_helper("checkfreq", "main", os.fspath(freq_gp))
            else:
                res.say(f"{material.mpid}: {freq_gp.name} not present")
            if (calc / "lambda.out").is_file() and not (calc / "freq.plot").is_file():
                write_lines(work / "freq.dat", ["0 -100"])
            if (work / "freq.dat").is_file():
                ensure_dir(target.parent)
                move_file(work / "freq.dat", target)
                res.say(f"  imaginary frequencies recorded in {target.name}")
        return res

    # ======================================================================== #
    #  8 - 12  post-phonon steps (q2r / matdyn / matdyn-dos / lambda / phonband)
    # ======================================================================== #
    def q2r_scan(self, start: int, end: int, track: str | None = None, *_) -> list[Result]:
        """``mainprogram 8`` -- real-space force constants."""
        self.banner("Submitting q2r.x calculations for force constants in real space")
        results = self.map(self._q2r_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _q2r_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        folder = material.phonon_folder("phonon")
        info = self.elph_status(material, folder)
        # bash overrode dync/netdyn with 1 == 1 two lines later, which disabled
        # this gate completely; it is a real check again.
        if info["state"] != self.ELPH_DONE:
            res.say("electron-phonon calculation not finished")
            res.status = "skipped"
            return res
        source = self.root / "q2r_dir" / f"q2r-{material.mpid}-{material.compound}.in"
        target = ensure_dir(material.sub(folder)) / "q2r.in"
        if folder == "phonon":
            write_lines(target, QEText.drop(read_lines(source), "la2F"))
        elif not copy_file(source, target):
            res.status = "skipped"
            res.message = f"{source.name} missing"
            return res
        return self.stage_and_submit(material, folder, "run-q2r.sh", "q2r", result=res)

    def matdyn_scan(self, start: int, end: int, track: str | None = None,
                    *_) -> list[Result]:
        """``mainprogram 9`` -- phonon bands along the high-symmetry path."""
        self.banner("Computing phonon bandstructure")
        results = self.map(self._matdyn_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _matdyn_one(self, material: Material) -> Result:
        return self._matdyn_generic(material, "matdyn.in", "matdyn", "run-matdyn.sh",
                                    f"matdyn-{material.mpid}-{material.compound}.in")

    def matdyn_dos_scan(self, start: int, end: int, track: str | None = None,
                        *_) -> list[Result]:
        """``mainprogram 10`` -- phonon DOS and linewidths."""
        self.banner("Submitting phonon dos and phonon linewidth calculations")
        results = self.map(self._matdyn_dos_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _matdyn_dos_one(self, material: Material) -> Result:
        return self._matdyn_generic(
            material, "matdyn-dos.in", "matdyn-dos", "run-matdyn-dos.sh",
            f"matdyn-{material.mpid}-{material.compound}-dos.in")

    def _matdyn_generic(self, material: Material, target_name: str, tag: str,
                        script: str, source_name: str) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        if not material.dir.is_dir():
            res.status = "skipped"
            return res
        folder = material.phonon_folder("phonon")
        source = self.root / "matdyn_dir" / source_name
        target = ensure_dir(material.sub(folder)) / target_name
        if folder == "phonon":
            write_lines(target, QEText.drop(read_lines(source), "la2F"))
        elif not copy_file(source, target):
            res.status = "skipped"
            res.message = f"{source_name} missing"
            return res
        return self.stage_and_submit(material, folder, script, tag, result=res)

    def dynmat_scan(self, start: int, end: int, track: str | None = None,
                    *_) -> list[Result]:
        """``mainprogram 23`` -- diagonalise the dynamical matrix (``dynmat.axsf``)."""
        self.banner("Diagonalizing phonon dynamical matrices")
        results = self.map(self._dynmat_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _dynmat_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        folder = material.phonon_folder("calc")
        target = material.sub(folder)
        comp = material.compound
        if not ((target / f"{comp}.dyn1").is_file() or (target / f"{comp}.dyn").is_file()):
            res.say("electron-phonon calculation not finished")
            res.status = "skipped"
            return res

        prefix = self.prefix_of(material)
        with self.scratch(material, "dynmat", local_scf_dir=False) as work:
            run_helper("q2r", "q2r_in", material.mpid, comp, prefix)
            move_file(work / f"dynmat-{material.mpid}-{comp}.in", target / "dynmat.in")
            # bash left q2r-<id>-<comp>.in lying around in the project root
            move_file(work / f"q2r-{material.mpid}-{comp}.in",
                      self.root / "q2r_dir" / f"q2r-{material.mpid}-{comp}.in")
        copy_file(target / f"{comp}.dyn1", target / f"{comp}.dyn")
        return self.stage_and_submit(material, folder, "run-dynmat.sh", "dynmat",
                                     result=res)

    def lambda_scan(self, start: int, end: int, track: str | None = None,
                    qgauss: float = 0.12, smearing: int = 0, mustar: float = 0.16,
                    *_) -> list[Result]:
        """``mainprogram 11`` -- run ``lambda.x`` on the login node."""
        self.banner("Running lambda.x command")
        self.args = {"qgauss": qgauss, "smearing": smearing, "mustar": mustar}
        results = self.map(self._lambda_one, self.materials(start, end, track))
        LOG.info("all done")
        LOG.info("************* Other parameters used *****************")
        LOG.info("Smearing for q-mesh: %s", qgauss)
        LOG.info("Smearing type: %s (0 = gauss, 1 = MP)", smearing)
        LOG.info("Coloumb potential (mu_star): %s", mustar)
        LOG.info("*****************************************************")
        return results

    def _lambda_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        calc = material.sub("calc")
        if not calc.is_dir():
            res.status = "skipped"
            res.message = f"no {material.name}/calc directory"
            return res
        emax = QEText.max_phonon_freq(calc / "elph.out") + 5
        res.say(f"Max phonon freq: {emax}")
        with pushd(calc):
            run_helper("lambda_in", "main", material.compound, emax,
                       self.args["qgauss"], self.args["smearing"], self.args["mustar"])
            if self.dry_run:
                res.say("[dry-run] lambda.x < lambda.in > lambda.out")
                return res
            try:
                with open("lambda.in") as stdin, open("lambda.out", "w") as stdout:
                    subprocess.run(["lambda.x"], stdin=stdin, stdout=stdout,
                                   stderr=subprocess.STDOUT, check=False)
            except FileNotFoundError:
                res.status = "failed"
                res.message = "lambda.x not found on PATH"
        return res

    def phonband_scan(self, start: int, end: int, track: str | None = None,
                      *_) -> list[Result]:
        """``mainprogram 12`` -- process phonon bands and mode-resolved lambda."""
        self.banner("Processing phonon bandstructure")
        results = self.map(self._phonband_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _phonband_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        if not material.dir.is_dir():
            res.status = "skipped"
            return res
        calc = ensure_dir(material.sub("calc"))
        copy_file(self.root / "scf_dir" /
                  f"phonband-{material.mpid}-{material.compound}.in",
                  calc / "phonband.in")
        with pushd(calc):
            run_helper("freq_process", "freq_process", material.compound)
            write_lines(Path("gammaband.in"),
                        ["elph.gamma.2", "0 5000", "gamma.plot", "gamma.ps",
                         "0.0", "100 0"])
            if self.dry_run:
                res.say("[dry-run] plotband.x < gammaband.in > gammaband.out")
                return res
            try:
                with open("gammaband.in") as stdin, open("gammaband.out", "w") as stdout:
                    subprocess.run(["plotband.x"], stdin=stdin, stdout=stdout,
                                   stderr=subprocess.STDOUT, check=False)
            except FileNotFoundError:
                res.status = "failed"
                res.message = "plotband.x not found on PATH"
        return res

    # ======================================================================== #
    #  13 - 18  bands and DOS
    # ======================================================================== #
    def bandscf_scan(self, start: int, end: int, track: str | None = None,
                     *_) -> list[Result]:
        """``mainprogram 13`` -- stage ``bands/`` and submit the scf step."""
        self.banner("Preparing and submitting BAND and DOS calculations")
        self.args = {"nkpt": self.input.nkpt, "kcut": self.input.kcut}
        results = self.map(self._bandscf_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _bandscf_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        bands = ensure_dir(material.sub("bands"))
        relax = material.sub("relax")

        if not self.is_vasp:
            mpid, comp = material.mpid, material.compound
            return self.stage_and_submit(
                material, "bands", "run-scf.sh", "scf",
                files={"scf.in": f"scf_dir/scf-{mpid}-{comp}.in",
                       "scf-band.in": f"scf_dir/scf-{mpid}-{comp}-band.in",
                       "band.in": f"scf_dir/band-{mpid}-{comp}.in"},
                result=res)

        # ---- VASP ------------------------------------------------------------
        copy_file(self.root / "run-vasp.sh", bands / "run.sh")
        if not copy_file(relax / "CONTCAR", bands / "POSCAR"):
            copy_file(relax / "POSCAR", bands / "POSCAR")
        copy_file(relax / "POTCAR", bands / "POTCAR")
        copy_file(relax / "INCAR", bands / "INCAR")
        copy_file(self.root / "vdw_kernel.bindat", bands / "vdw_kernel.bindat")

        with pushd(bands):
            res.say(f"Creating soft links to CHGCAR of {material.name}/relax/")
            with contextlib.suppress(OSError):
                Path("CHGCAR").symlink_to(Path("../relax/CHGCAR"))

            nelec = to_int(last_match_field("NELECT", relax / "OUTCAR", 3), 0)
            template = self.root / "vasp-band.in"
            if template.is_file():
                copy_file(template, Path("vasp.in"))
                soc = grep_count("LSORBIT .TRUE.", template) > 0
                nbnd = to_int(last_match_field("NBANDS", "vasp.in", 2), 0)
                res.say(f"Number of band used: {nbnd}")
            else:
                soc = grep_count("LSORBIT = .TRUE.", "INCAR") > 0
                nbnd = self.nbnd_heuristic(nelec, soc=soc, style="vasp")
                res.say(f"vasp-band.in not present; creating vasp.in with NBANDS={nbnd}")
                write_lines(Path("vasp.in"),
                            ["LORBIT 11", "LCHARG .False.", "LWAVE .False.",
                             "NSW 0", "ISIF 2", "ICHARG 11", f"NBANDS {nbnd}",
                             "EDIFFG", "IBRION", "NELM"])

            metagga = grep_count("METAGGA", "vasp.in") > 0
            hybrid = grep_count("LHFCALC", "vasp.in") > 0
            if metagga:
                res.say("METAGGA tag detected")
                write_lines(Path("vasp.in"), QEText.drop(read_lines("vasp.in"), "ICHARG"))
                copy_file(relax / "IBZKPT", Path("KPOINTS"))
            elif hybrid:
                res.say("LHFCALC tag detected -- linking WAVECAR of relax/")
                lines = QEText.drop(read_lines("vasp.in"), "ICHARG") + ["NELMIN 3"]
                write_lines(Path("vasp.in"), lines)
                with contextlib.suppress(OSError):
                    Path("WAVECAR").symlink_to(Path("../relax/WAVECAR"))
                copy_file(relax / "IBZKPT", Path("KPOINTS"))
            else:
                res.say("GGA or LDA band calculations")

            self.run_vasp_process("POSCAR")
            nkpt, kcut = self.args["nkpt"], self.args["kcut"]

            if metagga or hybrid:
                res.say("Remember to turn on LASPH = .True. in the GGA run")
                copy_file(self.root / "input.in", Path("input.in"))
                self.run_kpath("point", "POSCAR", nkpt, kcut, 0)
                kpath_lines = read_lines("scf_dir/kpathlines.dat")[2:]
                if Path("KPT_OPT").is_file():
                    write_lines(Path("INCAR"), QEText.drop(read_lines("INCAR"), "NPAR"))
                    copy_file(relax / "KPOINTS", Path("KPOINTS"))
                    self.run_kpath("line", 50 if hybrid else 100)
                else:
                    copy_file(relax / "IBZKPT", Path("KPOINTS"))
                    kpoints = read_lines("KPOINTS")
                    old = to_int(kpoints[1] if len(kpoints) > 1 else "0", 0)
                    merged = [kpoints[0], str(old + nkpt)] + kpoints[2:] + kpath_lines
                    write_lines(Path("KPOINTS"), merged)
                remove(Path("input.in"))
            else:
                self.run_kpath("line", 100)

            res.say(f" submitting jobs in {material.name}")
            res.job = self.scheduler.submit("run.sh", Path.cwd(), material.mpid,
                                            material.compound, None)
        return res

    def band_scan(self, start: int, end: int, track: str | None = None, *_) -> list[Result]:
        """``mainprogram 14`` -- run the QE ``bands`` calculation."""
        self.banner("Submitting electronic bandstructure calculations")
        if self.is_vasp:
            LOG.info("VASP band structure is produced by process 13 -- nothing to do")
            return []
        results = self.map(self._band_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _band_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        ensure_dir(material.sub("bands"))
        return self.stage_and_submit(material, "bands", "run-band.sh", "band", result=res)

    def bandp_scan(self, start: int, end: int, track: str | None = None,
                   *_) -> list[Result]:
        """``mainprogram 15`` -- band post-processing (``bands.x`` / EIGENVAL)."""
        self.banner("Band structure processing")
        results = self.map(self._bandp_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _bandp_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        bands = ensure_dir(material.sub("bands"))

        if self.is_vasp:
            with pushd(bands):
                move_file(Path("vasp.in"), Path("vasp-old.in"))
                write_lines(Path("vasp.in"), ["NSW 0"])
                self.run_vasp_process("POSCAR")
                self.run_vasp_process("eigen")
                move_file(Path("band.dat.gnu"), Path(f"{material.compound}.dat.gnu"))
            return res

        return self.stage_and_submit(
            material, "bands", "run-bandp.sh", "bandp",
            files={"projwfc.in":
                   f"scf_dir/bandproj-{material.mpid}-{material.compound}.in"},
            result=res)

    def dos_scan(self, start: int, end: int, track: str | None = None, *_) -> list[Result]:
        """``mainprogram 16`` -- stage and submit the DOS run."""
        self.banner("Submitting DOS calculations")
        results = self.map(self._dos_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _dos_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        dos = ensure_dir(material.sub("dos"))
        relax = material.sub("relax")
        copy_file(self.root / "vdw_kernel.bindat", dos / "vdw_kernel.bindat")

        if self.is_vasp:
            copy_file(self.root / "run-vasp.sh", dos / "run.sh")
            if not copy_file(relax / "CONTCAR", dos / "POSCAR"):
                copy_file(relax / "POSCAR", dos / "POSCAR")
            for name in ("POTCAR", "INCAR", "KPOINTS"):
                copy_file(relax / name, dos / name)
            with pushd(dos):
                if not Path("vasp.in").is_file():
                    res.say("vasp.in not present -- creating one")
                    write_lines(Path("vasp.in"),
                                ["LORBIT 11", "LCHARG .True.", "LWAVE .True.",
                                 "NSW 0", "NELM 200", "ISIF 2", "ISMEAR -5",
                                 "NEDOS 3000", "EDIFFG", "SIGMA"])
                self.run_vasp_process("POSCAR")
                kpoints = read_lines("KPOINTS")
                if len(kpoints) > 3:
                    mesh = [to_int(t, 1) * 2 for t in kpoints[3].split()[:3]]
                    kpoints[3] = " ".join(str(v) for v in mesh)
                    write_lines(Path("KPOINTS"), kpoints)
                res.say(f" submitting jobs in {material.name}")
                res.job = self.scheduler.submit("run.sh", Path.cwd(), material.mpid,
                                                material.compound, None)
            return res

        mpid, comp = material.mpid, material.compound
        # bash copied the multi-GB <prefix>.save tree; a symlink does the job
        # (ph-scan already used ln -s for exactly this).
        save = material.sub("bands") / f"{comp}.save"
        if save.is_dir() and not (dos / f"{comp}.save").exists():
            with contextlib.suppress(OSError):
                (dos / f"{comp}.save").symlink_to(save.resolve(), target_is_directory=True)
        return self.stage_and_submit(
            material, "dos", "run-dos.sh", "dos",
            files={"scf-dos.in": f"scf_dir/scf-{mpid}-{comp}-dos.in",
                   "dos.in": f"scf_dir/dos-{mpid}-{comp}.in",
                   "pdos.in": f"scf_dir/pdos-{mpid}-{comp}.in"},
            result=res)

    def dosp_scan(self, start: int, end: int, track: str | None = None,
                  *_) -> list[Result]:
        """``mainprogram 17`` -- ``dos.x`` post-processing."""
        self.banner("Calculations for DOS processing")
        results = self.map(self._dosp_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _dosp_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        ensure_dir(material.sub("dos"))
        return self.stage_and_submit(material, "dos", "run-dosp.sh", "dosp", result=res)

    def pdos_scan(self, start: int, end: int, track: str | None = None,
                  *_) -> list[Result]:
        """``mainprogram 18`` -- ``projwfc.x`` partial DOS."""
        self.banner("Submitting partial DOS calculations")
        results = self.map(self._pdos_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _pdos_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        ensure_dir(material.sub("dos"))
        return self.stage_and_submit(material, "dos", "run-pdos.sh", "pdos", result=res)

    # ======================================================================== #
    #  19 - 21  plots, clean-up, extraction
    # ======================================================================== #
    def plot_scan(self, start: int, end: int, track: str | None = None,
                  nkpt: int | None = None, plot: str | None = None,
                  *_) -> list[Result]:
        """``mainprogram 19`` -- one PDF per material per requested plot type."""
        ensure_dir(self.root / "plots")
        types = [plot] if plot else list(self.input.plot_types)
        materials = self.materials(start, end, track)
        results: list[Result] = []
        for kind in types:
            self.args = {"nkpt": int(nkpt or self.input.nkpt),
                         "kcut": self.input.kcut, "plot": kind}
            self.banner(f"Plotting: {kind}")
            results += self.map(self._plot_one, materials)
        LOG.info("all done. Check inside plots/ folder for output pdf")
        return results

    def _plot_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        kind = self.args["plot"]
        kcut, nkpt = self.args["kcut"], self.args["nkpt"]
        mpid, comp = material.mpid, material.compound
        plots = ensure_dir(self.root / "plots")

        def in_dir(sub: str):
            target = material.sub(sub)
            return target if target.is_dir() else None

        if kind in ("eband", "vasp-line"):
            target = in_dir("bands")
            if target is None:
                res.status = "skipped"
                return res
            with pushd(target):
                self.run_plot("band", mpid, comp, nkpt, kcut)
                move_file(Path(f"{comp}-band.pdf"), Path(f"{mpid}-{comp}-band.pdf"))
                copy_file(Path(f"{mpid}-{comp}-band.pdf"), plots)

        elif kind == "phband":
            for sub in ("phonon", "calc"):
                target = in_dir(sub)
                if target is None:
                    continue
                copy_file(self.root / "config.json", target / "config.json")
                with pushd(target):
                    if sub == "phonon":
                        run_helper("freq_process", "freq_process", comp)
                    self.run_plot("phonband", mpid, comp, nkpt, kcut)
                    move_file(Path(f"{comp}-phonon.pdf"), plots / f"{mpid}-{comp}-phonon.pdf")
                    if sub == "calc":
                        self.run_plot("a2f", mpid, comp)
                        move_file(Path(f"{comp}-a2f.pdf"), plots / f"{mpid}-{comp}-a2f.pdf")

        elif kind == "gammaband":
            target = in_dir("calc")
            if target is None:
                res.status = "skipped"
                return res
            copy_file(self.root / "config.json", target / "config.json")
            with pushd(target):
                self.run_plot("gammaband", mpid, comp, nkpt, kcut)
                self.run_plot("a2f", mpid, comp)
                move_file(Path(f"{comp}-gamma.pdf"), plots / f"{mpid}-{comp}-gamma.pdf")
                move_file(Path(f"{comp}-a2f.pdf"), plots / f"{mpid}-{comp}-a2f.pdf")

        elif kind == "pdos":
            target = in_dir("dos")
            if target is None:
                res.status = "skipped"
                return res
            with pushd(target):
                copy_file(Path("scf-dos.in"), Path("scf.in"))
                self.run_plot("pdos", mpid, comp, kcut)
                copy_file(Path(f"pdos-{comp}.pdf"), plots / f"plot-pdos-{mpid}-{comp}.pdf")
                copy_file(Path("pdos-spin-resolved.pdf"),
                          plots / f"plot-pdos-{mpid}-{comp}-spin-resolved.pdf")

        elif kind == "wann_band":
            target = in_dir("epw")
            if target is None:
                res.status = "skipped"
                return res
            with pushd(target):
                self.run_plot("wann_band", mpid, comp, kcut)
                move_file(Path("plot.pdf"), plots / f"plot-band-wann-{mpid}-{comp}.pdf")

        elif kind == "bandproj":
            target = in_dir("bands")
            if target is None:
                res.status = "skipped"
                return res
            with pushd(target):
                run_helper("plot_bandproj", "main", comp)
                for pdf in sorted(Path(".").glob("*.pdf")):
                    copy_file(pdf, plots / pdf.name)

        elif kind == "phonproj":
            # bash put `continue` *before* `cd ../../` here, so the shell stayed
            # inside phonon/ for every later material.  pushd makes that a non-issue.
            for sub in ("phonon", "calc"):
                target = in_dir(sub)
                if target is None:
                    continue
                with pushd(target):
                    run_helper("projection_phband", "main", comp, nkpt)
                    self.run_plot("phonproj", mpid, comp, nkpt)
                    move_file(Path(f"plot-proj-{mpid}-{comp}.pdf"), plots)
                break
        else:
            res.status = "skipped"
            res.message = f"unknown plot type {kind!r}"
        return res

    def clean_scan(self, start: int, end: int, track: str | None = None,
                   *_) -> list[Result]:
        """``mainprogram 20`` -- drop heavy files, copy finished runs to ``completed/``."""
        ensure_dir(self.root / "completed")
        self.banner("Cleaning heavy files and copying to 'completed' folder")
        results = self.map(self._clean_one, self.materials(start, end, track))
        write_track_file(self.root / "mpid-finished.in",
                         [(r.mpid, r.compound) for r in results if r.rows.get("finished")])
        write_track_file(self.root / "mpid-notfinished.in",
                         [(r.mpid, r.compound) for r in results
                          if r.rows.get("notfinished")])
        LOG.info("all done")
        return results

    def _clean_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        relax, calc = material.sub("relax"), material.sub("calc")
        comp = material.compound
        if not material.dir.is_dir():
            res.status = "skipped"
            return res

        if self.is_vasp:
            if relax.is_dir():
                for name in ("CHG", "CONTCAR", "DOSCAR", "EIGENVAL", "ELFCAR",
                             "IBZKPT", "LOCPOT", "OSZICAR", "OUTCAR", "PCDAT",
                             "PROCAR", "PROOUT", "XDATCAR", "REPORT"):
                    self.remove(relax / name)
            return res

        # --- QE: every rm is anchored at the material directory ---------------
        if relax.is_dir():
            self.remove(relax / f"{comp}.save", recursive=True)
            self.remove(relax / f"{comp}.xml")
            self.remove_glob(os.fspath(relax / "slurm*"))
        if not calc.is_dir():
            return res

        lambda_out = calc / "lambda.out"
        nlines = len(read_lines(lambda_out)) if lambda_out.is_file() else 0
        if nlines > 0:
            res.say(f"lambda.out: {nlines}")
            self.remove_glob(os.fspath(calc / f"{comp}.wfc*"))
            self.remove(calc / f"{comp}.save", recursive=True)
            self.remove(calc / f"{comp}.xml")
            self.remove(calc / "_ph0", recursive=True)
            self.remove_glob(os.fspath(calc / "slurm*"))
            if self.dry_run:
                LOG.info("[dry-run] would copy %s into completed/", material.name)
            else:
                copy_file(material.dir, self.root / "completed" / material.name)
            res.rows["finished"] = True
        else:
            res.say(f"{material.mpid}-{comp}: not finished")
            res.rows["notfinished"] = True
        return res

    def extract_scan(self, start: int, end: int, track: str | None = None,
                     *_) -> list[Result]:
        """``mainprogram 21`` -- Tc / lambda summary into ``result.csv`` + cifs."""
        self.banner("Extracting superconducting critical temperatures, "
                    "stored in result.csv")
        materials = self.materials(start, end, track)
        self.phcheck_scan(start, end, track)
        self.checkfreq_scan(start, end, track)

        ensure_dir(self.root / "cif")
        results = self.map(self._extract_one, materials)

        # one writer, after the loop -> deterministic row order
        with open(self.root / "result.csv", "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["ID", "compound", "degauss", "lambda", "logomega",
                             "Tc", "Phonon_freq", "qbreasym", "fft_check"])
            for res in results:
                if res.rows.get("row"):
                    writer.writerow(res.rows["row"])

        LOG.info("Checking convergence according to the recipe of Nepal et al. (ML paper)")
        for stale in glob.glob(os.fspath(self.root / "fitting_params_Tc_*")):
            remove(stale)
        remove(self.root / "problem_in_fit.in")
        with pushd(self.root):
            try:
                run_helper("fitting_elph_smearing", "main", "result.csv", 13, 0.005, 10)
            except Exception as exc:                        # noqa: BLE001
                LOG.warning("fitting step skipped: %s", exc)
        LOG.info("all done")
        return results

    def _extract_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        mpid, comp = material.mpid, material.compound
        calc, relax = material.sub("calc"), material.sub("relax")

        degauss = QEText.value(read_text(calc / "scf.in"), "degauss", "")
        freq_tag = "negative_freq" if (self.root / "elph_dir" /
                                       f"{mpid}-{comp}-freq.dat").is_file() \
            else "nonegative_freq"
        qb_tag = "qbreaksymerror" if (self.root / "scf_dir" /
                                      f"{mpid}-qbreaksym").is_file() else "noqbreaksymerror"
        fft_tag = "FFTcompatibleerror" if (self.root / "scf_dir" /
                                           f"{mpid}-fft-grid").is_file() else "noFFTerror"

        lambda_out = calc / "lambda.out"
        lam = logom = tc = "NaN"
        if lambda_out.is_file():
            lines = read_lines(lambda_out)
            row = lines[12] if len(lines) > 12 else ""
            lam, logom, tc = (field_of(row, i) for i in (1, 2, 3))
            lam = lam if is_number(lam) else "NaN"
            logom = logom if is_number(logom) else "NaN"
            tc = tc if is_number(tc) else "NaN"
            res.rows["row"] = [mpid, comp, degauss, lam, logom, tc,
                               freq_tag, qb_tag, fft_tag]

        if (relax / "scf.out").is_file():
            # bash overwrote relax/scf.in with the relaxed input here, destroying
            # the provenance of the run; the cif is produced in a scratch dir.
            with self.scratch(material, "extract", local_scf_dir=False) as work:
                if copy_file(material.scf_relaxed, work / "scf.in"):
                    self.run_scftocif()
                    copy_file(work / "relax.cif", self.root / "cif" / f"{mpid}.cif")
        return res

    # ======================================================================== #
    #  download / info / k-mesh / magmom / fermisurface / substitution
    # ======================================================================== #
    def download_input(self, start: int | None = None, end: int | None = None,
                       track: str | None = None, *_) -> list[Result]:
        """``mainprogram download`` -- pull MP structures and build the inputs.

        The bash script ignored its own ``$1 $2 $3`` and re-read the bounds from
        ``input.in``; here the arguments win and ``input.in`` is only the
        fallback, so the method behaves like every other stage.
        """
        start = self.input.start if start is None else int(start)
        end = self.input.end if end is None else int(end)
        self.banner(" INPUTS DOWNLOADING .......... ")
        if self.is_vasp and not (self.root / "vasp.in").is_file():
            LOG.info("vasp.in not present.  Format: '<VASP keyword> <value>' to replace "
                     "a value from the Materials Project INCAR, '<VASP keyword>' alone "
                     "to remove it.  List replacements first, removals last.")
        ensure_dir(self.root / "input_cif")
        # MP web requests are rate limited: keep the download pool small.
        saved, self.workers = self.workers, min(self.workers, 4)
        try:
            results = self.map(self._download_one, self.materials(start, end, track))
        finally:
            self.workers = saved
        LOG.info("all done")
        return results

    def _download_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        relax = material.sub("relax")

        if self.is_vasp:
            poscar_lines = len(read_lines(relax / "POSCAR"))
            lsorbit = grep_count("LSORBIT .TRUE.", self.root / "vasp.in") > 0
            if poscar_lines < 1:
                res.say(f"Downloading: {material.mpid} {material.compound}")
                with pushd(self.root):
                    run_helper("vasp_input", "main", material.mpid, material.compound)
            else:
                res.say(f"{material.name}/relax/POSCAR present")
            if not (relax / "NSW_0_DETECTED").is_file() or lsorbit:
                copy_file(self.root / "vasp.in", relax / "vasp.in")
                copy_file(self.root / "config.json", relax / "config.json")
                remove(relax / "EIGENVAL")
                with pushd(relax):
                    copy_file(Path("INCAR"), Path("INCAR_backup"))
                    self.run_vasp_process("POSCAR")
            return res

        with pushd(self.root):
            self.run_qe_input(material.mpid)
            move_file(self.root / f"{material.mpid}.cif",
                      self.root / "input_cif" / f"{material.mpid}.cif")
        return res

    def info_scan(self, start: int, end: int, track: str | None = None,
                  *_) -> list[Result]:
        """``mainprogram compound`` -- print structural/electronic information."""
        self.banner("Printing info about compounds. Redirect to a file to keep it.")
        results = self.map(self._info_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _info_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        relax = material.sub("relax")
        res.say("*" * 54)
        res.say(f"*               Compound: {material.compound}")
        res.say("*" * 54)

        if self.is_vasp:
            res.say("Structural parameters before relaxation")
            source = relax / ("POSCAR1" if (relax / "POSCAR1").is_file() else "POSCAR")
            if source.is_file():
                run_helper("crystal", "main", os.fspath(source))
            if (relax / "CONTCAR").is_file():
                res.say("Structural parameters after relaxation")
                run_helper("crystal", "main", os.fspath(relax / "CONTCAR"))
            outcar = relax / "OUTCAR"
            if outcar.is_file():
                res.say(f"*  Valence Electrons: "
                        f"{last_match_field('NELECT', outcar, 3)}")
                res.say(f"*  Fermi Energy: "
                        f"{last_match_field('Fermi energy', outcar, 3)} eV")
                res.say("*  K-mesh info: " +
                        " ".join(read_lines(relax / "KPOINTS")[2:4]))
            return res

        if not material.scf_template.is_file():
            with pushd(self.root):
                self.run_qe_input(material.mpid)
        if material.scf_template.is_file():
            res.say(f"Looking at scf-{material.mpid}.in (before relaxation)")
            run_helper("crystal", "main", os.fspath(material.scf_template))
        if material.scf_relaxed.is_file():
            res.say(f"Looking at {material.scf_relaxed.name} (after relaxation)")
            run_helper("crystal", "main", os.fspath(material.scf_relaxed))

        scf_out = material.sub("calc") / "scf.out"
        if scf_out.is_file():
            res.say(f"*  Valence Electrons: {QEText.nelec(scf_out)}")
            res.say(f"*  Fermi Energy: "
                    f"{last_match_field('Fermi', scf_out, 5)} eV")
        res.say(f"*  KEcutoff: {QEText.value(read_text(material.scf_template), 'ecutwfc')} Ry")
        # bash ran `grep ... $comp1` with comp1 unset, so grep read stdin and the
        # command hung on a terminal; the mesh is simply parsed here.
        source = material.scf_relaxed if material.scf_relaxed.is_file() \
            else material.scf_template
        mesh, shift = QEText.kmesh(read_text(source))
        res.say(f"*  K-mesh info: (automatic) {mesh} shift {shift}")
        elph_in = self.root / "elph_dir" / \
            f"elph-{material.mpid}-{material.compound}.in"
        if elph_in.is_file():
            res.say("*  q-mesh info: " + "; ".join(grep("nq", elph_in)))
        res.say("*" * 54)
        return res

    def double_kmesh(self, start: int, end: int, track: str | None = None,
                     *_) -> list[Result]:
        """``mainprogram change_k`` -- rewrite the k-mesh per ``kpoint.in``."""
        self.banner("Modifying K-point mesh")
        results = self.map(self._double_kmesh_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _double_kmesh_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)

        if self.is_vasp:
            kpoints_path = material.sub("relax") / "KPOINTS"
            lines = read_lines(kpoints_path)
            if len(lines) < 4:
                res.status = "skipped"
                res.message = "KPOINTS not found"
                return res
            base = [to_int(t, 1) for t in lines[3].split()[:3]]
            mesh, _ = read_mesh_file(self.root / "kpoint.in", base, default_divisor=1) \
                if (self.root / "kpoint.in").is_file() else ([2 * b for b in base], None)
            copy_file(kpoints_path, material.sub("relax") / "KPOINTS_old")
            lines[3] = " ".join(str(v) for v in mesh)
            write_lines(kpoints_path, lines)
            res.say(f"  k-mesh {base} -> {mesh}")
            return res

        template = material.scf_template
        if not template.is_file():
            res.status = "skipped"
            res.message = f"scf_dir/scf-{material.mpid}.in not present"
            return res
        text = read_text(template)
        base, shift = QEText.kmesh(text)
        if (self.root / "kpoint.in").is_file():
            mesh, new_shift = read_mesh_file(self.root / "kpoint.in", base,
                                             default_divisor=1)
        else:
            res.say("kpoint.in not present. Default 2*k-mesh along each direction")
            mesh, new_shift = [2 * b for b in base], [0, 0, 0]
        # single backup slot in bash; keep a numbered history instead
        backup = self.root / "scf_dir" / f"scf-{material.mpid}-old.in"
        n = 0
        while backup.exists():
            n += 1
            backup = self.root / "scf_dir" / f"scf-{material.mpid}-old{n}.in"
        copy_file(template, backup)
        lines = read_lines(template)
        for i, line in enumerate(lines):
            if line.strip().startswith("K_POINTS") and i + 1 < len(lines):
                lines[i + 1] = " ".join(str(v) for v in list(mesh) + list(new_shift))
                break
        write_lines(template, lines)
        res.say(f"  k-mesh {base} -> {mesh} (backup: {backup.name})")
        return res

    def magmom_extract(self, start: int, end: int, track: str | None = None,
                       *_) -> list[Result]:
        """``mainprogram magmom_extract`` -- dump the MAGMOM block of each OUTCAR."""
        self.banner("Extracting magnetic moments")
        if not self.is_vasp:
            LOG.info("Available only for VASP")
            return []
        ensure_dir(self.root / "magmom_file")
        results = self.map(self._magmom_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _magmom_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        relax = material.sub("relax")
        outcar, incar = relax / "OUTCAR", relax / "INCAR"
        if not outcar.is_file():
            res.status = "skipped"
            return res
        if not (grep_count("ISPIN = 2", incar) and grep_count("LORBIT", incar)):
            res.status = "skipped"
            res.message = "ISPIN = 2 / LORBIT not set"
            return res
        res.say("ISPIN = 2 and LORBIT found -- extracting MAGMOM from OUTCAR")
        nion = to_int(last_match_field("NIONS", outcar, 12), 0)
        lines = read_lines(outcar)
        block: list[str] = []
        for i, line in enumerate(lines):
            if "magnetization (x)" in line:
                block = lines[i:i + nion + 4]
        write_lines(self.root / "magmom_file" /
                    f"magmom-{material.mpid}-{material.compound}.txt", block)
        return res

    def ifermi_scan(self, start: int, end: int, track: str | None = None,
                    *_) -> list[Result]:
        """``mainprogram fermisurface`` -- ``ifermi`` info + a plotting job."""
        self.banner("Submitting Fermisurface calculations "
                    "(the 'ifermi' package must be installed)")
        if not self.is_vasp:
            LOG.info("ifermi plot works only for DFT = vasp")
            return []
        ensure_dir(self.root / "IFERMI")
        results = self.map(self._ifermi_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _ifermi_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        relax = ensure_dir(material.sub("relax"))
        mpid, comp = material.mpid, material.compound

        if not copy_file(self.root / "run-ifermi.sh", relax / "run-ifermi.sh"):
            body = QEText.drop(read_lines(self.root / "run-vasp.sh"), "mpirun")
            body.append("ifermi plot --property velocity -a 90 --interpolation-factor 10 "
                        "--output fermi-surface.jpg --property-colormap bwr --hide-labels")
            write_lines(relax / "run-ifermi.sh", body)
        # bash appended these three lines on *every* invocation
        lines = [ln for ln in read_lines(relax / "run-ifermi.sh")
                 if "../../IFERMI/" not in ln]
        lines += [f"cp fermi-surface.html ../../IFERMI/{mpid}-{comp}-fs.html",
                  f"cp fermi_info.dat ../../IFERMI/{mpid}-{comp}-fs.dat",
                  f"cp fermi-surface.jpg ../../IFERMI/{mpid}-{comp}-fs.jpg"]
        write_lines(relax / "run-ifermi.sh", lines)

        with pushd(relax):
            if not self.dry_run:
                try:
                    with open("fermi_info.dat", "w") as out:
                        subprocess.run(["ifermi", "info", "--property", "velocity"],
                                       stdout=out, stderr=subprocess.STDOUT, check=False)
                except FileNotFoundError:
                    res.say("ifermi not found on PATH -- info step skipped")
            res.job = self.scheduler.submit("run-ifermi.sh", Path.cwd(), mpid, comp,
                                            "ifermi")
        return res

    def sitesub_scan(self, start: int, end: int, track: str | None = None,
                     *_) -> list[Result]:
        """``mainprogram 29`` -- site-substitution inputs (needs ``bsym``)."""
        self.banner("Preparing input files for site substitutions "
                    "(the 'bsym' package must be installed)")
        results = self.map(self._sitesub_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _sitesub_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        with pushd(self.root):
            if self.is_vasp and not material.dir.is_dir():
                run_helper("vasp_input", "main", material.mpid, material.compound)
            run_helper("site_subs", "main", material.mpid, self.dft, material.compound)
        return res

    def atom_scan(self, start: int, end: int, track: str | None = None,
                  *_) -> list[Result]:
        """Isolated-atom QE inputs in a 23x24x25 A box (called manually)."""
        self.banner("Creating QE input for isolated ion")
        results = self.map(self._atom_one, self.materials(start, end, track))
        write_track_file(self.root / "mpid-atom-list.in",
                         [(f"{r.mpid}-atom", r.rows["element"]) for r in results
                          if r.rows.get("element")])
        LOG.info("all done")
        return results

    def _atom_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        template = material.scf_template
        if not template.is_file():
            res.status = "skipped"
            res.message = f"scf-{material.mpid}.in not present"
            return res
        text = read_text(template)
        species = QEText.card(text, "ATOMIC_SPECIES")
        element = field_of(species[-1], 1) if len(species) > 1 else material.compound

        header = QEText.drop_namelist(QEText.header(text), "IONS")
        header = QEText.drop_namelist(header, "CELL")
        header = QEText.drop(header, "nat", "prefix", "calculation")
        header = QEText.insert_after(header, "&SYSTEM", "  nat = 1,")
        header = QEText.insert_after(header, "&CONTROL", f"  prefix = '{element}',")

        body = species + ["K_POINTS automatic", "1 1 1 0 0 0",
                          "CELL_PARAMETERS angstrom", "23 0 0", "0 24 0", "0 0 25",
                          "ATOMIC_POSITIONS crystal", f"{element} 0 0 0"]
        write_lines(self.root / "scf_dir" /
                    f"scf-{material.mpid}-atom-{element}.in", header + body)
        res.rows["element"] = element
        return res

    def charge_input(self, start: int, end: int, track: str | None = None,
                     *_) -> list[Result]:
        """``mainprogram charge-input`` -- inputs for charged cells."""
        self.banner("Generating input files for systems with non-zero net charge")
        charge_file = self.root / "charge.in"
        if not charge_file.is_file():
            LOG.info("'charge.in' not found -- creating one extra electron and one hole")
            write_lines(charge_file, ["v1 1", "v2 -1"])
        self.args = {"charges": read_indexed_file(charge_file)}
        results = self.map(self._charge_one, self.materials(start, end, track))

        rows = [entry for res in results for entry in res.rows.get("entries", [])]
        old = self.root / "mpid-charge.in"
        if old.is_file():
            move_file(old, self.root / "mpid-charge-2.in")
        write_track_file(old, rows)
        LOG.info("*" * 98)
        LOG.info("Update 'input.in' with mpid-charge.in and run 'mainprogram 1'")
        return results

    def _charge_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        entries: list[tuple[str, str]] = []
        mpid, comp = material.mpid, material.compound
        relax = material.sub("relax")
        nelect = to_int(last_match_field("NELECT", relax / "OUTCAR", 3), 0)

        for _, delta in self.args["charges"]:
            if not self.is_vasp:
                source = material.scf_relaxed if material.scf_relaxed.is_file() \
                    else material.scf_template
                res.say("Utilizing relaxed structure" if material.scf_relaxed.is_file()
                        else "Utilizing unrelaxed structure")
                lines = QEText.insert_after(read_lines(source), "&SYSTEM",
                                            f"  tot_charge = {delta},")
                write_lines(self.root / "scf_dir" / f"scf-{mpid}-{delta}.in", lines)
                entries.append((f"{mpid}-{delta}", comp))
            else:
                total = nelect + to_float(delta, 0.0)
                target = ensure_dir(self.root / f"R{mpid}-{total}-{comp}" / "relax")
                for name in ("KPOINTS", "POTCAR", "run.sh", "POSCAR"):
                    copy_file(relax / name, target / name)
                write_lines(target / "INCAR",
                            [f"NELECT = {total}"] + read_lines(relax / "INCAR"))
                entries.append((f"{mpid}-{total}", comp))
        res.rows["entries"] = entries
        return res

    # ======================================================================== #
    #  generic inner-loop parallelism (modes, pressures, supercells)
    # ======================================================================== #
    def map_any(self, method: Callable[[Any], Any], items: Sequence[Any],
                parallel: bool = True) -> list[Any]:
        """Like :meth:`map` but for the inner loops over modes/pressures/cells."""
        if not items:
            return []
        workers = min(self.workers, len(items)) if parallel else 1
        if workers <= 1:
            return [method(item) for item in items]
        ctx = mp.get_context("fork" if hasattr(os, "fork") else "spawn")
        with ctx.Pool(processes=workers) as pool:
            return pool.map(_pool_entry_any,
                            [(self, method.__func__.__name__, item) for item in items])

    # ======================================================================== #
    #  pressure / volume series
    # ======================================================================== #
    def pressure_input(self, start: int, end: int, track: str | None = None,
                       *_) -> list[Result]:
        """``mainprogram pressure-input`` -- inputs for a pressure/volume series."""
        self.banner("Generating input files for different pressure")
        pressure_file = self.root / "pressure.in"
        if not pressure_file.is_file():
            if self.is_vasp:
                LOG.info("'pressure.in' not found -- isotropic volume scaling of 0.94")
                write_lines(pressure_file, ["v1 0.94"])
            else:
                LOG.info("'pressure.in' not found -- 10 GPa with cell_dofree = 'all'")
                write_lines(pressure_file, ["all", "v1 100"])
        self.args = {"points": read_indexed_file(pressure_file),
                     "relax_type": (read_lines(pressure_file) or ["all"])[0].split()[0]
                     if not self.is_vasp else "all"}
        results = self.map(self._pressure_input_one, self.materials(start, end, track))

        flat, per_material = [], []
        for res in results:
            rows = res.rows.get("entries", [])
            per_material += [f"v{i} {mpid} {comp}" for i, (mpid, comp)
                             in enumerate(rows, 1)]
            flat += rows
        # kept for backwards compatibility: per-material numbering restarts at v1
        write_lines(self.root / "mpid-pressure-1.in", per_material)
        write_track_file(self.root / "mpid-pressure-2.in", flat)
        LOG.info("*" * 98)
        LOG.info("For energy-volume calculations use 'mainprogram 26' "
                 "(reads mpid-pressure-1.in); for anything else update input.in "
                 "with mpid-pressure-2.in and run 'mainprogram 1'.")
        return results

    def _pressure_input_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        mpid, comp = material.mpid, material.compound
        entries: list[tuple[str, str]] = []
        relax_type = self.args["relax_type"]

        for _, raw in self.args["points"]:
            value = to_float(raw, 0.0)
            scaling = "." in str(raw)          # a fraction scales the cell
            tag = raw

            if self.is_vasp:
                target = ensure_dir(self.root / f"R{mpid}-{tag}-{comp}" / "relax")
                for name in ("KPOINTS", "POTCAR", "run.sh"):
                    copy_file(material.sub("relax") / name, target / name)
                incar = QEText.drop(read_lines(material.sub("relax") / "INCAR"),
                                    "ISIF", "NSW")
                if scaling:
                    write_lines(target / "INCAR",
                                ["ISIF = 4", "NSW = 200"] + incar)
                    poscar = read_lines(material.sub("relax") / "POSCAR")
                    poscar[1] = str(value)
                    write_lines(target / "POSCAR", poscar)
                else:
                    write_lines(target / "INCAR",
                                ["ISIF = 3", "NSW = 200"] + incar + [f"PSTRESS = {value}"])
                    copy_file(material.sub("relax") / "POSCAR", target / "POSCAR")
            else:
                source = material.scf_relaxed if material.scf_relaxed.is_file() \
                    else material.scf_template
                if not material.scf_relaxed.is_file():
                    res.say("Structure is not relaxed -- using the unrelaxed one "
                            "(run processes 1-3 first)")
                if scaling:
                    text = read_text(source)
                    header = QEText.drop(QEText.header(text), "pseudo_dir")
                    header = QEText.insert_after(header, "&CONTROL",
                                                 "pseudo_dir = '../../../../pp/',")
                    header = QEText.insert_after(header, "&CELL",
                                                 f"  cell_dofree = '{relax_type}',")
                    species = QEText.card(text, "ATOMIC_SPECIES")
                    positions = QEText.card(text, "ATOMIC_POSITIONS")
                    kpoints = QEText.card(text, "K_POINTS")
                    cell = QEText.card(text, "CELL_PARAMETERS")[1:]
                    scaled = ["CELL_PARAMETERS angstrom"]
                    for line in cell:
                        if line.split():
                            scaled.append(" ".join(f"{to_float(t, 0.0) * value:.8f}"
                                                   for t in line.split()))
                    write_lines(self.root / "scf_dir" / f"scf-{mpid}-{tag}.in",
                                header + species + positions + kpoints + scaled)
                else:
                    press = f"  press = {value},"
                    if relax_type != "all":
                        press = f"  press = {value}, cell_dofree = '{relax_type}',"
                    write_lines(self.root / "scf_dir" / f"scf-{mpid}-{tag}.in",
                                QEText.insert_after(read_lines(source), "&CELL", press))
            entries.append((f"{mpid}-{tag}", comp))
        res.rows["entries"] = entries
        return res

    def pressure_relax_scan(self, start: int, end: int, track: str | None = None,
                            *_) -> list[Result]:
        """``mainprogram 26`` -- relaxation of every point of the pressure series."""
        self.banner("Performing relaxation of systems under pressure")
        self.args = {"points": read_indexed_file(self.root / "pressure.in")}
        results = self.map(self._pressure_relax_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _pressure_relax_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        mpid, comp = material.mpid, material.compound
        pressure_dir = ensure_dir(material.sub("pressure"))
        points = self.args["points"]

        # the per-material track file is built from the entries themselves, not
        # by `grep "$B" ... | grep "$A"` (which matched mp-1 inside mp-10).
        write_track_file(pressure_dir / "mpid-pressure.in",
                         [(f"{mpid}-{raw}", comp) for _, raw in points])
        input_lines = read_lines(self.root / "input.in")
        while len(input_lines) < 4:
            input_lines.append("")
        input_lines[0] = "1"
        input_lines[1] = str(len(points) + 1)
        input_lines[3] = "mpid-pressure.in"
        write_lines(pressure_dir / "input.in", input_lines)

        marker = pressure_dir / "CALC_ALREADY_STARTED"
        if not marker.is_file():
            marker.touch()
            copy_file(self.root / "scf_dir", pressure_dir / "scf_dir")
        else:
            nested = HTESPWorkflow(pressure_dir, workers=1, dry_run=self.dry_run,
                                   submit_command=self.submit_command)
            nested.further_relax_input(1, len(points) + 1, "mpid-pressure.in")

        for index, (_, raw) in enumerate(points, 1):
            name = f"R{mpid}-{raw}-{comp}"
            if self.is_vasp:
                source = self.root / name
                if source.is_dir() and not (pressure_dir / name).is_dir():
                    move_file(source, pressure_dir / name)
                target = ensure_dir(pressure_dir / name / "relax")
                if (target / "NSW_0_DETECTED").is_file():
                    res.say("Already NSW = 0 found")
                    continue
                copy_file(self.root / "run-vasp.sh", target / "run.sh")
                res.job = self.scheduler.submit("run.sh", target, mpid, comp, str(index))
                continue

            target = ensure_dir(pressure_dir / name / "relax")
            scf_out = target / "scf.out"
            if scf_out.is_file() and grep_count("Error", scf_out):
                # bash deleted the *whole* pressure tree for every point here
                res.say(f"  {name}: scf.out reports an error -- skipped "
                        f"(nothing removed)")
                continue
            if not (target / "scf.in").is_file():
                copy_file(self.root / "scf_dir" / f"scf-{mpid}-{raw}.in",
                          target / "scf.in")
            niter = QEText.iterations(scf_out) if scf_out.is_file() else 100
            if niter and niter < 3:
                res.say(f"  {name}: structure already relaxed")
                continue
            copy_file(self.root / "run-scf.sh", target / "run-scf.sh")
            res.say(f" submitting jobs in {name}/relax")
            res.job = self.scheduler.submit("run-scf.sh", target, mpid, comp, str(index))
        return res

    def pressure_ph_scan(self, start: int, end: int, track: str | None = None,
                         *_) -> list[Result]:
        """``mainprogram 27`` -- el-ph for every point of the pressure series."""
        self.banner("Submitting phonon calculations for different pressure")
        self.args = {"points": read_indexed_file(self.root / "pressure.in")}
        results = self.map(self._pressure_ph_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _pressure_ph_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        mpid, comp = material.mpid, material.compound
        pressure_dir = material.sub("pressure")
        if not pressure_dir.is_dir():
            res.status = "skipped"
            return res
        master = self.root / "elph_dir" / f"ph-{mpid}-{comp}.in"

        for _, raw in self.args["points"]:
            # bash looked in R$kk here but pressure-relax-scan creates
            # R<id>-<P>-<comp>; the two never matched.
            point = pressure_dir / f"R{mpid}-{raw}-{comp}" / "relax"
            if not point.is_dir():
                res.say(f"  {point.parent.name}: not present")
                continue

            scf_in, scf_out = point / "scf.in", point / "scf.out"
            block = QEText.final_coordinates(scf_out)
            if block:
                header = QEText.drop_namelist(
                    QEText.drop_namelist(QEText.header(read_text(scf_in)), "IONS"), "CELL")
                header = QEText.replace(header, "'vc-relax'", "'scf'")
                species = QEText.card(read_text(scf_in), "ATOMIC_SPECIES")
                write_lines(scf_in, header + species + block)

            elph_out = point / "elph.out"
            converged = grep_count("Convergence has been achieved", elph_out) > 0
            unconverged = grep_count("No convergence has been achieved", elph_out) > 0
            elph_in = point / "elph.in"

            if converged:
                res.say(f"  {point.parent.name}: el-ph completed, nothing to do")
                continue
            if elph_in.is_file() and not unconverged:
                res.say(f"  {point.parent.name}: walltime -- recover=.true.")
                write_lines(elph_in, QEText.insert_after(
                    QEText.drop(read_lines(elph_in), "recover"),
                    "&inputph", "  recover=.true.,"))
            elif not elph_in.is_file():
                if not copy_file(master, elph_in):
                    res.say(f"  {master.name} missing -- run 'mainprogram epw1' first")
                    continue
            else:
                res.say(f"  {point.parent.name}: not converged, restarting")
                lines = read_lines(master)
                if any("alpha_mix" in ln for ln in lines):
                    lines = QEText.drop(lines, "alpha_mix")
                    write_lines(master, lines)
                    lines = QEText.insert_after(lines, "&inputph",
                                                "  alpha_mix=0.3, nmix_ph=8,")
                else:
                    lines = QEText.insert_after(lines, "&inputph", "  alpha_mix=0.3,")
                write_lines(elph_in, lines)

            copy_file(self.root / "run-elph.sh", point / "run-elph.sh")
            res.job = self.scheduler.submit("run-elph.sh", point, mpid, comp, "elph")
        return res

    def pressure_reset(self, start: int, end: int, track: str | None = None,
                       *_) -> list[Result]:
        """``mainprogram 28`` -- remove the ``pressure/`` tree of each material."""
        self.banner("Removing pressure folder")
        results = self.map(self._pressure_reset_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _pressure_reset_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        if material.dir.is_dir():
            self.remove(material.sub("pressure"), recursive=True)
        return res

    # ======================================================================== #
    #  24 / 25  soft-mode distortions
    # ======================================================================== #
    def distortion_relax_scan(self, start: int, end: int, track: str | None = None,
                              *_) -> list[Result]:
        """``mainprogram 24`` -- relax the structure distorted along each mode."""
        self.banner("Submitting relaxation of structures corresponding to "
                    "different mode displacement.")
        results = self.map(self._distortion_relax_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _distortion_relax_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        if not material.dir.is_dir():
            res.status = "skipped"
            return res
        scf_in = material.sub("relax") / "scf.in"
        if not scf_in.is_file():
            res.status = "skipped"
            res.message = f"scf.in not present inside {material.name}/relax"
            return res

        text = read_text(scf_in)
        natoms = to_int(QEText.value(text, "nat"), 0)
        nmodes = natoms * 3
        header = QEText.header(text) + QEText.card(text, "ATOMIC_SPECIES") + \
            QEText.card(text, "K_POINTS")
        header = QEText.insert_after(header, "&CONTROL", "  disk_io = 'nowf',",
                                     "  nstep = 300,")
        write_lines(material.dir / "scf-header.in", header)

        folder = material.phonon_folder("calc")
        if not copy_file(material.sub(folder) / "dynmat.axsf",
                         material.dir / "dynmat.axsf"):
            res.status = "skipped"
            res.message = "dynmat.axsf not present"
            return res
        copy_file(self.root / "run-scf.sh", material.dir / "run-scf.sh")

        res.say(f"Submitting distortion calculations for {nmodes} modes")
        self.args = {"material": material.name}
        self.map_any(self._band_distort_one,
                     [(material, mode) for mode in range(1, nmodes + 1)])
        return res

    def band_distort_scan(self, start: int, end: int, material_dir: str | None = None,
                          *_) -> list[Any]:
        """``band-distort-scan <first-mode> <last-mode>`` (inside a material dir).

        Kept as a public entry point because ``distortion.sh`` documents it as a
        stand-alone command.
        """
        base = Path(material_dir or Path.cwd()).resolve()
        material = Material(0, base.name.lstrip("R").split("-")[0],
                            "-".join(base.name.lstrip("R").split("-")[1:]), base.parent)
        return self.map_any(self._band_distort_one,
                            [(material, m) for m in range(int(start), int(end))])

    def _band_distort_one(self, item: tuple[Material, int]) -> Result:
        material, mode = item
        res = Result(mode, material.mpid, material.compound)
        base = material.dir
        mode_dir = ensure_dir(base / f"R{mode}")
        scf_out = mode_dir / "scf.out"
        header = base / "scf-header.in"

        relaxed = QEText.is_relaxed(scf_out) if scf_out.is_file() else False
        niter = QEText.iterations(scf_out) if scf_out.is_file() else 100

        if scf_out.is_file() and relaxed:
            if niter < 3:
                res.say(f"mode {mode}: structure already relaxed")
                res.status = "skipped"
                return res
            res.say(f"mode {mode}: resubmitting from the relaxed coordinates")
            concat(mode_dir / "scf.in", header, QEText.final_coordinates(scf_out))
        elif scf_out.is_file() and not relaxed:
            res.say(f"mode {mode}: time out -- restarting from the last cell")
            natoms = QEText.natoms(scf_out)
            block = QEText.last_cell_block(scf_out, natoms)
            move_file(mode_dir / "scf.in", mode_dir / "scf-initial.in")
            move_file(scf_out, mode_dir / "scf-initial.out")
            concat(mode_dir / "scf.in", header, block)
        else:
            with pushd(base):
                cell = run_helper("qe_axsf2cellpos", "main", "dynmat.axsf", mode, 1.0,
                                  capture=True)
            concat(mode_dir / "scf.in", header, cell.splitlines())

        if niter > 2:
            copy_file(base / "run-scf.sh", mode_dir / "run-scf.sh")
            res.say(f" submitting jobs in mode {mode}")
            res.job = self.scheduler.submit("run-scf.sh", mode_dir, material.mpid,
                                            material.compound, f"mode{mode}",
                                            rename=False)
        return res

    def distortion_energy_scan(self, start: int, end: int, track: str | None = None,
                               *_) -> list[Result]:
        """``mainprogram 25`` -- collect the per-mode total energies."""
        self.banner("Extracting the total energies of structures "
                    "corresponding to different modes")
        results = self.map(self._distortion_energy_one,
                           self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _distortion_energy_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        if not material.dir.is_dir():
            res.status = "skipped"
            return res
        copy_file(self.root / "config.json", material.dir / "config.json")
        scf_in = material.sub("relax") / "scf.in"
        if not scf_in.is_file():
            res.status = "skipped"
            res.message = f"scf.in not present inside {material.name}/relax"
            return res
        natoms = to_int(QEText.value(read_text(scf_in), "nat"), 0)
        rows = self.energy_distort_scan(1, natoms * 3 + 1, material.dir)
        with open(material.dir / "Energy-mode.csv", "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["mode", "a", "b", "c", "alpha", "beta", "gamma",
                             "kmesh", "Energy(eV)", "iteration"])
            writer.writerows([row for row in rows if row])
        return res

    def energy_distort_scan(self, start: int, end: int,
                            material_dir: str | None = None, *_) -> list[list]:
        """``energy-distort-scan <first-mode> <last-mode>`` -> csv rows."""
        base = Path(material_dir or Path.cwd()).resolve()
        self.args = {"base": os.fspath(base)}
        ensure_dir(base / "cif")
        return self.map_any(self._energy_distort_one,
                            list(range(int(start), int(end))))

    def _energy_distort_one(self, mode: int) -> list | None:
        base = Path(self.args["base"])
        mode_dir = base / f"R{mode}"
        scf_out = mode_dir / "scf.out"
        if not mode_dir.is_dir() or not scf_out.is_file():
            return None
        iterations = QEText.iterations(scf_out)
        if not QEText.is_relaxed(scf_out):
            LOG.info("%s: not relaxed, iteration: %s", mode, iterations)
            return None

        natoms = QEText.natoms(scf_out) or 1
        nkpt = last_match_field("number of k points=", scf_out, 5)
        energy = QEText.total_energy(scf_out) * RY_TO_EV / natoms

        work = base / SCRATCH_ROOT / f"mode-{mode}"
        ensure_dir(work)
        try:
            concat(work / "scf.in", base / "scf-header.in",
                   QEText.final_coordinates(scf_out))
            copy_file(base / "config.json", work / "config.json")
            with pushd(work):
                self.run_scftocif()
            cellpar = read_text(work / "cellpar.in").split()
            move_file(work / "relax.cif", base / "cif" / f"{mode}.cif")
        finally:
            if not self.keep_scratch:
                shutil.rmtree(work, ignore_errors=True)

        cellpar += [""] * (6 - len(cellpar))
        LOG.info("%s: relaxed, iteration: %s", mode, iterations)
        return [mode, *cellpar[:6], nkpt, f"{energy:.8f}", iterations]

    # ======================================================================== #
    #  phonopy
    # ======================================================================== #
    def _phonopy_dim(self) -> list[int]:
        """``DIM`` from ``setting.conf`` -- tolerant of ``DIM = 2 2 2``.

        bash used ``grep -oP 'DIM=\\K.*'`` which requires no spaces, so with
        phonopy's own spelling the supercell k-mesh became ``0 0 0``.
        """
        for line in read_lines(self.root / "setting.conf"):
            if line.strip().startswith("#") or "DIM" not in line.upper():
                continue
            values = [to_int(t, 0) for t in line.split("=", 1)[-1].split()]
            if len(values) >= 3 and all(values[:3]):
                return values[:3]
        return [2, 2, 2]

    def phonopy_scan(self, start: int, end: int, track: str | None = None,
                     step: Any = 1, *_) -> list[Result]:
        """``mainprogram e0 / phono1..5 / eos-* / ev-collect / phono*-pressure``."""
        step = str(step)
        self.args = {"step": step, "dim": self._phonopy_dim()}
        titles = {"0": "Extracting total energies and .cif structures",
                  "1": "Submitting scf calculations with different displacements",
                  "2": "Computing force constants",
                  "3": "Computing thermodynamic properties",
                  "4": "Phonon bandstructure calculation",
                  "5": "Phonon symmetry analysis",
                  "vp-ph-qha": "phonopy-qha: temperature-dependent properties",
                  "ev-collect": "Extracting energy-volume data into e-v.dat",
                  "eos-bm": "Birch-Murnaghan EOS fit",
                  "eos-vinet": "Vinet EOS fit"}
        self.banner(titles.get(step, f"Phonopy step {step}"))

        materials = self.materials(start, end, track)
        results = self.map(self._phonopy_one, materials)

        if step == "0":
            # bash truncated econv.csv on *every* step, wiping what e0 collected
            with open(self.root / "econv.csv", "w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["ID", "comp", "NIONS", "energy", "niteration"])
                for res in results:
                    if res.rows.get("econv"):
                        writer.writerow(res.rows["econv"])
        if step in ("eos-bm", "eos-vinet"):
            with open(self.root / "eos-fit.dat", "a") as handle:
                for res in results:
                    handle.write(res.rows.get("eos", ""))
        LOG.info("all done")
        return results

    def _phonopy_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        step = self.args["step"]
        handler = {
            "0": self._phonopy_energies, "1": self._phonopy_displacements,
            "2": self._phonopy_forces, "3": self._phonopy_thermal,
            "4": self._phonopy_bands, "5": self._phonopy_irreps,
            "vp-ph-qha": self._phonopy_qha, "ev-collect": self._phonopy_ev_collect,
            "eos-bm": self._phonopy_eos, "eos-vinet": self._phonopy_eos,
        }.get(step)

        if handler is None and step.endswith("-pressure"):
            # bash called `mainprogram vp-phN`, which is not a valid process
            # name, so phono1-pressure .. phono4-pressure only printed
            # "Bad input".  Both spellings are accepted here: the documented
            # `phonoN-pressure` and the legacy internal `vp-phN-pressure`,
            # whose numbering was offset by one (vp-ph1 == e0, vp-ph2 == phono1).
            base = step[: -len("-pressure")]
            if base.startswith("vp-ph"):
                inner = "phono" + str(max(0, to_int(base[len("vp-ph"):], 1) - 1))
            else:
                inner = base
            pressure_dir = material.sub("pressure")
            if not pressure_dir.is_dir():
                res.status = "skipped"
                return res
            nested = HTESPWorkflow(pressure_dir, workers=self.workers,
                                   dry_run=self.dry_run,
                                   submit_command=self.submit_command)
            points = len(read_lines(pressure_dir / "mpid-pressure.in"))
            nested.phonopy_scan(1, points + 1, "mpid-pressure.in",
                                inner.replace("phono", ""))
            return res

        if handler is None:
            res.status = "skipped"
            res.message = f"unknown phonopy step {step!r}"
            return res
        return handler(material, res)

    def _phonopy(self, *args: Any, cwd: Path, generates_input: bool = False) -> None:
        """Run ``phonopy``; a no-op under ``--dry-run`` unless it writes inputs.

        FIX: ``--dry-run`` suppressed *every* phonopy call, including
        ``phonopy -d``.  But ``-d`` is input generation -- a symmetry analysis
        of the relaxed cell that writes the displaced supercells and
        ``phonopy_disp.yaml`` -- and no part of it is a DFT calculation.
        Suppressing it made ``mainprogram phono1 --dry-run`` print "Number of
        supercells: 0" and produce nothing, which contradicts what
        ``--dry-run`` promises: build every input file, never call the
        scheduler.  The submission that follows is still suppressed, by
        :class:`Scheduler`.

        Every other phonopy call *consumes* results -- ``-f`` reads the forces
        of finished runs, ``-t``/``-p`` read FORCE_SETS -- so those stay
        suppressed: there is nothing for them to read here, and phonopy's own
        error would be less clear than saying so.
        """
        if self.dry_run and not generates_input:
            LOG.info("[dry-run] phonopy %s (in %s)", " ".join(str(a) for a in args), cwd)
            return
        subprocess.run([self._phonopy_command(generates_input),
                        *[str(a) for a in args]],
                       cwd=os.fspath(cwd), check=False)

    @staticmethod
    @functools.lru_cache(maxsize=None)
    def _phonopy_command(generates_input: bool = False) -> str:
        """``phonopy`` or ``phonopy-init``, whichever this version wants.

        FIX: phonopy 4 moved the setup operations out of the ``phonopy``
        command::

            phonopy: error: '--dim' is a setup operation that moved to
            'phonopy-init' in v4.

        so ``mainprogram phono1`` produced "Number of supercells: 0" and no
        ``phonopy_disp.yaml`` against any phonopy >= 4 -- silently, because
        phonopy exits 0 after printing that to stderr.  The displacement
        generation goes to ``phonopy-init`` when that executable exists; every
        other call (``-f``, ``-t``, ``-p``) still belongs to ``phonopy``, in
        v4 as before.  On phonopy 3 and earlier there is no ``phonopy-init``
        and the original command is used, so both versions work.
        """
        if generates_input and shutil.which("phonopy-init"):
            return "phonopy-init"
        return "phonopy"

    def _phonopy_energies(self, material: Material, res: Result) -> Result:
        ensure_dir(self.root / "cif")
        relax = material.sub("relax")
        if not self.is_vasp:
            scf_out = relax / "scf.out"
            natoms = QEText.natoms(scf_out) or 1
            energy = QEText.total_energy(scf_out) * RY_TO_EV / natoms
            res.rows["econv"] = [material.mpid, material.compound, natoms,
                                 f"{energy:.10f}", QEText.iterations(scf_out)]
            with pushd(relax):
                self.run_scftocif()
            copy_file(relax / "relax.cif", self.root / "cif" / f"{material.mpid}.cif")
            return res

        outcar, incar = relax / "OUTCAR", relax / "INCAR"
        converged = grep_count(
            "reached required accuracy - stopping structural energy minimisation",
            outcar) > 0
        nsw_zero = grep_count("NSW = 0", incar) > 0
        nsw_absent = grep_count("NSW", incar) < 1
        if not (converged or nsw_zero or nsw_absent):
            res.say(f"not converged: {material.mpid} {material.compound}")
            res.status = "skipped"
            return res
        energy = to_float(last_match_field("TOTEN", outcar, 5), 0.0)
        natoms = to_int(last_match_field("NIONS", outcar, 12), 1) or 1
        res.rows["econv"] = [material.mpid, material.compound, natoms,
                             f"{energy / natoms:.6f}",
                             grep_count("y  w", outcar)]
        with pushd(relax):
            self.run_scftocif("POSCAR")
        copy_file(relax / "relax.cif", self.root / "cif" / f"{material.mpid}.cif")
        return res

    def _phonopy_displacements(self, material: Material, res: Result) -> Result:
        if not material.dir.is_dir():
            res.status = "skipped"
            return res
        phonopy_dir = ensure_dir(material.sub("phonopy"))
        setting = self.root / "setting.conf"
        copy_file(setting, phonopy_dir / "setting.conf")
        dim = self.args["dim"]

        if not self.is_vasp:
            if not copy_file(self.root / "run-scf.sh", phonopy_dir / "run-scf.sh"):
                copy_file(material.sub("relax") / "run-scf.sh",
                          phonopy_dir / "run-scf.sh")
            copy_file(material.sub("relax") / "scf.in", phonopy_dir / "scf.in")
            with pushd(phonopy_dir):
                if setting.is_file():
                    self._phonopy("--qe", "-d", "setting.conf", "-c", "scf.in",
                                  cwd=Path.cwd(), generates_input=True)
                else:
                    self._phonopy("--qe", "-d", f"--dim={' '.join(map(str, dim))}",
                                  "-c", "scf.in", cwd=Path.cwd(),
                                  generates_input=True)
                cells = sorted(Path(".").glob("supercell-*.in"))
                base = read_text("scf.in")
                supercell = read_text("supercell.in")
                natoms = len(QEText.card(supercell, "ATOMIC_POSITIONS")) - 1

                header = QEText.drop_namelist(
                    QEText.drop_namelist(QEText.header(base), "IONS"), "CELL")
                header = QEText.replace(header, "'vc-relax'", "'scf'")
                header = QEText.drop(header, "nat", "pseudo_dir")
                header = QEText.insert_after(header, "&CONTROL",
                                             "pseudo_dir = '../../../pp/',")
                header = QEText.insert_after(header, "&SYSTEM", f"  nat = {natoms},")
                species = QEText.card(base, "ATOMIC_SPECIES")
                mesh, _ = QEText.kmesh(base)
                sub_mesh = [max(1, m // d) for m, d in zip(mesh, dim)]
                kpoints = ["K_POINTS automatic",
                           " ".join(str(v) for v in sub_mesh) + " 0 0 0"]
                cell = [ln for ln in QEText.card(supercell, "CELL_PARAMETERS")[1:]
                        if ln.split()]
                scaled = ["CELL_PARAMETERS angstrom"] + [
                    " ".join(f"{to_float(t, 0.0) * BOHR_TO_ANG:.8f}"
                             for t in line.split()) for line in cell]

                res.say(f"Number of supercells: {len(cells)}")
                for number, cell_file in enumerate(cells, 1):
                    target = ensure_dir(Path(f"R{number}"))
                    positions = QEText.card(read_text(cell_file), "ATOMIC_POSITIONS")
                    write_lines(target / "scf.in",
                                header + species + positions + kpoints + scaled)
                    copy_file("run-scf.sh", target / "run-scf.sh")
                    res.job = self.scheduler.submit("run-scf.sh", target,
                                                    material.mpid, material.compound,
                                                    f"cell{number}", rename=False)
            return res

        # ---- VASP -------------------------------------------------------------
        if not copy_file(self.root / "run-vasp.sh", phonopy_dir / "run.sh"):
            copy_file(material.sub("relax") / "run.sh", phonopy_dir / "run.sh")
        if not (self.root / "vasp-phonopy.in").is_file():
            self.vasp_phonopy_template()
        copy_file(self.root / "vasp-phonopy.in", phonopy_dir / "vasp.in")
        for name in ("INCAR", "POSCAR", "POTCAR", "KPOINTS"):
            copy_file(material.sub("relax") / name, phonopy_dir / name)
        copy_file(self.root / "config.json", phonopy_dir / "config.json")
        with pushd(phonopy_dir):
            self.run_vasp_process("symmetrize")
            if setting.is_file():
                self._phonopy("-d", "setting.conf", cwd=Path.cwd(),
                              generates_input=True)
            else:
                self._phonopy("-d", f"--dim={' '.join(map(str, dim))}",
                              cwd=Path.cwd(), generates_input=True)
            cells = sorted(Path(".").glob("POSCAR-[0-9]*"))
            res.say(f"Number of supercells: {len(cells)}")
            for number, cell_file in enumerate(cells, 1):
                target = ensure_dir(Path(f"R{number}"))
                for name in ("INCAR", "KPOINTS", "POTCAR", "vasp.in", "config.json"):
                    copy_file(name, target / name)
                copy_file(cell_file, target / "POSCAR")
                copy_file("run.sh", target / f"run-{number}.sh")
                with pushd(target):
                    self.run_vasp_process("POSCAR")
                res.job = self.scheduler.submit(f"run-{number}.sh", target,
                                                material.mpid, material.compound,
                                                f"cell{number}", rename=False)
        return res

    def _phonopy_forces(self, material: Material, res: Result) -> Result:
        phonopy_dir = material.sub("phonopy")
        if not phonopy_dir.is_dir():
            res.status = "skipped"
            return res
        pattern = "R*/scf.out" if not self.is_vasp else "R*/vasprun.xml"
        files = sorted(os.fspath(p) for p in phonopy_dir.glob(pattern))
        res.say(f"Number of force calculations: {len(files)}")
        self._phonopy("-f", *[Path(f).relative_to(phonopy_dir) for f in files],
                      cwd=phonopy_dir)
        return res

    def _atom_names(self, material: Material) -> str:
        phonopy_dir = material.sub("phonopy")
        if not self.is_vasp:
            species = QEText.card(read_text(phonopy_dir / "scf.in"), "ATOMIC_SPECIES")
            return " ".join(field_of(ln, 1) for ln in species[1:] if ln.split())
        lines = read_lines(phonopy_dir / "POSCAR")
        return lines[5].strip() if len(lines) > 5 else ""

    def _phonopy_thermal(self, material: Material, res: Result) -> Result:
        phonopy_dir = material.sub("phonopy")
        if not phonopy_dir.is_dir():
            res.status = "skipped"
            return res
        dim = self.args["dim"]
        write_lines(phonopy_dir / "mesh.conf",
                    [f"ATOM_NAME = {self._atom_names(material)}",
                     f"DIM = {' '.join(map(str, dim))}",
                     "MP = 48 48 48", "TPROP = .TRUE.", "TMAX = 2100"])
        self._phonopy("-t", "-p", "-s", "mesh.conf", cwd=phonopy_dir)
        return res

    def _phonopy_bands(self, material: Material, res: Result) -> Result:
        phonopy_dir = material.sub("phonopy")
        if not phonopy_dir.is_dir():
            res.status = "skipped"
            return res
        with pushd(phonopy_dir):
            self.run_vasp_process("scf.in" if not self.is_vasp else "POSCAR")
            dim = self.args["dim"]
            head = [f"ATOM_NAME = {self._atom_names(material)}",
                    f"DIM = {' '.join(map(str, dim))}",
                    "FC_SYMMETRY = .TRUE.", "PRIMITIVE_AXES = AUTO",
                    "EIGENVECTORS=.TRUE.", "BAND_POINTS = 50"]
            write_lines(Path("band.conf"), head + read_lines("band_phonopy.in"))
            self._phonopy("-p", "-s", "band.conf", cwd=Path.cwd())
        return res

    def _phonopy_irreps(self, material: Material, res: Result) -> Result:
        phonopy_dir = material.sub("phonopy")
        if not phonopy_dir.is_dir():
            res.status = "skipped"
            return res
        high_symm = self.root / "elph_dir" / \
            f"high_symm-{material.mpid}-{material.compound}.in"
        copy_file(high_symm, phonopy_dir / "high_symm.in")
        setting = (phonopy_dir / "setting.conf").is_file()
        dim = self.args["dim"]
        with pushd(phonopy_dir):
            analysis = []
            for number, line in enumerate(read_lines("high_symm.in"), 1):
                if not line.strip():
                    continue
                res.say(line)
                args = ["setting.conf"] if setting else \
                    [f"--dim={' '.join(map(str, dim))}"]
                if self.dry_run:
                    continue
                proc = subprocess.run(["phonopy", *args, f"--irreps={line}"],
                                      capture_output=True, text=True, check=False)
                analysis.append(proc.stdout)
                copy_file("irreps.yaml", f"irreps-{number}.yaml")
            write_lines(Path("symmetry_analysis.in"), analysis)
        return res

    def _phonopy_qha(self, material: Material, res: Result) -> Result:
        pressure_dir = material.sub("pressure")
        if not pressure_dir.is_dir():
            res.status = "skipped"
            return res
        if self.dry_run:
            res.say("[dry-run] phonopy-qha --tmax 2000 e-v.dat R*/phonopy/...")
            return res
        yaml_files = sorted(os.fspath(p.relative_to(pressure_dir))
                            for p in pressure_dir.glob("R*/phonopy/thermal_properties.yaml"))
        with open(pressure_dir / "thermo_qha.dat", "a") as handle:
            subprocess.run(["phonopy-qha", "--tmax", "2000", "e-v.dat", *yaml_files],
                           cwd=os.fspath(pressure_dir), stdout=handle,
                           stderr=subprocess.STDOUT, check=False)
        return res

    def _phonopy_ev_collect(self, material: Material, res: Result) -> Result:
        pressure_dir = material.sub("pressure")
        if not pressure_dir.is_dir():
            res.status = "skipped"
            return res
        points = read_track_file(pressure_dir / "mpid-pressure.in", 1, 10_000,
                                 pressure_dir)
        ev_rows, csv_rows = [], []
        for point in points:
            relax = point.dir / "relax"
            if not self.is_vasp:
                scf_out = relax / "scf.out"
                energy = QEText.total_energy(scf_out) * RY_TO_EV
                volume = QEText.volume(scf_out) * AU_TO_ANG3
                press = QEText.pressure(scf_out)
                pv = volume * press * PV_UNIT_CONV
                csv_rows.append([f"{press * 0.1:.3f}", f"{volume:.10f}", f"{pv:.10f}",
                                 f"{energy:.10f}"])
            else:
                outcar = relax / "OUTCAR"
                energy = to_float(last_match_field("free  energy   TOTEN  =", outcar, 5))
                volume = to_float(last_match_field("  volume of cell :", outcar, 5))
                csv_rows.append(["", f"{volume}", "", f"{energy}"])
            if energy == energy and volume == volume:        # not NaN
                ev_rows.append(f"{volume} {energy}")
        write_lines(pressure_dir / "e-v.dat", ev_rows)
        with open(pressure_dir / "e_p_v.csv", "w", newline="") as handle:
            csv.writer(handle).writerows(csv_rows)
        res.say(f"  {len(ev_rows)} converged points written to e-v.dat")
        return res

    def _phonopy_eos(self, material: Material, res: Result) -> Result:
        pressure_dir = material.sub("pressure")
        if not pressure_dir.is_dir():
            res.status = "skipped"
            return res
        eos = "birch_murnaghan" if self.args["step"] == "eos-bm" else "vinet"
        rows = [ln for ln in read_lines(pressure_dir / "e-v.dat") if len(ln.split()) == 2]
        removed = len(read_lines(pressure_dir / "e-v.dat")) - len(rows)
        write_lines(pressure_dir / "e-v.dat", rows)
        res.say(f"{removed} calculations not converged. removing them")
        output = ""
        if not self.dry_run:
            proc = subprocess.run(["phonopy-qha", f"--eos={eos}", "-b", "e-v.dat"],
                                  cwd=os.fspath(pressure_dir), capture_output=True,
                                  text=True, check=False)
            output = proc.stdout
        res.rows["eos"] = (f"Materials id: {material.mpid}, "
                           f"Compound: {material.compound}\n{output}"
                           "------------------------------------------\n")
        return res

    # ======================================================================== #
    #  EPW / WANNIER90 / WannierTools
    # ======================================================================== #
    def epw_bash_scripts(self, start: int, end: int, track: str | None = None,
                         process: str = "epw1", projection: str | None = None,
                         *_) -> list[Result]:
        """``mainprogram epw1..5 / qe-ph / wann-* / epw-*`` -- the EPW pipeline."""
        self.banner("EPW and WANNIER90 calculations", f"step: {process}")
        self.args = {"process": process, "projection": projection,
                     "nkpt": self.input.nkpt, "kcut": self.input.kcut}
        results = self.map(self._epw_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _epw_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        process = self.args["process"]

        if not material.dir.is_dir():
            res.status = "skipped"
            res.message = f"{material.name} folder doesn't exist"
            return res
        if not self.is_vasp and not QEText.final_coordinates(
                material.sub("relax") / "scf.out"):
            res.status = "skipped"
            res.message = "Perform relaxation using processes 1 to 4 first"
            return res

        handler = {"epw1": self._epw1, "epw2": self._epw2, "epw3": self._epw3,
                   "cpcharge": self._epw_cpcharge, "ciftoxsf": self._epw_ciftoxsf,
                   "epw4": self._epw4, "proj": self._epw_proj,
                   "band_wann": self._epw_band_wann,
                   "band_wann2": self._epw_band_wann2, "epw": self._epw_final
                   }.get(process)
        if handler is None:
            res.status = "skipped"
            res.message = f"unknown EPW step {process!r}"
            return res
        return handler(material, res)

    def _epw1(self, material: Material, res: Result) -> Result:
        """Build the scf/nscf/ph inputs for Wannier90, EPW and WannierTools."""
        mpid, comp = material.mpid, material.compound
        if self.is_vasp:
            epw = ensure_dir(material.sub("epw"))
            relax = material.sub("relax")
            nelec = to_int(last_match_field("NELECT", relax / "OUTCAR", 3), 0)
            # bash grepped INCAR before it was copied into epw/
            soc = grep_count("ISPIN = 2", relax / "INCAR") > 0
            nbnd = self.nbnd_heuristic(nelec, soc=soc, style="vasp")
            res.say(f"Number of bands used: {nbnd}")
            copy_file(relax / "KPOINTS", epw / "KPOINTS")
            if not copy_file(relax / "CONTCAR", epw / "POSCAR"):
                copy_file(relax / "POSCAR", epw / "POSCAR")
            for name in ("POTCAR", "INCAR"):
                copy_file(relax / name, epw / name)
            if not copy_file(relax / "run.sh", epw / "run.sh"):
                copy_file(self.root / "run-vasp.sh", epw / "run.sh")
            remove(epw / "EIGENVAL")
            if not (epw / "vasp.in").is_file():
                res.say("vasp.in not present inside epw/ -- creating one")
                write_lines(epw / "vasp.in",
                            ["LWAVE .True.", "NSW 0", "NELM 200", "ISIF 2",
                             f"NBANDS {nbnd}", "EDIFFG", "LORBIT"])
            with pushd(epw):
                self.run_vasp_process("POSCAR")
                res.say("Submitting vasp scf calculation")
                res.job = self.scheduler.submit("run.sh", Path.cwd(), mpid, comp, None)
            return res

        template = material.scf_template
        if not template.is_file():
            res.status = "skipped"
            res.message = f"scf-{mpid}.in not found inside scf_dir"
            return res

        text = read_text(template)
        prefix = QEText.prefix(text)
        block = QEText.final_coordinates(material.sub("relax") / "scf.out")
        nelec = QEText.nelec(material.sub("relax") / "scf.out")
        soc = "lspinorb = .true." in read_text(material.sub("relax") / "scf.in")
        nbnd = self.nbnd_heuristic(nelec, soc=soc, style="epw")
        nbnd_proj = nelec * 2 + 120 if soc else nelec + 120

        kmesh, kshift = QEText.kmesh(text)
        qmesh = self.qmesh_for(kmesh)
        species = QEText.card(text, "ATOMIC_SPECIES")
        masses = [field_of(ln, 2) for ln in species[1:] if ln.split()]

        base = QEText.drop_namelist(QEText.drop_namelist(QEText.header(text),
                                                         "IONS"), "CELL")
        header_scf = QEText.replace(base, "'vc-relax'", "'scf'")
        header_nscf = QEText.replace(base, "'vc-relax'", "'nscf'")
        header_nscf = QEText.drop(header_nscf, "conv_thr")
        header_band = QEText.insert_after(header_nscf, "&SYSTEM", f"  nbnd={nbnd},")
        header_band = QEText.insert_after(header_band, "&ELECTRONS",
                                          "  conv_thr = 1d-10,")
        header_proj = QEText.insert_after(header_nscf, "&SYSTEM",
                                          f"  nbnd={nbnd_proj},")
        header_proj = QEText.insert_after(header_proj, "&ELECTRONS",
                                          "  conv_thr = 1d-10,")
        header_proj = QEText.replace(header_proj, "smearing = 'gauss'",
                                     "smearing = 'cold'")
        kpoint_card = ["K_POINTS automatic",
                       " ".join(str(v) for v in list(kmesh) + list(kshift))]

        with self.scratch(material, "epw1", copy=[f"scf_dir/scf-{mpid}.in"]) as work:
            sdir = work / "scf_dir"
            write_lines(sdir / f"scf-{mpid}-{comp}.in",
                        header_scf + species + kpoint_card + block)
            write_lines(work / "mass.dat", masses)
            write_lines(work / "qpoint.dat", [" ".join(str(q) for q in qmesh)])
            write_lines(work / "kpoint.dat", [" ".join(str(k) for k in kmesh)])

            # coarse (q) grid, then the dense (k) grid
            write_lines(work / "kmesh.grid", [" ".join(str(q) for q in qmesh)])
            self.run_create_epw(mpid, comp, prefix, "nscf")
            write_lines(sdir / f"epw-{mpid}-{comp}-nscf.in",
                        header_band + species + read_lines(work / "nscf_grid.out") + block)

            write_lines(work / "kmesh.grid", [" ".join(str(k) for k in kmesh)])
            self.run_create_epw(mpid, comp, prefix, "nscf")
            grid = read_lines(work / "nscf_grid.out")
            write_lines(sdir / f"scf-{mpid}-{comp}-nscf.in",
                        header_band + species + grid + block)
            write_lines(sdir / f"scf-{mpid}-{comp}-nscf-proj.in",
                        header_proj + species + grid + block)

            self.run_create_epw(mpid, comp, prefix, "ph")
            move_file(work / f"ph-{mpid}-{comp}.in",
                      self.root / "elph_dir" / f"ph-{mpid}-{comp}.in")
            for name in (f"scf-{mpid}-{comp}.in", f"epw-{mpid}-{comp}-nscf.in",
                         f"scf-{mpid}-{comp}-nscf.in",
                         f"scf-{mpid}-{comp}-nscf-proj.in"):
                copy_file(sdir / name, self.root / "scf_dir" / name)
        return res

    def _epw2(self, material: Material, res: Result) -> Result:
        """scf + phonon run that feeds EPW."""
        mpid, comp = material.mpid, material.compound
        return self.stage_and_submit(
            material, "phonon", "run-ph.sh", "ph",
            files={"scf.in": f"scf_dir/scf-{mpid}-{comp}.in",
                   "elph.in": f"elph_dir/ph-{mpid}-{comp}.in"},
            result=res)

    def _epw3(self, material: Material, res: Result) -> Result:
        res.say("epw step 3: copy <QE>/EPW/bin/pp.py next to this file and run "
                "`pp.py <prefix>` inside the phonon/ folder to populate save/.")
        res.say(f"  prefix = {self.prefix_of(material)}")
        return res

    def _epw_cpcharge(self, material: Material, res: Result) -> Result:
        prefix = self.prefix_of(material).strip("'")
        epw = ensure_dir(material.sub("epw"))
        save = ensure_dir(epw / f"{prefix}.save")
        for name in ("charge-density.dat", "data-file-schema.xml"):
            copy_file(material.sub("phonon") / f"{prefix}.save" / name, save / name)
        res.say(f"Saved charge density into {save}")
        return res

    def _epw_ciftoxsf(self, material: Material, res: Result) -> Result:
        with pushd(self.root):
            self.run_create_epw(material.mpid, material.compound,
                                self.prefix_of(material), "ciftoxsf")
        return res

    def _epw4(self, material: Material, res: Result) -> Result:
        """Projectability run used to determine the SCDM parameters."""
        mpid, comp = material.mpid, material.compound
        prefix = self.prefix_of(material)
        epw = ensure_dir(material.sub("epw"))
        ensure_dir(self.root / "epw_dir")
        copy_file(self.root / "scf_dir" / f"scf-{mpid}-{comp}-nscf-proj.in",
                  epw / "nscf-proj.in")
        write_lines(epw / "scf.in",
                    QEText.replace(read_lines(self.root / "scf_dir" /
                                              f"scf-{mpid}-{comp}.in"),
                                   "smearing = 'gauss'", "smearing = 'cold'"))
        with self.scratch(material, "epw4", local_scf_dir=False) as work:
            self.run_create_epw(mpid, comp, prefix, "projection")
            move_file(work / f"projwfc-{mpid}-{comp}.in", epw / "projwfc.in")
        return self.stage_and_submit(material, "epw", "run-proj.sh", "proj", result=res)

    def _epw_proj(self, material: Material, res: Result) -> Result:
        """Extract the projectability curve and fit the SCDM parameters."""
        mpid, comp = material.mpid, material.compound
        prefix = self.prefix_of(material)
        epw = material.sub("epw")
        projwfc = epw / "projwfc.out"
        if not projwfc.is_file():
            res.status = "failed"
            res.message = "projwfc.out not found -- aborting"
            return res

        energies = [field_of(ln, 5) for ln in grep("==", projwfc)]
        weights = [field_of(ln, 3) for ln in grep("|psi|^2", projwfc)]
        rows = sorted(zip(energies, weights), key=lambda r: to_float(r[0], 0.0))
        write_lines(epw / "p_vs_e.dat", [f"{e}\t{p}" for e, p in rows])
        res.say("Install lmfit to automate the SCDM initialisation "
                "(conda install -c conda-forge lmfit)")

        fermi = last_match_field("Fermi", epw / "scf.out", 5, "0.0")
        scdm_dir = ensure_dir(self.root / "scdm_dir")
        with pushd(epw):
            self.run_create_epw(mpid, comp, prefix, "scdm", fermi)
        for pattern in (f"scdm-{mpid}-{comp}", f"scdm-proj-{mpid}-{comp}.png",
                        "p_vs_e_fit*"):
            for hit in sorted(epw.glob(pattern)) + sorted(self.root.glob(pattern)):
                move_file(hit, scdm_dir / hit.name)
        return res

    def _epw_band_wann(self, material: Material, res: Result) -> Result:
        """Stage the Wannier90 run (VASP: LWANNIER90; QE: pw2wannier90)."""
        mpid, comp = material.mpid, material.compound
        projection = self.args["projection"]
        epw = ensure_dir(material.sub("epw"))

        if self.is_vasp:
            copy_file(self.root / "wannier90.json", epw / "wannier90.json")
            with pushd(epw):
                copy_file(Path("INCAR"), Path("INCAR-scf"))
                incar = QEText.drop(read_lines("INCAR"), "NPAR", "LWANNIER90")
                incar.append("LWANNIER90 = .TRUE.")
                write_lines(Path("INCAR"), incar)
                self.run_create_epw(mpid, comp, "0", "epw_band", projection)
                copy_file(Path("INCAR"), Path("INCAR_backup"))
                concat(Path("INCAR_wannier"), Path("INCAR"), Path("wannier-vasp.in"))
                copy_file(Path("INCAR_wannier"), Path("INCAR"))
                if not Path("run.sh").is_file():
                    copy_file(self.root / "run-wannier.sh", Path("run.sh"))
                res.job = self.scheduler.submit("run.sh", Path.cwd(), mpid, comp, None)
            return res

        prefix = self.prefix_of(material)
        nkpt, kcut = self.args["nkpt"], self.args["kcut"]
        ensure_dir(self.root / "epw_dir")
        with self.scratch(material, "band_wann",
                          copy=[f"scf_dir/scf-{mpid}-{comp}.in",
                                f"scf_dir/scf-{mpid}-{comp}-band.in"]) as work:
            write_lines(work / "band.dat",
                        grep("nbnd", work / "scf_dir" / f"scf-{mpid}-{comp}-band.in"))
            ensure_dir(work / "epw_dir")
            self.run_create_epw(mpid, comp, prefix, "pw2wan", projection)
            self.run_create_epw(mpid, comp, prefix, "kpathwan", nkpt, kcut)
            self.run_create_epw(mpid, comp, prefix, "wankmesh")
            self.run_create_epw(mpid, comp, prefix, "epw_band", projection)
            self.collect_scratch(work, {"epw_dir/*": "epw_dir"})
        copy_file(self.root / "scf_dir" / f"scf-{mpid}-{comp}.in", epw / "scf.in")
        copy_file(self.root / "scf_dir" / f"scf-{mpid}-{comp}-nscf.in", epw / "nscf.in")
        copy_file(self.root / "epw_dir" / f"pw2wan-{mpid}-{comp}.in", epw / "pw2wan.in")
        copy_file(self.root / "epw_dir" / f"ex-{mpid}-{comp}.win", epw / "ex.win")
        copy_file(self.root / "run-wannier_band.sh", epw / "run-wannier_band.sh")
        res.say(f"run sbatch run-wannier_band.sh inside {material.name}/epw/")
        return res

    def _epw_band_wann2(self, material: Material, res: Result) -> Result:
        """Band run on the k-path Wannier90 produced (``*_band.kpt``)."""
        mpid, comp = material.mpid, material.compound
        prefix = self.prefix_of(material)
        epw = ensure_dir(material.sub("epw"))
        with self.scratch(material, "band_wann2",
                          copy=[f"scf_dir/scf-{mpid}.in"]) as work:
            self.run_create_epw(mpid, comp, prefix, "band")
            band_in = work / "scf_dir" / f"{mpid}-{comp}-band.in"
            lines = [ln for ln in read_lines(band_in)
                     if "K_POINTS crystal" not in ln]
            lines = QEText.drop(lines, "CELL_PARAMETERS")
            lines.append("K_POINTS crystal")
            write_lines(band_in, lines)
            copy_file(band_in, self.root / "scf_dir" / band_in.name)
        kpt = sorted(epw.glob("*_band.kpt"))
        concat(epw / "scf-band.in", self.root / "scf_dir" / f"{mpid}-{comp}-band.in",
               *kpt)
        copy_file(self.root / "scf_dir" / f"band-{mpid}-{comp}.in", epw / "band.in")
        copy_file(self.root / "wannier_band2.sh", epw / "wannier_band2.sh")
        res.say(f"run sbatch wannier_band2.sh inside {material.name}/epw/")
        return res

    def _epw_final(self, material: Material, res: Result) -> Result:
        """Write ``epw.in`` and stage the ``EPW/`` folder."""
        mpid, comp = material.mpid, material.compound
        prefix = self.prefix_of(material)
        projection = self.args["projection"]
        ensure_dir(self.root / "epw_dir")

        text = read_text(material.scf_template)
        kmesh, _ = QEText.kmesh(text)
        qmesh = self.qmesh_for(kmesh)
        with self.scratch(material, "epw", copy=[f"scf_dir/scf-{mpid}.in"]) as work:
            write_lines(work / "qpoint.dat", [" ".join(str(q) for q in qmesh)])
            write_lines(work / "kpoint.dat", [" ".join(str(k) for k in kmesh)])
            self.run_create_epw(mpid, comp, prefix, "epw", projection)
            move_file(work / "epw.in", self.root / "epw_dir" / f"epw-{mpid}-{comp}.in")

        epw_dir = ensure_dir(material.sub("EPW"))
        copy_file(self.root / "scf_dir" / f"scf-{mpid}-{comp}-nscf.in",
                  epw_dir / "nscf_epw.in")
        copy_file(self.root / "epw_dir" / f"epw-{mpid}-{comp}.in", epw_dir / "epw.in")
        script = self.scheduler.job_name(mpid, comp, "epw")
        copy_file(self.root / "run-epw.sh", epw_dir / script)
        res.say(f"run sbatch {script} inside {material.name}/EPW/")

        bare = prefix.strip("'")
        save = ensure_dir(epw_dir / f"{bare}.save")
        for name in ("charge-density.dat", "data-file-schema.xml"):
            copy_file(material.sub("phonon") / f"{bare}.save" / name, save / name)
        return res

    def wt_bash_scripts(self, start: int, end: int, track: str | None = None,
                        process: str = "wt1-b", *_) -> list[Result]:
        """``mainprogram wt1 / wt2`` -- WannierTools inputs (bulk / with surface)."""
        self.banner("WANNIERTOOLS CALCULATION INPUTS")
        self.args = {"process": process, "nkpt": self.input.nkpt,
                     "kcut": self.input.kcut}
        ensure_dir(self.root / "WT_dir")
        results = self.map(self._wt_one, self.materials(start, end, track))
        LOG.info("all done")
        return results

    def _wt_one(self, material: Material) -> Result:
        res = Result(material.index, material.mpid, material.compound)
        res.say(material.mpid, material.compound)
        mpid, comp = material.mpid, material.compound
        if not material.dir.is_dir():
            res.status = "skipped"
            res.message = f"{material.name} folder doesn't exist"
            return res
        if not material.scf_template.is_file():
            res.status = "skipped"
            res.message = f"scf-{mpid}.in not found inside scf_dir"
            return res
        # bash let scfcheck carry the previous material's value here
        if not QEText.final_coordinates(material.sub("relax") / "scf.out"):
            res.status = "skipped"
            res.message = "Perform relaxation using processes 1 to 4 first"
            return res

        prefix = QEText.prefix(read_text(material.scf_template))
        surface = self.args["process"] == "wt1-s"
        nkpt, kcut = self.args["nkpt"], self.args["kcut"]
        with self.scratch(material, "wt",
                          copy=[f"scf_dir/scf-{mpid}.in",
                                f"scf_dir/scf-{mpid}-{comp}.in"]) as work:
            if surface:
                res.say("wannier tools step 1: inputs including the surface")
                copy_file(material.sub("epw") / "POSCAR-slab", work / "POSCAR-slab")
            else:
                res.say("wannier tools step 1: bulk inputs")
                self.run_create_wt(mpid, comp, prefix, "initialize")
            self.run_create_wt(mpid, comp, prefix, "kpathwan", nkpt, kcut)
            self.run_create_wt(mpid, comp, prefix, "body", "T" if surface else "F")
            move_file(work / f"wt-{mpid}-{comp}.in",
                      self.root / "WT_dir" / f"wt-{mpid}-{comp}.in")
        return res

    # ======================================================================== #
    #  small utility scripts
    # ======================================================================== #
    def cancel_job(self, start_number: int, num_points: int,
                   cancel_command: str = "scancel", *_) -> None:
        """``cancel_job <first id> <count> <command>`` (count is exclusive now)."""
        start_number, num_points = int(start_number), int(num_points)
        for job in range(start_number, start_number + num_points):
            LOG.info("%s %s", cancel_command, job)
            if not self.dry_run:
                subprocess.run([cancel_command, str(job)], check=False)

    def check_calc(self, queue_command: str = "squeue", account_id: str = "",
                   *_) -> None:
        """``check_calc <queue command> <account>`` -- show the slurm logs of jobs."""
        proc = subprocess.run([queue_command], capture_output=True, text=True,
                              check=False)
        for line in proc.stdout.splitlines():
            if account_id and account_id not in line:
                continue
            job = field_of(line, 1)
            if not job.isdigit():
                continue
            LOG.info(job)
            for hit in sorted(self.root.rglob("slurm*")):
                if job in hit.name:
                    LOG.info("  %s", hit.relative_to(self.root))

    def history(self, count: int = 10, *_) -> None:
        """``history.sh`` -- the last ``mainprogram`` commands of the shell history."""
        LOG.info("-" * 53)
        LOG.info("Showing latest %d mainprogram commands", count)
        LOG.info("-" * 53)
        lines = [ln for ln in read_lines(Path.home() / ".bash_history")
                 if "mainprogram " in ln and "history" not in ln and not ln.startswith("vi")]
        for line in lines[-int(count):]:
            LOG.info(line)

    @staticmethod
    def distortion_help() -> str:
        """``distortion.sh`` -- the soft-mode tutorial text."""
        text = (
            "Tutorial: http://www.fisica.uniud.it/~giannozz/QE-Tutorial/handson_phon.html\n"
            "1. Perform the ground-state calculation (ionic and electronic relaxation).\n"
            "2. Phonon calculation at a single q point, e.g. Gamma.\n"
            "3. Diagonalise the dynamical matrix with dynmat.x -> dynmat.axsf.\n"
            "4. Write a header file 'scf-header.in' without lattice parameters.\n"
            "5. Relax each mode:  band_distort_scan(first_mode, last_mode)"
            "  (total modes = 3 * number of ions).\n"
            "6. Run with 'scf' first to check that symmetry-equivalent modes share "
            "the same total energy.\n"
            "7. Then switch scf-header.in to 'vc-relax' for the full relaxation.\n"
            "8. Collect with energy_distort_scan(first_mode, last_mode): relaxed "
            "structures land in cif/ and energies in Energy-mode.csv.\n"
            "9. Provide a submission script named 'run.sh' in the working directory.\n")
        LOG.info(text)
        return text

    def sumpdos(self, atom: str, orbital: str, cwd: str | None = None) -> None:
        """``sumpdos.sh <atom> <orbital>`` -- sum the projected DOS files."""
        where = Path(cwd or Path.cwd())
        jobs = {f"{atom}-tot.dat": f"*({atom})*",
                f"{orbital}-tot.dat": f"*({orbital})*",
                f"{atom}-{orbital}.dat": f"*({atom})*({orbital})*"}
        for output, pattern in jobs.items():
            matches = sorted(glob.glob(os.fspath(where / pattern)))
            if not matches:
                continue
            with open(where / output, "w") as handle:
                subprocess.run(["sumpdos.x", *matches], stdout=handle,
                               stderr=subprocess.STDOUT, cwd=os.fspath(where),
                               check=False)

    def vasp_phonopy_template(self) -> Path:
        """``vasp-phonopy.sh`` -- write the default ``vasp-phonopy.in``."""
        target = self.root / "vasp-phonopy.in"
        write_lines(target, ["PREC Accurate", "IBRION -1", "ISMEAR 0", "SIGMA 0.05",
                             "IALGO 38", "LREAL Auto", "LWAVE .FALSE.",
                             "LCHARG .FALSE.", "NSW 0", "ISIF", "EDIFFG"])
        return target

    def generate_submission_file(self, which_calc: str, parallel_command: str = "mpirun",
                                 nproc: int = 1, *_) -> list[Path]:
        """``generate_submission_file.sh`` -- build ``run-*.sh`` from ``batch.header``.

        The legacy bash version emitted ``run.sh``, ``q2r.sh``, ``matdyn.sh`` ...
        while every scan script expects ``run-scf.sh``, ``run-q2r.sh``,
        ``run-matdyn.sh`` ...  The names produced here are the ones the rest of
        the workflow actually looks for.
        """
        header = self.root / "batch.header"
        if not header.is_file():
            LOG.error("batch.header not found in %s.  Copy one from "
                      "examples/QE/batch.header or examples/VASP/batch.header.",
                      self.root)
            return []
        par, n = parallel_command, int(nproc)
        cmd = {
            "scf": f"{par} -np {n} pw.x < scf.in > scf.out",
            "band": f"{par} -np {n} pw.x < scf-band.in > scf-band.out",
            "bandp": f"{par} -np {n} bands.x < band.in > band.out",
            "dos": f"{par} -np {n} pw.x < scf-dos.in > scf-dos.out",
            "dosp": f"{par} -np {n} dos.x < dos.in > dos.out",
            "pdos": f"{par} -np {n} projwfc.x < pdos.in > pdos.out",
            "elph": f"{par} -np {n} ph.x < elph.in > elph.out",
            "q2r": f"{par} -np {n} q2r.x < q2r.in > q2r.out",
            "dynmat": f"{par} -np {n} dynmat.x < dynmat.in > dynmat.out",
            "matdyn": f"{par} -np {n} matdyn.x < matdyn.in > matdyn.out",
            "matdyn-dos": f"{par} -np {n} matdyn.x < matdyn-dos.in > matdyn-dos.out",
            "lambda": f"{par} -np {n} lambda.x < lambda.in > lambda.out",
            "ph": f"{par} -np {n} pw.x < scf.in > scf.out\n"
                  f"{par} -np {n} ph.x < elph.in > elph.out",
            "proj": f"{par} -np {n} pw.x < nscf-proj.in > nscf-proj.out\n"
                    f"{par} -np {n} projwfc.x -in projwfc.in > projwfc.out",
            "epw": f"{par} -np {n} pw.x < nscf_epw.in > nscf_epw.out\n"
                   f"{par} -np {n} epw.x -npools {n} -i epw.in > epw.out",
            "vasp": f"{par} -np {n} vasp_std",
            "ifermi": "ifermi plot --property velocity --interpolation-factor 10 "
                      "--property-colormap bwr",
        }
        groups = {"qe-elph": ["scf", "band", "bandp", "dos", "dosp", "pdos", "elph",
                              "q2r", "dynmat", "matdyn", "matdyn-dos", "lambda"],
                  "epw-elph": ["scf", "ph", "proj", "epw"],
                  "vasp": ["vasp", "ifermi"]}
        if which_calc == "help":
            LOG.info("Usage: generate_submission_file(which_calc, parallel_command, "
                     "nproc);  which_calc in %s", list(groups))
            LOG.info("batch.header may contain a line reading 'submission here', "
                     "which is replaced by the command; if it does not, the "
                     "command is appended to the header instead.")
            return []
        if which_calc not in groups:
            LOG.error("unknown target %r (expected one of %s)", which_calc, list(groups))
            return []

        template = read_text(header)
        # FIX: the placeholder was mandatory here while `mainprogram jobscript`
        # (generate_submission.py) appends instead -- and not one shipped
        # batch.header contains it.  The result was a run-*.sh that allocated
        # the job and ran nothing.  Both conventions now work.
        has_placeholder = PLACEHOLDER in template
        if not has_placeholder:
            LOG.info("batch.header has no %r line; appending the command instead",
                     PLACEHOLDER)
        written = []
        for name in groups[which_calc]:
            if has_placeholder:
                body = template.replace(PLACEHOLDER, cmd[name])
            else:
                body = template.rstrip("\n") + "\n\n" + cmd[name] + "\n"
            out = self.root / f"run-{name}.sh"
            Path(out).write_text(body)
            os.chmod(out, 0o755)
            written.append(out)
        LOG.info("wrote %s", ", ".join(p.name for p in written))
        return written

    def elph_finished_not_copied(self, start: int, end: int,
                                 track: str | None = None, *_) -> list[Material]:
        """Finished el-ph runs that are not in ``completed/`` yet."""
        pending = []
        for material in self.materials(start, end, track):
            info = self.elph_status(material)
            if info["state"] == self.ELPH_DONE and \
                    not (self.root / "completed" / material.name).is_dir():
                LOG.info("%s %s", material.mpid, material.compound)
                pending.append(material)
        LOG.info("all done")
        return pending


def _pool_entry(payload: tuple["HTESPWorkflow", str, Material]) -> Result:
    """Module-level trampoline so ``Pool.map`` can pickle the work item."""
    workflow, method_name, material = payload
    return workflow._guard(getattr(workflow, method_name), material)


def _pool_entry_any(payload: tuple["HTESPWorkflow", str, Any]) -> Any:
    """Trampoline for :meth:`HTESPWorkflow.map_any`."""
    workflow, method_name, item = payload
    return getattr(workflow, method_name)(item)


# --------------------------------------------------------------------------- #
#  command line interface
# --------------------------------------------------------------------------- #
#: bash script name -> method name.  ``mainprogram`` calls the scripts by these
#: names, so ``python workflow.py <script-name> <args...>`` is a drop-in
#: replacement for running the script itself.
SCRIPTS: dict[str, str] = {
    "relax-scan": "relax_scan",
    "further-relax-input": "further_relax_input",
    "further-relax-scan": "further_relax_scan",
    "create-inputs": "create_inputs",
    "fine-scan": "fine_scan",
    "coarse-scan": "coarse_scan",
    "ph-scan": "ph_scan",
    "q2r-scan": "q2r_scan",
    "matdyn-scan": "matdyn_scan",
    "matdyn-dos-scan": "matdyn_dos_scan",
    "lambda-scan": "lambda_scan",
    "phonband-scan": "phonband_scan",
    "bandscf-scan": "bandscf_scan",
    "band-scan": "band_scan",
    "bandp-scan": "bandp_scan",
    "dos-scan": "dos_scan",
    "dosp-scan": "dosp_scan",
    "pdos-scan": "pdos_scan",
    "plot-scan": "plot_scan",
    "clean-scan": "clean_scan",
    "extract-scan": "extract_scan",
    "dynmat-scan": "dynmat_scan",
    "distortion-relax-scan": "distortion_relax_scan",
    "band-distort-scan": "band_distort_scan",
    "distortion-energy-scan": "distortion_energy_scan",
    "energy-distort-scan": "energy_distort_scan",
    "pressure-input": "pressure_input",
    "pressure-relax-scan": "pressure_relax_scan",
    "pressure-ph-scan": "pressure_ph_scan",
    "pressure-reset": "pressure_reset",
    "charge-input": "charge_input",
    "phonopy-scan": "phonopy_scan",
    "epw-bash-scripts": "epw_bash_scripts",
    "wt-bash-scripts": "wt_bash_scripts",
    "download-input": "download_input",
    "info-scan": "info_scan",
    "phcheck-scan": "phcheck_scan",
    "checkfreq-scan": "checkfreq_scan",
    "double-kmesh": "double_kmesh",
    "magmom-extract": "magmom_extract",
    "ifermi-scan": "ifermi_scan",
    "sitesub-scan": "sitesub_scan",
    "atom-scan": "atom_scan",
    "cancel_job": "cancel_job",
    "check_calc": "check_calc",
    "history.sh": "history",
    "distortion.sh": "distortion_help",
    "sumpdos.sh": "sumpdos",
    "vasp-phonopy.sh": "vasp_phonopy_template",
    "generate_submission_file.sh": "generate_submission_file",
    "elph_finished_but_not_copied_to_completed_folder": "elph_finished_not_copied",
}

#: ``mainprogram`` process number -> method name (documentation / convenience).
PROCESSES: dict[str, str] = {
    "1": "relax_scan", "2": "further_relax_input", "3": "further_relax_scan",
    "4": "create_inputs", "5": "fine_scan", "6": "coarse_scan", "7": "ph_scan",
    "8": "q2r_scan", "9": "matdyn_scan", "10": "matdyn_dos_scan",
    "11": "lambda_scan", "12": "phonband_scan", "13": "bandscf_scan",
    "14": "band_scan", "15": "bandp_scan", "16": "dos_scan", "17": "dosp_scan",
    "18": "pdos_scan", "19": "plot_scan", "20": "clean_scan", "21": "extract_scan",
    "23": "dynmat_scan", "24": "distortion_relax_scan",
    "25": "distortion_energy_scan", "26": "pressure_relax_scan",
    "27": "pressure_ph_scan", "28": "pressure_reset", "29": "sitesub_scan",
    "checkph": "phcheck_scan", "checkfreq": "checkfreq_scan",
    "change_k": "double_kmesh", "compound": "info_scan", "download": "download_input",
    "magmom_extract": "magmom_extract", "fermisurface": "ifermi_scan",
    "pressure-input": "pressure_input", "charge-input": "charge_input",
}


def _coerce(token: str):
    """CLI arguments arrive as strings; make the numeric ones numbers."""
    try:
        return int(token)
    except ValueError:
        pass
    try:
        return float(token)
    except ValueError:
        return token


def main(argv_in: Sequence[str] | None = None) -> int:
    args = list(argv_in if argv_in is not None else sys.argv[1:])
    flags = {"workers": None, "dry_run": False, "keep_scratch": False,
             "root": ".", "verbose": False}
    positional: list[str] = []
    index = 0
    while index < len(args):
        token = args[index]
        if token in ("-h", "--help"):
            print(__doc__)
            print("known commands:\n  " + "\n  ".join(sorted(SCRIPTS)))
            return 0
        if token == "--workers":
            index += 1
            flags["workers"] = int(args[index])
        elif token.startswith("--workers="):
            flags["workers"] = int(token.split("=", 1)[1])
        elif token in ("--dry-run", "-n"):
            flags["dry_run"] = True
        elif token == "--keep-scratch":
            flags["keep_scratch"] = True
        elif token == "--root":
            index += 1
            flags["root"] = args[index]
        elif token.startswith("--root="):
            flags["root"] = token.split("=", 1)[1]
        elif token in ("--verbose", "-v"):
            flags["verbose"] = True
        else:
            positional.append(token)
        index += 1

    if not positional:
        print(__doc__)
        print("known commands:\n  " + "\n  ".join(sorted(SCRIPTS)))
        return 1

    command, rest = positional[0], [_coerce(t) for t in positional[1:]]
    method_name = SCRIPTS.get(command) or PROCESSES.get(command) or \
        command.replace("-", "_")

    workflow = HTESPWorkflow(root=flags["root"], workers=flags["workers"],
                             dry_run=flags["dry_run"],
                             keep_scratch=flags["keep_scratch"],
                             log_level=logging.DEBUG if flags["verbose"] else logging.INFO)
    method = getattr(workflow, method_name, None)
    if method is None or not callable(method):
        LOG.error("unknown command %r.  Known commands: %s",
                  command, ", ".join(sorted(SCRIPTS)))
        return 1
    method(*rest)
    if workflow.failed_count:
        LOG.error("%s: %d material(s) failed\n%s", command,
                  workflow.failed_count, workflow.failure_summary())
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
