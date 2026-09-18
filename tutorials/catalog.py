#!/usr/bin/env python
"""The HTESP tutorial catalogue, expressed as data.

``examples/`` ships 42 worked tutorials -- ``examples/QE/tutorial1..21`` and
``examples/VASP/tutorial1..21``.  Their ``README`` files describe, in prose,
which ``mainprogram`` commands to run in which order, which files have to be
copied in from an earlier tutorial first, and what each step is supposed to
leave behind.  This module turns that prose into records so that
:mod:`tutorials.runner` can execute it and :mod:`tutorials.report` can say
exactly which step of which tutorial stopped.

Nothing here imports :mod:`htesp`; the catalogue is plain data and can be
inspected on a machine with no scientific stack at all.

The QE/VASP numbering offset
----------------------------
The two example trees cover the same topics in the same order up to tutorial
10.  QE tutorial 11 (DFPT electron-phonon coupling and the superconducting
Tc) has no VASP counterpart, so from there on the VASP tutorial *n* is the QE
tutorial *n + 1*; VASP then adds a 21st tutorial (3D Fermi surface with the
IFermi package) that QE does not have.  Rather than keeping two hand-written
lists that can drift apart, both trees are generated from one ordered list of
topics per code -- see :data:`QE_TOPICS`, :data:`VASP_TOPICS` and
:func:`vasp_number_to_qe_number`.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Sequence

#: repository root (``tutorials/`` lives directly below it)
PACKAGE_ROOT = Path(__file__).resolve().parent.parent

def find_examples() -> Path:
    """Locate the read-only example tree the tutorials are copied out of.

    ``examples/`` is 185 MB of reference data and is deliberately not shipped
    inside the wheel, so ``PACKAGE_ROOT / "examples"`` exists only in a source
    checkout or an editable install.  After an ordinary ``pip install .`` the
    package lives in ``site-packages`` and that path points at nothing, which
    is what ``the example tree .../site-packages/examples is missing`` means.

    Looked for, in order: ``$HTESP_EXAMPLES``, ``./examples`` under the current
    directory, then the package-relative path.  ``--examples DIR`` overrides
    all three.
    """
    candidates = []
    from_env = os.environ.get("HTESP_EXAMPLES")
    if from_env:
        candidates.append(Path(from_env).expanduser())
    candidates.append(Path.cwd() / "examples")
    candidates.append(PACKAGE_ROOT / "examples")
    for candidate in candidates:
        if (candidate / "QE").is_dir() or (candidate / "VASP").is_dir():
            return candidate.resolve()
    return PACKAGE_ROOT / "examples"


def searched_for_examples() -> list:
    """The candidate paths, in order, for a message that has to explain itself."""
    out = []
    from_env = os.environ.get("HTESP_EXAMPLES")
    if from_env:
        out.append(f"$HTESP_EXAMPLES -> {Path(from_env).expanduser()}")
    out.append(f"./examples -> {Path.cwd() / 'examples'}")
    out.append(f"beside the package -> {PACKAGE_ROOT / 'examples'}")
    return out


#: the read-only example tree the tutorials are copied out of
EXAMPLES = find_examples()

#: QE tutorial 11 (DFPT el-ph) is the one the VASP tree does not have.
VASP_OFFSET_FROM = 11


# --------------------------------------------------------------------------- #
#  records
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class InputPatch:
    """Edits applied to ``input.in`` before a step runs.

    ``mainprogram`` reads the material range, the tracking file and the plot
    types from ``input.in``; several tutorials (the elastic-constant one most
    obviously) switch tracking file mid-way.  ``None`` means "leave alone".
    """

    start: int | None = None
    end: int | None = None
    track: str | None = None
    plot: str | None = None


@dataclass(frozen=True)
class Step:
    """One ``mainprogram`` invocation inside a tutorial.

    Parameters
    ----------
    id
        Short, unique-within-the-tutorial identifier; also the name of the log
        file and what ``--from`` matches against.
    label
        One line of English for the report.
    command
        A ``mainprogram`` process number (``"1"``, ``"29"``) or command name
        (``"e0"``, ``"pressure-input"``) -- or a callable taking
        ``(workdir, step)`` and returning an exit code, for the few steps that
        are plumbing rather than a command.
    args
        Extra arguments appended after the process token.
    submits
        True when the step hands work to the cluster, so the runner has to wait
        for the jobs before the next step can see their output.
    artifacts
        Glob patterns, relative to the work directory, that must match at least
        one path for the step to count as done.  "Exited 0 and produced
        nothing" is this package's most common silent failure, so a step with
        no artifacts declared is the exception, not the rule.
    job_dirs
        Glob patterns for the stage directories in which the workflow layer
        records job ids (``<stage dir>/.htesp_job.json``).
    timeout
        Seconds before the subprocess is killed.  ``None`` means the runner
        default.
    input_patch
        ``input.in`` edits to apply first.
    needs_api_key
        The step talks to the Materials Project and needs ``MP_API_KEY``.
    needs_dft_output
        The step reads the *output* of a real Quantum ESPRESSO / VASP run.  In
        ``--dry-run`` there is no such output, so the runner records the step as
        skipped with that reason instead of pretending it failed.
    needs_potcar
        The step writes VASP inputs, which are only usable with a POTCAR.  They
        are licensed, so they are never in ``examples/`` and a machine may have
        none at all.  When pymatgen cannot produce one the step is recorded as
        *skipped with that reason* rather than done: the inputs it would write
        are incomplete, and calling that a pass is a false green.
    needs_enumlib
        The step enumerates orderings through pymatgen's ``EnumlibAdaptor``,
        which shells out to ``enum.x``/``multienum.x`` and ``makestr.x`` from
        enumlib.  Those are a separate C/Fortran package, not a Python
        dependency, so preflight warns when they are absent rather than letting
        the step die three minutes in with a RuntimeError.
    after
        Ids of earlier steps *in this tutorial* whose output this step consumes.
        When one of them was skipped there is nothing for this step to read, so
        it is skipped too rather than failing for a reason that is not its own.
    note
        Anything a reader of the report needs to know about this step.
    """

    id: str
    label: str
    command: str | Callable[[Path, "Step"], int]
    args: tuple[str, ...] = ()
    submits: bool = False
    artifacts: tuple[str, ...] = ()
    job_dirs: tuple[str, ...] = ()
    timeout: int | None = None
    input_patch: InputPatch | None = None
    needs_api_key: bool = False
    needs_dft_output: bool = False
    needs_enumlib: bool = False
    needs_potcar: bool = False
    after: tuple[str, ...] = ()
    note: str = ""

    @property
    def is_callable(self) -> bool:
        return not isinstance(self.command, str)

    def describe(self) -> str:
        if self.is_callable:
            return f"<python: {getattr(self.command, '__name__', 'callable')}>"
        return " ".join(("mainprogram", str(self.command), *self.args))


@dataclass(frozen=True)
class Seed:
    """One copy instruction used to populate a tutorial's work directory.

    ``source`` is one of

    ``"code"``
        ``examples/<QE|VASP>/`` -- the shared ``batch.header``, ``config.json``,
        ``input.in`` / ``vasp.in`` and the QE ``pp/`` pseudopotentials.
    ``"self"``
        the tutorial's own example directory.
    ``"archive"``
        members of the tutorial's ``reference*.tar.gz``, flattened into the work
        directory (tutorial 8 keeps its ``.cif`` inputs there).
    anything else
        a tutorial code such as ``"QE/9"``; its *work* directory is the source,
        which is how the relaxation hub feeds everything downstream.
    """

    source: str
    patterns: tuple[str, ...]
    overwrite: bool = True
    link: bool = False
    rename: tuple[tuple[str, str], ...] = ()
    required: bool = False


@dataclass(frozen=True)
class Loop:
    """A block of steps repeated until a convergence probe is satisfied.

    Only the relaxation tutorial needs one: ``mainprogram 2`` / ``3`` / ``e0``
    are repeated until ``econv.csv`` reports ``niteration < 3``.
    """

    steps: tuple[str, ...]
    max_cycles: int = 3
    probe: str = "econv"


@dataclass(frozen=True)
class Tutorial:
    """One tutorial: where it lives, what it needs and what it runs."""

    code: str
    dft: str
    number: int
    topic: str
    title: str
    directory: Path
    depends_on: tuple[str, ...] = ()
    seeds: tuple[Seed, ...] = ()
    steps: tuple[Step, ...] = ()
    loop: Loop | None = None
    stub: bool = False
    note: str = ""

    @property
    def workdir_name(self) -> str:
        return self.code.replace("/", "-")

    def step(self, step_id: str) -> Step | None:
        for step in self.steps:
            if step.id == step_id:
                return step
        return None


# --------------------------------------------------------------------------- #
#  topics
# --------------------------------------------------------------------------- #
#: topic key -> human title (shared by both codes; the DFT name is filled in)
TITLES: dict[str, str] = {
    "jobscript": "Generate submission scripts from batch.header",
    "mp-element": "Materials Project search, element mode",
    "mp-chemsys": "Materials Project search, chemsys mode",
    "oqmd": "OQMD search and input generation",
    "aflow": "AFLOW search and input generation",
    "magnetic": "Input generation in magnetic configuration",
    "combine": "Combine data from several databases",
    "fromcif": "Input generation from .cif files",
    "relax": "Structural relaxation (the hub every later tutorial starts from)",
    "convergence": "Cutoff and k-point convergence tests",
    "elph": "Electron-phonon coupling and superconducting Tc",
    "bands": "Band structure and density of states",
    "pressure": "Input files for different pressures or volumes",
    "substitution": "Input files with site substitutions",
    "elastic": "Elastic constants",
    "hull": "Thermodynamic stability (convex hull)",
    "phonopy": "Phonon band structure with phonopy",
    "eos": "Equation of state",
    "wannier": "Wannier-interpolated band structure",
    "charge": "Input files for a non-zero net charge",
    "magorder": "Enumeration of magnetic orderings",
    "fermisurface": "3D Fermi surface with IFermi",
}

#: the QE tutorial order, tutorial 1 first.
QE_TOPICS: tuple[str, ...] = (
    "jobscript", "mp-element", "mp-chemsys", "oqmd", "aflow", "magnetic",
    "combine", "fromcif", "relax", "convergence", "elph", "bands", "pressure",
    "substitution", "elastic", "hull", "phonopy", "eos", "wannier", "charge",
    "magorder",
)

#: VASP follows the same order with the DFPT el-ph tutorial removed (VASP has
#: no counterpart) and the IFermi Fermi-surface tutorial appended.
VASP_TOPICS: tuple[str, ...] = (
    tuple(t for t in QE_TOPICS if t != "elph") + ("fermisurface",)
)


def vasp_number_to_qe_number(number: int) -> int | None:
    """Return the QE tutorial number covering the same topic as VASP *number*.

    ``None`` for VASP 21, whose topic (IFermi Fermi surfaces) the QE tree does
    not cover.
    """
    if number < VASP_OFFSET_FROM:
        return number
    if number == len(VASP_TOPICS):
        return None
    return number + 1


# --------------------------------------------------------------------------- #
#  seeding
# --------------------------------------------------------------------------- #
#: files every work directory gets from ``examples/<code>/`` before anything
#: the tutorial itself ships is laid on top.
BASE_SEED_FILES = ("batch.header", "config.json", "input.in", "vasp.in")

#: what a downstream tutorial takes from the relaxation hub's work directory.
HUB_OUTPUTS = ("R*-*", "scf_dir", "mpid.in", "econv.csv", "pp")

#: files never copied out of an example directory: they are documentation or
#: the *expected* answer, not input.
SEED_EXCLUDE = ("README", "README.txt", "log", "reference", "reference.tar.gz",
                "reference_*", "reference*", "*.tar.gz", "Y-C_references.tar.gz")


def _base_seeds(dft: str) -> tuple[Seed, ...]:
    seeds = [Seed("code", BASE_SEED_FILES, overwrite=True)]
    if dft == "QE":
        seeds.append(Seed("code", ("pp",), overwrite=False, link=True))
    return tuple(seeds)


def _hub_seed(hub: str) -> Seed:
    return Seed(hub, HUB_OUTPUTS, overwrite=True)


def _self_seed() -> Seed:
    """Everything the tutorial itself ships, minus documentation/references."""
    return Seed("self", ("*",), overwrite=True)


# --------------------------------------------------------------------------- #
#  building the catalogue
# --------------------------------------------------------------------------- #
def _topics(dft: str) -> tuple[str, ...]:
    return QE_TOPICS if dft == "QE" else VASP_TOPICS


def _number_of(topic: str, dft: str) -> int:
    """1-based tutorial number of *topic* in the *dft* tree."""
    return _topics(dft).index(topic) + 1


def code_for(topic: str, dft: str) -> str:
    """``"QE/9"`` for ``("relax", "QE")``; raises when the tree lacks the topic."""
    return f"{dft}/{_number_of(topic, dft)}"


def _resolve_seed(seed: Seed, dft: str) -> Seed:
    """Turn a ``"@topic"`` seed source into a concrete tutorial code."""
    if not seed.source.startswith("@"):
        return seed
    return Seed(code_for(seed.source[1:], dft), seed.patterns, seed.overwrite,
                seed.link, seed.rename, seed.required)


def _build_one(topic: str, dft: str) -> Tutorial:
    from tutorials.steps import spec_for       # imported late: steps imports us

    spec = spec_for(topic, dft)
    number = _number_of(topic, dft)
    code = f"{dft}/{number}"
    directory = EXAMPLES / dft / f"tutorial{number}"
    depends = tuple(code_for(t, dft) for t in spec.depends_topics
                    if t in _topics(dft))
    seeds = (*_base_seeds(dft),
             *(_resolve_seed(s, dft) for s in spec.seeds
               if not s.source.startswith("@") or s.source[1:] in _topics(dft)),
             _self_seed())
    stub = spec.stub or not (directory / "config.json").is_file()
    note = spec.note
    if stub and not spec.stub:
        note = (note + "  " if note else "") + (
            f"examples/{dft}/tutorial{number}/ ships no config.json; the work "
            "directory is seeded from the shared example configuration and from "
            "the tutorials this one depends on.")
    return Tutorial(code=code, dft=dft, number=number, topic=topic,
                    title=f"{TITLES[topic]} ({dft})", directory=directory,
                    depends_on=depends, seeds=seeds, steps=spec.steps,
                    loop=spec.loop, stub=spec.stub, note=note)


def build_catalog() -> dict[str, Tutorial]:
    """Build the whole 42-tutorial catalogue, keyed by code (``"QE/9"``)."""
    out: dict[str, Tutorial] = {}
    for dft in ("QE", "VASP"):
        for topic in _topics(dft):
            tutorial = _build_one(topic, dft)
            out[tutorial.code] = tutorial
    return out


#: the catalogue, built once at import time.
CATALOG: dict[str, Tutorial] = build_catalog()


def use_examples(root) -> dict:
    """Point the catalogue at a different example tree and rebuild it.

    Rebuilding matters: ``_build_one`` reads ``config.json`` from each tutorial
    directory to decide whether it is a stub, so a catalogue built against the
    wrong root marks every tutorial as one.
    """
    global EXAMPLES, CATALOG                         # noqa: PLW0603
    EXAMPLES = Path(root).expanduser().resolve()
    CATALOG = build_catalog()
    return CATALOG


# --------------------------------------------------------------------------- #
#  selection and ordering
# --------------------------------------------------------------------------- #
def normalise_code(text: str) -> str:
    """Accept ``qe/9``, ``QE-9``, ``QE/tutorial9`` and return ``"QE/9"``."""
    token = text.strip().replace("\\", "/").replace("-", "/")
    if "/" not in token:
        raise ValueError(f"{text!r} is not a tutorial code such as 'QE/9'")
    dft, number = token.split("/", 1)
    number = number.lower().removeprefix("tutorial")
    if not number.isdigit():
        raise ValueError(f"{text!r} is not a tutorial code such as 'QE/9'")
    return f"{dft.upper()}/{int(number)}"


def parse_codes(text: str | None) -> list[str]:
    """Split a ``--only`` / ``--skip`` value into normalised codes."""
    if not text:
        return []
    return [normalise_code(part) for part in text.split(",") if part.strip()]


def topological_order(codes: Iterable[str],
                      catalog: dict[str, Tutorial] | None = None) -> list[str]:
    """Order *codes* so that every dependency comes before its dependants.

    Dependencies that are not part of *codes* are left out (the runner reports
    them as blocked only when they were selected and failed).
    """
    catalog = CATALOG if catalog is None else catalog
    wanted = [c for c in codes]
    seen: set[str] = set()
    out: list[str] = []

    def visit(code: str, trail: tuple[str, ...] = ()) -> None:
        if code in seen or code not in wanted:
            return
        if code in trail:
            raise ValueError(f"dependency cycle: {' -> '.join((*trail, code))}")
        for dep in catalog[code].depends_on:
            visit(dep, (*trail, code))
        seen.add(code)
        out.append(code)

    for code in wanted:
        visit(code)
    return out


def select(catalog: dict[str, Tutorial] | None = None, *, code: str = "both",
           only: Sequence[str] = (), skip: Sequence[str] = (),
           include_stubs: bool = True) -> list[str]:
    """Return the ordered list of tutorial codes a run should cover.

    Parameters
    ----------
    code
        ``"QE"``, ``"VASP"`` or ``"both"``.
    only, skip
        Normalised tutorial codes to keep / drop.
    include_stubs
        Keep tutorials the example tree cannot actually run.
    """
    catalog = CATALOG if catalog is None else catalog
    chosen = [c for c, t in catalog.items()
              if code in ("both", t.dft) and (include_stubs or not t.stub)]
    if only:
        chosen = [c for c in chosen if c in set(only)]
    if skip:
        chosen = [c for c in chosen if c not in set(skip)]
    return topological_order(chosen, catalog)


def missing_directories(catalog: dict[str, Tutorial] | None = None) -> list[str]:
    """Tutorial codes whose example directory is not on disk."""
    catalog = CATALOG if catalog is None else catalog
    return [c for c, t in catalog.items() if not t.directory.is_dir()]


def format_catalog(catalog: dict[str, Tutorial] | None = None) -> str:
    """The ``--list`` output: one block per tutorial, steps indented."""
    catalog = CATALOG if catalog is None else catalog
    lines: list[str] = []
    for dft in ("QE", "VASP"):
        codes = [c for c, t in catalog.items() if t.dft == dft]
        lines.append(f"=== {dft} ({len(codes)} tutorials) " + "=" * 40)
        for code in codes:
            tut = catalog[code]
            flag = "  [STUB]" if tut.stub else ""
            lines.append(f"{code:<9} {tut.title}{flag}")
            if tut.depends_on:
                lines.append(f"{'':<9}   after: {', '.join(tut.depends_on)}")
            for number, step in enumerate(tut.steps, 1):
                marks = []
                if step.submits:
                    marks.append("submits")
                if step.needs_api_key:
                    marks.append("MP_API_KEY")
                if step.needs_dft_output:
                    marks.append("needs DFT output")
                suffix = f"   ({', '.join(marks)})" if marks else ""
                lines.append(f"{'':<9}   {number:>2}. {step.id:<18} "
                             f"{step.describe()}{suffix}")
        lines.append("")
    return "\n".join(lines)
