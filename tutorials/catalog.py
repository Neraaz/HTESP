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
    timeout
        Seconds before the subprocess is killed.  ``None`` means the runner
        default.
    input_patch
        ``input.in`` edits to apply first.
    needs_api_key
        The step talks to the Materials Project and needs ``MP_API_KEY``.
    needs_dft_output
        The step reads the *output* of a real Quantum ESPRESSO / VASP run.  In
        This runner never performs one, so the step is recorded as skipped
        with that reason -- and with the tutorial's own instructions for
        running it for real -- instead of being pretended to have failed.
    needs_relax_output
        The step reads a finished *relaxation* and nothing else -- the total
        energy, the relaxed structure, the cell to build its own inputs from.
        It needs no DFT code, only that output, and the tutorial's reference
        has it (:data:`REFERENCE_OUTPUT`), so it runs here.  This is the
        narrow sibling of ``needs_dft_output``: a step needing a *later*
        stage -- FORCE_SETS, ph.x output, computed bands -- keeps the wider
        flag and is still skipped, because nothing here can produce those.
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
    timeout: int | None = None
    input_patch: InputPatch | None = None
    needs_api_key: bool = False
    needs_dft_output: bool = False
    needs_relax_output: bool = False
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
    ``"reference-output"``
        the *relaxation output* from the tutorial's reference, placed back into
        ``R<mpid>-<compound>/relax/``.  See :data:`REFERENCE_OUTPUT`.
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
    #: how many times to run this tutorial before calling it failed.  More
    #: than one only where the failure is likely to be someone else's server
    #: -- see :data:`ATTEMPTS`.
    attempts: int = 1
    #: True when a timeout here is the service's doing, not HTESP's, and
    #: should be reported as skipped -- see :data:`FLAKY_SERVICE_TOPICS`.
    flaky_service: bool = False
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


#: the names a tutorial's written instructions go by
README_NAMES = ("README", "README.txt")


def readme_for(tutorial: "Tutorial",
               catalog: dict | None = None) -> Path | None:
    """The written instructions for running *tutorial* for real.

    Twelve of the forty-two tutorials ship no ``README`` of their own --
    QE/15 and eleven VASP ones.  The two trees cover the same topics in the
    same order, though, so the counterpart's instructions are the right ones
    to read: ``VASP/9`` has none, ``QE/9`` describes the same relaxation.
    The QE/VASP numbering diverges from 11 onward, which
    :func:`vasp_number_to_qe_number` already knows about.

    Returns None only when neither the tutorial nor its counterpart has one
    -- VASP/21 (IFermi), whose topic the QE tree does not cover at all.
    Pointing at a path that does not exist would be worse than saying
    nothing.
    """
    catalog = CATALOG if catalog is None else catalog

    def own(candidate) -> Path | None:
        for name in README_NAMES:
            path = candidate.directory / name
            if path.is_file():
                return path
        return None

    found = own(tutorial)
    if found is not None:
        return found

    if tutorial.dft == "VASP":
        number = vasp_number_to_qe_number(tutorial.number)
        other = catalog.get(f"QE/{number}") if number else None
    else:
        other = next((o for o in catalog.values()
                      if o.dft == "VASP"
                      and vasp_number_to_qe_number(o.number) == tutorial.number),
                     None)
    return own(other) if other is not None else None


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

#: The relaxation output the reference carries, per code -- and *only* that.
#:
#: Seventeen steps across the catalogue do nothing but read a finished
#: relaxation: `e0` wants the total energy, `mainprogram 2` the relaxed
#: structure, `mainprogram 4` and `pressure-input` and `phono1` the relaxed
#: cell to build their own inputs from.  None of them needs a DFT code, only
#: its output, and the references have it -- so with these files in place they
#: run for real instead of being skipped.
#:
#: The whitelist is deliberately narrow.  The reference also contains
#: ``econv.csv`` and ``scf_dir/scf-relax-*.in``, which are the *answers* those
#: steps must produce; seeding them would let a step pass by finding an
#: artefact it never wrote.  The numbered copies (``scf.out1``, ``OUTCAR1``,
#: ``POSCAR1``) are left out for the same reason -- ``mainprogram 2`` makes
#: them, counting the existing ones to pick the next number.  And
#: ``NSW_0_DETECTED`` is a state marker saying the relaxation already
#: finished, which would make the step skip its own work.
REFERENCE_OUTPUT = {
    "QE": ("scf.out",),
    "VASP": ("OUTCAR", "CONTCAR", "OSZICAR"),
}

#: topics worth attempting more than once, and how many attempts in total.
#:
#: OQMD is the least reliable service the tutorials touch: in one session it
#: was unresponsive long enough to fail a sweep outright, in another its
#: `oqmd-download` hung for nineteen minutes with two open sockets (its
#: client, `qmpy_rester`, builds a bare `requests.Session()` with no timeout,
#: so a stalled connection blocks for ever).  A search that normally takes 35
#: seconds is not broken because one call stalled, so the tutorial gets a
#: second run at it before being called a failure.
ATTEMPTS = {"oqmd": 2}

#: topics where running out of time says more about the service than about
#: HTESP, so the run is not marked failed for it.
#:
#: OQMD alone qualifies.  Its searches have finished in 35 seconds and in 100;
#: it has been unresponsive for an entire sweep; and its client has no request
#: timeout, so a stalled connection once ran for nineteen minutes.  None of
#: that is something a reader can act on, and turning a whole sweep red for it
#: hides the failures that *are* actionable -- the same reasoning that already
#: records an absent POTCAR or API key as skipped rather than failed.
#:
#: Only a *timeout* is forgiven.  An OQMD query that answers with the wrong
#: data still fails, because that is a defect somebody can fix.
FLAKY_SERVICE_TOPICS = ("oqmd",)

#: the tutorial whose reference supplies that output, per code
REFERENCE_SOURCE = {"QE": "QE/9", "VASP": "VASP/9"}

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


def _reference_output_seed(dft: str) -> Seed:
    """Put the reference relaxation output back into ``R*/relax/``."""
    return Seed("reference-output", REFERENCE_OUTPUT[dft], overwrite=False)


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
    # The relaxation hub is where the reference output belongs; everything
    # downstream inherits it through HUB_OUTPUTS' copy of R*-*.
    #
    # The hub's *own* steps are not unblocked by it.  QE/9 and VASP/9 exist to
    # run a relaxation; feeding them the answer and calling it a pass would
    # describe the wrong thing.  They stay skipped, pointed at their README
    # like every other step that needs a real run -- while the relaxed
    # structure they would have produced is made available to the eleven
    # downstream steps that only need to *read* one.
    if topic == "relax":
        seeds = (*seeds, _reference_output_seed(dft))
    stub = spec.stub or not (directory / "config.json").is_file()
    note = spec.note
    if stub and not spec.stub:
        note = (note + "  " if note else "") + (
            f"examples/{dft}/tutorial{number}/ ships no config.json; the work "
            "directory is seeded from the shared example configuration and from "
            "the tutorials this one depends on.")
    return Tutorial(code=code, dft=dft, number=number, topic=topic,
                    title=f"{TITLES[topic]} ({dft})", directory=directory,
                    depends_on=depends, seeds=seeds,
                    steps=_with_job_scripts(spec.steps),
                    loop=spec.loop, stub=spec.stub, note=note,
                    attempts=ATTEMPTS.get(topic, 1),

                    flaky_service=topic in FLAKY_SERVICE_TOPICS)


def _with_job_scripts(steps: tuple["Step", ...]) -> tuple["Step", ...]:
    """Put a ``jobscript`` step in front of any tutorial that submits.

    ``HTESPWorkflow.stage_and_submit`` copies ``run-<stage>.sh`` from the
    project root into the stage directory and submits *that*; when the script
    is not there it returns ``status="skipped"`` with "run-scf.sh not found in
    project root" and carries on.  Thirteen tutorials declared submitting steps
    without ever running ``mainprogram jobscript``, so in real mode they
    skipped every submission, exited 0, and left the relaxation -- and
    everything that depends on it -- undone while the run looked healthy.

    The scripts are built from ``batch.header`` and ``job_script.command_list``,
    so they are per-work-directory and cannot simply be seeded once.
    """
    from tutorials.steps import JOBSCRIPT_STEP

    if not any(step.submits for step in steps):
        return steps
    if any(step.command == "jobscript" for step in steps):
        return steps
    return (JOBSCRIPT_STEP, *steps)


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
#: the two example trees, as ``--only`` will accept them
DFT_TREES = ("QE", "VASP")


def normalise_code(text: str) -> str:
    """Accept ``qe/9``, ``QE-9``, ``QE/tutorial9`` and return ``"QE/9"``.

    A bare tree name -- ``QE`` or ``vasp`` -- comes back as ``"QE/*"``, which
    :func:`select` expands to every tutorial in that tree.  It replaces the
    old ``--code`` flag: two ways of narrowing the same selection was one more
    than this needed, and ``--only QE`` reads the same as ``--only QE/9``.
    """
    token = text.strip().replace("\\", "/").replace("-", "/")
    if "/" not in token:
        if token.upper() in DFT_TREES:
            return f"{token.upper()}/*"
        raise ValueError(
            f"{text!r} is not a tutorial code such as 'QE/9', "
            "nor a tree name ('QE' or 'VASP')")
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


def iters(codes: Iterable[str],
          catalog: dict[str, Tutorial] | None = None) -> list[list[str]]:
    """Group *codes* into dependency iters: iteration *n* may start once *n-1* ends.

    A real campaign cannot be one pass.  Iteration 0 is everything that depends on
    nothing -- input generation, the database front ends, and the relaxation
    hubs -- and it ends with jobs sitting in the queue.  Nothing in iteration 1 can
    honestly start until those relaxations have finished and converged, because
    iteration 1 *is* the tutorials that read the relaxed structure.  Running the
    whole catalogue in one pass either blocks for days inside `squeue` polling
    or, worse, proceeds on a structure that is not relaxed yet.

    The grouping is derived from ``depends_on`` alone, so it stays correct as
    the catalogue changes; a tutorial whose dependency was not selected counts
    as a root, matching :func:`topological_order`, which drops unselected
    dependencies rather than pulling them in.

    Returns
    -------
    list[list[str]]
        Iteration 0 first.  Each iteration is in topological order, so it can be run as
        it stands.  An empty selection gives an empty list.
    """
    catalog = CATALOG if catalog is None else catalog
    ordered = topological_order(codes, catalog)
    selected = set(ordered)
    depth: dict[str, int] = {}
    for code in ordered:                      # topological order: deps first
        deps = [d for d in catalog[code].depends_on if d in selected]
        depth[code] = 1 + max((depth[d] for d in deps), default=-1)
    out: list[list[str]] = [[] for _ in range(max(depth.values(), default=-1) + 1)]
    for code in ordered:
        out[depth[code]].append(code)
    return out


def _expand(codes: Sequence[str], catalog: dict[str, Tutorial]) -> set[str]:
    """Turn ``{"QE/*", "VASP/14"}`` into the concrete codes it names."""
    out: set[str] = set()
    for code in codes:
        tree, _, number = code.partition("/")
        if number == "*":
            out.update(c for c, t in catalog.items() if t.dft == tree)
        else:
            out.add(code)
    return out


def select(catalog: dict[str, Tutorial] | None = None, *,
           only: Sequence[str] = (), skip: Sequence[str] = ()) -> list[str]:
    """Return the ordered list of tutorial codes a run should cover.

    Parameters
    ----------
    only, skip
        Normalised tutorial codes to keep / drop.  ``"QE/*"`` -- what
        :func:`normalise_code` returns for a bare ``QE`` -- means the whole
        tree, which is how ``--only QE`` replaced the old ``--code`` flag.

    Stub tutorials (ones the example tree cannot run as shipped) are always
    included: preflight warns about them and the report names them, which is
    more useful than a flag for hiding them.
    """
    catalog = CATALOG if catalog is None else catalog
    chosen = list(catalog)
    if only:
        chosen = [c for c in chosen if c in _expand(only, catalog)]
    if skip:
        chosen = [c for c in chosen if c not in _expand(skip, catalog)]
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
