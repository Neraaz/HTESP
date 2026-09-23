#!/usr/bin/env python
"""What a run produced, and what each file is for.

The tutorial runner has always answered "did it work?".  It could not answer
"what did it give me?", and that gap hid three real defects: two steps
declared artefacts they never write (``pressure-input`` claimed a directory
that ``mainprogram 26`` creates later; ``update-input`` claimed a file only
written when a structure is *not* relaxed), and one tutorial's declared
artefact is shipped inside ``examples/`` itself, so the glob matched whether
or not the step ran.

Every step now records the paths that appeared or changed while it ran
(:attr:`tutorials.state.StepState.produced`).  This module turns that into a
table: the file, the step that wrote it, and a line of English saying what it
is for.  A step that produced nothing is shown as such -- that is the signal
those three defects would have given.
"""
from __future__ import annotations

import fnmatch

#: glob -> what that file is, most specific pattern first.
#:
#: The step's own label says what the *step* did; this says what the *file*
#: is, which is the question someone reading a work directory actually has.
DESCRIPTIONS: tuple[tuple[str, str], ...] = (
    ("batch.header", "scheduler directives and module loads, copied into every "
                     "submission script"),
    ("run-*.sh", "submission script: batch.header with one stage's run command "
                 "appended"),
    ("CALC_VISIBLE_WITH_*", "marker choosing how submitted jobs are named"),
    ("log", "mainprogram's own log of this run"),
    ("input.in", "the six-line control file: material range, k-mesh, tracking "
                 "file, plot types, code"),
    ("config.json", "the configuration this tutorial ran with"),
    ("mpid-list.in", "search results: one 'v<N> <id> <compound>' line per "
                     "material found"),
    ("mpid.in", "tracking file naming the materials this campaign works on"),
    ("mpid-pressure-*.in", "tracking lists for the volume series"),
    ("mpid-list-not-relaxed.in", "materials whose relaxation has not converged "
                                 "yet"),
    ("download/*.csv", "raw search results from one database, before filtering"),
    ("scf_dir/scf-relax-*.in", "Quantum ESPRESSO input rebuilt from the relaxed "
                               "structure"),
    ("scf_dir/*.in", "Quantum ESPRESSO input template, one per material"),
    ("*/relax/scf.in", "Quantum ESPRESSO relaxation input"),
    ("*/relax/scf.out", "Quantum ESPRESSO relaxation output"),
    ("*/relax/INCAR*", "VASP calculation parameters"),
    ("*/relax/KPOINTS", "VASP k-point mesh"),
    ("*/relax/POSCAR*", "VASP structure"),
    ("*/relax/CONTCAR", "VASP relaxed structure"),
    ("*/relax/POTCAR", "VASP pseudopotentials, concatenated per element"),
    ("*/relax/OUTCAR*", "VASP output"),
    ("*/phonopy/phonopy_disp.yaml", "phonopy's summary of the displacements it "
                                    "generated"),
    ("*/phonopy/supercell-*.in", "one displaced supercell, as a QE input"),
    ("*/phonopy/POSCAR-*", "one displaced supercell, as a VASP structure"),
    ("*/phonopy/SPOSCAR", "the undisplaced supercell"),
    ("*/phonopy/supercell.in", "the undisplaced supercell, as a QE input"),
    ("*/phonopy/R*/*", "one displacement's own run directory"),
    ("*/phonopy/scf.in", "the relaxed cell phonopy displaces"),
    ("*/pressure/*", "one isotropically scaled cell of the volume series"),
    ("econv.csv", "total energy and iteration count per material"),
    ("elastic.csv", "elastic tensor and the moduli derived from it"),
    ("eos-fit.dat", "equation-of-state fit"),
    ("e-v.dat", "energy against volume, for the equation-of-state fit"),
    ("convergence_result", "cutoff / k-mesh convergence curve"),
    ("cellpar.dat", "lattice parameters of the scaled cells"),
    ("pressure.in", "the pressures or volumes the series was built for"),
    ("kpath/*.dat", "high-symmetry k-path through the Brillouin zone"),
    ("plots/*", "generated figure"),
    ("*.cif", "structure in CIF form"),
    ("*/relax/*", "input for this material's relaxation"),
)


#: files every `mainprogram` invocation touches whatever it does.
#:
#: `log` is appended by the command itself on every run, so counting it as
#: output means no step ever looks like it produced nothing -- and "this step
#: passed without writing anything" is the whole signal this module exists to
#: give.  It is still listed; it just does not count.
AMBIENT = ("log",)


def is_ambient(path: str) -> bool:
    """True for files that say nothing about what a step accomplished."""
    return path in AMBIENT


def real_output(step) -> list[str]:
    """What a step wrote, excluding the files every command touches."""
    return [path for path in step.produced if not is_ambient(path)]


def describe(path: str) -> str:
    """One line of English for a produced file, or ``""`` when unknown.

    Matching is tried three ways, most specific first: the whole path, the
    path with a leading directory implied (so ``*/relax/scf.in`` matches a
    top-level ``relax/scf.in``), and -- only for patterns that name no
    directory -- the bare file name, so ``run-*.sh`` still describes
    ``R<id>/phonopy/R1/run-scf.sh``.  ``fnmatch``'s ``*`` crosses ``/``, so
    the whole-path attempt alone would miss every nested copy.
    """
    name = path.rsplit("/", 1)[-1]
    for pattern, text in DESCRIPTIONS:
        if fnmatch.fnmatch(path, pattern):
            return text
        if fnmatch.fnmatch("*/" + path, pattern):
            return text
        if "/" not in pattern and fnmatch.fnmatch(name, pattern):
            return text
    return ""


def rows(state) -> list[tuple[str, str, str, str]]:
    """``(tutorial, step, path, description)`` for everything a run wrote."""
    out = []
    for code, tutorial in sorted(state.tutorials.items()):
        for step in sorted(tutorial.steps.values(),
                           key=lambda s: (s.index, s.cycle)):
            for path in step.produced:
                out.append((code, step.step_id, path, describe(path)))
    return out


def render(state, limit_per_step: int = 12) -> str:
    """The manifest, as text.

    Long lists are trimmed: a phonopy tutorial writes a supercell per
    displacement and a reader does not need all of them enumerated.
    """
    lines = ["", "WHAT THIS RUN PRODUCED", ""]
    total = 0
    for code, tutorial in sorted(state.tutorials.items()):
        steps = sorted(tutorial.steps.values(), key=lambda s: (s.index, s.cycle))
        wrote = [s for s in steps if real_output(s)]
        if not wrote and tutorial.status in ("skipped", "blocked"):
            continue
        lines.append(f"{code}  {tutorial.title}")
        for step in steps:
            if step.status in ("skipped", "blocked"):
                continue
            if not real_output(step):
                lines.append(f"    {step.step_id:<18s} (wrote nothing"
                             + (" but the log)" if step.produced else ")"))
                continue
            shown = step.produced[:limit_per_step]
            for number, path in enumerate(shown):
                label = step.step_id if number == 0 else ""
                text = describe(path)
                lines.append(f"    {label:<18s} {path}"
                             + (f"\n    {'':<18s}     -- {text}" if text else ""))
            hidden = len(step.produced) - len(shown)
            if hidden > 0:
                lines.append(f"    {'':<18s} ... and {hidden} more")
            total += len(step.produced)
        lines.append("")
    lines.append(f"{total} file(s) written under the work directories.")
    return "\n".join(lines)
