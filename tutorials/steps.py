#!/usr/bin/env python
"""Per-topic step lists for :mod:`tutorials.catalog`.

One function per tutorial *topic* (not per tutorial number), because the QE and
the VASP tree cover the same topics under different numbers -- see the offset
note in :mod:`tutorials.catalog`.  Each builder takes the DFT code and returns
a :class:`TopicSpec`: the ordered steps, an optional convergence loop, the
topics this one has to run after, extra seeding instructions, and whether the
topic is a stub (shipped without enough input to run).

The command sequences come from ``examples/QE/README.txt``,
``examples/VASP/README.txt`` and the 42 per-tutorial ``README`` files,
cross-checked against the command tables in :mod:`htesp.cli`
(``SPECIAL_COMMANDS``, ``WORKFLOW_COMMANDS`` and ``NUMBERED``).

The artefact globs are what the step has to leave behind in the work directory
for it to count as done.  They were read off the ``reference*.tar.gz`` archives
that ship next to each tutorial and off the writers in
:mod:`htesp.workflow`, so "exited 0 and wrote nothing" is caught.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from tutorials.catalog import InputPatch, Loop, Seed, Step

#: material folders are ``R<mpid>-<compound>`` in both codes
# FIX: VASP inputs land in R<mpid>-<compound>/**relax**/, not directly
# under R<mpid>-<compound>/.  Two artifact globs omitted that level, so a
# download that had written INCAR/KPOINTS/POSCAR/POTCAR was still reported
# FAILED under the runner's "exited 0 and produced nothing" rule.
MAT = "R*-*"


@dataclass(frozen=True)
class TopicSpec:
    """Everything a topic contributes to a :class:`~tutorials.catalog.Tutorial`."""

    steps: tuple[Step, ...]
    loop: Loop | None = None
    depends_topics: tuple[str, ...] = ()
    seeds: tuple[Seed, ...] = ()
    stub: bool = False
    note: str = ""


def _inputs_artifacts(dft: str) -> tuple[str, ...]:
    """What ``download`` leaves behind: QE writes ``scf_dir``, VASP writes folders."""
    return ("scf_dir/scf-*.in",) if dft == "QE" else (f"{MAT}/relax/POSCAR",)


def _relax_artifacts(dft: str) -> tuple[str, ...]:
    return (f"{MAT}/relax/scf.in",) if dft == "QE" else (f"{MAT}/relax/INCAR",)


# --------------------------------------------------------------------------- #
#  topics 1-8: submission scripts and the database front ends
# --------------------------------------------------------------------------- #
#: Building the submission scripts is its own tutorial *and* a prerequisite of
#: every tutorial that submits, so it is defined once here and prepended by
#: `catalog._with_job_scripts`.
JOBSCRIPT_STEP = Step(
    "jobscript", "Build the run-*.sh submission scripts from batch.header",
    "jobscript", artifacts=("run-*.sh",),
    note="edit the job_script dictionary in config.json first")


def t_jobscript(dft: str) -> TopicSpec:
    return TopicSpec(steps=(JOBSCRIPT_STEP,))


def _search_download(dft: str, search: str, download: str, mode: str) -> tuple[Step, ...]:
    return (
        Step("search", f"Search the {mode} database", search,
             artifacts=("mpid-list.in",), needs_api_key=(search == "search"),
             timeout=3600),
        Step("download", f"Build {dft} inputs for the search hits", download,
             needs_potcar=(dft != "QE"),
             artifacts=_inputs_artifacts(dft) + ("mpid.in",),
             needs_api_key=(download == "download"), timeout=7200,
             input_patch=InputPatch(start=1, end=3, track="mpid-list.in"),
             note="only the first two hits are built, to keep the run short"),
    )


def t_mp_element(dft: str) -> TopicSpec:
    return TopicSpec(steps=_search_download(dft, "search", "download",
                                            "Materials Project (element mode)"),
                     note="needs MP_API_KEY in the environment")


def t_mp_chemsys(dft: str) -> TopicSpec:
    return TopicSpec(steps=_search_download(dft, "search", "download",
                                            "Materials Project (chemsys mode)"),
                     note="config.json sets download.info.mode = 'chemsys'")


def t_oqmd(dft: str) -> TopicSpec:
    return TopicSpec(steps=_search_download(dft, "oqmd-search", "oqmd-download", "OQMD"),
                     note="needs the optional qmpy-rester package")


def t_aflow(dft: str) -> TopicSpec:
    return TopicSpec(steps=_search_download(dft, "aflow-search", "aflow-download", "AFLOW"))


def t_magnetic(dft: str) -> TopicSpec:
    return TopicSpec(
        steps=_search_download(dft, "search", "download",
                               "Materials Project, magnetic configuration"),
        note="same as tutorials 2-5 with pwscf_in.magnetic = true in config.json")


def t_combine(dft: str) -> TopicSpec:
    return TopicSpec(
        steps=(Step("data-combine", "Merge the per-database downloads into one set",
                    "data-combine", artifacts=("mpid.in",)),),
        depends_topics=("mp-element", "oqmd", "aflow"),
        seeds=(Seed("@mp-element", ("scf_dir", MAT, "mpid.in"), overwrite=False),
               Seed("@oqmd", ("scf_dir", MAT, "mpid.in"), overwrite=False),
               Seed("@aflow", ("scf_dir", MAT, "mpid.in"), overwrite=False)),
        note="combines whatever the three database tutorials left in this folder")


def t_fromcif(dft: str) -> TopicSpec:
    return TopicSpec(
        steps=(Step("download", f"Build {dft} inputs from the .cif files", "download",
                    needs_potcar=(dft != "QE"),
                    artifacts=_inputs_artifacts(dft), timeout=3600),),
        seeds=(Seed("archive", ("*.cif",)),),
        note="config.json sets download.mode = 'fromcif'; the .cif files are "
             "unpacked from the tutorial's reference archive")


# --------------------------------------------------------------------------- #
#  topic 9: the relaxation hub
# --------------------------------------------------------------------------- #
def t_relax(dft: str) -> TopicSpec:
    relax = _relax_artifacts(dft)
    updated = ("scf_dir/scf-relax-*.in",) if dft == "QE" else ("mpid-list-not-relaxed.in",)
    return TopicSpec(
        steps=(
            Step("relax-submit", "Submit the first structural relaxation", "1",
                 submits=True, artifacts=relax, job_dirs=(f"{MAT}/relax",),
                 check_converged=True),
            Step("energy", "Collect total energies into econv.csv", "e0",
                 artifacts=("econv.csv",), needs_dft_output=True),
            Step("update-input", "Harvest the relaxed structure into a new input", "2",
                 artifacts=updated, needs_dft_output=True),
            Step("resubmit", "Resubmit the relaxation with the updated input", "3",
                 submits=True, artifacts=relax, job_dirs=(f"{MAT}/relax",),
                 check_converged=True),
            Step("energy-2", "Re-collect energies and check niteration", "e0",
                 artifacts=("econv.csv",), needs_dft_output=True),
        ),
        loop=Loop(("update-input", "resubmit", "energy-2"), max_cycles=4, probe="econv"),
        note="repeat 2 -> 3 -> e0 until econv.csv reports niteration < 3")


def t_convergence(dft: str) -> TopicSpec:
    return TopicSpec(steps=(
        Step("convtest", "Submit the cutoff / k-point convergence series", "convtest",
             needs_potcar=(dft != "QE"),
             submits=True, artifacts=(f"{MAT}/*/R*",), job_dirs=(f"{MAT}/*/R*",)),
        Step("extract", "Extract the convergence curve", "22",
             artifacts=("convergence_result",), needs_dft_output=True),
    ), note="config.json conv_test.param selects 'ecut' or 'kpoint'")


# --------------------------------------------------------------------------- #
#  topic 11: DFPT electron-phonon coupling -- QE only
# --------------------------------------------------------------------------- #
def t_elph(dft: str) -> TopicSpec:
    calc = f"{MAT}/calc"
    return TopicSpec(
        steps=(
            Step("even-kmesh", "Make the k-mesh even so the q-mesh divides it",
                 "change_k", artifacts=("kpoint.in",)),
            Step("create-inputs", "Build every downstream scf / el-ph input", "4",
                 artifacts=("elph_dir/elph-*.in", "kpath/kpath-*.dat"),
                 needs_dft_output=True),
            Step("fine-scf", "scf on the fine (doubled) k-mesh", "5", submits=True, after=("create-inputs",),
                 artifacts=(f"{calc}/scf.in",), job_dirs=(calc,)),
            Step("coarse-scf", "scf on the coarse mesh", "6", submits=True, after=("create-inputs",),
                 artifacts=(f"{calc}/scf.in",), job_dirs=(calc,)),
            Step("elph", "Electron-phonon coupling (ph.x)", "7", submits=True,
                 after=("create-inputs",),
                 artifacts=(f"{calc}/*.in",), job_dirs=(calc,), timeout=None),
            Step("q2r", "Force constants in real space (q2r.x)", "8",
                 artifacts=("q2r_dir/q2r-*.in",), needs_dft_output=True),
            Step("matdyn", "Phonon dispersion (matdyn.x)", "9",
                 artifacts=("matdyn_dir/matdyn-*.in",), needs_dft_output=True),
            Step("matdyn-dos", "Phonon DOS and EPC quantities (matdyn.x)", "10",
                 artifacts=("matdyn_dir/matdyn-*-dos.in",), needs_dft_output=True),
            Step("lambda", "Superconducting properties (lambda.x)", "11",
                 artifacts=("lambda*",), needs_dft_output=True),
            Step("phonband", "Process the phonon dispersion", "12",
                 artifacts=(f"{calc}/*",), needs_dft_output=True),
            Step("plot", "Plot alpha^2F and the phonon band structure", "19",
                 artifacts=("plots/*",), needs_dft_output=True,
                 input_patch=InputPatch(plot="phband")),
        ),
        depends_topics=("relax",),
        seeds=(Seed("@relax", ("R*-*", "scf_dir", "mpid.in", "input.in", "config.json"),
                    overwrite=False),),
        note="QE only; VASP has no DFPT counterpart, which is where the two "
             "tutorial numberings diverge")


def t_bands(dft: str) -> TopicSpec:
    if dft == "QE":
        bands = f"{MAT}/bands"
        steps = (
            Step("create-inputs", "Build the band / DOS inputs from the relaxed cell",
                 "4", artifacts=("kpath/kpath-*.dat",), needs_dft_output=True),
            Step("band-scf", "scf for the band structure", "13", submits=True, after=("create-inputs",),
                 artifacts=(f"{bands}/scf.in",), job_dirs=(bands,)),
            Step("band-nscf", "nscf along the high-symmetry path", "14", submits=True, after=("create-inputs",),
                 artifacts=(f"{bands}/band.in",), job_dirs=(bands,)),
            Step("band-post", "bands.x post-processing", "15", submits=True, after=("create-inputs",),
                 artifacts=(f"{bands}/*",), job_dirs=(bands,)),
            Step("dos-scf", "nscf on a dense mesh for the DOS", "16", submits=True, after=("create-inputs",),
                 artifacts=(f"{MAT}/dos/*",), job_dirs=(f"{MAT}/dos",)),
            Step("dos-post", "dos.x post-processing", "17", submits=True, after=("create-inputs",),
                 artifacts=(f"{MAT}/dos/*",), job_dirs=(f"{MAT}/dos",)),
            Step("pdos", "projwfc.x projected DOS", "18", submits=True, after=("create-inputs",),
                 artifacts=(f"{MAT}/dos/*",), job_dirs=(f"{MAT}/dos",)),
            Step("plot-bands", "Plot the band structure", "19",
                 artifacts=("plots/*",), needs_dft_output=True,
                 input_patch=InputPatch(plot="eband")),
            Step("plot-dos", "Plot the projected DOS", "19",
                 artifacts=("plots/*",), needs_dft_output=True,
                 input_patch=InputPatch(plot="pdos")),
        )
    else:
        steps = (
            Step("update-incar", "Rewrite INCAR with NSW = 0 and LCHARG = .TRUE.",
                 "download", artifacts=(f"{MAT}/relax/INCAR",)),
            Step("charge-density", "Charge-density run", "3", submits=True,
                 artifacts=(f"{MAT}/relax/INCAR",), job_dirs=(f"{MAT}/relax",)),
            Step("band-scf", "Band-structure run", "13", submits=True,
                 artifacts=(f"{MAT}/bands/*",), job_dirs=(f"{MAT}/bands",)),
            # FIX: vasp_process writes KPOINTS_band only when EIGENVAL exists,
            # i.e. only after a real VASP run, and `eigen` then moves it.  This
            # step therefore reads DFT output and cannot run under --dry-run;
            # it was failing with "No such file or directory: 'KPOINTS_band'"
            # instead of being skipped for the reason that is actually its own.
            Step("band-post", "Band post-processing", "15", submits=True,
                 artifacts=(f"{MAT}/bands/*",), job_dirs=(f"{MAT}/bands",),
                 needs_dft_output=True),
            Step("plot-bands", "Plot the band structure", "19",
                 artifacts=("plots/*",), needs_dft_output=True,
                 input_patch=InputPatch(plot="vasp-line")),
        )
    return TopicSpec(steps=steps, depends_topics=("relax",),
                     seeds=(Seed("@relax", ("R*-*", "scf_dir", "mpid.in", "econv.csv"),
                                 overwrite=True),))


# --------------------------------------------------------------------------- #
#  topics 13-21: the post-relaxation utilities
# --------------------------------------------------------------------------- #
def _from_relax(*extra: Seed) -> tuple[Seed, ...]:
    return (Seed("@relax", ("R*-*", "scf_dir", "mpid.in", "econv.csv"), overwrite=True),
            *extra)


def t_pressure(dft: str) -> TopicSpec:
    return TopicSpec(
        steps=(Step("pressure-input", "Build inputs for each pressure / volume point",
                    "pressure-input",
                    artifacts=("mpid-pressure-*.in", f"{MAT}/pressure"),
                    needs_dft_output=True),),
        depends_topics=("relax",), seeds=_from_relax(),
        note="the scale factors or pressures come from pressure.in")


def t_substitution(dft: str) -> TopicSpec:
    extra = ("scf_dir/scf-*-1.in",) if dft == "QE" else (f"{MAT}",)
    return TopicSpec(steps=(
        Step("substitute", "Enumerate the site substitutions", "29",
             needs_potcar=(dft != "QE"),
             artifacts=("mpid-substitute.in",) + extra),
    ), note="needs the optional bsym package")


def t_elastic(dft: str) -> TopicSpec:
    relax = _relax_artifacts(dft)
    return TopicSpec(
        steps=(
            Step("elastic-input", "Build the 24 deformed cells", "elastic-input",
                 needs_potcar=(dft != "QE"),
                 artifacts=("mpid-deformed.in",),
                 input_patch=InputPatch(start=1, end=2, track="mpid.in")),
            Step("relax-deformed", "Relax every deformed cell", "1", submits=True,
                 artifacts=relax, job_dirs=(f"{MAT}/relax",), check_converged=True,
                 input_patch=InputPatch(start=1, end=25, track="mpid-deformed.in")),
            Step("update-deformed", "Harvest the relaxed deformed cells", "2",
                 artifacts=("mpid-list-not-relaxed.in",), needs_dft_output=True,
                 input_patch=InputPatch(start=1, end=25, track="mpid-deformed.in")),
            Step("resubmit-deformed", "Resubmit the deformed relaxations", "3",
                 submits=True, artifacts=relax, job_dirs=(f"{MAT}/relax",),
                 input_patch=InputPatch(start=1, end=25, track="mpid-deformed.in")),
            Step("compute-elastic", "Fit the elastic tensor", "compute-elastic",
                 artifacts=("elastic.csv",), needs_dft_output=True,
                 input_patch=InputPatch(start=1, end=2, track="mpid.in")),
        ),
        depends_topics=("relax",), seeds=_from_relax(),
        note="input.in switches between mpid.in and mpid-deformed.in mid-tutorial; "
             "the runner rewrites it per step")


def t_hull(dft: str) -> TopicSpec:
    return TopicSpec(
        steps=(
            Step("energy", "Collect the formation energies", "e0",
                 artifacts=("econv.csv",), needs_dft_output=True),
            Step("phase-diagram", "Build the convex hull", "pd",
                 artifacts=("*.pdf",), needs_dft_output=True, needs_api_key=True),
        ),
        depends_topics=("relax",), seeds=_from_relax(),
        note="the hull needs the competing phases from the Materials Project")


def t_phonopy(dft: str) -> TopicSpec:
    ph = f"{MAT}/phonopy"
    return TopicSpec(
        steps=(
            Step("displacements", "Build the displaced supercells and submit them",
                 "phono1", submits=True, needs_dft_output=True,
                 artifacts=(f"{ph}/phonopy_disp.yaml",), job_dirs=(f"{ph}/R*",)),
            Step("force-constants", "Compute the force constants", "phono2",
                 artifacts=(f"{ph}/FORCE_SETS",), needs_dft_output=True),
            Step("thermal", "Thermodynamic properties", "phono3",
                 artifacts=(f"{ph}/mesh.conf",), needs_dft_output=True),
            Step("phonon-bands", "Phonon band structure", "phono4",
                 artifacts=(f"{ph}/band.conf",), needs_dft_output=True),
        ),
        depends_topics=("relax",), seeds=_from_relax(),
        note="needs the phonopy executable on PATH")


def t_eos(dft: str) -> TopicSpec:
    return TopicSpec(
        steps=(
            Step("pressure-input", "Build the isotropically scaled cells",
                 "pressure-input", artifacts=(f"{MAT}/pressure",),
                 needs_dft_output=True),
            Step("relax-volumes", "Relax at fixed volume", "26", submits=True,
                 artifacts=(f"{MAT}/pressure/*",),
                 job_dirs=(f"{MAT}/pressure/R*/relax",), check_converged=True),
            Step("ev-collect", "Collect the energy-volume curve", "ev-collect",
                 artifacts=(f"{MAT}/pressure/e-v.dat",), needs_dft_output=True),
            Step("eos-bm", "Birch-Murnaghan fit", "eos-bm",
                 artifacts=("eos-fit.dat",), needs_dft_output=True),
            Step("eos-vinet", "Vinet fit", "eos-vinet",
                 artifacts=("eos-fit.dat",), needs_dft_output=True),
        ),
        depends_topics=("relax",), seeds=_from_relax(),
        note="mainprogram 26 has to be repeated until every volume is converged")


def t_wannier(dft: str) -> TopicSpec:
    if dft == "QE":
        prepare = (Step("create-inputs", "Build the scf / nscf inputs", "4",
                        artifacts=("kpath/kpath-*.dat",), needs_dft_output=True),)
        epw_art = ("scf_dir/*-nscf.in",)
    else:
        prepare = ()
        epw_art = (f"{MAT}/epw/POSCAR",)
    return TopicSpec(
        steps=prepare + (
            Step("epw1", "Build the Wannier90 / EPW inputs", "epw1",
                 artifacts=epw_art, needs_dft_output=True),
            Step("jobscript", "Build the Wannier submission script", "jobscript",
                 artifacts=("run-*.sh",)),
            Step("wann-file", "Wannierise with the projections from projection.in",
                 "wann-file", artifacts=("epw_dir",), needs_dft_output=True,
                 after=("epw1",)),
        ),
        depends_topics=("relax",), seeds=_from_relax(),
        note="the disentanglement windows come from the band/pDOS tutorial and "
             "from wannier90.json")


def t_charge(dft: str) -> TopicSpec:
    return TopicSpec(
        steps=(Step("charge-input", "Build one input per net charge", "charge-input",
                    artifacts=("mpid-charge.in",)),),
        depends_topics=("relax",), seeds=_from_relax(),
        note="the charges come from charge.in; +ve removes electrons in QE")


def t_magorder(dft: str) -> TopicSpec:
    return TopicSpec(
        steps=(Step("magenum", "Enumerate the magnetic orderings", "magenum",
                    needs_potcar=(dft != "QE"),
                    artifacts=("mpid-magnetic.in",), needs_enumlib=True),),
        depends_topics=("relax",), seeds=_from_relax(),
        note="config.json magmom.type must be 'ordering'")


def t_fermisurface(dft: str) -> TopicSpec:
    return TopicSpec(
        steps=(Step("fermisurface", "3D Fermi surface with IFermi", "fermisurface",
                    artifacts=(f"{MAT}/*",), needs_dft_output=True),),
        depends_topics=("relax",), seeds=_from_relax(), stub=True,
        note="examples/VASP/tutorial21 ships only ifermi.tar.gz -- no config.json, "
             "no input.in and no vasprun.xml, so this tutorial cannot run as "
             "shipped; the step list is a best-effort reconstruction")


#: topic key -> builder
BUILDERS = {
    "jobscript": t_jobscript, "mp-element": t_mp_element, "mp-chemsys": t_mp_chemsys,
    "oqmd": t_oqmd, "aflow": t_aflow, "magnetic": t_magnetic, "combine": t_combine,
    "fromcif": t_fromcif, "relax": t_relax, "convergence": t_convergence,
    "elph": t_elph, "bands": t_bands, "pressure": t_pressure,
    "substitution": t_substitution, "elastic": t_elastic, "hull": t_hull,
    "phonopy": t_phonopy, "eos": t_eos, "wannier": t_wannier, "charge": t_charge,
    "magorder": t_magorder, "fermisurface": t_fermisurface,
}


def spec_for(topic: str, dft: str) -> TopicSpec:
    """Return the :class:`TopicSpec` for *topic* under *dft* (``"QE"``/``"VASP"``)."""
    try:
        builder = BUILDERS[topic]
    except KeyError as exc:                     # pragma: no cover - guarded by tests
        raise KeyError(f"no step list for topic {topic!r}") from exc
    return builder(dft)
