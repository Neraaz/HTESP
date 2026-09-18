#!/usr/bin/env python
"""``mainprogram`` -- the HTESP command line.

This replaces the 469-line ``if/elif`` chain of the original
``src/mainprogram.py``.  The behaviour the user sees is unchanged:

    mainprogram 4
    mainprogram download
    mainprogram phono1

still read ``input.in``, still write the banner to ``log`` and still act on
``[start, end)`` of the tracking file.  What changed underneath:

* commands are a dictionary, not a chain of 90 ``elif`` branches;
* the numbered processes call :class:`htesp.workflow.HTESPWorkflow` methods
  **in this process** instead of ``os.system("<bash script> ...")``, so a
  failure is an exception with a traceback rather than a discarded exit code;
* ``input.in`` is parsed by :class:`htesp.inputin.InputIn`, which reports what
  is wrong with a short file instead of raising ``IndexError`` on ``lines[4]``
  or leaving ``start``/``end`` unbound;
* ``mainprogram 19`` iterates over the *list* of plot types.  The old code
  created ``input.in`` with the string ``'phband'`` and then iterated over it,
  launching six plot jobs named ``p``, ``h``, ``b``, ``a``, ``n``, ``d``;
* the four long help blocks live in :mod:`htesp.help_text`, which
  ``docs/command.rst`` is generated from;
* ``mainprogram`` with no argument prints the usage instead of raising
  ``IndexError`` on ``sys.argv[1]``.

New, additive options (every existing invocation still works):

    --workers N     size of the per-material process pool (default: min(cpu, 8))
    --dry-run       build every input file but never call the scheduler
    --root DIR      run against DIR instead of the working directory
    --config FILE   use FILE instead of searching for config.json
    -v/--verbose    debug logging
    --list          list every command with a one-line description
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Any, Callable

from htesp.banner import __version__, banner
from htesp.help_text import HELP, PROCESS_SUMMARY, SUMMARY
from htesp.inputin import InputIn, InputInError

LOG = logging.getLogger("htesp")

#: commands that take ``(start, end, track)`` and an extra positional argument.
#: name -> (workflow method, extra args appended after the track file)
WORKFLOW_COMMANDS: dict[str, tuple[str, tuple[Any, ...]]] = {
    # named commands
    "compound": ("info_scan", ()),
    "pressure-input": ("pressure_input", ()),
    "charge-input": ("charge_input", ()),
    "fermisurface": ("ifermi_scan", ()),
    "checkph": ("phcheck_scan", ()),
    "checkfreq": ("checkfreq_scan", ()),
    "change_k": ("double_kmesh", ()),
    "magmom_extract": ("magmom_extract", ()),
    # phonopy family
    "e0": ("phonopy_scan", (0,)),
    "phono1": ("phonopy_scan", (1,)),
    "phono2": ("phonopy_scan", (2,)),
    "phono3": ("phonopy_scan", (3,)),
    "phono4": ("phonopy_scan", (4,)),
    "phono5": ("phonopy_scan", (5,)),
    "phono-qha": ("phonopy_scan", ("vp-ph-qha",)),
    "eos-bm": ("phonopy_scan", ("eos-bm",)),
    "eos-vinet": ("phonopy_scan", ("eos-vinet",)),
    "ev-collect": ("phonopy_scan", ("ev-collect",)),
    # FIX: the original dispatched these to 'vp-ph2' ... 'vp-ph5', which
    # phonopy-scan does not implement, so all four printed "Bad input" and
    # did nothing.  They are the pressure variants of phono1..phono4.
    "phono1-pressure": ("phonopy_scan", ("phono1-pressure",)),
    "phono2-pressure": ("phonopy_scan", ("phono2-pressure",)),
    "phono3-pressure": ("phonopy_scan", ("phono3-pressure",)),
    "phono4-pressure": ("phonopy_scan", ("phono4-pressure",)),
    # EPW / Wannier90 family
    "epw1": ("epw_bash_scripts", ("epw1",)),
    "qe-ph": ("epw_bash_scripts", ("epw2",)),
    "epw2": ("epw_bash_scripts", ("epw3",)),
    "epw3": ("epw_bash_scripts", ("epw4",)),
    "epw4": ("epw_bash_scripts", ("proj",)),
    "epw5": ("epw_bash_scripts", ("band_wann2",)),
    "wann-scdm": ("epw_bash_scripts", ("band_wann", "scdm")),
    "wann-file": ("epw_bash_scripts", ("band_wann", "fromfile")),
    "wann-random": ("epw_bash_scripts", ("band_wann", "random")),
    "epw-scdm": ("epw_bash_scripts", ("epw", "scdm")),
    "epw-file": ("epw_bash_scripts", ("epw", "fromfile")),
    "epw-random": ("epw_bash_scripts", ("epw", "random")),
    # WannierTools
    "wt1": ("wt_bash_scripts", ("wt1-b",)),
    "wt2": ("wt_bash_scripts", ("wt1-s",)),
    # single-mode helper
    "singlemode": ("distortion_help", None),   # None == takes no range
    "history": ("history", None),
}

#: numbered processes -> (workflow method, extra args)
NUMBERED: dict[int, tuple[str, tuple[Any, ...]]] = {
    1: ("relax_scan", ()),
    2: ("further_relax_input", ("first",)),
    3: ("further_relax_scan", ()),
    4: ("create_inputs", ()),          # nkpt appended from input.in
    5: ("fine_scan", ()),
    6: ("coarse_scan", ()),
    7: ("ph_scan", ()),
    8: ("q2r_scan", ()),
    9: ("matdyn_scan", ()),
    10: ("matdyn_dos_scan", ()),
    11: ("lambda_scan", ()),
    12: ("phonband_scan", ()),
    13: ("bandscf_scan", ()),
    14: ("band_scan", ()),
    15: ("bandp_scan", ()),
    16: ("dos_scan", ()),
    17: ("dosp_scan", ()),
    18: ("pdos_scan", ()),
    19: ("plot_scan", ()),             # handled specially (plot types)
    20: ("clean_scan", ()),
    21: ("extract_scan", ()),
    23: ("dynmat_scan", ()),
    24: ("distortion_relax_scan", ()),
    25: ("distortion_energy_scan", ()),
    26: ("pressure_relax_scan", ()),
    27: ("pressure_ph_scan", ()),
    28: ("pressure_reset", ()),
    29: ("sitesub_scan", ()),
}

#: processes that are not workflow scans
SPECIAL_NUMBERED = {0, 22}


class CommandError(RuntimeError):
    """A command could not be carried out; the message is for the user."""


# --------------------------------------------------------------------------- #
#  context
# --------------------------------------------------------------------------- #
class Context:
    """Everything a command needs: the parsed ``input.in``, config and options."""

    def __init__(self, args: argparse.Namespace):
        self.root = Path(args.root).resolve()
        self.workers = args.workers
        self.dry_run = args.dry_run
        self.force = getattr(args, "force", False)
        init = getattr(args, "init_header", None)
        self.init_header = init.lower() if init else None
        self.verbose = args.verbose
        if args.config:
            os.environ["HTESP_CONFIG"] = str(Path(args.config).resolve())
        self._workflow = None
        self._input = None
        self._config = None

    @property
    def config(self) -> dict:
        if self._config is None:
            from htesp.config import config as _config
            self._config = _config(self.root)
        return self._config

    @property
    def input(self) -> InputIn:
        if self._input is None:
            calc = str(self.config["download"]["inp"]["calc"]).upper()
            self._input = InputIn.load_or_create(self.root / "input.in", dft=calc)
        return self._input

    @property
    def has_workflow(self) -> bool:
        """True once the workflow layer has actually been constructed."""
        return self._workflow is not None

    @property
    def workflow(self):
        if self._workflow is None:
            from htesp.workflow import HTESPWorkflow
            self._workflow = HTESPWorkflow(
                root=self.root, workers=self.workers, dry_run=self.dry_run,
                log_level=logging.DEBUG if self.verbose else logging.INFO)
        return self._workflow

    @property
    def range(self) -> tuple[int, int, str]:
        inp = self.input
        return inp.start, inp.end, inp.track

    def write_log(self) -> None:
        """Write the banner, the resolved range and the configuration to ``log``.

        Which ``config.json`` was in force is the first thing anyone asks when
        a run produced unexpected numbers: the file is searched for up the
        directory tree, so a stage running in ``R<id>-<comp>/relax/`` may be
        reading a different one from the one next to ``input.in``. Recording
        the resolved path makes the answer part of the run rather than
        something to reconstruct afterwards.
        """
        from htesp.config import config_path
        inp = self.input
        resolved = config_path(self.root)
        source = str(resolved) if resolved else "(packaged default)"
        LOG.info("configuration: %s", source)
        try:
            with open(self.root / "log", "w") as handle:
                handle.write(banner())
                handle.write("\n")
                handle.write("#" * 131 + "\n")
                handle.write("# " + inp.summary() + "\n")
                handle.write("# configuration: " + source + "\n")
                handle.write("#" * 131 + "\n")
        except OSError as exc:            # a read-only directory must not abort the run
            LOG.warning("could not write the log file: %s", exc)

    def warn_about_range(self) -> None:
        for problem in self.input.check_track_file(self.root):
            LOG.warning("input.in: %s", problem)


# --------------------------------------------------------------------------- #
#  non-workflow commands
# --------------------------------------------------------------------------- #
def cmd_jobscript(ctx: Context, rest: list[str]) -> int:
    """``mainprogram jobscript`` -- build the submission scripts.

    ``--init-header qe|vasp`` instead writes a starting ``batch.header`` by
    asking SLURM and Lmod what this machine has.  The shipped examples say
    ``--partition=dense``, which exists on one cluster and nowhere else, so
    copying one produces a job the scheduler rejects before anything runs.
    """
    if ctx.init_header:
        from htesp import batch_header
        from htesp.config import config as _config

        name = _config(ctx.root).get("job_script", {}).get("batch", "batch.header")
        return batch_header.write(ctx.root / name, ctx.init_header,
                                  force=ctx.force)
    from htesp import generate_submission
    generate_submission.main()
    return 0


def cmd_search(ctx: Context, rest: list[str]) -> int:
    from htesp import element_extract
    from htesp.config import config_path

    mode = ctx.config["download"]["mode"] if config_path(ctx.root) else "element"
    if mode not in ("fromcif", "fromvasp"):
        LOG.info("searching the Materials Project database in %r mode", mode)
    existing = ctx.root / "mpid-list.in"
    if existing.is_file():
        backup = ctx.root / "mpid-list-2.in"
        LOG.info("an mpid-list.in was found; renaming it to %s", backup.name)
        existing.replace(backup)
    element_extract.main()
    if mode not in ("fromcif", "fromvasp"):
        LOG.info("adjust the indices in input.in, then run 'mainprogram download'")
    else:
        LOG.info("run 'mainprogram download' to build the inputs from .cif/.vasp files")
    return 0


def cmd_download(ctx: Context, rest: list[str]) -> int:
    download = ctx.config["download"]
    mode = download.get("mode", "")
    calc = download["inp"]["calc"]
    if mode == "fromcif":
        from htesp import cif_to_gsinput
        LOG.info("building %s inputs from the .cif files", calc)
        if download["inp"].get("use_cif2cell"):
            LOG.info("this mode needs the cif2cell package (pip install cif2cell)")
        cif_to_gsinput.main(calc)
        LOG.info("inputs are in %s",
                 "Rmp-<id>-<compound>/" if str(calc).lower() == "vasp" else "scf_dir/")
        return 0
    if mode == "fromvasp":
        from htesp import poscar_to_vasp
        LOG.info("building %s inputs from the .vasp files", calc)
        poscar_to_vasp.main()
        return 0
    start, end, track = ctx.range
    ctx.workflow.download_input(start, end, track)
    return 0


def cmd_convtest(ctx: Context, rest: list[str]) -> int:
    from htesp import convergence_test
    convergence_test.main(["calculate"])
    return 0


def cmd_convextract(ctx: Context, rest: list[str]) -> int:
    from htesp import convergence_test
    convergence_test.main(["extract"])
    return 0


def cmd_oqmd_search(ctx: Context, rest: list[str]) -> int:
    from htesp import oqmd_extract
    oqmd_extract.main(["search"])
    return 0


def cmd_oqmd_download(ctx: Context, rest: list[str]) -> int:
    from htesp import oqmd_extract
    oqmd_extract.main(["download"])
    return 0


def cmd_aflow_search(ctx: Context, rest: list[str]) -> int:
    from htesp import aflow_extract
    aflow_extract.main(["search"])
    return 0


def cmd_aflow_download(ctx: Context, rest: list[str]) -> int:
    from htesp import aflow_extract
    aflow_extract.main(["download"])
    return 0


def cmd_data_combine(ctx: Context, rest: list[str]) -> int:
    from htesp import structure_group
    structure_group.main()
    return 0


def cmd_elastic_input(ctx: Context, rest: list[str]) -> int:
    from htesp import elastic
    elastic.main(["input"])
    return 0


def cmd_compute_elastic(ctx: Context, rest: list[str]) -> int:
    from htesp import elastic
    elastic.main(["compute_elastic"])
    return 0


def cmd_magenum(ctx: Context, rest: list[str]) -> int:
    from htesp import magnetic
    magnetic.main()
    return 0


def cmd_phase_diagram(ctx: Context, rest: list[str]) -> int:
    from htesp import pymatgen_phase_diagram
    pymatgen_phase_diagram.main()
    return 0


def cmd_primtoconv(ctx: Context, rest: list[str]) -> int:
    """Convert every relaxed cell in the range to its conventional setting."""
    import shutil
    import subprocess

    start, end, track = ctx.range
    materials = ctx.workflow.materials(start, end, track)
    for material in materials:
        relax = material.sub("relax")
        if not relax.is_dir():
            LOG.warning("%s: no relax/ directory, skipping", material.name)
            continue
        backup = material.dir / "relax_prim"
        if not backup.exists():
            shutil.copytree(relax, backup)
        subprocess.run([sys.executable, "-m", "htesp.vasp_process", "conventional"],
                       cwd=relax, check=True)
    return 0


def cmd_mkdirs(ctx: Context, rest: list[str]) -> int:
    ctx.workflow.ensure_project_dirs("scf_dir", "elph_dir", "matdyn_dir", "q2r_dir")
    return 0


def cmd_plot(ctx: Context, rest: list[str]) -> int:
    """Process 19 -- one plot-scan per plot type listed on input.in line 5."""
    start, end, track = ctx.range
    inp = ctx.input
    # FIX: the original wrote plot_type as the *string* 'phband' when it
    # created input.in and then iterated over it character by character.
    for kind in inp.plot_types:
        LOG.info("plotting: %s", kind)
        ctx.workflow.plot_scan(start, end, track, inp.nkpt, kind)
    return 0


def cmd_config_init(ctx: Context, rest: list[str]) -> int:
    """Write the packaged default ``config.json`` into the project directory.

    Without a ``config.json`` of its own a project silently runs on the
    packaged default -- every key has a value, so nothing fails, and the
    cutoffs and k-point density are whatever the package ships rather than
    what the study needs.  This is the counterpart to ``config-validate``:
    it produces the file that command then reports on.
    """
    import shutil
    from htesp.config import CONFIG_FILENAME, DEFAULT_CONFIG_PATH

    force = ctx.force
    targets = [item for item in rest if not item.startswith("-")]
    target = Path(targets[0]) if targets else ctx.root / CONFIG_FILENAME
    if target.is_dir():
        target = target / CONFIG_FILENAME
    if target.exists() and not force:
        LOG.error("%s already exists; pass --force to overwrite it", target)
        return 2
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(DEFAULT_CONFIG_PATH, target)
    print(f"wrote {target}")
    print("  the Materials Project key is NOT read from this file: "
          "set $MP_API_KEY or ~/.config/htesp/credentials")
    print(f"  check it with: mainprogram config-validate")
    return 0


def cmd_config_validate(ctx: Context, rest: list[str]) -> int:
    from htesp.config import config_path, validate
    from htesp.config import CONFIG_FILENAME, DEFAULT_CONFIG_PATH, SEARCH_DEPTH
    path = config_path(ctx.root)
    print(f"configuration: {path or '(packaged default)'}")
    print(f"  searched:    $HTESP_CONFIG, then {CONFIG_FILENAME} in {ctx.root} "
          f"and up to {SEARCH_DEPTH} parent directories")
    print(f"  merged over: {DEFAULT_CONFIG_PATH}")
    if "-v" in rest or "--verbose" in rest or ctx.verbose:
        import json as _json
        print(_json.dumps(ctx.config, indent=2))
    problems = validate(ctx.config)
    if not problems:
        print("no problems found")
        return 0
    print(f"{len(problems)} problem(s):")
    for problem in problems:
        print(f"  - {problem}")
    return 1


def cmd_help_block(name: str) -> Callable[[Context, list[str]], int]:
    def run(ctx: Context, rest: list[str]) -> int:
        print(HELP[name])
        return 0
    return run


SPECIAL_COMMANDS: dict[str, Callable[[Context, list[str]], int]] = {
    "jobscript": cmd_jobscript,
    "search": cmd_search,
    "download": cmd_download,
    "convtest": cmd_convtest,
    "oqmd-search": cmd_oqmd_search,
    "oqmd-download": cmd_oqmd_download,
    "aflow-search": cmd_aflow_search,
    "aflow-download": cmd_aflow_download,
    "data-combine": cmd_data_combine,
    "elastic-input": cmd_elastic_input,
    "compute-elastic": cmd_compute_elastic,
    "magenum": cmd_magenum,
    "pd": cmd_phase_diagram,
    "primtoconv": cmd_primtoconv,
    "config-init": cmd_config_init,
    "config-validate": cmd_config_validate,
    "basicinfo": cmd_help_block("basicinfo"),
    "process-info": cmd_help_block("process-info"),
    "epw-info": cmd_help_block("epw-info"),
    "wt-info": cmd_help_block("wt-info"),
}

#: commands that never need input.in, a tracking file or the workflow layer
NO_INPUT_NEEDED = {"basicinfo", "process-info", "epw-info", "wt-info",
                   "config-init", "config-validate", "jobscript"}


# --------------------------------------------------------------------------- #
#  dispatch
# --------------------------------------------------------------------------- #
def all_commands() -> list[str]:
    names = set(SPECIAL_COMMANDS) | set(WORKFLOW_COMMANDS)
    return sorted(names)


def print_usage() -> None:
    print(banner())
    print("usage: mainprogram <process> [options]\n")
    print("  <process> is a number 0-29 or one of the named commands below.")
    print("  Run 'mainprogram basicinfo' for the introduction,")
    print("  'mainprogram process-info' for the numbered processes,")
    print("  'mainprogram --list' for a one-line summary of everything.\n")
    print("named commands:")
    for name in all_commands():
        print(f"  {name:<20} {SUMMARY.get(name, '')}")


def print_list() -> None:
    print("numbered processes")
    for number, text in sorted(PROCESS_SUMMARY.items(), key=lambda kv: int(kv[0])):
        print(f"  {number:>3}  {text}")
    print("\nnamed commands")
    for name in all_commands():
        print(f"  {name:<20} {SUMMARY.get(name, '')}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="mainprogram", add_help=False,
        description="HTESP -- High Throughput Electronic Structure Package")
    parser.add_argument("process", nargs="?", help="process number or command name")
    parser.add_argument("rest", nargs="*", help="extra arguments for the command")
    parser.add_argument("--workers", type=int, default=None,
                        help="per-material process pool size (default: min(cpu, 8))")
    parser.add_argument("--dry-run", action="store_true",
                        help="prepare every input file but do not submit anything")
    parser.add_argument("--root", default=".", help="project directory (default: .)")
    parser.add_argument("--config", default=None,
                        help="config.json to use instead of searching for one")
    parser.add_argument("--init-header", metavar="CODE", default=None,
                        choices=("qe", "vasp", "QE", "VASP"),
                        help="with 'jobscript': write a starting batch.header "
                             "for this code, filled in from SLURM and Lmod")
    parser.add_argument("--force", action="store_true",
                        help="overwrite a file the command would otherwise "
                             "refuse to replace (config-init)")
    parser.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    parser.add_argument("--list", action="store_true",
                        help="list every command and exit")
    parser.add_argument("--version", action="store_true")
    parser.add_argument("-h", "--help", action="store_true")
    return parser


def run_command(ctx: Context, name: str, rest: list[str]) -> int:
    if name in SPECIAL_COMMANDS:
        return SPECIAL_COMMANDS[name](ctx, rest)

    if name in WORKFLOW_COMMANDS:
        method_name, extra = WORKFLOW_COMMANDS[name]
        method = getattr(ctx.workflow, method_name)
        if extra is None:                      # takes no material range
            result = method()
            if isinstance(result, str):
                print(result)
            return 0
        start, end, track = ctx.range
        method(start, end, track, *extra)
        return 0

    raise CommandError(
        f"unknown command {name!r}.  Run 'mainprogram --list' for the full list, "
        "or 'mainprogram basicinfo' for the introduction."
    )


def run_numbered(ctx: Context, number: int) -> int:
    if number == 0:
        return cmd_mkdirs(ctx, [])
    if number == 19:
        return cmd_plot(ctx, [])
    if number == 22:
        return cmd_convextract(ctx, [])
    if number not in NUMBERED:
        raise CommandError(
            f"process {number} does not exist.  Valid numbers are 0-29 "
            "(22 needs 'mainprogram convtest' first); run 'mainprogram process-info'."
        )
    method_name, extra = NUMBERED[number]
    start, end, track = ctx.range
    method = getattr(ctx.workflow, method_name)
    if number == 4:                            # create-inputs takes nkpt
        method(start, end, track, ctx.input.nkpt)
    else:
        method(start, end, track, *extra)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(message)s")

    if args.version:
        print(f"HTESP {__version__}")
        return 0
    if args.list:
        print_list()
        return 0
    # FIX: the original did sys.argv[1] unconditionally -> IndexError
    if args.help or not args.process:
        print_usage()
        return 0 if args.help else 1

    ctx = Context(args)
    name = args.process

    try:
        if name not in NO_INPUT_NEEDED:
            ctx.write_log()
            ctx.warn_about_range()
        if name.lstrip("-").isdigit():
            status = run_numbered(ctx, int(name))
        else:
            status = run_command(ctx, name, list(args.rest))
        # FIX: a stage in which every material failed used to exit 0, because
        # the bash layer discarded every exit code.  The next stage then ran
        # on nothing and reported success too.
        if status == 0 and ctx.has_workflow and ctx.workflow.failed_count:
            LOG.error("%d material(s) failed:\n%s",
                      ctx.workflow.failed_count, ctx.workflow.failure_summary())
            return 1
        return status
    except InputInError as exc:
        LOG.error("%s", exc)
        return 2
    except CommandError as exc:
        LOG.error("%s", exc)
        return 2
    except FileNotFoundError as exc:
        LOG.error("%s", exc)
        return 2
    except KeyboardInterrupt:
        LOG.error("interrupted")
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
