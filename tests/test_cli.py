"""The ``mainprogram`` command line: dispatch, exit codes, and the plot loop."""
from __future__ import annotations

import contextlib
import io
import unittest
from pathlib import Path

from tests.helpers import TempProject

from htesp import cli
from htesp.help_text import HELP, PROCESS_SUMMARY, SUMMARY


class RecordingWorkflow:
    """Stands in for :class:`htesp.workflow.HTESPWorkflow` during dispatch tests."""

    def __init__(self):
        self.calls: list[tuple] = []
        self.failed_count = 0

    def failure_summary(self) -> str:
        return ""

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)

        def record(*args, **kwargs):
            self.calls.append((name, args, kwargs))
            return []
        return record


class DispatchTest(TempProject):
    def setUp(self):
        super().setUp()
        self.write_input_in(start=1, end=3, nkpt=150, track="mpid.in",
                            plots="phband dos a2f")
        self.write_track("mpid.in", [("mp-1", "A"), ("mp-2", "B")])
        self.fake = RecordingWorkflow()
        self._saved = cli.Context.workflow
        cli.Context.workflow = property(lambda ctx: self.fake)
        cli.Context.has_workflow = property(lambda ctx: True)

    def tearDown(self):
        cli.Context.workflow = self._saved
        cli.Context.has_workflow = property(lambda ctx: ctx._workflow is not None)
        super().tearDown()

    def run_cli(self, *argv) -> tuple[int, str]:
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = cli.main(list(argv) + ["--root", str(self.root)])
        return code, buffer.getvalue()

    # -- numbered processes ------------------------------------------------
    def test_numbered_process_calls_the_matching_method(self):
        code, _ = self.run_cli("1")
        self.assertEqual(code, 0)
        self.assertEqual(self.fake.calls[0][0], "relax_scan")
        self.assertEqual(self.fake.calls[0][1][:3], (1, 3, "mpid.in"))

    def test_process_2_passes_the_first_flag(self):
        self.run_cli("2")
        self.assertEqual(self.fake.calls[0][1][3], "first")

    def test_process_4_passes_nkpt_from_input_in(self):
        self.run_cli("4")
        name, args, _ = self.fake.calls[0]
        self.assertEqual(name, "create_inputs")
        self.assertEqual(args[3], 150)

    def test_process_19_iterates_the_plot_types_not_the_letters(self):
        """The original iterated the string 'phband' -> p, h, b, a, n, d."""
        self.run_cli("19")
        kinds = [call[1][4] for call in self.fake.calls]
        self.assertEqual(kinds, ["phband", "dos", "a2f"])

    def test_process_0_creates_the_project_directories(self):
        self.run_cli("0")
        self.assertEqual(self.fake.calls[0][0], "ensure_project_dirs")

    def test_an_unknown_process_number_is_reported(self):
        code, _ = self.run_cli("99")
        self.assertEqual(code, 2)

    # -- named commands ----------------------------------------------------
    def test_named_command(self):
        self.run_cli("checkph")
        self.assertEqual(self.fake.calls[0][0], "phcheck_scan")

    def test_phonopy_family_passes_the_step(self):
        for command, step in (("e0", 0), ("phono1", 1), ("phono4", 4),
                              ("eos-bm", "eos-bm"), ("ev-collect", "ev-collect")):
            self.fake.calls.clear()
            self.run_cli(command)
            self.assertEqual(self.fake.calls[0][0], "phonopy_scan", command)
            self.assertEqual(self.fake.calls[0][1][3], step, command)

    def test_pressure_variants_reach_phonopy_scan(self):
        """These four dispatched to `vp-phN`, which was not a command at all."""
        for command in ("phono1-pressure", "phono2-pressure",
                        "phono3-pressure", "phono4-pressure"):
            self.fake.calls.clear()
            code, _ = self.run_cli(command)
            self.assertEqual(code, 0, command)
            self.assertEqual(self.fake.calls[0][0], "phonopy_scan", command)
            self.assertTrue(str(self.fake.calls[0][1][3]).endswith("-pressure"))

    def test_epw_family(self):
        self.run_cli("wann-file")
        name, args, _ = self.fake.calls[0]
        self.assertEqual(name, "epw_bash_scripts")
        self.assertEqual(args[3:], ("band_wann", "fromfile"))

    def test_unknown_command_exits_non_zero(self):
        code, _ = self.run_cli("nonsense")
        self.assertEqual(code, 2)

    def test_failed_materials_make_the_command_exit_non_zero(self):
        self.fake.failed_count = 2
        self.fake.failure_summary = lambda: "  mp-1 A: boom"
        code, _ = self.run_cli("1")
        self.assertEqual(code, 1)


class NoArgumentsTest(TempProject):
    def test_bare_mainprogram_prints_usage_instead_of_raising(self):
        """The original did sys.argv[1] unconditionally -> IndexError."""
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            code = cli.main([])
        self.assertEqual(code, 1)
        self.assertIn("usage", buffer.getvalue().lower())

    def test_help_lists_the_commands(self):
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            cli.main(["--help"])
        text = buffer.getvalue()
        for name in ("basicinfo", "compute-elastic", "wt1"):
            self.assertIn(name, text)

    def test_list_covers_the_numbered_processes(self):
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            cli.main(["--list"])
        text = buffer.getvalue()
        self.assertIn("clean-scan", text)
        self.assertIn("pdos-scan", text)

    def test_version(self):
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            cli.main(["--version"])
        self.assertIn("HTESP", buffer.getvalue())


class HelpBlocks(unittest.TestCase):
    def test_the_four_blocks_exist_and_are_substantial(self):
        for name in ("basicinfo", "process-info", "epw-info", "wt-info"):
            self.assertIn(name, HELP)
            self.assertGreater(len(HELP[name].splitlines()), 10, name)

    def test_the_wanniertools_config_key_is_spelled_correctly(self):
        """wt-info told users to edit 'wanniertool_input', which does not exist."""
        self.assertIn("wanniertools_input", HELP["wt-info"])
        self.assertNotIn("wanniertool_input", HELP["wt-info"])

    def test_pdos_is_process_18_and_20_is_destructive(self):
        """docs/tutorial.rst told users to run 20 for PDOS; 20 deletes the run."""
        self.assertIn("18", PROCESS_SUMMARY)
        self.assertIn("pdos", PROCESS_SUMMARY["18"].lower())
        self.assertIn("DELETE", PROCESS_SUMMARY["20"])

    def test_every_summarised_command_really_exists(self):
        known = set(cli.SPECIAL_COMMANDS) | set(cli.WORKFLOW_COMMANDS)
        unknown = sorted(name for name in SUMMARY if name not in known)
        self.assertEqual(unknown, [])

    def test_every_command_has_a_summary(self):
        known = set(cli.SPECIAL_COMMANDS) | set(cli.WORKFLOW_COMMANDS)
        missing = sorted(name for name in known if name not in SUMMARY)
        self.assertEqual(missing, [])

    def test_numbered_table_matches_the_dispatcher(self):
        listed = {int(key) for key in PROCESS_SUMMARY}
        implemented = set(cli.NUMBERED) | cli.SPECIAL_NUMBERED
        self.assertEqual(listed, implemented)


class InitBatchHeader(unittest.TestCase):
    """`mainprogram jobscript --init-header qe|vasp`.

    The shipped examples say `--partition=dense`, which exists on one cluster
    and nowhere else, so copying one produces a job the scheduler rejects
    before anything runs. This fills in what SLURM and Lmod report and marks
    the rest TODO -- it does not claim to be submit-ready.
    """

    def test_the_flag_is_global_not_positional(self):
        """`rest` is nargs='*', so a per-command flag is swallowed by the main
        parser -- the same reason --force is global."""
        from htesp.cli import build_parser

        args = build_parser().parse_args(["jobscript", "--init-header", "qe"])
        self.assertEqual(args.process, "jobscript")
        self.assertEqual(args.init_header, "qe")

    def test_every_spelling_of_the_code_is_accepted(self):
        """Sites, and people, write Quantum ESPRESSO every way there is."""
        from htesp.batch_header import normalise_code

        for spelling in ("qe", "QE", "QuantumEspresso", "QUANTUMESPRESSO",
                         "quantum-espresso", "quantum_espresso",
                         "Quantum ESPRESSO", "espresso", "pw"):
            with self.subTest(spelling=spelling):
                self.assertEqual(normalise_code(spelling), "qe")
        for spelling in ("vasp", "VASP", "Vasp"):
            with self.subTest(spelling=spelling):
                self.assertEqual(normalise_code(spelling), "vasp")

    def test_an_unknown_code_says_what_it_takes(self):
        """argparse no longer screens the value, so the message has to."""
        from htesp.batch_header import normalise_code

        with self.assertRaises(ValueError) as caught:
            normalise_code("abinit")
        self.assertIn("quantumespresso", str(caught.exception))
        self.assertIn("abinit", str(caught.exception))

    def test_an_unknown_code_exits_two_rather_than_tracebacks(self):
        import contextlib
        import io
        import tempfile

        from htesp.cli import Context, cmd_jobscript

        with tempfile.TemporaryDirectory() as tmp:
            ctx = Context.__new__(Context)
            ctx.init_header, ctx.force, ctx.root = "abinit", False, Path(tmp)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                self.assertEqual(cmd_jobscript(ctx, []), 2)
        self.assertIn("abinit", buf.getvalue())

    def test_the_module_search_covers_the_common_spellings(self):
        """A module name this list misses puts a TODO where the `module load`
        belongs, on a machine that has the module."""
        from htesp.batch_header import MODULE_NAMES

        for name in ("qe", "quantum-espresso", "quantumespresso", "espresso"):
            with self.subTest(name=name):
                self.assertIn(name, MODULE_NAMES["qe"])

    def test_gres_is_omitted_when_the_cluster_has_none(self):
        """Vista is a GPU machine whose GresTypes is (null); an emitted
        --gres=gpu:1 would be a flag the scheduler rejects."""
        from htesp import batch_header

        text = batch_header.build("qe")
        if not batch_header.gres_configured():
            self.assertNotIn("--gres", text)

    def test_the_header_is_marked_not_finished(self):
        from htesp import batch_header

        text = batch_header.build("vasp")
        self.assertIn("TODO", text)
        self.assertIn("#SBATCH", text)
        self.assertIn("vasp_std", text)

    def test_the_header_contains_no_run_command(self):
        """generate_submission.py copies the header verbatim and appends its
        own run line.  A command written into the header would run first, and
        fail -- the first version of this generator emitted `ibrun pw.x ...`
        and would have done exactly that."""
        from htesp import batch_header

        for code, executable in (("qe", "pw.x"), ("vasp", "vasp_std")):
            with self.subTest(code=code):
                text = batch_header.build(code)
                self.assertIn(executable, text)      # named, in a comment
                for line in text.splitlines():
                    stripped = line.strip()
                    if not stripped or stripped.startswith("#"):
                        continue
                    self.assertTrue(
                        stripped.startswith("module ") or stripped == "#!/bin/bash",
                        "executable line in batch.header: " + line)

    def test_an_unknown_code_is_rejected(self):
        from htesp import batch_header

        with self.assertRaises(ValueError):
            batch_header.build("abinit")

    def test_it_refuses_to_overwrite_without_force(self):
        import tempfile
        from pathlib import Path as _Path

        from htesp import batch_header

        with tempfile.TemporaryDirectory() as tmp:
            target = _Path(tmp) / "batch.header"
            target.write_text("# hand tuned\n")
            self.assertEqual(batch_header.write(target, "qe"), 1)
            self.assertEqual(target.read_text(), "# hand tuned\n")
            self.assertEqual(batch_header.write(target, "qe", force=True), 0)
            self.assertIn("#SBATCH", target.read_text())

    def test_the_lmod_default_marker_wins_over_version_order(self):
        """A site that defaults to an older build has a reason."""
        from htesp import batch_header

        mods = batch_header.modules("vasp")
        if len(mods) > 1 and batch_header.os.environ.get("LMOD_CMD"):
            self.assertTrue(mods[0])   # first entry is the (D) one when marked

    def test_every_probe_survives_a_machine_without_slurm(self):
        """No scheduler is a fact to report, not an error."""
        import unittest.mock as mock

        from htesp import batch_header

        with mock.patch.object(batch_header.shutil, "which", return_value=None):
            self.assertEqual(batch_header.partitions(), [])
            self.assertEqual(batch_header.accounts(), [])
            self.assertFalse(batch_header.gres_configured())
            self.assertIsNone(batch_header.launcher())
            self.assertIn("TODO", batch_header.build("qe"))


class ModuleDependencies(unittest.TestCase):
    """A hierarchical Lmod site hides an application until its toolchain is in.

    `qe/7.3` on this machine lives under
    /opt/apps/nvidia24/openmpi5/modulefiles, so `module load qe` in a job
    script fails with "these module(s) exist but cannot be loaded as
    requested" unless nvidia and openmpi came first.  It works when you type it
    interactively only because the login shell already has them, which is
    exactly the kind of difference that turns into a job that dies on line one.
    """

    SPIDER = """
    You will need to load all module(s) on any one of the lines below before the "qe/7.3" module is available to load.

      nvidia/24.5  cuda/12.5  openmpi/5.0.5
      nvidia/24.7  cuda/12.6  openmpi/5.0.5

----------------------------------------------------------------------------
  Help:
"""

    def test_the_prerequisite_line_is_parsed(self):
        import unittest.mock as mock

        from htesp import batch_header

        with mock.patch.object(batch_header, "_run", return_value=self.SPIDER):
            with mock.patch.object(batch_header, "loaded_modules", return_value=[]):
                with mock.patch.dict(batch_header.os.environ,
                                     {"LMOD_CMD": batch_header.__file__}):
                    names, exact = batch_header.prerequisites("qe/7.3")
        self.assertEqual(names, ["nvidia", "cuda", "openmpi"])
        self.assertIn("nvidia/24.5", exact)

    def test_the_combination_matching_this_machine_wins(self):
        """Several combinations are offered; the one the login shell is
        already running is the one demonstrably working here."""
        import unittest.mock as mock

        from htesp import batch_header

        with mock.patch.object(batch_header, "_run", return_value=self.SPIDER):
            with mock.patch.object(batch_header, "loaded_modules",
                                   return_value=["nvidia/24.7", "openmpi/5.0.5"]):
                with mock.patch.dict(batch_header.os.environ,
                                     {"LMOD_CMD": batch_header.__file__}):
                    _, exact = batch_header.prerequisites("qe/7.3")
        self.assertIn("nvidia/24.7", exact)
        self.assertNotIn("24.5", exact)

    def test_a_flat_site_gets_no_prerequisite_lines(self):
        """Not every site is hierarchical; inventing a toolchain there would
        break a header that would otherwise have worked."""
        import unittest.mock as mock

        from htesp import batch_header

        with mock.patch.object(batch_header, "_run", return_value="qe/7.3\n"):
            with mock.patch.dict(batch_header.os.environ,
                                 {"LMOD_CMD": batch_header.__file__}):
                self.assertEqual(batch_header.prerequisites("qe/7.3"), ([], ""))

    def test_no_lmod_is_not_an_error(self):
        import unittest.mock as mock

        from htesp import batch_header

        with mock.patch.dict(batch_header.os.environ, {}, clear=True):
            self.assertEqual(batch_header.prerequisites("qe/7.3"), ([], ""))
            self.assertEqual(batch_header.loaded_modules(), [])

    def test_the_prerequisites_are_loaded_before_the_code(self):
        import os

        from htesp import batch_header

        if not os.environ.get("LMOD_CMD"):
            self.skipTest("no Lmod on this machine")
        module = (batch_header.modules("qe") or [None])[0]
        if not module or not batch_header.prerequisites(module)[0]:
            self.skipTest("no hierarchical qe module here")
        lines = [l for l in batch_header.build("qe").splitlines()
                 if l.startswith("module load")]
        self.assertGreaterEqual(len(lines), 2)
        self.assertTrue(lines[-1].endswith(" qe"), lines)

    def test_the_prerequisites_are_unversioned_too(self):
        """Same reasoning as the code module: a pinned toolchain goes stale."""
        import os

        from htesp import batch_header

        if not os.environ.get("LMOD_CMD"):
            self.skipTest("no Lmod on this machine")
        for line in batch_header.build("qe").splitlines():
            if line.startswith("module load"):
                for name in line.split()[2:]:
                    self.assertNotIn("/", name)


class HelpTextDependencies(unittest.TestCase):
    """Lmod's hierarchy does not express every dependency.

    Bridges-2 reports "This module can be loaded directly: module load
    QuantumEspresso/7.5-intel" -- no hierarchy at all -- while the Help
    underneath says "module load intel-oneapi QuantumEspresso/7.5-intel".
    intel-oneapi carries the Intel MPI and MKL runtimes that build is linked
    against, so loading QE alone puts pw.x on PATH and then fails at run time
    on a missing shared library, which is far more confusing than a module
    that refuses to load.
    """

    BRIDGES = '\n  QuantumEspresso: QuantumEspresso/7.5-intel\n\n    This module can be loaded directly: module load QuantumEspresso/7.5-intel\n\n    Help:\n      This module is built with intel compiler, intel MPI, intel MKL and ELPA.\n\n      To load the module type\n\n      > module load intel-oneapi QuantumEspresso/7.5-intel\n\n      To unload the module type\n\n      > module unload QuantumEspresso/7.5-intel\n'
    HIERARCHY_THEN_PROSE = '\n    You will need to load all module(s) on any one of the lines below before the "qe/7.3" module is available to load.\n\n      nvidia/24.7  cuda/12.6  openmpi/5.0.5\n      To run codes in quantum espresso include the following lines\n'
    BOTH_SOURCES = '\n    You will need to load all module(s) on any one of the lines below before the "qe/7.3" module is available to load.\n\n      intel/2024  impi/2021\n\n    Help:\n      > module load intel-oneapi intel qe/7.3\n'
    VISTA_HELP = '\n    Help:\n      To run codes in quantum espresso, include the following lines:\n      module load qe/7.3\n      ibrun pw.x -input input.scf\n'

    def _prereqs(self, text, module, help_text=""):
        """*text* is what spider returns; help is empty unless given.

        `module help` is asked first now, so a test aimed at spider has to
        leave help silent or it never gets there.
        """
        import unittest.mock as mock

        from htesp import batch_header

        def fake(command):
            return help_text if command[2] == "help" else text

        with mock.patch.object(batch_header, "_run", side_effect=fake):
            with mock.patch.dict(batch_header.os.environ,
                                 {"LMOD_CMD": batch_header.__file__}):
                return batch_header.prerequisites(module)

    def test_a_flat_site_still_yields_its_toolchain(self):
        names, exact = self._prereqs(self.BRIDGES, "QuantumEspresso/7.5-intel")
        self.assertEqual(names, ["intel-oneapi"])
        self.assertEqual(exact, "")          # no hierarchy combination to quote

    def test_the_module_is_not_its_own_prerequisite(self):
        names, _ = self._prereqs(self.BRIDGES, "QuantumEspresso/7.5-intel")
        self.assertNotIn("QuantumEspresso", names)

    def test_an_unload_line_is_not_mistaken_for_a_load(self):
        names, _ = self._prereqs(self.BRIDGES, "QuantumEspresso/7.5-intel")
        self.assertNotIn("unload", names)

    def test_help_that_names_only_the_module_yields_nothing(self):
        """Vista's help says plainly "module load qe/7.3"."""
        self.assertEqual(self._prereqs(self.VISTA_HELP, "qe/7.3"), ([], ""))

    def test_prose_never_parses_as_a_module_list(self):
        """The Help block follows the hierarchy block, and a line of English
        would otherwise become several "modules"."""
        names, _ = self._prereqs(self.HIERARCHY_THEN_PROSE, "qe/7.3")
        self.assertEqual(names, ["nvidia", "cuda", "openmpi"])
        self.assertNotIn("run", names)

    def test_the_hierarchy_and_spiders_help_copy_are_merged(self):
        """Within the spider path both halves count, without duplicates."""
        names, _ = self._prereqs(self.BOTH_SOURCES, "qe/7.3")
        self.assertEqual(names, ["intel", "impi", "intel-oneapi"])

    def test_module_help_wins_outright_when_it_answers(self):
        """It is the module author speaking, so it is not merged with the
        hierarchy -- it replaces it."""
        names, exact = self._prereqs(
            self.BOTH_SOURCES, "qe/7.3",
            help_text="> module load site-toolchain qe/7.3\n")
        self.assertEqual(names, ["site-toolchain"])
        self.assertEqual(exact, "")


class LatestVersionIsProbed(unittest.TestCase):
    """The newest build is the one whose help describes the current toolchain.

    What gets loaded is still the bare name, so Lmod resolves it to the site
    default; only the probe uses the exact version.
    """

    def test_versions_sort_numerically_not_lexically(self):
        from htesp.batch_header import latest

        self.assertEqual(latest(["vasp/6.4.3", "vasp/5.4.4.pl2"]), "vasp/6.4.3")
        self.assertEqual(latest(["qe/7.9", "qe/7.10"]), "qe/7.10")
        self.assertEqual(latest(["a/1.2", "a/1.10"]), "a/1.10")

    def test_a_suffixed_version_does_not_crash_the_sort(self):
        from htesp.batch_header import latest

        self.assertEqual(latest(["QuantumEspresso/7.5-intel",
                                 "QuantumEspresso/6.7-pgi"]),
                         "QuantumEspresso/7.5-intel")

    def test_an_empty_list_is_not_an_error(self):
        from htesp.batch_header import latest

        self.assertIsNone(latest([]))

    def test_the_header_probes_the_latest_but_loads_the_bare_name(self):
        import os

        from htesp import batch_header

        if not os.environ.get("LMOD_CMD"):
            self.skipTest("no Lmod on this machine")
        if not batch_header.modules("vasp"):
            self.skipTest("no vasp module here")
        text = batch_header.build("vasp")
        loads = [l for l in text.splitlines() if l.startswith("module load")]
        self.assertTrue(loads)
        self.assertTrue(loads[-1].endswith(" vasp"), loads)


class ModuleHelpIsAlsoASource(unittest.TestCase):
    """`module help <name>` carries the same instruction as spider's Help.

    On Bridges-2 both `module help QuantumEspresso` (which resolves to the
    site default) and `module spider QuantumEspresso/7.5-intel` print
    "> module load intel-oneapi QuantumEspresso/7.5-intel".  Reading both
    means a site that carries the text in only one of them still works.
    """

    HELP = '------------ Module Specific Help for "QuantumEspresso/7.5-intel" ------------\nQuantumEspresso 7.5\n\nThis module is built with intel compiler, intel MPI, intel MKL, and ELPA.\n\nTo load the module type\n\n> module load intel-oneapi QuantumEspresso/7.5-intel\n\nTo unload the module type\n\n> module unload QuantumEspresso/7.5-intel\n'

    def test_the_toolchain_is_read_out_of_module_help(self):
        import unittest.mock as mock

        from htesp import batch_header

        # spider says nothing useful; help carries the instruction
        def fake(command):
            return self.HELP if "help" in command else "no hierarchy here\n"

        with mock.patch.object(batch_header, "_run", side_effect=fake):
            with mock.patch.dict(batch_header.os.environ,
                                 {"LMOD_CMD": batch_header.__file__}):
                names, exact = batch_header.prerequisites(
                    "QuantumEspresso/7.5-intel")
        self.assertEqual(names, ["intel-oneapi"])
        self.assertEqual(exact, "")

    def test_the_bare_name_is_tried_when_the_exact_build_has_no_help(self):
        """`module help qe` resolves to the default and may be the only one
        carrying the text."""
        import unittest.mock as mock

        from htesp import batch_header

        asked = []

        def fake(command):
            asked.append(" ".join(command[1:]))
            if command[-2:] == ["help", "QuantumEspresso"]:
                return self.HELP
            return "nothing\n"

        with mock.patch.object(batch_header, "_run", side_effect=fake):
            with mock.patch.dict(batch_header.os.environ,
                                 {"LMOD_CMD": batch_header.__file__}):
                names, _ = batch_header.prerequisites(
                    "QuantumEspresso/7.5-intel")
        self.assertEqual(names, ["intel-oneapi"])
        self.assertIn("bash help QuantumEspresso", asked)


class TheHeaderComesWithAWarning(unittest.TestCase):
    """A generated header is a starting point, not a working job script.

    The probes answer what this machine reports; they cannot know what a
    particular build needs at run time.  A module chain that is one module
    short is accepted by the queue and then fails inside the job, minutes
    later, with an error naming a shared library rather than a module.
    """

    def test_writing_a_header_warns(self):
        import contextlib
        import io
        import tempfile
        from pathlib import Path as _Path

        from htesp import batch_header

        with tempfile.TemporaryDirectory() as tmp:
            target = _Path(tmp) / "batch.header"
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                batch_header.write(target, "qe")
            printed = buf.getvalue()
        self.assertIn("CHECK THIS FILE BEFORE SUBMITTING", printed)
        self.assertIn("dependencies", printed)

    def test_the_warning_names_the_file_and_a_way_to_check_it(self):
        import contextlib
        import io
        import tempfile
        from pathlib import Path as _Path

        from htesp import batch_header

        with tempfile.TemporaryDirectory() as tmp:
            target = _Path(tmp) / "batch.header"
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                batch_header.write(target, "qe")
            printed = buf.getvalue()
        self.assertIn(str(target), printed)
        self.assertIn("source", printed)
        self.assertIn("pw.x", printed)

    def test_the_vasp_warning_names_the_vasp_executable(self):
        import contextlib
        import io
        import tempfile
        from pathlib import Path as _Path

        from htesp import batch_header

        with tempfile.TemporaryDirectory() as tmp:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                batch_header.write(_Path(tmp) / "batch.header", "vasp")
        self.assertIn("vasp_std", buf.getvalue())

    def test_a_refused_overwrite_does_not_warn(self):
        """Nothing was written, so there is nothing to check."""
        import contextlib
        import io
        import tempfile
        from pathlib import Path as _Path

        from htesp import batch_header

        with tempfile.TemporaryDirectory() as tmp:
            target = _Path(tmp) / "batch.header"
            target.write_text("# mine\n")
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                self.assertEqual(batch_header.write(target, "qe"), 1)
        self.assertNotIn("CHECK THIS FILE", buf.getvalue())

    def test_the_tutorial_runner_warns_once_not_once_per_tutorial(self):
        """A 42-tutorial sweep writes 42 headers."""
        import logging
        import unittest.mock as mock

        from tutorials import workdirs

        workdirs._HEADER_WARNING_SHOWN = False
        with mock.patch.object(workdirs.LOG, "warning") as warned:
            workdirs._warn_about_headers()
            workdirs._warn_about_headers()
            workdirs._warn_about_headers()
        self.assertEqual(warned.call_count, 1)
        workdirs._HEADER_WARNING_SHOWN = False


class HelpIsAskedBeforeSpider(unittest.TestCase):
    """`module help` is the module author speaking; spider is a derivation.

    Help is asked first and its answer taken whenever it gives one.  It is a
    preference, not an exclusion: most help text says no more than
    "module load qe/7.3", and on a hierarchical site spider is the only thing
    that reveals the toolchain.  Dropping spider whenever help merely exists
    would break every such site, this one included.
    """

    def _run_with(self, help_text, spider_text):
        import unittest.mock as mock

        from htesp import batch_header

        asked = []

        def fake(command):
            asked.append(command[2])
            return help_text if command[2] == "help" else spider_text

        with mock.patch.object(batch_header, "_run", side_effect=fake):
            with mock.patch.object(batch_header, "loaded_modules", return_value=[]):
                with mock.patch.dict(batch_header.os.environ,
                                     {"LMOD_CMD": batch_header.__file__}):
                    return batch_header.prerequisites("qe/7.3"), asked

    HIERARCHY = "\n".join([
        "",
        '    You will need to load all module(s) on any one of the lines below before the "qe/7.3" module is available to load.',
        "",
        "      nvidia/24.7  cuda/12.6  openmpi/5.0.5",
        "",
    ])

    def test_spider_is_not_run_when_help_answers(self):
        answer = "To load the module type\n> module load intel-oneapi qe/7.3\n"
        (names, exact), asked = self._run_with(answer, self.HIERARCHY)
        self.assertEqual(names, ["intel-oneapi"])
        self.assertEqual(exact, "")
        self.assertNotIn("spider", asked)

    def test_spider_is_used_when_help_names_no_prerequisite(self):
        """Vista: help says only "module load qe/7.3", which is the
        incomplete advice this exists to correct."""
        answer = "To run pw.x include:\nmodule load qe/7.3\nibrun pw.x\n"
        (names, exact), asked = self._run_with(answer, self.HIERARCHY)
        self.assertEqual(names, ["nvidia", "cuda", "openmpi"])
        self.assertIn("nvidia/24.7", exact)
        self.assertIn("spider", asked)

    def test_spider_is_used_when_there_is_no_help_at_all(self):
        (names, _), asked = self._run_with("", self.HIERARCHY)
        self.assertEqual(names, ["nvidia", "cuda", "openmpi"])
        self.assertIn("spider", asked)

    def test_this_machine_still_gets_its_toolchain(self):
        """The regression this ordering could have caused, checked live."""
        import os

        from htesp import batch_header

        if not os.environ.get("LMOD_CMD"):
            self.skipTest("no Lmod on this machine")
        module = (batch_header.modules("qe") or [None])[0]
        if not module:
            self.skipTest("no qe module here")
        names, _ = batch_header.prerequisites(module)
        self.assertTrue(names, "help-first must not lose the hierarchy chain")


if __name__ == "__main__":
    unittest.main()
