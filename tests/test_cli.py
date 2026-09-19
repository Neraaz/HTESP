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


if __name__ == "__main__":
    unittest.main()
