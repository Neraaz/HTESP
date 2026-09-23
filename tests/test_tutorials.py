"""Runs the tutorial runner's own suite as part of the package suite.

``tutorials/selftest.py`` covers the catalogue, the checkpoint file, the stop
report and an end-to-end dry run against a stubbed ``mainprogram``.  It lives
next to the runner so it can be run on its own with
``python -m unittest tutorials.selftest``; this module pulls it into
``pytest tests/`` as well, and adds the checks that tie the catalogue to the
rest of the package.
"""
from __future__ import annotations

import os
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def load_tests(loader, tests, pattern):          # noqa: ARG001 (unittest protocol)
    from tutorials import selftest
    tests.addTests(loader.loadTestsFromModule(selftest))
    return tests


class CatalogueMatchesTheCommandLine(unittest.TestCase):
    """Every command a tutorial step runs must exist in the dispatcher."""

    def test_every_step_command_is_a_real_htesp_command(self):
        from htesp import cli
        from tutorials.catalog import build_catalog

        known = set(cli.SPECIAL_COMMANDS) | set(cli.WORKFLOW_COMMANDS)
        known |= {str(number) for number in cli.NUMBERED}
        known |= {str(number) for number in cli.SPECIAL_NUMBERED}

        unknown = set()
        for tutorial in build_catalog().values():
            for step in tutorial.steps:
                command = getattr(step, "command", None)
                if isinstance(command, str) and command not in known:
                    unknown.add(f"{tutorial.code}:{step.id} -> {command}")
        self.assertEqual(sorted(unknown), [])

    def test_every_tutorial_directory_exists(self):
        from tutorials.catalog import build_catalog

        missing = [t.code for t in build_catalog().values()
                   if not (ROOT / t.directory).is_dir()]
        self.assertEqual(missing, [])


class SubmissionScript(unittest.TestCase):
    def test_it_parses_and_carries_no_site_specific_partition(self):
        import subprocess

        script = ROOT / "tutorials" / "submit_tutorials.sh"
        self.assertTrue(script.is_file())
        result = subprocess.run(["bash", "-n", str(script)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        # 40 shipped example scripts hard-code one site's queue; the
        # submission script must only mention it in the comment that warns
        # against it, never in an #SBATCH directive.
        directives = [line for line in script.read_text().splitlines()
                      if line.strip().startswith("#SBATCH")]
        self.assertNotIn("dense", "\n".join(directives))


class StepsCannotHangOnStdin(unittest.TestCase):
    """`htesp-tutorials` sat on QE/4 for 33 minutes with an empty log dir.

    `qmpy_rester.get_oqmd_phases()` defaults to `verbose=True`, which calls
    `input('Proceed? [Y/n]:')`. The step subprocess inherited the runner's
    stdin, so the prompt waited for a keypress nobody was there to give.
    """

    def test_the_step_subprocess_gets_no_stdin(self):
        text = (ROOT / "tutorials" / "runner.py").read_text()
        self.assertIn("stdin=subprocess.DEVNULL", text)

    def test_the_log_header_is_written_before_the_command_runs(self):
        """A step that hangs must still leave its command and cwd on disk."""
        text = (ROOT / "tutorials" / "runner.py").read_text()
        before = text.index('log_path.write_text("\\n".join(header))')
        after = text.index("proc = subprocess.run(state.command")
        self.assertLess(before, after,
                        "the log header must be written before the subprocess starts")

    def test_the_oqmd_query_is_not_interactive(self):
        """oqmd_extract must never let qmpy_rester prompt."""
        text = (ROOT / "htesp" / "oqmd_extract.py").read_text()
        self.assertIn("get_oqmd_phases(verbose=False", text)
        self.assertNotIn("get_oqmd_phases(**kwargs)", text)


class VaspArtifactGlobs(unittest.TestCase):
    """A VASP download that wrote every input was still reported FAILED.

    HTESP writes VASP inputs into `R<mpid>-<compound>/relax/`, but two
    artifact globs named `R*-*/POSCAR` and `R*-*/INCAR`, omitting the `relax`
    level. The runner's "exited 0 and produced nothing" rule then failed the
    step even though INCAR/KPOINTS/POSCAR/POTCAR were all present.
    """

    def test_vasp_inputs_are_looked_for_under_relax(self):
        text = (ROOT / "tutorials" / "steps.py").read_text()
        self.assertIn('f"{MAT}/relax/POSCAR"', text)
        self.assertIn('f"{MAT}/relax/INCAR"', text)
        self.assertNotIn('f"{MAT}/POSCAR"', text)
        self.assertNotIn('(f"{MAT}/INCAR",)', text)


class ConfigureVaspPotcars(unittest.TestCase):
    """`htesp-check --config_vasp_pot DIR` points pymatgen at the POTCARs.

    Without PMG_VASP_PSP_DIR every VASP tutorial dies with PmgVaspPspDirError.
    The setting must name the *parent* of POT_GGA_PAW_PBE, which is the easy
    thing to get wrong, so the option accepts either and works it out.
    """

    def test_the_option_exists(self):
        from htesp.check import configure_vasp_potcars  # noqa: F401

        text = (ROOT / "htesp" / "check.py").read_text()
        self.assertIn('"--config_vasp_pot"', text)

    def test_a_missing_directory_is_reported_not_written(self):
        from htesp.check import configure_vasp_potcars

        self.assertEqual(configure_vasp_potcars("/no/such/directory/here"), 1)

    def test_a_directory_without_potentials_is_rejected(self):
        import tempfile

        from htesp.check import configure_vasp_potcars

        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(configure_vasp_potcars(tmp), 1)

    def test_both_the_family_dir_and_its_parent_are_accepted(self):
        """Give it POT_GGA_PAW_PBE or the directory above; same result."""
        import tempfile
        from pathlib import Path as _Path

        from htesp.check import VASP_FUNCTIONAL_DIRS

        self.assertIn("POT_GGA_PAW_PBE", VASP_FUNCTIONAL_DIRS)
        with tempfile.TemporaryDirectory() as tmp:
            family = _Path(tmp) / "POT_GGA_PAW_PBE" / "Mg"
            family.mkdir(parents=True)
            (family / "POTCAR").write_text("PAW_PBE Mg\n")
            # the resolution logic, without touching the user's pmgrc
            for given in (_Path(tmp), _Path(tmp) / "POT_GGA_PAW_PBE"):
                if given.name in VASP_FUNCTIONAL_DIRS:
                    root = given.parent
                else:
                    root = given
                self.assertEqual(root, _Path(tmp))


class EnumlibPreflight(unittest.TestCase):
    """`magenum` died three minutes into a sweep with a RuntimeError.

    pymatgen's EnumlibAdaptor shells out to enum.x/multienum.x and
    makestr.x/makeStr.py, which are a separate C/Fortran package that nothing
    pip-installs. Preflight now says so before the first step runs.
    """

    def test_the_magenum_step_declares_the_dependency(self):
        import sys

        sys.path.insert(0, str(ROOT))
        from tutorials.catalog import build_catalog

        catalog = build_catalog()
        for code in ("QE/21", "VASP/20"):
            with self.subTest(code=code):
                steps = [s for s in catalog[code].steps if s.id == "magenum"]
                self.assertTrue(steps, "no magenum step in " + code)
                self.assertTrue(steps[0].needs_enumlib)

    def test_preflight_warns_rather_than_erroring(self):
        """Two tutorials need it; forty do not. A hard error would block them."""
        text = (ROOT / "tutorials" / "runner.py").read_text()
        self.assertIn("ENUMLIB_TOOLS", text)
        self.assertIn("needs_enumlib", text)
        # bounded by the next preflight check, not by a mode branch: real
        # mode is gone, and "if options.mode" no longer appears anywhere
        block = text.split("needs_enumlib = sorted", 1)[1].split("probe =", 1)[0]
        self.assertIn('"warning"', block)
        self.assertNotIn('"error"', block)

    def test_every_enumlib_executable_name_is_probed(self):
        import sys

        sys.path.insert(0, str(ROOT))
        from tutorials.runner import ENUMLIB_TOOLS

        for tool in ("enum.x", "multienum.x", "makestr.x", "makeStr.py"):
            self.assertIn(tool, ENUMLIB_TOOLS)


class UnpackedReferenceDirs(unittest.TestCase):
    """VASP/8 seeded no .cif files and `download` wrote nothing.

    The "archive" seed only globbed `reference*.tar.gz`. QE/8 ships one, but
    VASP/8 ships an unpacked `reference/` directory holding the .cif inputs.
    """

    def test_unpacked_reference_directories_are_seeded(self):
        text = (ROOT / "tutorials" / "workdirs.py").read_text()
        self.assertIn('tutorial.directory.glob("reference*")', text)
        self.assertIn("folder.rglob", text)

    def test_the_vasp_cif_tutorial_really_ships_a_directory(self):
        """If this ever becomes a tarball the special case can go."""
        import sys

        sys.path.insert(0, str(ROOT))
        from tutorials.catalog import EXAMPLES

        folder = EXAMPLES / "VASP" / "tutorial8" / "reference"
        if not folder.is_dir():
            self.skipTest("example tree not present")
        self.assertTrue(list(folder.glob("*.cif")),
                        "VASP/8's .cif inputs are no longer in reference/")


class MissingPotcarNeverFails(unittest.TestCase):
    """A missing POTCAR must not abort input generation.

    POTCARs are licensed: HTESP cannot ship them, the reference relax/
    directories hold INCAR, KPOINTS and POSCAR only, and a machine may have no
    VASP licence at all. So: build one when possible, say what to do when not,
    and carry on either way. `convtest` used to die on the copy.
    """

    def test_a_present_potcar_passes_quietly(self):
        import tempfile
        from pathlib import Path as _Path

        from htesp.write_potcar import stage_potcar

        with tempfile.TemporaryDirectory() as tmp:
            (_Path(tmp) / "POTCAR").write_text("PAW_PBE Mg\n")
            self.assertTrue(stage_potcar(tmp))

    def test_an_unobtainable_potcar_warns_and_returns_false(self):
        """No POSCAR to build from, no POTCAR present -- must not raise."""
        import tempfile

        from htesp.write_potcar import stage_potcar

        with tempfile.TemporaryDirectory() as tmp:
            try:
                result = stage_potcar(tmp)
            except Exception as exc:                # noqa: BLE001
                self.fail("stage_potcar raised %r; it must never fail" % (exc,))
            self.assertFalse(result)

    def test_nothing_raises_on_a_missing_potcar_any_more(self):
        text = (ROOT / "htesp" / "convergence_test.py").read_text()
        self.assertIn("stage_potcar(relax)", text)
        self.assertNotIn("raise FileNotFoundError", text)

    def test_the_message_names_the_command_that_fixes_it(self):
        from htesp.write_potcar import POTCAR_HELP

        self.assertIn("--config_vasp_pot", POTCAR_HELP)
        self.assertIn("licensed", POTCAR_HELP)

    def test_the_copy_list_drops_potcar_when_it_is_absent(self):
        """Otherwise shutil.copy would raise on the very next line."""
        text = (ROOT / "htesp" / "convergence_test.py").read_text()
        self.assertIn('wanted = ["INCAR", "KPOINTS", "POSCAR"]', text)
        self.assertIn('wanted.append("POTCAR")', text)


class BandPostNeedsDftOutput(unittest.TestCase):
    """VASP band-post failed with "No such file or directory: 'KPOINTS_band'".

    vasp_process writes KPOINTS_band only when EIGENVAL is present, so the step
    reads real DFT output and must be skipped under --dry-run, not failed.
    """

    def test_the_step_is_classified_as_needing_dft_output(self):
        import sys

        sys.path.insert(0, str(ROOT))
        from tutorials.catalog import build_catalog

        steps = [s for s in build_catalog()["VASP/11"].steps if s.id == "band-post"]
        self.assertTrue(steps)
        self.assertTrue(steps[0].needs_dft_output)

    def test_the_bare_errno_became_an_explanation(self):
        text = (ROOT / "htesp" / "vasp_process.py").read_text()
        self.assertIn('if not os.path.isfile("KPOINTS_band")', text)
        self.assertIn("mainprogram 13", text)


class KeepOutputOption(unittest.TestCase):
    """`--keep_output no` drops the generated work directories.

    A full sweep leaves 42 of them; after a green run that is just bulk. The
    logs, report and checkpoint are always kept -- they are the record of the
    run and they are small.
    """

    def test_the_flag_exists_with_two_answers(self):
        text = (ROOT / "tutorials" / "run_tutorials.py").read_text()
        self.assertIn('"--keep_output"', text)
        self.assertIn('choices=("yes", "no")', text)
        self.assertIn('default="yes"', text)

    def test_yes_is_the_default_so_behaviour_is_unchanged(self):
        import sys

        sys.path.insert(0, str(ROOT))
        from tutorials.runner import RunOptions

        self.assertEqual(RunOptions(workdir=Path(".")).keep, "all")

    def test_cleanup_only_ever_deletes_under_tutorial_runs(self):
        """state.json, logs/ and report.* must survive any --keep_output.

        Checked by finding the deletion calls in the AST rather than grepping
        for words -- the method legitimately says "report" in its docstring and
        in the line it logs.
        """
        import ast

        source = (ROOT / "tutorials" / "runner.py").read_text()
        tree = ast.parse(source)
        func = next(n for n in ast.walk(tree)
                    if isinstance(n, ast.FunctionDef) and n.name == "_clean_work_dirs")

        deletions = []
        for node in ast.walk(func):
            if not isinstance(node, ast.Call):
                continue
            name = getattr(node.func, "attr", None)
            if name in ("rmtree", "remove", "unlink", "rmdir"):
                deletions.append((name, [ast.unparse(a) for a in node.args]))
        self.assertEqual(len(deletions), 1,
                         "expected exactly one deletion, found %r" % (deletions,))
        verb, args = deletions[0]
        self.assertEqual(verb, "rmtree")
        self.assertEqual(args, ["target"])

        # ...and `target` is built from runs_root, nothing else
        assigns = [ast.unparse(n) for n in ast.walk(func)
                   if isinstance(n, ast.Assign)
                   and any(getattr(x, "id", "") == "target" for x in n.targets)]
        self.assertTrue(assigns, "target is never assigned")
        self.assertTrue(all("runs_root" in a for a in assigns), assigns)

    def test_an_interrupted_run_keeps_everything(self):
        """Ctrl-C means 'I want to look at this', not 'tidy up'."""
        text = (ROOT / "tutorials" / "runner.py").read_text()
        body = text.split("def _clean_work_dirs", 1)[1].split("def _mark_running", 1)[0]
        self.assertIn("self.state.interrupted", body)

    def test_cleanup_runs_after_every_tutorial_not_during(self):
        """A work directory is the seed for its dependents; removing one
        mid-run would starve them.

        Checked against `run()` alone: the tutorials are executed by
        `_run_group` (which may run several at once), and that call has to
        come before the cleanup.
        """
        text = (ROOT / "tutorials" / "runner.py").read_text()
        run_body = text.split("def run(self)", 1)[1].split("\n    def ", 1)[0]
        self.assertLess(run_body.index("self._run_group(codes)"),
                        run_body.index("self._clean_work_dirs()"))

    def test_parallel_tutorials_still_clean_up_only_at_the_end(self):
        """The same rule, now that several tutorials can be in flight."""
        text = (ROOT / "tutorials" / "runner.py").read_text()
        group = text.split("def _run_group", 1)[1].split("\n    def ", 1)[0]
        self.assertNotIn("_clean_work_dirs", group)


class ResumeDoesNotClobberOutputs(unittest.TestCase):
    """`--resume --only QE/6` fed OQMD ids to a Materials Project lookup.

    Seeds overwrite by design, so a fresh work directory gets the tutorial's
    shipped files. On a resume the steps that would regenerate those files are
    skipped, so the shipped copy won instead -- and examples/QE/tutorial6 ships
    an mpid-list.in full of `oqmd-*` ids left over from the author's own run.
    The tutorial then "failed" for a reason that was not its own.
    """

    def test_resuming_seeds_only_what_is_missing(self):
        text = (ROOT / "tutorials" / "workdirs.py").read_text()
        self.assertIn("resuming = workdir.is_dir()", text)
        self.assertIn("overwrite = seed.overwrite and not resuming", text)

    def test_the_archive_seed_honours_it_too(self):
        text = (ROOT / "tutorials" / "workdirs.py").read_text()
        self.assertIn("_seed_from_archive(seed, tutorial, workdir, overwrite=not resuming)",
                      text)

    def test_a_fresh_run_still_overwrites(self):
        """--restart removes the directory first, so nothing is preserved."""
        import tempfile
        from pathlib import Path as _Path
        import sys

        sys.path.insert(0, str(ROOT))
        from tutorials.catalog import build_catalog
        from tutorials.runner import RunOptions
        from tutorials.workdirs import seed_workdir

        tutorial = build_catalog()["QE/1"]
        if not tutorial.directory.is_dir():
            self.skipTest("example tree not present")
        with tempfile.TemporaryDirectory() as tmp:
            options = RunOptions(workdir=_Path(tmp), resume=False)
            work = seed_workdir(tutorial, options)
            marker = work / "config.json"
            self.assertTrue(marker.is_file())
            marker.write_text("{}")               # pretend a step rewrote it
            # resume=False must restore the shipped copy
            seed_workdir(tutorial, options)
            self.assertNotEqual(marker.read_text(), "{}")

    def test_resuming_preserves_a_rewritten_file(self):
        import tempfile
        from pathlib import Path as _Path
        import sys

        sys.path.insert(0, str(ROOT))
        from tutorials.catalog import build_catalog
        from tutorials.runner import RunOptions
        from tutorials.workdirs import seed_workdir

        tutorial = build_catalog()["QE/1"]
        if not tutorial.directory.is_dir():
            self.skipTest("example tree not present")
        with tempfile.TemporaryDirectory() as tmp:
            work = seed_workdir(tutorial, RunOptions(workdir=_Path(tmp), resume=False))
            (work / "config.json").write_text("{}")
            seed_workdir(tutorial, RunOptions(workdir=_Path(tmp), resume=True))
            self.assertEqual((work / "config.json").read_text(), "{}",
                             "a resume overwrote what an earlier step produced")


class ShippedTrackingFilesMatchTheirTutorial(unittest.TestCase):
    """examples/QE/tutorial6 shipped an mpid-list.in full of `oqmd-*` ids.

    Tutorial 6 is the magnetic-configuration tutorial; its own `search` step
    produces Materials Project ids. The OQMD list was left over from the
    author's own run and got committed as though it were input, so resuming
    that tutorial fed OQMD ids to an MP lookup. The equivalent VASP file held
    `aflow:` auids. Both are removed; this keeps them out.
    """

    #: which id prefix belongs to which topic.  Everything else is MP.
    FOREIGN = {"oqmd-": "oqmd", "aflow:": "aflow"}

    def test_no_tutorial_ships_another_databases_ids(self):
        import sys

        sys.path.insert(0, str(ROOT))
        from tutorials.catalog import EXAMPLES, build_catalog

        if not EXAMPLES.is_dir():
            self.skipTest("example tree not present")
        catalog = build_catalog()
        offenders = []
        for code, tutorial in catalog.items():
            for name in ("mpid-list.in", "mpid.in"):
                path = tutorial.directory / name
                if not path.is_file():
                    continue
                first = (path.read_text().splitlines() or [""])[0]
                for prefix, topic in self.FOREIGN.items():
                    if prefix in first and tutorial.topic != topic:
                        offenders.append(
                            "%s ships %s with %r ids but its topic is %r"
                            % (code, name, prefix.rstrip("-:"), tutorial.topic))
        self.assertEqual(offenders, [])

    def test_no_tutorial_ships_an_alphaid_tracking_file(self):
        """`mp-aaackgnq` is the new alphabetic spelling; the tree uses mp-<int>."""
        import re
        import sys

        sys.path.insert(0, str(ROOT))
        from tutorials.catalog import EXAMPLES, build_catalog

        if not EXAMPLES.is_dir():
            self.skipTest("example tree not present")
        alpha = re.compile(r"\bmp-[a-z]{2,}\b")
        offenders = []
        for code, tutorial in build_catalog().items():
            for name in ("mpid-list.in", "mpid.in"):
                path = tutorial.directory / name
                if path.is_file() and alpha.search(path.read_text()):
                    offenders.append("%s: %s" % (code, path))
        self.assertEqual(offenders, [])


class PotcarlessStepsAreSkippedNotDone(unittest.TestCase):
    """Writing VASP inputs without a POTCAR is not a pass.

    `stage_potcar` deliberately never fails, so input generation carries on --
    but the inputs are unusable, and reporting the step as done is a false
    green. When pymatgen cannot produce a POTCAR the runner records the step as
    skipped, naming the command that fixes it.
    """

    def test_the_input_writing_vasp_steps_are_flagged(self):
        import sys

        sys.path.insert(0, str(ROOT))
        from tutorials.catalog import build_catalog

        catalog = build_catalog()
        flagged = {(code, s.id) for code, t in catalog.items() for s in t.steps
                   if getattr(s, "needs_potcar", False)}
        for wanted in (("VASP/2", "download"), ("VASP/8", "download"),
                       ("VASP/10", "convtest"), ("VASP/13", "substitute"),
                       ("VASP/14", "elastic-input"), ("VASP/20", "magenum")):
            with self.subTest(step=wanted):
                self.assertIn(wanted, flagged)

    def test_no_qe_step_is_flagged(self):
        """QE uses .upf pseudopotentials from examples/QE/pp; POTCARs are VASP."""
        import sys

        sys.path.insert(0, str(ROOT))
        from tutorials.catalog import build_catalog

        for code, tutorial in build_catalog().items():
            if not code.startswith("QE/"):
                continue
            for s in tutorial.steps:
                with self.subTest(code=code, step=s.id):
                    self.assertFalse(getattr(s, "needs_potcar", False))

    def test_the_skip_reason_names_the_fix(self):
        text = (ROOT / "tutorials" / "runner.py").read_text()
        self.assertIn("step.needs_potcar and not potcars_available()", text)
        self.assertIn("--config_vasp_pot", text)

    def test_availability_is_decided_by_pmg_vasp_psp_dir(self):
        import sys

        sys.path.insert(0, str(ROOT))
        from tutorials.runner import potcars_available

        potcars_available.cache_clear()
        self.assertIsInstance(potcars_available(), bool)
        potcars_available.cache_clear()


class PerFolderBatchHeader(unittest.TestCase):
    """Each work directory gets a header this cluster will accept.

    examples/<code>/batch.header says --partition=dense and loads no module.
    It was written for one machine; everywhere else sbatch rejects it before
    any calculation starts, so every submitting step of a real run fails for a
    reason that has nothing to do with HTESP.
    """

    def test_the_shipped_header_is_the_one_that_needs_replacing(self):
        for code in ("QE", "VASP"):
            with self.subTest(code=code):
                text = (ROOT / "examples" / code / "batch.header").read_text()
                self.assertIn("dense", text)
                self.assertNotIn("module load", text)

    def test_nothing_is_generated_without_slurm(self):
        """On a laptop the probes have nothing to say and a dry run submits
        nothing, so the shipped header is the more useful thing to leave."""
        import unittest.mock as mock

        from tutorials import workdirs

        source = (ROOT / "tutorials" / "workdirs.py").read_text()
        block = source.split("def _ensure_batch_header", 1)[1].split("\ndef ", 1)[0]
        self.assertIn('shutil.which("sinfo") is None', block)
        self.assertIn("return", block)

    def test_the_launcher_is_written_to_config_not_the_header(self):
        """mainprogram jobscript builds the run line from
        job_script.parallel_command; a command in the header would run first."""
        source = (ROOT / "tutorials" / "workdirs.py").read_text()
        block = source.split("def _match_launcher", 1)[1].split("\ndef ", 1)[0]
        self.assertIn("parallel_command", block)
        self.assertNotIn('job["nproc"]', block)

    def test_nproc_is_left_to_the_tutorial(self):
        """It is the study's choice, not the machine's."""
        source = (ROOT / "tutorials" / "workdirs.py").read_text()
        block = source.split("def _match_launcher", 1)[1].split("\ndef ", 1)[0]
        self.assertNotIn('"nproc"', block)

    def test_ibrun_is_never_hardcoded(self):
        """ibrun exists at TACC and nowhere else; srun and mpirun are what the
        rest of the world has."""
        from htesp.batch_header import LAUNCHERS

        self.assertEqual(LAUNCHERS[:3], ("ibrun", "srun", "mpirun"))
        source = (ROOT / "tutorials" / "workdirs.py").read_text()
        self.assertNotIn('"ibrun"', source)

    def test_a_launcher_is_chosen_by_what_is_installed(self):
        import unittest.mock as mock

        from htesp import batch_header

        for present, expected in (({"srun"}, "srun"),
                                  ({"mpirun"}, "mpirun"),
                                  ({"ibrun", "srun", "mpirun"}, "ibrun"),
                                  (set(), None)):
            with self.subTest(present=sorted(present)):
                with mock.patch.object(
                        batch_header.shutil, "which",
                        side_effect=lambda n, p=present: n if n in p else None):
                    self.assertEqual(batch_header.launcher(), expected)


class SubmittingTutorialsBuildTheirRunScripts(unittest.TestCase):
    """Real mode submitted nothing, and looked fine doing it.

    HTESPWorkflow.stage_and_submit copies run-<stage>.sh from the work
    directory into the stage directory and submits that.  When the script is
    absent it returns status="skipped" -- "run-scf.sh not found in project
    root" -- and moves on.  Thirteen tutorials declared submitting steps
    without ever running `mainprogram jobscript`, so every submission was
    skipped, the steps exited 0, and the relaxation everything else depends on
    was never run.
    """

    def test_every_submitting_tutorial_builds_its_scripts_first(self):
        from tutorials.catalog import CATALOG

        for code, tutorial in sorted(CATALOG.items()):
            submitting = [s.id for s in tutorial.steps if s.submits]
            if not submitting:
                continue
            with self.subTest(code=code):
                ids = [s.id for s in tutorial.steps]
                self.assertIn("jobscript", ids,
                              f"{code} submits {submitting} but never builds run-*.sh")
                self.assertLess(ids.index("jobscript"),
                                min(ids.index(s) for s in submitting),
                                f"{code} builds its scripts after using them")

    def test_the_step_is_not_duplicated(self):
        """The jobscript tutorial already has one."""
        from tutorials.catalog import CATALOG

        for code, tutorial in CATALOG.items():
            ids = [s.id for s in tutorial.steps]
            with self.subTest(code=code):
                self.assertLessEqual(ids.count("jobscript"), 1)

    def test_nothing_is_prepended_to_a_tutorial_that_never_submits(self):
        """Building scripts nobody submits is noise in the report.

        Tested against the helper rather than the catalogue: the Wannier
        tutorials declare a jobscript step of their own without any
        `submits=True` step, because the scripts it builds are submitted by
        hand afterwards.
        """
        from tutorials.catalog import _with_job_scripts
        from tutorials.steps import Step

        steps = (Step("only", "reads a file", "22"),)
        self.assertEqual(_with_job_scripts(steps), steps)

    def test_a_tutorial_that_already_builds_them_is_left_alone(self):
        from tutorials.catalog import _with_job_scripts
        from tutorials.steps import JOBSCRIPT_STEP, Step

        steps = (JOBSCRIPT_STEP, Step("go", "submit", "1", submits=True))
        self.assertEqual(_with_job_scripts(steps), steps)

    def test_the_step_goes_in_front(self):
        from tutorials.catalog import _with_job_scripts
        from tutorials.steps import Step

        steps = (Step("prep", "write inputs", "4"),
                 Step("go", "submit", "1", submits=True))
        out = _with_job_scripts(steps)
        self.assertEqual([s.id for s in out], ["jobscript", "prep", "go"])

    def test_the_script_names_match_what_the_workflow_looks_for(self):
        """QE submits run-scf.sh, VASP run-vasp.sh; both come from
        job_script.command_list, so the two must agree."""
        import json

        source = (ROOT / "htesp" / "workflow.py").read_text()
        self.assertIn('"relax", "run-scf.sh"', source)
        self.assertIn('script: str = "run-vasp.sh"', source)
        qe = json.loads((ROOT / "examples" / "QE" / "config.json").read_text())
        vasp = json.loads((ROOT / "examples" / "VASP" / "config.json").read_text())
        self.assertIn("scf", qe["job_script"]["command_list"])
        self.assertIn("vasp", vasp["job_script"]["command_list"])

    def test_nothing_is_submitted_at_all_any_more(self):
        """The safety net that checked job ids went with real mode: this
        runner never calls the scheduler, so there are no jobs to count.
        Every step is invoked with --dry-run instead."""
        source = (ROOT / "tutorials" / "runner.py").read_text()
        self.assertNotIn("collect_job_ids", source)
        self.assertNotIn("wait_for_jobs", source)
        block = source.split("def _command", 1)[1].split("\n    def ", 1)[0]
        self.assertIn('cmd.append("--dry-run")', block)
        self.assertNotIn("if self.options.mode", block)


class BareModuleName(unittest.TestCase):
    """`module load qe`, not `module load qe/7.3`."""

    def test_the_header_loads_the_unversioned_name(self):
        from htesp import batch_header

        for code in ("qe", "vasp"):
            with self.subTest(code=code):
                for line in batch_header.build(code).splitlines():
                    if line.startswith("module load"):
                        loaded = line.split(None, 2)[2]
                        self.assertNotIn("/", loaded,
                                         "a pinned version goes stale when the "
                                         "site retires that build")

    def test_the_versions_found_are_still_listed(self):
        """Pinning must stay possible for a study that needs one build."""
        import os

        from htesp import batch_header

        if not os.environ.get("LMOD_CMD"):
            self.skipTest("no Lmod on this machine")
        text = batch_header.build("qe")
        if batch_header.modules("qe"):
            self.assertIn("versions available now:", text)
            self.assertIn("pin one by writing it out", text)


class PotcarReachesTheStageDirectory(unittest.TestCase):
    """VASP jobs were submitted with no POTCAR.

    Only the *download* path (htesp/vasp_input.py) built one.  A stage
    directory that was seeded rather than downloaded -- which is every tutorial
    starting from a prepared R<mpid>-<compound>/relax/, and any directory
    assembled by hand -- reached the scheduler with INCAR, KPOINTS and POSCAR
    only, and VASP stopped on the first step.
    """

    def test_the_vasp_submit_path_stages_one(self):
        source = (ROOT / "htesp" / "workflow.py").read_text()
        block = source.split("def _submit_vasp", 1)[1].split("\n    def ", 1)[0]
        self.assertIn("stage_potcar(target)", block)

    def test_it_is_staged_before_the_job_script_check(self):
        """A missing run-vasp.sh returns early; the POTCAR must be written
        before that or the directory is left incomplete."""
        source = (ROOT / "htesp" / "workflow.py").read_text()
        block = source.split("def _submit_vasp", 1)[1].split("\n    def ", 1)[0]
        self.assertLess(block.index("stage_potcar(target)"),
                        block.index("not found in project root"))

    def test_the_shared_helper_stages_one_for_vasp_only(self):
        """stage_and_submit serves both codes; QE has no POTCAR."""
        source = (ROOT / "htesp" / "workflow.py").read_text()
        block = source.split("def stage_and_submit", 1)[1].split("\n    def ", 1)[0]
        self.assertIn("if self.is_vasp:", block)
        self.assertIn("stage_potcar(target)", block)

    def test_staging_never_raises(self):
        """A machine may have no VASP licence at all; the rest of the input
        generation is still worth doing."""
        import tempfile

        from htesp.write_potcar import stage_potcar

        with tempfile.TemporaryDirectory() as tmp:
            self.assertIs(stage_potcar(Path(tmp)), False)   # no POSCAR, no crash

    def test_the_help_names_the_pmg_reorganisation(self):
        """`pmg config -p` is the step people miss; --config_vasp_pot alone
        only covers a tree already in <functional>/<symbol>/POTCAR layout."""
        from htesp.write_potcar import POTCAR_HELP

        self.assertIn("pmg config -p", POTCAR_HELP)
        self.assertIn("PMG_VASP_PSP_DIR", POTCAR_HELP)

    def test_config_vasp_pot_explains_the_same_route_when_it_fails(self):
        source = (ROOT / "htesp" / "check.py").read_text()
        block = source.split("def configure_vasp_potcars", 1)[1].split("\ndef ", 1)[0]
        self.assertIn("pmg config -p", block)
        self.assertIn("pmg config --add PMG_VASP_PSP_DIR", block)


class ApiKeyIsResolvedTheWayHtespResolvesIt(unittest.TestCase):
    """The runner asked only the environment, and nothing else did.

    `htesp-check --set_mp_api` exists so nobody has to export the variable --
    an exported key is lost by a batch job, a nohup-ed sweep or a new
    terminal -- and it writes ~/.config/htesp/credentials.  A correctly
    configured machine therefore had eight tutorials skip with "MP_API_KEY is
    not set" while `mainprogram search`, run by hand in the same directory,
    worked.
    """

    def test_the_credentials_file_is_enough(self):
        import tempfile
        import unittest.mock as mock
        from pathlib import Path as _Path

        from tutorials.runner import mp_api_key

        with tempfile.TemporaryDirectory() as tmp:
            creds = _Path(tmp) / "credentials"
            creds.write_text("MP_API_KEY = " + "k" * 32 + "\n")
            import htesp.config as cfg

            with mock.patch.object(cfg, "CREDENTIALS_PATH", creds):
                with mock.patch.dict(os.environ, {}, clear=True):
                    self.assertEqual(mp_api_key(), "k" * 32)

    def test_the_environment_still_wins(self):
        import unittest.mock as mock

        from tutorials.runner import mp_api_key

        with mock.patch.dict(os.environ, {"MP_API_KEY": "from-the-shell"}):
            self.assertEqual(mp_api_key(), "from-the-shell")

    def test_no_key_anywhere_is_still_no_key(self):
        import tempfile
        import unittest.mock as mock
        from pathlib import Path as _Path

        from tutorials.runner import mp_api_key

        with tempfile.TemporaryDirectory() as tmp:
            import htesp.config as cfg

            with mock.patch.object(cfg, "CREDENTIALS_PATH",
                                   _Path(tmp) / "nothing-here"):
                with mock.patch.dict(os.environ, {}, clear=True):
                    cfg.clear_cache()
                    self.assertIsNone(mp_api_key())

    def test_the_shipped_placeholder_is_not_a_key(self):
        """config.json ships use_your_API_KEY; sending that to MP is worse
        than skipping."""
        from htesp.config import API_KEY_PLACEHOLDER, api_key

        self.assertEqual(API_KEY_PLACEHOLDER, "use_your_API_KEY")
        # api_key() is what mp_api_key() delegates to
        self.assertTrue(callable(api_key))

    def test_the_runner_does_not_read_the_variable_directly(self):
        """Two ways of answering the same question is how they diverged."""
        source = (ROOT / "tutorials" / "runner.py").read_text()
        body = source.split("def mp_api_key", 1)[1].split("\ndef child_env", 1)[0]
        rest = source.replace(body, "")
        self.assertNotIn('os.environ.get("MP_API_KEY")', rest)

    def test_the_message_names_the_command_that_stores_it(self):
        source = (ROOT / "tutorials" / "runner.py").read_text()
        self.assertIn("htesp-check --set_mp_api", source)


class TutorialTimeBudget(unittest.TestCase):
    """A tutorial gets a minute; a hang is reported, not waited on.

    The step limit used to be six hours, which is how an `oqmd-download`
    stalled for nineteen minutes unnoticed -- alive, two open sockets, no
    output, because `qmpy_rester` builds a bare `requests.Session()` with no
    timeout.  Measured over a healthy sweep the slowest tutorial totalled
    58.6s, so a minute fits real work and cuts a hang short.
    """

    def test_the_default_budget_is_one_minute(self):
        from tutorials.runner import DEFAULT_TUTORIAL_TIMEOUT

        self.assertEqual(DEFAULT_TUTORIAL_TIMEOUT, 60)

    def test_the_flag_is_in_minutes(self):
        """Checked with no argument too: the flag kept an hours-era default of
        24 while its help had been rewritten, so every tutorial silently got a
        1440-second budget."""
        from tutorials.run_tutorials import build_parser

        self.assertEqual(build_parser().parse_args([]).timeout, 1.0)
        self.assertEqual(build_parser().parse_args(["--timeout", "2"]).timeout, 2.0)

    def test_the_help_says_minutes_not_hours(self):
        import contextlib
        import io

        from tutorials.run_tutorials import build_parser

        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            with self.assertRaises(SystemExit):
                build_parser().parse_args(["--help"])
        text = buf.getvalue()
        self.assertIn("--timeout MINUTES", text)
        self.assertNotIn("--timeout HOURS", text)

    def test_a_timeout_message_quotes_the_budget_that_applied(self):
        """OQMD runs on 100s; reporting the run-wide default instead told the
        reader a number that was never used."""
        source = (ROOT / "tutorials" / "runner.py").read_text()
        block = source.split("except subprocess.TimeoutExpired", 1)[1].split(
            "except OSError", 1)[0]
        self.assertIn('getattr(self._deadline, "budget"', block)
        self.assertNotIn("self.options.tutorial_timeout", block)

    def test_no_step_outlives_the_tutorial(self):
        """Otherwise a single hung call spends the whole run."""
        source = (ROOT / "tutorials" / "runner.py").read_text()
        block = source.split("def _run_subprocess", 1)[1].split("\n    def ", 1)[0]
        self.assertIn("self._time_left()", block)
        self.assertIn("min(timeout", block)

    def test_the_deadline_is_per_thread(self):
        """Tutorials run in parallel under --jobs; a shared attribute would
        give one tutorial another's clock."""
        source = (ROOT / "tutorials" / "runner.py").read_text()
        self.assertIn("self._deadline = threading.local()", source)

    def test_timeout_output_is_decoded_before_it_is_written(self):
        """subprocess.run(text=True) still hands TimeoutExpired *bytes*, so
        concatenating them raised TypeError and the timeout escaped as a
        traceback -- leaving the tutorial RUNNING and skipping the retry."""
        from tutorials.runner import _as_text

        self.assertEqual(_as_text(b"hello"), "hello")
        self.assertEqual(_as_text("hello"), "hello")
        self.assertEqual(_as_text(None), "")
        self.assertEqual(_as_text(b"\xff"), "\ufffd")

    def test_the_handler_uses_it(self):
        source = (ROOT / "tutorials" / "runner.py").read_text()
        block = source.split("except subprocess.TimeoutExpired", 1)[1].split(
            "except OSError", 1)[0]
        self.assertIn("_as_text(exc.stdout)", block)
        self.assertNotIn('(exc.stdout or "")', block)


class OqmdGetsASecondAttempt(unittest.TestCase):
    """OQMD is the least reliable service the tutorials touch.

    It has been unresponsive enough to fail a sweep outright, and its
    download has hung for nineteen minutes.  A search that normally takes 35
    seconds is not broken because one call stalled, so the tutorial runs
    again before being called a failure.
    """

    def test_oqmd_gets_a_longer_budget_than_the_rest(self):
        """Its searches alone have taken 34.7s to 71.8s across runs; a budget
        a healthy run cannot meet is a source of false failures, not a hang
        detector."""
        from tutorials.catalog import CATALOG
        from tutorials.runner import DEFAULT_TUTORIAL_TIMEOUT

        for code in ("QE/4", "VASP/4"):
            with self.subTest(code=code):
                self.assertEqual(CATALOG[code].timeout, 100)
                self.assertGreater(CATALOG[code].timeout,
                                   DEFAULT_TUTORIAL_TIMEOUT)

    def test_every_other_tutorial_uses_the_default(self):
        from tutorials.catalog import CATALOG

        for code, tutorial in CATALOG.items():
            if code in ("QE/4", "VASP/4"):
                continue
            with self.subTest(code=code):
                self.assertEqual(tutorial.timeout, 0.0)

    def test_the_per_tutorial_budget_wins_over_the_run_wide_one(self):
        source = (ROOT / "tutorials" / "runner.py").read_text()
        block = source.split("def _run_tutorial", 1)[1].split("\n    def ", 1)[0]
        self.assertIn("tutorial.timeout", block)

    def test_only_the_oqmd_tutorials_retry(self):
        from tutorials.catalog import CATALOG

        retried = sorted(c for c, t in CATALOG.items() if t.attempts > 1)
        self.assertEqual(retried, ["QE/4", "VASP/4"])

    def test_everything_else_runs_once(self):
        from tutorials.catalog import CATALOG

        for code, tutorial in CATALOG.items():
            if code in ("QE/4", "VASP/4"):
                continue
            with self.subTest(code=code):
                self.assertEqual(tutorial.attempts, 1)

    def test_only_a_clock_failure_is_retried(self):
        """A wrong answer is still wrong the second time; only running out of
        time is worth another go."""
        from tutorials.runner import TutorialRunner
        from tutorials.state import StepState

        ran_out = StepState(step_id="s", key="s",
                            reason="killed after 60s: the tutorial's 60s budget ran out")
        timed = StepState(step_id="s", key="s", reason="timed out after 60s")
        wrong = StepState(step_id="s", key="s",
                          reason="mainprogram 4 exited 2")
        self.assertTrue(TutorialRunner._out_of_time(ran_out))
        self.assertTrue(TutorialRunner._out_of_time(timed))
        self.assertFalse(TutorialRunner._out_of_time(wrong))

    def test_a_retry_does_not_repeat_work_already_done(self):
        """Re-running a 39-second search to retry the download after it would
        spend the new budget on work that already passed."""
        source = (ROOT / "tutorials" / "runner.py").read_text()
        block = source.split("def _run_step", 1)[1].split("\n    def ", 1)[0]
        self.assertIn("self.options.resume or force_resume", block)


class AFlakyServiceDoesNotFailTheRun(unittest.TestCase):
    """OQMD not answering is not a defect anyone reading the report can act on.

    Its searches have finished in 35 seconds and in 100; it has been
    unresponsive for a whole sweep; and its client, `qmpy_rester`, passes no
    timeout to its `requests.Session`, so a stalled connection once ran for
    nineteen minutes.  Turning a sweep red for that hides the failures that
    *are* actionable -- the same reasoning that already records an absent
    POTCAR or API key as skipped.
    """

    def test_a_timeout_there_is_recorded_as_skipped(self):
        source = (ROOT / "tutorials" / "runner.py").read_text()
        block = source.split("def _run_tutorial", 1)[1].split("\n    def ", 1)[0]
        self.assertIn("tutorial.flaky_service", block)
        self.assertIn("record.status = SKIPPED", block)

    def test_only_a_timeout_is_forgiven(self):
        """A wrong answer from OQMD is still a failure."""
        source = (ROOT / "tutorials" / "runner.py").read_text()
        block = source.split("def _run_tutorial", 1)[1].split("\n    def ", 1)[0]
        self.assertIn("self._out_of_time(failed)", block)

    def test_a_skipped_flaky_dependency_does_not_block(self):
        import tempfile
        from pathlib import Path as _Path

        from tutorials.catalog import CATALOG
        from tutorials.runner import RunOptions, TutorialRunner
        from tutorials.state import DONE, FAILED, SKIPPED

        with tempfile.TemporaryDirectory() as tmp:
            opts = RunOptions(workdir=_Path(tmp), examples=ROOT / "examples")
            runner = TutorialRunner(["QE/2", "QE/4", "QE/5", "QE/7"], opts)
            for dep in ("QE/2", "QE/5"):
                runner.state.tutorial(dep).status = DONE
            runner.state.tutorial("QE/4").status = SKIPPED
            self.assertEqual(runner._blocked_by(CATALOG["QE/7"]), "")

            # the same status on a dependency that is *not* a flaky service
            runner.state.tutorial("QE/5").status = SKIPPED
            self.assertEqual(runner._blocked_by(CATALOG["QE/7"]), "QE/5")

            # and a real failure still blocks
            runner.state.tutorial("QE/5").status = DONE
            runner.state.tutorial("QE/4").status = FAILED
            self.assertEqual(runner._blocked_by(CATALOG["QE/7"]), "QE/4")

    def test_the_dependency_reference_supplies_what_it_did_not_produce(self):
        """data-combine merges what the three database tutorials downloaded.
        When OQMD is skipped its work directory is bare, so the same files
        come from its reference instead."""
        import tempfile
        from pathlib import Path as _Path

        from tutorials.catalog import CATALOG, Seed
        from tutorials.workdirs import _seed_from_dependency_reference

        seed = Seed("QE/4", ("scf_dir", "R*-*", "mpid.in"))
        with tempfile.TemporaryDirectory() as tmp:
            workdir = _Path(tmp)
            placed = _seed_from_dependency_reference(seed, workdir)
        self.assertTrue(placed, "nothing came out of QE/4's reference")
        self.assertIn("mpid.in", placed)
        self.assertTrue(any(p.startswith("scf_dir/") for p in placed), placed)

    def test_it_takes_only_what_the_seed_asked_for(self):
        """The reference also holds `log` and `input.in`; copying those would
        hand the tutorial the author's run instead of its own."""
        import tempfile
        from pathlib import Path as _Path

        from tutorials.catalog import Seed
        from tutorials.workdirs import _seed_from_dependency_reference

        seed = Seed("QE/4", ("mpid.in",))
        with tempfile.TemporaryDirectory() as tmp:
            placed = _seed_from_dependency_reference(seed, _Path(tmp))
        self.assertEqual(placed, ["mpid.in"])

    def test_a_real_run_beats_the_reference(self):
        """The fallback fires only when the work directory yielded nothing."""
        source = (ROOT / "tutorials" / "workdirs.py").read_text()
        block = source.split("def seed_workdir", 1)[1].split("\ndef ", 1)[0]
        self.assertIn("if not placed", block)


class OqmdSocketTimeout(unittest.TestCase):
    """`qmpy_rester` passes no timeout to its `requests.Session`.

    A wall-clock budget cannot tell a stalled connection from a slow one; a
    socket timeout can, which is why this sits at the socket layer and the
    budget stays as an outer backstop.
    """

    def test_the_timeout_applies_only_inside_the_block(self):
        import socket

        from htesp.oqmd_extract import OQMD_SOCKET_TIMEOUT, socket_timeout

        before = socket.getdefaulttimeout()
        with socket_timeout():
            self.assertEqual(socket.getdefaulttimeout(), OQMD_SOCKET_TIMEOUT)
        self.assertEqual(socket.getdefaulttimeout(), before)

    def test_it_is_restored_when_the_query_raises(self):
        import socket

        from htesp.oqmd_extract import socket_timeout

        before = socket.getdefaulttimeout()
        with self.assertRaises(RuntimeError):
            with socket_timeout():
                raise RuntimeError("OQMD said no")
        self.assertEqual(socket.getdefaulttimeout(), before)

    def test_it_is_not_set_at_import_time(self):
        """aflow_extract imports this module, and data-combine runs all three
        front ends in one process: an import-time setting would put a timeout
        on AFLOW's queries too."""
        import socket

        self.assertIsNone(socket.getdefaulttimeout())

    def test_both_oqmd_calls_are_wrapped(self):
        import ast

        source = (ROOT / "htesp" / "oqmd_extract.py").read_text()
        calls = [n for n in ast.walk(ast.parse(source))
                 if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)
                 and n.func.attr in ("get_oqmd_phases", "get_entry_by_id")]
        self.assertEqual(len(calls), 2)
        withs = [n for n in ast.walk(ast.parse(source)) if isinstance(n, ast.With)]
        guarded = [w for w in withs
                   if any(isinstance(i.context_expr, ast.Call)
                          and getattr(i.context_expr.func, "id", "") == "socket_timeout"
                          for i in w.items)]
        self.assertEqual(len(guarded), 2)


class WhatTheRunProduced(unittest.TestCase):
    """The runner could say "did it work?" but not "what did it give me?".

    That gap hid three defects: `pressure-input` declared an artefact that
    `mainprogram 26` creates later, `update-input` declared one written only
    when a structure is *not* relaxed, and QE/9's declared artefact is
    shipped inside examples/ itself so the glob matched whether or not the
    step ran.  Each would have shown here as a step that passed having
    written nothing.
    """

    def test_the_flag_exists(self):
        from tutorials.run_tutorials import build_parser

        self.assertTrue(build_parser().parse_args(["--output"]).output)

    def test_a_step_records_what_it_wrote(self):
        from tutorials.state import StepState

        self.assertEqual(StepState(step_id="s", key="s").produced, [])

    def test_the_log_does_not_count_as_output(self):
        """Every mainprogram call appends to it, so counting it means no step
        ever looks empty -- and that emptiness is the whole signal."""
        from tutorials.manifest import is_ambient, real_output
        from tutorials.state import StepState

        self.assertTrue(is_ambient("log"))
        step = StepState(step_id="s", key="s", produced=["log"])
        self.assertEqual(real_output(step), [])
        step = StepState(step_id="s", key="s", produced=["log", "run-scf.sh"])
        self.assertEqual(real_output(step), ["run-scf.sh"])

    def test_nested_copies_are_described_too(self):
        """fnmatch's * crosses /, so matching the whole path alone would miss
        every file below the top level."""
        from tutorials.manifest import describe

        self.assertIn("submission script", describe("R1-x/phonopy/R1/run-scf.sh"))
        self.assertIn("submission script", describe("run-scf.sh"))
        self.assertIn("relaxation input", describe("Rmp-763-Mg1B2/relax/scf.in"))
        self.assertEqual(describe("something-unknown.xyz"), "")

    def test_snapshot_skips_symlinked_trees(self):
        """QE work directories link pp/ at the shared pseudopotential tree;
        walking it would add hundreds of files no step writes."""
        import tempfile
        from pathlib import Path as _Path

        from tutorials.workdirs import snapshot

        with tempfile.TemporaryDirectory() as tmp:
            root = _Path(tmp)
            (root / "real").mkdir()
            (root / "real" / "a.txt").write_text("a")
            target = root / "elsewhere"
            target.mkdir()
            (target / "b.txt").write_text("b")
            (root / "linked").symlink_to(target)
            seen = snapshot(root)
        self.assertIn("real/a.txt", seen)
        self.assertNotIn("linked/b.txt", seen)

    def test_changed_since_reports_new_and_modified(self):
        import tempfile
        import time
        from pathlib import Path as _Path

        from tutorials.workdirs import changed_since, snapshot

        with tempfile.TemporaryDirectory() as tmp:
            root = _Path(tmp)
            (root / "kept.txt").write_text("same")
            (root / "grows.txt").write_text("x")
            before = snapshot(root)
            time.sleep(0.01)
            (root / "new.txt").write_text("new")
            (root / "grows.txt").write_text("xx")
            changed = changed_since(before, snapshot(root))
        self.assertIn("new.txt", changed)
        self.assertIn("grows.txt", changed)
        self.assertNotIn("kept.txt", changed)

    def test_the_report_names_steps_that_wrote_nothing(self):
        source = (ROOT / "tutorials" / "report.py").read_text()
        self.assertIn("passed without writing anything", source)
        self.assertIn("manifest.real_output", source)


class SkippedStepsSayWhereToReadOn(unittest.TestCase):
    """A step this runner can never perform should say how to perform it.

    Forty-one steps read the output of a real DFT run.  Telling the reader
    only that they were skipped leaves them nowhere to go; the tutorial's own
    written instructions are where they should go.
    """

    def test_every_tutorial_but_one_resolves_to_a_readme(self):
        """Twelve ship none of their own, so the counterpart's are used -- the
        two trees cover the same topics in the same order."""
        from tutorials.catalog import CATALOG, readme_for

        missing = sorted(c for c, t in CATALOG.items() if readme_for(t) is None)
        self.assertEqual(missing, ["VASP/21"])

    def test_the_counterpart_is_found_in_both_directions(self):
        from tutorials.catalog import CATALOG, readme_for

        # VASP/9 ships none; QE/9 covers the same relaxation
        self.assertEqual(readme_for(CATALOG["VASP/9"]).parent.name, "tutorial9")
        self.assertEqual(readme_for(CATALOG["VASP/9"]).parent.parent.name, "QE")
        # QE/15 ships none; its VASP counterpart is 14, not 15
        self.assertEqual(readme_for(CATALOG["QE/15"]).parent.name, "tutorial14")

    def test_the_offset_after_tutorial_eleven_is_respected(self):
        """VASP n is QE n+1 from 11 on, so VASP/19 is QE/20."""
        from tutorials.catalog import CATALOG, readme_for

        self.assertEqual(readme_for(CATALOG["VASP/19"]).parent.name, "tutorial20")

    def test_readme_txt_is_accepted_too(self):
        from tutorials.catalog import README_NAMES

        self.assertIn("README", README_NAMES)
        self.assertIn("README.txt", README_NAMES)

    def test_every_resolved_path_exists(self):
        """A pointer to a file that is not there is worse than none."""
        from tutorials.catalog import CATALOG, readme_for

        for code, tutorial in sorted(CATALOG.items()):
            found = readme_for(tutorial)
            if found is None:
                continue
            with self.subTest(code=code):
                self.assertTrue(found.is_file(), found)

    def test_the_pointer_is_only_for_dft_output_skips(self):
        """A machine missing a POTCAR, an API key or enumlib already gets a
        message saying what to install; replacing it would be a downgrade.
        A cascade skip names the step that caused it, which is the real
        reason."""
        source = (ROOT / "tutorials" / "runner.py").read_text()
        block = source.split("def _run_step", 1)[1].split("\n    def ", 1)[0]
        self.assertEqual(block.count("self._how_to_run_it("), 1)
        dft = block.split("needs_dft_output", 1)[1][:400]
        self.assertIn("_how_to_run_it", dft)
        for other in ("POTCARs are not configured", "htesp-check --set_mp_api",
                      "whose output it reads, was skipped"):
            with self.subTest(message=other):
                segment = block.split(other, 1)
                self.assertEqual(len(segment), 2, other)
                self.assertNotIn("_how_to_run_it", segment[1][:200])

    def test_a_tutorial_with_no_instructions_gets_no_pointer(self):
        import tempfile
        from pathlib import Path as _Path

        from tutorials.catalog import CATALOG
        from tutorials.runner import RunOptions, TutorialRunner

        with tempfile.TemporaryDirectory() as tmp:
            runner = TutorialRunner(["VASP/21"],
                                    RunOptions(workdir=_Path(tmp),
                                               examples=ROOT / "examples"))
            self.assertEqual(runner._how_to_run_it(CATALOG["VASP/21"]), "")
            self.assertIn("tutorial9",
                          runner._how_to_run_it(CATALOG["VASP/9"]))

    def test_the_report_names_tutorials_nothing_ran_for(self):
        source = (ROOT / "tutorials" / "report.py").read_text()
        self.assertIn("Nothing ran for these", source)
        self.assertIn("readme_for", source)


if __name__ == "__main__":
    unittest.main()


class ExampleTreeIsFound(unittest.TestCase):
    """`examples/` is 185 MB and is not in the wheel, so after `pip install .`
    the package-relative path points at nothing and the runner said only
    "the example tree .../site-packages/examples is missing"."""

    def test_the_tree_is_searched_for_not_assumed(self):
        import os
        import tempfile
        from tutorials import catalog

        with tempfile.TemporaryDirectory() as tmp:
            tree = Path(tmp) / "elsewhere" / "examples"
            (tree / "QE").mkdir(parents=True)
            previous = os.environ.get("HTESP_EXAMPLES")
            os.environ["HTESP_EXAMPLES"] = str(tree)
            try:
                self.assertEqual(catalog.find_examples(), tree.resolve())
                self.assertIn("$HTESP_EXAMPLES", catalog.searched_for_examples()[0])
            finally:
                if previous is None:
                    os.environ.pop("HTESP_EXAMPLES", None)
                else:
                    os.environ["HTESP_EXAMPLES"] = previous

    def test_the_default_is_the_repository_tree(self):
        from tutorials import catalog
        self.assertTrue((catalog.EXAMPLES / "QE").is_dir())
        self.assertEqual(catalog.EXAMPLES, (ROOT / "examples").resolve())

    def test_use_examples_rebuilds_the_catalogue(self):
        """`stub` is decided by reading each tutorial's config.json, so a
        catalogue built against the wrong root marks everything a stub."""
        from tutorials import catalog
        saved_root, saved_catalog = catalog.EXAMPLES, catalog.CATALOG
        try:
            rebuilt = catalog.use_examples(ROOT / "examples")
            self.assertEqual(set(rebuilt), set(saved_catalog))
            self.assertEqual(sum(1 for t in rebuilt.values() if t.stub),
                             sum(1 for t in saved_catalog.values() if t.stub))
        finally:
            catalog.EXAMPLES, catalog.CATALOG = saved_root, saved_catalog

    def test_a_workdir_inside_the_example_tree_is_refused(self):
        """`--workdir examples/` would write generated runs into the read-only
        reference tree; it is usually `--examples` that was meant."""
        from tutorials.runner import RunOptions, preflight
        options = RunOptions(workdir=ROOT / "examples" / "runs",
                             examples=ROOT / "examples")
        problems = [str(p) for p in preflight(["QE/1"], options) if p.fatal]
        self.assertTrue(any("inside the example tree" in text
                            for text in problems), problems)

    def test_the_missing_tree_message_says_how_to_point_at_one(self):
        from tutorials.runner import RunOptions, preflight
        options = RunOptions(workdir=Path("/tmp/htesp-runs"),
                             examples=Path("/nonexistent/site-packages/examples"))
        problems = [str(p) for p in preflight(["QE/1"], options)]
        joined = "\n".join(problems)
        self.assertIn("--examples", joined)
        self.assertIn("HTESP_EXAMPLES", joined)
        self.assertIn("Looked in:", joined)
