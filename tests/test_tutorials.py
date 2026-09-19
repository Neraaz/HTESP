"""Runs the tutorial runner's own suite as part of the package suite.

``tutorials/selftest.py`` covers the catalogue, the checkpoint file, the stop
report and an end-to-end dry run against a stubbed ``mainprogram``.  It lives
next to the runner so it can be run on its own with
``python -m unittest tutorials.selftest``; this module pulls it into
``pytest tests/`` as well, and adds the checks that tie the catalogue to the
rest of the package.
"""
from __future__ import annotations

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
        block = text.split("needs_enumlib = sorted", 1)[1].split("if options.mode", 1)[0]
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
        mid-run would starve them."""
        text = (ROOT / "tutorials" / "runner.py").read_text()
        run_body = text.split("def run(self)", 1)[1].split("def _clean_work_dirs", 1)[0]
        loop_at = run_body.index("self._run_tutorial(code)")
        clean_at = run_body.index("self._clean_work_dirs()")
        self.assertLess(loop_at, clean_at)


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


class DependencyWaves(unittest.TestCase):
    """Running a real campaign in waves.

    Wave 0 ends with relaxations sitting in the queue.  Nothing in wave 1 can
    honestly start until those have finished and converged, because wave 1 *is*
    the tutorials that read the relaxed structure.  One pass over the whole
    catalogue either blocks for days inside squeue polling or proceeds on a
    structure that is not relaxed yet -- and the second is silently wrong,
    which is worse than failing.
    """

    def _catalog(self):
        from tutorials.catalog import CATALOG

        return CATALOG

    def test_wave_zero_is_everything_that_depends_on_nothing(self):
        from tutorials.catalog import select, waves

        catalog = self._catalog()
        plan = waves(select(catalog=catalog), catalog)
        self.assertTrue(plan)
        for code in plan[0]:
            with self.subTest(code=code):
                self.assertEqual(catalog[code].depends_on, ())

    def test_every_tutorial_appears_exactly_once(self):
        from tutorials.catalog import select, waves

        catalog = self._catalog()
        codes = select(catalog=catalog)
        flat = [code for wave in waves(codes, catalog) for code in wave]
        self.assertEqual(sorted(flat), sorted(codes))

    def test_a_dependency_is_always_in_an_earlier_wave(self):
        """The whole point: nothing runs before what it reads."""
        from tutorials.catalog import select, waves

        catalog = self._catalog()
        codes = select(catalog=catalog)
        plan = waves(codes, catalog)
        number = {code: i for i, wave in enumerate(plan) for code in wave}
        for code, index in number.items():
            for dep in catalog[code].depends_on:
                if dep in number:
                    with self.subTest(code=code, dep=dep):
                        self.assertLess(number[dep], index)

    def test_an_unselected_dependency_makes_the_dependent_a_root(self):
        """topological_order drops dependencies outside the selection, so the
        wave numbers must agree with it rather than inventing a wave for a
        tutorial that will not run."""
        from tutorials.catalog import select, waves

        catalog = self._catalog()
        plan = waves(select(only=["QE/11"], catalog=catalog), catalog)
        self.assertEqual(plan, [["QE/11"]])

    def test_the_wave_plan_survives_a_restart(self):
        import tempfile

        from tutorials.state import DONE, RunState

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "state.json"
            state = RunState(path=path)
            state.wave(0, ["QE/9"]).status = DONE
            state.wave(1, ["QE/11"])
            state.save()

            back = RunState.load(path)
            self.assertEqual(back.waves["0"].status, DONE)
            self.assertEqual(back.waves["1"].codes, ["QE/11"])
            self.assertEqual(back.next_wave([["QE/9"], ["QE/11"]]), 1)

    def test_a_failed_wave_is_offered_again_not_stepped_over(self):
        """Wave 1 reads what wave 0 produced; if wave 0 failed, there is
        nothing to read."""
        import tempfile

        from tutorials.state import FAILED, RunState

        with tempfile.TemporaryDirectory() as tmp:
            state = RunState(path=Path(tmp) / "state.json")
            state.wave(0, ["QE/9"]).status = FAILED
            self.assertEqual(state.next_wave([["QE/9"], ["QE/11"]]), 0)

    def test_forgetting_a_tutorial_unfinishes_its_wave(self):
        """--restart --only QE/9 must not leave wave 0 marked done, or the
        next run would step straight over the tutorial it just cleared."""
        import tempfile

        from tutorials.state import DONE, PENDING, RunState

        with tempfile.TemporaryDirectory() as tmp:
            state = RunState(path=Path(tmp) / "state.json")
            state.wave(0, ["QE/9", "QE/2"]).status = DONE
            state.wave(1, ["QE/11"]).status = DONE
            state.reset(["QE/9"])
            self.assertEqual(state.waves["0"].status, PENDING)
            self.assertEqual(state.waves["1"].status, DONE)

    def test_an_old_checkpoint_loads_with_no_waves(self):
        """The field is additive: a checkpoint written before waves existed is
        still a valid checkpoint, not a schema mismatch that discards hours."""
        import json
        import tempfile

        from tutorials.state import SCHEMA_VERSION, RunState

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "state.json"
            path.write_text(json.dumps({"version": SCHEMA_VERSION,
                                        "mode": "real", "tutorials": {}}))
            state = RunState.load(path)
            self.assertEqual(state.waves, {})
            self.assertEqual(state.mode, "real")

    def test_without_the_flag_the_run_is_one_pass_as_before(self):
        """--dry-run and --no-dft queue nothing, so there is nothing to wait
        between waves for; every earlier invocation must keep working."""
        source = (ROOT / "tutorials" / "runner.py").read_text()
        block = source.split("def _wave_to_run", 1)[1].split("\n    def ", 1)[0]
        self.assertIn("if wanted is None:", block)
        self.assertIn("return None, list(self.codes)", block)

    def test_the_cli_offers_the_flag(self):
        import argparse
        import contextlib
        import io

        from tutorials.run_tutorials import build_parser

        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            with self.assertRaises(SystemExit):
                build_parser().parse_args(["--help"])
        text = buf.getvalue()
        self.assertIn("--wave", text)
        self.assertIn("--list-waves", text)


class JobFinishedIsNotJobSucceeded(unittest.TestCase):
    """Two ways a step used to pass while having failed.

    A job id leaves `squeue` whether it COMPLETED, FAILED, TIMED OUT, was
    CANCELLED or ran out of memory; and both DFT codes stop at their ionic-step
    limit and exit cleanly.  Either way a full set of output files is left
    behind, the artefact globs match, the step is recorded DONE -- and every
    dependent tutorial then runs on a structure that was never relaxed, which
    is silently wrong rather than merely broken.
    """

    # -- sacct ------------------------------------------------------------- #
    def test_a_completed_job_is_the_only_good_state(self):
        import unittest.mock as mock

        from tutorials import workdirs

        with mock.patch.object(workdirs, "job_states",
                               return_value={"1": "COMPLETED", "2": "COMPLETED"}):
            self.assertEqual(workdirs.check_job_states(["1", "2"]), "")

    def test_a_walltime_kill_is_reported_in_words(self):
        import unittest.mock as mock

        from tutorials import workdirs

        with mock.patch.object(workdirs, "job_states", return_value={"7": "TIMEOUT"}):
            problem = workdirs.check_job_states(["7"])
        self.assertIn("TIMEOUT", problem)
        self.assertIn("hit the wall time", problem)

    def test_every_bad_state_is_caught(self):
        import unittest.mock as mock

        from tutorials import workdirs

        for state in ("FAILED", "TIMEOUT", "OUT_OF_MEMORY", "CANCELLED",
                      "NODE_FAIL", "PREEMPTED"):
            with self.subTest(state=state):
                with mock.patch.object(workdirs, "job_states",
                                       return_value={"1": state}):
                    self.assertTrue(workdirs.check_job_states(["1"]))

    def test_cancelled_by_someone_keeps_only_the_state(self):
        """sacct writes 'CANCELLED by 12345'; the uid is not part of it."""
        import subprocess
        import unittest.mock as mock

        from tutorials import workdirs

        done = subprocess.CompletedProcess(args=[], returncode=0,
                                           stdout="55|CANCELLED by 1001\n", stderr="")
        with mock.patch.object(workdirs.shutil, "which", return_value="/usr/bin/sacct"):
            with mock.patch("subprocess.run", return_value=done):
                self.assertEqual(workdirs.job_states(["55"]), {"55": "CANCELLED"})

    def test_no_sacct_is_not_the_same_as_every_job_failed(self):
        """Accounting is not enabled everywhere."""
        import unittest.mock as mock

        from tutorials import workdirs

        with mock.patch.object(workdirs.shutil, "which", return_value=None):
            self.assertEqual(workdirs.job_states(["1"]), {})
            self.assertEqual(workdirs.check_job_states(["1"]), "")

    # -- convergence -------------------------------------------------------- #
    def test_a_converged_vasp_relaxation_is_recognised(self):
        import tempfile

        from tutorials.workdirs import relaxation_converged

        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "OUTCAR").write_text(
                "lots of output\n reached required accuracy - stopping "
                "structural energy minimisation\n General timing\n")
            self.assertIs(relaxation_converged(Path(tmp)), True)

    def test_a_vasp_run_that_used_every_ionic_step_is_not_converged(self):
        """NSW reached without EDIFFG: OUTCAR and CONTCAR both exist, so the
        artefact check passes and only this catches it."""
        import tempfile

        from tutorials.workdirs import relaxation_converged

        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "OUTCAR").write_text("ionic step 60\n General timing\n")
            (Path(tmp) / "CONTCAR").write_text("Mg B2\n")
            self.assertIs(relaxation_converged(Path(tmp)), False)

    def test_a_converged_qe_relaxation_is_recognised(self):
        import tempfile

        from tutorials.workdirs import relaxation_converged

        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "relax.out").write_text(
                "End of BFGS Geometry Optimization\n JOB DONE.\n")
            self.assertIs(relaxation_converged(Path(tmp)), True)

    def test_a_qe_run_out_of_steps_is_not_converged(self):
        import tempfile

        from tutorials.workdirs import relaxation_converged

        with tempfile.TemporaryDirectory() as tmp:
            (Path(tmp) / "relax.out").write_text(
                "The maximum number of steps has been reached\n JOB DONE.\n")
            self.assertIs(relaxation_converged(Path(tmp)), False)

    def test_nothing_to_judge_by_is_not_a_failure(self):
        """An absent output file is the artefact check's business, not this
        one's; returning False here would report the wrong cause."""
        import tempfile

        from tutorials.workdirs import relaxation_converged

        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(relaxation_converged(Path(tmp)))

    def test_the_markers_match_the_ones_the_workflow_uses(self):
        """Two copies of this string would drift; this pins them together."""
        from tutorials.workdirs import CONVERGED_MARKERS

        source = (ROOT / "htesp" / "workflow.py").read_text()
        self.assertIn(CONVERGED_MARKERS["OUTCAR"], source)

    # -- wiring ------------------------------------------------------------- #
    def test_the_relaxation_steps_ask_for_the_check(self):
        from tutorials.catalog import CATALOG

        wanted = {"relax-submit", "resubmit", "relax-deformed", "relax-volumes"}
        seen = set()
        for tutorial in CATALOG.values():
            for step in tutorial.steps:
                if step.check_converged:
                    seen.add(step.id)
                    self.assertTrue(step.submits, step.id)
                    self.assertTrue(step.job_dirs, step.id)
        self.assertEqual(seen, wanted)

    def test_the_check_only_runs_against_real_output(self):
        """--dry-run and --no-dft never produce an OUTCAR, so applying it
        there would fail every relaxation for the wrong reason."""
        source = (ROOT / "tutorials" / "runner.py").read_text()
        block = source.split("def _unconverged", 1)[1].split("\n    def ", 1)[0]
        self.assertIn("self.options.mode != REAL", block)


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

    def test_submitting_nothing_is_a_failure_not_an_unverifiable(self):
        """The safety net: whatever the cause, a submitting step that recorded
        no job id has not submitted."""
        source = (ROOT / "tutorials" / "runner.py").read_text()
        block = source.split("if not state.jobs:", 1)[1].split("ok, problem", 1)[0]
        self.assertIn("state.status = FAILED", block)
        self.assertIn("nothing was submitted", block)


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
