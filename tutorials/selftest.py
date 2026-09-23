#!/usr/bin/env python
"""Self-tests for the tutorial runner.

Plain :mod:`unittest`, stdlib only: ``pytest`` cannot be installed on every
machine this has to run on, and ``unittest`` classes are collected by both::

    python -m unittest tutorials.selftest -v
    pytest tutorials/selftest.py

Nothing here needs ``pymatgen``, Quantum ESPRESSO, VASP or SLURM.  The
end-to-end test replaces ``mainprogram`` with a stub script that writes the
artefacts the catalogue asks for, so the runner's control flow -- seeding,
execution, verification, checkpointing, resuming, blocking and reporting -- is
proven without the scientific stack.

Coverage:

* ``CatalogTest``       -- catalogue integrity: codes, dependencies, step ids,
                           the QE/VASP numbering offset, directories on disk.
* ``StateTest``         -- checkpoint save/load/resume round-trip.
* ``WorkdirTest``       -- ``input.in`` patching and artefact verification.
* ``JobWaitTest``       -- ``squeue`` polling and the "no squeue" path.
* ``ReportTest``        -- what the stop report says about a synthesised failure.
* ``EndToEndTest``      -- a whole dry run against a fake ``mainprogram``.
"""
from __future__ import annotations

import json
import logging
import os
import shutil
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

from tutorials import report as report_mod
from tutorials.catalog import (CATALOG, QE_TOPICS, VASP_TOPICS, Loop, Seed, Step,
                               Tutorial, normalise_code, parse_codes, select,
                               topological_order, vasp_number_to_qe_number)
from tutorials.runner import RunOptions, TutorialRunner, preflight
from tutorials.state import (BLOCKED, DONE, FAILED, SKIPPED, RunState, StepState,
                             TutorialState, step_key)
from tutorials.workdirs import missing_artifacts, patch_input_in


# --------------------------------------------------------------------------- #
class CatalogTest(unittest.TestCase):
    """The catalogue has to be internally consistent before anything runs."""

    def test_size_and_codes(self):
        self.assertEqual(len(CATALOG), 42)
        for dft, count in (("QE", 21), ("VASP", 21)):
            codes = [c for c, t in CATALOG.items() if t.dft == dft]
            self.assertEqual(len(codes), count)
            self.assertEqual(sorted(codes, key=lambda c: int(c.split("/")[1])),
                             [f"{dft}/{n}" for n in range(1, count + 1)])

    def test_every_dependency_exists(self):
        for code, tut in CATALOG.items():
            for dep in tut.depends_on:
                self.assertIn(dep, CATALOG, f"{code} depends on unknown {dep}")
                self.assertEqual(CATALOG[dep].dft, tut.dft,
                                 f"{code} depends across DFT codes")
                self.assertNotEqual(dep, code)

    def test_no_dependency_cycles(self):
        order = topological_order(list(CATALOG))
        self.assertEqual(len(order), len(CATALOG))
        position = {code: i for i, code in enumerate(order)}
        for code, tut in CATALOG.items():
            for dep in tut.depends_on:
                self.assertLess(position[dep], position[code])

    def test_step_ids_unique_within_a_tutorial(self):
        for code, tut in CATALOG.items():
            ids = [s.id for s in tut.steps]
            self.assertEqual(len(ids), len(set(ids)), f"{code} repeats a step id")
            self.assertTrue(ids, f"{code} has no steps")

    def test_step_after_refers_to_earlier_steps(self):
        for code, tut in CATALOG.items():
            ids = [s.id for s in tut.steps]
            for position, step in enumerate(tut.steps):
                for other in step.after:
                    self.assertIn(other, ids, f"{code}.{step.id} waits on {other}")
                    self.assertLess(ids.index(other), position)

    def test_loop_steps_exist_and_are_contiguous(self):
        for code, tut in CATALOG.items():
            if tut.loop is None:
                continue
            positions = [i for i, s in enumerate(tut.steps) if s.id in tut.loop.steps]
            self.assertEqual(len(positions), len(tut.loop.steps), code)
            self.assertEqual(positions, list(range(positions[0], positions[-1] + 1)),
                             f"{code}: the loop block is not contiguous")

    def test_example_directories_exist(self):
        missing = [c for c, t in CATALOG.items() if not t.directory.is_dir()]
        self.assertEqual(missing, [])

    def test_vasp_offset(self):
        self.assertEqual(len(VASP_TOPICS), len(QE_TOPICS))      # 21 each
        self.assertNotIn("elph", VASP_TOPICS)
        for number in range(1, 11):
            self.assertEqual(vasp_number_to_qe_number(number), number)
            self.assertEqual(VASP_TOPICS[number - 1], QE_TOPICS[number - 1])
        for number in range(11, 21):
            self.assertEqual(vasp_number_to_qe_number(number), number + 1)
            self.assertEqual(VASP_TOPICS[number - 1], QE_TOPICS[number])
        self.assertIsNone(vasp_number_to_qe_number(21))
        self.assertEqual(VASP_TOPICS[20], "fermisurface")
        self.assertEqual(CATALOG["QE/11"].topic, "elph")
        self.assertEqual(CATALOG["VASP/11"].topic, CATALOG["QE/12"].topic)
        self.assertEqual(CATALOG["VASP/14"].topic, CATALOG["QE/15"].topic)

    def test_code_parsing(self):
        for text in ("qe/9", "QE-9", "QE/tutorial9", " QE/9 "):
            self.assertEqual(normalise_code(text), "QE/9")
        self.assertEqual(parse_codes("QE/9,VASP/14"), ["QE/9", "VASP/14"])
        with self.assertRaises(ValueError):
            normalise_code("nine")

    def test_selection(self):
        # a bare tree name means the whole tree; it replaced the --code flag
        self.assertEqual(len(select(only=["QE/*"])), 21)
        self.assertEqual(len(select(only=["VASP/*"])), 21)
        self.assertEqual(len(select()), 42)
        self.assertEqual(select(only=["QE/9"]), ["QE/9"])
        self.assertNotIn("QE/9", select(skip=["QE/9"]))
        self.assertEqual(select(only=["QE/12", "QE/9"]), ["QE/9", "QE/12"])
        # trees and single codes mix, and skip understands a tree too
        self.assertEqual(len(select(only=["QE/*", "VASP/14"])), 22)
        self.assertEqual(len(select(skip=["VASP/*"])), 21)

    def test_a_bare_tree_name_is_accepted(self):
        self.assertEqual(normalise_code("QE"), "QE/*")
        self.assertEqual(normalise_code("vasp"), "VASP/*")
        with self.assertRaises(ValueError):
            normalise_code("nonsense")


# --------------------------------------------------------------------------- #
class StateTest(unittest.TestCase):
    """The checkpoint has to survive a round trip and drive ``--resume``."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def test_round_trip(self):
        state = RunState(path=self.tmp / "state.json")
        tut = state.tutorial("QE/9", "relaxation")
        tut.status = FAILED
        tut.steps["relax-submit"] = StepState(
            step_id="relax-submit", key="relax-submit", index=1, status=DONE,
            command=["python", "-m", "htesp", "1"], duration=1.5)
        tut.steps["energy"] = StepState(
            step_id="energy", key="energy", index=2, status=FAILED,
            exit_code=2, reason="boom", missing=["econv.csv"])
        state.save()

        again = RunState.load(self.tmp / "state.json")
        self.assertEqual(again.mode, "dry-run")
        self.assertEqual(again.tutorials["QE/9"].status, FAILED)
        self.assertEqual(again.tutorials["QE/9"].steps["energy"].exit_code, 2)
        self.assertTrue(again.is_done("QE/9", "relax-submit"))
        self.assertFalse(again.is_done("QE/9", "energy"))
        self.assertEqual(again.tutorials["QE/9"].first_failure().step_id, "energy")

    def test_missing_and_corrupt_files_give_a_fresh_state(self):
        self.assertEqual(RunState.load(self.tmp / "nope.json").tutorials, {})
        bad = self.tmp / "bad.json"
        bad.write_text("{not json")
        self.assertEqual(RunState.load(bad).tutorials, {})

    def test_reset_and_step_keys(self):
        state = RunState(path=self.tmp / "s.json")
        state.tutorial("QE/9").steps["a"] = StepState(step_id="a", key="a", status=DONE)
        state.reset(["QE/9"])
        self.assertNotIn("QE/9", state.tutorials)
        self.assertEqual(step_key("resubmit"), "resubmit")
        self.assertEqual(step_key("resubmit", 3), "resubmit#3")


# --------------------------------------------------------------------------- #
class WorkdirTest(unittest.TestCase):
    """``input.in`` editing and the artefact check that catches silent failures."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def test_patch_creates_and_edits(self):
        from tutorials.catalog import InputPatch

        lines = patch_input_in(self.tmp, None, "VASP")
        self.assertEqual(lines[5], "DFT = VASP")
        (self.tmp / "mpid-deformed.in").write_text("v1 mp-763 Mg1B2\n")
        lines = patch_input_in(self.tmp, InputPatch(start=1, end=25,
                                                    track="mpid-deformed.in"), "QE")
        self.assertEqual((lines[0], lines[1], lines[3]), ("1", "25", "mpid-deformed.in"))
        self.assertEqual((self.tmp / "input.in").read_text().splitlines()[3],
                         "mpid-deformed.in")

    def test_missing_track_file_falls_back(self):
        from tutorials.catalog import InputPatch

        (self.tmp / "mpid.in").write_text("v1 mp-763 Mg1B2\n")
        lines = patch_input_in(self.tmp, InputPatch(track="mpid-list.in"), "QE")
        self.assertEqual(lines[3], "mpid.in")

    def test_artifact_verification(self):
        self.assertEqual(missing_artifacts(self.tmp, ["econv.csv"]), ["econv.csv"])
        (self.tmp / "econv.csv").write_text("ID\n")
        (self.tmp / "Rmp-763-Mg1B2" / "relax").mkdir(parents=True)
        (self.tmp / "Rmp-763-Mg1B2" / "relax" / "scf.in").write_text("&CONTROL\n")
        self.assertEqual(missing_artifacts(self.tmp, ["econv.csv", "R*-*/relax/scf.in"]),
                         [])
        self.assertEqual(missing_artifacts(self.tmp, ["R*-*/bands/scf.in"]),
                         ["R*-*/bands/scf.in"])


# --------------------------------------------------------------------------- #
class ReportTest(unittest.TestCase):
    """The stop report is the deliverable; check it says all of it."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        log = self.tmp / "step.log"
        log.write_text("\n".join(f"line {i}" for i in range(80)))
        self.state = RunState(path=self.tmp / "state.json",
                              argv=["htesp-tutorials", "--dry-run"])
        tut = self.state.tutorial("QE/9", "Structural relaxation (QE)")
        tut.status, tut.started, tut.workdir = FAILED, 1.0, str(self.tmp)
        tut.steps["relax-submit"] = StepState(
            step_id="relax-submit", key="relax-submit", index=1, status=DONE)
        tut.steps["energy"] = StepState(
            step_id="energy", key="energy", index=2, status=FAILED, exit_code=2,
            command=[sys.executable, "-m", "htesp", "e0"], workdir=str(self.tmp),
            duration=12.5, log=str(log), expected=["econv.csv"],
            missing=["econv.csv"], reason="the command exited 0 but produced nothing")
        blocked = self.state.tutorial("QE/17", "phonopy (QE)")
        blocked.status, blocked.blocked_by = BLOCKED, "QE/9"
        blocked.reason = "QE/9 did not finish"

    def test_stopped_at_is_the_first_failure(self):
        tut, step = self.state.stopped_at()
        self.assertEqual((tut.code, step.step_id), ("QE/9", "energy"))

    def test_markdown_has_everything_needed_to_restart(self):
        text = report_mod.build_markdown(self.state, ["QE/9", "QE/17"], self.tmp)
        for needle in ("QE/9", "energy", "exit code", "2", "econv.csv", "12.5 s",
                       "-m htesp e0", str(self.tmp), "line 79",
                       "htesp-tutorials --resume --only QE/9 --from energy",
                       "## Blocked", "QE/17", "blocked by"):
            self.assertIn(needle, text, f"the report never mentions {needle!r}")
        self.assertNotIn("line 39", text)            # only the last 40 lines

    def test_json_mirrors_the_markdown(self):
        blob = report_mod.build_json(self.state, ["QE/9", "QE/17"], self.tmp)
        self.assertEqual(blob["stopped_at"]["tutorial"], "QE/9")
        self.assertEqual(blob["stopped_at"]["step_id"], "energy")
        self.assertEqual(blob["stopped_at"]["missing"], ["econv.csv"])
        self.assertEqual(len(blob["stopped_at"]["log_tail"]), 40)
        self.assertEqual(blob["blocked"][0]["blocked_by"], "QE/9")
        self.assertIn("--from energy", blob["stopped_at"]["retry"])

    def test_summary_table_counts_every_tutorial(self):
        table = report_mod.summary_table(self.state, ["QE/9", "QE/17"])
        self.assertIn("2 tutorials:", table)
        self.assertIn("FAILED", table)
        self.assertIn("blocked", table)

    def test_write_report_makes_both_files(self):
        md, js = report_mod.write_report(self.state, self.tmp, ["QE/9", "QE/17"])
        self.assertTrue(md.is_file() and js.is_file())
        json.loads(js.read_text())
        self.assertIn("WHERE IT STOPPED",
                      report_mod.console_report(self.state, self.tmp,
                                                ["QE/9", "QE/17"]))


# --------------------------------------------------------------------------- #
#  end to end, against a fake mainprogram
# --------------------------------------------------------------------------- #
#: a stand-in for ``python -m htesp``: writes the artefact named after the
#: process, counts its invocations, and honours a "fail" / "silent" process.
STUB_MAIN = '''\
import os, sys
from pathlib import Path

process = sys.argv[1] if len(sys.argv) > 1 else ""
Path("calls.log").open("a").write(" ".join(sys.argv[1:]) + "\\n")
print("stub mainprogram: process", process)
if process == "fail":
    sys.stderr.write("stub: deliberate failure\\n")
    raise SystemExit(3)
if process == "silent":
    raise SystemExit(0)                      # exits 0, writes nothing
if process == "submit":
    stage = Path("Rmp-763-Mg1B2") / "relax"
    stage.mkdir(parents=True, exist_ok=True)
    (stage / "scf.in").write_text("&CONTROL\\n")
    raise SystemExit(0)
Path(process + ".txt").write_text("written by the stub\\n")
'''


def _make_fake_htesp(root: Path) -> Path:
    """Build an importable package whose ``python -m htesp`` is :data:`STUB_MAIN`."""
    package = root / "fakeroot" / "htesp"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("")
    (package / "__main__.py").write_text(STUB_MAIN)
    return root / "fakeroot"


class EndToEndTest(unittest.TestCase):
    """Seeding, execution, verification, checkpointing, blocking and reporting."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)
        self.fakeroot = _make_fake_htesp(self.tmp)

        log = logging.getLogger("htesp.tutorials")      # keep the test output clean
        self.addCleanup(log.setLevel, log.level)
        log.setLevel(logging.CRITICAL)

        import tutorials.runner as runner_mod
        original = runner_mod.child_env
        runner_mod.child_env = lambda root=None, **kw: {     # noqa: ARG005
            **os.environ, "PYTHONPATH": str(self.fakeroot)}
        self.addCleanup(setattr, runner_mod, "child_env", original)

        self.examples = self.tmp / "examples"
        code_dir = self.examples / "QE"
        (code_dir / "pp").mkdir(parents=True)
        (code_dir / "pp" / "Mg.upf").write_text("pseudo\n")
        (code_dir / "batch.header").write_text("#!/bin/bash\n")
        (code_dir / "config.json").write_text('{"download": {"inp": {"calc": "QE"}}}')
        (code_dir / "input.in").write_text("1\n2\n200 0\nmpid.in\nphband\nDFT = QE\n")
        for number in (1, 2, 3):
            tutorial_dir = code_dir / f"tutorial{number}"
            tutorial_dir.mkdir()
            (tutorial_dir / "config.json").write_text('{"download": {"inp": '
                                                      '{"calc": "QE"}}}')
            (tutorial_dir / "mpid.in").write_text("v1 mp-763 Mg1B2\n")
        self.catalog = self._catalog()

    def _catalog(self) -> dict[str, Tutorial]:
        base = (Seed("code", ("batch.header", "config.json", "input.in")),
                Seed("code", ("pp",), overwrite=False, link=True),
                Seed("self", ("*",)))
        first = Tutorial(
            code="QE/1", dft="QE", number=1, topic="jobscript",
            title="the one that works", directory=self.examples / "QE" / "tutorial1",
            seeds=base,
            steps=(Step("make", "write an artefact", "alpha", artifacts=("alpha.txt",)),
                   Step("submit", "submit something", "submit", submits=True,
                        artifacts=("R*-*/relax/scf.in",))))
        second = Tutorial(
            code="QE/2", dft="QE", number=2, topic="relax",
            title="the one that stops", directory=self.examples / "QE" / "tutorial2",
            depends_on=("QE/1",),
            seeds=base + (Seed("QE/1", ("alpha.txt",)),),
            steps=(Step("copy", "reuse what QE/1 produced", "beta",
                        artifacts=("beta.txt", "alpha.txt")),
                   Step("nothing", "exits 0 and writes nothing", "silent",
                        artifacts=("never.txt",)),
                   Step("after", "never reached", "gamma", artifacts=("gamma.txt",))))
        third = Tutorial(
            code="QE/3", dft="QE", number=3, topic="bands",
            title="the one that is blocked", directory=self.examples / "QE" / "tutorial3",
            depends_on=("QE/2",), seeds=base,
            steps=(Step("never", "never runs", "delta", artifacts=("delta.txt",)),))
        return {"QE/1": first, "QE/2": second, "QE/3": third}

    def _options(self, **kwargs) -> RunOptions:
        return RunOptions(workdir=self.tmp / "run",
                          examples=self.examples, step_timeout=120, **kwargs)

    def test_full_run_stops_where_it_should(self):
        options = self._options()
        state = TutorialRunner(["QE/1", "QE/2", "QE/3"], options,
                               catalog=self.catalog).run()

        self.assertEqual(state.tutorials["QE/1"].status, DONE)
        self.assertEqual(state.tutorials["QE/2"].status, FAILED)
        self.assertEqual(state.tutorials["QE/3"].status, BLOCKED)
        self.assertEqual(state.tutorials["QE/3"].blocked_by, "QE/2")

        # seeding: the dependency's output and the shared pseudopotentials arrived
        workdir = options.runs_root / "QE-2"
        self.assertTrue((workdir / "alpha.txt").is_file())
        self.assertTrue((workdir / "config.json").is_file())
        self.assertTrue((options.runs_root / "QE-1" / "pp" / "Mg.upf").is_file())

        # the failure is "exited 0 and produced nothing", and it is recorded
        failed = state.tutorials["QE/2"].first_failure()
        self.assertEqual(failed.step_id, "nothing")
        self.assertEqual(failed.exit_code, 0)
        self.assertEqual(failed.missing, ["never.txt"])
        self.assertIn("produced nothing", failed.reason)
        self.assertTrue(Path(failed.log).is_file())

        # the step after the failure was never attempted
        self.assertNotIn("after", state.tutorials["QE/2"].steps)

        # the checkpoint on disk says the same thing
        reloaded = RunState.load(options.workdir / "state.json")
        self.assertEqual(reloaded.tutorials["QE/2"].first_failure().step_id, "nothing")

        # and the report points at it
        tut, step = state.stopped_at()
        self.assertEqual((tut.code, step.step_id), ("QE/2", "nothing"))
        md, js = report_mod.write_report(state, options.workdir,
                                         ["QE/1", "QE/2", "QE/3"])
        text = md.read_text()
        self.assertIn("QE/2", text)
        self.assertIn("never.txt", text)
        self.assertIn("htesp-tutorials --resume --only QE/2 --from nothing", text)
        self.assertEqual(json.loads(js.read_text())["stopped_at"]["step_id"], "nothing")

    def test_resume_does_not_rerun_finished_steps(self):
        options = self._options()
        TutorialRunner(["QE/1"], options, catalog=self.catalog).run()
        calls = (options.runs_root / "QE-1" / "calls.log").read_text().splitlines()
        self.assertEqual(len(calls), 2)

        state = RunState.load(options.workdir / "state.json")
        TutorialRunner(["QE/1"], options, catalog=self.catalog, state=state).run()
        calls_again = (options.runs_root / "QE-1" / "calls.log").read_text().splitlines()
        self.assertEqual(calls_again, calls, "--resume re-ran a finished step")

    def test_restart_reruns_everything(self):
        options = self._options()
        TutorialRunner(["QE/1"], options, catalog=self.catalog).run()
        fresh = RunState(path=options.workdir / "state.json")
        TutorialRunner(["QE/1"], self._options(resume=False), catalog=self.catalog,
                       state=fresh).run()
        calls = (options.runs_root / "QE-1" / "calls.log").read_text().splitlines()
        self.assertEqual(len(calls), 4)

    def test_from_step_skips_the_earlier_ones(self):
        options = self._options(from_step="submit")
        state = TutorialRunner(["QE/1"], options, catalog=self.catalog).run()
        steps = state.tutorials["QE/1"].steps
        self.assertEqual(steps["make"].status, SKIPPED)
        self.assertIn("--from", steps["make"].reason)
        self.assertEqual(steps["submit"].status, DONE)

    def test_dry_run_passes_the_flag_through(self):
        options = self._options()
        TutorialRunner(["QE/1"], options, catalog=self.catalog).run()
        calls = (options.runs_root / "QE-1" / "calls.log").read_text()
        self.assertIn("--dry-run", calls)

    def test_preflight_reports_problems_without_running(self):
        options = self._options()
        options.workdir.mkdir(parents=True, exist_ok=True)
        problems = preflight(["QE/1", "QE/404"], options, catalog=self.catalog)
        messages = " ".join(p.message for p in problems)
        self.assertIn("QE/404", messages)
        self.assertTrue(any(p.fatal for p in problems))


if __name__ == "__main__":                              # pragma: no cover
    unittest.main(verbosity=2)
