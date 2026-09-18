"""The workflow layer's shared machinery: tracking files, meshes, scratch, jobs."""
from __future__ import annotations

import json
import os
import unittest
from pathlib import Path

from tests.helpers import TempProject

from htesp.workflow import (PROCESSES, SCRIPTS, HTESPWorkflow, Material,
                            Scheduler, pushd, read_mesh_file, read_track_file)


class TrackFile(TempProject):
    def setUp(self):
        super().setUp()
        self.write_track("mpid.in", [("mp-1", "A"), ("mp-2", "B"), ("mp-3", "C")])

    def test_end_is_exclusive(self):
        found = read_track_file(self.root / "mpid.in", 1, 3, self.root)
        self.assertEqual([m.mpid for m in found], ["mp-1", "mp-2"])

    def test_a_missing_index_is_skipped_not_turned_into_R_dash(self):
        """Bash built `R-/relax` and then ran `rm -r` from the project root."""
        self.write("sparse.in", "v1 mp-1 A\nv3 mp-3 C\n")
        found = read_track_file(self.root / "sparse.in", 1, 4, self.root)
        self.assertEqual([m.index for m in found], [1, 3])
        self.assertTrue(all(m.mpid for m in found))

    def test_blank_and_comment_lines_are_ignored(self):
        self.write("messy.in", "\nv1 mp-1 A\n   \nnot-a-row\nv2 mp-2 B\n")
        found = read_track_file(self.root / "messy.in", 1, 3, self.root)
        self.assertEqual(len(found), 2)


class MaterialPaths(TempProject):
    def test_names_and_directories(self):
        m = Material(1, "mp-763", "Mg1B2", self.root)
        self.assertEqual(m.name, "Rmp-763-Mg1B2")
        self.assertEqual(m.dir, self.root / "Rmp-763-Mg1B2")
        self.assertEqual(m.sub("relax"), self.root / "Rmp-763-Mg1B2" / "relax")
        self.assertEqual(m.scf_template.name, "scf-mp-763.in")
        self.assertEqual(m.scf_relaxed.name, "scf-relax-mp-763-Mg1B2.in")

    def test_phonon_folder_falls_back_to_calc(self):
        m = Material(1, "mp-763", "Mg1B2", self.root)
        self.assertEqual(m.phonon_folder(), "calc")
        m.sub("phonon").mkdir(parents=True)
        self.assertEqual(m.phonon_folder(), "phonon")


class MeshFile(TempProject):
    """One parser replaced five copies of the qpoint.in/kpoint.in logic."""

    def test_missing_file_halves_the_base_mesh(self):
        mesh, shift = read_mesh_file(self.root / "absent.in", [12, 12, 8])
        self.assertEqual(mesh, [6, 6, 4])
        self.assertEqual(shift, [0, 0, 0])

    def test_single_divisor(self):
        self.write("q.in", "4\n")
        self.assertEqual(read_mesh_file(self.root / "q.in", [12, 12, 8])[0], [3, 3, 2])

    def test_explicit_mesh(self):
        self.write("q.in", "3 3 2\n")
        self.assertEqual(read_mesh_file(self.root / "q.in", [12, 12, 8])[0], [3, 3, 2])

    def test_divisor_with_shift(self):
        mesh, shift = read_mesh_file(self.write("q.in", "2 1 1 1\n"), [12, 12, 8])
        self.assertEqual(mesh, [6, 6, 4])
        self.assertEqual(shift, [1, 1, 1])

    def test_mesh_with_shift(self):
        mesh, shift = read_mesh_file(self.write("q.in", "4 4 4 1 0 1\n"), [12, 12, 8])
        self.assertEqual(mesh, [4, 4, 4])
        self.assertEqual(shift, [1, 0, 1])

    def test_never_returns_zero(self):
        """A zero in K_POINTS automatic is an invalid QE input."""
        mesh, _ = read_mesh_file(self.root / "absent.in", [1, 1, 1])
        self.assertTrue(all(value >= 1 for value in mesh), mesh)


class Pushd(TempProject):
    def test_cwd_is_restored_even_on_an_exception(self):
        (self.root / "sub").mkdir()
        before = Path.cwd()
        with self.assertRaises(ValueError):
            with pushd(self.root / "sub"):
                raise ValueError("boom")
        self.assertEqual(Path.cwd(), before)

    def test_missing_directory_raises_instead_of_running_in_the_project_root(self):
        """`cd R$A-$B/... || :` is how bash came to `rm -r` in the root."""
        with self.assertRaises(OSError):
            with pushd(self.root / "does-not-exist"):
                pass


class JobRecords(TempProject):
    def test_ids_round_trip(self):
        stage = self.root / "Rmp-1-A" / "relax"
        stage.mkdir(parents=True)
        Scheduler.record(stage, "12345", "scf")
        Scheduler.record(stage, "12346", "scf")
        self.assertEqual(Scheduler.job_ids(stage), ["12345", "12346"])
        stored = json.loads((stage / ".htesp_job.json").read_text())
        self.assertIn("scf", stored)

    def test_an_empty_id_is_not_recorded(self):
        stage = self.root / "s"
        stage.mkdir()
        Scheduler.record(stage, "", "scf")
        self.assertEqual(Scheduler.job_ids(stage), [])

    def test_missing_store_is_empty(self):
        self.assertEqual(Scheduler.job_ids(self.root), [])


class CommandTables(unittest.TestCase):
    def test_every_legacy_bash_script_has_a_workflow_method(self):
        missing = [name for name, method in SCRIPTS.items()
                   if not hasattr(HTESPWorkflow, method)]
        self.assertEqual(missing, [])

    def test_every_numbered_process_has_a_workflow_method(self):
        missing = [key for key, method in PROCESSES.items()
                   if not hasattr(HTESPWorkflow, method)]
        self.assertEqual(missing, [])

    def test_the_legacy_scripts_are_all_covered(self):
        legacy = Path(__file__).resolve().parent.parent / "legacy" / "bash"
        if not legacy.is_dir():
            self.skipTest("legacy/bash not present")
        names = {p.name for p in legacy.iterdir() if p.is_file()}
        # jobscript.sh is sourced, not dispatched; it has a shim of its own
        uncovered = names - set(SCRIPTS) - {"jobscript.sh"}
        self.assertEqual(uncovered, set())


class WorkflowConstruction(TempProject):
    def test_it_reads_input_in_and_the_config(self):
        self.write_input_in(dft="VASP")
        wf = HTESPWorkflow(root=self.root, workers=1, dry_run=True)
        self.assertTrue(wf.is_vasp)
        self.assertEqual(wf.input.nkpt, 200)
        self.assertIn("pseudo", wf.config)

    def test_workers_is_at_least_one(self):
        wf = HTESPWorkflow(root=self.root, workers=0, dry_run=True)
        self.assertGreaterEqual(wf.workers, 1)

    def test_environment_override(self):
        os.environ["HTESP_WORKERS"] = "3"
        try:
            self.assertEqual(HTESPWorkflow(root=self.root, dry_run=True).workers, 3)
        finally:
            os.environ.pop("HTESP_WORKERS")

    def test_a_missing_track_file_is_a_clear_error(self):
        wf = HTESPWorkflow(root=self.root, workers=1, dry_run=True)
        with self.assertRaises(FileNotFoundError):
            wf.materials(1, 3, "nowhere.in")

    def test_failure_accounting_starts_empty(self):
        wf = HTESPWorkflow(root=self.root, workers=1, dry_run=True)
        self.assertEqual(wf.failed_count, 0)

    def test_a_failing_body_is_counted_not_swallowed(self):
        """Bash discarded every exit code: a stage where all materials blew
        up still reported success and the next stage ran on nothing."""
        self.write_track("mpid.in", [("mp-1", "A"), ("mp-2", "B")])
        wf = HTESPWorkflow(root=self.root, workers=1, dry_run=True)
        materials = wf.materials(1, 3, "mpid.in")

        def explode(material):
            raise RuntimeError("no scf.out")

        results = wf.map(explode, materials, parallel=False)
        self.assertEqual(len(results), 2)
        self.assertTrue(all(r.status == "failed" for r in results))
        self.assertEqual(wf.failed_count, 2)
        self.assertIn("no scf.out", wf.failure_summary())

    def test_results_come_back_in_material_order(self):
        self.write_track("mpid.in", [(f"mp-{i}", f"C{i}") for i in range(1, 6)])
        wf = HTESPWorkflow(root=self.root, workers=1, dry_run=True)
        materials = wf.materials(1, 6, "mpid.in")
        from htesp.workflow import Result

        def body(material):
            return Result(material.index, material.mpid, material.compound)

        results = wf.map(body, materials, parallel=False)
        self.assertEqual([r.index for r in results], [1, 2, 3, 4, 5])


if __name__ == "__main__":
    unittest.main()
