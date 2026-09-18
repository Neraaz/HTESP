"""``--dry-run`` must be safe, and batch.header must not need a magic line.

Both were found by running the tutorial runner against the real package:
``mainprogram 20 --dry-run`` deleted exactly what the user was checking, and
``generate_submission_file`` required a ``submission here`` line that no
shipped ``batch.header`` contains, so it wrote job scripts that allocated the
job and ran nothing.
"""
from __future__ import annotations

import unittest

from tests.helpers import TempProject

from htesp.workflow import PLACEHOLDER, HTESPWorkflow


class DryRunIsNonDestructive(TempProject):
    def setUp(self):
        super().setUp()
        self.write_input_in(start=1, end=2, track="mpid.in")
        self.write_track("mpid.in", [("mp-763", "Mg1B2")])
        self.relax = self.root / "Rmp-763-Mg1B2" / "relax"
        self.calc = self.root / "Rmp-763-Mg1B2" / "calc"
        self.relax.mkdir(parents=True)
        self.calc.mkdir(parents=True)
        (self.relax / "Mg1B2.xml").write_text("x")
        (self.relax / "Mg1B2.save").mkdir()
        (self.calc / "lambda.out").write_text("lambda line\n")
        (self.calc / "Mg1B2.wfc1").write_text("wave function")
        (self.calc / "_ph0").mkdir()
        self.pressure = self.root / "Rmp-763-Mg1B2" / "pressure"
        self.pressure.mkdir()
        (self.pressure / "keep.txt").write_text("results")

    def test_clean_scan_deletes_nothing_under_dry_run(self):
        wf = HTESPWorkflow(root=self.root, workers=1, dry_run=True)
        wf.clean_scan(1, 2, "mpid.in")
        self.assertTrue((self.calc / "Mg1B2.wfc1").is_file())
        self.assertTrue((self.calc / "_ph0").is_dir())
        self.assertTrue((self.relax / "Mg1B2.xml").is_file())

    def test_clean_scan_does_delete_for_real(self):
        wf = HTESPWorkflow(root=self.root, workers=1, dry_run=False)
        wf.clean_scan(1, 2, "mpid.in")
        self.assertFalse((self.calc / "Mg1B2.wfc1").exists())
        self.assertFalse((self.calc / "_ph0").exists())

    def test_pressure_reset_deletes_nothing_under_dry_run(self):
        wf = HTESPWorkflow(root=self.root, workers=1, dry_run=True)
        wf.pressure_reset(1, 2, "mpid.in")
        self.assertTrue((self.pressure / "keep.txt").is_file())

    def test_pressure_reset_does_delete_for_real(self):
        wf = HTESPWorkflow(root=self.root, workers=1, dry_run=False)
        wf.pressure_reset(1, 2, "mpid.in")
        self.assertFalse(self.pressure.exists())

    def test_a_wildcard_never_matches_a_neighbouring_compound(self):
        """`rm scf_dir/kpoint-$A-$B*` also matched MgB2O when cleaning MgB2."""
        neighbour = self.root / "Rmp-763-Mg1B2O" / "calc"
        neighbour.mkdir(parents=True)
        (neighbour / "Mg1B2O.wfc1").write_text("keep me")
        wf = HTESPWorkflow(root=self.root, workers=1, dry_run=False)
        wf.clean_scan(1, 2, "mpid.in")
        self.assertTrue((neighbour / "Mg1B2O.wfc1").is_file())


class SubmissionScripts(TempProject):
    HEADER_WITHOUT_PLACEHOLDER = (
        "#!/bin/bash\n#SBATCH --job-name=htesp\n#SBATCH --time=04:00:00\n"
        "#SBATCH --nodes=1\n")

    def setUp(self):
        super().setUp()
        self.write_input_in()

    def test_a_header_without_the_placeholder_still_gets_the_command(self):
        self.write("batch.header", self.HEADER_WITHOUT_PLACEHOLDER)
        wf = HTESPWorkflow(root=self.root, workers=1, dry_run=True)
        written = wf.generate_submission_file("qe-elph", "mpirun", 24)
        self.assertTrue(written)
        body = (self.root / "run-scf.sh").read_text()
        self.assertIn("pw.x < scf.in > scf.out", body)
        self.assertIn("#SBATCH --job-name=htesp", body)

    def test_the_placeholder_is_still_honoured_when_present(self):
        self.write("batch.header",
                   self.HEADER_WITHOUT_PLACEHOLDER + PLACEHOLDER + "\n")
        wf = HTESPWorkflow(root=self.root, workers=1, dry_run=True)
        wf.generate_submission_file("qe-elph", "mpirun", 24)
        body = (self.root / "run-q2r.sh").read_text()
        self.assertIn("q2r.x < q2r.in > q2r.out", body)
        self.assertNotIn(PLACEHOLDER, body)

    def test_the_names_are_the_ones_the_scans_look_for(self):
        """The legacy generator wrote run.sh / q2r.sh; the scans want run-*.sh."""
        self.write("batch.header", self.HEADER_WITHOUT_PLACEHOLDER)
        wf = HTESPWorkflow(root=self.root, workers=1, dry_run=True)
        names = {path.name for path in wf.generate_submission_file("qe-elph")}
        self.assertIn("run-scf.sh", names)
        self.assertIn("run-q2r.sh", names)
        self.assertNotIn("q2r.sh", names)

    def test_a_missing_header_is_reported_not_crashed_on(self):
        wf = HTESPWorkflow(root=self.root, workers=1, dry_run=True)
        self.assertEqual(wf.generate_submission_file("qe-elph"), [])

    def test_an_unknown_target_is_reported(self):
        self.write("batch.header", self.HEADER_WITHOUT_PLACEHOLDER)
        wf = HTESPWorkflow(root=self.root, workers=1, dry_run=True)
        self.assertEqual(wf.generate_submission_file("nonsense"), [])


if __name__ == "__main__":
    unittest.main()
