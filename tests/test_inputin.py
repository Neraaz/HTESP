"""``input.in`` parsing -- the four historic failure modes plus round-tripping."""
from __future__ import annotations

import unittest

from tests.helpers import TempProject

from htesp.inputin import InputIn, InputInError


class Parsing(TempProject):
    def test_full_six_line_file(self):
        inp = InputIn.parse("2\n30\n200 1\nmpid-list.in\nphband dos\nDFT = VASP\n")
        self.assertEqual((inp.start, inp.end, inp.nkpt, inp.kcut), (2, 30, 200, 1))
        self.assertEqual(inp.track, "mpid-list.in")
        self.assertEqual(inp.plot_types, ["phband", "dos"])
        self.assertTrue(inp.is_vasp)
        self.assertFalse(inp.is_qe)

    def test_four_line_file_does_not_raise_index_error(self):
        """`if len(lines) >= 4:` and then `lines[4]` -- the original crashed."""
        inp = InputIn.parse("1\n5\n120\nmpid.in\n")
        self.assertEqual(inp.plot_types, ["phband"])
        self.assertEqual(inp.nkpt, 120)

    def test_three_line_file_reports_what_is_wrong(self):
        """The original left start/end/element unbound and raised NameError."""
        with self.assertRaises(InputInError) as caught:
            InputIn.parse("1\n2\n3\n")
        self.assertIn("at least four", str(caught.exception))

    def test_missing_dft_line_defaults_to_qe_rather_than_an_empty_string(self):
        """31 bash sites did `[ $dft == 'vasp' ]`, which aborts on an empty $dft."""
        inp = InputIn.parse("1\n5\n120\nmpid.in\nphband\n")
        self.assertEqual(inp.dft, "QE")
        self.assertFalse(inp.is_vasp)

    def test_end_must_be_greater_than_start(self):
        with self.assertRaises(InputInError):
            InputIn.parse("5\n5\n120\nmpid.in\n")

    def test_non_numeric_start_is_reported_by_name(self):
        with self.assertRaises(InputInError) as caught:
            InputIn.parse("first\n5\n120\nmpid.in\n")
        self.assertIn("line 1", str(caught.exception))

    def test_end_is_exclusive(self):
        self.assertEqual(InputIn.parse("1\n4\n10\nt.in\n").count, 3)


class Creation(TempProject):
    def test_created_file_has_a_list_of_plot_types_not_a_string(self):
        """The original wrote plot_type='phband' and then iterated its letters."""
        inp = InputIn.load_or_create(self.root / "input.in", dft="QE")
        self.assertIsInstance(inp.plot_types, list)
        self.assertEqual(inp.plot_types, ["phband"])
        reread = InputIn.read(self.root / "input.in")
        self.assertEqual(reread.plot_types, ["phband"])

    def test_created_file_always_carries_the_dft_line(self):
        InputIn.load_or_create(self.root / "input.in", dft="VASP")
        self.assertIn("DFT = VASP", (self.root / "input.in").read_text())

    def test_round_trip(self):
        original = InputIn(start=3, end=9, nkpt=150, kcut=2, track="mpid-2.in",
                           plot_types=["band", "dos", "a2f"], dft="VASP")
        original.write(self.root / "input.in")
        again = InputIn.read(self.root / "input.in")
        for field in ("start", "end", "nkpt", "kcut", "track", "plot_types"):
            self.assertEqual(getattr(again, field), getattr(original, field), field)
        self.assertTrue(again.is_vasp)


class TolerantLoad(TempProject):
    """`InputIn.load` is what the workflow layer uses: it must never raise."""

    def test_missing_file(self):
        inp = InputIn.load(self.root / "nope.in")
        self.assertEqual(inp.dft, "qe")
        self.assertEqual(inp.track, "mpid.in")

    def test_garbage_lines_fall_back_to_the_defaults(self):
        self.write("input.in", "not a number\nalso not\n\n\n")
        inp = InputIn.load(self.root / "input.in")
        self.assertEqual(inp.start, 1)
        self.assertEqual(inp.end, 2)

    def test_dft_comes_from_the_config_when_the_line_is_absent(self):
        self.write("input.in", "1\n3\n100\nmpid.in\n")
        inp = InputIn.load(self.root / "input.in",
                           {"download": {"inp": {"calc": "VASP"}}})
        self.assertTrue(inp.is_vasp)

    def test_load_lowercases_dft_for_the_workflow_layer(self):
        self.write("input.in", "1\n3\n100\nmpid.in\nphband\nDFT = VASP\n")
        self.assertEqual(InputIn.load(self.root / "input.in").dft, "vasp")


class TrackFileChecks(TempProject):
    def test_missing_track_file_is_reported(self):
        inp = InputIn.parse("1\n3\n100\nmpid.in\n")
        problems = inp.check_track_file(self.root)
        self.assertTrue(any("not found" in p for p in problems))

    def test_range_past_the_end_is_reported(self):
        self.write_track("mpid.in", [("mp-1", "A"), ("mp-2", "B")])
        inp = InputIn.parse("1\n9\n100\nmpid.in\n")
        problems = inp.check_track_file(self.root)
        self.assertTrue(any("past the" in p for p in problems), problems)

    def test_the_documented_off_by_one_is_not_flagged(self):
        """`end` is exclusive, so end == len(entries) + 1 is correct."""
        self.write_track("mpid.in", [("mp-1", "A"), ("mp-2", "B")])
        inp = InputIn.parse("1\n3\n100\nmpid.in\n")
        self.assertEqual(inp.check_track_file(self.root), [])


if __name__ == "__main__":
    unittest.main()
