"""``QEText`` -- the Quantum ESPRESSO reader/writer that replaced the sed pipelines.

Nine separate ``sed '1,4d' | sed '5d'`` copies in the bash layer assumed the
exact line layout of a converged ``vc-relax`` run.  These tests pin the three
cases that layout got wrong.
"""
from __future__ import annotations

import unittest

from tests.helpers import (RELAX_OUT, SCF_IN, TempProject, UNCONVERGED_OUT,
                           VC_RELAX_OUT)

from htesp.workflow import QEText


class FinalCoordinates(TempProject):
    def test_vc_relax_keeps_the_cell_and_every_atom(self):
        path = self.write("scf.out", VC_RELAX_OUT)
        lines = QEText.final_coordinates(path)
        text = "\n".join(lines)
        self.assertIn("CELL_PARAMETERS", text)
        self.assertIn("ATOMIC_POSITIONS", text)
        self.assertEqual(text.count("B "), 2, text)
        self.assertIn("Mg", text)
        self.assertNotIn("Begin final coordinates", text)
        self.assertNotIn("End final coordinates", text)

    def test_relax_output_keeps_every_atom(self):
        """No cell block: the fixed `sed` pipeline ate the header *and* Mg."""
        path = self.write("scf.out", RELAX_OUT)
        text = "\n".join(QEText.final_coordinates(path))
        self.assertIn("ATOMIC_POSITIONS", text)
        self.assertIn("Mg", text)
        self.assertEqual(text.count("B "), 2, text)

    def test_unconverged_output_does_not_produce_a_bogus_block(self):
        path = self.write("scf.out", UNCONVERGED_OUT)
        lines = QEText.final_coordinates(path)
        # either empty, or the last positions block -- never a truncated cell
        if lines:
            self.assertIn("ATOMIC_POSITIONS", "\n".join(lines))
            self.assertNotIn("CELL_PARAMETERS", "\n".join(lines))

    def test_missing_file_is_empty_not_an_exception(self):
        self.assertEqual(QEText.final_coordinates(self.root / "absent.out"), [])


class Scrapers(TempProject):
    def test_nelec(self):
        self.assertEqual(QEText.nelec(self.write("o", VC_RELAX_OUT)), 26)

    def test_natoms(self):
        self.assertEqual(QEText.natoms(self.write("o", VC_RELAX_OUT)), 3)

    def test_is_relaxed(self):
        self.assertTrue(QEText.is_relaxed(self.write("a", VC_RELAX_OUT)))
        self.assertFalse(QEText.is_relaxed(self.write("b", UNCONVERGED_OUT)))


class Sections(unittest.TestCase):
    def test_prefix(self):
        self.assertEqual(QEText.prefix(SCF_IN), "'Mg1B2'")

    def test_kmesh(self):
        mesh, shift = QEText.kmesh(SCF_IN)
        self.assertEqual(mesh, [12, 12, 8])
        self.assertEqual(shift, [0, 0, 0])

    def test_value(self):
        self.assertIn("vc-relax", QEText.value(SCF_IN, "calculation"))
        self.assertEqual(QEText.value(SCF_IN, "nat").strip(" ,"), "3")

    def test_card_returns_only_that_card(self):
        species = QEText.card(SCF_IN, "ATOMIC_SPECIES")
        body = "\n".join(species)
        self.assertIn("Mg.upf", body)
        self.assertIn("B.upf", body)
        self.assertNotIn("K_POINTS", body)
        self.assertNotIn("0.333333", body)

    def test_set_key_replaces_in_place(self):
        lines = QEText.set_key(SCF_IN.splitlines(), "SYSTEM", "ecutwfc", 90.0)
        joined = "\n".join(lines)
        self.assertIn("90.0", joined)
        self.assertEqual(joined.count("ecutwfc"), 1)

    def test_set_key_adds_a_missing_key(self):
        lines = QEText.set_key(SCF_IN.splitlines(), "SYSTEM", "nbnd", 24)
        self.assertIn("nbnd", "\n".join(lines))

    def test_drop_removes_matching_lines_only(self):
        lines = QEText.drop(SCF_IN.splitlines(), "pseudo_dir")
        joined = "\n".join(lines)
        self.assertNotIn("pseudo_dir", joined)
        self.assertIn("outdir", joined)


if __name__ == "__main__":
    unittest.main()
