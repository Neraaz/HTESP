"""One test per defect the 2026 review found.

Most of these modules need pymatgen, ase or scipy, which may not be installed.
Where the fix can be checked by executing the function it is, and where it
cannot the test asserts on the source, so the fix cannot be silently reverted.
Every test names the original symptom.
"""
from __future__ import annotations

import ast
import tempfile
import unittest
from pathlib import Path

from tests.helpers import have, skip_without

PKG = Path(__file__).resolve().parent.parent / "htesp"
ROOT = PKG.parent


def source(module: str) -> str:
    """The module's full text, comments and docstrings included."""
    return (PKG / f"{module}.py").read_text()


def code_only(module: str) -> str:
    """The module's text with comments and string literals removed.

    Every fix carries a ``# FIX(n): <what was wrong>`` comment naming the old
    pattern, so a naive "the old pattern is gone" assertion would match the
    comment describing it.  Negative assertions run against this instead.
    """
    import io
    import tokenize

    pieces = []
    text = source(module)
    try:
        tokens = tokenize.generate_tokens(io.StringIO(text).readline)
        for token in tokens:
            if token.type in (tokenize.COMMENT, tokenize.STRING):
                continue
            pieces.append(token.string)
    except (tokenize.TokenError, IndentationError):      # pragma: no cover
        return text
    return "\n".join(pieces)


class ElasticConstants(unittest.TestCase):
    """The moduli were 10x too large (VASP) and had the wrong sign (QE)."""

    def test_vasp_stress_is_converted_kbar_to_gpa(self):
        text = source("elastic")
        self.assertIn("KBAR_TO_GPA = 0.1", text)
        self.assertNotIn("-1.0*np.array", code_only("elastic"))

    def test_the_last_ionic_step_is_used_not_the_first(self):
        text = source("elastic")
        self.assertIn("ionic_steps[-1]", text)
        self.assertNotIn("ionic_steps[0]", code_only("elastic"))

    def test_the_qe_path_reads_the_kbar_columns(self):
        """It took the Ry/bohr^3 columns and applied the Ry/angstrom^3 factor."""
        text = source("elastic")
        self.assertNotIn("21798.7", code_only("elastic"))

    def test_deformed_cells_are_not_restandardised(self):
        """primitive=True snapped the 0.5-1% strains away and rotated the frame."""
        self.assertIn("primitive=False", source("elastic"))

    def test_csv_header_and_rows_have_the_same_width(self):
        text = source("elastic")
        self.assertNotIn("'structure',", text.split("PROPNAME")[1][:400]
                         if "PROPNAME" in text else "")


class QuantumEspressoInputs(unittest.TestCase):
    def test_matdyn_asks_for_crystal_q_coordinates(self):
        """create_matdyn omitted q_in_cryst_coord; wrong for every non-cubic cell."""
        self.assertIn("q_in_cryst_coord", source("htepc"))

    def test_generate_kpath_writes_crystal_k_points(self):
        """Fractional ASE k-points read as tpiba give a wrong band path."""
        self.assertIn("crystal", source("htepc"))
        self.assertIn("qualifier", source("htepc"))

    def test_one_symprec_for_the_mesh_and_the_cell(self):
        """getkpt used 0.1, setting_qeinput 0.01: two different cells."""
        self.assertIn("SYMPREC = 0.01", source("htepc"))
        self.assertNotIn("symprec=0.1", code_only("htepc").replace(" ", ""))

    def test_the_sssp_table_covers_oxygen(self):
        """The shipped table lacked O, Se and Au, so every oxide raised KeyError."""
        import json
        table = json.loads((PKG / "data" / "config.json").read_text())
        for element in ("O", "Se", "Au"):
            self.assertIn(element, table["pseudo"]["PSEUDO"])


class ImportTimeFailures(unittest.TestCase):
    def test_the_wanniertools_key_is_spelled_with_an_s(self):
        """`config['wanniertool_input']` was a KeyError at import: wt1/wt2 dead."""
        text = source("create_wt_inputs")
        self.assertIn("wanniertools_input", text)
        self.assertNotIn("wanniertool_input", code_only("create_wt_inputs"))

    def test_the_epw_writer_imports_ase_read(self):
        """`read(...)` in ciftoxsf was never imported -> NameError."""
        self.assertIn("from ase.io import read", source("create_epw_inputs"))

    def test_the_epw_writer_does_not_write_to_a_closed_file(self):
        """epw_write.write ran after its `with open(...)` block had closed."""
        tree = ast.parse(source("create_epw_inputs"))
        for node in ast.walk(tree):
            if not isinstance(node, ast.FunctionDef):
                continue
            with_blocks = [n for n in node.body if isinstance(n, ast.With)]
            if not with_blocks:
                continue
            after = node.body[node.body.index(with_blocks[-1]) + 1:]
            for statement in after:
                for inner in ast.walk(statement):
                    if (isinstance(inner, ast.Attribute) and inner.attr == "write"
                            and isinstance(inner.value, ast.Name)
                            and "write" in inner.value.id):
                        self.fail(f"{node.name}: writes after the with-block closed")

    def test_standard_module_was_retired(self):
        """Dead and broken: Cell.get_bravais_lattice(Atoms), BravaisLattice[i][j]."""
        self.assertFalse((PKG / "standard.py").exists())
        self.assertTrue((ROOT / "_removed" / "standard.py").exists())


class Plotting(unittest.TestCase):
    def test_lsorbit_is_parsed_not_counted(self):
        """`grep LSORBIT INCAR | wc -l` counted `LSORBIT = .FALSE.` as on."""
        text = source("plot")
        self.assertNotIn("grep LSORBIT", code_only("plot"))
        self.assertIn("read_lsorbit", text)

    def test_the_fermi_level_is_the_last_match(self):
        """QE: the first match is the first ionic step of a relaxation."""
        text = source("plot")
        self.assertIn("fermi_from_qe_out", text)
        self.assertIn("fermi_from_outcar", text)

    def test_both_vasp_wordings_are_accepted(self):
        """One branch matched VASP 5 (`E-fermi`), the other VASP 6."""
        text = source("plot")
        self.assertIn("E-fermi", text)
        self.assertIn("Fermi energy", text)

    def test_legend_columns_are_at_least_one(self):
        compact = code_only("plot").replace(" ", "").replace("\n", "")
        self.assertIn("max(1,int(ndos_data/4))", compact)
        self.assertNotIn("ncol=int(", compact)

    def test_band_projection_variable_is_assigned(self):
        """procar_jband was used but its assignment was commented out."""
        text = source("plot_bandproj")
        self.assertIn("procar_jband =", text)


class Scheduling(unittest.TestCase):
    def test_which_calc_is_normalised_in_one_place(self):
        """`which_calc: "WANNIER"` left submission_files unbound -> NameError."""
        self.assertIn("def normalise_calc", source("generate_submission"))

    def test_submission_uses_the_shared_config_loader(self):
        text = source("generate_submission")
        self.assertIn("from htesp.check_json import config", text)
        self.assertNotIn("./config.json", code_only("generate_submission"))


class Physics(unittest.TestCase):
    @skip_without("numpy")
    def test_allen_dynes_returns_nan_below_the_physical_threshold(self):
        """np.abs() in the denominator manufactured a finite Tc for lambda<0.19."""
        if not have("mp_api") or not have("pymatgen"):
            # allen_dynes is dependency-free; read and exec just that function
            import math
            import numpy as np
            text = source("ml_processing")
            tree = ast.parse(text)
            func = next(n for n in tree.body
                        if isinstance(n, ast.FunctionDef) and n.name == "allen_dynes")
            namespace = {"np": np, "math": math, "MUSTAR": 0.1,
                         "warnings": __import__("warnings")}
            exec(compile(ast.Module([func], []), "<allen>", "exec"), namespace)
            allen_dynes = namespace["allen_dynes"]
        else:
            from htesp.ml_processing import allen_dynes
        import math
        self.assertTrue(math.isnan(allen_dynes(500.0, 0.10)))
        self.assertGreater(allen_dynes(500.0, 1.0), 0.0)

    def test_precision_and_recall_labels_are_not_swapped(self):
        text = source("ml_processing")
        index_recall = text.find("recall_score:")
        index_precision = text.find("precision_score:")
        if index_recall == -1 or index_precision == -1:
            self.skipTest("labels rewritten")
        window = text[min(index_recall, index_precision):
                      max(index_recall, index_precision) + 200]
        self.assertNotIn("recall_score:\", precision", window)

    def test_phonon_displacements_use_arrays(self):
        """A nested list from YAML divided by a float raised TypeError."""
        self.assertIn("np.array(self.eigen[i]", source("displace_phonopy"))


class VaspInputHandling(unittest.TestCase):
    def test_incar_keys_are_matched_whole(self):
        """`sed -i '/ENCUT/d'` also deleted ENCUTGW; NELM ate NELMIN/NELMDL."""
        text = source("vasp_process")
        self.assertIn("drop_incar_keys", text)
        self.assertNotIn("sed -i", code_only("vasp_process"))

    def test_vasp_in_is_parsed_as_pairs(self):
        """Keys and values in two independent lists shifted on a delete line."""
        self.assertIn("def parse_vasp_in", source("vasp_process"))

    def test_phonopy_band_path_keeps_discontinuities(self):
        """Stripping ',' created spurious segments in the phonopy BAND= line."""
        self.assertIn("split_path_labels", source("vasp_process"))

    @skip_without("ase", "pymatgen")
    def test_drop_incar_keys_behaviour(self):
        """`sed -i '/ENCUT/d'` also deleted ENCUTGW, and `/NELM/d` NELMIN.

        The signature is ``drop_incar_keys(keys, path="INCAR")``: it rewrites
        the file in place and returns how many lines it removed.  The previous
        version of this test passed a list of lines as ``keys`` and a list of
        keys as ``path``, so it exercised nothing.
        """
        from htesp.vasp_process import drop_incar_keys
        with tempfile.TemporaryDirectory() as tmp:
            incar = Path(tmp) / "INCAR"
            incar.write_text("ENCUT = 500\nENCUTGW = 300\nNELM = 60\nNELMIN = 4\n")
            removed = drop_incar_keys(["ENCUT", "NELM"], str(incar))
            kept = incar.read_text()
        self.assertEqual(removed, 2)
        self.assertIn("ENCUTGW = 300", kept)
        self.assertIn("NELMIN = 4", kept)
        self.assertNotIn("ENCUT = 500", kept)
        self.assertNotIn("NELM = 60", kept)


class PathLabels(unittest.TestCase):
    def test_kcutoff_slices_labels_not_characters(self):
        """`G1`/`K1` were cut in half by a character slice of the path string."""
        text = source("kpath")
        self.assertIn("path_tokens", text)

    def test_discontinuities_do_not_create_a_bogus_segment(self):
        """`K|U` was emitted as a consecutive pair -> a K->U segment."""
        self.assertIn("split_path_segments", source("kpoint_path"))


class SafetyNets(unittest.TestCase):
    def test_the_workflow_never_deletes_with_a_bare_glob(self):
        """`rm scf_dir/kpoint-$A-$B*` also matched MgB2O when cleaning MgB2."""
        text = source("workflow")
        self.assertIn("def remove_glob", text)

    def test_the_workflow_counts_failures(self):
        self.assertIn("failed_count", source("workflow"))

    def test_job_ids_are_captured(self):
        """`squeue | grep "$B"` matched every job containing B, C, Si..."""
        text = source("workflow")
        self.assertIn("--parsable", text)
        self.assertIn(".htesp_job.json", text)


class MagneticOrderingFilter(unittest.TestCase):
    """`mainprogram search` found 465 compounds and wrote 2 to mpid-list.in.

    ``emmet.core`` declares ``Ordering`` as a plain ``Enum``, so
    ``str(Ordering.NM)`` is ``'Ordering.NM'`` and ``Ordering.NM == 'NM'`` is
    False.  ``element_extract.download()`` wrote ``str(prop)`` into
    ``download/data-<elm>.csv``, so the ordering column held ``Ordering.NM``
    and the ``extract()`` filter ``data['ordering'] == 'NM'`` -- the value
    config.json documents -- matched nothing.  The chemsys path compared the
    raw member, so ``mag_logic`` was False for every entry.
    """

    def test_plain_value_unwraps_an_enum(self):
        from enum import Enum

        from htesp.element_extract import plain_value

        class Ordering(Enum):        # exactly how emmet declares it
            NM = "NM"
            FM = "FM"

        self.assertNotEqual(str(Ordering.NM), "NM")      # the trap
        self.assertEqual(plain_value(Ordering.NM), "NM")
        self.assertEqual(str(plain_value(Ordering.NM)), "NM")
        # anything that is not an Enum is returned untouched
        for value in ("NM", 3, 0.5, None, ["a"]):
            self.assertEqual(plain_value(value), value)

    def test_emmet_ordering_really_is_a_bare_enum(self):
        """If emmet ever makes Ordering a StrEnum this test says so."""
        if not have("emmet"):
            self.skipTest("needs emmet-core")
        try:
            from emmet.core.electronic_structure import Ordering
        except ImportError:
            self.skipTest("emmet moved Ordering")
        self.assertNotEqual(
            str(Ordering.NM), "NM",
            "emmet's Ordering now stringifies to its value; plain_value() is "
            "still correct, but the comment explaining why can be simplified")

    def test_both_extraction_paths_normalise_the_field(self):
        """The CSV column and the chemsys comparison must both go through it."""
        # code_only() splits into tokens, so match against the raw source;
        # none of these call expressions appear in a comment or docstring.
        text = source("element_extract")
        self.assertIn("plain_value(search.dict()[propty])", text)
        self.assertIn("plain_value(obj.data['ordering'])", text)
        self.assertIn("plain_value(obj.data[propty])", text)
        self.assertIn("plain_value(obj.data['ordering'])", source("crystal"))


    def test_material_ids_keep_the_legacy_mp_integer_spelling(self):
        """mpid-list.in switched from `mp-149` to `mp-aaaaft`.

        Materials Project is migrating ids from MPID to AlphaID. 1.x wrote
        ``search['material_id'].string``; the rewrite wrote
        ``search.dict()['material_id']`` and stringified it, which for an
        AlphaID is the new alphabetic spelling. Every R<mpid>-<compound>
        directory is named from this, so a running campaign stopped matching
        its own directories.
        """
        from htesp.element_extract import legacy_mpid

        self.assertEqual(legacy_mpid("mp-763"), "mp-763")

        class FakeAlphaID(str):                 # str() is the alpha spelling
            string = "mp-149"                   # .string is the legacy one

        self.assertEqual(legacy_mpid(FakeAlphaID("mp-aaaaft")), "mp-149")

        text = source("element_extract")
        # the *attribute*, not .dict()[...]: see the docstring of legacy_mpid
        self.assertIn("legacy_mpid(getattr(search, propty", text)
        self.assertIn("legacy_mpid(entries[i].data['material_id'])", text)

    def test_legacy_mpid_decodes_a_bare_alphabetic_string(self):
        """`.dict()` hands back a plain 'mp-bdj' str with no `.string` left.

        The first version of this fix only consulted `.string`, so it returned
        'mp-bdj' unchanged on the one path that actually matters.
        """
        from htesp.element_extract import legacy_mpid

        self.assertEqual(legacy_mpid("mp-bdj"), "mp-763")     # MgB2
        self.assertEqual(legacy_mpid("mp-763"), "mp-763")     # already legacy
        self.assertEqual(legacy_mpid("mvc-bvm"), "mvc-1234")  # prefix preserved
        # identifiers that are not Materials Project ids pass through untouched
        self.assertEqual(legacy_mpid("12345"), "12345")
        self.assertEqual(legacy_mpid("aflow:abc123"), "aflow:abc123")

    @unittest.skipUnless(have("emmet"), "needs emmet-core")
    def test_against_the_real_alphaid(self):
        try:
            from emmet.core.mpid import MPID, AlphaID
        except ImportError:
            self.skipTest("this emmet-core has no AlphaID")
        from htesp.element_extract import legacy_mpid

        alpha = AlphaID(149, padlen=6, prefix="mp")
        self.assertEqual(str(alpha), "mp-aaaaft")       # the trap
        self.assertEqual(legacy_mpid(alpha), "mp-149")
        self.assertEqual(legacy_mpid(MPID("mp-763")), "mp-763")


    @skip_without("pandas")
    def test_ordering_filter_is_skippable_and_accepts_a_list(self):
        """465 hits -> 2 rows, because MP now says `Unknown` for 202 of them.

        `ordering` was the only filter in extract() with no on/off switch,
        so a configuration could not say "any ordering" or "NM or Unknown".
        """
        import pandas as pd

        from htesp.element_extract import filter_ordering

        data = pd.DataFrame({
            "material_id": ["a", "b", "c", "d"],
            "ordering": ["NM", "Unknown", "FM", "NM"],
        })
        # a bare string keeps the historical behaviour
        self.assertEqual(list(filter_ordering(data, "NM")["material_id"]),
                         ["a", "d"])
        # null switches the filter off
        self.assertEqual(len(filter_ordering(data, None)), 4)
        # a list keeps any of its values -- the case MP's data now needs
        self.assertEqual(
            list(filter_ordering(data, ["NM", "Unknown"])["material_id"]),
            ["a", "b", "d"])
        # an empty list means "no filter", never "match nothing"
        self.assertEqual(len(filter_ordering(data, [])), 4)
        # the index is reset, as the other filters in extract() do
        self.assertEqual(list(filter_ordering(data, ["Unknown", "FM"]).index),
                         [0, 1])

    @skip_without("pandas")
    def test_ordering_filter_tolerates_an_enum(self):
        """Belt and braces: a member must not become 'Ordering.NM' here."""
        from enum import Enum

        import pandas as pd

        from htesp.element_extract import filter_ordering

        class Ordering(Enum):
            NM = "NM"

        data = pd.DataFrame({"ordering": ["NM", "FM"]})
        self.assertEqual(len(filter_ordering(data, Ordering.NM)), 1)
        self.assertEqual(len(filter_ordering(data, [Ordering.NM])), 1)

    def test_extract_routes_through_the_helper(self):
        self.assertIn("filter_ordering(data, ordering)",
                      source("element_extract"))


    def test_the_cutoff_tables_cover_the_sssp_elements(self):
        """He, Ne, Ar, Kr, Xe and every lanthanide were absent from both
        tables, so any compound containing one raised KeyError in
        getecut_sssp() -- and after FIX(2) merely warned and fell back.
        Filled from SSSP 1.3.0 PBE efficiency.
        """
        import json

        from htesp.htepc import SSSP_EFFICIENCY

        shipped = json.loads(
            (PKG / "data" / "config.json").read_text())["pseudo"]["PSEUDO"]
        for element in ("He", "Ne", "Ar", "Kr", "Xe",
                        "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
                        "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"):
            with self.subTest(element=element):
                self.assertIn(element, shipped)
                self.assertIn(element, SSSP_EFFICIENCY)
        # the packaged default and the in-code fallback must not drift apart
        self.assertEqual(sorted(shipped), sorted(SSSP_EFFICIENCY))
        for element, value in SSSP_EFFICIENCY.items():
            with self.subTest(element=element):
                self.assertEqual(shipped[element], value)

    def test_every_vasp_potential_name_is_plausible(self):
        """pseudo.pot names are POTCAR directory names; a typo is only found
        when write_potcar() cannot find the directory, material by material."""
        import json
        import re

        pot = json.loads(
            (PKG / "data" / "config.json").read_text())["pseudo"]["pot"]
        for element, name in pot.items():
            with self.subTest(element=element):
                self.assertTrue(
                    re.fullmatch(r"[A-Z][a-z]?(_[a-z0-9]+)*", name),
                    "%s -> %r is not a POTCAR name" % (element, name))
                self.assertTrue(name.split("_")[0] == element,
                                "%s -> %r names a different element"
                                % (element, name))

class MagneticPseudoKeys(unittest.TestCase):
    """Every magnetic QE input died with
    `PWInputError: Missing Mg,spin=1.0 in pseudo specification!`

    `PWInput.__init__` validates with
        for species in self.structure.composition:
            if str(species) not in pseudo: raise PWInputError(...)
    so for a spin-decorated structure the key is the decorated string.
    FIX(1b) had changed the dict to bare symbols, believing PWInput used
    `site.specie.symbol`. It does not. htesp-tutorials QE/6 was the symptom.
    """

    def test_the_pseudo_dict_is_keyed_the_way_pwinput_looks_it_up(self):
        text = source("htepc")
        self.assertIn("pseudo1[str(species)]", text)
        self.assertIn("for species in self.structure.composition:", text)

    @skip_without("pymatgen")
    def test_pwinput_accepts_what_we_build(self):
        """Build the dict the way setting_qeinput does and hand it to PWInput."""
        from pymatgen.core import Lattice, Structure
        from pymatgen.io.pwscf import PWInput

        from htesp.htepc import bare_element

        structure = Structure(Lattice.cubic(3.2), ["Mg", "Mg"],
                              [[0, 0, 0], [0.5, 0.5, 0.5]])
        structure.add_spin_by_element({"Mg": 1.0})

        pseudo = {}
        for species in structure.composition:
            element = bare_element(getattr(species, "symbol", species))
            pseudo[str(species)] = element + ".upf"
            pseudo.setdefault(element, element + ".upf")

        self.assertIn("Mg,spin=1.0", pseudo)          # the key that was missing
        PWInput(structure, pseudo=pseudo)             # raised PWInputError before


class AfluxRangeOperators(unittest.TestCase):
    """The AFLOW search returned ternaries for a two-element query.

    The AFLUX docs say ',' is OR and ':' is AND inside a property's
    parentheses, which reads as though a bounded range must be
    nspecies(1*:*2). Against the live API that is wrong:

        species('Mg','B'),nspecies(1*,*2)  -> all binary  (B3Mg1, B4Mg2 ...)
        species('Mg','B'),nspecies(1*:*2)  -> all TERNARY (Ag2B1Mg2 ...)

    examples/QE/tutorial5's reference is 64 binary Mg-B entries; the ':'
    spelling produced 64 ternaries. 1.x used ',' and it is correct.
    """

    def test_ranges_use_the_comma_spelling(self):
        text = source("aflow_extract")
        for wanted in ("nspecies(1*,*{}),", "natoms(1*,*{}),",
                       "enthalpy_formation_atom(-100*,*0),"):
            with self.subTest(pattern=wanted):
                self.assertIn(wanted, text)

    def test_the_colon_spelling_is_not_reintroduced(self):
        """It looks more correct than it is; keep it out of the query builder."""
        text = code_only("aflow_extract")
        for forbidden in ("nspecies(1*:*", "natoms(1*:*",
                          "enthalpy_formation_atom(-100*:*"):
            with self.subTest(pattern=forbidden):
                self.assertNotIn(forbidden, text)


class NumpyTwoSpinLabels(unittest.TestCase):
    """`mainprogram magenum` died with
    `ValueError: could not convert string to float: 'np.float64(5.0)'`

    `Species.__str__` renders the spin with `repr`, and numpy 2 changed a
    scalar's repr from `5.0` to `np.float64(5.0)`. The parser did
    `float(item.split('=')[-1])`, which had been fine under numpy 1.
    """

    def test_both_label_spellings_parse(self):
        from htesp.htepc import parse_spin

        for label, wanted in (("Fe,spin=5", 5.0),
                              ("Fe,spin=-5", -5.0),
                              ("Fe,spin=5.0", 5.0),
                              ("Fe,spin=np.float64(5.0)", 5.0),
                              ("Fe,spin=np.float64(-5.0)", -5.0),
                              ("Fe,spin=1e-8", 1e-8),
                              ("Pd", 0.0)):
            with self.subTest(label=label):
                self.assertEqual(parse_spin(label), wanted)

    def test_the_float64_in_the_type_name_is_not_read_as_the_value(self):
        """The first number in 'np.float64(5.0)' is the 64 of the type name."""
        from htesp.htepc import parse_spin

        self.assertNotEqual(parse_spin("Fe,spin=np.float64(5.0)"), 64.0)

    def test_species_lists_come_out_the_same_either_way(self):
        from htesp.htepc import convert_species_list

        plain = convert_species_list(["Fe,spin=5", "Fe,spin=-5", "Pd"])
        numpy2 = convert_species_list(
            ["Fe,spin=np.float64(5.0)", "Fe,spin=np.float64(-5.0)", "Pd"])
        self.assertEqual(plain, numpy2)
        self.assertEqual(plain[0], ["Fe1", "Fe2", "Pd"])


class SummaryFieldRequest(unittest.TestCase):
    """`mainprogram download` raised ValidationError for almost every material.

    `MpConnect.setting()` requested `available_fields[:-29]`, which still
    contains the nested `bandstructure` and `dos` documents. emmet-core
    validates every requested sub-document, and the API's payload for those no
    longer carries the fields the model declares:

        dos.elemental.Ru.total.1.task_id   Field required
        bandstructure.setyawan_curtarolo.equivalent_labels   Field required

    One unusable sub-document rejects the whole SummaryDoc, so only materials
    with neither band structure nor DOS downloaded at all -- 1 of 64 in the
    reported case.
    """

    def test_nested_documents_are_not_requested(self):
        text = source("htepc")
        self.assertIn("UNREQUESTABLE_FIELDS", text)
        from htesp.htepc import UNREQUESTABLE_FIELDS

        self.assertIn("bandstructure", UNREQUESTABLE_FIELDS)
        self.assertIn("dos", UNREQUESTABLE_FIELDS)

    def test_the_request_filters_them_out(self):
        """The slice alone is not enough; the filter has to be applied.

        code_only() tokenises, so match the raw source; this expression
        appears nowhere in a comment or docstring.
        """
        self.assertIn("if field not in UNREQUESTABLE_FIELDS", source("htepc"))

    def test_nothing_reads_the_fields_that_were_dropped(self):
        """Safe to drop only because HTESP computes bands/DOS, never fetches."""
        for module in ("htepc", "element_extract", "crystal"):
            text = code_only(module)
            for pattern in ("data['bandstructure']", 'data["bandstructure"]',
                            "data['dos']", 'data["dos"]'):
                with self.subTest(module=module, pattern=pattern):
                    self.assertNotIn(pattern, text)


class LauncherProcessCountFlag(unittest.TestCase):
    """The run line spelled the process count `-np` whatever the launcher.

    `-np` is mpirun syntax.  srun spells it `-n` and rejects `-np`; ibrun takes
    no count at all and rejects both.  On a TACC or Cray machine every
    generated run-*.sh therefore died on its first line with a launcher usage
    message -- before the DFT code was reached, so the output gave no hint
    where the problem was.
    """

    def test_srun_gets_minus_n(self):
        from htesp.generate_submission import launch

        self.assertEqual(launch("srun", 4), "srun -n 4")

    def test_ibrun_gets_no_count_at_all(self):
        """ibrun runs the whole SLURM allocation; a count is an error."""
        from htesp.generate_submission import launch

        self.assertEqual(launch("ibrun", 4), "ibrun")

    def test_mpirun_is_unchanged(self):
        """The historical behaviour was right for the mpirun family."""
        from htesp.generate_submission import launch

        self.assertEqual(launch("mpirun", 4), "mpirun -np 4")
        self.assertEqual(launch("mpiexec", 4), "mpiexec -np 4")

    def test_an_unknown_launcher_keeps_the_historical_flag(self):
        """Site wrappers are nearly always mpirun-like, and this is what they
        have been running; guessing otherwise would break working setups."""
        from htesp.generate_submission import launch

        self.assertEqual(launch("SiteWrapper", 4), "SiteWrapper -np 4")

    def test_a_full_path_is_matched_on_its_basename(self):
        from htesp.generate_submission import launch

        self.assertEqual(launch("/opt/apps/bin/srun", 8), "/opt/apps/bin/srun -n 8")

    def test_a_launcher_that_already_has_flags_is_left_alone(self):
        """A site that wrote out its launcher means it; appending a second
        count flag would be the same class of bug in reverse."""
        from htesp.generate_submission import launch

        for spelled in ("srun --cpu-bind=cores -n 8", "mpirun -np 8 --bind-to core"):
            with self.subTest(command=spelled):
                self.assertEqual(launch(spelled, 4), spelled)

    def test_no_launcher_runs_the_executable_directly(self):
        """An empty parallel_command used to produce ' -np 4 pw.x ...', which
        the shell read as a command named ''."""
        from htesp.generate_submission import generate_submission_files, launch

        self.assertEqual(launch("", 4), "")
        line = generate_submission_files(
            "qe", "", "4", {"scf": ("pw.x", "scf.in", "scf.out")})["scf"]
        self.assertTrue(line.startswith("pw.x "), line)

    def test_every_branch_uses_the_helper(self):
        """qe, epw, wannier and vasp each built their own run line, so a fix
        applied to one branch would leave the other three broken.

        The FIX comment quotes the old pattern, so comment lines are dropped
        before looking for it.
        """
        text = "\n".join(line for line in source("generate_submission").splitlines()
                         if not line.strip().startswith("#"))
        self.assertNotIn("-np {", text)
        self.assertGreaterEqual(text.count("{run}"), 6)

    def test_the_epw_pools_flag_is_not_a_launcher_flag(self):
        """-npools is passed to epw.x and must survive the rewrite."""
        from htesp.generate_submission import generate_submission_files

        line = generate_submission_files(
            "epw", "ibrun", "4", {"epw": ("epw.x", "epw.in", "epw.out")})["epw"]
        self.assertEqual(line, "ibrun epw.x -npools 4 -i epw.in > epw.out")


if __name__ == "__main__":
    unittest.main()
