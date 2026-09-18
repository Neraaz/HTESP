"""Configuration loading -- the behaviour that replaced ``check_json.config``."""
from __future__ import annotations

import json
import os
import unittest
from pathlib import Path

from tests.helpers import TempProject

from htesp import config as cfgmod


class ConfigSearch(TempProject):
    def test_packaged_default_is_returned_when_nothing_is_found(self):
        # the old loader returned None here and every caller raised TypeError
        cfg = cfgmod.config(self.root)
        self.assertIsInstance(cfg, dict)
        self.assertIn("pseudo", cfg)
        self.assertIsNone(cfgmod.config_path(self.root))

    def test_config_in_the_working_directory_wins(self):
        self.write("config.json", json.dumps({"kptden": 0.111}))
        cfg = cfgmod.config(self.root)
        self.assertEqual(cfg["kptden"], 0.111)
        self.assertEqual(cfgmod.config_path(self.root), self.root / "config.json")

    def test_parent_directories_are_searched(self):
        self.write("config.json", json.dumps({"kptden": 0.222}))
        deep = self.root / "Rmp-763-Mg1B2" / "relax"
        deep.mkdir(parents=True)
        self.assertEqual(cfgmod.config(deep)["kptden"], 0.222)

    def test_search_stops_after_the_documented_depth(self):
        self.write("config.json", json.dumps({"kptden": 0.333}))
        deep = self.root.joinpath(*["level"] * (cfgmod.SEARCH_DEPTH + 2))
        deep.mkdir(parents=True)
        self.assertIsNone(cfgmod.find_config(deep))

    def test_htesp_config_environment_variable_wins(self):
        self.write("config.json", json.dumps({"kptden": 0.1}))
        other = self.write("elsewhere/config.json", json.dumps({"kptden": 0.9}))
        os.environ["HTESP_CONFIG"] = str(other)
        cfgmod.clear_cache()
        self.assertEqual(cfgmod.config(self.root)["kptden"], 0.9)


class ConfigMerge(TempProject):
    def test_a_partial_config_still_has_every_key(self):
        """An example config written for the 2023 schema is missing 54 keys."""
        self.write("config.json", json.dumps({"download": {"mode": "chemsys"}}))
        cfg = cfgmod.config(self.root)
        self.assertEqual(cfg["download"]["mode"], "chemsys")
        for key in ("chull_cutoff", "kpt_opt", "elph_mode", "wanniertools_input"):
            self.assertIn(key, cfg, f"{key} should come from the packaged default")
        self.assertIn("calc", cfg["download"]["inp"])

    def test_merge_is_deep_not_shallow(self):
        self.write("config.json", json.dumps({"pseudo": {"pot": "mine"}}))
        cfg = cfgmod.config(self.root)
        self.assertEqual(cfg["pseudo"]["pot"], "mine")
        self.assertIn("PSEUDO", cfg["pseudo"], "a shallow merge would drop this")

    def test_the_default_table_covers_the_elements_the_fallback_knew(self):
        """O, Se and Au were in the hard-coded fallback but not in config.json."""
        table = cfgmod.config(self.root)["pseudo"]["PSEUDO"]
        for element in ("O", "Se", "Au", "Mg", "B"):
            self.assertIn(element, table)

    def test_returned_dict_is_a_copy(self):
        first = cfgmod.config(self.root)
        first["kptden"] = "mutated"
        self.assertNotEqual(cfgmod.config(self.root)["kptden"], "mutated")


class ApiKey(TempProject):
    def test_placeholder_counts_as_absent(self):
        self.assertIsNone(cfgmod.api_key(cfgmod.config(self.root)))

    def test_environment_variable_wins_over_the_file(self):
        self.write("config.json", json.dumps(
            {"mpi_key": {"API_KEY": {"key": "from-the-file"}}}))
        os.environ["MP_API_KEY"] = "from-the-environment"
        cfgmod.clear_cache()
        self.assertEqual(cfgmod.api_key(cfgmod.config(self.root)),
                         "from-the-environment")

    def test_the_file_is_still_honoured(self):
        self.write("config.json", json.dumps(
            {"mpi_key": {"API_KEY": {"key": "from-the-file"}}}))
        self.assertEqual(cfgmod.api_key(cfgmod.config(self.root)), "from-the-file")

    def test_require_api_key_explains_how_to_set_it(self):
        with self.assertRaises(cfgmod.ConfigError) as caught:
            cfgmod.require_api_key(cfgmod.config(self.root))
        self.assertIn("MP_API_KEY", str(caught.exception))


class Validation(TempProject):
    def test_the_packaged_default_only_lacks_a_key(self):
        problems = cfgmod.validate(cfgmod.config(self.root))
        self.assertEqual(len(problems), 1)
        self.assertIn("API key", problems[0])

    def test_a_bad_calc_is_reported(self):
        self.write("config.json", json.dumps({"download": {"inp": {"calc": "ABINIT"}}}))
        problems = cfgmod.validate(cfgmod.config(self.root))
        self.assertTrue(any("calc" in p for p in problems), problems)

    def test_a_bad_plot_limit_is_reported(self):
        self.write("config.json", json.dumps({"plot": {"xlim": [1, 2, 3]}}))
        problems = cfgmod.validate(cfgmod.config(self.root))
        self.assertTrue(any("plot.xlim" in p for p in problems), problems)


class PackagedDefault(unittest.TestCase):
    def test_it_is_valid_json_and_carries_no_key(self):
        text = Path(cfgmod.DEFAULT_CONFIG_PATH).read_text()
        data = json.loads(text)
        self.assertEqual(data["mpi_key"]["API_KEY"]["key"],
                         cfgmod.API_KEY_PLACEHOLDER)
        for section in cfgmod.REQUIRED_SECTIONS:
            self.assertIn(section, data)

    def test_the_wanniertools_key_is_spelled_correctly(self):
        """`wanniertool_input` (no s) was a KeyError at import time."""
        data = json.loads(Path(cfgmod.DEFAULT_CONFIG_PATH).read_text())
        self.assertIn("wanniertools_input", data)
        self.assertNotIn("wanniertool_input", data)


if __name__ == "__main__":
    unittest.main()


class WhichConfigWasUsed(unittest.TestCase):
    """"Which config.json did that run read?" must be answerable after the
    fact, not only by re-deriving the search."""

    def test_config_validate_names_the_file_and_the_search(self):
        import io
        import contextlib
        from htesp import cli
        import tempfile
        import os
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "config.json").write_text('{"kptden": 0.5}')
            saved = os.getcwd()
            os.chdir(root)
            try:
                from htesp import config as config_module
                config_module.clear_cache()
                args = cli.build_parser().parse_args(["config-validate"])
                ctx = cli.Context(args)
                buffer = io.StringIO()
                with contextlib.redirect_stdout(buffer):
                    cli.cmd_config_validate(ctx, [])
            finally:
                os.chdir(saved)
                from htesp import config as config_module
                config_module.clear_cache()
        text = buffer.getvalue()
        self.assertIn("configuration:", text)
        self.assertIn("config.json", text)
        self.assertIn("searched:", text)
        self.assertIn("merged over:", text)

    def test_the_log_file_records_the_configuration(self):
        import os
        import tempfile
        from htesp import cli
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "config.json").write_text('{"kptden": 0.5}')
            (root / "input.in").write_text(
                "1\n2\n200 0\nmpid.in\nphband\nDFT = QE\n")
            (root / "mpid.in").write_text("v1 mp-763 Mg1B2\n")
            saved = os.getcwd()
            os.chdir(root)
            try:
                from htesp import config as config_module
                config_module.clear_cache()
                args = cli.build_parser().parse_args(["4"])
                ctx = cli.Context(args)
                ctx.write_log()
                written = (root / "log").read_text()
            finally:
                os.chdir(saved)
                from htesp import config as config_module
                config_module.clear_cache()
        self.assertIn("# configuration: ", written)
        self.assertIn("config.json", written)


class ConfigInit(unittest.TestCase):
    """Without a config.json of its own a project runs silently on the packaged
    default: nothing fails, but the cutoffs are whatever the package ships."""

    def _run(self, argv, cwd):
        import contextlib
        import io
        import os
        from htesp import cli
        from htesp import config as config_module
        saved = os.getcwd()
        os.chdir(cwd)
        try:
            config_module.clear_cache()
            args = cli.build_parser().parse_args(argv)
            ctx = cli.Context(args)
            buffer = io.StringIO()
            with contextlib.redirect_stdout(buffer):
                status = cli.SPECIAL_COMMANDS[args.process](ctx, list(args.rest))
            return status, buffer.getvalue()
        finally:
            os.chdir(saved)
            config_module.clear_cache()

    def test_it_writes_the_packaged_default_and_validate_then_finds_it(self):
        import json
        import tempfile
        from htesp.config import DEFAULT_CONFIG_PATH
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            status, out = self._run(["config-init"], root)
            self.assertEqual(status, 0)
            written = root / "config.json"
            self.assertTrue(written.is_file())
            self.assertEqual(json.loads(written.read_text()),
                             json.loads(DEFAULT_CONFIG_PATH.read_text()))
            self.assertIn("MP_API_KEY", out)         # the key is not in the file
            _status, validated = self._run(["config-validate"], root)
            self.assertIn(str(written), validated)

    def test_it_refuses_to_overwrite_without_force(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "config.json").write_text('{"mine": true}')
            status, _out = self._run(["config-init"], root)
            self.assertEqual(status, 2)
            self.assertEqual((root / "config.json").read_text(), '{"mine": true}')
            status, _out = self._run(["config-init", "--force"], root)
            self.assertEqual(status, 0)
            self.assertNotEqual((root / "config.json").read_text(), '{"mine": true}')

    def test_it_accepts_an_explicit_target_and_makes_parents(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            status, _out = self._run(["config-init", "sub/dir/my.json"], root)
            self.assertEqual(status, 0)
            self.assertTrue((root / "sub" / "dir" / "my.json").is_file())

    def test_the_written_file_carries_no_api_key(self):
        import json
        import tempfile
        from htesp.config import API_KEY_PLACEHOLDER
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._run(["config-init"], root)
            data = json.loads((root / "config.json").read_text())
            self.assertEqual(data["mpi_key"]["API_KEY"]["key"],
                             API_KEY_PLACEHOLDER)
