"""The package must behave the same on x86_64 and on arm64/aarch64.

Nothing in HTESP is architecture-specific -- no assembly, no ``ctypes``, no
``struct`` packing, no endianness assumption, no ``platform.machine()`` branch,
and every shim is POSIX ``sh``.  These tests keep it that way, because the one
thing that *does* differ between architectures is the compiled wheels
underneath, and a machine-specific workaround in HTESP would hide that rather
than fix it.
"""
from __future__ import annotations

import ast
import re
import subprocess
import unittest
from pathlib import Path

from tests.helpers import python_sources

ROOT = Path(__file__).resolve().parent.parent
PKG = ROOT / "htesp"

#: things that make a source behave differently on one architecture.
#: ``htesp/check.py`` is the deliberate exception: reporting the machine is
#: the whole point of it.
ARCH_TOKENS = (
    "platform.machine", "platform.processor", "platform.architecture",
    "os.uname", "/proc/cpuinfo", "uname -m", "x86_64", "aarch64",
    "arm64", "amd64", "/usr/local/", "/opt/homebrew", "/opt/intel",
    "float128", "longdouble", "longfloat", "float96",
)
ARCH_EXEMPT = {"check.py"}


class NoArchitectureAssumptions(unittest.TestCase):
    def test_no_module_branches_on_the_machine_it_runs_on(self):
        offenders = []
        for path in python_sources(PKG):
            if path.name in ARCH_EXEMPT:
                continue
            text = path.read_text()
            for number, line in enumerate(text.splitlines(), start=1):
                stripped = line.strip()
                if stripped.startswith("#"):
                    continue          # a comment naming an arch is fine
                for token in ARCH_TOKENS:
                    if token in stripped:
                        offenders.append(f"{path.name}:{number} {token}")
        self.assertEqual(offenders, [])

    def test_no_shim_or_script_hardcodes_an_architecture(self):
        offenders = []
        scripts = sorted(p for p in (ROOT / "bin").iterdir() if p.is_file())
        scripts += sorted((ROOT / "tutorials").glob("*.sh"))
        for path in scripts:
            text = path.read_text(errors="replace")
            for number, line in enumerate(text.splitlines(), start=1):
                if line.strip().startswith("#"):
                    continue
                for token in ("x86_64", "aarch64", "arm64", "amd64",
                              "/proc/cpuinfo", "uname -m", "/usr/local/",
                              "/opt/homebrew", "/opt/intel"):
                    if token in line:
                        offenders.append(f"{path.name}:{number} {token}")
        self.assertEqual(offenders, [])

    def test_the_package_is_pure_python(self):
        """A compiled file in the package would make the wheel arch-specific."""
        binaries = sorted(p.name for p in PKG.rglob("*")
                          if p.suffix in (".so", ".pyd", ".dylib", ".c",
                                          ".pyx", ".o", ".a"))
        self.assertEqual(binaries, [])

    def test_no_module_imports_ctypes_or_packs_native_structs(self):
        offenders = []
        for path in python_sources(PKG):
            tree = ast.parse(path.read_text())
            for node in ast.walk(tree):
                roots = []
                if isinstance(node, ast.ImportFrom) and node.level == 0:
                    roots.append((node.module or "").split(".")[0])
                elif isinstance(node, ast.Import):
                    roots += [a.name.split(".")[0] for a in node.names]
                for root in roots:
                    if root in ("ctypes", "cffi", "mmap", "sysconfig"):
                        offenders.append(f"{path.name}:{node.lineno} {root}")
        self.assertEqual(offenders, [])


class PosixShims(unittest.TestCase):
    """The shims run under dash on Linux and bash-as-sh on macOS."""

    def _scripts(self):
        return sorted(p for p in (ROOT / "bin").iterdir() if p.is_file())

    def test_every_sh_shim_parses_as_posix_sh(self):
        for path in self._scripts():
            first = path.read_text(errors="replace").splitlines()[0]
            if "/bin/sh" not in first:
                continue
            with self.subTest(shim=path.name):
                result = subprocess.run(["sh", "-n", str(path)],
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)

    def test_no_sh_shim_uses_a_bash_only_construct(self):
        bashisms = (r"\[\[", r"\$\{[A-Za-z_][A-Za-z0-9_]*,,\}",
                    r"\$\{[A-Za-z_][A-Za-z0-9_]*\^\^\}", r"\bdeclare -A\b",
                    r"\bmapfile\b", r"\breadarray\b", r"\bsource\b", r"&>")
        offenders = []
        for path in self._scripts():
            text = path.read_text(errors="replace")
            if "/bin/sh" not in text.splitlines()[0]:
                continue
            for number, line in enumerate(text.splitlines(), start=1):
                if line.strip().startswith("#"):
                    continue
                for pattern in bashisms:
                    if re.search(pattern, line):
                        offenders.append(f"{path.name}:{number} {pattern}")
        self.assertEqual(offenders, [])


class Check(unittest.TestCase):
    """`htesp-check` is how a new machine gets checked; it must not itself
    need anything that might be broken on that machine."""

    def test_check_imports_with_the_standard_library_alone(self):
        import importlib
        importlib.import_module("htesp.check")

    def test_probe_reports_rather_than_raises(self):
        from htesp.check import probe
        self.assertEqual(probe("json")["status"], "ok")
        self.assertEqual(probe("definitely_not_installed_xyz")["status"],
                         "missing")

    def test_probe_survives_a_dependency_that_aborts_the_interpreter(self):
        """deltalake did exactly this: SIGABRT, which no except can catch."""
        import sys
        import tempfile
        from htesp.check import probe
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "wobbly_probe.py").write_text("import os\nos.abort()\n")
            saved = sys.path[:]
            sys.path.insert(0, tmp)
            try:
                import os as _os
                previous = _os.environ.get("PYTHONPATH")
                _os.environ["PYTHONPATH"] = tmp
                try:
                    result = probe("wobbly_probe")
                finally:
                    if previous is None:
                        _os.environ.pop("PYTHONPATH", None)
                    else:
                        _os.environ["PYTHONPATH"] = previous
            finally:
                sys.path[:] = saved
        self.assertEqual(result["status"], "aborted")
        self.assertIn("SIGABRT", result["detail"])

    def test_the_environment_block_names_the_machine(self):
        from htesp.check import environment
        env = environment()
        for key in ("python", "system", "machine", "byteorder", "cpus"):
            self.assertIn(key, env)
        self.assertIn(env["byteorder"], ("little", "big"))

    def test_mp_api_is_probed_through_the_submodule_that_uses_deltalake(self):
        """`import mp_api` runs only __init__.py and reported a false "ok" on a
        machine where `mp_api.client` aborted."""
        from htesp.check import IMPORT_TARGET
        self.assertEqual(IMPORT_TARGET["mp_api"], "mp_api.client")

    def test_classify_separates_a_page_size_mismatch_from_a_wrong_wheel(self):
        """The first real abort was a *correct* aarch64 wheel whose bundled
        jemalloc assumed 4 KiB pages on a 64 KiB-page kernel.  Telling the user
        to reinstall it would have been useless advice."""
        from htesp.check import classify
        self.assertEqual(classify("<jemalloc>: Unsupported system page size"),
                         "page-size mismatch")
        self.assertEqual(classify("memory allocation of 16 bytes failed"),
                         "page-size mismatch")
        self.assertEqual(classify("incompatible architecture"),
                         "wrong architecture")
        self.assertEqual(classify("Illegal instruction"), "wrong CPU baseline")
        self.assertEqual(classify("something else entirely"), "")

    def test_the_environment_reports_the_memory_page_size(self):
        from htesp.check import environment, page_size
        self.assertGreater(page_size(), 0)
        self.assertEqual(environment()["pagesize"], page_size())

    def test_page_size_advice_names_the_extension_and_the_build_flag(self):
        """The importer (mp_api) is pure Python; the pip command must name the
        extension (deltalake) and the lg_page for this kernel."""
        import contextlib
        import io
        from htesp import check
        aborted = {"status": "aborted", "cause": "page-size mismatch",
                   "detail": "killed by SIGABRT; memory allocation of 16 bytes failed"}
        result = {
            "environment": {"python": "3.11.16", "implementation": "CPython",
                            "executable": "python3.11", "system": "Linux",
                            "release": "5.14.0.aarch64+64k", "machine": "aarch64",
                            "byteorder": "little", "pagesize": 65536, "cpus": 144},
            "required": [dict(aborted, module="mp_api")],
            "transitive": [dict(aborted, module="deltalake")],
            "extras": {}, "executables": {},
        }
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            check.render(result)
        text = buffer.getvalue()
        self.assertIn("JEMALLOC_SYS_WITH_LG_PAGE=16", text)
        self.assertIn("--no-binary deltalake deltalake", text)
        self.assertNotIn("--no-binary mp_api", text)
        self.assertIn("64 KiB", text)

    def test_every_probe_name_maps_to_an_installable_distribution(self):
        """"No module named 'qmpy_rester'" does not tell you to install
        `qmpy-rester`.  Every import name htesp-check probes must resolve to a
        distribution actually declared in pyproject.toml."""
        from htesp.check import EXTRAS, REQUIRED, distribution

        def normalise(name):
            return name.lower().replace("_", "-")

        text = (ROOT / "pyproject.toml").read_text()
        required_block = text.split("dependencies = [", 1)[1].split("]", 1)[0]
        declared = {normalise(name)
                    for name in re.findall(r'"([A-Za-z0-9_.-]+)', required_block)}
        for module in REQUIRED:
            with self.subTest(module=module):
                self.assertIn(normalise(distribution(module)), declared)

        extras_block = text.split("[project.optional-dependencies]", 1)[1]
        extras_block = extras_block.split("[project.urls]", 1)[0]
        for extra, modules in EXTRAS.items():
            line = extras_block.split(f"{extra} = [", 1)[1].split("]", 1)[0]
            names = {normalise(name)
                     for name in re.findall(r'"([A-Za-z0-9_.-]+)', line)}
            for module in modules:
                with self.subTest(extra=extra, module=module):
                    self.assertIn(normalise(distribution(module)), names)

    def test_a_missing_package_is_reported_with_how_to_install_it(self):
        from htesp.check import install_hint
        # required: the distribution name, which is not the import name
        self.assertEqual(install_hint("qmpy_rester"), "pip install qmpy-rester")
        self.assertEqual(install_hint("yaml"), "pip install PyYAML")
        # an extra: the extra, plus the distribution when it differs
        hint = install_hint("sklearn")
        self.assertIn('pip install "htesp[ml]"', hint)
        self.assertIn("scikit-learn", hint)

    def test_a_transitive_abort_is_not_escalated_when_nothing_requires_it(self):
        """After `pip install mp-api==0.45.0` deltalake is no longer imported by
        anything, so an aborting deltalake is a leftover, not a fault."""
        import contextlib
        import io
        from htesp import check
        result = {
            "environment": {"python": "3.11.16", "implementation": "CPython",
                            "executable": "py", "system": "Linux",
                            "release": "5.14.0.aarch64+64k", "machine": "aarch64",
                            "byteorder": "little", "pagesize": 65536, "cpus": 144},
            "required": [{"module": "mp_api", "status": "ok", "detail": "0.45.0"}],
            "transitive": [{"module": "deltalake", "status": "aborted",
                            "cause": "page-size mismatch",
                            "detail": "killed by SIGABRT"}],
            "extras": {}, "executables": {},
        }
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            check.render(result)
        text = buffer.getvalue()
        self.assertIn("Nothing required is broken", text)
        self.assertIn("Every required dependency imports", text)
        self.assertNotIn("JEMALLOC_SYS_WITH_LG_PAGE", text)

    def test_installation_reports_where_htesp_is_running_from(self):
        """"I reinstalled and nothing changed" is usually a second copy on
        sys.path, or the 1.x `HTESP` distribution still present."""
        from htesp.check import installation
        info = installation()
        self.assertTrue(info["package"].endswith("htesp"))
        self.assertIsInstance(info["distributions"], dict)
        self.assertIsInstance(info["warnings"], list)

    def test_the_version_column_holds_versions_not_commands(self):
        """An `ok` row shows the installed version, so a MISSING row must show
        version information too -- what this installation asks for -- and the
        pip commands belong in one block at the end."""
        import contextlib
        import io
        from importlib import metadata
        from unittest import mock
        from htesp import check

        declared = ["mp-api>=0.33", 'matminer>=0.9; extra == "ml"']
        result = {
            "environment": {"python": "3.11", "implementation": "CPython",
                            "executable": "py", "system": "Linux", "release": "r",
                            "machine": "aarch64", "byteorder": "little",
                            "pagesize": 4096, "cpus": 8},
            "installation": {"package": "/x/htesp",
                             "distributions": {"htesp": "2.0.0"}, "warnings": []},
            "required": [{"module": "mp_api", "status": "missing", "detail": "x"},
                         {"module": "numpy", "status": "ok", "detail": "2.4.6"}],
            "transitive": [],
            "extras": {"ml": [{"module": "matminer", "status": "missing",
                               "detail": "x"}]},
            "executables": {},
        }
        buffer = io.StringIO()
        with mock.patch.object(metadata, "requires", return_value=declared):
            check.required_versions.cache_clear()
            try:
                with contextlib.redirect_stdout(buffer):
                    check.render(result)
            finally:
                check.required_versions.cache_clear()
        text = buffer.getvalue()
        rows = [line for line in text.splitlines() if "mp_api" in line
                and line.startswith("  ")]
        self.assertTrue(rows)
        self.assertIn("needs >=0.33", rows[0])
        self.assertNotIn("pip install", rows[0])       # not in the row
        self.assertIn("To install what is missing:", text)
        self.assertIn('pip install "mp-api>=0.33"', text)
        self.assertIn('pip install "htesp[ml]"', text)

    def test_required_versions_parses_the_installed_metadata(self):
        from importlib import metadata
        from unittest import mock
        from htesp import check
        declared = ["numpy>=1.23", "PyYAML>=6.0", "mp_api>=0.33",
                    'scikit-learn>=1.2; extra == "ml"']
        with mock.patch.object(metadata, "requires", return_value=declared):
            check.required_versions.cache_clear()
            try:
                self.assertEqual(check.wanted_version("mp_api"), ">=0.33")
                self.assertEqual(check.wanted_version("yaml"), ">=6.0")
                self.assertEqual(check.wanted_version("sklearn"), ">=1.2")
                self.assertEqual(check.wanted_version("nothing_here"), "")
            finally:
                check.required_versions.cache_clear()

    def test_installed_but_broken_is_not_reported_as_missing(self):
        """`from mp_api.client import ...` raised `ImportError: cannot import
        name 'BSPathType'` -- mp-api was installed, emmet-core was too new.
        Reporting that as MISSING would have advised reinstalling mp-api,
        which fixes nothing."""
        import os
        import tempfile
        from htesp.check import probe
        with tempfile.TemporaryDirectory() as tmp:
            Path(tmp, "brokenprobe.py").write_text(
                "from json import ThisNameDoesNotExist\n")
            Path(tmp, "needsdepprobe.py").write_text("import absent_dep_xyz\n")
            previous = os.environ.get("PYTHONPATH")
            os.environ["PYTHONPATH"] = tmp
            try:
                broken = probe("brokenprobe")
                needs = probe("needsdepprobe")
                absent = probe("no_such_module_at_all_xyz")
            finally:
                if previous is None:
                    os.environ.pop("PYTHONPATH", None)
                else:
                    os.environ["PYTHONPATH"] = previous
        self.assertEqual(broken["status"], "broken")
        self.assertIn("ThisNameDoesNotExist", broken["detail"])
        self.assertEqual(needs["status"], "broken")      # its dependency, not it
        self.assertIn("absent_dep_xyz", needs["detail"])
        self.assertEqual(absent["status"], "missing")    # genuinely not there
        self.assertIn("install", absent)

    def test_a_broken_package_is_not_offered_for_installation(self):
        import contextlib
        import io
        from htesp import check
        result = {
            "environment": {"python": "3.11", "implementation": "CPython",
                            "executable": "py", "system": "Linux", "release": "r",
                            "machine": "aarch64", "byteorder": "little",
                            "pagesize": 65536, "cpus": 144},
            "installation": {"package": "/x/htesp", "distributions": {},
                             "warnings": []},
            "required": [{"module": "mp_api", "status": "broken",
                          "detail": "ImportError: cannot import name 'BSPathType'"}],
            "transitive": [], "extras": {}, "executables": {},
        }
        buffer = io.StringIO()
        with contextlib.redirect_stdout(buffer):
            check.render(result)
        text = buffer.getvalue()
        self.assertIn("version conflict", text)
        self.assertIn("BSPathType", text)
        self.assertNotIn("To install what is missing", text)

    def test_required_list_matches_pyproject(self):
        """htesp-check must probe exactly the declared dependencies."""
        from htesp.check import REQUIRED
        text = (ROOT / "pyproject.toml").read_text()
        block = text.split("dependencies = [", 1)[1].split("]", 1)[0]
        declared = re.findall(r'"([A-Za-z0-9_.-]+)', block)
        #: distribution name -> import name, where they differ
        as_import = {"PyYAML": "yaml", "mp-api": "mp_api",
                     "qmpy-rester": "qmpy_rester"}
        expected = {as_import.get(name, name).lower() for name in declared}
        self.assertEqual({name.lower() for name in REQUIRED}, expected)


class SetMpApiKey(unittest.TestCase):
    """`htesp-check --set_mp_api KEY` stores the key where it survives.

    `$MP_API_KEY` lives only in the shell that exported it, so a nohup-ed
    sweep or a new terminal loses it and every database command starts
    skipping. The credentials file survives, and keeps the key out of
    config.json -- which is how a real key ended up in 218 tracked files.
    """

    def test_the_option_exists(self):
        from htesp.check import set_mp_api_key  # noqa: F401

        text = (ROOT / "htesp" / "check.py").read_text()
        self.assertIn('"--set_mp_api"', text)

    def test_an_empty_key_is_an_error_not_a_fallthrough(self):
        """`if args.set_mp_api:` treated '' as absent and printed the report."""
        from htesp.check import set_mp_api_key

        self.assertEqual(set_mp_api_key("", verify=False), 1)
        text = (ROOT / "htesp" / "check.py").read_text()
        self.assertIn("if args.set_mp_api is not None:", text)

    def test_the_shipped_placeholder_is_rejected(self):
        from htesp.check import set_mp_api_key

        self.assertEqual(set_mp_api_key("use_your_API_KEY", verify=False), 1)

    def test_a_rejected_key_is_never_written(self):
        """The first version wrote the file and then asked MP, so a typo
        replaced a working key with a broken one."""
        text = (ROOT / "htesp" / "check.py").read_text()
        verify_at = text.index("the Materials Project rejected it")
        write_at = text.index("path.write_text(body)")
        self.assertLess(verify_at, write_at,
                        "the key must be verified before it is stored")
        self.assertIn("the existing credentials are untouched", text)

    def test_the_key_is_never_printed_in_full(self):
        from htesp.check import _mask

        self.assertEqual(_mask("abcdefghijklmnopqrstuvwxyz"), "abcd...wxyz")
        self.assertNotIn("efghijklmnopqrstuv", _mask("abcdefghijklmnopqrstuvwxyz"))
        self.assertEqual(_mask("short"), "*****")

    def test_it_writes_the_name_config_py_reads(self):
        """config.py accepts MP_API_KEY, api_key or key; write the first."""
        from htesp.config import API_KEY_ENV

        text = (ROOT / "htesp" / "check.py").read_text()
        self.assertIn('"{}={}".format(API_KEY_ENV, key)', text)
        self.assertEqual(API_KEY_ENV, "MP_API_KEY")


class EnumlibIsReported(unittest.TestCase):
    """A successful `--install-enumlib` still looked like it had done nothing.

    `htesp-check --executables` never probed for enum.x or makestr.x, so the
    binaries could be built, installed and on PATH and the report stayed
    silent about them.
    """

    def test_the_enumlib_binaries_are_probed(self):
        from htesp.check import ENUMLIB_BINARIES, EXECUTABLES

        for name in ENUMLIB_BINARIES:
            with self.subTest(name=name):
                self.assertIn(name, EXECUTABLES,
                              "--executables does not probe for %s" % name)

    def test_what_is_installed_is_what_is_reported(self):
        """The list the installer writes and the list the report reads are
        the same object, so they cannot drift apart."""
        text = (ROOT / "htesp" / "check.py").read_text()
        self.assertIn('+ ENUMLIB_BINARIES', text)
        # and only one definition of it
        self.assertEqual(text.count("ENUMLIB_BINARIES = ("), 1)


class CleanSourceTree(unittest.TestCase):
    """`htesp-check --clean` removes build artifacts, and nothing else."""

    def test_it_refuses_an_installed_copy(self):
        from htesp.check import clean_source_tree

        self.assertEqual(
            clean_source_tree("/usr/lib/python3/site-packages"), 1)

    def test_it_refuses_a_directory_that_is_not_a_checkout(self):
        import tempfile

        from htesp.check import clean_source_tree

        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(clean_source_tree(tmp), 1)

    def test_the_deliberate_archive_is_never_swept(self):
        """_removed/ holds .pyc files the 2.0 work moved there on purpose;
        treating them as build output destroys the record."""
        from htesp.check import CLEAN_EXCLUDE

        self.assertIn("_removed", CLEAN_EXCLUDE)
        self.assertIn(".git", CLEAN_EXCLUDE)

    def test_dry_run_removes_nothing(self):
        import tempfile
        from pathlib import Path as _Path

        from htesp.check import clean_source_tree

        with tempfile.TemporaryDirectory() as tmp:
            root = _Path(tmp)
            (root / "pyproject.toml").write_text("[project]\n")
            (root / "htesp").mkdir()
            (root / "build").mkdir()
            (root / "_removed").mkdir()
            (root / "_removed" / "kept.pyc").write_bytes(b"x")
            self.assertEqual(clean_source_tree(str(root), dry_run=True), 0)
            self.assertTrue((root / "build").is_dir(), "dry run deleted something")

    def test_it_removes_build_output_but_spares_the_archive(self):
        import tempfile
        from pathlib import Path as _Path

        from htesp.check import clean_source_tree

        with tempfile.TemporaryDirectory() as tmp:
            root = _Path(tmp)
            (root / "pyproject.toml").write_text("[project]\n")
            (root / "htesp" / "__pycache__").mkdir(parents=True)
            (root / "htesp" / "__pycache__" / "a.pyc").write_bytes(b"x")
            (root / "build").mkdir()
            (root / "htesp.egg-info").mkdir()
            (root / "_removed").mkdir()
            (root / "_removed" / "kept.pyc").write_bytes(b"x")
            (root / "htesp" / "real_source.py").write_text("x = 1\n")

            self.assertEqual(clean_source_tree(str(root)), 0)
            self.assertFalse((root / "build").exists())
            self.assertFalse((root / "htesp.egg-info").exists())
            self.assertFalse((root / "htesp" / "__pycache__").exists())
            self.assertTrue((root / "_removed" / "kept.pyc").is_file(),
                            "the deliberate archive was swept")
            self.assertTrue((root / "htesp" / "real_source.py").is_file(),
                            "a source file was removed")

    def test_it_leaves_enumlib_and_configuration_alone(self):
        """--clean is about this tree, not about the machine."""
        import ast

        source = (ROOT / "htesp" / "check.py").read_text()
        func = next(n for n in ast.walk(ast.parse(source))
                    if isinstance(n, ast.FunctionDef) and n.name == "clean_source_tree")
        if ast.get_docstring(func):
            func.body = func.body[1:]      # the docstring names them to say it does NOT
        body = ast.unparse(func)
        for untouched in ("ENUMLIB", "CREDENTIALS", "PMG_VASP_PSP_DIR"):
            with self.subTest(untouched=untouched):
                self.assertNotIn(untouched, body)


if __name__ == "__main__":
    unittest.main()
