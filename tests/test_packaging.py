"""Packaging, shims and hygiene: the things that made `pip install .` unusable."""
from __future__ import annotations

import ast
import json
import re
import subprocess
import sys
import unittest
from pathlib import Path

from tests.helpers import python_sources, stray_sidecars

ROOT = Path(__file__).resolve().parent.parent
PKG = ROOT / "htesp"


def _normalise(name: str) -> str:
    """PyPI treats `mp_api`, `mp-api` and `MP-API` as the same project."""
    return name.lower().replace("_", "-")


def _declared_dependencies() -> set:
    """The names in `[project] dependencies` of pyproject.toml."""
    text = (ROOT / "pyproject.toml").read_text()
    block = text.split("dependencies = [", 1)[1].split("]", 1)[0]
    return {_normalise(name)
            for name in re.findall(r'"([A-Za-z0-9_.-]+)', block)}


def _declared_extras() -> set:
    """The extra names in `[project.optional-dependencies]` of pyproject.toml."""
    text = (ROOT / "pyproject.toml").read_text()
    block = text.split("[project.optional-dependencies]", 1)[1]
    block = block.split("[project.urls]", 1)[0]
    return set(re.findall(r"^([A-Za-z0-9_.-]+)\s*=\s*\[", block, re.MULTILINE))


def _pinned_names(path) -> set:
    """The package names in a requirements file, ignoring versions."""
    names = set()
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        name = re.split(r"[=<>!~\[]", line, maxsplit=1)[0].strip()
        if name:
            names.add(_normalise(name))
    return names


class PackageLayout(unittest.TestCase):
    def test_every_module_compiles(self):
        result = subprocess.run(
            [sys.executable, "-m", "py_compile", *[str(p) for p in python_sources(PKG)]],
            capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_no_module_imports_a_flat_sibling(self):
        """`from kpath import ...` only worked if the cwd happened to be src/."""
        siblings = {p.stem for p in python_sources(PKG)} - {"__init__"}
        offenders = []
        for path in python_sources(PKG):
            tree = ast.parse(path.read_text())
            for node in ast.walk(tree):
                if isinstance(node, ast.ImportFrom) and node.level == 0:
                    root = (node.module or "").split(".")[0]
                    if root in siblings:
                        offenders.append(f"{path.name}: from {node.module}")
                elif isinstance(node, ast.Import):
                    for alias in node.names:
                        if alias.name.split(".")[0] in siblings:
                            offenders.append(f"{path.name}: import {alias.name}")
        self.assertEqual(offenders, [])

    def test_the_default_config_is_shipped_inside_the_package(self):
        data = ROOT / "htesp" / "data" / "config.json"
        self.assertTrue(data.is_file())
        json.loads(data.read_text())

    def test_entry_points_resolve(self):
        text = (ROOT / "pyproject.toml").read_text()
        for target in ("htesp.cli:main", "htesp.workflow:main",
                       "tutorials.run_tutorials:main", "htesp.check:main"):
            self.assertIn(target, text)
        module, _, function = "htesp.cli:main".partition(":")
        __import__(module)
        self.assertTrue(callable(getattr(sys.modules[module], function)))


class Shims(unittest.TestCase):
    def test_every_legacy_script_has_a_shim(self):
        legacy = ROOT / "legacy" / "bash"
        if not legacy.is_dir():
            self.skipTest("legacy/bash not present")
        legacy_names = {p.name for p in legacy.iterdir() if p.is_file()}
        shims = {p.name for p in (ROOT / "bin").iterdir() if p.is_file()}
        self.assertEqual(legacy_names - shims, set())

    def test_every_shim_is_executable_and_parses(self):
        for shim in sorted((ROOT / "bin").iterdir()):
            with self.subTest(shim=shim.name):
                self.assertTrue(shim.stat().st_mode & 0o111, "not executable")
                result = subprocess.run(["bash", "-n", str(shim)],
                                        capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)

    def test_every_shim_is_listed_in_pyproject(self):
        text = (ROOT / "pyproject.toml").read_text()
        block = text[text.index("script-files = ["):]
        block = block[:block.index("]")]
        listed = set(re.findall(r'"bin/([^"]+)"', block))
        actual = {p.name for p in (ROOT / "bin").iterdir() if p.is_file()}
        self.assertEqual(actual - listed, set(), "shim not listed in pyproject.toml")
        self.assertEqual(listed - actual, set(), "pyproject.toml lists a missing shim")

    def test_shims_delegate_to_the_python_workflow(self):
        sample = (ROOT / "bin" / "create-inputs").read_text()
        self.assertIn("htesp.workflow", sample)
        self.assertNotIn("sbatch", sample)


class Hygiene(unittest.TestCase):
    """Regressions on the whole-repo problems the review found."""

    def test_no_committed_api_key(self):
        """A real-format MP key sat in 218 tracked JSON files."""
        pattern = re.compile(r'"key"\s*:\s*"([A-Za-z0-9]{28,40})"')
        offenders = []
        for path in ROOT.rglob("*.json"):
            if "_removed" in path.parts or ".git" in path.parts:
                continue
            try:
                text = path.read_text(errors="ignore")
            except OSError:
                continue
            for match in pattern.finditer(text):
                if match.group(1) != "use_your_API_KEY":
                    offenders.append(str(path.relative_to(ROOT)))
        self.assertEqual(offenders[:10], [])

    def test_no_os_system_calls_in_the_package(self):
        """310 unchecked os.system calls: a failed stage was invisible."""
        offenders = []
        for path in python_sources(PKG):
            tree = ast.parse(path.read_text())
            for node in ast.walk(tree):
                if (isinstance(node, ast.Call)
                        and isinstance(node.func, ast.Attribute)
                        and node.func.attr == "system"
                        and isinstance(node.func.value, ast.Name)
                        and node.func.value.id == "os"):
                    offenders.append(f"{path.name}:{node.lineno}")
        self.assertEqual(offenders, [])

    def test_no_bare_except_in_the_package(self):
        offenders = []
        for path in python_sources(PKG):
            tree = ast.parse(path.read_text())
            for node in ast.walk(tree):
                if isinstance(node, ast.ExceptHandler) and node.type is None:
                    offenders.append(f"{path.name}:{node.lineno}")
        self.assertEqual(offenders, [])

    def test_no_eval_or_exec(self):
        offenders = []
        for path in python_sources(PKG):
            tree = ast.parse(path.read_text())
            for node in ast.walk(tree):
                if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                        and node.func.id in ("eval", "exec")):
                    offenders.append(f"{path.name}:{node.lineno}")
        self.assertEqual(offenders, [])

    def test_no_shell_true_subprocesses(self):
        offenders = []
        for path in python_sources(PKG) + python_sources(ROOT / "tutorials"):
            tree = ast.parse(path.read_text())
            for node in ast.walk(tree):
                if isinstance(node, ast.Call):
                    for keyword in node.keywords:
                        if (keyword.arg == "shell"
                                and isinstance(keyword.value, ast.Constant)
                                and keyword.value.value is True):
                            offenders.append(f"{path.name}:{node.lineno}")
        self.assertEqual(offenders, [])


    def test_no_appledouble_sidecars(self):
        """`._module.py` blobs from copying the tree off macOS break every scan.

        scp/rsync/tar from macOS to a non-HFS filesystem writes one `._x`
        beside each `x`, holding the extended attributes.  They match
        `*.py`, so `ast.parse` sees null bytes and `py_compile` refuses them.
        Copy with `rsync -a --exclude='._*'` or `COPYFILE_DISABLE=1 tar`, and
        clear the existing ones with `find . -name '._*' -delete`.
        """
        strays = [str(p.relative_to(ROOT)) for p in stray_sidecars(ROOT)]
        self.assertEqual(strays[:20], [])

    def test_dependency_names_are_distributions_not_modules(self):
        """`mp_api` is the module; the distribution on PyPI is `mp-api`. Both
        install, because pip normalises the name, but writing the module name
        invites the mistake that `pymatgen-core` already made in reverse --
        a plausible-looking name that is not the package you meant."""
        text = (ROOT / "pyproject.toml").read_text()
        blocks = [text.split("dependencies = [", 1)[1].split("]", 1)[0]]
        extras = text.split("[project.optional-dependencies]", 1)[1]
        blocks.append(extras.split("[project.urls]", 1)[0])
        offenders = []
        for block in blocks:
            for name in re.findall(r'"([A-Za-z0-9_.-]+)', block):
                if name.startswith("htesp["):
                    continue
                if "_" in name:
                    offenders.append(name)
        self.assertEqual(offenders, [],
                         "use the PyPI distribution name, with hyphens")

    def test_the_requirements_files_use_distribution_names_too(self):
        """The same rule, applied to the two requirements files."""
        for relative in ("requirements.txt", "INSTALL/requirements1.txt"):
            raw = [re.split(r"[=<>!~]", line.strip(), maxsplit=1)[0].strip()
                   for line in (ROOT / relative).read_text().splitlines()
                   if line.strip() and not line.startswith("#")]
            with self.subTest(file=relative):
                self.assertEqual([name for name in raw if "_" in name], [],
                                 "use the PyPI distribution name, with hyphens")

    #: pinned although no extra declares them, each explained in the file
    TRANSITIVE_PINS = {"emmet-core"}

    def test_requirements_txt_covers_the_required_set_and_no_extras(self):
        """It drifted the other way once: `qmpy-rester` sat in
        `INSTALL/requirements1.txt` and not here, which reads like an omission
        until you know the extras are deliberately excluded. Pin the rule."""
        declared = _declared_dependencies()
        listed = _pinned_names(ROOT / "requirements.txt")
        self.assertEqual(sorted(declared - listed), [], "required package absent")
        self.assertEqual(sorted(listed - declared - self.TRANSITIVE_PINS), [],
                         "an extra, or an undocumented transitive pin")

    def test_both_requirements_files_are_pinned_and_agree(self):
        """They are meant to describe one environment; a version that differs
        between them means one of the two was edited alone."""
        def pins(relative):
            out = {}
            for line in (ROOT / relative).read_text().splitlines():
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                self.assertIn("==", line, f"{relative}: {line} is not pinned")
                name, _, version = line.partition("==")
                out[_normalise(name)] = version
            return out
        first = pins("requirements.txt")
        second = pins("INSTALL/requirements1.txt")
        disagree = {name: (first[name], second[name])
                    for name in first.keys() & second.keys()
                    if first[name] != second[name]}
        self.assertEqual(disagree, {})

    def test_the_pinned_environment_covers_every_required_package(self):
        """requirements1.txt pinned ifermi and qmpy-rester (both extras) while
        omitting spglib and PyYAML (both required), so the environment it
        produced could not import `crystal` or `displace_phonopy`."""
        pinned = _pinned_names(ROOT / "INSTALL" / "requirements1.txt")
        self.assertEqual(sorted(_declared_dependencies() - pinned), [])

    def test_everything_pinned_beyond_the_required_set_is_declared_somewhere(self):
        """A pin for a package no extra provides is an untracked dependency."""
        text = (ROOT / "pyproject.toml").read_text()
        extras = text.split("[project.optional-dependencies]", 1)[1]
        extras = extras.split("[project.urls]", 1)[0]
        known = _declared_dependencies() | {
            _normalise(name) for name in re.findall(r'"([A-Za-z0-9_.-]+)', extras)}
        known.add("emmet-core")     # transitive pin of mp-api, documented as such
        extra_pins = _pinned_names(ROOT / "INSTALL" / "requirements1.txt") - known
        self.assertEqual(sorted(extra_pins), [])

    def test_every_pin_is_exact(self):
        for relative in ("requirements.txt", "INSTALL/requirements1.txt"):
            text = (ROOT / relative).read_text()
            loose = [line.strip() for line in text.splitlines()
                     if line.strip() and not line.startswith("#")
                     and "==" not in line]
            with self.subTest(file=relative):
                self.assertEqual(loose, [])

    def test_no_document_offers_an_extra_that_does_not_exist(self):
        """`htesp[oqmd]` was offered by README.md, docs/usage.rst and a comment
        in check.py, but no such extra was ever declared: `qmpy-rester` is a
        required dependency. `pip install ".[oqmd]"` fails outright, so the
        one instruction a reader follows when `import qmpy_rester` fails was
        the one that could not work."""
        declared = _declared_extras()
        sources = [ROOT / "README.md", ROOT / "CHANGELOG.md",
                   ROOT / "INSTALL" / "README"]
        sources += sorted((ROOT / "docs").glob("*.rst"))
        for directory in (PKG, ROOT / "tutorials", ROOT / "tools"):
            sources += python_sources(directory)

        # Two independent passes: one line can carry both spellings, and a
        # single alternation lets the `pip install ...` branch swallow an
        # `htesp[...]` that appears earlier on the same line.
        patterns = (re.compile(r"htesp\[([^\]]+)\]"),
                    re.compile(r"pip install[^\n]*?\.\[([^\]]+)\]"))
        offenders = []
        for path in sources:
            for number, line in enumerate(path.read_text().splitlines(), 1):
                groups = [match.group(1) for pattern in patterns
                          for match in pattern.finditer(line)]
                for group in groups:
                    for name in group.split(","):
                        name = name.strip().strip("\"'`")
                        # `htesp[{extra}]` in an f-string is a placeholder;
                        # the values it takes come from EXTRAS, checked in
                        # tests/test_portability.py.
                        if not name or not re.fullmatch(r"[A-Za-z0-9_.-]+", name):
                            continue
                        if name not in declared:
                            offenders.append(
                                f"{path.relative_to(ROOT)}:{number}: "
                                f"htesp[{name}] is not declared in pyproject.toml")
        self.assertEqual(offenders, [])

    def test_gitignore_covers_the_generated_artefacts(self):
        text = (ROOT / ".gitignore").read_text()
        for pattern in ("__pycache__/", "*.egg-info/", ".DS_Store", "._*", "log",
                        "slurm-*.out", "result.csv", "econv.csv"):
            self.assertIn(pattern, text, pattern)


class OneDefinitionOfAPythonSource(unittest.TestCase):
    """`tools/check_names.py` died on AppleDouble sidecars; the suite did not.

    Both enumerate Python files, but the tool globbed `*.py` itself while the
    suite went through `tests.helpers.python_sources`. One copy from a Mac and
    the tool was unusable -- `UnicodeDecodeError` before it scanned anything --
    which is how it stayed broken for a whole session.
    """

    def test_the_tool_uses_the_shared_enumerator(self):
        text = (ROOT / "tools" / "check_names.py").read_text()
        self.assertIn("from tests.helpers import python_sources", text)
        self.assertIn("python_sources(path, recursive=True)", text)
        self.assertNotIn('path.rglob("*.py")', text)

    def test_python_sources_can_recurse(self):
        import tempfile
        from pathlib import Path as _Path

        from tests.helpers import python_sources

        with tempfile.TemporaryDirectory() as tmp:
            root = _Path(tmp)
            (root / "pkg").mkdir()
            (root / "top.py").write_text("x = 1\n")
            (root / "pkg" / "deep.py").write_text("y = 2\n")
            flat = [p.name for p in python_sources(root)]
            deep = [p.name for p in python_sources(root, recursive=True)]
            self.assertEqual(flat, ["top.py"])
            self.assertEqual(deep, ["deep.py", "top.py"])

    def test_a_sidecar_is_excluded_either_way(self):
        import tempfile
        from pathlib import Path as _Path

        from tests.helpers import python_sources

        with tempfile.TemporaryDirectory() as tmp:
            root = _Path(tmp)
            (root / "real.py").write_text("x = 1\n")
            (root / "._real.py").write_bytes(b"\x00\x05\x16\x07Mac OS X\x00")
            for recursive in (False, True):
                with self.subTest(recursive=recursive):
                    names = [p.name for p in python_sources(root, recursive=recursive)]
                    self.assertEqual(names, ["real.py"])


if __name__ == "__main__":
    unittest.main()
