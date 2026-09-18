"""Documentation checks -- the defects here had user-visible consequences.

The worst was ``docs/tutorial.rst`` telling the reader to run ``mainprogram 20``
for the partial density of states.  Process 20 is ``clean-scan``: it deletes the
wavefunctions and moves the run to ``completed/``.  Following the tutorial
destroyed the calculation.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs"


def rst_files():
    return sorted(DOCS.glob("*.rst"))


class GeneratedPages(unittest.TestCase):
    def test_command_rst_is_up_to_date(self):
        result = subprocess.run(
            [sys.executable, str(ROOT / "tools" / "gen_command_rst.py"), "--check"],
            capture_output=True, text=True, cwd=ROOT)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


class WrongNames(unittest.TestCase):
    """Four commands and two file names were documented under names that do
    not exist in the code."""

    FORBIDDEN = {
        "elastic-compute": "the command is compute-elastic",
        "epw6-file": "the command is wann-file",
        "epw8-file": "the command is epw-file",
        "distort-extract.py": "the file is distort_extract.py",
        "generate_submission.sh": "the file is generate_submission_file.sh",
        "econv_vasp.csv": "the file written is econv.csv",
        "wanniertool_input": "the config key is wanniertools_input",
        "mainprogram inputinfo": "the command is mainprogram basicinfo",
    }

    def test_no_document_uses_a_name_that_does_not_exist(self):
        offenders = []
        for path in rst_files():
            text = path.read_text()
            for wrong, correction in self.FORBIDDEN.items():
                if wrong in text:
                    line = next((i + 1 for i, ln in enumerate(text.splitlines())
                                 if wrong in ln), 0)
                    offenders.append(
                        f"{path.name}:{line}: {wrong!r} -- {correction}")
        self.assertEqual(offenders, [])


class DestructiveInstructions(unittest.TestCase):
    def test_pdos_is_not_documented_as_process_20(self):
        """Process 20 is clean-scan; it deletes wavefunctions."""
        pattern = re.compile(r"(?i)(partial\s+dos|pdos)[^\n]{0,120}mainprogram\s+20")
        offenders = []
        for path in rst_files():
            text = path.read_text()
            for match in pattern.finditer(text):
                line = text[:match.start()].count("\n") + 1
                offenders.append(f"{path.name}:{line}")
        self.assertEqual(offenders, [])

    def test_process_20_is_described_as_destructive_where_it_appears(self):
        text = "\n".join(p.read_text() for p in rst_files())
        if "mainprogram 20" not in text and "process = 20" not in text:
            self.skipTest("process 20 is not mentioned")
        lowered = text.lower()
        self.assertTrue(any(word in lowered for word in
                            ("delete", "removes", "removing", "clean")),
                        "process 20 must be described as destructive")


class CodeBlocks(unittest.TestCase):
    def test_every_json_code_block_parses(self):
        """The flagship config block in param.rst had a stray ] and two
        unclosed braces, so anyone copy-pasting it got a parse error."""
        offenders = []
        for path in rst_files():
            lines = path.read_text().splitlines()
            index = 0
            while index < len(lines):
                if lines[index].strip().startswith(".. code-block:: json"):
                    index += 1
                    while index < len(lines) and not lines[index].strip():
                        index += 1
                    body, indent = [], None
                    while index < len(lines):
                        line = lines[index]
                        if not line.strip():
                            body.append("")
                            index += 1
                            continue
                        current = len(line) - len(line.lstrip())
                        if indent is None:
                            indent = current
                        if current < indent:
                            break
                        body.append(line[indent:])
                        index += 1
                    snippet = "\n".join(body).strip()
                    if not snippet or not snippet.startswith("{"):
                        continue          # a documented fragment, not a document
                    try:
                        json.loads(snippet)
                    except ValueError as exc:
                        offenders.append(f"{path.name}: {exc}")
                else:
                    index += 1
        self.assertEqual(offenders, [])


class References(unittest.TestCase):
    def test_no_duplicate_labels(self):
        """`_pressure-label` was defined twice, so its two :ref:s were random."""
        seen: dict[str, str] = {}
        duplicates = []
        for path in rst_files():
            for number, line in enumerate(path.read_text().splitlines(), 1):
                match = re.match(r"\.\.\s+_([A-Za-z0-9_-]+):\s*$", line.strip())
                if not match:
                    continue
                label = match.group(1)
                where = f"{path.name}:{number}"
                if label in seen:
                    duplicates.append(f"{label} ({seen[label]} and {where})")
                seen[label] = where
        self.assertEqual(duplicates, [])

    def test_every_ref_resolves(self):
        labels = set()
        for path in rst_files():
            labels.update(re.findall(r"\.\.\s+_([A-Za-z0-9_-]+):", path.read_text()))
        broken = []
        for path in rst_files():
            for target in re.findall(r":ref:`[^<`]*<([A-Za-z0-9_-]+)>`",
                                     path.read_text()):
                if target not in labels:
                    broken.append(f"{path.name}: {target}")
        self.assertEqual(broken, [])


class CommandLineReferenceIsComplete(unittest.TestCase):
    """The reference must list every installed command and every option.

    It was written by reading the parsers, and parsers change. These assertions
    fail when a new entry point or option is added without documenting it.
    """

    def _usage(self):
        return (DOCS / "usage.rst").read_text()

    def test_every_console_script_is_documented(self):
        import re

        pyproject = (ROOT / "pyproject.toml").read_text()
        block = pyproject.split("[project.scripts]", 1)[1].split("\n[", 1)[0]
        names = re.findall(r"^([A-Za-z0-9_-]+)\s*=", block, re.MULTILINE)
        self.assertTrue(names)
        text = self._usage()
        for name in names:
            if name == "htesp":          # documented as an alias of mainprogram
                continue
            with self.subTest(command=name):
                self.assertIn(name, text)

    def test_every_htesp_check_option_is_documented(self):
        import argparse
        import io
        import contextlib
        import sys

        sys.path.insert(0, str(ROOT))
        from htesp.check import main as check_main

        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            try:
                check_main(["--help"])
            except SystemExit:
                pass
        import re

        flags = set(re.findall(r"(--[a-z0-9_-]+)", buf.getvalue()))
        text = self._usage()
        for flag in sorted(flags):
            with self.subTest(flag=flag):
                self.assertIn(flag, text)

    def test_every_mainprogram_global_option_is_documented(self):
        """`mainprogram`'s own flags drift the same way `htesp-check`'s do."""
        import re
        import sys

        sys.path.insert(0, str(ROOT))
        from htesp.cli import build_parser

        flags = set()
        for action in build_parser()._actions:
            flags.update(f for f in action.option_strings if f.startswith("--"))
        self.assertIn("--init-header", flags)
        text = self._usage()
        for flag in sorted(flags):
            with self.subTest(flag=flag):
                self.assertIn(flag, text)

    def test_the_environment_variables_are_listed(self):
        text = self._usage()
        for var in ("MP_API_KEY", "HTESP_CONFIG", "HTESP_WORKERS",
                    "HTESP_EXAMPLES", "HTESP_PYTHON", "PMG_VASP_PSP_DIR"):
            with self.subTest(var=var):
                self.assertIn(var, text)

    def test_the_readme_covers_the_post_install_steps(self):
        text = (ROOT / "README.md").read_text()
        for step in ("--set_mp_api", "--config_vasp_pot", "--install-enumlib",
                     "--clean", "config-init", "htesp-tutorials",
                     "--init-header"):
            with self.subTest(step=step):
                self.assertIn(step, text)

    def test_the_readme_does_not_call_required_packages_optional(self):
        """lmfit and bsym are in [project] dependencies, not extras."""
        text = (ROOT / "README.md").read_text()
        self.assertNotIn("Extra packages", text)


if __name__ == "__main__":
    unittest.main()
