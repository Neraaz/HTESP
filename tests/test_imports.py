"""Import smoke test.

`create_epw_inputs`, `create_wt_inputs` and `element_extract` used to fail at
*import* time -- `create_wt_inputs` read a mis-spelled config key at module
scope, so `wt1`/`wt2` could never run.  Any module that cannot be imported on a
complete installation is a bug, so this test imports every one of them and
skips only the ones whose third-party dependency is genuinely absent here.

The sweep runs in a **child interpreter**.  A third-party extension module that
fails to load does not always raise: ``mp_api`` -> ``deltalake`` aborts the
process outright (``Fatal Python error: Aborted``, SIGABRT) on a machine whose
compiled ``deltalake`` does not match its ``pyarrow``.  Imported in-process
that killed the whole pytest run mid-file, with no result for the 100-odd tests
that had not run yet.  The child prints its progress line by line, so when it
dies the parent still knows which module was being imported and turns it into
one ordinary test failure.
"""
from __future__ import annotations

import importlib
import os
import signal
import subprocess
import sys
import unittest
from pathlib import Path

from tests.helpers import TempProject, have, python_sources

REPO_ROOT = Path(__file__).resolve().parent.parent
PKG = REPO_ROOT / "htesp"

#: module -> the third-party packages it needs.  A module missing from this
#: table must import with nothing but the standard library and numpy.
OPTIONAL_DEPENDENCIES = {
    # aflow_extract imports oqmd_extract, which imports qmpy_rester
    "aflow_extract": ("pymatgen", "qmpy_rester"),
    "cif_to_gsinput": ("scipy", "pymatgen"),
    "convergence_test": ("pymatgen",), "create_epw_inputs": ("lmfit", "ase"),
    "create_wt_inputs": ("ase",), "crystal": ("pymatgen", "spglib"),
    "displace_phonopy": ("pymatgen", "yaml"), "elastic": ("pymatgen",),
    "element_extract": ("scipy", "pymatgen", "ase"), "elph": ("pymatgen",),
    "htepc": ("scipy", "pymatgen", "ase"), "kpath": ("ase",),
    "kpoint_path": ("ase",), "magnetic": ("pymatgen",),
    "ml_processing": ("scipy", "pymatgen", "ase"),
    "oqmd_extract": ("ase", "pymatgen", "qmpy_rester"),
    "plot": ("pymatgen",), "plot_bandproj": ("pymatgen",),
    "poscar_to_vasp": ("pymatgen",), "projection_phband": ("ase",),
    "pymatgen_phase_diagram": ("pymatgen",), "qe_input": ("scipy", "pymatgen"),
    "scftocif": ("ase",), "site_subs": ("bsym", "pymatgen"),
    "structure_group": ("pymatgen",), "vasp_input": ("pymatgen",),
    "vasp_process": ("ase", "pymatgen"), "wannier90": ("ase",),
    "write_potcar": ("pymatgen",),
}

#: the modules the CLI and the workflow layer need: these must import with the
#: standard library alone, so `mainprogram --help` works on a bare machine
CORE_MODULES = ("banner", "cli", "config", "check_json", "help_text",
                "inputin", "mainprogram", "workflow")

#: ``mp_api`` pulls in ``deltalake``; see the module docstring.  Import it
#: lazily, in the function that contacts Materials Project, never at module
#: scope -- otherwise writing a QE input file loads a stack of Arrow/Rust
#: extensions it has no use for.
FORBIDDEN_AT_MODULE_SCOPE = ("mp_api", "matminer", "sklearn", "ifermi", "plotly")

#: run in a child interpreter: imports the named modules one at a time and
#: reports each outcome on its own flushed line, so a crash is attributable.
_SWEEP = r"""
import importlib, sys

for name in sys.argv[1:]:
    print("TRY " + name, flush=True)
    try:
        importlib.import_module("htesp." + name)
    except BaseException as exc:
        print("FAIL {0} {1}: {2}".format(name, type(exc).__name__, exc), flush=True)
    else:
        print("OK " + name, flush=True)
print("DONE", flush=True)
"""


class CoreImports(unittest.TestCase):
    def test_the_core_needs_nothing_but_the_standard_library(self):
        for name in CORE_MODULES:
            with self.subTest(module=name):
                importlib.import_module(f"htesp.{name}")


class AllImports(TempProject):
    def test_every_module_imports_or_is_honestly_missing_a_dependency(self):
        wanted, skipped = [], []
        for path in python_sources(PKG):
            name = path.stem
            if name.startswith("__"):
                continue
            missing = [dep for dep in OPTIONAL_DEPENDENCIES.get(name, ())
                       if not have(dep)]
            if missing:
                skipped.append(f"{name} (needs {', '.join(missing)})")
            else:
                wanted.append(name)

        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join(
            [str(REPO_ROOT)] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
        # matplotlib must not try to open a display from the child
        env.setdefault("MPLBACKEND", "Agg")

        failures = []
        remaining = list(wanted)
        while remaining:
            child = subprocess.run(
                [sys.executable, "-c", _SWEEP, *remaining],
                capture_output=True, text=True, env=env, timeout=900)
            lines = child.stdout.splitlines()
            failures += [line[len("FAIL "):] for line in lines
                         if line.startswith("FAIL ")]
            if "DONE" in lines:
                break
            # the child died: the last TRY without a verdict names the culprit.
            # Record it and restart on what is left, so one module that takes
            # the interpreter down with it does not hide the others.
            settled = {line.split()[1] for line in lines
                       if line.startswith(("OK ", "FAIL "))}
            attempted = [line[len("TRY "):] for line in lines
                         if line.startswith("TRY ")]
            culprit = next((m for m in reversed(attempted) if m not in settled),
                           None)
            failures.append(
                f"{culprit or remaining[0]} killed the interpreter "
                f"({_exit_reason(child.returncode)}). "
                f"stderr: {child.stderr.strip()[-2000:]}")
            cut = remaining.index(culprit) if culprit in remaining else 0
            remaining = remaining[cut + 1:]

        if skipped:
            print(f"\n  skipped {len(skipped)} module(s) with absent dependencies:")
            for line in skipped:
                print(f"    {line}")

        self.assertEqual(failures, [])

    def test_no_module_reads_a_file_at_import_time(self):
        """A module-scope `config()` call makes the import depend on the cwd."""
        import ast
        offenders = []
        for path in python_sources(PKG):
            tree = ast.parse(path.read_text())
            for node in tree.body:                   # module scope only
                for inner in ast.walk(node):
                    if not isinstance(inner, ast.Call):
                        continue
                    func = inner.func
                    called = getattr(func, "id", None) or getattr(func, "attr", None)
                    if called in ("config", "load", "loads") and isinstance(node, ast.Assign):
                        offenders.append(f"{path.name}:{inner.lineno} {called}()")
        self.assertEqual(offenders, [])

    def test_no_module_imports_a_heavy_optional_package_at_module_scope(self):
        """`import mp_api` at module scope loaded deltalake for input writing."""
        import ast
        offenders = []
        for path in python_sources(PKG):
            tree = ast.parse(path.read_text())
            for node in tree.body:                   # module scope only
                roots = []
                if isinstance(node, ast.ImportFrom) and node.level == 0:
                    roots.append((node.module or "").split(".")[0])
                elif isinstance(node, ast.Import):
                    roots += [alias.name.split(".")[0] for alias in node.names]
                for root in roots:
                    if root in FORBIDDEN_AT_MODULE_SCOPE:
                        offenders.append(f"{path.name}:{node.lineno} {root}")
        self.assertEqual(offenders, [])


def _exit_reason(returncode: int) -> str:
    if returncode < 0:
        try:
            return f"killed by {signal.Signals(-returncode).name}"
        except ValueError:                           # pragma: no cover
            return f"killed by signal {-returncode}"
    return f"exit status {returncode}"


if __name__ == "__main__":
    unittest.main()
