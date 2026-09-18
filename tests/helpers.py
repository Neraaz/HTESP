"""Shared fixtures and guards for the HTESP test suite."""
from __future__ import annotations

import importlib.util
import os
import shutil
import sys
import tempfile
import unittest
from contextlib import contextmanager
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


#: names that look like Python sources but are not: macOS AppleDouble
#: sidecars (``._module.py``), editor swap files and dotfiles.  Copying the
#: tree to a non-HFS filesystem -- scp/rsync/tar to an HPC login node --
#: creates one ``._x.py`` per ``x.py``, each a small binary blob.  Globbing
#: ``*.py`` then hands ast.parse and py_compile a file full of null bytes.
def is_python_source(path) -> bool:
    """True for a real module file, False for AppleDouble and other sidecars."""
    name = Path(path).name
    return name.endswith(".py") and not name.startswith((".", "._"))


def python_sources(directory, recursive: bool = False) -> list:
    """Every real ``*.py`` in ``directory``, sorted, sidecars excluded.

    ``recursive=True`` walks sub-directories as well, which is what
    ``tools/check_names.py`` needs.  It lives here rather than in the tool so
    that "what counts as a Python source" has exactly one definition: the tool
    used to glob ``*.py`` itself, and so died on the first AppleDouble sidecar
    it met while the test suite sailed past them.
    """
    walk = Path(directory).rglob if recursive else Path(directory).glob
    return sorted(p for p in walk("*.py") if is_python_source(p))


def stray_sidecars(root) -> list:
    """AppleDouble/resource-fork files anywhere under ``root`` (excluding .git)."""
    return sorted(p for p in Path(root).rglob("._*")
                  if ".git" not in p.parts and "_removed" not in p.parts)


def have(module: str) -> bool:
    """True when ``module`` can be imported without importing it."""
    try:
        return importlib.util.find_spec(module) is not None
    except (ImportError, ValueError):
        return False


def skip_without(*modules: str):
    """Skip a test that needs optional third-party packages."""
    missing = [m for m in modules if not have(m)]
    return unittest.skipIf(
        bool(missing), f"needs {', '.join(missing)}")


class TempProject(unittest.TestCase):
    """A throw-away project directory, made current for the duration of a test.

    Every HTESP command is relative to "the project directory" (the one holding
    ``input.in``), so almost every test wants one of these.
    """

    def setUp(self) -> None:
        super().setUp()
        self._previous = Path.cwd()
        self._tmp = tempfile.mkdtemp(prefix="htesp-test-")
        self.root = Path(self._tmp)
        os.chdir(self.root)
        # a stale $HTESP_CONFIG from the developer's shell would defeat the
        # config-search tests
        self._saved_env = {k: os.environ.get(k)
                           for k in ("HTESP_CONFIG", "MP_API_KEY", "HTESP_WORKERS")}
        for key in self._saved_env:
            os.environ.pop(key, None)
        from htesp import config as config_module
        # FIX: api_key() also reads ~/.config/htesp/credentials, which is
        # outside the temporary project.  Until this was isolated the suite
        # passed or failed depending on whether the developer happened to have
        # a real credentials file -- writing one with
        # `htesp-check --set_mp_api` broke four tests that had nothing to do
        # with it.  Point it at the throw-away directory instead.
        self._saved_credentials = config_module.CREDENTIALS_PATH
        config_module.CREDENTIALS_PATH = self.root / "no-credentials-here"
        config_module.clear_cache()

    def tearDown(self) -> None:
        os.chdir(self._previous)
        shutil.rmtree(self._tmp, ignore_errors=True)
        for key, value in self._saved_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
        from htesp import config as config_module
        config_module.CREDENTIALS_PATH = self._saved_credentials
        config_module.clear_cache()
        super().tearDown()

    # -- convenience -------------------------------------------------------
    def write(self, relative: str, text: str) -> Path:
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def write_input_in(self, start: int = 1, end: int = 3, nkpt: int = 200,
                       track: str = "mpid.in", plots: str = "phband",
                       dft: str = "QE") -> Path:
        return self.write(
            "input.in", f"{start}\n{end}\n{nkpt} 0\n{track}\n{plots}\nDFT = {dft}\n")

    def write_track(self, name: str = "mpid.in", entries=(("mp-763", "Mg1B2"),)) -> Path:
        lines = [f"v{i + 1} {mpid} {comp}" for i, (mpid, comp) in enumerate(entries)]
        return self.write(name, "\n".join(lines) + "\n")


@contextmanager
def pushd(path):
    """``cd`` for the duration of a block (tests only; the package has its own)."""
    previous = Path.cwd()
    os.chdir(path)
    try:
        yield Path(path)
    finally:
        os.chdir(previous)


#: a minimal but realistic `vc-relax` output, used by several parser tests
VC_RELAX_OUT = """
     Program PWSCF v.7.2 starts
     number of electrons       =        26.00
     number of Kohn-Sham states=           18
     number of atoms/cell      =            3

     the Fermi energy is     8.1234 ev

     total   stress  (Ry/bohr**3)                   (kbar)     P=      -12.34
   0.00001000   0.00000000   0.00000000           1.47        0.00        0.00
   0.00000000   0.00001000   0.00000000           0.00        1.47        0.00
   0.00000000   0.00000000   0.00002000           0.00        0.00        2.94

     the Fermi energy is     8.4321 ev

     A final scf calculation at the relaxed structure.
     Final scf calculation at the relaxed structure.

Begin final coordinates
     new unit-cell volume =     62.12345 a.u.^3 (    9.2050 Ang^3 )
     density =      2.6250 g/cm^3

CELL_PARAMETERS (angstrom)
   3.0850000000   0.0000000000   0.0000000000
  -1.5425000000   2.6716000000   0.0000000000
   0.0000000000   0.0000000000   3.5230000000

ATOMIC_POSITIONS (crystal)
Mg       0.0000000000        0.0000000000        0.0000000000
B        0.3333333333        0.6666666667        0.5000000000
B        0.6666666667        0.3333333333        0.5000000000
End final coordinates

     JOB DONE.
"""

#: a `relax` output: no cell block, so the fixed `sed '1,4d' | sed '5d'`
#: pipeline of the bash layer deleted the first atoms
RELAX_OUT = """
     number of electrons       =        26.00
     number of atoms/cell      =            3

     the Fermi energy is     7.7777 ev

Begin final coordinates
ATOMIC_POSITIONS (crystal)
Mg       0.0000000000        0.0000000000        0.0000000000
B        0.3333333333        0.6666666667        0.5000000000
B        0.6666666667        0.3333333333        0.5000000000
End final coordinates

     JOB DONE.
"""

#: an output that stopped on the wall clock: no final-coordinates block at all
UNCONVERGED_OUT = """
     number of electrons       =        26.00
     number of atoms/cell      =            3

ATOMIC_POSITIONS (crystal)
Mg       0.0100000000        0.0000000000        0.0000000000
B        0.3333333333        0.6666666667        0.5000000000
B        0.6666666667        0.3333333333        0.5000000000

     Maximum CPU time exceeded
"""

SCF_IN = """&CONTROL
  calculation = 'vc-relax',
  prefix = 'Mg1B2',
  pseudo_dir = '../../pp',
  outdir = './',
/
&SYSTEM
  ibrav = 0,
  nat = 3,
  ntyp = 2,
  ecutwfc = 45.0,
  ecutrho = 360.0,
  occupations = 'smearing',
  smearing = 'mv',
  degauss = 0.02,
/
&ELECTRONS
  conv_thr = 1e-08,
/
&IONS
/
&CELL
/
ATOMIC_SPECIES
  Mg 24.305 Mg.upf
  B  10.811 B.upf
ATOMIC_POSITIONS crystal
  Mg 0.000000 0.000000 0.000000
  B  0.333333 0.666667 0.500000
  B  0.666667 0.333333 0.500000
K_POINTS automatic
  12 12 8 0 0 0
CELL_PARAMETERS angstrom
  3.085000 0.000000 0.000000
 -1.542500 2.671600 0.000000
  0.000000 0.000000 3.523000
"""
