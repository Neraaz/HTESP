#!/usr/bin/env python
"""Report whether this machine can actually run HTESP.

``htesp-check``   (or ``python -m htesp.check``)

Why this exists
---------------
Nothing in HTESP itself depends on the processor architecture.  There is no
assembly, no ``ctypes``, no ``struct`` packing, no endianness assumption, no
``platform.machine()`` branch, and every shell shim is POSIX ``sh``.  The
package is the same on x86_64 and on arm64/aarch64.

What *does* differ is every compiled wheel underneath it: numpy, scipy,
pymatgen's Cython extensions, spglib, pyarrow, and -- the one that broke a real
run on an aarch64 login node -- ``deltalake``, which ``mp_api`` imports.  A
wheel that is missing, built for the wrong architecture, or built against a
mismatched sibling does not always raise ``ImportError``.  It can abort the
interpreter outright with SIGABRT, and no ``try``/``except`` catches that: the
process is simply gone, which is why it took a whole ``pytest`` run with it.

So every probe here runs in a **child interpreter**.  A dependency that aborts
is reported as one line of output rather than taking the report down with it.

Exit status is ``0`` when every required dependency imports, ``1`` otherwise.
"""
from __future__ import annotations

import argparse
import functools
import importlib.util
import json
import os
import platform
import re
import shutil
import signal
import subprocess
import sys
from importlib import metadata
from pathlib import Path

#: the dependencies declared in pyproject.toml ``[project] dependencies``,
#: as import names (``PyYAML`` imports as ``yaml``, ``mp_api`` as ``mp_api``).
REQUIRED = ("numpy", "scipy", "pandas", "matplotlib", "pymatgen", "mp_api",
            "emmet.core", "ase", "spglib", "yaml", "bsym", "lmfit",
            "qmpy_rester")

#: optional extras -> the import names they provide
EXTRAS = {
    "ml": ("matminer", "sklearn"),
    "fermisurface": ("ifermi", "plotly"),
    "test": ("pytest",),
    "docs": ("sphinx",),
}

#: transitive packages that are worth naming because they carry the compiled
#: pieces that differ between architectures
TRANSITIVE = ("pyarrow", "deltalake")

#: import name -> the name to type after ``pip install``, where they differ.
#: "No module named 'qmpy_rester'" is true but unhelpful: the distribution is
#: ``qmpy-rester``, and nothing in the error says so.
DISTRIBUTION = {
    "yaml": "PyYAML",
    "sklearn": "scikit-learn",
    "qmpy_rester": "qmpy-rester",
    "mp_api": "mp-api",
    "emmet.core": "emmet-core",
}

#: import name -> the extra that provides it, so the hint can be the whole
#: extra (``pip install "htesp[fermisurface]"``) rather than a loose
#: package.  Only the extras in EXTRAS appear here: ``qmpy_rester`` is a
#: required dependency, not an extra, so its hint is ``pip install
#: qmpy-rester`` -- see test_install_hint_names_the_distribution.
PROVIDED_BY = {}

#: the module to import when it is not the top-level package.  ``import
#: mp_api`` only executes ``mp_api/__init__.py`` and reports "ok" on a machine
#: where ``mp_api.client`` -- the thing HTESP imports, and the thing that pulls
#: in ``deltalake`` -- aborts.  Probe what is actually used.
IMPORT_TARGET = {"mp_api": "mp_api.client"}

for _extra, _modules in EXTRAS.items():
    for _module in _modules:
        PROVIDED_BY[_module] = _extra

#: external programs HTESP shells out to.  Absent is not an error -- a login
#: node that only prepares inputs needs none of them.
#: the enumlib executables pymatgen's EnumlibAdaptor shells out to, and the
#: two that ``--install-enumlib`` builds.  One source of truth: they are part
#: of :data:`EXECUTABLES` below, because ``--executables`` not listing them was
#: how a successful ``--install-enumlib`` still looked like it had done
#: nothing -- the report simply never probed for them.
ENUMLIB_BINARIES = ("enum.x", "makestr.x")

EXECUTABLES = ("pw.x", "ph.x", "q2r.x", "matdyn.x", "epw.x", "lambda.x",
               "plotband.x", "sumpdos.x", "wannier90.x", "wt.x",
               "vasp_std", "phonopy", "sbatch", "squeue") + ENUMLIB_BINARIES

#: run in a child interpreter: import one module, print one line, exit.
_PROBE = (
    "import importlib, sys\n"
    "name = sys.argv[1]\n"
    "try:\n"
    "    module = importlib.import_module(name)\n"
    "except ModuleNotFoundError as exc:\n"
    "    absent = getattr(exc, 'name', None)\n"
    "    own = absent and (name == absent or name.startswith(absent + '.'))\n"
    "    print(('MISSING ' if own else 'BROKEN ModuleNotFoundError: ') + str(exc))\n"
    "except ImportError as exc:\n"
    "    print('BROKEN ImportError: ' + str(exc))\n"
    "except BaseException as exc:\n"
    "    print('ERROR {0}: {1}'.format(type(exc).__name__, exc))\n"
    "else:\n"
    "    version = getattr(module, '__version__', '')\n"
    "    if not version:\n"
    "        import importlib.metadata as meta\n"
    "        try:\n"
    "            version = meta.version(name.split('.')[0])\n"
    "        except Exception:\n"
    "            version = ''\n"
    "    print('OK ' + str(version))\n"
)


def _exit_reason(returncode: int) -> str:
    """``-6`` -> ``killed by SIGABRT``; anything else -> the status."""
    if returncode < 0:
        try:
            return f"killed by {signal.Signals(-returncode).name}"
        except ValueError:                       # pragma: no cover
            return f"killed by signal {-returncode}"
    return f"exit status {returncode}"


#: substrings in a dying child's output -> (cause, what to do about it).
#: An aborting extension is not always a wrong-architecture wheel; on this
#: project the first real case was a *right*-architecture wheel built for the
#: wrong memory page size.
_CAUSES = (
    (("unsupported system page size", "memory allocation of"),
     "page-size mismatch",
     "The wheel is for this architecture but its bundled allocator (jemalloc)\n"
     "was built assuming 4 KiB pages, and this kernel uses {pagesize}.\n"
     "Reinstalling the same wheel will not help; build it here instead:\n"
     "  JEMALLOC_SYS_WITH_LG_PAGE={lg_page} pip install --force-reinstall \\\n"
     "      --no-cache-dir --no-binary {names_csv} {names}"),
    (("incompatible architecture", "wrong elf class", "cannot open shared object",
      "mach-o", "no matching architecture"),
     "wrong architecture",
     "The extension is built for a different architecture than {machine}.\n"
     "  pip install --force-reinstall --no-cache-dir {names}"),
    (("illegal instruction",),
     "wrong CPU baseline",
     "The extension uses instructions this CPU does not have.  Build it here:\n"
     "  pip install --force-reinstall --no-binary {names_csv} {names}"),
)


def distribution(module: str) -> str:
    """The name to install for an import name (``yaml`` -> ``PyYAML``)."""
    return DISTRIBUTION.get(module, module)


def install_hint(module: str) -> str:
    """What to type to get ``module``: the extra if it belongs to one.

    The distribution is named too whenever it differs from the import name --
    ``qmpy_rester`` comes from ``qmpy-rester``, ``sklearn`` from
    ``scikit-learn`` -- because the ImportError never says so.
    """
    dist = distribution(module)
    extra = PROVIDED_BY.get(module)
    if extra:
        suffix = f"   # {dist}" if dist != module else ""
        return f'pip install "htesp[{extra}]"{suffix}'
    return f"pip install {dist}"


def page_size() -> int:
    """The kernel's memory page size in bytes (4096, 16384 or 65536)."""
    try:
        return os.sysconf("SC_PAGESIZE")
    except (ValueError, AttributeError, OSError):   # pragma: no cover
        return 0


def classify(text: str) -> str:
    """Name the cause of a dying child from what it printed, or ``""``."""
    lowered = (text or "").lower()
    for needles, cause, _ in _CAUSES:
        if any(needle in lowered for needle in needles):
            return cause
    return ""


def probe(name: str, timeout: float = 120.0) -> dict:
    """Import ``name`` in a child interpreter; never raise, never abort.

    Returns ``{"module", "status", "detail"}`` where ``status`` is one of
    ``ok``, ``missing``, ``error`` or ``aborted``; an aborted probe also
    carries ``cause``.
    """
    target = IMPORT_TARGET.get(name, name)
    env = dict(os.environ)
    env.setdefault("MPLBACKEND", "Agg")
    try:
        child = subprocess.run([sys.executable, "-c", _PROBE, target],
                               capture_output=True, text=True,
                               env=env, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"module": name, "status": "error",
                "detail": f"import did not finish in {timeout:.0f}s"}

    line = (child.stdout or "").strip().splitlines()
    head = line[-1] if line else ""
    if head.startswith("OK"):
        return {"module": name, "status": "ok",
                "detail": head[len("OK"):].strip(),
                "distribution": distribution(name)}
    if head.startswith("MISSING "):
        return {"module": name, "status": "missing",
                "detail": head[len("MISSING "):],
                "distribution": distribution(name),
                "install": install_hint(name)}
    if head.startswith("BROKEN "):
        # installed, but its own import failed -- nearly always a version
        # mismatch with one of its dependencies, which `pip install` of this
        # package alone will not fix
        return {"module": name, "status": "broken",
                "detail": head[len("BROKEN "):],
                "distribution": distribution(name)}
    if head.startswith("ERROR "):
        return {"module": name, "status": "error",
                "detail": head[len("ERROR "):]}
    # no verdict: the child died before it could print one
    detail = _exit_reason(child.returncode)
    stderr = (child.stderr or "").strip()
    if stderr:
        detail = f"{detail}; {stderr.splitlines()[-1]}"
    return {"module": name, "status": "aborted", "detail": detail,
            "cause": classify(f"{child.stdout}\n{stderr}")}


@functools.lru_cache(maxsize=1)
def required_versions() -> dict:
    """distribution -> version specifier, from the installed htesp metadata.

    ``importlib.metadata`` is used rather than ``pyproject.toml`` because the
    latter is not shipped with the installed package, and what matters here is
    what *this* installation asks for.
    """
    try:
        requirements = metadata.requires("htesp") or []
    except metadata.PackageNotFoundError:
        return {}
    specifiers = {}
    for line in requirements:
        item = line.split(";")[0].strip()          # drop the extra marker
        match = re.match(r"^([A-Za-z0-9_.-]+)\s*(.*)$", item)
        if match:
            specifiers[match.group(1).lower().replace("_", "-")] = \
                match.group(2).strip()
    return specifiers


def wanted_version(module: str) -> str:
    """The specifier this installation declares for ``module``, or ""."""
    dist = distribution(module).lower().replace("_", "-")
    return required_versions().get(dist, "")


def installation() -> dict:
    """Which HTESP is running, from where, and any stale install beside it.

    "I reinstalled and nothing changed" is almost always one of two things: a
    second copy earlier on ``sys.path``, or the HTESP 1.x distribution still
    present.  1.x was named ``HTESP`` (not ``htesp``) and installed a top-level
    package literally called ``src`` plus 46 modules as bare commands on
    ``$PATH``, so the two can sit side by side without pip noticing.
    """
    found = {}
    for name in ("htesp", "HTESP"):
        try:
            found[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            continue
    warnings = []
    if len(found) > 1:
        warnings.append("both `htesp` and `HTESP` are installed -- the second "
                        "is HTESP 1.x; `pip uninstall -y HTESP` removes it")
    try:
        if importlib.util.find_spec("src") is not None:
            warnings.append("a top-level `src` package is importable; that is "
                            "HTESP 1.x's packaging bug -- uninstall it")
    except (ImportError, ValueError):                # pragma: no cover
        pass
    return {"package": str(Path(__file__).resolve().parent),
            "distributions": found, "warnings": warnings}


def environment() -> dict:
    """The facts that decide which wheels this interpreter needs."""
    return {
        "python": platform.python_version(),
        "implementation": platform.python_implementation(),
        "executable": sys.executable,
        "system": platform.system(),
        "release": platform.release(),
        "machine": platform.machine(),
        "byteorder": sys.byteorder,
        "pagesize": page_size(),
        "cpus": os.cpu_count(),
    }


def report(extras: bool = True) -> dict:
    """Probe everything and return the whole result as plain data."""
    result = {"environment": environment(),
              "installation": installation(),
              "required": [probe(name) for name in REQUIRED],
              "transitive": [probe(name) for name in TRANSITIVE],
              "extras": {}, "executables": {}}
    if extras:
        for extra, modules in EXTRAS.items():
            result["extras"][extra] = [probe(name) for name in modules]
    for name in EXECUTABLES:
        result["executables"][name] = shutil.which(name)
    return result


_SYMBOL = {"ok": "ok      ", "missing": "MISSING ", "broken": "BROKEN  ",
           "error": "ERROR   ", "aborted": "ABORTED "}


def _print_group(title: str, probes: list) -> None:
    print(f"\n{title}")
    for entry in probes:
        mark = _SYMBOL.get(entry["status"], "?       ")
        detail = entry["detail"]
        if entry["status"] == "missing":
            # the version column stays a version column: say what this
            # installation asks for, not how to get it -- the commands are
            # gathered under "To install what is missing" at the end
            specifier = wanted_version(entry["module"])
            detail = f"not installed, needs {specifier}" if specifier \
                else "not installed"
        print(f"  {mark} {entry['module']:<12} {detail}"[:110].rstrip())


def render(result: dict, show_executables: bool = False) -> None:
    """Print the human-readable form of :func:`report`."""
    env = result["environment"]
    print("HTESP environment report")
    pagesize = env.get("pagesize") or 0
    print(f"  {env['implementation']} {env['python']} on "
          f"{env['system']} {env['release']} / {env['machine']} "
          f"({env['byteorder']}-endian, {pagesize // 1024 or '?'} KiB pages, "
          f"{env['cpus']} CPUs)")
    print(f"  {env['executable']}")

    install = result.get("installation", {})
    versions = ", ".join(f"{name} {version}"
                         for name, version in install.get("distributions", {}).items())
    print(f"  HTESP {versions or '(not installed -- running from the source tree)'}")
    print(f"  running from {install.get('package', '?')}")
    for warning in install.get("warnings", []):
        print(f"  ! {warning}")

    _print_group("Required:", result["required"])
    _print_group("Compiled transitive dependencies "
                 "(informational -- they do not affect the exit status):",
                 result["transitive"])
    for extra, probes in result["extras"].items():
        _print_group(f"Extra htesp[{extra}]:", probes)

    if show_executables:
        print("\nExternal programs:")
        for name, path in result["executables"].items():
            print(f"  {'ok      ' if path else 'absent  '} {name:<12} "
                  f"{path or ''}")

    broken = [entry for entry in result["required"]
              if entry["status"] != "ok"]
    aborted = [entry for entry in result["required"] + result["transitive"]
               if entry["status"] == "aborted"]
    # if nothing required is broken, a transitive abort is no longer in the
    # way: mp-api 0.45.0 does not import deltalake at all, so an aborting
    # deltalake beside a healthy required set is a leftover, not a fault
    if aborted and not broken:
        print("Nothing required is broken.  "
              + ", ".join(entry["module"] for entry in aborted)
              + " still aborts, but no required package imports it here;")
        print("it can be left alone or uninstalled.")
        aborted = []
    print()
    if aborted:
        # a pure-Python package aborts only because of the compiled one it
        # imports: name the extension in the pip command, not its importer
        culprits = [entry for entry in aborted
                    if entry["module"] in TRANSITIVE] or aborted
        names = " ".join(entry["module"] for entry in culprits)
        names_csv = ",".join(entry["module"] for entry in culprits)
        causes = {entry.get("cause") for entry in aborted} - {""}
        pagesize = env.get("pagesize") or 4096
        print("A dependency ABORTED the interpreter instead of raising, so no")
        print("try/except could have caught it.")
        advice = [text for needles, cause, text in _CAUSES if cause in causes]
        if not advice:
            advice = ["Cause not recognised from the child's output.  Run it by\n"
                      "hand to see the message in full:\n"
                      "  python -c 'import {names}'"]
        for text in advice:
            print()
            print(text.format(names=names, names_csv=names_csv,
                              machine=env["machine"],
                              pagesize=f"{pagesize // 1024} KiB",
                              lg_page=max(pagesize.bit_length() - 1, 12)))
    if broken:
        print("Missing or broken required packages: "
              + ", ".join(entry["module"] for entry in broken))
    else:
        print("Every required dependency imports on this machine.")

    _print_install_block(result)


def _requirement(module: str) -> str:
    """``mp-api>=0.33``, quoted when a shell would mangle it."""
    dist = distribution(module)
    specifier = wanted_version(module)
    return f'"{dist}{specifier}"' if specifier else dist


def _print_install_block(result: dict) -> None:
    """One place for every command, after the report rather than inside it."""
    broken = [entry for entry in result["required"]
              if entry["status"] == "broken"]
    required = [entry["module"] for entry in result["required"]
                if entry["status"] == "missing"]
    extras = sorted(name for name, probes in result.get("extras", {}).items()
                    if any(entry["status"] == "missing" for entry in probes))
    if broken:
        print("\nInstalled but failing to import -- this is a version conflict,")
        print("not something `pip install <package>` will fix:")
        for entry in broken:
            print(f"  {entry['module']}: {entry['detail']}")
            print(f"    pip install \"{distribution(entry['module'])}"
                  f"{wanted_version(entry['module'])}\" --force-reinstall")
            print("    ...or pin the dependency it names to a version that "
                  "still provides it.")
    if not required and not extras:
        return
    print("\nTo install what is missing:")
    if required:
        print("  pip install " + " ".join(_requirement(name) for name in required))
    if extras:
        print(f'  pip install "htesp[{",".join(extras)}]"')



#: pymatgen's directory name for the PBE PAW set, and the other functionals it
#: recognises.  ``PMG_VASP_PSP_DIR`` must be the *parent* of one of these:
#: pymatgen looks for ``$PMG_VASP_PSP_DIR/<functional_dir>/<symbol>/POTCAR``.
VASP_FUNCTIONAL_DIRS = (
    "POT_GGA_PAW_PBE", "POT_GGA_PAW_PBE_52", "POT_GGA_PAW_PBE_54",
    "POT_GGA_PAW_PW91", "POT_LDA_PAW", "POT_LDA_PAW_52", "POT_LDA_PAW_54",
)


def configure_vasp_potcars(path: str) -> int:
    """Point pymatgen at a VASP POTCAR tree (``--config_vasp_pot``).

    VASP's POTCARs are licensed and cannot be shipped, so every command that
    writes VASP inputs fails with ``PmgVaspPspDirError: PMG_VASP_PSP_DIR is not
    set`` until pymatgen is told where they are.  Doing that by hand means
    knowing that the setting points at the *parent* of ``POT_GGA_PAW_PBE``
    rather than at it -- the mistake costs a confusing "POTCAR not found" for
    every material.  This accepts either and works out which was meant.

    Returns a process exit status: 0 when a POTCAR could actually be read
    afterwards, 1 otherwise.
    """
    given = Path(path).expanduser().resolve()
    if not given.is_dir():
        print("not a directory: {}".format(given))
        return 1

    # Accept the functional directory itself or its parent.
    if given.name in VASP_FUNCTIONAL_DIRS:
        root, functional = given.parent, given.name
    else:
        present = [name for name in VASP_FUNCTIONAL_DIRS if (given / name).is_dir()]
        if not present:
            print("{} contains no VASP potential directory.\n"
                  "Expected one of {} either at this path or inside it."
                  .format(given, ", ".join(VASP_FUNCTIONAL_DIRS)))
            return 1
        root, functional = given, present[0]

    # pymatgen accepts <functional>/<symbol>/POTCAR and <functional>/POTCAR.<symbol>
    family = root / functional
    symbols = [child.name for child in family.iterdir()
               if (child / "POTCAR").is_file()]
    flat = list(family.glob("POTCAR.*"))
    if not symbols and not flat:
        print("{} has neither <symbol>/POTCAR nor POTCAR.<symbol> entries.\n"
              "If this is a raw VASP tarball, unpack it first with:\n"
              "    pmg config -p {} <target>".format(family, family))
        return 1

    try:
        from pymatgen.cli.pmg_config import add_config_var
    except ImportError as exc:
        print("pymatgen is needed to write the configuration: {}".format(exc))
        return 1
    add_config_var(["PMG_VASP_PSP_DIR", str(root)], "")
    print("PMG_VASP_PSP_DIR = {}".format(root))
    print("  potential set   : {} ({} entries)".format(
        functional, len(symbols) or len(flat)))

    # Prove it: writing a POTCAR is the thing that was failing.
    try:
        import tempfile
        import warnings

        warnings.filterwarnings("ignore")
        from pymatgen.core import Lattice, Structure
        from pymatgen.io.vasp.sets import MPRelaxSet

        structure = Structure(Lattice.cubic(3.5), ["Mg"], [[0, 0, 0]])
        with tempfile.TemporaryDirectory() as tmp:
            MPRelaxSet(structure=structure).write_input(output_dir=tmp)
            first = (Path(tmp) / "POTCAR").read_text().splitlines()[0].strip()
        print("  verified        : wrote a POTCAR ({})".format(first))
        return 0
    except Exception as exc:                       # noqa: BLE001 - report anything
        print("  WARNING: the setting was written but a POTCAR still could not "
              "be produced: {}: {}".format(type(exc).__name__, exc))
        return 1



#: where enumlib is built from.  It is Fortran with a git submodule, has no
#: Python packaging of any kind, and is not on PyPI -- so it cannot be a
#: dependency in pyproject.toml and ``pip install -e .`` can never bring it in.
ENUMLIB_REPO = "https://github.com/msg-byu/enumlib.git"



def install_enumlib(prefix: str | None = None, compiler: str = "gfortran") -> int:
    """Build enumlib from source and install its executables (--install-enumlib).

    ``mainprogram magenum`` with ``magmom.type: "ordering"`` and an
    antiferromagnetic strategy goes through pymatgen's ``MagneticStructure
    Enumerator`` -> ``EnumerateStructureTransformation`` -> ``EnumlibAdaptor``,
    which shells out to ``enum.x`` and ``makestr.x``.  Those come from a
    Fortran codebase that pip cannot install and whose conda-forge build is
    linux-64/osx-64 only, so on aarch64 a source build is the only route.

    Builds in a temporary directory and installs into ``prefix`` (default: the
    ``bin`` of the running interpreter's environment, so it lands on PATH
    whenever that environment is active).

    Returns 0 when both executables are present and runnable afterwards.
    """
    import shutil as _shutil
    import subprocess as _subprocess
    import tempfile

    target = Path(prefix).expanduser().resolve() if prefix else Path(sys.prefix) / "bin"
    for tool in ("git", "make", compiler):
        if _shutil.which(tool) is None:
            print("{} is not on PATH; enumlib is Fortran and needs it to build."
                  .format(tool))
            if tool == compiler:
                print("  pass a different compiler with --fortran-compiler ifort")
            return 1
    try:
        target.mkdir(parents=True, exist_ok=True)
        probe = target / ".htesp-write-probe"
        probe.touch()
        probe.unlink()
    except OSError as exc:
        print("cannot write to {}: {}".format(target, exc))
        print("  choose somewhere else with --prefix ~/.local/bin")
        return 1

    print("building enumlib ({}), installing into {}".format(compiler, target))
    with tempfile.TemporaryDirectory(prefix="htesp-enumlib-") as tmp:
        root = Path(tmp) / "enumlib"
        steps = (
            (["git", "clone", "--recursive", "--depth", "1", ENUMLIB_REPO,
              str(root)], Path(tmp), "cloning"),
            (["make", "F90=" + compiler], root / "symlib" / "src",
             "building symlib"),
            (["make", "F90=" + compiler], root / "src", "building enumlib"),
            (["make", "F90=" + compiler, "enum.x"], root / "src", "linking enum.x"),
            (["make", "F90=" + compiler, "makestr.x"], root / "src",
             "linking makestr.x"),
        )
        for command, cwd, what in steps:
            print("  {}...".format(what))
            proc = _subprocess.run(command, cwd=os.fspath(cwd), capture_output=True,
                                   text=True, timeout=3600)
            if proc.returncode != 0:
                tail = (proc.stderr or proc.stdout).strip().splitlines()[-8:]
                print("  {} failed (exit {}):".format(what, proc.returncode))
                for line in tail:
                    print("    " + line)
                return 1

        installed = []
        for name in ENUMLIB_BINARIES:
            found = next((c for c in (root / "src" / name, root / name)
                          if c.is_file()), None)
            if found is None:
                hits = list(root.rglob(name))
                found = hits[0] if hits else None
            if found is None:
                print("  built, but {} was not produced".format(name))
                return 1
            _shutil.copy2(found, target / name)
            (target / name).chmod(0o755)
            installed.append(name)
        # makeStr.py is an accepted alternative to makestr.x and ships in aux_src
        aux = root / "aux_src" / "makeStr.py"
        if aux.is_file():
            _shutil.copy2(aux, target / "makeStr.py")
            (target / "makeStr.py").chmod(0o755)
            installed.append("makeStr.py")

    print("  installed: {}".format(", ".join(installed)))
    missing = [n for n in ENUMLIB_BINARIES if not (target / n).is_file()]
    if missing:
        print("  missing after install: {}".format(", ".join(missing)))
        return 1
    if _shutil.which("enum.x") is None:
        print("  NOTE: {} is not on your PATH; add it with".format(target))
        print("        export PATH=\"{}:$PATH\"".format(target))
    print("  verified: {}".format(target / "enum.x"))
    return 0



def _mask(key: str) -> str:
    """``'abcd...wxyz'`` -- never print a credential in full."""
    key = str(key)
    return key[:4] + "..." + key[-4:] if len(key) > 12 else "*" * len(key)


def set_mp_api_key(key: str, verify: bool = True) -> int:
    """Write the Materials Project key to the credentials file (--set_mp_api).

    ``htesp/config.py`` reads ``$MP_API_KEY`` first, then
    ``~/.config/htesp/credentials``, then ``config.json``.  The environment
    variable is the fragile one: it lives only in the shell that exported it,
    so a batch job, a ``nohup``-ed sweep or a new terminal silently loses it
    and every database command starts skipping or failing.  The credentials
    file survives all of those, and keeps the key out of ``config.json`` --
    which is what put a real key into 218 tracked files in the first place.

    Other lines in the file are preserved; only ``MP_API_KEY`` is replaced.
    The file is written ``0600``.  The key is never echoed in full.
    """
    from htesp.config import CREDENTIALS_PATH, API_KEY_ENV, API_KEY_PLACEHOLDER

    key = (key or "").strip()
    if not key:
        print("no key given")
        return 1
    if key == API_KEY_PLACEHOLDER:
        print("that is the placeholder shipped in the example configs, "
              "not a key.  Get one at https://next-gen.materialsproject.org/api")
        return 1

    # FIX: check the key before storing it.  The first version wrote the file
    # and *then* asked Materials Project, so a typo replaced a working key with
    # a broken one and the next command failed for a new reason.
    if verify:
        try:
            from htesp.htepc import mprester

            docs = mprester(key).materials.summary.search(
                material_ids=["mp-763"], fields=["material_id"])
            if not docs:
                print("the query returned nothing; the key may be wrong. "
                      "Nothing was written.")
                return 1
            verified = True
        except ImportError as exc:
            print("not verified (mp_api unavailable: {})".format(exc))
            verified = False
        except Exception as exc:                   # noqa: BLE001 - report anything
            print("the Materials Project rejected it: {}: {}"
                  .format(type(exc).__name__, str(exc)[:160]))
            print("Nothing was written; the existing credentials are untouched.")
            return 1
    else:
        verified = False

    path = Path(CREDENTIALS_PATH)
    kept = []
    if path.is_file():
        for line in path.read_text().splitlines():
            name = line.split("=", 1)[0].strip()
            if name not in (API_KEY_ENV, "api_key", "key"):
                kept.append(line)
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "\n".join(kept + ["{}={}".format(API_KEY_ENV, key)]) + "\n"
    path.write_text(body)
    try:
        path.chmod(0o600)
    except OSError:                                # pragma: no cover
        pass
    print("wrote {} ({})".format(path, _mask(key)))
    if kept:
        print("  kept {} other line(s) in the file".format(len(kept)))

    # The environment wins over the file, so say so rather than let a stale
    # exported value quietly keep winning.
    from_env = os.environ.get(API_KEY_ENV, "").strip()
    if from_env and from_env != key:
        print("  NOTE: ${} is set to a different key ({}) and takes "
              "precedence.  unset it, or that is the key that will be used."
              .format(API_KEY_ENV, _mask(from_env)))

    if verified:
        print("  verified: the Materials Project accepted it")
    return 0



#: what `pip install .` / `setup.py build` leave behind in a source checkout.
#: Directories are removed whole; ``*.pyc`` is swept separately.
BUILD_ARTIFACTS = ("build", "dist", "htesp.egg-info", "HTESP.egg-info", ".pytest_cache")

#: never swept.  ``.git`` is obvious; ``_removed/`` is a deliberate archive --
#: the 2.0 work *moved* the tracked .pyc files and HTESP.egg-info there on
#: purpose, so a sweep that treats them as build output destroys the record it
#: was meant to preserve (13 .pyc files, in this tree).
CLEAN_EXCLUDE = (".git", "_removed")


def clean_source_tree(root: str | None = None, dry_run: bool = False) -> int:
    """Remove build artifacts, returning a checkout to its pre-build state.

    ``build/``, ``*.egg-info/``, ``__pycache__/`` and stray ``*.pyc`` are
    generated by building or importing the package, are all in ``.gitignore``,
    and go stale in ways that mislead: this tree carries a ``build/`` from an
    aarch64 build and an ``htesp.egg-info/`` that no longer matches
    ``pyproject.toml``.  Nothing here is a source file, and nothing outside the
    checkout is touched.

    Deliberately does *not* remove the enumlib binaries, ``PMG_VASP_PSP_DIR`` or
    the stored API key: those are configuration of the machine, not artifacts
    of this tree, and rebuilding them is expensive.

    Refuses to run against an installed copy -- ``site-packages`` has no
    ``pyproject.toml``, and deleting ``__pycache__`` from under a live
    installation is never what anyone means by "clean the source tree".
    """
    here = Path(root).expanduser().resolve() if root else Path(__file__).resolve().parent.parent
    if "site-packages" in here.parts or "dist-packages" in here.parts:
        print("{} is an installed copy, not a source checkout; refusing to clean it."
              .format(here))
        print("  run this from the repository, or pass --root /path/to/HTESP")
        return 1
    if not (here / "pyproject.toml").is_file() or not (here / "htesp").is_dir():
        print("{} does not look like an HTESP checkout "
              "(no pyproject.toml beside an htesp/ package).".format(here))
        return 1

    verb = "would remove" if dry_run else "removed"
    removed = 0
    for name in BUILD_ARTIFACTS:
        target = here / name
        if not target.exists():
            continue
        print("  {} {}/".format(verb, name))
        if not dry_run:
            shutil.rmtree(target, ignore_errors=True)
        removed += 1

    caches = [d for d in here.rglob("__pycache__")
              if d.is_dir() and not set(d.parts) & set(CLEAN_EXCLUDE)]
    for cache in caches:
        if not dry_run:
            shutil.rmtree(cache, ignore_errors=True)
    if caches:
        print("  {} {} __pycache__/ director{}".format(
            verb, len(caches), "y" if len(caches) == 1 else "ies"))
        removed += len(caches)

    strays = [f for f in here.rglob("*.pyc")
              if not set(f.parts) & set(CLEAN_EXCLUDE)]
    for stray in strays:
        if not dry_run:
            try:
                stray.unlink()
            except OSError:                        # pragma: no cover
                pass
    if strays:
        print("  {} {} stray *.pyc".format(verb, len(strays)))
        removed += len(strays)

    if not removed:
        print("  nothing to remove; the tree is already clean")

    # Reported, never removed: these are not build artifacts, they are copy
    # damage, and deleting files is not something to do as a side effect.
    sidecars = [f for f in here.rglob("._*")
                if not set(f.parts) & set(CLEAN_EXCLUDE)]
    if sidecars:
        print("\n  note: {} macOS AppleDouble file(s) (._*) are also present. "
              "They are not\n  build output, so they are left alone; remove them "
              "with\n      find {} -name '._*' -delete".format(len(sidecars), here))
    return 0


def main(argv: list | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="htesp-check",
        description="Check that this machine can run HTESP "
                    "(architecture, interpreter, compiled dependencies).")
    parser.add_argument("--json", action="store_true",
                        help="print the report as JSON")
    parser.add_argument("--no-extras", action="store_true",
                        help="probe only the required dependencies")
    parser.add_argument("--executables", action="store_true",
                        help="also look for pw.x, vasp_std, sbatch, ...")
    parser.add_argument("--config_vasp_pot", metavar="DIR", default=None,
                        help="point pymatgen at a VASP POTCAR tree and exit; "
                             "give either POT_GGA_PAW_PBE or its parent")
    parser.add_argument("--clean", action="store_true",
                        help="remove build artifacts (build/, *.egg-info/, "
                             "__pycache__/, *.pyc) from the source checkout, "
                             "returning it to its pre-build state; leaves "
                             "enumlib and all configuration alone")
    parser.add_argument("--root", metavar="DIR", default=None,
                        help="the checkout --clean acts on (default: the one "
                             "this htesp package lives in)")
    parser.add_argument("--set_mp_api", metavar="KEY", default=None,
                        help="write the Materials Project API key to "
                             "~/.config/htesp/credentials and verify it; "
                             "survives new shells and batch jobs, unlike an "
                             "exported MP_API_KEY")
    parser.add_argument("--dry-run", action="store_true",
                        help="with --clean, list what would be removed and "
                             "remove nothing")
    parser.add_argument("--no-verify", action="store_true",
                        help="with --set_mp_api, skip the live check")
    parser.add_argument("--install-enumlib", action="store_true",
                        help="build enumlib from source and install enum.x and "
                             "makestr.x (needed by 'mainprogram magenum'); it "
                             "is not on PyPI, so pip cannot install it")
    parser.add_argument("--prefix", metavar="DIR", default=None,
                        help="where --install-enumlib puts the executables "
                             "(default: the bin/ of this Python environment)")
    parser.add_argument("--fortran-compiler", metavar="FC", default="gfortran",
                        help="compiler for --install-enumlib (default: gfortran)")
    args = parser.parse_args(argv)

    if args.clean:
        return clean_source_tree(args.root, dry_run=args.dry_run)
    if args.set_mp_api is not None:
        return set_mp_api_key(args.set_mp_api, verify=not args.no_verify)
    if args.config_vasp_pot:
        return configure_vasp_potcars(args.config_vasp_pot)
    if args.install_enumlib:
        return install_enumlib(args.prefix, args.fortran_compiler)

    result = report(extras=not args.no_extras)
    if args.json:
        print(json.dumps(result, indent=2))
    else:
        render(result, show_executables=args.executables)

    return 0 if all(entry["status"] == "ok"
                    for entry in result["required"]) else 1


if __name__ == "__main__":
    sys.exit(main())
