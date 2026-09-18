#!/usr/bin/env python
# coding: utf-8
# Written by Niraj K. Nepal, Ph.D.
"""
Writing POTCAR file from POSCAR.

FIX(23): every third-party import in this module (``pymatgen``) is a core
dependency of the package, so nothing needs deferring here; the module is
imported by eight other modules and must stay cheap.
"""
import os
from pymatgen.core import structure
from pymatgen.io.vasp.sets import Potcar
from htesp.check_json import config
def poscar2potcar(poscar="POSCAR", outfile="POTCAR"):
    """
    Function to generate the POTCAR file based on the element potentials provided in
    config.json file.

    If config.json file is present in the current directory or its parent directory, the function retrieves
    the potential dictionary from the 'pseudo' key in the config.json file.

    If the config.json file is not found, a message indicating the absence of the 'pseudo' keyword is displayed.

    The function reads the POSCAR file to determine the elements present in the structure and generates the
    corresponding POTCAR file based on the potential dictionary.

    """
    input_data = config()
    # check_json.config() always returns a populated dict now, so the guard
    # that used to leave ``pot1`` unbound (and raise NameError one line later)
    # is gone; a configuration without a 'pseudo' section is reported instead.
    try:
        pot2 = input_data['pseudo']['pot']
    except (KeyError, TypeError) as exc:
        raise KeyError(
            "config.json has no 'pseudo'.'pot' mapping, so no POTCAR can be "
            "assembled") from exc
    struc = structure.Structure.from_file(poscar)
    element_list = []
    for elm in struc.composition.elements:
        try:
            element_list.append(pot2[str(elm)])
        except KeyError as exc:
            raise KeyError(
                "no pseudo.pot entry for element {!r} in config.json".format(
                    str(elm))) from exc
    potcar = Potcar(element_list)
    potcar.write_file(outfile)


def ensure_potcar(directory=".", poscar="POSCAR", outfile="POTCAR"):
    """Make sure ``directory`` has a POTCAR, building one from its POSCAR.

    FIX(38): POTCARs are licensed, so ``examples/`` cannot ship them -- the
    reference ``R<mpid>-<compound>/relax/`` directories contain INCAR, KPOINTS
    and POSCAR but no POTCAR.  Anything that stages those four files into a
    working directory (``mainprogram convtest`` most visibly) then died with
    ``FileNotFoundError: 'R<...>/relax/POTCAR'`` even though HTESP can build
    one itself from the POSCAR and ``pseudo.pot``.

    Does nothing when a POTCAR is already there.  Returns True when one exists
    afterwards, False when it could not be produced (no POSCAR, or pymatgen
    has no PMG_VASP_PSP_DIR -- the caller decides whether that is fatal).
    """
    directory = os.fspath(directory)
    target = os.path.join(directory, outfile)
    if os.path.isfile(target):
        return True
    source = os.path.join(directory, poscar)
    if not os.path.isfile(source):
        return False
    previous = os.getcwd()
    try:
        os.chdir(directory)
        poscar2potcar(poscar=poscar, outfile=outfile)
    except Exception as exc:                      # noqa: BLE001 - caller decides
        print("could not build {} from {}: {}: {}".format(
            target, source, type(exc).__name__, exc))
        return False
    finally:
        os.chdir(previous)
    return os.path.isfile(target)


#: printed once per run when a POTCAR is missing and cannot be built
POTCAR_HELP = (
    "VASP POTCARs are licensed, so HTESP cannot ship them and the reference\n"
    "  R<mpid>-<compound>/relax/ directories contain INCAR, KPOINTS and POSCAR\n"
    "  only.  To have HTESP build them from your own set:\n"
    "      htesp-check --config_vasp_pot /path/to/POT_GGA_PAW_PBE\n"
    "  (give either that directory or its parent).  Until then the inputs are\n"
    "  written without a POTCAR and VASP will not run on them."
)

_POTCAR_HELP_SHOWN = False


def stage_potcar(relax_dir, work_dir=None):
    """Put a POTCAR in ``relax_dir`` if one can be had; never fail if not.

    FIX(40): a missing POTCAR is not an error.  It cannot be shipped, it is not
    always configurable (a machine may simply have no VASP licence), and the
    rest of the input generation is still useful without it -- so this reports
    what to do and lets the caller carry on rather than aborting the scan.

    Returns True when ``relax_dir/POTCAR`` exists afterwards.  The explanatory
    message is printed at most once per process, so a 60-material scan does not
    repeat it 60 times.
    """
    global _POTCAR_HELP_SHOWN
    if ensure_potcar(relax_dir):
        return True
    print("no POTCAR in {}".format(relax_dir))
    if not _POTCAR_HELP_SHOWN:
        print("  " + POTCAR_HELP)
        _POTCAR_HELP_SHOWN = True
    return False
