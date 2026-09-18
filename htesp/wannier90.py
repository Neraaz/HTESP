#!/usr/bin/env python
#"""Written by Niraj K. Nepal, Ph.D"""
"""Module to prepare wannier90 input files"""
import os
import json
import shutil
import numpy as np
from ase.io import espresso
from htesp.kpoint_path import kpoint_path

# Number of Wannier functions contributed by one projection of each type.
_ORBITAL_SIZE = {"s": 1, "p": 3, "d": 5, "f": 7,
                 "sp": 2, "sp2": 3, "sp3": 4, "sp3d": 5, "sp3d2": 6}

# Fallback used when `num_wann` is absent and cannot be derived at all.
DEFAULT_NUM_WANN = None


def _as_bool(value):
    """
    Interpret a wannier90/JSON flag as a Python bool.

    FIX(9): ``plot_settings['bands_plot'] == '.true.'`` is a *string* compare.
    The shipped JSON happens to hold the string ``".true."``, but a JSON
    boolean ``true`` - which the schema also allows - silently failed the test
    and the whole ``Begin Kpoint_Path`` block was skipped.

    Parameters
    ----------
    value : bool or str or int or None
        Flag as read from the JSON file.

    Returns
    -------
    bool
        True for ``True``, ``'.true.'``, ``'true'``, ``'t'``, ``'.t.'``,
        ``'yes'``, ``1``; False otherwise.
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        return value.strip().strip(".").lower() in ("true", "t", "yes", "1")
    return False


def _count_projection_wannier(projection_file="projection.in"):
    """
    Count the Wannier functions requested by a ``projection.in`` file.

    Each non-empty line looks like ``X:s;p`` or ``f=0.25,0.25,0.25:s``; the
    orbitals after the final ``:`` are separated by ``;``.

    Parameters
    ----------
    projection_file : str
        Path to the projection file.

    Returns
    -------
    int or None
        The number of Wannier functions, or None if the file is absent or
        holds an orbital this table does not know.
    """
    if not os.path.isfile(projection_file):
        return None
    total = 0
    with open(projection_file, "r") as read_proj:
        for line in read_proj:
            line = line.split("!")[0].split("#")[0].strip()
            if not line or ":" not in line:
                continue
            for orb in line.rsplit(":", 1)[1].split(";"):
                orb = orb.strip().lower()
                if not orb:
                    continue
                if orb not in _ORBITAL_SIZE:
                    return None
                total += _ORBITAL_SIZE[orb]
    return total or None


def resolve_num_wann(config_settings, projection_file="projection.in"):
    """
    Work out ``num_wann`` for the wannier90 input.

    FIX(8): ``config_settings["num_wann"]`` raised KeyError with the shipped
    ``utility/input_files/wannier90.json``, which does not define the key, so
    every VASP wannier90 input failed.  The value is now looked up with
    ``.get`` and derived when it is missing.

    Resolution order:
      1. ``config_settings['num_wann']`` if present;
      2. the number of orbitals in ``projection_file``;
      3. ``config_settings['num_bands']`` (a crude upper bound - a warning is
         printed);
      4. otherwise a ValueError naming exactly what to add to the JSON file.

    Parameters
    ----------
    config_settings : dict
        The ``config_settings`` block of the wannier90 JSON file.
    projection_file : str
        Projection file used for step 2. Default: 'projection.in'.

    Returns
    -------
    int
        The number of Wannier functions.
    """
    num_wann = config_settings.get("num_wann", DEFAULT_NUM_WANN)
    if num_wann is not None:
        return int(num_wann)
    derived = _count_projection_wannier(projection_file)
    if derived is not None:
        print("num_wann not set in the wannier90 JSON file; "
              "derived {} from {}\n".format(derived, projection_file))
        return derived
    num_bands = config_settings.get("num_bands")
    if num_bands is not None:
        print("WARNING: num_wann not set and no usable {}; "
              "falling back to num_bands = {}\n".format(projection_file, num_bands))
        return int(num_bands)
    raise ValueError(
        "cannot determine 'num_wann': add \"num_wann\": <n> to the "
        "\"config_settings\" block of the wannier90 JSON file, or provide a "
        "'{}' file listing the projections".format(projection_file))


def epw_bandcheck(infile='scf.in', out="ex.win", proj=" ", json_file="wannier90.json", dft="QE"):
    """
    function to write input file for wannier90 calculations. Default: 'ex.win'

    parameters
    -------------
    infile : input file for scf calculation
    out : output file. Default: 'ex.win'
    proj : 'scdm' if used SCDM projection, '' otherwise
    json_file: JSON file containing element-value pairs for configurational settings
    dft: DFT package name, can be 'QE' (Quantum ESPRESSO), 'VASP', or 'other'
    """

    data = None
    if os.path.isfile("scf_dir/{}".format(infile)):
        data = espresso.read_espresso_in("scf_dir/{}".format(infile))
    else:
        print("Creates wannier-vasp.in from json file\n")

    # Load settings from combined JSON file
    with open(json_file, 'r') as read_json_file:
        all_settings = json.load(read_json_file)

    config_settings = all_settings["config_settings"]
    plot_settings = all_settings["plot_settings"]


    if dft in ("VASP","vasp"):
        # Write the WANNIER90_WIN file for VASP
        with open("wannier-vasp.in", 'w') as wannier_vasp_file:
            settings_str = " ".join(["{}={}\n".format(k, v) for k, v in config_settings.items()])

            # Include projection settings for VASP
            settings_str += "Begin Projections\n"

            if proj == 'scdm':
                settings_str += "auto_projections=.true.\n"
            elif proj == "fromfile":
                # FIX: os.system("cp ...") replaced by shutil.copy.
                if os.path.isfile("../../projection.in"):
                    shutil.copy("../../projection.in", "projection.in")
                if os.path.isfile("projection.in"):
                    with open("projection.in", "r") as gfile:
                        lines = gfile.readlines()
                    len_l = len(lines)
                    projections = ""
                    for line in lines:
                        projections += line
                    #projections = " ".join([line.strip() for line in lines])
                    settings_str += "{}".format(projections)
                else:
                    print("projection.in file not found\n")
                    print("write projections in different line, 'X:s', 'Y:pz', ... so on\n")
            else:
                settings_str += "random\n"
            settings_str += "End Projections\n"
            settings_str += " ".join(["{}={}\n".format(k, v) for k, v in plot_settings.items()])
            # FIX(9): accept a JSON boolean as well as '.true.'/'true'.
            if _as_bool(plot_settings.get('bands_plot')):
                kpoint_path(infile)
                settings_str += "Begin Kpoint_Path\n"
                with open("wannier_kpath.in", "r") as gfile:
                    lines = gfile.readlines()
                for line in lines:
                    settings_str += line
                settings_str += "END Kpoint_Path\n"
            # FIX(8): was config_settings["num_wann"] -> KeyError with the
            # shipped utility/input_files/wannier90.json.
            num_wann = resolve_num_wann(config_settings)
            wannier_vasp_file.write(f"NUM_WANN = {num_wann}\n")
            wannier_vasp_file.write("WANNIER90_WIN = \" \n")
            wannier_vasp_file.write("{}".format(settings_str))
            wannier_vasp_file.write("\"")

    elif dft in ("QE","qe"):
        if data is None:
            raise FileNotFoundError(
                "scf_dir/{} not found; the QE branch needs it to read the "
                "cell and atomic positions".format(infile))
        # Extract relevant data from the input file
        cell = data.cell
        symbol = data.get_chemical_symbols()
        pos = data.get_scaled_positions()
        kmesh = np.loadtxt('kmesh.grid')
        # Write the standard EPW input file
        with open(out, 'w') as epw_write:
            # Write configuration settings from JSON file
            for key, value in config_settings.items():
                epw_write.write("{} = {}\n".format(key, value))

            # Write Kpoint Path
            epw_write.write("Begin Kpoint_Path\n")
            with open("wannier_kpath.in", "r") as gfile:
                lines = gfile.readlines()
            for line in lines:
                epw_write.write(line)
            epw_write.write("End Kpoint_Path\n\n")

            # Write projections based on the 'proj' parameter
            if proj == 'scdm':
                epw_write.write("auto_projections = .true.\n")
            elif proj == "fromfile":
                if os.path.isfile("projection.in"):
                    with open("projection.in", "r") as gfile:
                        lines = gfile.readlines()
                    len_l = len(lines)
                    epw_write.write("begin projections\n")
                    for i in range(len_l):
                        epw_write.write("{}\n".format(lines[i].split("\n")[0]))
                    epw_write.write("end projections\n")
                else:
                    print("projection.in file not found\n")
                    print("write projections in different line, 'X:s', 'Y:pz', ... so on\n")
            else:
                epw_write.write("begin projections\n")
                epw_write.write("random\n")
                epw_write.write("end projections\n")

            # Write remaining data
            for key, value in plot_settings.items():
                epw_write.write("{} = {}\n".format(key, value))

            epw_write.write("begin unit_cell_cart\n")
            epw_write.write("Ang\n")
            for lat in cell:
                epw_write.write("{} {} {}\n".format(lat[0], lat[1], lat[2]))
            epw_write.write("end unit_cell_cart\n\n")
            epw_write.write("begin atoms_frac\n")
            for sym, pos in zip(symbol, pos):
                epw_write.write("{} {} {} {}\n".format(sym, pos[0], pos[1], pos[2]))
            epw_write.write("end atoms_frac\n\n")
            epw_write.write("mp_grid : {} {} {}\n".format(int(kmesh[0]), int(kmesh[1]), int(kmesh[2])))
            epw_write.write("begin kpoints\n")
            with open('wann_grid.out', 'r') as gfile:
                lines = gfile.readlines()
            for line in lines:
                epw_write.write(line)
            epw_write.write("end kpoints\n")

    else:
        print("Allowed values for 'dft' parameter are 'QE' or 'VASP'.")

if __name__ == "__main__":
    epw_bandcheck(dft="VASP")
