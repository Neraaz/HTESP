#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to write wanniertools input files"""
import sys
import os
import re
import numpy as np
import pandas as pd
from ase.io import espresso,vasp
from htesp.check_json import config
# FIX(7): this module used to carry its own copy of ``kpoint_path`` that
# duplicated htesp/kpoint_path.py (and silently lacked the VASP fallback and
# the discontinuity handling).  The shared one is imported instead; its
# ``slab=True`` flag covers the extra 2D path this module needs.
from htesp.kpoint_path import kpoint_path

# FIX(6): the default block below used to be built at *import* time from
# config['wanniertool_input'] - a key that does not exist (the real one is
# ``wanniertools_input``) - and the surrounding ``except FileNotFoundError``
# does not catch the resulting KeyError, so importing this module blew up and
# `wt1`/`wt2` were unusable.  The config is now read inside wt_body().
DEFAULT_WT_INPUT = {
    'tb_file': {'Hrfile': "'ex_hr.dat'", "Package": "'QE'"},
    'control': {"BulkBand_calc": "T", "BulkFS_calc": "F", "BulkGap_cube_calc": "F",
                "BulkGap_plane_calc": "F", "FindNodes_calc": "F", "SlabBand_calc": "F",
                "WireBand_calc": "F", "Dos_calc": "F", "JDos_calc": "F",
                "SlabSS_calc": "F", "SlabArc_calc": "F", "SlabQPI_calc": "F",
                "SlabSpintexture_calc": "F", "wanniercenter_calc": "F",
                "Z2_3D_calc": "F", "Chern_3D_calc": "F", "WeylChirality_calc": "F",
                "BerryPhase_calc": "F", "BerryCurvature_calc": "F", "AHC_calc": "F"},
    'parameters': {"E_arc": "0.0", "Eta_Arc": "0.001", "OmegaMin": "0.0",
                   "OmegaMax": "0.0", "OmegaNum": "100", "Nk1": "10", "Nk2": "10",
                   "Nk3": "10", "NP": "2", "Gap_threshold": "0.1"},
    'system': {"NSlab": "10", "NSlab1": "1", "NSlab2": "1", "NumOccupied": "1",
               "SOC": "1", "E_FERMI": "0.0", "Bx": "0.0", "By": "0.0", "Bz": "0.0",
               "surf_onsite": "0.0"},
    'surface': {"surface": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "KPATH_SLAB": {},
                "KPLANE_SLAB": {}, "EFFECTIVE_MASS": "0.0", "SELECTED_ATOMS": {}},
}


def read_wt_input():
    """
    Read the ``wanniertools_input`` block of config.json.

    FIX(6): the key was misspelled ``wanniertool_input`` and the read happened
    at module scope, so a missing/renamed key broke the *import* for every
    caller rather than the one function that needs the data.

    Returns
    -------
    dict
        The wanniertools input block, falling back to DEFAULT_WT_INPUT when
        config.json does not define one.
    """
    try:
        input_data = config()
    except (OSError, ValueError) as err:
        print("config.json could not be read ({}), using default "
              "wanniertools_input values\n".format(err))
        return DEFAULT_WT_INPUT
    wt_input = input_data.get('wanniertools_input')
    if not wt_input:
        print("'wanniertools_input' not found in config.json, "
              "using default values\n")
        return DEFAULT_WT_INPUT
    merged = {key: value for key, value in DEFAULT_WT_INPUT.items()}
    merged.update(wt_input)
    return merged


_WF_CENTRE = re.compile(
    r"\(\s*([-+0-9.eEdD]+)\s*,\s*([-+0-9.eEdD]+)\s*,\s*([-+0-9.eEdD]+)\s*\)")


def write_wannier_centres(wout_file, out="wannier.csv"):
    """
    Extract the final Wannier centres from a wannier90 ``.wout`` file.

    Replaces the ``sed ... | awk '{print $7 $8 $9}'`` pipeline that used to
    build ``wannier.csv``.

    Parameters
    ----------
    wout_file : str
        wannier90 output file, e.g. ``R<mpid>-<compound>/epw/ex.wout``.
    out : str
        CSV file to write. Default ``'wannier.csv'``.

    Returns
    -------
    int
        Number of Wannier centres written.
    """
    with open(wout_file, "r") as read_wout:
        lines = read_wout.readlines()
    centres = []
    inside = False
    for line in lines:
        if "Final State" in line:
            inside = True
            centres = []
            continue
        if not inside:
            continue
        if "Sum of centres and spreads" in line:
            inside = False
            continue
        match = _WF_CENTRE.search(line)
        if match:
            centres.append(match.groups())
    if not centres:
        raise ValueError(
            "no 'Final State' Wannier centres found in {}".format(wout_file))
    with open(out, "w") as write_csv:
        write_csv.write("x,y,z\n")
        for centre in centres:
            write_csv.write(",".join(centre) + "\n")
    return len(centres)


def wt_body(mpid,compound,surface=False):
    """
    Writes wt.in file in the format of 'wt-mpid-compound.in' inside WT_dir, based on information
    given in the wanniertools_input keyword within config.json file.

    Parameters:
    - mpid (str): Materials ID.
    - compound (str): Compound name used in mpid-list.in.
    - surface (bool): If True, includes cards printed for surface calculations. Default is False.

    Returns: None
    Example:
    ---------
    >>> wt_body("mp-1234567", "Si2") for bulk properties calculation
    """
    # FIX(6): the wanniertools block is read here, not at import time.
    wt_input = read_wt_input()
    # Check if the input file exists for either QE or VASP and load to structure
    if os.path.isfile("R{}-{}/relax/scf.in".format(mpid,compound)):
        input_file="R{}-{}/relax/scf.in".format(mpid,compound)
        data = espresso.read_espresso_in(input_file)
    elif os.path.isfile("R{}-{}/relax/POSCAR".format(mpid,compound)):
        input_file="R{}-{}/relax/POSCAR".format(mpid,compound)
        data = vasp.read_vasp(input_file)
    else:
        # FIX: `data` used to stay unbound here, so the next statement raised
        # NameError instead of reporting the missing structure file.
        raise FileNotFoundError(
            "no R{}-{}/relax/scf.in or POSCAR found".format(mpid, compound))
    # Extract wannier center data from ex.wout and save to wannier.csv
    # FIX: three chained `os.system` sed/awk calls that interpolated the
    # compound name into a shell string, replaced by a pure-Python parse.
    write_wannier_centres("R{}-{}/epw/ex.wout".format(mpid, compound))
    # Extract cell, symbols, and positions data
    cell=data.cell
    symbol=data.get_chemical_symbols()
    pos=data.get_scaled_positions()
    # Dictionary mapping orbital symbols to orbitals
    dict_orb = {'s':'s', 'p': 'pz px py', 'd': 'dz2 dxz dyz dx2-y2 dxy'}
    # Write data to wt-mpid-compound.in file
    with open("wt-{}-{}.in".format(mpid,compound), "w") as wtfile_write:
        # Write tb_file section from wanniertool_input
        wtfile_write.write("&TB_FILE\n")
        keys = list(wt_input['tb_file'].keys())
        values = list(wt_input['tb_file'].values())
        for tb_i,_ in enumerate(wt_input['tb_file']):
            wtfile_write.write(keys[tb_i] + "=" + values[tb_i] + "\n")
        wtfile_write.write("/\n")
        wtfile_write.write("\n")
        # Write control section from wanniertool_input
        wtfile_write.write("&CONTROL\n")
        keys = list(wt_input['control'].keys())
        values = list(wt_input['control'].values())
        for control_i,_ in enumerate(wt_input['control']):
            wtfile_write.write(keys[control_i] + "=" + values[control_i] + "\n")
        wtfile_write.write("/\n")
        wtfile_write.write("\n")
        # Write system section from wanniertool_input
        wtfile_write.write("&SYSTEM\n")
        keys = list(wt_input['system'].keys())
        values = list(wt_input['system'].values())
        for system_i,_ in enumerate(wt_input['system']):
            wtfile_write.write(keys[system_i] + "=" + values[system_i] + "\n")
        wtfile_write.write("/\n")
        wtfile_write.write("\n")
        # Write parameters section from wanniertool_input
        wtfile_write.write("&PARAMETERS\n")
        keys = list(wt_input['parameters'].keys())
        values = list(wt_input['parameters'].values())
        for param_i,_ in enumerate(wt_input['parameters']):
            wtfile_write.write(keys[param_i] + "=" + values[param_i] + "\n")
        wtfile_write.write("/\n")
        wtfile_write.write("\n")
        # Write lattice section
        wtfile_write.write("LATTICE\n")
        wtfile_write.write("Angstrom\n")
        for lat in cell:
            wtfile_write.write(str(lat[0]) + " " + str(lat[1]) + " " + str(lat[2]) + "\n")
        wtfile_write.write("\n")
        # Write atom positions section
        wtfile_write.write("ATOM_POSITIONS\n")
        wtfile_write.write(str(len(symbol)) + "\n")
        wtfile_write.write('Direct\n')
        for s_i,pos_i in enumerate(pos):
            wtfile_write.write(symbol[s_i] + " " + str(np.round(pos_i[0],10)) + " ")
            wtfile_write.write(str(np.round(pos_i[1],10)) + " " + str(np.round(pos_i[2],10)) + "\n")
        wtfile_write.write("\n")
        # Write projectors section
        wtfile_write.write("PROJECTORS\n")
        # Check if proj-wt.in file exists
        if os.path.isfile('proj-wt.in'):
            print("proj-wt.in file found, Reading .....")
            # Parse data and write to wt-mpid-compound.in file
            with open("proj-wt.in", "r") as gfile:
                lines = gfile.readlines()
            dict_elm = {}
            for line in lines:
                dict_elm[line.split('\n')[0].split()[0]] = line.split('\n')[0].split()[1:]
            for sym in symbol:
                norb = 0
                sorb = dict_elm[sym]
                for orb in sorb:
                    if orb in dict_orb.keys():
                        norb += len(dict_orb[orb].split())
                    else:
                        norb += 1
                wtfile_write.write(str(norb) + " ")
            wtfile_write.write("!number of projectors\n")
            for sym in symbol:
                orb_str = ''
                sorb = dict_elm[sym]
                for orb in sorb:
                    if orb in dict_orb.keys():
                        orb_str = orb_str + " " +  dict_orb[orb]
                    else:
                        orb_str = orb_str + " " + orb
                wtfile_write.write(sym + " " + orb_str + "\n")
        else:
            # Write default s orbital if proj-wt.in file not found
            print("proj-wt.in file not found, writing only s orbital .....")
            for sym in symbol:
                wtfile_write.write("1 ")
            wtfile_write.write("!number of projectors\n")
            for sym in symbol:
                wtfile_write.write(sym + " s" + "\n")
        wtfile_write.write("\n")
        # Write surface section
        if surface:
            wtfile_write.write("SURFACE\n")
            # Check if config.json or surface-wt.in file exists
            if os.path.isfile("config.json") or os.path.isfile("../../config.json"):
                print("Taking surface from config.json\n")
                for surf_lat in wt_input['surface']['surface']:
                    wtfile_write.write(str(surf_lat[0]) + " " + str(surf_lat[1]) + " " + str(surf_lat[2]) + "\n")
            # Read surface data from surface-wt.in file
            elif os.path.isfile("surface-wt.in"):
                print("surface-wt.in file found")
                with open('surface-wt.in', 'r') as gfile:
                    lines = gfile.readlines()
                for line in lines:
                    wtfile_write.write(line)
            else:
                print("surface-wt.in and surface inside wanniertools_input not found, printing default values")
                wtfile_write.write("1 0 0 ! write vectors determining surfaces\n")
                wtfile_write.write("0 1 0 \n")
        wtfile_write.write("\n")
        # Write bulk k-path section
        wtfile_write.write("KPATH_BULK\n")
        with open("R{}-{}/epw/ex.win".format(mpid,compound), "r") as hfile:
            lines = hfile.readlines()
        # FIX: `os.system("sed -n '/Begin/,/End/p' ... | sed 1d | sed $d > ...")`
        # replaced by a Python slice; the old test also only matched when the
        # markers had no surrounding whitespace.
        begin = end = None
        for iline, line in enumerate(lines):
            if begin is None and "Begin Kpoint_Path" in line:
                begin = iline
            elif begin is not None and "End Kpoint_Path" in line:
                end = iline
                break
        if begin is not None and end is not None:
            with open("wannier_kpath.in", "w") as kpath_write:
                for line in lines[begin + 1:end]:
                    kpath_write.write(line)
        with open("wannier_kpath.in", "r") as gfile:
            lines = gfile.readlines()
        wtfile_write.write(str(len(lines)) + "\n")
        for i,line in enumerate(lines):
            wtfile_write.write(line)
        wtfile_write.write("\n")
        # Write slab k-path section if surface is True
        if surface:
            if os.path.isfile("wannier_kpath_slab.in"):
                with open("wannier_kpath_slab.in", "r") as wanpath:
                    kpath_2d = wanpath.readlines()
                wtfile_write.write("KPATH_SLAB\n")
                wtfile_write.write(str(len(kpath_2d)) + "\n")
                for line in kpath_2d:
                    wtfile_write.write(line)
                wtfile_write.write("\n")
        # Write wannier centers section
        wtfile_write.write("WANNIER_CENTRES\n")
        wtfile_write.write("Cartesian\n")
        # Read wannier.csv file and write to wt-mpid-compound.in file
        wannier = pd.read_csv("wannier.csv")
        for i in range(wannier.shape[0]):
            wtfile_write.write(str(wannier.x[i]) + " " + str(wannier.y[i]) + " " + str(wannier.z[i]) + "\n")
def main():
    """
    Main function to execute different operations based on command-line arguments.

    Command-line arguments:
    - mpid (str): Material ID.
    - compound (str): Compound name.
    - prefix (str): Prefix value.
    - condition (str): Condition for different operations.
    - surface_prop (str): Surface property indicator.

    Returns: None
    """
    # Extract command-line arguments
    mpid = sys.argv[1]
    compound = sys.argv[2]
    prefix = sys.argv[3]
    condition = sys.argv[4]
    if condition == "kpathwan":
        # Check if scf file exists for the specified material and compound
        file_name = "scf_dir/scf-{}-{}.in".format(mpid, compound)
        # FIX: `file_name` used to stay unbound on the miss -> NameError.
        if not os.path.isfile(file_name):
            print("scf-{}-{}.in not found in scf_dir/\n".format(mpid, compound))
            print("Complete mainprogram.py from 1 to 4\n")
            return
        # FIX(7): the shared htesp.kpoint_path.kpoint_path, with slab=True so
        # it still writes wannier_kpath_slab.in when POSCAR-slab is present.
        kpoint_path(file_name, slab=True)
    elif condition == "body":
        # Determine if surface property should be included
        surface_prop = sys.argv[5]
        if surface_prop == 'T':
            wt_body(mpid,compound,surface=True)
        else:
            wt_body(mpid,compound,surface=False)
    else:
        print("bad inputs")
if __name__ == "__main__":
    main()
