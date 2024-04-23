#!/usr/bin/env python
#"""Written by Niraj K. Nepal, Ph.D"""
"""Module to prepare wannier90 input files"""
import os
import json
import numpy as np
from ase.io import espresso
from kpoint_path import kpoint_path

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
                if os.path.isfile("../../projection.in"):
                    os.system("cp ../../projection.in .")
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
            if plot_settings['bands_plot'] == '.true.':
                kpoint_path(infile)
                settings_str += "Begin Kpoint_Path\n"
                with open("wannier_kpath.in", "r") as gfile:
                    lines = gfile.readlines()
                for line in lines:
                    settings_str += line
                settings_str += "END Kpoint_Path\n"
            num_wann = config_settings["num_wann"]
            wannier_vasp_file.write(f"NUM_WANN = {num_wann}\n")
            wannier_vasp_file.write("WANNIER90_WIN = \" \n")
            wannier_vasp_file.write("{}".format(settings_str))
            wannier_vasp_file.write("\"")

    elif dft in ("QE","qe"):
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
