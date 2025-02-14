#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D.
"""This script will download vasp input files from materials project database.
Script is run within 'download-input' bash script."""
import os
import sys
import warnings
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from cif_to_gsinput import pos_to_kpt
from write_potcar import poscar2potcar
from htepc import MpConnect
from check_json import config
# To ignore all warnings
warnings.filterwarnings("ignore")
def vasp_input(mpid,compound):
    """
    Download VASP input files from the Materials Project.
    Downloads INCAR, POSCAR, POTCAR, and KPOINTS inside Rmpid-compound/relax folder.
    Updates 'mpid.in' with entry number, mpid, and compound name.

    Parameters:
    - mpid (str): Materials ID.
    - compound (str): Compound name.
    """
    input_data = config()
    obj = MpConnect()
    obj.setting(mpid)
    obj.download()
    if not os.path.isdir("input_cif"):
        os.mkdir("input_cif")
    if os.path.isfile("{}.cif".format(mpid)):
        os.system("mv {}.cif input_cif/".format(mpid))
    #obtain vasp inputs from MPRelaxSet
    structure = SpacegroupAnalyzer(obj.structure, symprec=0.1).get_primitive_standard_structure()
    relax_set = MPRelaxSet(structure=structure)
    #create a folder with structures
    if not os.path.isdir("R{}-{}".format(mpid,compound)):
        os.mkdir("R{}-{}".format(mpid,compound))
    if not os.path.isdir("R{}-{}/relax".format(mpid,compound)):
        os.mkdir("R{}-{}/relax".format(mpid,compound))
    relax_set.write_input(output_dir="R{}-{}/relax".format(mpid,compound))
    relax_set.poscar.write_file("R{}-{}/relax/POSCAR".format(mpid,compound))
    if os.path.isfile("config.json"):
        os.system("cp config.json R{}-{}/relax/".format(mpid,compound))
    os.chdir("R{}-{}/relax/".format(mpid,compound))
    poscar2potcar()
    os.chdir("../../")
    if os.path.isfile("config.json") or os.path.isfile("../../config.json"):
        d = input_data['download']
    evenkpt = d['inp']['evenkpt']
    kptden = input_data['kptden']
    if evenkpt:
        print("Even kpoint mesh is utilized\n")
        pos_to_kpt("R{}-{}/relax/POSCAR".format(mpid,compound),kptden,True)
    else:
        pos_to_kpt("R{}-{}/relax/POSCAR".format(mpid,compound),kptden)
    os.system("mv KPOINTS R{}-{}/relax/".format(mpid,compound))
    if not os.path.isfile('mpid.in'):
        entry = 0
        with open("mpid.in", "w") as write_mpid:
            write_mpid.write("v{}".format(entry+1) + " " + obj.mpid + " " + compound + "\n")
    else:
        with open('mpid.in', 'r') as read_mpid:
            lines = read_mpid.readlines()
        entry = len(lines)
        new_mpid = mpid + " "
        if not any(new_mpid in line for line in lines):
            with open("mpid.in", "a") as write_mpid:
                write_mpid.write("v{}".format(entry+1) + " " + obj.mpid + " " + compound + "\n")
def main():
    """
    main function
    """
    mpid = sys.argv[1]
    compound = sys.argv[2]
    vasp_input(mpid,compound)
if __name__ == "__main__":
    main()
