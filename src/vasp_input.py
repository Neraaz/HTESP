#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D.
This script will download vasp input files from materials project database.
Script is run within 'download-input' bash script."""
import os
import sys
import json
import warnings
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from cif_to_gsinput import pos_to_kpt
from write_potcar import poscar2potcar
from htepc import MpConnect
# To ignore all warnings
warnings.filterwarnings("ignore")
try:
    PWD = os.getcwd()
    if os.path.isfile(PWD+"/htepc.json"):
        JSONFILE = PWD+"/htepc.json"
    else:
        JSONFILE = "../../htepc.json"
    with open(JSONFILE, "r") as readjson:
        input_data = json.load(readjson)
except FileNotFoundError:
    print("htepc.json file not found\n")
def vasp_input(mpid,compound):
    """
    Download VASP input files from the Materials Project.
    Downloads INCAR, POSCAR, POTCAR, and KPOINTS inside Rmpid-compound/relax folder.
    Updates 'mpid.in' with entry number, mpid, and compound name.

    Parameters:
    - mpid (str): Materials ID.
    - compound (str): Compound name.
    """
    obj = MpConnect()
    obj.setting(mpid)
    obj.download()
    if not os.path.isdir("input_cif"):
        os.mkdir("input_cif")
    if os.path.isfile("{}.cif".format(mpid)):
        os.system("mv {}.cif input_cif/".format(mpid))
    #Create poscar separately
    #from pymatgen.io.vasp.inputs import Poscar
    #poscar = Poscar(structure=struct,comment="MgB2")
    #poscar = Poscar(structure=obj.structure,comment="MgB2")
    #poscar
    #poscar.write_file(filename="POSCAR")
    #obtain vasp inputs from MPRelaxSet
    structure = SpacegroupAnalyzer(obj.structure, symprec=0.1).get_primitive_standard_structure()
    #relax_set = MPRelaxSet(structure=obj.structure)
    relax_set = MPRelaxSet(structure=structure)
    #process 1
    #For POTCAR. Suppose we have POTCARS as  POT_GGA_PAW_PBE/Mg_p/POTCAR
    #Mg_p can be found in POTCAR.spec
    #pmg config -p /home/nnepal/bin/POT_GGA_PAW_PBE PBE52
    # After that add path to .pmgrc.yaml
    #pmg config --add PMG_VASP_PSP_DIR PBE52
    #write INCAR and POSCAR separately
    #incar = relax_set.incar
    #incar.write_file('INCAR')
    #create a folder with structures
    if not os.path.isdir("R{}-{}".format(mpid,compound)):
        os.mkdir("R{}-{}".format(mpid,compound))
    if not os.path.isdir("R{}-{}/relax".format(mpid,compound)):
        os.mkdir("R{}-{}/relax".format(mpid,compound))
    relax_set.write_input(output_dir="R{}-{}/relax".format(mpid,compound))
    #once we complete process 1, we can write POTCAR file as
    #relax_set.potcar.write_file("R{}-{}/relax/POTCAR".format(mpid,compound))
    relax_set.poscar.write_file("R{}-{}/relax/POSCAR".format(mpid,compound))
    if os.path.isfile("htepc.json"):
        os.system("cp htepc.json R{}-{}/relax/".format(mpid,compound))
    os.chdir("R{}-{}/relax/".format(mpid,compound))
    poscar2potcar()
    os.chdir("../../")
    if os.path.isfile("htepc.json") or os.path.isfile("../../htepc.json"):
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
    #MPID_LIST = MPID.split("-")[0:2]
    #MPID = MPID_LIST[0] + "-" + MPID_LIST[1]
    compound = sys.argv[2]
    vasp_input(mpid,compound)
if __name__ == "__main__":
    main()
