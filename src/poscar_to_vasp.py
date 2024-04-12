#!/usr/bin/env python
# coding: utf-8
# Written by Niraj K. Nepal, Ph.D.
"""
Program to write vasp inputfiles from poscars in mpid.vasp format
"""
import os
import glob
import json
from pymatgen.core import structure
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from cif_to_gsinput import pos_to_kpt
from write_potcar import poscar2potcar
from htepc import MpConnect
try:
    PWD = os.getcwd()
    if os.path.isfile(PWD+"/config.json"):
        JSONFILE = PWD+"/config.json"
    else:
        JSONFILE = "../../config.json"
    with open(JSONFILE, "r") as readjson:
        input_data = json.load(readjson)
except FileNotFoundError:
    print("config.json file not found\n")
def main():
    """
    Main function to prepare VASP input files and directories for relaxation.

    This function performs the following tasks:
    - Searches for VASP structure files in the current directory.
    - Converts the structure to the primitive standard form and generates VASP input files.
    - Sets up directories for relaxation calculations and organizes input files accordingly.
    - Appends the information of each structure to 'mpid.in' file.

    Note:
    -----
    Ensure that the necessary VASP input files (POSCAR, INCAR, KPOINTS, POTCAR) are available.
    The 'vasp_process.py' script is assumed to be available for processing the POSCAR file.

    Returns:
    --------
    None
    """
    structures = glob.glob("*.vasp",recursive=True)
    kptden = input_data['kptden']
    dft = input_data['download']['inp']['calc']
    for struc in structures:
        mpid = struc.split(".")[0]
        struc_poscar = structure.Structure.from_file(struc)
        structure_standard = SpacegroupAnalyzer(structure=struc_poscar,symprec=0.1).get_primitive_standard_structure()
        structure_standard.to("POSCAR")
        if dft in ('VASP','vasp'):
            relax_set = MPRelaxSet(structure_standard)
            pos_to_kpt("POSCAR",kptden)
            poscar2potcar()
            relax_set.incar.write_file("INCAR")
            compound = str(struc_poscar.composition).replace(" ","")
            if not os.path.isdir("R{}-{}".format(mpid,compound)):
                os.mkdir("R{}-{}".format(mpid,compound))
            if not os.path.isdir("R{}-{}/relax".format(mpid,compound)):
                os.mkdir("R{}-{}/relax".format(mpid,compound))
            os.system("mv POSCAR KPOINTS POTCAR INCAR R{}-{}/relax".format(mpid,compound))
            if os.path.isfile("vasp.in"):
                os.system("cp vasp.in R{}-{}/relax".format(mpid,compound))
            os.chdir("R{}-{}/relax".format(mpid,compound))
            os.system("vasp_process.py POSCAR")
            os.chdir("../../")
        else:
            obj = MpConnect()
            struc = structure_standard
            magnetic = input_data['pwscf_in']['magnetic']
            if magnetic:
                default_magmoms = input_data['magmom']['magmom']
                struc.add_spin_by_element(default_magmoms)
                obj.structure = struc
            else:
                obj.structure = struc
            comp_list = []
            for composition in struc.composition.elements:
                comp_list.append(str(composition))
            obj.comp_list = comp_list
            obj.getkpt()
            evenkpt = input_data['download']['inp']['evenkpt']
            if evenkpt:
                print("Utilizing even kpoint mesh\n")
                obj.getevenkpt()
            obj.maxecut_sssp_for_subs()
            obj.prefix = struc.composition.alphabetical_formula.replace(" ","")
            if magnetic:
                obj.setting_qeinput(magnetic=True,pseudo_dir='../../pp')
            else:
                obj.setting_qeinput(pseudo_dir='../../pp')
            os.system("""mv scf-None.in""" + f""" scf_dir/scf-{mpid}.in""")
            print(mpid,obj.prefix)
        if not os.path.isfile('mpid.in'):
            entry = 0
            with open("mpid.in", "w") as write_mpid:
                write_mpid.write("v{}".format(entry+1) + " " + mpid + " " + compound + "\n")
        else:
            with open('mpid.in', 'r') as read_mpid:
                lines = read_mpid.readlines()
            entry = len(lines)
            new_mpid = mpid + " "
            if not any(new_mpid in line for line in lines):
                with open("mpid.in", "a") as write_mpid:
                    write_mpid.write("v{}".format(entry+1) + " " + mpid + " " + compound + "\n")
if __name__ == "__main__":
    main()
