#!/usr/bin/env python
# coding: utf-8
# Written by Niraj K. Nepal, Ph.D.
"""
Program to write vasp inputfiles from poscars in mpid.vasp format
"""
import glob
import os
import shutil
import subprocess
import sys
from pymatgen.core import structure
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from htesp.cif_to_gsinput import pos_to_kpt
from htesp.write_potcar import poscar2potcar
from htesp.htepc import MpConnect
from htesp.check_json import config


def run_vasp_process(what, cwd=None):
    """Run ``vasp_process.py <what>`` without a shell, optionally in ``cwd``."""
    # FIX(36): ``vasp_process.py`` was an executable only while 1.x put src/ on
    # $PATH.  In 2.0 it is a module, so the bare name is a FileNotFoundError --
    # "[Errno 2] No such file or directory: 'vasp_process.py'" -- and every
    # command routed through here died.  workflow.py and oqmd_extract.py were
    # converted; these four call sites were missed.
    subprocess.run([sys.executable, "-m", "htesp.vasp_process", what],
                   check=True, cwd=cwd)


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
    input_data = config()
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
            relax_dir = "R{}-{}/relax".format(mpid,compound)
            os.makedirs(relax_dir, exist_ok=True)
            for name in ("POSCAR", "KPOINTS", "POTCAR", "INCAR"):
                shutil.move(name, os.path.join(relax_dir, name))
            if os.path.isfile("vasp.in"):
                shutil.copy("vasp.in", relax_dir)
            run_vasp_process("POSCAR", cwd=relax_dir)
        else:
            obj = MpConnect()
            struc = structure_standard
            compound = str(struc.composition).replace(" ","")
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
            # FIX(27): scf_dir is created by the VASP branch of every other
            # module but never on this one, so the move failed and the input
            # was left behind as scf-None.in in the working directory.
            os.makedirs("scf_dir", exist_ok=True)
            shutil.move("scf-{}.in".format(obj.mpid), f"scf_dir/scf-{mpid}.in")
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
