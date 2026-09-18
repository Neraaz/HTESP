#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to generate input files for substitutions"""
import os
import shutil
import subprocess
import sys
import warnings
from bsym.interface.pymatgen import unique_structure_substitutions as us
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io import pwscf
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import structure
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


def substitution(mpid,obj):
    """
    Substitute elements using the bsym package.

    This function substitutes elements using the bsym package. It requires a 'substitute' key in the config.json file.
    Inside that, the user needs to define 'mode', 'elm', 'sub', and 'new_sub' as explained below:

    Parameters:
    - mpid (str): Materials ID.
    - obj: Object of the MpConnect class.

    'substitute' dictionary:
    - 'mode' (int): 1 for unique_structure_substitutions, 2 for substituting each key by its value pair.
    - 'elm' (str): Element to be substituted.
    - 'sub' (dict or tuple): Dictionary for single substitution or tuple for multiple substitutions.
    - 'new_sub' (dict): Dictionary for new substitutions.

    If 'mode' is set to 1:
    - 'sub' (dict): Dictionary of elements to be substituted and their corresponding count.

    If 'mode' is set to 2:
    - 'new_sub' (dict): Dictionary of substitutions.

    Returns:
    None
    """
    input_data = config()
    if mpid in ('help', 'h'):
        msg="""required package: bsym 'pip install bsym'
               put substitute dictionary file in the config.json file. Inside that define
               'mode':2 # 'mode':1 for unique_structure_substitutions
               # 'mode':3 for new_structure_from_substitution function of bsym.interface.pymatgen
               'elm':'X', (element to be substituted)
               'sub':{'X':n1, 'Y':n2} (substitution)
               'new_sub':{'A1':'A2', 'B1':'B2', 'C1':'C2'}
               for multiple substitution, put 'sub':sub1,sub2,....
               suppose total number of X element is n
               then n2 of them is replaced by element Y"""
        print(msg)
    else:
        dft = sys.argv[2]
        orig_prefix=sys.argv[3]
        sb = input_data['substitute']
        mode = sb['mode']
        try:
            obj.structure = structure.Structure.from_file("R{}-{}/relax/POSCAR".format(mpid,orig_prefix))
        except (OSError, ValueError):
            obj.structure = pwscf.PWInput.from_file("scf_dir/scf-{}.in".format(mpid)).structure
        structure_sym = SpacegroupAnalyzer(obj.structure, symprec=0.1).get_primitive_standard_structure()
        list_str = []
        if mode == 1:
            if isinstance(sb['sub'],dict):
                list_str = us(structure_sym,sb['elm'],sb['sub'])
            elif isinstance(sb['sub'],list):
                list_str = []
                for sub_i in sb['sub']:
                    list_str += us(structure_sym,sb['elm'],sub_i)
        elif mode == 2:
            new_sub = sb['new_sub']
            # FIX(18): the old loop assigned to ``structure_sym[i].specie.symbol``.
            # ``specie`` is a pymatgen ``Element``, an enum singleton shared by
            # the whole process, so that renamed the element itself -- every
            # other structure, and every later material handled by the same
            # process, saw the substituted name.  ``replace_species`` changes
            # this structure and nothing else.
            present = {str(site.specie.symbol) for site in structure_sym}
            mapping = {old: new for old, new in new_sub.items() if old in present}
            missing = sorted(set(new_sub) - set(mapping))
            if missing:
                warnings.warn(
                    "substitute.new_sub names {} which are not in {}; "
                    "they are ignored".format(missing, structure_sym.composition),
                    RuntimeWarning)
            if mapping:
                structure_sym.replace_species(mapping)
            list_str = [structure_sym]
        else:
            print("wrong mode. Use either mode = 1 or 2 \n")
        os.makedirs("scf_dir", exist_ok=True)
        if not os.path.isfile('mpid-substitute.in'):
            entry = 0
        else:
            with open('mpid-substitute.in', 'r') as mpid_read:
                lines = mpid_read.readlines()
            entry = len(lines)
        for i,struc in enumerate(list_str):
            if dft in ('vasp', 'VASP'):
                obj.prefix = struc.composition.alphabetical_formula.replace(" ","")
                poscar = Poscar(structure=struc,comment=obj.prefix)
                relax_dir = "R{}-{}-{}/relax".format(mpid,i+1,obj.prefix)
                os.makedirs(relax_dir, exist_ok=True)
                shutil.copy("R{}-{}/relax/INCAR".format(mpid,orig_prefix), relax_dir)
                poscar.write_file(filename="{}/POSCAR".format(relax_dir))
                d = input_data['download']
                evenkpt = d['inp']['evenkpt']
                kptden = input_data['kptden']
                if evenkpt:
                    print("Even kpoint mesh is utilized\n")
                    pos_to_kpt("{}/POSCAR".format(relax_dir),kptden,True)
                else:
                    pos_to_kpt("{}/POSCAR".format(relax_dir),kptden)
                shutil.move("KPOINTS", "{}/KPOINTS".format(relax_dir))
                structure_file = structure.Structure.from_file("{}/POSCAR".format(relax_dir))
                relax_set = MPRelaxSet(structure=structure_file)
                if os.path.isfile("vasp.in"):
                    shutil.copy("vasp.in", relax_dir)
                if os.path.isfile("config.json"):
                    shutil.copy("config.json", relax_dir)
                # poscar2potcar() works on the current directory; the chdir is
                # wrapped so an exception cannot leave the process inside the
                # substituted directory.
                pwd = os.getcwd()
                try:
                    os.chdir(relax_dir)
                    poscar2potcar()
                    run_vasp_process("POSCAR")
                finally:
                    os.chdir(pwd)
            else:
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
                shutil.move("scf-{}.in".format(obj.mpid),
                            "scf_dir/scf-{}-{}.in".format(mpid,i+1))
            print(mpid,obj.prefix)
            with open("mpid-substitute.in", "a") as mpid_append:
                mpid_append.write("v{}".format(entry+1+i) + " " + mpid + "-{}".format(i+1) + " " + obj.prefix + "\n")


def main():
    """
    main function
    """
    mpid = sys.argv[1]
    obj = MpConnect()
    substitution(mpid,obj)


if __name__ == "__main__":
    main()
