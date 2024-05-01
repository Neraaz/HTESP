#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to generate input files for substitutions"""
import sys
import os
from bsym.interface.pymatgen import unique_structure_substitutions as us
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io import pwscf
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import structure
from cif_to_gsinput import pos_to_kpt
from write_potcar import poscar2potcar
from htepc import MpConnect
from check_json import config

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
        except:
            obj.structure = pwscf.PWInput.from_file("scf_dir/scf-{}.in".format(mpid)).structure
        structure_sym = SpacegroupAnalyzer(obj.structure, symprec=0.1).get_primitive_standard_structure()
        if mode == 1:
            if isinstance(sb['sub'],dict):
                list_str = us(structure_sym,sb['elm'],sb['sub'])
            elif isinstance(sb['sub'],list):
                list_str = []
                for sub_i in sb['sub']:
                    list_str += us(structure_sym,sb['elm'],sub_i)
        elif mode == 2:
            new_sub = sb['new_sub']
            nelement = len(structure_sym)
            for i in range(nelement):
                try:
                    structure_sym[i].specie.symbol = new_sub[str(structure_sym[i].specie.symbol)]
                except:
                    continue
            #structure_sym.merge_sites()
            list_str = [structure_sym]
        else:
            print("wrong mode. Use either mode = 1 or 2 \n")
        if not os.path.isdir('scf_dir'):
            os.mkdir("scf_dir")
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
                if not os.path.isdir("R{}-{}-{}".format(mpid,i+1,obj.prefix)):
                    os.mkdir("R{}-{}-{}".format(mpid,i+1,obj.prefix))
                if not os.path.isdir("R{}-{}-{}/relax".format(mpid,i+1,obj.prefix)):
                    os.mkdir("R{}-{}-{}/relax".format(mpid,i+1,obj.prefix))
                pwd=os.getcwd()
                os.system("cp R{}-{}/relax/INCAR {}/R{}-{}-{}/relax/".format(mpid,orig_prefix,pwd,mpid,i+1,obj.prefix))
                poscar.write_file(filename="R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix))
                if os.path.isfile("config.json") or os.path.isfile("../../config.json"):
                    d = input_data['download']
                evenkpt = d['inp']['evenkpt']
                kptden = input_data['kptden']
                if evenkpt:
                    print("Even kpoint mesh is utilized\n")
                    pos_to_kpt("R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix),kptden,True)
                else:
                    pos_to_kpt("R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix),kptden)
                os.system("mv KPOINTS R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
                structure_file = structure.Structure.from_file("R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix))
                relax_set = MPRelaxSet(structure=structure_file)
                if os.path.isfile("vasp.in"):
                    os.system("cp vasp.in R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
                if os.path.isfile("config.json"):
                    os.system("cp config.json R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
                os.chdir("R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
                poscar2potcar()
                os.system("vasp_process.py POSCAR")
                os.chdir("../../")
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
                os.system("""mv scf-None.in""" + """ scf_dir/scf-{}-{}.in""".format(mpid,i+1))
            print(mpid,obj.prefix)
            with open("mpid-substitute.in", "a") as mpid_append:
                mpid_append.write("v{}".format(entry+1+i) + " " + mpid + "-{}".format(i+1) + " " + obj.prefix + "\n")
def main():
    """
    main function
    """
    mpid = sys.argv[1]
    #MPID_LIST = MPID.split("-")[0:2]
    #MPID = MPID_LIST[0] + "-" + MPID_LIST[1]
    #print(MPID)
    obj = MpConnect()
    substitution(mpid,obj)
if __name__ == "__main__":
    main()
