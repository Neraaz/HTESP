#!/usr/bin/env python
# coding: utf-8
"""
Creates different magnetic ordering according to MagneticStructureEnumerator functions.
"""
import os
import json
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core import structure
from pymatgen.io.pwscf import PWInput
from pymatgen.analysis.magnetism import MagneticStructureEnumerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from cif_to_gsinput import pos_to_kpt
from write_potcar import poscar2potcar
from htepc import MpConnect
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
def magnetic_structure(obj,mpid,compound,magconfig,dft):
    """
    Functions to generate possible magnetic structures.

    Parameters:
    obj (MpConnect): An object representing an instance of the MpConnect class.
    mpid (str): Materials ID.
    compound (str): Compound name.
    magconfig (list): List of possible collinear magnetic configurations, e.g., ['ferromagnetic', 'antiferromagnetic'].
    dft (str): Type of DFT calculation, e.g., 'VASP' or 'QE'.

    Look for pymatgen.analysis.magnetism.MagneticStructureEnumerator class for more details.
    Install Enumlib library: https://github.com/msg-byu/enumlib to run this module.
    """
    try:
        strucinit = structure.Structure.from_file("R{}-{}/relax/POSCAR".format(mpid,compound))
    except FileNotFoundError:
        strucinit = PWInput.from_file("R{}-{}/relax/scf.in".format(mpid,compound)).structure
    if os.path.isfile("htepc.json"):
        #from magmom import magmom, order
        default_magmoms = input_data['magmom']['magmom']
        order = input_data['magmom']['order']
        #print(magmom,order)
        newstructure = MagneticStructureEnumerator(strucinit,default_magmoms=default_magmoms,strategies=order,truncate_by_symmetry=True).ordered_structures
    else:
        newstructure = MagneticStructureEnumerator(strucinit,strategies=magconfig,truncate_by_symmetry=True).ordered_structures
    print(len(newstructure))
    if not os.path.isfile('mpid-magnetic.in'):
        entry = 0
    else:
        with open('mpid-magnetic.in', 'r') as mpid_read:
            lines = mpid_read.readlines()
        entry = len(lines)
    if len(newstructure) > 0:
        for i,struc in enumerate(newstructure):
            struc = SpacegroupAnalyzer(struc,symprec=0.1).get_refined_structure()
            if dft in ('vasp', 'VASP'):
                obj.prefix = struc.composition.alphabetical_formula.replace(" ","")
                poscar = Poscar(structure=struc,comment=obj.prefix)
                if not os.path.isdir("R{}-{}-{}".format(mpid,i+1,obj.prefix)):
                    os.mkdir("R{}-{}-{}".format(mpid,i+1,obj.prefix))
                if not os.path.isdir("R{}-{}-{}/relax".format(mpid,i+1,obj.prefix)):
                    os.mkdir("R{}-{}-{}/relax".format(mpid,i+1,obj.prefix))
                pwd=os.getcwd()
                orig_prefix = compound
                os.system("cp R{}-{}/relax/INCAR {}/R{}-{}-{}/relax/".format(mpid,orig_prefix,pwd,mpid,i+1,obj.prefix))
                poscar.write_file(filename="R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix))
                #os.system("cp R{}-{}/relax/KPOINTS {}/R{}-{}-{}/relax/".format(mpid,orig_prefix,pwd,mpid,i+1,obj.prefix))
                #if os.path.isfile("download.py"):
                #    import download as d
                #d = input_data['download']
                evenkpt = input_data['download']['inp']['evenkpt']
                if evenkpt:
                    print("Even kpoint mesh is utilized\n")
                    pos_to_kpt("R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix),0.025,True)
                else:
                    pos_to_kpt("R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix),0.025)
                os.system("mv KPOINTS R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
                #structure_file = structure.Structure.from_file("R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix))
                #relax_set = MPRelaxSet(structure=structure_file)
                #relax_set.potcar.write_file("R{}-{}-{}/relax/POTCAR".format(mpid,i+1,obj.prefix))
                if os.path.isfile("vasp.in"):
                    os.system("cp vasp.in R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
                #if os.path.isfile("magmom.py"):
                #    os.system("cp magmom.py R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
                if os.path.isfile("htepc.json"):
                    os.system("cp htepc.json R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
                os.chdir("R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
                poscar2potcar()
                os.system("vasp_process.py POSCAR")
                with open("INCAR", "r") as read_incar:
                    incar = read_incar.readlines()
                for j,line in enumerate(incar):
                    if 'MAGMOM' in line:
                        os.system("""sed -i '{}d' INCAR""".format(j+1))
                        os.system("""sed -i '/ISPIN/d' INCAR""")
                maglist = ""
                for _,specie in enumerate(struc.species):
                    if specie.spin is not None:
                        maglist += str(specie.spin) + " "
                    else:
                        maglist += "0 "
                os.system("""echo "ISPIN = 2" >> INCAR""")
                os.system("""echo "MAGMOM = {}" >> INCAR""".format(maglist))
                os.chdir("../../")
            else:
                obj.structure = struc
                comp_list = []
                for composition in struc.composition.elements:
                    comp_list.append(composition.name)
                obj.comp_list = comp_list
                obj.getkpt()
                evenkpt = input_data['download']['inp']['evenkpt']
                if evenkpt:
                    print("Utilizing even kpoint mesh\n")
                    obj.getevenkpt()
                obj.maxecut_sssp_for_subs()
                obj.prefix = struc.composition.alphabetical_formula.replace(" ","")
                obj.setting_qeinput(magnetic=True,pseudo_dir='../../pp')
                if not os.path.isdir("scf_dir"):
                    os.system("mkdir scf_dir")
                os.system("""mv scf-None.in""" + """ scf_dir/scf-{}-{}.in""".format(mpid,i+1))
            #print(mpid,obj.prefix)
            with open("mpid-magnetic.in", "a") as mpid_append:
                mpid_append.write("v{}".format(entry+1+i) + " " + mpid + "-{}".format(i+1) + " " + obj.prefix + "\n")
                #mpid_append.write("v{}".format(entry+1+i) + " " + obj.mpid  + " " + obj.prefix + "\n")
    else:
        print("Structures not created\n")
def changeaxis(mpid,comp):
    """
    Function to change magnetic axis to compute magnetic anisotropy.

    Parameters:
    mpid (str): Material ID.
    comp (str): Compound name.

    The function iterates over magnetic axis configurations specified in the saxis variable,
    updates the necessary input files for VASP calculations, and performs the VASP process.
    """
    if not os.path.isfile('mpid-magnetic.in'):
        entry = 0
    else:
        with open('mpid-magnetic.in', 'r') as mpid_read:
            lines = mpid_read.readlines()
        entry = len(lines)
    saxis = input_data['magmom']['saxis']
    for i,axis in enumerate(saxis):
        sx = int(axis[0])
        sy = int(axis[1])
        sz = int(axis[2])
        if not os.path.isdir("R{}-saxis-{}{}{}-{}".format(mpid,sx,sy,sz,comp)):
            os.mkdir("R{}-saxis-{}{}{}-{}".format(mpid,sx,sy,sz,comp))
        if not os.path.isdir("R{}-saxis-{}{}{}-{}/relax".format(mpid,sx,sy,sz,comp)):
            os.mkdir("R{}-saxis-{}{}{}-{}/relax".format(mpid,sx,sy,sz,comp))
        os.chdir("R{}-{}/relax/".format(mpid,comp))
        os.system("cp INCAR POTCAR POSCAR KPOINTS ../../R{}-saxis-{}{}{}-{}/relax/".format(mpid,sx,sy,sz,comp))
        os.chdir("../../")
        os.system("""sed -i '/SAXIS/d' R{}-saxis-{}{}{}-{}/relax/INCAR""".format(mpid,sx,sy,sz,comp))
        os.system("""echo 'SAXIS = {} {} {}' >> R{}-saxis-{}{}{}-{}/relax/INCAR""".format(sx,sy,sz,mpid,sx,sy,sz,comp))
        with open("mpid-magnetic.in", "a") as mpid_append:
            mpid_append.write("v{}".format(entry+1+i) + " " + mpid + "-saxis-{}{}{}".format(sx,sy,sz) + " " + comp + "\n")
        if os.path.isfile("htepc.json"):
            os.system("""cp htepc.json R{}-saxis-{}{}{}-{}/relax/""".format(mpid,sx,sy,sz,comp))
        if os.path.isfile("vasp.in"):
            os.system("""cp vasp.in R{}-saxis-{}{}{}-{}/relax/""".format(mpid,sx,sy,sz,comp))
        os.chdir("R{}-saxis-{}{}{}-{}/relax/""".format(mpid,sx,sy,sz,comp))
        os.system("vasp_process.py POSCAR")
        os.chdir("../../")
def main():
    """
    Main function to control the workflow.

    Reads input data, determines the type of calculation and magnetic configuration,
    and performs the corresponding operations.

    Parameters:
    None

    Returns:
    None
    """
    #mpid = sys.argv[1]
    #comp = sys.argv[2]
    #dft = sys.argv[3]
    dft = input_data['download']['inp']['calc']
    mag_type = input_data['magmom']['type']
    with open("input.in","r") as read_in:
        lines = read_in.readlines()
    start = int(lines[0].split("\n")[0])
    end = int(lines[1].split("\n")[0])
    filename = lines[3].split("\n")[0]
    with open(filename,'r') as read_mpid:
        lines = read_mpid.readlines()
    lines = lines[start-1:end-1]
    obj = MpConnect()
    config = ['ferromagnetic','antiferromagnetic']
    for line in lines:
        mpid = line.split(" ")[1]
        comp = line.split(" ")[2].split("\n")[0]
        if mag_type == "ordering":
            magnetic_structure(obj,mpid,comp,config,dft)
        elif mag_type == "anisotropy":
            changeaxis(mpid,comp)
        else:
            print("Only ordering and anisotropy allowed\n")
if __name__ == "__main__":
    main()
