#!/usr/bin/env python
# coding: utf-8
# Written by Niraj K. Nepal, Ph.D.
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
    if os.path.isfile(PWD+"/config.json"):
        JSONFILE = PWD+"/config.json"
    else:
        JSONFILE = "../../config.json"
    with open(JSONFILE, "r") as readjson:
        input_data = json.load(readjson)
except FileNotFoundError:
    print("config.json file not found\n")
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

    Example:
    >>> obj = MpConnect()
    >>> magnetic_structure(obj, 'mp-123', 'MnO', ['ferromagnetic', 'antiferromagnetic'], 'VASP')
    """
    # load the structure
    try:
        strucinit = structure.Structure.from_file("R{}-{}/relax/POSCAR".format(mpid,compound))
    except FileNotFoundError:
        strucinit = PWInput.from_file("R{}-{}/relax/scf.in".format(mpid,compound)).structure
    # Obtain magnetic structures based on given magnetic configurations
    if os.path.isfile("config.json"):
        default_magmoms = input_data['magmom']['magmom']
        order = input_data['magmom']['order']
        newstructure = MagneticStructureEnumerator(strucinit,default_magmoms=default_magmoms,strategies=order,truncate_by_symmetry=True).ordered_structures
    else:
        newstructure = MagneticStructureEnumerator(strucinit,strategies=magconfig,truncate_by_symmetry=True).ordered_structures
    print(len(newstructure))
    # Check the entry number
    if not os.path.isfile('mpid-magnetic.in'):
        entry = 0
    else:
        with open('mpid-magnetic.in', 'r') as mpid_read:
            lines = mpid_read.readlines()
        entry = len(lines)
    # Process each generated structure
    if len(newstructure) > 0:
        for i,struc in enumerate(newstructure):
            # Refine the structure
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
                evenkpt = input_data['download']['inp']['evenkpt']
                if evenkpt:
                    print("Even kpoint mesh is utilized\n")
                    pos_to_kpt("R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix),0.025,True)
                else:
                    pos_to_kpt("R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix),0.025)
                os.system("mv KPOINTS R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
                if os.path.isfile("vasp.in"):
                    os.system("cp vasp.in R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
                if os.path.isfile("config.json"):
                    os.system("cp config.json R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
                os.chdir("R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
                # Process POTCAR and INCAR
                poscar2potcar()
                os.system("vasp_process.py POSCAR")
                # Modify INCAR for magnetic calculations
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
                # Convert IStructure to Structure
                struc = structure.Structure(struc.lattice, struc.species, struc.cart_coords, coords_are_cartesian=True)
                obj.structure = struc
                comp_list = []
                for composition in struc.composition.elements:
                    comp_list.append(composition.name)
                obj.comp_list = comp_list
                obj.getkpt(primitive=False)
                evenkpt = input_data['download']['inp']['evenkpt']
                if evenkpt:
                    print("Utilizing even kpoint mesh\n")
                    obj.getevenkpt()
                obj.maxecut_sssp_for_subs()
                obj.prefix = struc.composition.alphabetical_formula.replace(" ","")
                obj.setting_qeinput(magnetic=True,monoclinic=False,pseudo_dir='../../pp')
                if not os.path.isdir("scf_dir"):
                    os.system("mkdir scf_dir")
                os.system("""mv scf-None.in""" + """ scf_dir/scf-{}-{}.in""".format(mpid,i+1))
            with open("mpid-magnetic.in", "a") as mpid_append:
                mpid_append.write("v{}".format(entry+1+i) + " " + mpid + "-{}".format(i+1) + " " + obj.prefix + "\n")
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
    # Define magnetic axis configurations
    saxis = input_data['magmom']['saxis']
    # Process each magnetic axis configuration
    for i,axis in enumerate(saxis):
        sx = int(axis[0])
        sy = int(axis[1])
        sz = int(axis[2])
         # Create directory for the current magnetic axis configuration
        if not os.path.isdir("R{}-saxis-{}{}{}-{}".format(mpid,sx,sy,sz,comp)):
            os.mkdir("R{}-saxis-{}{}{}-{}".format(mpid,sx,sy,sz,comp))
        if not os.path.isdir("R{}-saxis-{}{}{}-{}/relax".format(mpid,sx,sy,sz,comp)):
            os.mkdir("R{}-saxis-{}{}{}-{}/relax".format(mpid,sx,sy,sz,comp))
        os.chdir("R{}-{}/relax/".format(mpid,comp))
        # Copy necessary input files
        os.system("cp INCAR POTCAR POSCAR KPOINTS ../../R{}-saxis-{}{}{}-{}/relax/".format(mpid,sx,sy,sz,comp))
        if os.path.isfile("CHGCAR"):
            os.system("cp CHGCAR ../../R{}-saxis-{}{}{}-{}/relax/".format(mpid,sx,sy,sz,comp))
        os.chdir("../../")
        # Uncomment following if calculations read from previous CHGCAR
        #os.system(f"""echo 'ICHARG = 11' >> R{}-saxis-{}{}{}-{}/relax/INCAR""".format(mpid,sx,sy,sz,comp))
        # Update INCAR with the new magnetic axis
        os.system("""sed -i '/SAXIS/d' R{}-saxis-{}{}{}-{}/relax/INCAR""".format(mpid,sx,sy,sz,comp))
        os.system("""echo 'SAXIS = {} {} {}' >> R{}-saxis-{}{}{}-{}/relax/INCAR""".format(sx,sy,sz,mpid,sx,sy,sz,comp))
        # Append to the mpid-magnetic.in file
        with open("mpid-magnetic.in", "a") as mpid_append:
            mpid_append.write("v{}".format(entry+1+i) + " " + mpid + "-saxis-{}{}{}".format(sx,sy,sz) + " " + comp + "\n")
        if os.path.isfile("config.json"):
            os.system("""cp config.json R{}-saxis-{}{}{}-{}/relax/""".format(mpid,sx,sy,sz,comp))
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
    # DFT calculation type
    dft = input_data['download']['inp']['calc']
    # Ordering or anisotropy ?
    mag_type = input_data['magmom']['type']
    # Read input.in to obtain list of materials
    with open("input.in","r") as read_in:
        lines = read_in.readlines()
    start = int(lines[0].split("\n")[0])
    end = int(lines[1].split("\n")[0])
    filename = lines[3].split("\n")[0]
    with open(filename,'r') as read_mpid:
        lines = read_mpid.readlines()
    lines = lines[start-1:end-1]
    # Initiate MpConnect object
    obj = MpConnect()
    config = ['ferromagnetic','antiferromagnetic']
    # Loop over materials
    for line in lines:
        mpid = line.split(" ")[1]
        comp = line.split(" ")[2].split("\n")[0]
        # Creating different magnetic ordering
        if mag_type == "ordering":
            magnetic_structure(obj,mpid,comp,config,dft)
        # Creating input files with different SAXIS
        elif mag_type == "anisotropy":
            changeaxis(mpid,comp)
        else:
            print("Only ordering and anisotropy allowed\n")
if __name__ == "__main__":
    main()