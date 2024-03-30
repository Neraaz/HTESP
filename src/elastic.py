#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D. (tug11655@temple.edu)"""
"""Module to perform elastic calculations"""
import sys
import os
import json
import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io import pwscf
from pymatgen.io.vasp.sets import Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.magnetism import MagneticStructureEnumerator
from pymatgen.core import structure
from pymatgen.analysis.elasticity import diff_fit,ElasticTensor,Stress
from pymatgen.analysis.elasticity import DeformedStructureSet,find_eq_stress
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

def deformation(mpid,obj,dft,orig_prefix,deformed_struc):
    """
    Function to create deformed structures utilizing pymatgen.analysis.elastic class.

    Parameters:
    -----------
    mpid : str
        Materials id.

    obj : object
        Object of MpConnect class.

    dft : str
        Density Functional Theory (DFT) method used, e.g., 'vasp', 'qe'.

    orig_prefix : str
        Prefix for the original undeformed structure.

    deformed_struc : object
        Object containing deformed structures.

    Returns:
    --------
    None

    Example:
    --------
    >>> from mpconnect import MpConnect
    >>> from pymatgen import Structure
    >>> from pymatgen.io.vasp import Poscar
    >>> obj = MpConnect()
    >>> mpid = "mp-1234"
    >>> orig_prefix = "Si2"
    >>> deformed_struc = ...  # Object containing deformed structures
    >>> # obtained with pymatgen.analysis.elasticity.DeformedStructureSet
    >>> deformation(mpid, obj, "vasp", orig_prefix, deformed_struc)
    """
    # Extracting deformed structures
    list_str = deformed_struc.deformed_structures
    if not os.path.isdir('scf_dir'):
        os.mkdir("scf_dir")
    # Initialize index when mpid-deformed.in file not found
    if not os.path.isfile('mpid-deformed.in'):
        entry = 0
    else:
        with open('mpid-deformed.in', 'r') as mpid_read:
            lines = mpid_read.readlines()
        entry = len(lines)
    # Loop over deformed structures
    for i,struc in enumerate(list_str):
        if dft in ('vasp', 'VASP'):
            # Write VASP input files inside R{mpid}-{i+1}-{name}/relax
            obj.prefix = struc.composition.alphabetical_formula.replace(" ","")
            poscar = Poscar(structure=struc,comment=obj.prefix)
            if not os.path.isdir("R{}-{}-{}".format(mpid,i+1,obj.prefix)):
                os.mkdir("R{}-{}-{}".format(mpid,i+1,obj.prefix))
            if not os.path.isdir("R{}-{}-{}/relax".format(mpid,i+1,obj.prefix)):
                os.mkdir("R{}-{}-{}/relax".format(mpid,i+1,obj.prefix))
            pwd=os.getcwd()
            os.system("cp R{}-{}/relax/INCAR {}/R{}-{}-{}/relax/".format(mpid,orig_prefix,pwd,mpid,i+1,obj.prefix))
            poscar.write_file(filename="R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix))
            #os.system("cp R{}-{}/relax/KPOINTS {}/R{}-{}-{}/relax/".format(mpid,orig_prefix,pwd,mpid,i+1,obj.prefix))
            if os.path.isfile("config.json") or os.path.isfile("../../config.json"):
                dwn = input_data['download']
            evenkpt = dwn['inp']['evenkpt']
            kptden = input_data['kptden']
            if evenkpt:
                print("Even kpoint mesh is utilized\n")
                pos_to_kpt("R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix),kptden,True)
            else:
                pos_to_kpt("R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix),kptden)
            os.system("mv KPOINTS R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
            structure_file = structure.Structure.from_file("R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix))
            relax_set = MPRelaxSet(structure=structure_file)
            #relax_set.potcar.write_file("R{}-{}-{}/relax/POTCAR".format(mpid,i+1,obj.prefix))
            if os.path.isfile("vasp.in"):
                os.system("cp vasp.in R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
            if os.path.isfile("config.json"):
                os.system("cp config.json R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
            os.chdir("R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
            # Write POTCAR from POSCAR
            poscar2potcar()
            # process INCAR file with vasp.in file
            os.system("vasp_process.py POSCAR")
            os.system("""sed -i '/ISIF/d' INCAR""")
            os.system("""echo "ISIF = 2" >> INCAR""")
            os.chdir("../../")
        else:
            # Read magnetic flag from pwscf_in dictionary
            magnetic = input_data['pwscf_in']['magnetic']
            # Create input files with FM ordering if magnetic flag is true
            if magnetic:
                default_magmoms = input_data['magmom']['magmom']
                struc = MagneticStructureEnumerator(struc,default_magmoms=default_magmoms,strategies=['ferromagnetic'],truncate_by_symmetry=True).ordered_structures
                obj.structure = struc[0]
            else:
                obj.structure = struc
            # Create scf.in file
            comp_list = []
            for composition in obj.structure.composition.elements:
                comp_list.append(composition.name)
            obj.comp_list = comp_list
            obj.getkpt()
            evenkpt = input_data['download']['inp']['evenkpt']
            if evenkpt:
                print("Utilizing even kpoint mesh\n")
                obj.getevenkpt()
            obj.maxecut_sssp_for_subs()
            obj.prefix = struc.composition.alphabetical_formula.replace(" ","")
            if magnetic:
                obj.setting_qeinput(calculation='relax',magnetic=True,pseudo_dir='../../pp')
            else:
                obj.setting_qeinput(calculation='relax',pseudo_dir='../../pp')
            os.system("""mv scf-None.in""" + """ scf_dir/scf-{}-{}.in""".format(mpid,i+1))
        print(mpid,obj.prefix)
        with open("mpid-deformed.in", "a") as mpid_append:
            mpid_append.write("v{}".format(entry+1+i) + " " + mpid + "-{}".format(i+1) + " " + obj.prefix + "\n")
def main(mpid,orig_prefix):
    """
    Main function to handle deformation and elastic computations.

    Parameters:
    -----------
    mpid : str
        Materials id.

    orig_prefix : str
        Prefix for the undeformed structure.

    Returns:
    --------
    None
    """
    mode = sys.argv[1]
    dft = input_data['download']['inp']['calc']
    obj = MpConnect()
    # Determine strain based on file existence
    if os.path.isfile("config.json") or os.path.isfile("../../config.json"):
        strain = input_data['strain']
    else:
        strain = [-0.01,-0.005,0.005,0.01]
    try:
        # Try loading the structure from relaxed POSCAR
        obj.structure = structure.Structure.from_file("R{}-{}/relax/POSCAR".format(mpid,orig_prefix))
    except:
        # Load structure from SCF input file if relaxed POSCAR doesn't exist
        obj.structure = pwscf.PWInput.from_file("scf_dir/scf-{}.in".format(mpid)).structure
    # Get conventional standardized structure
    structure_sym = SpacegroupAnalyzer(obj.structure, symprec=0.1).get_conventional_standard_structure()
    prefix_conv = structure_sym.composition.alphabetical_formula.replace(" ","")
    # Generate deformed structure set
    deformed_struc = DeformedStructureSet(structure_sym,norm_strains=strain)
    # Handle different modes of operation
    if mode == "input":
        # Prepare inputs and submit calculations
        deformation(mpid,obj,dft,orig_prefix,deformed_struc)
    elif mode == "compute_elastic":
        # Compute elastic constants
        nstruc = 6*len(strain)
        deformation_mat = deformed_struc.deformations
        deform_list = []
        stress = []
        #dataeq = Vasprun("R{}-{}/relax/vasprun.xml".format(mpid,orig_prefix))
        #stresseq = -1.0*np.array(dataeq.as_dict()['output']['ionic_steps'][0]['stress'])
        #print(stresseq)
        # Loop over inputs and extract stress from deformed structures after relaxation
        for istruc in range(nstruc):
            if dft in ('VASP','vasp'):
                data = Vasprun("R{}-{}-{}/relax/vasprun.xml".format(mpid,istruc+1,prefix_conv))
                calc_stress = -1.0*np.array(data.as_dict()['output']['ionic_steps'][0]['stress'])
            else:
                os.system(f"""grep -A 3 'total   stress' R{mpid}-{istruc+1}-{prefix_conv}/relax/scf.out | tail -n 3 > stress-qe""")
                data = np.loadtxt("stress-qe")[:,:3]
                calc_stress = data*21798.7 #1 Ry/angstrom^3 to 1 kbar conversion
                #print("Currently supporting only for VASP.\n")
            # Creates stress object
            stress_obj = Stress(calc_stress)
            pk2stress = stress_obj.piola_kirchoff_2(deformation_mat[istruc])
            #pk2stress = calc_stress
            stress.append(pk2stress)
            # Extracts green lagrange strain
            deform_list.append(deformation_mat[istruc].green_lagrange_strain)
        stress = np.array(stress)
        strain = np.array(deform_list)
        # Fitting to obtain elastic tensor
        elastic_tens = diff_fit(strain,stress,order=2)[0]
        elastic_obj = ElasticTensor(elastic_tens)
        #print(ELASTIC_OBJ.property_dict)
        properties = elastic_obj.get_structure_property_dict(structure_sym)
        propname = properties.keys()
        print(propname)
        # Save elastic properties in elastic.csv file
        with open("elastic.csv","a") as write_elastic:
            write_elastic.write("{},{}".format(mpid,orig_prefix))
            for prop in properties:
                if prop != "structure":
                    write_elastic.write("," + str(round(properties[prop],2)))
            write_elastic.write("\n")
    else:
        print("Only 2 mode is available, either input or compute_elastic\n")
if __name__ == "__main__":
    with open("input.in","r") as read_in:
        lines1 = read_in.readlines()
    START = int(lines1[0].split("\n")[0])
    END = int(lines1[1].split("\n")[0])
    FILENAME = lines1[3].split("\n")[0]
    with open(FILENAME,'r') as read_mpid:
        LINES = read_mpid.readlines()
    LINES = LINES[START-1:END-1]
    PROPNAME = ['trans_v', 'long_v', 'snyder_ac', 'snyder_opt', 'snyder_total', 'clarke_thermalcond', 'cahill_thermalcond', 'debye_temperature', 'structure', 'k_voigt', 'k_reuss', 'k_vrh', 'g_voigt', 'g_reuss', 'g_vrh', 'universal_anisotropy', 'homogeneous_poisson', 'y_mod']
    if not os.path.isfile("elastic.csv"):
        with open("elastic.csv","w") as write_elastic1:
            write_elastic1.write("materials_id,compound")
            for prop1 in PROPNAME:
                write_elastic1.write("," + prop1)
            write_elastic1.write("\n")
    for LINE in LINES:
        MPID = LINE.split(" ")[1]
        COMP = LINE.split(" ")[2].split("\n")[0]
        main(MPID,COMP)
