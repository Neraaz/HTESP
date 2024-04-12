#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to write QE/VASP input files from .cif files"""
import sys
import os
import glob
import json
import warnings
import scipy.linalg as alg
from ase.io import vasp
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import structure
import numpy as np
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
def pos_to_kpt(structure_filename,kpoint_density,evenkpt=False):
    """
    Function to obtain a k-point mesh from a structure file.

    Parameters:
    - structure_filename (str): Path to the structure file (QE scf.in or VASP POSCAR).
    - kpoint_density (float): Desired k-point density.
    - evenkpt (bool): Flag indicating whether to enforce an even number of k-points along each axis.
      Default is False.

    Returns:
    - kmesh (list): K-point mesh according to the k-point density.
    Example:
        >>> # Generating a k-point mesh for a VASP POSCAR with a k-point density of 0.05
        >>> pos_to_kpt("POSCAR", 0.05, evenkpt=True)
    """
    kptsp = kpoint_density
    # Read the structure file
    with open(structure_filename,"r") as read_struc:
        lines = read_struc.readlines()
    # Get cell vectors and compute reciprocal lattice vectors
    line = lines[1].split()
    latscl = float(line[0])
    ain = np.zeros((3, 3))
    amat = np.zeros((3, 3))
    bmat = np.zeros((3, 3))
    anorm = np.zeros(3)
    for i in range(3):
        line = lines[2 + i].split()
        for j in range(3):
            ain[i][j] = float(line[j]) * latscl
            amat[j][i] = ain[i][j]
        anorm[i] = np.sqrt(ain[i][0] ** 2 + ain[i][1] ** 2 + ain[i][2] ** 2)
    bmat = alg.inv(amat)
    #bnorm = np.zeros(3)
    # Compute the length of reciprocal lattice vectors
    bnorm = alg.norm(bmat,axis=1)
    kratio = [bnorm[i] / bnorm[0] for i in range(3)]
    # Compute the k-point mesh
    klat = bnorm[0] / kptsp
    kmesh = [int(kratio[i] * klat + 0.5) if int(kratio[i] * klat + 0.5) != 0 else 1 for i in range(3)]
    # Enforce even number of kpoints if evenkpt flag is true
    if evenkpt:
        for i in range(3):
            if kmesh[i]%2 == 0:
                kmesh[i] = kmesh[i]
            else:
                kmesh[i] = kmesh[i] + 1
    # Write KPOINTS file
    with open("KPOINTS", "w") as write_kpt:
        write_kpt.write("KPOINTS automatic" + "\n")
        write_kpt.write("0"+"\n")
        write_kpt.write("Gamma\n")
        write_kpt.write(f"{kmesh[0]} {kmesh[1]} {kmesh[2]}"+"\n")
        write_kpt.write("0 0 0\n")
    return kmesh
def pymatgen_cif(infile):
    """
    Function to convert (ICSD) CIF files to pymatgen format.

    Parameters:
    - infile (str): Path to the input CIF file.

    Returns:
    - structure (Structure): Pymatgen Structure object representing the CIF structure.
    Example:
        >>> # Convert a CIF file named 'example.cif' to pymatgen format
        >>> pymatgen_cif("example.cif")
    """
    # Parse the CIF file using CifParser
    cif_parser = CifParser(infile)
    structure = cif_parser.get_structures()[0]  # Assuming there's only one structure in the CIF
    oxidation = False
    # Check the oxidation states and merge sites if present
    for elem in structure.elements:
        if '+' in str(elem):
            oxidation = True
    if oxidation:
        structure.merge_sites()
        structure = structure.remove_oxidation_states()
    # Write the CIF file using CifWriter
    cif_writer = CifWriter(structure, symprec=0.1)
    cif_writer.write_file(infile)
    # Get the primitive standard structure using SpacegroupAnalyzer
    structure = SpacegroupAnalyzer(structure=structure,symprec=0.1).get_primitive_standard_structure()
    return structure
def ciftoscf(calc_type,mpid,cif2cell=True,keven=False):
    """
    Function to convert a CIF file to a Quantum ESPRESSO (QE) input file 'scf.in'.

    Parameters:
    - calc_type (str): The type of calculation, either 'VASP' or 'QE'.
    - mpid (str): The Material ID.
    - cif2cell (bool): If True, uses cif2cell to process the CIF file. Default is True.
    - keven (bool): If True, even k-mesh is used. Default is False.

    Returns:
    - compound (str): The compound name.

    Writes the 'scf-mpid.in' file.
    Example:
        >>> # Convert a CIF file to a Quantum ESPRESSO input file
        >>> ciftoscf("QE", "mp-1234", cif2cell=False, keven=False)
    """
    # Get current directory
    pwd = os.getcwd()
    # Check if CIF file exists in the current directory
    check_pwd = len(glob.glob("{}.cif".format(mpid)))
    try:
        if check_pwd > 0:
            file_path = glob.glob(pwd+"/{}.cif".format(mpid), recursive=True)[0]
    except FileNotFoundError:
        print("File {}.cif not found".format(mpid) + "\n")
    # Process CIF file using cif2cell if use_cif2cell flag is true
    # Otherwise, it uses pymatgen
    if not cif2cell:
        # Use pymatgen to process CIF file and obtain structure
        struc_poscar = pymatgen_cif(file_path)
        cell = struc_poscar.lattice.matrix
        compound = str(struc_poscar.composition).replace(" ", "")
        pos = struc_poscar.frac_coords
        symbol = []
        for site in struc_poscar.sites:
            symbol.append(str(site.specie))
        # Write POSCAR file
        with open("POSCAR", "w") as write_poscar:
            write_poscar.write("Structure for {}".format(mpid) + "\n")
            write_poscar.write("1.0\n")
            for i in range(3):
                write_poscar.write("{} {} {}".format(cell[i][0], cell[i][1], cell[i][2]) + "\n")
            dict_symbol = {}
            for i,sym in enumerate(symbol):
                if sym not in dict_symbol:
                    dict_symbol[sym] = 1
                else:
                    dict_symbol[sym] += 1
            for symb in dict_symbol:
                write_poscar.write(symb + " ")
            write_poscar.write("\n")
            for symb in dict_symbol:
                write_poscar.write(str(dict_symbol[symb]) + " ")
            write_poscar.write("\n")
            write_poscar.write("Direct\n")
            for i,_ in enumerate(symbol):
                write_poscar.write(str(pos[i][0]) + " ")
                write_poscar.write(str(pos[i][1]) + " " + str(pos[i][2]) + " " + symbol[i] +"\n")
    else:
        # Run cif2cell to generate POSCAR
        os.system("""cif2cell {} -p vasp --vasp-cartesian-lattice-vectors""".format(file_path))
        with open("POSCAR", "r") as read_poscar:
            lines = read_poscar.readlines()
        list_first = lines[0].split("\n")[0].split(" ")
        for i,element in enumerate(list_first):
            if "order:" in element:
                index = i
        insert_text = ""
        ion_name = list_first[index+1:-1]
        for ion in ion_name:
            insert_text += ion + " "
        os.system("""sed '5 a {}' POSCAR > POSCAR_new""".format(insert_text))
        os.system("""mv POSCAR_new POSCAR""")
        symbol = ion_name
        data = vasp.read_vasp('POSCAR')
        compound = str(data.symbols)
    # Obtain k-point mesh
    evenkpt = input_data['download']['inp']['evenkpt']
    kptden = input_data['kptden']
    if evenkpt:
        k_mesh = pos_to_kpt("POSCAR",kptden,True)
    else:
        k_mesh = pos_to_kpt("POSCAR",kptden)
    # Generate VASP input set
    if calc_type in ('VASP','vasp'):
        if not os.path.isdir('R{}-{}'.format(mpid,compound)):
            os.system("mkdir R{}-{}".format(mpid,compound))
        if not os.path.isdir('R{}-{}/relax/'.format(mpid,compound)):
            os.system("mkdir R{}-{}/relax/".format(mpid,compound))
        structure_file = structure.Structure.from_file("POSCAR")
        # Generate MPRelaxSet object
        relax_set = MPRelaxSet(structure=structure_file)
        relax_set.poscar.write_file("POSCAR")
        # Generate POTCAR from POSCAR file
        poscar2potcar()
        relax_set.incar.write_file("INCAR")
        os.system("mv KPOINTS POSCAR INCAR POTCAR R{}-{}/relax/".format(mpid,compound))
        print(compound)
        # Copy vasp.in file if it exists
        if os.path.isfile("vasp.in"):
            os.system("cp vasp.in R{}-{}/relax".format(mpid,compound))
            pwd = os.getcwd()
            relax_folder = pwd + "/R{}-{}/relax".format(mpid,compound)
            os.chdir(relax_folder)
            # Process VASP input files reading vasp.in file if it exists
            # Otherwise, it will simply print INCAR from MPRelaxSet
            os.system("vasp_process.py POSCAR")
            os.chdir(pwd)
    else:
        os.system("rm KPOINTS")
    if keven:
        for i in range(3):
            if k_mesh[i]%2 == 0:
                k_mesh[i] = k_mesh[i]
            else:
                k_mesh[i] = k_mesh[i] + 1
    # QE calculation setup
    if calc_type in ('QE','qe'):
        obj = MpConnect()
        structure_file = structure.Structure.from_file("POSCAR")
        magnetic = input_data['pwscf_in']['magnetic']
        # Magnetic (FM) structure generation if magnetic flag is true
        if magnetic:
            # Reading magnetic moments from input file
            default_magmoms = input_data['magmom']['magmom']
            # Obtaining Ferromagnetic structure
            structure_file.add_spin_by_element(default_magmoms)
            obj.structure = structure_file
        else:
            obj.structure = structure_file
        # Get composition list
        comp_list = []
        for composition in obj.structure.composition.elements:
            comp_list.append(composition.name)
        obj.comp_list = comp_list
        # Getting k-point mesh
        obj.getkpt()
        evenkpt = input_data['download']['inp']['evenkpt']
        if evenkpt:
            print("Utilizing even kpoint mesh\n")
            obj.getevenkpt()
        # Getting maximum kinetic energy cutoff among elements
        obj.maxecut_sssp()
        obj.prefix = compound
        obj.mpid = mpid
        # Writing QE input scf.in file
        if magnetic:
            obj.setting_qeinput(magnetic=True,pseudo_dir='../../pp')
        else:
            obj.setting_qeinput(pseudo_dir='../../pp')
        if not os.path.isdir("scf_dir"):
            os.system("mkdir scf_dir")
        os.system("mv scf-{}.in scf_dir".format(mpid))
    return compound
def main():
    """
    The main function of the script.

    This function performs the following operations:
    - Retrieves a list of CIF files in the current directory.
    - Sets the calculation type.
    - Sets the flag for using cif2cell.
    - Iterates over CIF files, converts them to QE input files, and writes data to 'mpid.in'.

    Parameters:
    - calc_type (str): Calculation type (e.g., 'VASP' or 'QE').

    Returns: None
    """
    warnings.filterwarnings('ignore')
    # Get a list of CIF files in the current directory
    list_cif = glob.glob("*.cif",recursive=True)
    # Initialize index for 'mpid.in' entries
    #k_ind = 1
    # Read calculation type from command line argument
    calc_type = sys.argv[1]
    # CIF2CELL if True, uses cif2cell package to create POSCAR from given .cif files.
    # if False, It uses pymatgen cifparser to read and produce cif output, which then explicitely
    # read and uses to write POSCAR.
    cif2cell = input_data['download']['inp']['use_cif2cell']
    if cif2cell:
        print("CIF2CELL is True. Using cif2cell package....\n")
    # Iterate over CIF files
    for cif in list_cif:
        mpid = cif.split('.')[0] # Extract Material ID from file name
        print(mpid)
        # Convert CIF to QE input file
        compound = ciftoscf(calc_type,mpid,cif2cell,False)
        # If mpid.in does not exist, initialize index
        if not os.path.isfile('mpid.in'):
            k_ind = 0
            with open("mpid.in", "w") as write_mpid:
                write_mpid.write("v{}".format(k_ind+1) + " " + mpid + " " + compound + "\n")
        else:
            # If mpid.in exists, append data with new index
            with open('mpid.in', 'r') as read_mpid:
                lines = read_mpid.readlines()
            k_ind = len(lines)
            if not any(mpid in line for line in lines):
                with open("mpid.in", "a") as write_mpid:
                    write_mpid.write("v{}".format(k_ind+1) + " " + mpid + " " + compound + "\n")
if __name__ == '__main__':
    main()
