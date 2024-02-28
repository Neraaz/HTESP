#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D."""
import sys
import os
import glob
import json
#import re
import warnings
import scipy.linalg as alg
#from ase import Atoms
from ase.io import vasp
from ase.data import atomic_masses,chemical_symbols
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import structure
from pymatgen.io import pwscf
import numpy as np
from write_potcar import poscar2potcar
from htepc import INPUTscf
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
#try:
#    with open("htepc.json", "r") as readjson:
#        input_data = json.load(readjson)
#except FileNotFoundError:
#    print("htepc.json file not found\n")

def pos_to_kpt(structure_filename,kpoint_density,evenkpt=False):
    """
    Function to obtain a k-point mesh from a structure file.

    Parameters:
    - structure_filename (str): Path to the structure file (QE scf.in or VASP POSCAR).
    - kpoint_density (float): Desired k-point density.
    - evenkpt (bool): Flag indicating whether to enforce an even number of k-points along each axis. Default is False.

    Returns:
    - kmesh (list): K-point mesh according to the k-point density.
    """
    kptsp = kpoint_density
    with open(structure_filename,"r") as read_struc:
        lines = read_struc.readlines()
    # Get cell vectors
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
    bnorm = np.zeros(3)
    bnorm = alg.norm(bmat,axis=1)
    kratio = [bnorm[i] / bnorm[0] for i in range(3)]
    klat = bnorm[0] / kptsp
    kmesh = [int(kratio[i] * klat + 0.5) if int(kratio[i] * klat + 0.5) != 0 else 1 for i in range(3)]
    kratio = [bnorm[i] / bnorm[0] for i in range(3)]
    klat = bnorm[0] / kptsp
    kmesh = [int(kratio[i] * klat + 0.5) if int(kratio[i] * klat + 0.5) != 0 else 1 for i in range(3)]
    if evenkpt:
        for i in range(3):
            if kmesh[i]%2 == 0:
                kmesh[i] = kmesh[i]
            else:
                kmesh[i] = kmesh[i] + 1
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
    """
    cif_parser = CifParser(infile)
    structure = cif_parser.get_structures()[0]  # Assuming there's only one structure in the CIF
    cif_writer = CifWriter(structure, symprec=0.1)
    cif_writer.write_file(infile)
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
    """
    pwd = os.getcwd()
    check_pwd = len(glob.glob("{}.cif".format(mpid)))
    try:
        if check_pwd > 0:
            file_path = glob.glob(pwd+"/{}.cif".format(mpid), recursive=True)[0]
    except FileNotFoundError:
        print("File {}.cif not found".format(mpid) + "\n")
    if not cif2cell:
        struc_poscar = pymatgen_cif(file_path)
        cell = struc_poscar.lattice.matrix
        compound = str(struc_poscar.composition).replace(" ", "")
        pos = struc_poscar.frac_coords
        symbol = []
        for site in struc_poscar.sites:
            symbol.append(str(site.specie))
        #with open(file_path, 'r') as file:
        #    lines = file.readlines()
        #cell = []
        #pos=[]
        #symbol = []
        #for line in lines:
        #    if 'length' in line or 'angle' in line:
        #        cell.append(float(line.split()[1]))
        #    #if '_chemical_formula_structural' in line:
        #    if "_chemical_formula_sum" in line:
        #        compound = line.split("\n")[0].split("  ")[1].replace(" ","").replace("'","")
        #    #    compound = line.split('\n')[0].split()[1]
        #index = lines.index(' _atom_site_occupancy\n')
        #lines = lines[index+1:]
        #for line in lines:
        #    sym = line.split()[0]
        #    sym = re.sub(r'[^A-Za-z\s]', '', sym)
        #    symbol.append(sym)
        #    pos.append([float(line.split()[3]), float(line.split()[4]),float(line.split()[5])])
        #crystal = Atoms(symbols=symbol, positions=pos, cell=cell,pbc=True)
        #cell = crystal.cell
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
    #if os.path.isfile("download.py"):
    #    import download as d
    evenkpt = input_data['download']['inp']['evenkpt']
    kptden = input_data['kptden']
    if evenkpt:
        k_mesh = pos_to_kpt("POSCAR",kptden,True)
    else:
        k_mesh = pos_to_kpt("POSCAR",kptden)
    if calc_type in ('VASP','vasp'):
        if not os.path.isdir('R{}-{}'.format(mpid,compound)):
            os.system("mkdir R{}-{}".format(mpid,compound))
        if not os.path.isdir('R{}-{}/relax/'.format(mpid,compound)):
            os.system("mkdir R{}-{}/relax/".format(mpid,compound))
        structure_file = structure.Structure.from_file("POSCAR")
        relax_set = MPRelaxSet(structure=structure_file)
        #relax_set.potcar.write_file("POTCAR")
        relax_set.poscar.write_file("POSCAR")
        poscar2potcar()
        relax_set.incar.write_file("INCAR")
        os.system("mv KPOINTS POSCAR INCAR POTCAR R{}-{}/relax/".format(mpid,compound))
        print(compound)
        if os.path.isfile("vasp.in"):
            os.system("cp vasp.in R{}-{}/relax".format(mpid,compound))
            pwd = os.getcwd()
            relax_folder = pwd + "/R{}-{}/relax".format(mpid,compound)
            os.chdir(relax_folder)
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
    if calc_type in ('QE','qe'):
        structure_file = structure.Structure.from_file("POSCAR")
        #structure_file = SpacegroupAnalyzer(structure_file, symprec=0.1).get_primitive_standard_structure()
        if os.path.isfile("htepc.json") or os.path.isfile("../../htepc.json"):
            psd_data = input_data['pseudo']
            dict_element = psd_data.PSEUDO
        else:
            print("htepc.json file not found\n")
            print("Creating one with default values from SSSP efficiency set\n")
            print("https://www.materialscloud.org/discover/sssp/table/efficiency\n")
            dict_element = {'H': 60, 'Li': 40, 'Be': 40, 'N': 60, 'F': 45, 'Na': 40, 'Mg': 30, 'Al': 30, 'Si': 30, 'P': 30, 'S': 35, 'Cl': 40, 'K': 60, 'Ca': 30, 'Sc': 40, 'Ti': 35, 'V': 35, 'Cr': 40, 'Mn': 65, 'Fe': 90, 'Co': 45, 'Ni': 45, 'Cu': 55, 'Zn': 40, 'Ga': 70, 'Ge': 40, 'As': 35, 'Br': 30, 'Rb': 30, 'Sr': 30, 'Y': 35, 'Zr': 30, 'Nb': 40, 'Mo': 35, 'Tc': 30, 'Ru': 35, 'Rh': 35, 'Pd': 45, 'Ag': 50, 'Cd': 60, 'In': 50, 'Sn': 60, 'Sb': 40, 'Te': 30, 'I': 35, 'Cs': 30, 'Ba': 30, 'La': 40, 'Hf': 50, 'Ta': 45, 'W': 30, 'Re': 30, 'Os': 40, 'Ir': 55, 'Pt': 35, 'Hg': 50, 'Tl': 50, 'Pb': 40, 'Bi': 45, 'B': 35, 'C': 45, 'Au': 45, 'Se': 30, 'O': 60}
            with open("pseudo.txt", "w") as write_pseudo:
                write_pseudo.write("PSEUDO={}".format(dict_element))
        #pseudo = pd.read_csv('pseudo.csv')
        #print(pseudo)
        #element_list = pseudo['Element'].tolist()
        #SSSP_cutoff_list = pseudo['SSSPcutoff'].tolist()
        #element_list = list(pseudo.keys())
        #SSSP_cutoff_list = list(pseudo.values())
        #nlist = len(element_list)
        #dict_element = {}
        #for i in range(nlist):
        #    dict_element[element_list[i]] = SSSP_cutoff_list[i]
        #
        dict_mass = {}
        nel = len(atomic_masses)
        for i in range(nel):
            dict_mass[chemical_symbols[i]] = atomic_masses[i]
        cutoff_list = [dict_element[el] for el in symbol]
        cutoff = max(cutoff_list)
        rhocutoff = 8 * cutoff
        prefix = compound
        pseudo1 = {el:el+'.upf' for el in set(symbol)}
        if check_pwd > 0:
            if os.path.isfile("htepc.json") or os.path.isfile("../../htepc.json"):
                pwscf_in = input_data['pwscf_in']
                control = pwscf_in['control']
                system = pwscf_in['system']
                electrons = pwscf_in['electrons']
            else:
                control = {'calculation':'vc-relax', 'nstep':300, 'restart_mode':'from_scratch', 'pseudo_dir':'../../pp', 'outdir':'./', 'tprnfor':'.true.','tstress':'.true.', 'etot_conv_thr':0.00001, 'forc_conv_thr':0.0001}
                system = {'smearing':'gauss', 'occupations':'smearing', 'degauss':0.02}
                electrons = {'diagonalization':'david', 'mixing_mode':'plain', 'mixing_beta':0.7, 'conv_thr': 1e-16, 'electron_maxstep':300}
                with open("pwscf_in.txt", "w") as pw_in:
                    pw_in.write("control={}".format(control) + "\n")
                    pw_in.write("system={}".format(system) + "\n")
                    pw_in.write("electrons={}".format(electrons) + "\n")
            system['ecutwfc'] = cutoff
            system['ecutrho'] = rhocutoff
            control['prefix'] = prefix
            filename = pwscf.PWInput(structure_file, pseudo=pseudo1, control=control, system=system,electrons=electrons, kpoints_grid=tuple(k_mesh), ions={'ion_dynamics':'bfgs'}, cell={'cell_dynamics':'bfgs','press_conv_thr':0.05})
            filename.write_file("temp.in")
            os.system("""sed "s/'.true.'/.true./" {} > {}""".format("temp.in","scf-{}.in".format(mpid)))
            os.system("sed -n '/K_POINTS automatic/,/CELL_PARAMETERS angstrom/p' scf-{}.in | sed '$d' > kpoint-{}.dat".format(mpid,mpid))
            obj = INPUTscf("scf-{}.in".format(mpid))
            obj.standardize(mpid,output="temp.dat")
            os.system("sed -n '/&CONTROL/,/ATOMIC_SPECIES/p' scf-{}.in | sed '$d' > scf.header".format(mpid))
            os.system("sed -n '/ATOMIC_SPECIES/,/ATOMIC_POSITIONS crystal/p' scf-{}.in | sed '$d' > species".format(mpid))
            os.system("rm temp.in")
            os.system("cat scf.header species temp.dat > scf-{}.in".format(mpid))
            os.system("rm POSCAR scf.header species temp.dat kpoint*")
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
    list_cif = glob.glob("*.cif",recursive=True)
    k_ind = 1
    calc_type = sys.argv[1]
    # CIF2CELL if True, uses cif2cell package to create POSCAR from given .cif files.
    # if False, It uses pymatgen cifparser to read and produce cif output, which then explicitely
    # read and uses to write POSCAR.
    cif2cell = input_data['download']['inp']['use_cif2cell']
    if cif2cell:
        print("CIF2CELL is True. Using cif2cell package....\n")
    for cif in list_cif:
        mpid = cif.split('.')[0]
        compound = ciftoscf(calc_type,mpid,cif2cell,False)
        if not os.path.isfile('mpid.in'):
            k_ind = 0
            with open("mpid.in", "w") as write_mpid:
                write_mpid.write("v{}".format(k_ind+1) + " " + mpid + " " + compound + "\n")
        else:
            with open('mpid.in', 'r') as read_mpid:
                lines = read_mpid.readlines()
            k_ind = len(lines)
            if not any(mpid in line for line in lines):
                with open("mpid.in", "a") as write_mpid:
                    write_mpid.write("v{}".format(k_ind+1) + " " + mpid + " " + compound + "\n")
if __name__ == '__main__':
    main()
