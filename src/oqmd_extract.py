#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to extract data from OQMD database"""
import sys
import os
import re
from collections import Counter
from ase.io.vasp import read_vasp
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.core import structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import qmpy_rester as qr
from cif_to_gsinput import pos_to_kpt
from write_potcar import poscar2potcar
from htepc import MpConnect
from check_json import config

def poscar_to_input(calc_type,mpid,compound,keven):
    """
    Function to convert POSCAR into ground-state input files.

    Parameters
    ----------
    calc_type : str
        Calculation type. "QE" or "VASP".
    mpid : str
        Material id.
    compound : str
        Compound name.
    keven : bool
        If even k-mesh to use.

    Returns
    -------
    str
        Compound name.

    Notes
    -----
    This function reads a POSCAR file and generates input files
    required for ground-state calculations
    using either Quantum Espresso (QE) or VASP software.

    It checks if 'config.json' exists, and if so, it retrieves settings such as k-point density and
    download configurations.

    Based on the calculation type, it creates input files and directories accordingly. For VASP,
    it prepares INCAR, POSCAR, and POTCAR files and
    potentially runs VASP processing scripts if available.
    For QE, it creates input files with appropriate settings.

    Finally, it moves input files to the designated directories.

    """
    input_data = config()
    if os.path.isfile("config.json"):
        d = input_data['download']
    evenkpt = d['inp']['evenkpt']
    kptden = input_data['kptden']
    # Check if even kpoint mesh is required
    # Calculate k-mesh grid from POSCAR
    if evenkpt:
        print("Even kpoint mesh is utilized\n")
        k_mesh = pos_to_kpt("POSCAR",kptden,True)
    else:
        k_mesh = pos_to_kpt("POSCAR",kptden)
    # Read POSCAR with ASE
    data = read_vasp("POSCAR")
    symbol = list(data.symbols)
    # Creates input for VASP
    if calc_type in ('VASP','vasp'):
        if not os.path.isdir('R{}-{}'.format(mpid,compound)):
            os.system("mkdir R{}-{}".format(mpid,compound))
        if not os.path.isdir('R{}-{}/relax/'.format(mpid,compound)):
            os.system("mkdir R{}-{}/relax/".format(mpid,compound))
        # Write VASP input files
        structure_file = structure.Structure.from_file("POSCAR")
        structure_file = SpacegroupAnalyzer(structure_file, symprec=0.1).get_primitive_standard_structure()
        relax_set = MPRelaxSet(structure=structure_file)
        relax_set.poscar.write_file("POSCAR")
        poscar2potcar()
        relax_set.incar.write_file("INCAR")
        # Copy files to R{mpid}-{compound}/relax/ folder
        os.system("mv KPOINTS POSCAR INCAR POTCAR R{}-{}/relax/".format(mpid,compound))
        # Update INCAR with vasp.in
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
def search(kwargs,properties):
    """
    Function to search compounds.

    Parameters
    ----------
    kwargs : dict
        Dictionary object for query to the OQMD website.
    properties : list
        Properties to extract.

    Notes
    -----
    This function uses QMPYRester to retrieve data from the OQMD (Open Quantum Materials Database) website.
    It attempts to get data with the provided query parameters. If the limit is reached, it decreases the limit
    and retries until successful.

    It then writes the retrieved data to 'download.csv' file,
    including IDs and properties specified in the list
    'properties'. It also creates 'mpid-list.in' file containing
    IDs, names, and compounds for each entry.
    Example
    -------
    >>> kwargs = {'element_set': '(Fe-Mn),B', 'stability': '<-0.1', 'natom': '<10', 'limit': 100}
    >>> properties = ['entry_id', 'name', 'enthalpy', 'composition']
    >>> search(kwargs, properties)
    """
    # Initialize QMPYRester object
    obj = qr.QMPYRester()
    limit = kwargs['limit']
    orig_limit = limit
    # Loop to retry with a decreased limit if the initial attempt fails
    while limit > 0:
        try:
            # Attempt to retrieve data from OQMD with provided query parameters
            data = obj.get_oqmd_phases(**kwargs)
            # Extract data from response
            data = data['data']
            # Exit loop if data retrieval is successful
            break
        except:
            # If retrieval fails, decrease the limit and retry
            limit -= int(orig_limit*0.1)
            kwargs['limit'] = limit
            continue
    # Counter for entries
    ind = 1
    # Create 'download' directory if it doesn't exist
    if not os.path.isdir("download"):
        os.mkdir("download")
    # Write data to CSV file and create 'mpid-list.in' file
    with open("download/download-oqmd.csv", "w") as write_download:
        write_download.write("ID,")
        for i,prop in enumerate(properties):
            if prop == 'composition':
                write_download.write("compound,")
            else:
                if i < len(properties) - 1:
                    write_download.write(prop + ",")
                else:
                    write_download.write(prop)
        write_download.write("\n")
    with open("mpid-list.in", "w") as write_mpid:
        for subdata in data:
            strings = subdata['name']
            strings = re.findall(r'[A-Z][a-z]*', strings)
            elements = kwargs['element_set']
            elements = re.findall(r'[A-Za-z]+', elements)
            exists_in_string = any(string not in elements for string in strings)
            # Write entry to 'mpid-list.in' and CSV file if it meets the criteria
            if not exists_in_string:
                write_mpid.write("v{}".format(ind) + " " + "oqmd-"+str(subdata['entry_id']) + " " + subdata['name'] + "\n")
                ind += 1
                with open("download/download-oqmd.csv", "a") as write_download:
                    write_download.write("oqmd-" + str(subdata['entry_id']) + ",")
                    for i,prop in enumerate(properties):
                        if prop == 'composition':
                            write_download.write(subdata[prop].replace(" ", "") + ",")
                        else:
                            if i < len(properties) - 1:
                                write_download.write(str(subdata[prop]) + ",")
                            else:
                                write_download.write(str(subdata[prop]))
                    write_download.write("\n")
def download(calc_type,start,end):
    """
    Function to create input files.

    Parameters
    ----------
    calc_type : str
        Type of calculations, QE or VASP.
    start : int
        Start index of the compounds to download.
    end : int
        End index of the compounds to download.

    Notes
    -----
    This function utilizes QMPYRester to retrieve data from the Open Quantum Materials Database (OQMD).
    It reads material IDs and compounds from 'mpid-list.in' file and downloads the corresponding data.
    For each compound, it creates a POSCAR file containing the compound structure.
    Then, it invokes the 'poscar_to_input' function to generate input files for QE or VASP calculations.
    Finally, it updates 'mpid.in' file with downloaded compounds.

    """
    #kwargs = {
    #    'element_set': '(Fe-Mn),B',      # composition include (Fe OR Mn) AND O
    #    'stability': '<-0.1',            # hull distance smaller than -0.1 eV
    #    'natom': '<10',                  # number of atoms less than 10
    #    }
    obj = qr.QMPYRester()
    # Read the material IDs and compounds from 'mpid-list.in' file
    with open("mpid-list.in","r") as read_mpid:
        mpid_data = read_mpid.readlines()
    mpid_data = mpid_data[start-1:end-1]
    # Loop through the specified range of compounds
    for mpid in mpid_data:
        oqmd_id = int(mpid.split(" ")[1].split("-")[1])
        subd = obj.get_entry_by_id(oqmd_id)
        try:
            # Extract necessary data from the OQMD entry
            oqmd_id = subd['id']
            oqmd_id = "oqmd-" + str(oqmd_id)
            compound = subd['name']
            # Write POSCAR file with the compound structure
            with open('POSCAR', 'w') as write_poscar:
                write_poscar.write(compound + '\n')
                write_poscar.write("1.0 \n")
                unit_cell = subd['unit_cell']
                for lattice in unit_cell:
                    write_poscar.write(str(lattice[0]) + " " + str(lattice[1]) + " " + str(lattice[2]) + "\n")
                basis = []
                elements = []
                for i,elm in enumerate(subd['sites']):
                    elements.append(subd['sites'][i].split("@")[0].split(" ")[0])
                    basis.append(subd['sites'][i].split("@")[1].split(" ")[1:])
                elem_count = Counter(elements)
                for elm in elem_count:
                    write_poscar.write(elm + " ")
                write_poscar.write("\n")
                for elm in elem_count:
                    write_poscar.write(str(elem_count[elm]) + " ")
                write_poscar.write("\n")
                write_poscar.write("Direct\n")
                for i,bas in enumerate(basis):
                    write_poscar.write(bas[0] + " " + bas[1] + " " + bas[2] + " " + elements[i] + "\n")
            # Generate input files for QE or VASP calculations
            poscar_to_input(calc_type,oqmd_id,compound,False)
            if os.path.isfile("POSCAR"):
                os.system("mv POSCAR POSCAR-{}".format(oqmd_id))
            # Update 'mpid.in' file with downloaded compounds
            if not os.path.isfile('mpid.in'):
                k_ind = 0
                with open("mpid.in", "w") as write_mpid:
                    write_mpid.write("v{}".format(k_ind+1) + " " + oqmd_id + " " + compound + "\n")
            else:
                with open('mpid.in', 'r') as read_mpid:
                    lines = read_mpid.readlines()
                k_ind = len(lines)
                if not any(oqmd_id in line for line in lines):
                    with open("mpid.in", "a") as write_mpid:
                        write_mpid.write("v{}".format(k_ind+1) + " " + oqmd_id + " " + compound + "\n")
        except:
            print(f"Structure data not found for {oqmd_id}-{subd['name']}\n")
            continue
if __name__ == "__main__":
    input_data = config()
    CONDITION = sys.argv[1]
    # Read config.json file and creates KWARGS dictionary
    with open("input.in","r") as read_in:
        lines = read_in.readlines()
    START = int(lines[0].split("\n")[0])
    END = int(lines[1].split("\n")[0])
    if os.path.isfile("config.json"):
        d = input_data['download']
        #START = d['inp']['start']
        #END = d['inp']['end']
        doqmd = d['oqmd']
        LIMIT = doqmd['limit']
        NTYPE = doqmd['ntype_constraint']
        ELM_LIST = doqmd['entries']
        NELM = len(ELM_LIST)
        if NELM == 1:
            ELMS = ELM_LIST[0]
        else:
            ELMS = ''
            for i,ELM in enumerate(ELM_LIST):
                if i < NELM - 1:
                    ELMS = ELMS + ELM + "-"
                else:
                    ELMS = ELMS + ELM
        ELMS = '({})'.format(ELMS)
        MUST_INCLUDE = doqmd['must_include']
        if len(MUST_INCLUDE) > 0:
            ELMS = ELMS + ","
        for j, ELM in enumerate(MUST_INCLUDE):
            if j < len(MUST_INCLUDE) - 1 and j > 0:
                ELMS += ELM + ","
            else:
                ELMS += ELM
        PROPERTIES = doqmd['prop']
        #PROPERTIES = ['composition','spacegroup','volume','band_gap','stability']
        METAL=doqmd['metal']
        if METAL:
            BAND_GAP = 0.001
        else:
            BAND_GAP = None
        NEG_FE=doqmd['FE']
        if NEG_FE:
            FE = 0.0001
        else:
            FE = None
        THERMO_STABLE=doqmd['thermo_stable']
        if THERMO_STABLE:
            EBH = 0.001
        else:
            EBH = None
        NSITES=doqmd['size_constraint']
        SPACEGROUP=doqmd['spacegroup']
        KWARGS = {
                   'element_set': '{}'.format(ELMS),
                   'band_gap': '<{}'.format(BAND_GAP),
                   'delta_e': '<{}'.format(FE),
                   'stability': '<{}'.format(EBH),
                   'natoms': '<{}'.format(NSITES),
                   'ntypes': '<{}'.format(NTYPE),
                   'limit' : LIMIT,
                  }
        if not METAL:
            KWARGS.pop('band_gap')
        if not NEG_FE:
            KWARGS.pop('delta_e')
        if not THERMO_STABLE:
            KWARGS.pop('stability')
        CALC_TYPE = d['inp']['calc']
        print(CALC_TYPE)
    else:
        print("input file config.json not found\n")
        print("Utilizing default settings\n")
        NTYPE = '<5' #Number of different types of element in the compound.
        NSITES = '<10'
        ELM_LIST = ['B','C']
        EXCLUDE_EL = ['Lu']
        METAL = False
        NEG_FE = False
        THERMO_STABLE = False
        ORDERING = None
        SPACEGROUP = None
        CALC_TYPE = 'QE'
        PROPERTIES = ['composition','spacegroup','volume','band_gap','stability']
        KWARGS = {
                   'element_set': '(B-C)',
                   'band_gap': '<2.0',
                   'delta_e': '<0.0001',
                   'stability': '<-0.1',
                   'natoms': NSITES,
                   'ntypes': NTYPE,
                   'limit' : LIMIT,
                  }
    if CONDITION == 'search':
        search(KWARGS,PROPERTIES)
    elif CONDITION == 'download':
        download(CALC_TYPE,START,END)
    else:
        print("wrong options. Allowed options are either search or download\n")
