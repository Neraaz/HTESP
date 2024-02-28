#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D."""
import sys
import os
import json
import re
from collections import Counter
from ase.io.vasp import read_vasp
from ase.data import atomic_masses,chemical_symbols
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.core import structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io import pwscf
import qmpy_rester as qr
from cif_to_gsinput import pos_to_kpt
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
    This function reads a POSCAR file and generates input files required for ground-state calculations
    using either Quantum Espresso (QE) or VASP software.

    It checks if 'htepc.json' exists, and if so, it retrieves settings such as k-point density and
    download configurations.

    Based on the calculation type, it creates input files and directories accordingly. For VASP,
    it prepares INCAR, POSCAR, and POTCAR files and potentially runs VASP processing scripts if available.
    For QE, it creates input files with appropriate settings.

    Finally, it moves input files to the designated directories.

    """
    if os.path.isfile("htepc.json"):
        d = input_data['download']
    evenkpt = d['inp']['evenkpt']
    kptden = input_data['kptden']
    if evenkpt:
        print("Even kpoint mesh is utilized\n")
        k_mesh = pos_to_kpt("POSCAR",kptden,True)
    else:
        k_mesh = pos_to_kpt("POSCAR",kptden)
    data = read_vasp("POSCAR")
    symbol = list(data.symbols)
    #print(symbol)
    if calc_type in ('VASP','vasp'):
        if not os.path.isdir('R{}-{}'.format(mpid,compound)):
            os.system("mkdir R{}-{}".format(mpid,compound))
        if not os.path.isdir('R{}-{}/relax/'.format(mpid,compound)):
            os.system("mkdir R{}-{}/relax/".format(mpid,compound))
        structure_file = structure.Structure.from_file("POSCAR")
        structure_file = SpacegroupAnalyzer(structure_file, symprec=0.1).get_primitive_standard_structure()
        relax_set = MPRelaxSet(structure=structure_file)
        #relax_set.potcar.write_file("POTCAR")
        relax_set.poscar.write_file("POSCAR")
        poscar2potcar()
        relax_set.incar.write_file("INCAR")
        os.system("mv KPOINTS POSCAR INCAR POTCAR R{}-{}/relax/".format(mpid,compound))
        #os.system("mv KPOINTS POSCAR INCAR R{}-{}/relax/".format(mpid,compound))
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
        if os.path.isfile("htepc.json"):
            psd_data = input_data['pseudo']
            dict_element = psd_data['PSEUDO']
        else:
            print("pseudo not found\n")
            print("Creating one with default values from SSSP efficiency set\n")
            print("https://www.materialscloud.org/discover/sssp/table/efficiency\n")
            dict_element = {'H': 60, 'Li': 40, 'Be': 40, 'N': 60, 'F': 45, 'Na': 40, 'Mg': 30, 'Al': 30, 'Si': 30, 'P': 30, 'S': 35, 'Cl': 40, 'K': 60, 'Ca': 30, 'Sc': 40, 'Ti': 35, 'V': 35, 'Cr': 40, 'Mn': 65, 'Fe': 90, 'Co': 45, 'Ni': 45, 'Cu': 55, 'Zn': 40, 'Ga': 70, 'Ge': 40, 'As': 35, 'Br': 30, 'Rb': 30, 'Sr': 30, 'Y': 35, 'Zr': 30, 'Nb': 40, 'Mo': 35, 'Tc': 30, 'Ru': 35, 'Rh': 35, 'Pd': 45, 'Ag': 50, 'Cd': 60, 'In': 50, 'Sn': 60, 'Sb': 40, 'Te': 30, 'I': 35, 'Cs': 30, 'Ba': 30, 'La': 40, 'Hf': 50, 'Ta': 45, 'W': 30, 'Re': 30, 'Os': 40, 'Ir': 55, 'Pt': 35, 'Hg': 50, 'Tl': 50, 'Pb': 40, 'Bi': 45, 'B': 35, 'C': 45, 'Au': 45, 'Se': 30, 'O': 60}
            #with open("pseudo.py", "w") as write_pseudo:
            #    write_pseudo.write("PSEUDO={}".format(dict_element))
        dict_mass = {}
        nel = len(atomic_masses)
        for i in range(nel):
            dict_mass[chemical_symbols[i]] = atomic_masses[i]
        try:
            cutoff_list = [dict_element[el] for el in symbol]
        except:
            print("energy cutoff not provided for some of the elements. Please check htepc.json\n")
        cutoff = max(cutoff_list)
        rhocutoff = 8 * cutoff
        prefix = compound
        pseudo1 = {el:el+'.upf' for el in set(symbol)}
        check_pwd = 1
        if check_pwd > 0:
            try:
                pwscf_in = input_data['pwscf_in']
                control = pwscf_in['control']
                system = pwscf_in['system']
                electrons = pwscf_in['electrons']
            except ModuleNotFoundError:
                control = {'calculation':'vc-relax', 'nstep':300, 'restart_mode':'from_scratch', 'pseudo_dir':'../../pp', 'outdir':'./', 'tprnfor':'.true.','tstress':'.true.', 'etot_conv_thr':0.00001, 'forc_conv_thr':0.0001}
                system = {'smearing':'gauss', 'occupations':'smearing', 'degauss':0.02}
                electrons = {'diagonalization':'david', 'mixing_mode':'plain', 'mixing_beta':0.7, 'conv_thr': 1e-16, 'electron_maxstep':300}
                #with open("pwscf_in.py", "w") as pw_in:
                #    pw_in.write("control={}".format(control) + "\n")
                #    pw_in.write("system={}".format(system) + "\n")
                #    pw_in.write("electrons={}".format(electrons) + "\n")
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
    #return compound
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

    It then writes the retrieved data to 'download.csv' file, including IDs and properties specified in the list
    'properties'. It also creates 'mpid-list.in' file containing IDs, names, and compounds for each entry.
    """
    obj = qr.QMPYRester()
    limit = kwargs['limit']
    orig_limit = limit
    while limit > 0:
        try:
            data = obj.get_oqmd_phases(**kwargs)
            data = data['data']
            break
        except:
            limit -= int(orig_limit*0.1)
            kwargs['limit'] = limit
            continue
    ind = 1
    if not os.path.isdir("download"):
        os.mkdir("download")
    with open("download/download.csv", "w") as write_download:
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
            #elements = re.findall(r'[A-Za-z]+', elements)
            #print(strings,elements)
            exists_in_string = any(string not in elements for string in strings)
            if not exists_in_string:
                write_mpid.write("v{}".format(ind) + " " + "oqmd-"+str(subdata['entry_id']) + " " + subdata['name'] + "\n")
                ind += 1
                with open("download/download.csv", "a") as write_download:
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
    with open("mpid-list.in","r") as read_mpid:
        mpid_data = read_mpid.readlines()
    mpid_data = mpid_data[start-1:end-1]
    for mpid in mpid_data:
        #print(mpid)
        oqmd_id = int(mpid.split(" ")[1].split("-")[1])
        subd = obj.get_entry_by_id(oqmd_id)
    #data = obj.get_oqmd_phases(**kwargs)
    #for subd in data['data']:
        try:
            oqmd_id = subd['id']
            oqmd_id = "oqmd-" + str(oqmd_id)
            #print(oqmd_id + "\n")
            compound = subd['name']
            with open('POSCAR', 'w') as write_poscar:
                write_poscar.write(compound + '\n')
                write_poscar.write("1.0 \n")
                unit_cell = subd['unit_cell']
                for lattice in unit_cell:
                    write_poscar.write(str(lattice[0]) + " " + str(lattice[1]) + " " + str(lattice[2]) + "\n")
                #basis = d1['cartesian_site_positions']
                basis = []
                elements = []
                #print(subd['sites'])
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
            poscar_to_input(calc_type,oqmd_id,compound,False)
            if os.path.isfile("POSCAR"):
                os.system("mv POSCAR POSCAR-{}".format(oqmd_id))
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
    CONDITION = sys.argv[1]
    if os.path.isfile("htepc.json"):
        d = input_data['download']
        START = d['inp']['start']
        END = d['inp']['end']
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
        print("input file htepc.json not found\n")
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
        #with open("download.py", "w") as h:
        #    h.write("info="+str(default1) + "\n")
        #    h.write("inp="+str(default2) + "\n")
        #    h.write("chemsys="+str(chemsys) + "\n")
    if CONDITION == 'search':
        search(KWARGS,PROPERTIES)
    elif CONDITION == 'download':
        download(CALC_TYPE,START,END)
    else:
        print("wrong options. Allowed options are either search or download\n")
