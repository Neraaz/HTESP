#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to prepare QE input files"""
import os
from collections import OrderedDict
import numpy as np
import scipy.linalg as alg
from ase.io import espresso,cif
from ase.cell import Cell
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io import pwscf
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from mp_api.client import MPRester
from kpath import kpath
from check_json import config
input_data = config()
def pos_to_kpt(structure_filename,kpoint_density):
    """
    Obtain k-point mesh from a structure file.

    Parameters:
    ----------------------
    structure_filename : str
        Structure file (e.g., QE scf.in or VASP POSCAR).
    kpoint_density : float
        K-point density.

    Returns:
    ----------------------
    kmesh : list
        K-point mesh according to the k-point density.
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
    return kmesh
class MpConnect:
    """
    Class for connecting to Materials Project (MP) via MP API, extracting properties,
    and preparing input files for Quantum Espresso (QE) calculations.

    Parameters:
    --------------
    key : str, optional
        MP API KEY. Use your own key. If not provided, it should be available in config.json
        or ../../config.json file.

    Attributes:
    --------------
    prop : list
        List of available properties.
    mpid : str
        Materials ID.
    comp : str
        Compound name.
    data : dict
        Dictionary containing materials data.
    prefix : str
        Prefix for the compound.
    ecutwfc : float
        Kinetic energy cutoff for wavefunction.
    ecutrho : float
        Kinetic energy cutoff for charge density.
    kpt : list
        K-point grid.
    structure : Structure
        Structure object.
    pseudo_dir : str
        Path to pseudopotential files.
    outdir : str
        Directory for output files.
    evenkpt : tuple
        K-grid with an even number of points in all directions.
    kpshift : list
        K-point grid shifts.
    kptype : str
        K-point grid type.
    calc : str
        Type of calculation.
    smear : float
        Degauss value for smearing.
    smear_type : str
        Type of smearing.
    etot_conv_thr : float
        Total energy convergence threshold.
    forc_conv_thr : float
        Force convergence threshold.
    conv_thr : float
        Convergence threshold.
    dict_element : dict
        Dictionary containing kinetic energy cutoffs for different elements.
    comp_list : list
        List of elements in the compound.

    Example of using MpConnect class:
    >>> obj = MpConnect()  # Initialize MpConnect object
    >>> obj.setting('mp-763')  # Set Materials ID to 'mp-763'
    >>> obj.getkpt()  # Get k-points for the structure
    >>> obj.maxecut_sssp()  # Calculate maximum recommended energy cutoff using SSSP
    >>> obj.setting_qeinput()  # Set up Quantum Espresso input files based on the retrieved data
    """
    # Intialize the class
    def __init__(self):
        if os.path.isfile("config.json") or os.path.isfile("../../config.json"):
            input_data = config()
            key = input_data["mpi_key"]["API_KEY"]
        else:
            print("config.json file not found. Please provide with your materials project api key\n")
        try:
            self.key = key['key']
            self.mpr = MPRester(self.key)
        except:
            self.key = None
        self.prop = None
        self.mpid = None
        self.comp = None
        self.data = None
        self.prefix = None
        self.ecutwfc = 0
        self.ecutrho = 0
        self.kpt = None
        self.structure = None
        self.pseudo_dir = ""
        self.outdir = ""
        self.evenkpt = None
        self.kpshift = None
        self.kptype = ""
        self.calc = ""
        self.smear = None
        self.smear_type = ""
        self.etot_conv_thr = None
        self.forc_conv_thr = None
        self.conv_thr = None
        self.dict_element = {}
        self.comp_list = None
    def setting(self,comp):
        """
        Initialize the process with a Materials ID.

        Parameters:
        ---------------------
        comp : str
            Materials ID.

        Returns:
        ---------------------
        comp : str
            Materials ID.

        Notes:
        ---------------------
        This function initializes the process with a Materials ID. It retrieves the data related to the given ID from the Materials Project (MP) database, sets up necessary parameters, and prepares the structure for further calculations.
        """
        # Set the Material ID
        self.comp=comp
        # Retrieve available properties from Materials Project database
        # We only store properties except for last 29 elements
        # Mostly related to elastic properties, absent for many systems
        self.prop = self.mpr.materials.summary.available_fields[:-29]
        # Check if config.json file exists
        if os.path.isfile("config.json") or os.path.isfile("../../config.json"):
            # Load settings from input_data dictionary
            d = input_data["download"]
        # Append additional properties specified in the settings to the properties list
        for prop in d['element']['prop']:
            if prop not in self.prop:
                self.prop.append(prop)
        #self.data = self.mpr.materials.summary.get_data_by_id(self.comp,fields=self.prop).dict()
        # Fetch data for the given Materials ID from the Materials Project database
        self.data = self.mpr.materials.summary.search(material_ids=[self.comp],fields=self.prop)[0]
        self.structure = self.data.structure
        self.data = self.data.dict()
        # Extract relevant information from the fetched data
        self.mpid = self.data['material_id']
        symbol_comp = ""
        elm = list(self.data['composition'].keys())
        count = list(self.data['composition'].values())
        # Construct a prefix based on elemental composition
        for i,_ in enumerate(elm):
            symbol_comp += elm[i]+str(int(count[i]))
        self.prefix = symbol_comp
        #print("******************************************\n")
        #print("Use property method to extract the following info\n")
        #print("property('name') or get_properties(['name1', 'name2', ..])\n")
        #print("******************************************\n")
        return comp
    def get_prop_list(self):
        """
        Print the list of properties available when downloading the data.

        Returns:
        ---------------------
        prop_list : dict_keys
            A list of keys representing the available properties.

        Notes:
        ---------------------
        This function retrieves the list of properties available for download from the Materials Project (MP) database for the current Materials ID (mpid). It returns a list of keys that represent the available properties that can be downloaded and accessed for analysis or further processing.
        """
        #return self.mpr.summary.get_data_by_id(self.mpid).dict().keys()
        return self.prop
        #return self.mpr.materials.summary.search(material_ids=[self.mpid],fields=self.prop).dict().keys()
    def download(self,filetype='cif'):
        """
        Download the structure file.

        Parameters:
        ---------------
        filetype : str, optional
            File format for downloading. Default is 'cif'.

        Notes:
        ---------------
        This function downloads the structure file for the current Materials ID (mpid) from the Materials Project (MP) database. It saves the file in the specified format, with the default format being CIF ('.cif'). Other supported formats may be specified as needed.
        """
        if filetype == 'cif':
            unsym_struc = self.structure
            CifWriter(unsym_struc, symprec=0.1).write_file('{}.cif'.format(self.mpid))
    def property(self,name):
        """
         Extract a particular property.

         Parameters:
         ------------------
         name : str
             The name of the property available from the list obtained from the get_prop_list() function.

         Returns:
         ------------------
         value
             The value of the specified property.

         Notes:
         ------------------
         This function extracts the value of a specific property identified by its name from the data retrieved for the current Materials ID (mpid). The property name should be one of the properties listed in the output of the get_prop_list() function.
        """
        return self.data[name]
    def getkpt(self,primitive=True):
        """
        Compute k-points based on k-point density.
        Parameters:
        --------------------
        primitive: logical
                 Use primitive standard structure
        Returns:
        ---------------------
        kpt : list
            The k-point grid.
        kptype : str
            The type of k-point grid.
        kptshift : list
            The shifts in the k-point grid, returns [0,0,0].

        Notes:
        ---------------------
        This function computes the k-point grid based on a specified k-point density,
        defaulting to 0.05 if unspecified. It returns the k-point grid, type, and shifts.
        """
        if primitive:
            struc = SpacegroupAnalyzer(self.structure,symprec=0.1).get_primitive_standard_structure()
        else:
            struc = self.structure
        relax_set = MPRelaxSet(structure=struc)
        relax_set.poscar.write_file('POSCAR')
        try:
            kptden = input_data['kptden']
        except:
            print("Default kptden of 0.05 being utilized\n")
            kptden = 0.05
        self.kpt = pos_to_kpt("POSCAR",kptden)
        self.kptype = "automatic"
        self.kpshift = [0, 0, 0]
        os.system("rm POSCAR")
        print("*********************************\n")
        print("KPOINT with kpoint density of {} \n".format(kptden))
        print("*********************************\n")
        return self.kpt,self.kptype,self.kpshift

    def getevenkpt(self):
        """
        Make the k-point mesh even.

        Returns:
        ---------------
        evenkpt : tuple
            K-grid with an even number of points in all directions.

        Notes:
        ---------------
        Ensures the generated k-point grid has even points in each dimension
        by incrementing odd components to the next even number.
        Returns the resulting even k-point grid as a tuple.
        """
        kptsize = len(self.kpt)
        kpoint_list = self.kpt
        for i in range(kptsize):
            if kpoint_list[i]%2 == 0:
                kpoint_list[i] = kpoint_list[i]
            else:
                kpoint_list[i] = kpoint_list[i] + 1
        self.evenkpt = tuple(kpoint_list)
        return self.evenkpt

    def getecut_sssp(self,element):
        """
        Obtain the kinetic energy cutoff for a specific element.

        Parameters:
        ---------------------
        element : str
            The type of element for which the kinetic energy cutoff is required.

        Returns:
        ----------------------
        float
            Kinetic energy cutoff for the particular element.

        Notes:
        ----------------------
        Gets cutoff values for elements from config.json or defaults from SSSP efficiency set if not found.
        """
        if os.path.isfile("config.json") or os.path.isfile("../../config.json"):
            psd_data = input_data["pseudo"]
            self.dict_element = psd_data['PSEUDO']
        else:
            print("config.json file not found, using default values from SSSP efficiency set\n")
            #print("https://www.materialscloud.org/discover/sssp/table/efficiency\n")
            self.dict_element = {'H': 60, 'Li': 40, 'Be': 40, 'N': 60, 'F': 45, 'Na': 40, 'Mg': 30, 'Al': 30, 'Si': 30, 'P': 30, 'S': 35, 'Cl': 40, 'K': 60, 'Ca': 30, 'Sc': 40, 'Ti': 35, 'V': 35, 'Cr': 40, 'Mn': 65, 'Fe': 90, 'Co': 45, 'Ni': 45, 'Cu': 55, 'Zn': 40, 'Ga': 70, 'Ge': 40, 'As': 35, 'Br': 30, 'Rb': 30, 'Sr': 30, 'Y': 35, 'Zr': 30, 'Nb': 40, 'Mo': 35, 'Tc': 30, 'Ru': 35, 'Rh': 35, 'Pd': 45, 'Ag': 50, 'Cd': 60, 'In': 50, 'Sn': 60, 'Sb': 40, 'Te': 30, 'I': 35, 'Cs': 30, 'Ba': 30, 'La': 40, 'Hf': 50, 'Ta': 45, 'W': 30, 'Re': 30, 'Os': 40, 'Ir': 55, 'Pt': 35, 'Hg': 50, 'Tl': 50, 'Pb': 40, 'Bi': 45, 'B': 35, 'C': 45, 'Au': 45, 'Se': 30, 'O': 60}
        return self.dict_element[element]
    def maxecut_sssp(self):
        """
        Choose the maximum kinetic energy cutoff among elements in a compound.

        Returns:
        -----------------
        tuple
            A tuple containing the maximum kinetic energy cutoff for waveFunction (ecutwfc)
            and the maximum kinetic energy cutoff for charge density (ecutrho).

        Notes:
        -----------------
        This function finds the maximum kinetic energy cutoff for compound elements,
        then returs the cutoffs for waveFunction and density.
        """
        # Try to retrieve elements from the fetched data
        try:
            self.comp_list = self.data['elements']
            # Get the kinetic energy cutoff for each element in the compound
            cutoff_list = [self.getecut_sssp(el) for el in self.comp_list]
        except KeyError:
            # If 'elements' key is not found in the fetched data, extract elements from the structure
            self.comp_list = [str(el) for el in self.structure.elements]
            cutoff_list = [self.getecut_sssp(el) for el in self.comp_list]
        except:
            elements_spin = self.structure.elements
            cutoff_list = []
            for elm_spin in elements_spin:
                elm_spin = str(elm_spin)
                # Check if the element contains additional spin information (e.g., 'Mn, Spin=5')
                if "," in elm_spin:
                    # If additional spin information is present, extract the element name
                    cutoff_list.append(self.getecut_sssp(elm_spin.split(",")[0]))
                else:
                    cutoff_list.append(self.getecut_sssp(elm_spin))
        # Determine the maximum kinetic energy cutoff for waveFunction
        self.ecutwfc= max(cutoff_list)
        # Determine the maximum kinetic energy cutoff for charge density
        self.ecutrho = 8 * self.ecutwfc
        print("******************************************\n")
        print("SSSP tested K.E. cutoffs for waveFunction and density\n")
        print("******************************************\n")
        return self.ecutwfc,self.ecutrho
    def maxecut_sssp_for_subs(self):
        """
        Determine the maximum kinetic energy cutoff among elements in a compound during the substitution process.

        Returns:
        -----------------
        tuple
            A tuple containing the maximum kinetic energy cutoff for waveFunction (ecutwfc)
            and the maximum kinetic energy cutoff for charge density (ecutrho).

        Notes:
        -----------------
        Similar to maxecut_sssp but for substitution process.
        """
        try:
            cutoff_list = [self.getecut_sssp(el) for el in self.comp_list]
        except:
            elements_spin = self.structure.elements
            cutoff_list = []
            for elm_spin in elements_spin:
                elm_spin = str(elm_spin)
                if "," in elm_spin:
                    cutoff_list.append(self.getecut_sssp(elm_spin.split(",")[0]))
                else:
                    cutoff_list.append(self.getecut_sssp(elm_spin))
        self.ecutwfc= max(cutoff_list)
        self.ecutrho = 8 * self.ecutwfc
        return self.ecutwfc,self.ecutrho

    def ecut_set(self,ecutwfc=50.0,ecutrho=400.0):
        """
        Function to set kinetic energy cutoffs
        """
        self.ecutrho=ecutrho
        self.ecutwfc=ecutwfc

    def get_properties(self,property_name):
        """
        Function to extract multiple properties from MP.
        parameters
        -------------------
        property_name : list of properties: Default: ['material_id']
        Returns
        -------------------
        property_list : list of properties extracted.
        """
        property_list = []
        for prop in property_name:
            property_list.append(self.data[prop])
        property_list.insert(0,self.data['material_id'])
        with open('mpid.csv', 'a') as data:
            for prop in property_list:
                data.write(str(prop) + ",")
            data.write("\n")
        return property_list
    def setting_qeinput(self,calculation='vc-relax',occupations='smearing',restart_mode='from_scratch',pseudo_dir='./',smearing=0.02,smearing_type='gauss',etot_conv_thr=1e-05,forc_conv_thr=1e-04,conv_thr=1e-16,ion_dynamics='bfgs',cell_dynamics='bfgs',magnetic=False,primitive=True):
        """
        Function to create input file for QE ground-state calculations.

        Parameters:
        - calculation (str): Type of calculation. Default: 'vc-relax'. Other options are 'relax' (ionic only), 'bands' for band structure, 'scf' for SCF calculations.
        - occupations (str): Occupation. Default: 'smearing'. Other options could be 'tetrahedra' and so on.
        - restart_mode (str): How to start the calculations.
        - pseudo_dir (str): Path to pseudopotential files. Default: './' (current directory).
        - smearing (float): Degauss value. Default: 0.02.
        - smearing_type (str): Type of smearing. Default: 'gauss'.
        - etot_conv_thr, forc_conv_thr, conv_thr (float): Convergence parameters of QE calculations.
        - ion_dynamics, cell_dynamics (str): Algorithm to perform relaxation. Default: 'bfgs'.
        - magnetic (logical): Magnetic flag if magnetic input need to be created
        - primitive (logical): if True, then use primitive structure

        Returns:
        Creates input files in scf-mpid.in format inside scf_dir/.
        """
        # If magnetic calculation is enabled, handle pseudopotentials and magnetizations
        if magnetic:
            elements_spin = self.structure.elements
            pseudo1 = {}
            #magnetization = []
            magdict = OrderedDict()
            new_species = []
            new_coord = []
            for site in self.structure:
                new_species.append(site.label)
                new_coord.append(site.frac_coords)
            latest_species,_ = convert_species_list(new_species)
            pseudo_species = []
            for elm_spin in elements_spin:
                pseudo_species.append(str(elm_spin))
            pseudo_species_new,magdict = convert_species_list(pseudo_species)
            pseudo_change_dict = {}
            for i,_ in enumerate(pseudo_species):
                pseudo_change_dict[pseudo_species[i]] = pseudo_species_new[i]
            for i,elm_spin_i in enumerate(elements_spin):
                elm_spin = str(elm_spin_i)
                if "," in elm_spin:
                    elm_s = elm_spin.split(",")[0]
                    pseudo1[elm_spin] = elm_s + ".upf"
                else:
                    pseudo1[elm_spin] = elm_spin + ".upf"
        else:
            pseudo1 = {el:el+'.upf' for el in self.comp_list}
        # Set prefix for output files
        prefix = self.prefix
        #for i,_ in enumerate(latest_species):
        #    magdict[latest_species[i]] = magnetization[i]
        # Set up control, system, and electrons parameters
        if os.path.isfile("config.json") or os.path.isfile("../../config.json"):
            pwscf_in = input_data["pwscf_in"]
            control = pwscf_in['control']
            system = pwscf_in['system']
            electrons = pwscf_in['electrons']
        else:
            control = {'calculation':calculation, 'nstep':300, 'restart_mode':restart_mode, 'pseudo_dir':pseudo_dir, 'outdir':'./', 'tprnfor':'.true.','tstress':'.true.', 'etot_conv_thr':etot_conv_thr, 'forc_conv_thr':forc_conv_thr}
            system = {'smearing':smearing_type, 'occupations':occupations, 'degauss':smearing}
            electrons = {'diagonalization':'david', 'mixing_mode':'plain', 'mixing_beta':0.7, 'conv_thr': conv_thr, 'electron_maxstep':300}
        # Set kinetic energy cutoffs
        system['ecutwfc'] = self.ecutwfc
        system['ecutrho'] = self.ecutrho
        control['prefix'] = prefix
        # lW  2010_SC has monoclinic conventional cell with alpha<90, beta=gamma=90, different from international beta=90
        if primitive:
            tmpanalyzer = SpacegroupAnalyzer(self.structure)
            tmpstructure = tmpanalyzer.get_primitive_standard_structure(international_monoclinic=False)
            self.structure = tmpstructure
        # Generate the input file based on the chosen calculation type
        if calculation == 'vc-relax':
            filename = pwscf.PWInput(self.structure, pseudo=pseudo1, control=control, system=system,electrons=electrons, kpoints_grid=self.kpt, ions={'ion_dynamics':ion_dynamics}, cell={'cell_dynamics':cell_dynamics,'press_conv_thr':0.05})
        elif calculation == 'relax':
            filename = pwscf.PWInput(self.structure, pseudo=pseudo1, control=control, system=system,electrons=electrons, kpoints_grid=self.kpt, ions={'ion_dynamics':ion_dynamics})
        else:
            filename = pwscf.PWInput(self.structure, pseudo=pseudo1, control=control, system=system,electrons=electrons, kpoints_grid=self.kpt)
        # Write input file
        filename.write_file("temp.in")
        # Adjust the input file format and save it with the appropriate filename
        os.system("""sed "s/'.true.'/.true./" {} > {}""".format("temp.in","scf-{}.in".format(self.mpid)))
        os.system("sed -n '/K_POINTS automatic/,/CELL_PARAMETERS angstrom/p' scf-{}.in | sed '$d' > kpoint-{}.dat".format(self.mpid,self.mpid))
        obj = INPUTscf("scf-{}.in".format(self.mpid))
        # Use symmetric structure
        obj.standardize(self.mpid,output="temp.dat")
        os.system("sed -n '/&CONTROL/,/ATOMIC_SPECIES/p' scf-{}.in | sed '$d' > scf.header".format(self.mpid))
        os.system("sed -n '/ATOMIC_SPECIES/,/ATOMIC_POSITIONS crystal/p' scf-{}.in | sed '$d' > species".format(self.mpid))
        os.system("rm temp.in")
        if calculation == 'bands':
            os.system("sed -i '/K_POINTS/{N;d;}' temp.dat")
            obj = INPUTscf(filename="scf-{}.in".format(self.mpid))
            obj.generate_kpath(nqpoint=200,kcut=0,out="kpoint.dat")
            os.system("cat kpoint.dat temp.dat > temp2.dat")
            os.system("mv temp2.dat temp.dat")
            os.system("cat scf.header species temp.dat > scf-{}.in".format(self.mpid))
            print("Input for band calculation. Provide nband = <n> within SYSTEM section. Use sufficient <n>.")
        else:
            os.system("cat scf.header species temp.dat > scf-{}.in".format(self.mpid))
        # If magnetic calculation is enabled, modify the input file accordingly
        # to incorporate magnetic keywords in the input files
        if magnetic:
            maglist = list(magdict.values())
            os.system(f"""sed -i '/&SYSTEM/a nspin = 2,' scf-{self.mpid}.in""")
            nsites = len(maglist)
            with open(f"scf-{self.mpid}.in", "r") as readscf:
                scflines = readscf.readlines()
            old_species = []
            old_coord = []
            natom = 0
            for i,scfl in enumerate(scflines):
                if 'ATOMIC_POSITIONS' in scfl:
                    atom_start = i
                if 'ATOMIC_SPECIES' in scfl:
                    pseudo_start = i
            while natom < len(latest_species):
                old_species.append(scflines[atom_start+1+natom].split(" ")[0])
                old_x = float(scflines[atom_start+1+natom].split(" ")[1])
                old_y = float(scflines[atom_start+1+natom].split(" ")[2])
                old_z = float(scflines[atom_start+1+natom].split(" ")[3])
                old_coord.append([old_x,old_y,old_z])
                natom += 1
            new_list = []
            for k,_ in enumerate(pseudo_species):
                old = scflines[pseudo_start+1+k].split(" ")
                old = list(filter(lambda x: x != '', old))[0]
                new = pseudo_change_dict[old]
                new_list.append(new)
                os.system(f"""sed -i '{pseudo_start+2+k}s/{old}/{new}/' scf-{self.mpid}.in""")
            maglist = list(reorder_dictionary(magdict,new_list).values())
            for j,_ in enumerate(old_species):
                os.system(f"""sed -i '{atom_start+2+j}s/{old_species[j]}/{latest_species[j]}/' scf-{self.mpid}.in""")
                os.system(f"""sed -i '{atom_start+2+j}s/{old_coord[j][0]} {old_coord[j][1]} {old_coord[j][2]}/{new_coord[j][0]} {new_coord[j][1]} {new_coord[j][2]}/' scf-{self.mpid}.in""")
            for j,_ in enumerate(maglist):
                os.system(f"""sed -i '/&SYSTEM/a starting_magnetization({nsites-j}) = {maglist[nsites-j-1]}' scf-{self.mpid}.in""")
        # Clean up temporary files
        os.system("rm scf.header species temp.dat kpoint*")
def reorder_dictionary(original_dict, new_order):
    """
    Reorder the keys of a dictionary based on a new list of keywords and return the list of keys.

    Parameters:
    - original_dict (dict): The original dictionary whose keys are to be reordered.
    - new_order (list): The new list of keywords specifying the desired order of keys.

    Returns:
    - list: The list of keys of the reordered dictionary based on the new order.

    Example:
    >>> original_dict = {'a': 1, 'b': 2, 'c': 3, 'd': 4}
    >>> new_order = ['c', 'a', 'd', 'b']
    >>> reorder_dictionary(original_dict, new_order)
    """
    # Create a new OrderedDict with the keys ordered according to the new list
    ordered_dict = OrderedDict((key, original_dict[key]) for key in new_order if key in original_dict)

    # Return the list of keys in the new order
    return ordered_dict


def extract_elements_with_spin(original_list):
    """
    Extract elements with their spin values from the given list.

    Parameters:
    - original_list (list): List containing strings representing elements with spin values.

    Returns:
    - dict: A dictionary where keys are elements and values are lists of spin values.
    """
    elements_with_spin = {}
    for item in original_list:
        if ',' in item:
            element, spin = item.split(',')
            spin_value = spin.split('=')[1]
            if element not in elements_with_spin:
                elements_with_spin[element] = [spin_value]
            elif spin_value not in elements_with_spin[element]:
                elements_with_spin[element].append(spin_value)
    return elements_with_spin

def convert_species_list(original_list):
    """
    Rename elements in the original list with elements with different spins.

    Parameters:
    - original_list (list): List containing strings representing elements.

    Returns:
    - list: Updated list where elements are replaced with generic names for different spins.

    Example:
    >>> original_list = ['Fe,spin=5', 'Fe,spin=-5', 'pd', 'pd', 'I,spin=1', 'I,spin=-1']
    >>> convert_species_list(original_list)
    ['Fe1', 'Fe2', 'pd', 'pd', 'I1', 'I2']
    """
    elements_spin = extract_elements_with_spin(original_list)
    # Updated list
    updated_list = []
    magdict = {}
    # Iterate through the original list
    for item in original_list:
        if ',' in item:
            element = item.split(',')[0]
            spin = float(item.split('=')[-1])
            if len(elements_spin[element]) > 1 and spin > 0:
                updated_list.append(element + str(1))
                magdict[element + str(1)] = 1
            elif len(elements_spin[element]) > 1 and spin < 0:
                updated_list.append(element + str(2))
                magdict[element + str(2)] = -1
            else:
                updated_list.append(element)
                spin_s = spin
                mag_s = 1 if spin_s > 0 else -1 if spin_s < 0 else 0
                magdict[element] = mag_s
        else:
            updated_list.append(item)
            magdict[item] = 0
    return updated_list,magdict


def ase_cell_to_structure(ase_cell):
    # Extract lattice vectors from ASE Atoms object
    lattice_vectors = ase_cell.cell.tolist()
    
    # Extract atomic positions and species from ASE Atoms object
    atomic_positions = ase_cell.get_scaled_positions()
    species_symbols = ase_cell.get_chemical_symbols()
    
    # Create pymatgen Lattice object
    lattice = Lattice(lattice_vectors)
    
    # Create list of pymatgen Specie objects
    species = [Element(symbol) for symbol in species_symbols]
    
    # Create pymatgen Structure object
    structure = Structure(lattice, species, atomic_positions)
    
    return structure

class INPUTscf:
    """
    class to process QE input files
    parameters
    ------------------
    filename: (str) QE input file
    """
    def __init__(self,filename='scf.in'):
        self.filename = filename
        self.ase_cell = espresso.read_espresso_in(self.filename)
        # Convert the ASE cell.io.espresso object to a pymatgen Structure object
        structure = ase_cell_to_structure(self.ase_cell)
        # Use SpacegroupAnalyzer to get the symmetric structure
        symmetry_analyzer = SpacegroupAnalyzer(structure)
        symmetric_structure = symmetry_analyzer.get_symmetrized_structure()
        # Now, symmetric_structure contains the symmetrized version of the structure
        self.struc = symmetric_structure
        self.cell = symmetric_structure.lattice.matrix
        self.volume = symmetric_structure.lattice.volume
        #self.braiv_latt = Cell.get_bravais_lattice(self.cell)
        self.mpid = None
        self.comp = None
        self.prefix = None
        self.qpoint = []
        self.mass = []

    def scftocif(self,output='file.cif'):
        """
        Function to convert QE input file to structure file in .cif format
        parameters
        ---------------
        output : (str) name of the output file
        Returns
        ---------------
        file2 : output file object
        """
        self.struc.to(output)
        return self.struc

    def cellpar(self):
        """
        Function to calculate cell lengths and angles
        Returns
        -----------------------
        list of lenghts and angles
        """
        length = self.struc.lattice.abc
        angles = self.struc.lattice.angles
        return list(length) + list(angles)

    def standardize(self,mpid,output="standard.in"):
        """
        Function to get symmetrized structure.
        parameters
        -----------------
        mpid : (str) materials project ID
        output : (str) output file. Default: 'standard.in'
        """
        finalpos = self.struc.frac_coords
        finalcell = self.cell
        specieslist = self.struc.species
        with open("kpoint-{}.dat".format(mpid), "r") as kpoint:
            kplines = kpoint.readlines()
        with open(output, "w") as struc:
            struc.write("ATOMIC_POSITIONS crystal\n")
            for i,_ in enumerate(specieslist):
                struc.write(str(specieslist[i]) + " " + str(finalpos[i][0])+ " ")
                struc.write(str(finalpos[i][1]) + " " + str(finalpos[i][2]) + "\n")
            for i in range(2):
                struc.write(kplines[i])
            struc.write("CELL_PARAMETERS angstrom\n")
            for i in range(3):
                struc.write(str(finalcell[i][0]) + " " + str(finalcell[i][1]) + " " + str(finalcell[i][2]) + "\n")

    def generate_kpath(self,nqpoint=200,kcut=0,out='kpath.in'):
        """
        Function to generate and write k-point mesh for bandstructure calculation to a file
        parameters
        --------------------
        nqpoint : (int) size of the k-point mesh. Default: 200
        kcut : (int) cutoff to the high-symmetry path of the Brillouin zone. Default: 0 for full Brillouin zone
        out : (str) output file. Default: 'kpath.in'
        Returns
        --------------------
        kpts : numpy array of kpoints in linear axis after processing
        n : (int) size of kpts
        """
        kpts,_,_,_,_,_ = kpath(self.filename,nqpoint,kcut)
        nkpt = kpts.shape[0]
        with open(out, 'w') as kmesh:
            kmesh.write('K_POINTS\n')
            kmesh.write(str(nkpt) + '\n')
            for i in range(nkpt):
                kmesh.write(str(round(kpts[i][0],8)) + " " + str(round(kpts[i][1],8)) + " " + str(round(kpts[i][2],8)) + " " + str(0.0) + "\n")
        return nkpt,kpts

    def setting_input(self,comp,mass,qpoint):
        """
        Function to setup input
        parameters
        -------------
        comp : compound name
        mass: List of masses of elements
        qpoint : List of qpoint
        """
        self.comp=comp
        self.mass=mass
        self.qpoint=qpoint
        obj=MpConnect()
        obj.setting(self.comp)
        self.mpid = obj.mpid
        self.prefix = obj.prefix

    def create_elph(self,output='elph.in',tol=0.00000000000001,sigma=0.005):
        """
         Function similar to elph.py file inside src/
        """
        dynmat = self.prefix.replace("'", "") + ".dyn"
        nat = len(self.mass)
        with open(output, 'w') as elph:
            elph.write("electron phonon coupling \n")
            elph.write("&inputph" + "\n")
            elph.write("tr2_ph={},".format(tol) + "\n")
            elph.write("prefix={},".format(self.prefix) + "\n")
            elph.write("fildvscf='aldv'," + "\n")
            for i in range(1,nat+1):
                elph.write("amass({})={},".format(i,self.mass[i-1]) + "\n")
            elph.write("outdir='./'," + "\n")
            elph.write("fildyn='{}',".format(dynmat) + "\n")
            elph.write("electron_phonon='interpolated'," + "\n")
            elph.write("el_ph_sigma={},".format(sigma) + "\n")
            elph.write("el_ph_nsigma=10," + "\n")
            elph.write("trans=.true.," + "\n")
            elph.write("ldisp=.true." + "\n")
            elph.write("nq1={},nq2={},nq3={}".format(int(self.qpoint[0]),int(self.qpoint[1]),int(self.qpoint[2])) + "\n")
            elph.write("/" + "\n")

    def create_q2r(self,output='q2r.in'):
        """
         Function similar to q2r.py file inside src/
        """
        dynmat = self.prefix.replace("'", "") + ".dyn"
        frc = self.prefix.replace("'", "") + ".fc"
        with open(output, 'w') as q2r_write:
            q2r_write.write("&input" + "\n")
            q2r_write.write("zasr='simple'," + "\n")
            q2r_write.write("fildyn='{}',".format(dynmat) + "\n")
            q2r_write.write("flfrc='{}',".format(frc) + "\n")
            q2r_write.write("la2F=.true." + "\n")
            q2r_write.write("/" + "\n")

    def create_matdyn(self,out='matdyn.in',nqpt=60,kcut=0):
        """
         Function similar to matdyn.py file inside src/
        """
        freq = self.prefix.replace("'", "") + ".freq"
        frc = self.prefix.replace("'", "") + ".fc"
        eig = self.prefix.replace("'", "") + ".eig"
        nat = len(self.mass)
        nkpt,kpts = self.generate_kpath(nqpoint=nqpt,kcut=kcut)
        with open(out, 'w') as matdyn_write:
            matdyn_write.write("&input" + "\n")
            matdyn_write.write("asr='simple'," + "\n")
            for i in range(1,nat+1):
                matdyn_write.write("amass({})={},".format(i,self.mass[i-1]) + "\n")
            matdyn_write.write("flfrc='{}',".format(frc) + "\n")
            matdyn_write.write("flfrq='{}',".format(freq) + "\n")
            matdyn_write.write("fleig='{}',".format(eig) + "\n")
            matdyn_write.write("la2F=.true.," + "\n")
            matdyn_write.write("dos=.false." + "\n")
            matdyn_write.write("/" + "\n")
            matdyn_write.write(str(nkpt) + "\n")
            for i in range(nkpt):
                matdyn_write.write(str(round(kpts[i][0],8)) + " " + str(round(kpts[i][1],8)) + " ")
                matdyn_write.write(str(round(kpts[i][2],8)) + " " + str(0.0) + "\n")
    def create_phdos(self,kpts,out='phdos.in',ndos=200):
        """
         Function similar to matdyn_dos.py file inside src/
        """
        freq = self.prefix.replace("'", "") + "-dos.freq"
        frc = self.prefix.replace("'", "") + ".fc"
        nat = len(self.mass)
        with open(out, 'w') as phdos:
            phdos.write("&input" + "\n")
            phdos.write("asr='simple'," + "\n")
            for i in range(1,nat+1):
                phdos.write("amass({})={},".format(i,self.mass[i-1]) + "\n")
            phdos.write("flfrc='{}',".format(frc) + "\n")
            phdos.write("flfrq='{}',".format(freq) + "\n")
            phdos.write("la2F=.true.," + "\n")
            phdos.write("dos=.true.," + "\n")
            phdos.write("fldos='phonon.dos'," + "\n")
            phdos.write("nk1={},nk2={},nk3={},ndos={},".format(kpts[0],kpts[1],kpts[2],ndos) + "\n")
            phdos.write("/" + "\n")

    def create_dos(self,out1='dos.in',out2='pdos.in'):
        """
         Function similar to dos.py file inside src/
        """
        dynmat = self.prefix.replace("'", "") + ".dos"
        dynmat1 = self.prefix.replace("'", "") + ".pdos"
        with open(out1, 'w') as dos:
            dos.write("&dos" + "\n")
            dos.write("prefix={},".format(self.prefix) + "\n")
            dos.write("outdir='./'," + "\n")
            dos.write("fildos='{}',".format(dynmat) + "\n")
            dos.write("DeltaE=0.01" + "\n")
            dos.write("/" + "\n")
        with open(out2, 'w') as pdos:
            pdos.write("&projwfc" + "\n")
            pdos.write("prefix={},".format(self.prefix) + "\n")
            pdos.write("outdir='./'," + "\n")
            pdos.write("pfildos='{}',".format(dynmat1) + "\n")
            pdos.write("DeltaE=0.01" + "\n")
            pdos.write("/" + "\n")

    def post_band(self,out='band.in'):
        """
         Function similar to band.py file inside src/
        """
        dynmat = self.prefix.replace("'", "") + ".dat"
        with open(out, 'w') as band:
            band.write("&BANDS" + "\n")
            band.write("prefix={},".format(self.prefix) + "\n")
            band.write("outdir='./'," + "\n")
            band.write("filband='{}',".format(dynmat) + "\n")
            band.write("lsym=.true." + "\n")
            band.write("/" + "\n")

    def post_phband(self,out='phonband.in'):
        """
         Function similar to phonband.py file inside src/
        """
        freq = self.prefix.replace("'", "") + ".freq"
        with open(out, 'w') as phband:
            phband.write(freq + "\n")
            phband.write("0 5000" + "\n")
            phband.write("freq.plot" + "\n")
            phband.write("freq.ps" + "\n")
            phband.write("0.0" + "\n")
            phband.write("100.0 0.0" + "\n")

def scf_to_dos_scf(input_file='scf.in',out='scf-dos.in'):
    """
    Function to change QE input file from using smearing to 'tetrahedra' method, mainly for DOS calculation

    parameters
    ----------------------
    input_file : (str) input file for QE scf calculation
    out : (str) output file after modification

    """
    os.system("""sed -i '/calculation/d' {} """.format(input_file))
    insert(input_file=input_file,output='temp',keyword='&CONTROL',what="calculation = 'nscf',")
    os.system("""cat {} | sed "s/'smearing'/'tetrahedra'/" | sed '/degaus/d' | sed '/smearing/d' > {}""".format('temp',out))
    os.system("""rm temp""")
    print("Use denser k-mesh for dos calculations")

def insert(input_file='scf.in',output='scf-new.in',keyword=None,where="after",what=""):
    """
    Function to insert keyword in input file
    parameters
    -----------------------
    input_file : (str) input file. Default: 'scf.in'
    output : (str) output file. Default: 'scf-new.in'
    keyword : (str) keyword to look after
    where : (str) where to insert. Default: 'after', otherwise 'before'
    what : (str) what to insert. Default: ''
    """
    if where == "after":
        os.system("""sed "/{}/a   {}" {} > {}""".format(keyword,what,input_file,output))
    else:
        os.system("""sed "/{}/i   {}" {} > {}""".format(keyword,what,input_file,output))
# Extract relaxed structure and update QE input.
class OUTPUTscf:
    """
    class to extract relax structure from QE output file and update input file
    parameters
    ---------------------
    filename : (str) QE output file. Default: 'scf.out'
    """
    def __init__(self,filename='scf.out'):
        self.filename = filename

    def extract_relax(self,output='relax.dat'):
        """
        Function to extract cell and positions of crystal structures from scf output file and write to a file
        parameters
        ----------------
        filename : (str) output file. Default: 'relax.dat'
        """
        os.system("sed -n '/Begin final coordinates/,/End final coordinates/p' {} | sed '$d' | sed '1,4d'| sed '5d' > {}".format(self.filename,output))
    def update_scf(self,input_file='scf.in',output='scf-new.in',what='structure'):
        """
        Function to update QE input file and update with new structure
        parameters
        ------------
        input_file : (str) QE scf input file to update. Default: 'scf.in'
        output : (str) QE scf output file, updated with new structure. Default: 'scf-new.in'
        what : (str) what to update. Default: 'structure'. Other, not implemented yet !
        """
        os.system("sed -n '/K_POINTS automatic/,/CELL_PARAMETERS angstrom/p' {} | sed '$d' > kpoint.dat".format(input_file))
        os.system("sed -n '/&CONTROL/,/ATOMIC_POSITIONS crystal/p' {} | sed '$d'  > scf.header".format(input_file))
        os.system("sed -n '/ATOMIC_SPECIES/,/ATOMIC_POSITIONS crystal/p' {} | sed '$d' > species".format(input_file))
        if what == 'structure':
            self.extract_relax('relax.in')
            os.system("cat scf.header kpoint.dat relax.in > {}".format(output))
        else:
            print('do nothing\n')
        os.system("rm scf.header species kpoint.dat relax.in")

#def extract_energy(self,input='scf.out'):
#    os.system("""Rytoev=`echo "scale=6;13.605698" | bc`""")
#    os.system("""en=`grep "!    total energy              =   " {} | tail -n1 | awk '{print $5}'` """.format(input))
#    os.system("""e=`echo "scale=6; $en * $Rytoev  " | bc`""")
#    os.system("""echo "the total energy: $e" """)
#    return
#class for calculations. In progress........
#class calculation:
#    def __init__(self,batch_header='batch.header'):
#        self.batch_header = batch_header
#def create_submission_scripts(self):
#def relax_structure(self,input='scf.in',out='scf.out'):
#def calculate_Tc(self,mu=0.16,degaussq=0.12,ngaussq=0,phfile='elph.out',phsigma=0.005):
