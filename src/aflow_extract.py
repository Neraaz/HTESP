#!/usr/bin/env python
#
#Written by Niraj K. Nepal, PhD.
"""
Module to search and download input files from AFLOW database.
"""
# coding: utf-8
# Extracting Carbon compounds
import os
import sys
import ast
import json
from urllib.request import urlopen
import pandas as pd
#from ase.cell import Cell
from pymatgen.core.structure import IStructure
from pymatgen.core import structure, lattice
from pymatgen.core.composition import Composition
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from oqmd_extract import poscar_to_input
from cif_to_gsinput import pos_to_kpt
from write_potcar import poscar2potcar
from htepc import INPUTscf
try:
    PWD = os.getcwd()
    if os.path.isfile(PWD+"/config.json"):
        JSONFILE = PWD+"/config.json"
    else:
        JSONFILE = "../../config.json"
    with open(JSONFILE, "r") as readjson:
        input_data = json.load(readjson)['download']
except FileNotFoundError:
    print("config.json file not found\n")
def properties_string(dict_param):
    """
    Function to create a string for URL search from the dictionary provided in config.json.

    Parameters:
    dict_param (dict): A dictionary from the config.json file containing search criteria.

    Returns:
    str: A string representing the search criteria for URL.

    Example:
        >>> criteria_dict = {
        ...    "aflow": {
        ...        "filter": True,
        ...        "elm": ["Cu", "Zn"],
        ...        "nsites": 4,
        ...        "metal": True,
        ...        "FE": False,
        ...        "spacegroup": None,
        ...        "limit": 10,
        ...        "prop": ["energy", "bandgap"]
        ...    }
        ... }
        >>> criteria_string = properties_string(criteria_dict)
        >>> print(criteria_string)
        'species(Cu,Zn),natoms(1*,*4),Egap(0),enthalpy_formation_atom(*),spacegroup_relax(*),
        $paging(1,10),energy,bandgap'
    """
    criteria = ""
    search = dict_param['aflow'] # Extract aflow dictionary from dict_param
    apply_filter = search['filter'] # Check if filtering is applied
    if apply_filter:
        for key in search:
            if key == 'elm': # Element filter
                if len(tuple(search[key])) == 1:
                    criteria += "species({}),".format(tuple(search[key])[0])
                else:
                    criteria += "species{},".format(tuple(search[key]))
            elif key == 'nelm': # Number of elements filter
                criteria += "nspecies(1*,*{}),".format(search[key])
            elif key == 'nsites': # Number of sites filter
                criteria += "natoms(1*,*{}),".format(search[key])
            elif key == 'metal': # Metal filter
                if search[key]:
                    criteria += "Egap(0),"
                else:
                    criteria += "Egap(*),"
            elif key == 'FE': # Formation energy filter
                if search[key]:
                    criteria += "enthalpy_formation_atom(-100*,*0),"
                else:
                    criteria += "enthalpy_formation_atom(*),"
            elif key == 'spacegroup': # Spacegroup filter
                if search[key] is None:
                    criteria += "spacegroup_relax(*),"
                else:
                    criteria += "spacegroup_relax({}),".format(search[key])
            elif key == 'limit': # Limit filter for number of results filter
                criteria += "$paging(1,{}),".format(search[key])
            elif key == 'prop': # Additional properties to include
                prop_string = ""
                for i,_ in enumerate(search[key]):
                    if i < len(search[key]) - 1:
                        prop_string += search[key][i] + ","
                    else:
                        prop_string += search[key][i]
            else:
                print("Invalid key provided\n")
    else:
        for key in search:
            if key == 'elm':
                if len(tuple(search[key])) == 1:
                    criteria += "species({}),".format(tuple(search[key])[0])
                else:
                    criteria += "species{},".format(tuple(search[key]))
            elif key == 'nelm':
                criteria += "nspecies(1*,*{}),".format(search[key])
            elif key == 'prop':
                prop_string = ""
                for i,_ in enumerate(search[key]):
                    if i < len(search[key]) - 1:
                        prop_string += search[key][i] + ","
                    else:
                        prop_string += search[key][i]
            else:
                print("Invalid key provided\n")
    criteria = criteria + prop_string
    criteria = criteria.replace(" ","") # Remove any spaces in the criteria string
    return criteria
def search_data(dict_param):
    """
    Function to search data via AFLOW database.

    Parameters:
    dict_param (dict): A dictionary containing search criteria.

    Returns:
    None

    This function constructs a URL based on the search criteria provided
    in dict_param and retrieves data from the AFLOW database. It saves the
    retrieved data as a CSV file in the 'download' directory and also
    generates an 'mpid-list.in' file containing relevant data.

    Example:
        >>> search_criteria = {
        ...    "aflow": {
        ...        "filter": True,
        ...        "elm": ["Cu", "Zn"],
        ...        "nsites": 4,
        ...        "metal": True,
        ...        "FE": False,
        ...        "spacegroup": None,
        ...        "limit": 10,
        ...        "prop": ["energy", "bandgap"]
        ...    }
        ... }
        >>> search_data(search_criteria)
    """
    aflow_url = "http://www.aflowlib.org/API/aflux/"
    criteria = properties_string(dict_param) # Generating criteria string for URL
    search_url = aflow_url + "?" + criteria # Constructing the search URL
    success = False
    # Attempting to fetch data from the AFLOW database
    while not success:
        try:
            data = urlopen(search_url).read().decode('utf-8')
            success = True
        except:
            success = False
            continue
    data_processed = ast.literal_eval(data) # Processing retrieved data
    data_processed = pd.DataFrame(data_processed) # Converting to panda dataframe

    # Creating download directory if doesn't exists
    if not os.path.isdir("download"):
        os.mkdir("download")
    # Saving data as the CSV file in the 'download' directory
    if not os.path.isfile("download/download-aflow.csv"):
        data_processed.to_csv("download/download-aflow.csv",index=False)
    # Generating 'mpid-list.in' file
    if not os.path.isfile("mpid-list.in"):
        with open("mpid-list.in","w") as write_mpid:
            for i in range(data_processed.shape[0]):
                write_mpid.write("v{}".format(i+1) + " " + data_processed['auid'][i] + " " + data_processed['compound'][i] + "\n")
def poscar_create(data):
    """
    Function to create a POSCAR file from data obtained from a query.

    Parameters:
    data (dict): A dictionary containing information necessary for creating the POSCAR file.

    Returns:
    Structure: A pymatgen Structure object representing the structure.

    This function takes data including lattice parameters, compound information,
    fractional positions, and creates a pymatgen Structure object. It then returns
    the Structure object representing the structure. This function currently does
    not write the POSCAR file to disk but can be extended to do so.

    Example:
        >>> data = {
        ...    "geometry": [3.84, 3.84, 3.84, 90, 90, 90],
        ...    "compound": "Si",
        ...    "auid": "aflow:12345",
        ...    "positions_fractional":
               [[0.000000, 0.000000, 0.000000], [0.250000, 0.250000, 0.250000]]
        ... }
        >>> structure = poscar_create(data)
        >>> print(structure)
        Full Formula (Si2)
        Reduced Formula: Si
        abc   :   3.840198   3.840198   3.840198
        angles:  90.000000  90.000000  90.000000
        Sites (2)
          #  SP           auid
        ---  ----  -------------
          0  Si    aflow:12345
          1  Si    aflow:12345
    """
    # Create lattice from lattice parameters using ASE
    lat_data = data['geometry']
    a = lat_data[0]
    b = lat_data[1]
    c = lat_data[2]
    alpha = lat_data[3]
    beta = lat_data[4]
    gamma = lat_data[5]
    latt = lattice.Lattice.from_parameters(a,b,c,alpha,beta,gamma)
    lattice_matrix = latt.matrix
    comp = Composition(data['compound']) # Extract compound information
    aflow_id = data['auid'] # Extract AFLOW ID
    species = []
    comp = comp.as_dict()
    # Determine species from compound information
    for key in comp:
        for i in range(int(comp[key])):
            species.append(key)
    # Create structure object
    struc = IStructure(lattice_matrix,species=species,coords=data['positions_fractional'],to_unit_cell=True)
    struc = structure.Structure(struc.lattice,struc.species,struc.frac_coords,to_unit_cell=True)
    struc = SpacegroupAnalyzer(struc, symprec=0.1).get_primitive_standard_structure()
    # Add AFLOW ID as a site property
    site_prop = [aflow_id] * len(struc.sites)
    struc.add_site_property("auid",site_prop)
    return struc
def download(calc_type,start,end):
    """
    Function to create input files.

    Parameters:
    calc_type (str): Type of calculations, either 'QE' or 'VASP'.
    start (int): Starting index.
    end (int): Ending index.

    Returns:
    None

    This function retrieves data from the AFLOW database based on the specified
    range of indices, constructs input files for quantum mechanical calculations,
    and saves them. It also creates a file named 'mpid.in' containing information
    about the downloaded structures.

    Example:
        >>> download("QE", 1, 10)
    """
    aflow_url = "http://www.aflowlib.org/API/aflux/"
    # Read mpid-list.in to get AFLOW IDs
    with open("mpid-list.in","r") as read_mpid:
        mpid_data = read_mpid.readlines()
    # Extract data for the specified range of indices
    mpid_data = mpid_data[start-1:end-1]
    list_data = []
    # Retrieve data from AFLOW for each AFLOW ID
    for idx,mpid in enumerate(mpid_data):
        aflow_id = mpid.split(" ")[1]
        subd = aflow_url + "?auid('{}'),geometry,positions_cartesian,positions_fractional,species".format(aflow_id)
        data = urlopen(subd).read().decode('utf-8')
        data = ast.literal_eval(data)
        list_data.append(poscar_create(data[0]))
    # Group similar structures
    unique_structures = StructureMatcher().group_structures(list_data)
    # Create input files and save them
    for index,subd in enumerate(unique_structures):
        try:
            subd = subd[0]
            aflow_id = subd.site_properties['auid'][0]
            aflow_id = "aflow-" + aflow_id.split(":")[1]
            # Write POSCAR file
            poscar = subd.to("POSCAR")
            with open("POSCAR","w") as write_poscar:
                write_poscar.write(poscar)
            compound = subd.composition.formula.replace(" ","")
            print(compound)
            # Create input files
            poscar_to_input(calc_type,aflow_id,compound,False)
            # Move POSCAR file if it exists
            if os.path.isfile("POSCAR"):
                os.system("mv POSCAR POSCAR-{}".format(aflow_id))
            # Write mpid.in file
            if not os.path.isfile('mpid.in'):
                k_ind = 0
                with open("mpid.in", "w") as write_mpid:
                    write_mpid.write("v{}".format(k_ind+1) + " " + aflow_id + " " + compound + "\n")
            else:
                with open('mpid.in', 'r') as read_mpid:
                    lines = read_mpid.readlines()
                k_ind = len(lines)
                if not any(aflow_id in line for line in lines):
                    with open("mpid.in", "a") as write_mpid:
                        write_mpid.write("v{}".format(k_ind+1) + " " + aflow_id + " " + compound + "\n")
        except:
            continue
def main():
    """
    Main program.

    This function serves as the entry point to the program. It reads command-line
    arguments to determine the operation mode (search or download) and then performs
    the corresponding action based on the provided inputs.
    """
    # Read command-line argument to determine operation mode
    condition = sys.argv[1]
    # Extract input data
    with open("input.in","r") as read_in:
        lines = read_in.readlines()
    start = int(lines[0].split("\n")[0])
    end = int(lines[1].split("\n")[0])
    #start = input_data['inp']['start']
    #end = input_data['inp']['end']
    dft = input_data['inp']['calc']
    # Call search_data function if the mode is search
    if condition == 'search':
        search_data(input_data)
    elif condition == 'download':
        print(start,end,dft)
        # Call download function if the mode is download
        download(dft,start,end)
    else:
        print("Either search or download is allowed\n")
if __name__ == "__main__":
    main()
