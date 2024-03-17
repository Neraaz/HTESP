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
from ase.cell import Cell
from pymatgen.core.structure import IStructure
from pymatgen.core.composition import Composition
from pymatgen.analysis.structure_matcher import StructureMatcher
from oqmd_extract import poscar_to_input
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
        input_data = json.load(readjson)['download']
except FileNotFoundError:
    print("htepc.json file not found\n")
def properties_string(dict_param):
    """
    Function to create a string for URL search from the dictionary provided in htepc.json.

    Parameters:
    dict_param (dict): A dictionary from the htepc.json file containing search criteria.

    Returns:
    str: A string representing the search criteria for URL.
    """
    criteria = ""
    search = dict_param['aflow']
    apply_filter = search['filter']
    print(search,apply_filter)
    if apply_filter:
        for key in search:
            if key == 'elm':
                if len(tuple(search[key])) == 1:
                    criteria += "species({}),".format(tuple(search[key])[0])
                else:
                    criteria += "species{},".format(tuple(search[key]))
            elif key == 'nelm':
                criteria += "nspecies(1*,*{}),".format(search[key])
            elif key == 'nsites':
                criteria += "natoms(1*,*{}),".format(search[key])
            elif key == 'metal':
                if search[key]:
                    criteria += "Egap(0),"
                else:
                    criteria += "Egap(*),"
            elif key == 'FE':
                if search[key]:
                    criteria += "enthalpy_formation_atom(-100*,*0),"
                else:
                    criteria += "enthalpy_formation_atom(*),"
            elif key == 'spacegroup':
                if search[key] is None:
                    criteria += "spacegroup_relax(*),"
                else:
                    criteria += "spacegroup_relax({}),".format(search[key])
            elif key == 'limit':
                criteria += "$paging(1,{}),".format(search[key])
            elif key == 'prop':
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
    criteria = criteria.replace(" ","")
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
    """
    aflow_url = "http://www.aflowlib.org/API/aflux/"
    criteria = properties_string(dict_param)
    #print(criteria)
    search_url = aflow_url + "?" + criteria
    success = False
    while not success:
        try:
            data = urlopen(search_url).read().decode('utf-8')
            success = True
        except:
            success = False
            continue
    data_processed = ast.literal_eval(data)
    data_processed = pd.DataFrame(data_processed)
    #print(data_processed)
    if not os.path.isdir("download"):
        os.mkdir("download")
    if not os.path.isfile("download/download-aflow.csv"):
        data_processed.to_csv("download/download-aflow.csv",index=False)
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
    """
    lattice = Cell.fromcellpar(data['geometry']).array
    comp = Composition(data['compound'])
    aflow_id = data['auid']
    species = []
    comp = comp.as_dict()
    for key in comp:
        for i in range(int(comp[key])):
            species.append(key)
    struc = IStructure(lattice,species=species,coords=data['positions_fractional'],to_unit_cell=True)
    site_prop = [aflow_id] * len(struc.sites)
    struc.add_site_property("auid",site_prop)
    #poscar = struc.to("POSCAR")
    #with open("POSCAR","w") as write_poscar:
    #    write_poscar.write(poscar)
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
    """
    aflow_url = "http://www.aflowlib.org/API/aflux/"
    with open("mpid-list.in","r") as read_mpid:
        mpid_data = read_mpid.readlines()
    mpid_data = mpid_data[start-1:end-1]
    list_data = []
    for idx,mpid in enumerate(mpid_data):
        aflow_id = mpid.split(" ")[1]
        subd = aflow_url + "?auid('{}'),geometry,positions_cartesian,positions_fractional,species".format(aflow_id)
        data = urlopen(subd).read().decode('utf-8')
        data = ast.literal_eval(data)
        list_data.append(poscar_create(data[0]))
    unique_structures = StructureMatcher().group_structures(list_data)
    for index,subd in enumerate(unique_structures):
        try:
            subd = subd[0]
            #print(subd)
            #aflow_id = subd['id']
            aflow_id = subd.site_properties['auid'][0]
            aflow_id = "aflow-" + aflow_id.split(":")[1]
            #aflow_id = "aflow-{}".format(index+1)
            #oqmd_id = "aflow-" + str(oqmd_id)
            #print(oqmd_id + "\n")
            poscar = subd.to("POSCAR")
            with open("POSCAR","w") as write_poscar:
                write_poscar.write(poscar)
            compound = subd.composition.formula.replace(" ","")
            print(compound)
            poscar_to_input(calc_type,aflow_id,compound,False)
            if os.path.isfile("POSCAR"):
                os.system("mv POSCAR POSCAR-{}".format(aflow_id))
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
    condition = sys.argv[1]
    start = input_data['inp']['start']
    end = input_data['inp']['end']
    dft = input_data['inp']['calc']
    if condition == 'search':
        search_data(input_data)
    elif condition == 'download':
        print(start,end,dft)
        download(dft,start,end)
    else:
        print("Either search or download is allowed\n")
if __name__ == "__main__":
    main()
