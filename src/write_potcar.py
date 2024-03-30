#!/usr/bin/env python
# coding: utf-8
# Written by Niraj K. Nepal, Ph.D.
"""
Writing POTCAR file from POSCAR
"""
import os
import json
from pymatgen.core import structure
from pymatgen.io.vasp.sets import Potcar
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
def poscar2potcar():
    """
    Function to generate the POTCAR file based on the element potentials provided in
    config.json file.

    If config.json file is present in the current directory or its parent directory, the function retrieves
    the potential dictionary from the 'pseudo' key in the config.json file.

    If the config.json file is not found, a message indicating the absence of the 'pseudo' keyword is displayed.

    The function reads the POSCAR file to determine the elements present in the structure and generates the
    corresponding POTCAR file based on the potential dictionary.

    """
    if os.path.isfile("config.json") or os.path.isfile("../../config.json"):
        pot1 = input_data['pseudo']
    else:
        print("pseudo keyword not found\n")
    struc = structure.Structure.from_file("POSCAR")
    pot2 = pot1['pot']
    element_list = []
    for elm in struc.composition.elements:
        element_list.append(pot2[str(elm)])
    potcar = Potcar(element_list)
    potcar.write_file('POTCAR')
