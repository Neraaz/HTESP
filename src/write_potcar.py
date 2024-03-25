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
    if os.path.isfile(PWD+"/htepc.json"):
        JSONFILE = PWD+"/htepc.json"
    else:
        JSONFILE = "../../htepc.json"
    with open(JSONFILE, "r") as readjson:
        input_data = json.load(readjson)
except FileNotFoundError:
    print("htepc.json file not found\n")
def poscar2potcar():
    """
    Function to generate the POTCAR file based on the element potentials provided in
    htepc.json file.

    If htepc.json file is present in the current directory or its parent directory, the function retrieves
    the potential dictionary from the 'pseudo' key in the htepc.json file.

    If the htepc.json file is not found, a message indicating the absence of the 'pseudo' keyword is displayed.

    The function reads the POSCAR file to determine the elements present in the structure and generates the
    corresponding POTCAR file based on the potential dictionary.

    """
    if os.path.isfile("htepc.json") or os.path.isfile("../../htepc.json"):
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
