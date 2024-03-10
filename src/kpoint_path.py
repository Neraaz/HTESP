#!/usr/bin/env python
"""Written by Niraj K. Nepal, Ph.D"""
from ase.io import espresso, vasp
from ase.cell import Cell
from kpath import kpath
def kpoint_path(file_name):
    """
    Function to print Kpoint_path section in wannier input.
    parameters
    ----------------
    file_name : QE scf.in file
    """
    _,_,_,_,sym,_ = kpath(file_name,1,0)
    try:
        data = espresso.read_espresso_in(file_name)
    except:
        data = vasp.read_vasp(file_name)
    band = Cell.bandpath(data.cell)
    band_dict = band.special_points
    key = []
    value = []
    for i,subsym in enumerate(sym):
        key.append(subsym)
        value.append(band_dict[subsym])
    with open('wannier_kpath.in', 'w') as wan_kpath_write:
        for i,_ in enumerate(key):
            if i < len(key) - 1:
                wan_kpath_write.write("{} {} {} {} \t".format(key[i],round(value[i][0],5),round(value[i][1],5),round(value[i][2],5)))
                wan_kpath_write.write("{} {} {} {}".format(key[i+1],round(value[i+1][0],5),round(value[i+1][1],5),round(value[i+1][2],5)) + "\n")
