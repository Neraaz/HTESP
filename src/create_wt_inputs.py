#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to write wanniertools input files"""
import sys
import os
import json
import numpy as np
import pandas as pd
from ase.io import espresso,vasp
from ase.cell import Cell
from kpath import kpath

try:
    PWD = os.getcwd()
    if os.path.isfile(PWD+"/htepc.json"):
        JSONFILE = PWD+"/htepc.json"
    else:
        JSONFILE = "../../htepc.json"
    with open(JSONFILE, "r") as readjson:
        wt_input = json.load(readjson)['wanniertool_input']
except FileNotFoundError:
    print("htepc.json file not found, writting wanniertool_input.json file with default values\n")
    tb_file={'Hrfile':"'ex_hr.dat'", "Package":"'QE'"}
    control={"BulkBand_calc":"T","BulkFS_calc":"F","BulkGap_cube_calc":"F","BulkGap_plane_calc":"F","FindNodes_calc":"F","SlabBand_calc":"F","WireBand_calc":"F","Dos_calc":"F","JDos_calc":"F","SlabSS_calc":"F","SlabArc_calc":"F","SlabQPI_calc":"F","SlabSpintexture_calc":"F","wanniercenter_calc":"F","Z2_3D_calc":"F","Chern_3D_calc":"F","WeylChirality_calc":"F","BerryPhase_calc":"F","BerryCurvature_calc":"F","AHC_calc":"F"}
    parameters={"E_arc":"0.0","Eta_Arc":"0.001","OmegaMin":"0.0","OmegaMax":"0.0","OmegaNum":"100","Nk1":"10","Nk2":"10","Nk3":"10","NP":"2","Gap_threshold":"0.1"}
    system={"NSlab":"10","NSlab1":"1","NSlab2":"1","NumOccupied":"1","SOC":"1","E_FERMI":"0.0","Bx":"0.0","By":"0.0","Bz":"0.0","surf_onsite":"0.0"}
    surface_dict={"surface":[[1,0,0],[0,1,0],[0,0,1]],"KPATH_SLAB":{},"KPLANE_SLAB":{},"EFFECTIVE_MASS":"0.0","SELECTED_ATOMS":{}}
    wt_input = {
                'tb_file':tb_file,
                'control':control,
                'parameters':parameters,
                'system':system,
                'surface':surface_dict}
#try:
#    import wanniertool_input as wt_input
#except FileNotFoundError:
#    print("wanniertool_input.py not found, creating a default one\n")
#    with open("wanniertool_input.py", "w") as f:
#        f.write("tb_file={}".format(tb_file) + "\n")
#        f.write("control={}".format(control) + "\n")
#        f.write("parameters={}".format(parameters) + "\n")
#        f.write("system={}".format(system) + "\n")
#        f.write("surface={}".format(surface_dict) + "\n")
#    import wanniertool_input as wt_input

def kpoint_path(filename):
    """
    Extracts Brillouin zone high symmetry path using ASE modules for bulk and slab systems
    (if POSCAR-slab is available).

    Parameters:
    - filename (str): The filename of the QE input file 'scf.in' or VASP 'POSCAR' file.

    Returns: None
    """
    _,_,_,_,sym,_ = kpath(filename,1,0)
    data = espresso.read_espresso_in(filename)
    band = Cell.bandpath(data.cell)
    band_dict = band.special_points
    key = []
    value = []
    for i,sym_i in enumerate(sym):
        key.append(sym_i)
        value.append(band_dict[sym_i])
    with open('wannier_kpath.in', 'w') as wtfile_write:
        for i,_ in enumerate(key):
            if i < len(key) - 1:
                wtfile_write.write("{} {} {} {}  ".format(key[i],round(value[i][0],5),round(value[i][1],5),round(value[i][2],5)))
                wtfile_write.write("{} {} {} {}".format(key[i+1],round(value[i+1][0],5),round(value[i+1][1],5),round(value[i+1][2],5)) + "\n")
    if os.path.isfile("POSCAR-slab"):
        slab = vasp.read_vasp("POSCAR-slab")
        band_s = Cell.bandpath(slab.cell,pbc=[1,1,0])
        bands_dict = band_s.special_points
        sym_s = band_s.get_linear_kpoint_axis()[2]
        key_s = []
        value_s = []
        for i,sym_i in enumerate(sym_s):
            key_s.append(sym_i)
            value_s.append(bands_dict[sym_i])
        with open('wannier_kpath_slab.in', 'w') as wtfile_write:
            for i,_ in enumerate(key_s):
                if i < len(key_s) - 1:
                    wtfile_write.write("{} {} {}  ".format(key_s[i],round(value_s[i][0],5),round(value_s[i][1],5)))
                    wtfile_write.write("{} {} {}".format(key_s[i+1],round(value_s[i+1][0],5),round(value_s[i+1][1],5)) + "\n")
def wt_body(mpid,compound,surface=False):
    """
    Writes wt.in file in the format of 'wt-mpid-compound.in' inside WT_dir, based on information
    given in wanniertool_input keyword within htepc.json file.

    Parameters:
    - mpid (str): Materials ID.
    - compound (str): Compound name used in mpid-list.in.
    - surface (bool): If True, includes cards printed for surface calculations. Default is False.

    Returns: None
    """
    if os.path.isfile("R{}-{}/relax/scf.in".format(mpid,compound)):
        input_file="R{}-{}/relax/scf.in".format(mpid,compound)
        data = espresso.read_espresso_in(input_file)
    elif os.path.isfile("R{}-{}/relax/POSCAR".format(mpid,compound)):
        input_file="R{}-{}/relax/POSCAR".format(mpid,compound)
        data = vasp.read_vasp(input_file)
    else:
        print("No input file found, either vasp or QE" + "\n")
    os.system("sed -n '/Final State/,/Sum of centres and spreads/p'" + " R{}-{}/epw/ex.wout".format(mpid,compound) + "| awk '{ print $7 $8 $9}' > wannier.csv")
    os.system("sed -i '1s/^/x,y,z/' wannier.csv")
    os.system("sed -i '$d' wannier.csv")
    #l=data.cell.get_bravais_lattice()
    cell=data.cell
    symbol=data.get_chemical_symbols()
    pos=data.get_scaled_positions()
    dict_orb = {'s':'s', 'p': 'pz px py', 'd': 'dz2 dxz dyz dx2-y2 dxy'}
    with open("wt-{}-{}.in".format(mpid,compound), "w") as wtfile_write:
        wtfile_write.write("&TB_FILE\n")
        keys = list(wt_input['tb_file'].keys())
        values = list(wt_input['tb_file'].values())
        for tb_i,_ in enumerate(wt_input['tb_file']):
            wtfile_write.write(keys[tb_i] + "=" + values[tb_i] + "\n")
        wtfile_write.write("/\n")
        wtfile_write.write("\n")
        wtfile_write.write("&CONTROL\n")
        keys = list(wt_input['control'].keys())
        values = list(wt_input['control'].values())
        for control_i,_ in enumerate(wt_input['control']):
            wtfile_write.write(keys[control_i] + "=" + values[control_i] + "\n")
        wtfile_write.write("/\n")
        wtfile_write.write("\n")
        wtfile_write.write("&SYSTEM\n")
        keys = list(wt_input['system'].keys())
        values = list(wt_input['system'].values())
        for system_i,_ in enumerate(wt_input['system']):
            wtfile_write.write(keys[system_i] + "=" + values[system_i] + "\n")
        wtfile_write.write("/\n")
        wtfile_write.write("\n")
        wtfile_write.write("&PARAMETERS\n")
        keys = list(wt_input['parameters'].keys())
        values = list(wt_input['parameters'].values())
        for param_i,_ in enumerate(wt_input['parameters']):
            wtfile_write.write(keys[param_i] + "=" + values[param_i] + "\n")
        wtfile_write.write("/\n")
        wtfile_write.write("\n")
        wtfile_write.write("LATTICE\n")
        wtfile_write.write("Angstrom\n")
        for lat in cell:
            wtfile_write.write(str(lat[0]) + " " + str(lat[1]) + " " + str(lat[2]) + "\n")
        wtfile_write.write("\n")
        wtfile_write.write("ATOM_POSITIONS\n")
        wtfile_write.write(str(len(symbol)) + "\n")
        wtfile_write.write('Direct\n')
        for s_i,pos_i in enumerate(pos):
            wtfile_write.write(symbol[s_i] + " " + str(np.round(pos_i[0],10)) + " ")
            wtfile_write.write(str(np.round(pos_i[1],10)) + " " + str(np.round(pos_i[2],10)) + "\n")
        wtfile_write.write("\n")
        wtfile_write.write("PROJECTORS\n")
        if os.path.isfile('proj-wt.in'):
            print("proj-wt.in file found, Reading .....")
            with open("proj-wt.in", "r") as gfile:
                lines = gfile.readlines()
            dict_elm = {}
            for line in lines:
                dict_elm[line.split('\n')[0].split()[0]] = line.split('\n')[0].split()[1:]
            for sym in symbol:
                norb = 0
                sorb = dict_elm[sym]
                for orb in sorb:
                    if orb in dict_orb.keys():
                        norb += len(dict_orb[orb].split())
                    else:
                        norb += 1
                wtfile_write.write(str(norb) + " ")
            wtfile_write.write("!number of projectors\n")
            for sym in symbol:
                orb_str = ''
                sorb = dict_elm[sym]
                for orb in sorb:
                    if orb in dict_orb.keys():
                        orb_str = orb_str + " " +  dict_orb[orb]
                    else:
                        orb_str = orb_str + " " + orb
                wtfile_write.write(sym + " " + orb_str + "\n")
        else:
            print("proj-wt.in file not found, writing only s orbital .....")
            for sym in symbol:
                wtfile_write.write("1 ")
            wtfile_write.write("!number of projectors\n")
            for sym in symbol:
                wtfile_write.write(sym + " s" + "\n")
        wtfile_write.write("\n")
        if surface:
            wtfile_write.write("SURFACE\n")
            if os.path.isfile("htepc.json") or os.path.isfile("../../htepc.json"):
                print("Taking surface from htepc.json\n")
                for surf_lat in wt_input['surface']['surface']:
                    wtfile_write.write(str(surf_lat[0]) + " " + str(surf_lat[1]) + " " + str(surf_lat[2]) + "\n")
            elif os.path.isfile("surface-wt.in"):
                print("surface-wt.in file found")
                with open('surface-wt.in', 'r') as gfile:
                    lines = gfile.readlines()
                for line in lines:
                    wtfile_write.write(line)
            else:
                print("surface-wt.in and surface inside wanniertoo_input.py file not found, printing default values")
                wtfile_write.write("1 0 0 ! write vectors determining surfaces\n")
                wtfile_write.write("0 1 0 \n")
        wtfile_write.write("\n")
        wtfile_write.write("KPATH_BULK\n")
        with open("R{}-{}/epw/ex.win".format(mpid,compound), "r") as hfile:
            lines = hfile.readlines()
        if 'Begin Kpoint_Path\n' in lines and 'End Kpoint_Path\n' in lines:
            os.system("sed -n '/Begin Kpoint_Path/,/End Kpoint_Path/p' R{}-{}/epw/ex.win | sed '1d' | sed '$d' > wannier_kpath.in".format(mpid,compound))
        with open("wannier_kpath.in", "r") as gfile:
            lines = gfile.readlines()
        wtfile_write.write(str(len(lines)) + "\n")
        for i,line in enumerate(lines):
            wtfile_write.write(line)
        wtfile_write.write("\n")
        if surface:
            if os.path.isfile("wannier_kpath_slab.in"):
                with open("wannier_kpath_slab.in", "r") as wanpath:
                    kpath_2d = wanpath.readlines()
                wtfile_write.write("KPATH_SLAB\n")
                wtfile_write.write(str(len(kpath_2d)) + "\n")
                for line in kpath_2d:
                    wtfile_write.write(line)
                wtfile_write.write("\n")
        wtfile_write.write("WANNIER_CENTRES\n")
        wtfile_write.write("Cartesian\n")
        wannier = pd.read_csv("wannier.csv")
        for i in range(wannier.shape[0]):
            wtfile_write.write(str(wannier.x[i]) + " " + str(wannier.y[i]) + " " + str(wannier.z[i]) + "\n")
def main():
    """
    Main function to execute different operations based on command-line arguments.

    Command-line arguments:
    - mpid (str): Material ID.
    - compound (str): Compound name.
    - prefix (str): Prefix value.
    - condition (str): Condition for different operations.
    - surface_prop (str): Surface property indicator.

    Returns: None
    """
    mpid = sys.argv[1]
    compound = sys.argv[2]
    prefix = sys.argv[3]
    condition = sys.argv[4]
    if condition == "kpathwan":
        if os.path.isfile("scf_dir/scf-{}-{}.in".format(mpid,compound)):
            file_name="scf_dir/scf-{}-{}.in".format(mpid,compound)
        else:
            print("scf-{}-{}.in not found in scf_dir/".format(mpid,compound) + "\n")
            print("Complete mainprogram.py from 1 to 4\n")
        kpoint_path(file_name)
    elif condition == "body":
        surface_prop = sys.argv[5]
        if surface_prop == 'T':
            wt_body(mpid,compound,surface=True)
        else:
            wt_body(mpid,compound,surface=False)
    else:
        print("bad inputs")
if __name__ == "__main__":
    main()
