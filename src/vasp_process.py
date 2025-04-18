#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to process vasp input files"""
import sys
import os
import numpy as np
from ase.io import vasp,espresso
from ase.cell import Cell
from pymatgen.core import structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.sets import Incar
from cif_to_gsinput import pos_to_kpt
from check_json import config
def encut_check():
    """
    Function to check ENCUT in INCAR file and replace
    with 1.3xENCUTMAX of POTCAR, when ENCUT is less than that.
    """
    # Obtain ENMAX from POTCAR
    os.system("""grep ENMAX POTCAR | awk '{ print $3 }' | sed "s/;//" > encut.txt""")
    with open("encut.txt", "r") as read_encut:
        lines = read_encut.readlines()
    encut = []
    for line in lines:
        encut.append(float(line.split("\n")[0]))
    # Calculate 1.3 times of maximum ENCUT in found in POTCAR
    encutmax = max(encut)*1.3
    # Get ENCUT from current INCAR, set 0 if not found
    try:
        os.system("grep ENCUT INCAR | awk '{ print $3 }' > encut")
        with open("encut", "r") as read_ecut:
            lines = read_ecut.readlines()
        encut = float(lines[0].split("\n")[0])
    except:
        print("No ENCUT parameters found in INCAR or vasp.in\n")
        print("1.3 times ENMAX is used from POTCAR\n")
        encut = 0.0
    # If ENCUT in INCAR is less than encutmax, update it with ecutmax
    if encut < encutmax:
        print("------------------------------------------------------\n")
        print("ENCUT less than 1.3 times ENCUTMAX in POTCAR found\n")
        print("Adjusting ENCUT in INCAR .........................\n")
        print("------------------------------------------------------\n")
        os.system("sed -i '/ENCUT/d' INCAR")
        encut = encutmax
        os.system("echo ENCUT = {} >> INCAR".format(encut))
    os.system("rm encut encut.txt")
def vasp_process():
    """
    Processes the INCAR file using the 'vasp.in' file.
    Look for 'vasp.in' file in the utility folder as a demo.

    If 'vasp.in' is found, it reads the file and replaces corresponding keys in the 'INCAR' file with the provided values.
    Existing keys such as ISPIN, MAGMOM, and LORBIT in 'INCAR' are removed to avoid conflicts.
    If magnetic moment values are provided in 'config.json', it updates the 'MAGMOM' keyword accordingly.
    If 'METAGGA' is present, 'LMIXTAU = .TRUE.' is added to 'INCAR'.
    To remove unwanted keywords from 'INCAR', put those keywords at the bottom in 'vasp.in' file after the keyword that has values.

    Returns:
        None
    """
    input_data = config()
    # Get the type of magnetic enumeration
    magenum = input_data['magmom']['type']
    # Check if 'vasp.in' file exists
    if os.path.isfile("vasp.in"):
        if magenum == 'anisotropy':
            print("anisotropy type found in magmom dictionary\n")
            # If magnetic enumeration is anisotropy, remove NSW from 'vasp.in'
            os.system("""sed -i '/NSW/d' vasp.in""")
        # Read 'vasp.in' file
        with open("vasp.in", 'r') as read_vasp:
            lines = read_vasp.readlines()
        # Extract keys and values from 'vasp.in'
        key = []
        values = []
        for line in lines:
            key_element = line.split('\n')[0].split(' ')[0]
            if key_element != '':
                key.append(key_element)
            value_list = line.split('\n')[0].split(' ')[1:]
            value_list = list(filter(None,value_list))
            k_str = ''
            for ii,line_i in enumerate(value_list):
                if ii < len(value_list) - 1:
                    k_str += line_i + " "
                else:
                    k_str += line_i
            if k_str != '':
                values.append(k_str)
        # Remove all the keys in INCAR which is found in vasp.in
        for k_i in key:
            os.system("""sed -i '/{}/d' INCAR""".format(k_i))
        with open("INCAR", "r") as read_incar:
            newincar = read_incar.readlines()
        backupkey = {'ISPIN':1,'MAGMOM':False,'LORBIT':False}
        for _,oldkey in enumerate(newincar):
            keyx = oldkey.split(" ")[0]
            if keyx == 'ISPIN':
                backupkey[keyx] = int(oldkey.split(" ")[2].split("\n")[0])
            elif keyx == 'LORBIT':
                backupkey['LORBIT'] = True
            elif keyx == 'MAGMOM':
                backupkey['MAGMOM'] = True
            else:
                backupkey = backupkey
        # Remove MAGMOM keyword if exists and magenum is not 'anisotropy'
        if backupkey['MAGMOM'] and magenum != 'anisotropy':
            os.system("""sed -i '/MAGMOM/d' INCAR""")
        # Remove NSW keyword if magenum is 'anisotropy'
        if magenum == 'anisotropy':
            os.system("""sed -i '/NSW/d' INCAR""")
        if backupkey['LORBIT']:
            os.system("""sed -i '/LORBIT/d' INCAR""")
        # Write new keys and values to 'INCAR'
        metagga = False
        spinval = 1
        with open("INCAR", "a") as change_incar:
            for i,_ in enumerate(values):
                if key[i] == "ISPIN":
                    spinval = values[i]
                if key[i] == "METAGGA":
                    metagga = True
                try:
                    change_incar.write(key[i] + " = " + str(values[i]) + "\n")
                except IndexError:
                    continue
            # If ISPIN=2, update MAGMOM keyword if config.json exists
            if ("ISPIN" in key and int(spinval) == 2) or backupkey['ISPIN'] == 2:
                if os.path.isfile("config.json") or os.path.isfile("../../config.json"):
                    m=input_data['magmom']['magmom']
                    struc = structure.Structure.from_file("POSCAR")
                    sites = struc.sites
                    magmom_string = ""
                    for j,site in enumerate(sites):
                        element = str(site.specie)
                        if j < len(sites) - 1:
                            # if LSORBIT key found, rewrite MAGMOM in mx my mz format
                            if 'LSORBIT' not in key:
                                magmom_string += str(m[element]) + " "
                            else:
                                magmom_string += "0 0 " + str(m[element]) + "  "
                        else:
                            if 'LSORBIT' not in key:
                                magmom_string += str(m[element]) + "\n"
                            else:
                                magmom_string += "0 0 " + str(m[element]) + "\n"
                    if magenum != 'anisotropy':
                        change_incar.write("MAGMOM = {}".format(magmom_string))
                        change_incar.write("LORBIT = 11\n")
                    else:
                        print("type is anisotropy, therefore doesn't update MAGMOM keyword\n")
                else:
                    print("magmom values not provided in config.json\n")
                    print("Provide magnetic moment values as dictionary magmom={'A':2, 'B':3}\n")
            # Add LMIXTAU = .TRUE. if METAGGA is present in INCAR
            if metagga:
                change_incar.write("LMIXTAU = .TRUE.\n")
                change_incar.write("LASPH = .TRUE.\n")
    # Check and adjust ENCUT
    encut_check()
def eigen_process():
    """
    Process the EIGENVAL file generated by VASP.

    Reads the 'input.in' file to determine if the 'vasp-line' keyword is present, indicating
    a line-mode calculation.

    Reads the 'EIGENVAL' file to extract band structure data and writes it to 'band.dat'.

    If 'vasp-line' is present or the weight of a k-point is less than 0.000001, it writes the k-point index
    and corresponding band energies to 'band.dat'. Otherwise, it skips the k-point and its associated data.

    Additionally, it creates 'band.dat.gnu', a file suitable for plotting band structures in GNUPlot.

    Returns:
        None
    """
    with open("../../input.in","r") as read_inputin:
        inputline = read_inputin.readlines()
    vasp_line = False
    for line in inputline:
        if "vasp-line" in line:
            vasp_line = True
    if os.path.isfile("KPT_OPT"):
        vasp_line = True
    # TO DO: parse EIGENVAL for spin-polarized calculations
    incar = Incar.from_file("INCAR")
    if 'LSORBIT' not in incar.keys():
        incar['LSORBIT'] = False
    if 'ISPIN' not in incar.keys():
        incar['ISPIN'] = 1
    if 'LNONCOLLINEAR' not in incar.keys():
        incar['LNONCOLLINEAR'] = False
    if incar['LSORBIT'] or incar['LNONCOLLINEAR']:
        nspin = 4
    elif incar['ISPIN'] == 2 and not incar['LNONCOLLINEAR']:
        nspin = 2
    else:
        nspin = 1
    with open("EIGENVAL", "r") as read_eig:
        for i in range(5):
            read_eig.readline()
        _, nkpoints, nbands = [int(eig) for eig in read_eig.readline().split()]
        read_eig.readline()
        with open('band.dat', 'w') as write_eig:
            k_ind = 1
            for i in range(nkpoints):
                _, _, _, weight = [float(eig) for eig in read_eig.readline().split()]
                if weight < 0.000001 or vasp_line:
                    write_eig.write(str(k_ind) + " ")
                if weight < 0.000001 or vasp_line:
                    for j in range(nbands):
                        fields = read_eig.readline().split()
                        if j < nbands - 1:
                            if nspin == 2:
                                write_eig.write(str(fields[1]) + " " + str(fields[2]) + " ")
                            else:
                                write_eig.write(str(fields[1]) + " ")
                        else:
                            if nspin == 2:
                                write_eig.write(str(fields[1]) + " " + str(fields[2]) + "\n")
                            else:
                                write_eig.write(str(fields[1]) + "\n")
                    read_eig.readline()
                    k_ind += 1
                else:
                    for j in range(nbands):
                        fields = read_eig.readline().split()
                    read_eig.readline()
    with open("band.dat.gnu", "w") as write_dat_gnu:
        data = np.loadtxt('band.dat')
        nrow,_ = data.shape
        if nspin == 2:
            band_value = data[:,1:]
            band_value_1 = band_value[:, ::2]
            band_value_2 = band_value[:, 1::2]
            band_value = band_value_1
        else:
            band_value = data[:,1:]
        for i in range(band_value.shape[1]):
            for j in range(band_value.shape[0]):
                if j < nrow - 1:
                    if nspin == 2:
                        write_dat_gnu.write(str(j) + " " + str(band_value_1[j,i]) + " " + str(band_value_2[j,i]) + "\n")
                    else:
                        write_dat_gnu.write(str(j) + " " + str(band_value[j,i]) + "\n")
                else:
                    if nspin == 2:
                        write_dat_gnu.write(str(j) + " " + str(band_value_1[j,i]) + " " + str(band_value_2[j,i]) + "\n")
                        write_dat_gnu.write("\n")
                    else:
                        write_dat_gnu.write(str(j) + " " + str(band_value[j,i]) + "\n")
                        write_dat_gnu.write("\n")
def band_phonopy(filename):
    """
    Write k-points of high-symmetry points for Phonopy band structure plot using 'band.conf' file.

    Parameters:
    -----------
    filename : str
        The name of the file containing the crystal structure, either 'POSCAR' for VASP or 'scf.in' for Espresso.

    Returns:
    --------
    None
    """
    if filename == 'POSCAR':
        data = vasp.read_vasp('POSCAR')
    else:
        data = espresso.read_espresso_in('scf.in')
    cell_bp = Cell.bandpath(data.cell)
    res = cell_bp.path.replace(',','')
    #print(res)
    new_list = []
    ind = 0
    for _,s in enumerate(list(res)):
        if not s.isdigit():
            new_list.append(s)
            temp  = s
            index = ind
            ind += 1
        else:
            s1 = temp + str(s)
            new_list[index] = s1

    #res = "".join(filter(lambda x: not x.isdigit(), res))
    npt = len(new_list)
    #print(new_list,npt)
    band = data.cell.bandpath(path=cell_bp.path,npoints=npt)
    specialpt = band.special_points
    with open("band_phonopy.in", "w") as band_phon:
        band_phon.write("BAND = ")
        #for i in range(band.shape[0]):
        #    band_phon.write(str(band[i][0]) + " " + str(band[i][1]) + " " + str(band[i][2]) + "\t")
        for _,special in enumerate(new_list):
            band_phon.write(str(specialpt[special][0]) + " " + str(specialpt[special][1]) + " " + str(specialpt[special][2]) + "\t")
        band_phon.write("\n")
        band_phon.write("BAND_LABELS = ")
        for _,special in enumerate(new_list):
            band_phon.write(special + "\t")
    with open("high_symm.in", "w") as high_sym:
        special = list(Cell.bandpath(data.cell).special_points.values())
        for j,_ in enumerate(special):
            high_sym.write(str(special[j][0]) + " ")
            high_sym.write(str(special[j][1]) + " " + str(special[j][2]) + "\n")
def main():
    """
    main function
    """
    input_data = config()
    filename = sys.argv[1]
    kptden = input_data["kptden"]
    if filename == 'POSCAR':
        band_phonopy(filename)
        vasp_process()
        # remove EIGENVAL if you only update INCAR and KPOINTS
        if os.path.isfile("EIGENVAL"):
            os.system("cp KPOINTS KPOINTS_band")
            eigen_process()
        pos_to_kpt(filename,kptden)
    # Change POSCAR to conventional unit cell
    # and update KPOINTS
    elif filename == 'conventional':
        data = structure.Structure.from_file("POSCAR")
        data = SpacegroupAnalyzer(data,symprec=0.1).get_conventional_standard_structure()
        data.to("POSCAR")
        pos_to_kpt("POSCAR",kptden)
    elif filename == 'eigen':
        os.system("mv KPOINTS_band KPOINTS")
    # Symmetrize the POSCAR
    elif filename == 'symmetrize':
        print("Symmetrizing the primitive structure\n")
        data = structure.Structure.from_file("POSCAR")
        data = SpacegroupAnalyzer(data,symprec=0.1).get_primitive_standard_structure()
        data.to("POSCAR")
        pos_to_kpt("POSCAR",kptden)
    else:
        band_phonopy(filename)
if __name__ == "__main__":
    main()
