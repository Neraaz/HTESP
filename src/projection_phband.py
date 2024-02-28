#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D."""
import os
import sys
import numpy as np
from ase.io import espresso
np.set_printoptions(threshold=100000,linewidth=1000,suppress=True, formatter={'float_kind':'{:0.2f}'.format})
MSG="""
Format of phonproj.in file
_____________________________________________
plot.eig    #Eigenval file
2 111       #number of atom type and number of kpoints.
*
Mg 0        # first ion and its position in input file
O 1         # second ion and its position in input file
*
_____________________________________________
"""

def write_phonproj_in(comp,nkpoint,filename,outfile='phonproj.in'):
    """
    Function to write the 'phonproj.in' file required for atomic projection.

    This function generates the 'phonproj.in' file based on the provided compound name,
    the number of k-points, and the QE input scf.in file.

    Parameters:
    -----------
    comp : str
        Compound name.
    nkpoint : int
        Number of k-points.
    filename : str
        Path to the QE input scf.in file.
    outfile : str, optional
        Output file name for 'phonproj.in'. Default is 'phonproj.in'.

    Returns:
    --------
    None
    """
    struc = espresso.read_espresso_in(filename)
    symbol_dict = struc.symbols.indices()
    with open(outfile, "w") as write_phonproj:
        write_phonproj.write("{}.eig".format(comp)+ "\n")
        nat = len(symbol_dict.keys())
        write_phonproj.write("{} {}".format(nat,nkpoint) + "\n")
        write_phonproj.write("*\n")
        for _,elm in enumerate(symbol_dict.keys()):
            value = symbol_dict[elm]
            write_phonproj.write(elm + " ")
            for j,ind in enumerate(value):
                if j < len(value) - 1:
                    write_phonproj.write(str(ind) + " ")
                else:
                    write_phonproj.write(str(ind) + "\n")
        write_phonproj.write("*\n")
def reading_input(filename):
    """
    Function to read the 'phonproj.in' file.

    This function reads the 'phonproj.in' file and extracts relevant information
    such as the name of the .eig file, the number of atoms in the compound,
    the size of the q-point mesh, and a dictionary containing ions and their indices
    defining their position in the QE input file 'scf.in'.

    Parameters:
    -----------
    filename : str
        Path to the 'phonproj.in' file.

    Returns:
    --------
    eig_file : str
        Name of the .eig file.
    nat : int
        Number of atoms in the compound.
    nqpt : int
        Size of the q-point mesh (nq, 3).
    ions : dict
        Dictionary with ions and indices defining their position in the 'scf.in' file.
    """
    with open(filename, "r") as proj_file:
        lines = proj_file.readlines()
    ions = {}
    ind = []
    for i,line in enumerate(lines):
        if "*" in line:
            ind.append(i)
    ion_list = lines[ind[0]+1:ind[1]]
    na_t = len(ion_list)
    ions = {}
    nat = 0
    for i in range(na_t):
        if "#" in ion_list[i]:
            temp = ion_list[i].split("#")[0].split(" ")
        else:
            temp = ion_list[i].split("\n")[0].split(" ")
        name = temp[0]
        indices = [int(j) for j in temp[1:] if j != ""]
        ions[name] = indices
        nat += len(indices)
    for line in lines:
        if "#" in line:
            eig_file = lines[0].split("#")[0].split(" ")[0]
            nqpt = int(lines[1].split("#")[0].split(" ")[1])
        else:
            eig_file = lines[0].split("\n")[0].split(" ")[0]
            nqpt = int(lines[1].split("\n")[0].split(" ")[1])
    return eig_file, nat, nqpt, ions

def parse_input(filename,nat,nqpt):
    """
    Function to parse eigenvector file in .eig format obtained from matdyn.x using matdyn.in.

    This function reads the eigenvector file and extracts the q-points, frequencies, and
    eigenvectors for each mode.

    Parameters:
    -----------
    filename : str
        Path to the .eig file.
    nat : int
        Number of atoms.
    nqpt : int
        Number of q points.

    Returns:
    --------
    ql : numpy.ndarray
        Array of q points.
    freq : numpy.ndarray
        Array of frequencies of 3*nat modes.
    pol_q_mode : numpy.ndarray
        Array of eigenvectors of size (nqpt, 3*nat, nat, 3).
    """
    with open(filename, "r") as proj_file:
        lines = proj_file.readlines()
    qlist = []
    freq = []
    for line in lines:
        if "q =" in line:
            qpt = line.rstrip("\n").split()[2:]
            qlist.append([float(qpt[0]), float(qpt[1]), float(qpt[2])])
        if 'freq' in line:
            freq.append(float(line.split()[4]))
    qlist = np.array(qlist)
    freq = np.array(freq)
    modes = []
    for line in lines:
        if 'q' not in line and "*" not in line and 'diag' not in line and len(line) > 1:
            line_split=line.split()
            eigvec=[complex(float(line_split[1]),float(line_split[2])), complex(float(line_split[3]), float(line_split[4])), complex(float(line_split[5]),float(line_split[6]))]
            modes.append(eigvec)
    modes = np.array(modes)
    #nmodes = 3*nat
    pol_q_mode = []
    for iqpt in range(nqpt):
        modeq = modes[iqpt*3*nat*nat:(iqpt+1)*3*nat*nat]
        for imode in range(3*nat):
            pol_q_mode.append(modeq[imode*nat:(imode+1)*nat])
    pol_q_mode = np.array(pol_q_mode, dtype='object')
    pol_q_mode = pol_q_mode.reshape(nqpt,3*nat,nat,3)
    return qlist,freq,pol_q_mode

def projection(pol_q_mode,ion,indices):
    """
    Function to calculate atomic projection for mode nu and wavevector q.

    This function computes the atomic projection for a given mode nu and wavevector q.

    Parameters:
    -----------
    pol_q_mode : numpy.ndarray
        Array obtained from parse_input function containing eigenvectors.
    ion : str
        Type of atom.
    indices : list
        Indices of an ion in the scf.in input file.

    Returns:
    --------
    proj_array : numpy.ndarray
        Array of projections of size (nq, 3*nat).
    """
    nqpt = pol_q_mode.shape[0]
    nmode = pol_q_mode.shape[1]
    proj_array = np.zeros((nqpt,nmode))
    with open("phonon-{}.proj".format(ion), "w") as phon_proj:
        for iqpt in range(nqpt):
            for imode in range(nmode):
                eig=pol_q_mode[iqpt,imode]
                if eig.shape[0] > 0:
                    eigval=eig[indices]
                    proj = round(np.real(np.sum(np.conj(eigval) * eigval)),4)
                else:
                    proj = 0.
                proj_array[iqpt,imode] = proj
                phon_proj.write("{} {} {}".format(iqpt, imode, proj) + "\n")
    proj_array = proj_array.reshape(nqpt,nmode)
    return proj_array

def print_to_file():
    """
    Function to save mode resolved projection in phonon-compound.prog.gp format.

    This function reads the input file 'phonproj.in' to extract necessary information about the phonon modes and ions.
    It then calculates the mode-resolved projections for each ion and writes them to a file in the format
    'phonon-compound.proj.gp'.

    Returns:
    --------
    ions_list : list
        List of ion names for which projections are calculated and saved.
    """
    filename, nat, nqpoint, ions = reading_input("phonproj.in")
    _,_,modes = parse_input(filename,nat,nqpoint)
    with open("phonon-{}.proj.gp".format(filename.split(".")[0]), "w") as phon_proj:
        for _,ion in enumerate(ions.keys()):
            proj = projection(modes,ion,ions[ion])
            #f.write("#{}".format(x) + "\n")
            phon_proj.write(str(proj).replace(' [', '').replace('[', '').replace(']', '').replace(' ]','')+"\n")
            phon_proj.write("\n")
    return list(ions.keys())
def main():
    """
    Main function to generate phonon projections and print ion-color associations.

    This function serves as the main entry point for generating phonon projections and printing
    ion-color associations based on the projections.

    If 'phonproj.in' file is not found, it creates one.
    It then generates phonon projections for each ion
    and saves them to a file. Finally, it prints the association between ions and colors.

    """
    if not os.path.isfile("phonproj.in"):
        print("phonproj.in not found, creating one\n")
        compound = sys.argv[1]
        nkpoint = sys.argv[2]
        write_phonproj_in(compound,nkpoint,'scf.in')
    ions=print_to_file()
    color_list = ['red','blue','green','cyan','k']
    for i,_ in enumerate(ions):
        print("ION: {}, PROJECTED COLOR: {}".format(ions[i],color_list[i]))
    print("-------------------------------------------\n")
if __name__ == "__main__":
    main()
