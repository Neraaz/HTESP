#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to handle kpath for high-symmetry lines"""
import os
import sys
from ase.io import espresso,vasp
import pylab
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.core import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from check_json import config

def kpath(filename,npoint,kcutoff):
    """
    Function to write k-point mesh along the high-symmetry path of the Brillouin zone (BZ).

    Parameters:
    - filename (str): Input file to read, which contains the structure or VASP 'POSCAR' file.
    - npoint (int): Size of k-point mesh.
    - kcutoff (int): Cutoff to use for k-point path in ASE.
      '0' means full Brillouin zone,
      '1' means to remove one symmetry point from the BZ path.

    Returns:
    - kpoints (numpy.ndarray): K-mesh of size (npoint, 3).
    - sympoint (list): K-point in linear axis ready for plotting after processing.
    - symname (list): Naming for sympoint.
    - kpt (list): K-point in linear axis without processing.
    - spt (list): K-point in linear axis without processing at high-symmetry points.
    - sym (list): Naming for spt.
    Example:
    >>> kpoints, sympoint, symname, kpt, sym, spt = kpath("POSCAR", 100, 0)
    """
    input_data = config()
    try:
        # Attempt to read the input file as an espresso file.
        file_name = espresso.read_espresso_in(filename)
    except:
        # If reading as an espresso file fails, try VASP file.
        file_name = vasp.read_vasp(filename)
    inp_dict = input_data['download']['inp']
    if "kpath_pbc" in inp_dict.keys():
        pbc = inp_dict['kpath_pbc']
    else:
        pbc = [1, 1, 1]
    if pbc is None:
        pbc = [1, 1, 1]
    # Get the band path from the cell.
    bandpath = file_name.cell.bandpath(pbc=pbc)
    # Determine the path based on the provided cutoff.
    if kcutoff > 0:
        path=bandpath.path[:kcutoff]
    else:
        path=bandpath.path[:None]
    # Generate the band path with the specified number of points.
    bandpath = file_name.cell.bandpath(path,npoints=npoint,pbc=pbc)
    file_name.cell.bandpath().plot()
    pylab.savefig("BZ.pdf", format='pdf', bbox_inches='tight')
    # Retrieve k-points, linear axis, and symmetry names.
    kpoints = bandpath.kpts
    sympoint = bandpath.get_linear_kpoint_axis()[1]
    symname = bandpath.get_linear_kpoint_axis()[2]
    # Process the combined labels.
    sympoint2 = []
    idx = []
    for i in range(sympoint.shape[0]):
        # Check if the current point has the same position as the previous one.
        if sympoint[i-1] == sympoint[i]:
            # Append the previous label to the current point.
            symadd = symname[i-1] + "|" + symname[i]
            symname[i] = symadd
            # Store the index of the previous point to be removed.
            idx.append(i-1)
        else:
            sympoint2.append(round(sympoint[i],8))
    # Remove the redundant labels from the list.
    rmv = len(idx)
    while rmv > 0:
        symname.pop(idx[rmv-1])
        rmv = rmv - 1
    sympoint = sympoint2
    # Retrieve additional data.
    spt = bandpath.get_linear_kpoint_axis()[1]
    sym = bandpath.get_linear_kpoint_axis()[2]
    kpt = bandpath.get_linear_kpoint_axis()[0]
    return kpoints,sympoint,symname,kpt,sym,spt

def printk():
    """
    Print k-point mesh within high-symmetry points and between them.

    This function generates and prints k-points for quantum mechanics
    calculations, focusing on the Brillouin zone (BZ) and high-symmetry points.

    Parameters:
    None

    Returns:
    None

    Usage:
    The function expects command-line arguments in the following order:
    - sys.argv[2]: Filename containing the structure or VASP 'POSCAR' file.
    - sys.argv[3]: Number of k-points.
    - sys.argv[4]: Cutoff for the k-point path in ASE.
     '0' for the full Brillouin zone, '1' to remove one symmetry point from the BZ path.
    - sys.argv[5]: Weight of the k-point (optional).

    Output:
    The function generates two files in the 'scf_dir' directory:
    - 'kpathlines.dat': Contains the k-point mesh within high-symmetry points and between them.
    - 'kspecial-points.dat': Lists the high-symmetry points.

    """
    filename = sys.argv[2]
    nkpoint = int(sys.argv[3])
    kcutoff = int(sys.argv[4])
    weight = int(sys.argv[5])
    kpts,sympoint,symname,_,_,_ = kpath(filename,nkpoint,kcutoff)
    nkpt = kpts.shape[0]
    if not os.path.isdir("scf_dir"):
        os.system("mkdir scf_dir")
    with open('scf_dir/kpathlines.dat', 'w') as kpathlines:
        kpathlines.write('K_POINTS crystal\n')
        kpathlines.write(str(nkpt) + '\n')
        for i in range(nkpt):
            kpathlines.write(str(round(kpts[i][0],8)) + " " + str(round(kpts[i][1],8)))
            kpathlines.write(" " + str(round(kpts[i][2],8)) + " " + str(weight) + "\n")
    with open('scf_dir/kspecial-points.dat', 'w') as special_points:
        special_points.write(str(sympoint) + "\n")
        special_points.write(str(symname))

def make_line_kpt(filename="KPOINTS"):
    """
    Function to create high symmetry points for line mode using pymatgen.

    This function generates high symmetry points for line
    mode using the provided structure file ("POSCAR").
    The number of k-points for the line mode is specified
    by the second command-line argument.

    Parameters:
    None

    Returns:
    None

    Usage:
    The function expects a command-line argument specifying the number of k-points for line mode.

    Output:
    The function creates a KPOINTS file containing high symmetry points for line mode.

    """
    input_data = config()
    inp_dict = input_data['download']['inp']
    if "kpath_pbc" in inp_dict.keys():
        pbc = inp_dict['kpath_pbc']
    else:
        pbc = [1, 1, 1]
    if pbc is None:
        pbc = [1, 1, 1]
    nkpt=int(sys.argv[2])
    if pbc == [1, 1, 1]:
        struct = Structure.from_file("POSCAR")
        k_path = HighSymmKpath(struct)
        kpts = Kpoints.automatic_linemode(divisions=nkpt,ibz=k_path)
        kpts.write_file(filename)
    else:
        data = vasp.read_vasp("POSCAR")
        bandpath = data.cell.bandpath(pbc=pbc)
        path = bandpath.path
        special_points = bandpath.special_points
        with open(filename, "w") as write_kpt:
            write_kpt.write("Line_mode KPOINTS file\n")
            write_kpt.write(str(nkpt) + "\n")
            write_kpt.write("Line_mode\n")
            write_kpt.write("Reciprocal\n")
            for i in range(len(path) - 1):
                xpath = path[i]
                ypath = path[i+1]
                array1 = special_points[xpath]
                array2 = special_points[ypath]
                write_kpt.write(str(array1[0]) + " " + str(array1[1]) + " " + str(array1[2]) + " ! " + xpath + "\n")
                write_kpt.write(str(array2[0]) + " " + str(array2[1]) + " " + str(array2[2]) + " ! " + ypath +  "\n")
                write_kpt.write("\n")
        os.system(f"""sed -i '$d' {filename}""")
def main():
    """
    Main function for executing different modes of k-point generation.

    Parameters:
    None

    Returns:
    None
    """
    input_data = config()
    mode = sys.argv[1]
    if mode == "line":
        if os.path.isfile("KPT_OPT"):
            make_line_kpt("KPOINTS_OPT")
        else:
            make_line_kpt()
    elif mode == "point":
        printk()
    else:
        print("Either line or point mode available\n")
if __name__ == "__main__":
    input_data = config()
    if "kpt_opt" in input_data.keys():
        kpt_opt = input_data["kpt_opt"]
    else:
        kpt_opt = False
    if kpt_opt:
        os.system("touch KPT_OPT")
    main()
