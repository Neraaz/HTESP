#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to handle kpath for high-symmetry lines"""
import os
import re
import sys
from ase.io import espresso,vasp
# FIX(12): the BZ plot used to be drawn on *every* kpath() call with whatever
# matplotlib backend happened to be active (only plot.py selected 'Agg'), so a
# head-less batch run could block or fail in a GUI backend.  Select a
# non-interactive backend before pylab/pyplot is imported here.
import matplotlib
matplotlib.use("Agg")
import pylab
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.core import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from htesp.check_json import config

# FIX(11): a high-symmetry label can be longer than one character ("G1", "K1")
# and a break in the path arrives as a separator ("K|U").  Slicing the ASE path
# *string* by characters cut multi-character labels in half.
_LABEL = re.compile(r"[A-Za-z][0-9]*|[|,]")


def path_tokens(path):
    """
    Split an ASE band-path string into whole high-symmetry labels.

    ``'GXWK|UX'`` -> ``['G', 'X', 'W', 'K', '|', 'U', 'X']``.  Separators are
    kept so that ``''.join(tokens)`` reproduces the original path string.

    Parameters
    ----------
    path : str
        Band path as produced by ``cell.bandpath().path``.

    Returns
    -------
    list of str
        Labels and separators, in order.
    """
    return _LABEL.findall(str(path))


def cut_path(path, kcutoff):
    """
    Keep the first ``kcutoff`` high-symmetry *labels* of ``path``.

    ``kcutoff <= 0`` returns the full path.  Separators never count as labels
    and a trailing separator is dropped.

    Parameters
    ----------
    path : str
        Band path as produced by ``cell.bandpath().path``.
    kcutoff : int
        Number of labels to keep.

    Returns
    -------
    str
        The truncated path string.

    Example
    -------
    >>> cut_path("G1XWK|UX", 3)
    'G1XW'
    """
    if kcutoff is None or int(kcutoff) <= 0:
        return path
    kept = []
    nlabel = 0
    for token in path_tokens(path):
        if token in ("|", ","):
            kept.append(token)
            continue
        if nlabel >= int(kcutoff):
            break
        kept.append(token)
        nlabel += 1
    while kept and kept[-1] in ("|", ","):
        kept.pop()
    return "".join(kept)


def kpath(filename,npoint,kcutoff,plot_bz=False,bz_file="BZ.pdf"):
    """
    Function to write k-point mesh along the high-symmetry path of the Brillouin zone (BZ).

    Parameters:
    - filename (str): Input file to read, which contains the structure or VASP 'POSCAR' file.
    - npoint (int): Size of k-point mesh.
    - kcutoff (int): Cutoff to use for k-point path in ASE.
      '0' means full Brillouin zone,
      'n' > 0 keeps only the first n high-symmetry labels of the BZ path.
    - plot_bz (bool): Write the Brillouin-zone sketch to `bz_file`.
      Default: False (see FIX(12); only `printk` turns it on).
    - bz_file (str): Output file for the BZ sketch. Default: 'BZ.pdf'.

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
    # FIX: a bare ``except`` swallowed everything, KeyboardInterrupt included.
    except (OSError, ValueError, IndexError, KeyError):
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
    # FIX(11): tokenise the path into whole labels before truncating it, so a
    # multi-character label such as "G1" or "K1" is never cut in half.
    path = cut_path(bandpath.path, kcutoff)
    # Generate the band path with the specified number of points.
    bandpath = file_name.cell.bandpath(path,npoints=npoint,pbc=pbc)
    # FIX(12): the BZ sketch is now opt-in; it used to be redrawn and saved on
    # every single call (kpoint_path, create_wt_inputs, plot, htepc, ...).
    if plot_bz:
        file_name.cell.bandpath().plot()
        pylab.savefig(bz_file, format='pdf', bbox_inches='tight')
        pylab.close("all")
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

def printk(out_dir="scf_dir", plot_bz=True):
    """
    Print k-point mesh within high-symmetry points and between them.

    This function generates and prints k-points for quantum mechanics
    calculations, focusing on the Brillouin zone (BZ) and high-symmetry points.

    Parameters:
    - out_dir (str): Directory the two .dat files are written to.
      Default: 'scf_dir'. FIX(13): this used to be hard-coded, so a caller
      could not direct the output at a per-material directory.
    - plot_bz (bool): Write 'BZ.pdf' next to the output. Default: True,
      because the workflow harvests BZ.pdf from this call (see FIX(12)).

    Returns:
    None

    Usage:
    The function expects command-line arguments in the following order:
    - sys.argv[2]: Filename containing the structure or VASP 'POSCAR' file.
    - sys.argv[3]: Number of k-points.
    - sys.argv[4]: Cutoff for the k-point path in ASE.
     '0' for the full Brillouin zone, 'n' to keep only the first n labels.
    - sys.argv[5]: Weight of the k-point.
    - sys.argv[6]: Output directory (optional, defaults to 'scf_dir').

    Output:
    The function generates two files in `out_dir`:
    - 'kpathlines.dat': Contains the k-point mesh within high-symmetry points and between them.
    - 'kspecial-points.dat': Lists the high-symmetry points.

    """
    filename = sys.argv[2]
    nkpoint = int(sys.argv[3])
    kcutoff = int(sys.argv[4])
    weight = int(sys.argv[5])
    # FIX(13): keep the positional CLI contract - the optional 6th argument
    # overrides the default output directory.
    if len(sys.argv) > 6 and sys.argv[6]:
        out_dir = sys.argv[6]
    kpts,sympoint,symname,_,_,_ = kpath(filename,nkpoint,kcutoff,plot_bz=plot_bz)
    nkpt = kpts.shape[0]
    # FIX: isdir-then-"mkdir" through the shell replaced by makedirs(exist_ok).
    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, 'kpathlines.dat'), 'w') as kpathlines:
        kpathlines.write('K_POINTS crystal\n')
        kpathlines.write(str(nkpt) + '\n')
        for i in range(nkpt):
            kpathlines.write(str(round(kpts[i][0],8)) + " " + str(round(kpts[i][1],8)))
            kpathlines.write(" " + str(round(kpts[i][2],8)) + " " + str(weight) + "\n")
    with open(os.path.join(out_dir, 'kspecial-points.dat'), 'w') as special_points:
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
                if i < len(path) - 2:
                    write_kpt.write("\n")
        # FIX: the trailing blank line used to be removed by shelling out to
        # `sed -i '$d' <filename>`, interpolating a path into a shell string.
        # It is simply not written any more.
def main():
    """
    Main function for executing different modes of k-point generation.

    Parameters:
    None

    Returns:
    None
    """
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
    INPUT_DATA = config()
    KPT_OPT = INPUT_DATA.get("kpt_opt", False)
    if KPT_OPT:
        # FIX: `os.system("touch KPT_OPT")` replaced by a pure-Python touch.
        with open("KPT_OPT", "a", encoding="utf-8"):
            pass
    main()
