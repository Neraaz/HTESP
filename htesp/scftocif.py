#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to extract structures"""
import sys
import os
import warnings
from ase.io import cif,espresso,vasp
from htesp.cif_to_gsinput import pymatgen_cif
warnings.filterwarnings('ignore')
def scf_tocif(filename):
    """
    Function to convert QE scf.in file to .cif
    prameters
    ---------
    filename : QE input file
    """
    # The name, not just the basename, used to be compared, so a path such as
    # 'relax/scf.in' fell through with ``filename`` still a string.
    basename = os.path.basename(str(filename))
    if basename == "scf.in" or basename.endswith(".in"):
        filename = espresso.read_espresso_in(filename)
    elif basename in ("POSCAR", "CONTCAR"):
        filename = vasp.read_vasp(filename)
    else:
        raise ValueError(
            "Either scf.in or POSCAR expected, got {!r}".format(filename))
    cif.write_cif('relax.cif',filename)
    pymatgen_cif('relax.cif')
    return filename
def cellpar(filename):
    """
    Function to compute lattice parameters
    parameters
    ------------
    filename : QE scf input file
    """
    alat,blat,clat,alp,bet,gam = filename.get_cell_lengths_and_angles()
    with open("cellpar.in", "w") as cell_par:
        cell_par.write("{} {} {} ".format(round(alat,5),round(blat,5),round(clat,5)))
        cell_par.write("{} {} {}".format(round(alp,5),round(bet,5),round(gam,5)))
    return alat,blat,clat,alp,bet,gam
def main(filename=None, argv=None):
    """
    main function.

    FIX(21): ``main`` is the entry point htesp.workflow.run_helper() calls; the
    structure file may also be passed directly.
    """
    if filename is None:
        argv = list(sys.argv[1:] if argv is None else argv)
        filename = argv[0] if argv else "scf.in"
    filename=scf_tocif(filename)
    cellpar(filename)
if __name__ == "__main__":
    main()
