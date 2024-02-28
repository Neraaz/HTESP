#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D."""
import warnings
from ase.io import cif,espresso
warnings.filterwarnings('ignore')
def scf_tocif(filename):
    """
    Function to convert QE scf.in file to .cif
    prameters
    ---------
    filename : QE input file
    """
    filename = espresso.read_espresso_in(filename)
    cif.write_cif('relax.cif',filename)
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
def main():
    """
    main function
    """
    filename=scf_tocif('scf.in')
    cellpar(filename)
if __name__ == "__main__":
    main()
