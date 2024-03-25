#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to create matdyn.x input files"""
import sys
import numpy as np
def matdyn_dos():
    """
    Prepares input file for phonon DOS and superconductivity-related properties.

    Parameters:
    - mpid (str): Materials Project ID.
    - compound (str): Chemical compound.
    - prefix2 (str): Prefix for file names.

    Returns:
    Creates input file named matdyn-mpid-compound-dos.in.
    """
    # Get command-line arguments
    mpid = sys.argv[1]
    compound = sys.argv[2]
    prefix2=sys.argv[3]
    # File names
    freq = prefix2.replace("'", "") + "-dos.freq"
    frc = prefix2.replace("'", "") + ".fc"
    # Determine the number of atoms
    mass=np.loadtxt('mass.dat')
    if mass.ndim > 0:
        nat = mass.shape[0]
    else:
        nat = 1
        mass = [mass]
    # Write input file
    with open("matdyn-{}-{}-dos.in".format(mpid,compound), 'w') as matdos:
        matdos.write("&input" + "\n")
        matdos.write("asr='simple'," + "\n")
        for i in range(1,nat+1):
            matdos.write("amass({})={},".format(i,mass[i-1]) + "\n")
        matdos.write("flfrc='{}',".format(frc) + "\n")
        matdos.write("flfrq='{}',".format(freq) + "\n")
        matdos.write("dos=.true.," + "\n")
        matdos.write("fldos='phonon.dos'," + "\n")
        matdos.write("nk1=10,nk2=10,nk3=10,ndos=2000," + "\n")
        matdos.write("/" + "\n")

if __name__ == "__main__":
    matdyn_dos()
