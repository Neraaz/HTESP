#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to create matdyn.x input files"""
import sys
import numpy as np

#: default name of the atomic-mass scratch file written by the workflow layer
DEFAULT_MASS_FILE = "mass.dat"
def matdyn_dos(mpid=None, compound=None, prefix2=None, mass_file=None, argv=None):
    """
    Prepares input file for phonon DOS and superconductivity-related properties.

    Parameters:
    - mpid (str): Materials Project ID.
    - compound (str): Chemical compound.
    - prefix2 (str): Prefix for file names.

    Returns:
    Creates input file named matdyn-mpid-compound-dos.in.
    """
    # FIX(21): ``matdyn_dos`` is the entry point htesp.workflow.run_helper()
    # calls.  FIX(22): an optional fourth argument names the mass file, which
    # used to be the hard-coded 'mass.dat' in the current directory.
    argv = list(sys.argv[1:] if argv is None else argv)
    if mpid is None or compound is None or prefix2 is None:
        if len(argv) < 3:
            raise SystemExit("usage: matdyn_dos.py <mpid> <compound> <prefix> "
                             "[mass file]")
        mpid, compound, prefix2 = argv[0], argv[1], argv[2]
    if mass_file is None:
        mass_file = argv[3] if len(argv) > 3 else DEFAULT_MASS_FILE
    # File names
    freq = prefix2.replace("'", "") + "-dos.freq"
    frc = prefix2.replace("'", "") + ".fc"
    # Determine the number of atoms
    mass=np.loadtxt(mass_file)
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

main = matdyn_dos
if __name__ == "__main__":
    matdyn_dos()
