#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to create input file to compute phonon dispersion at high-symmetry BZ path"""
import sys
import numpy as np

#: default name of the atomic-mass scratch file written by the workflow layer
DEFAULT_MASS_FILE = "mass.dat"
def matdyn_in(mpid=None, compound=None, prefix2=None, mass_file=None, argv=None):
    """
    Prepares input file for phonon band structure and eigenvectors.

    Parameters:
    - mpid (str): Materials ID.
    - compound (str): Chemical compound.
    - prefix2 (str): Prefix for file names.

    Returns:
    Creates input file named matdyn-mpid-compound.in.
    """
    # Look for matdyn_dos.py
    # FIX(21): ``matdyn_in`` is the entry point htesp.workflow.run_helper()
    # calls.  FIX(22): an optional fourth argument names the mass file, which
    # used to be the hard-coded 'mass.dat' in the current directory.
    argv = list(sys.argv[1:] if argv is None else argv)
    if mpid is None or compound is None or prefix2 is None:
        if len(argv) < 3:
            raise SystemExit("usage: matdyn.py <mpid> <compound> <prefix> "
                             "[mass file]")
        mpid, compound, prefix2 = argv[0], argv[1], argv[2]
    if mass_file is None:
        mass_file = argv[3] if len(argv) > 3 else DEFAULT_MASS_FILE
    freq = prefix2.replace("'", "") + ".freq"
    frc = prefix2.replace("'", "") + ".fc"
    eig = prefix2.replace("'", "") + ".eig"
    mass=np.loadtxt(mass_file)
    if mass.ndim > 0:
        nat = mass.shape[0]
    else:
        nat = 1
        mass = [mass]
    with open("matdyn-{}-{}.in".format(mpid,compound), 'w') as matdyn:
        matdyn.write("&input" + "\n")
        matdyn.write("asr='simple'," + "\n")
        for i in range(1,nat+1):
            matdyn.write("amass({})={},".format(i,mass[i-1]) + "\n")
        matdyn.write("flfrc='{}',".format(frc) + "\n")
        matdyn.write("flfrq='{}',".format(freq) + "\n")
        matdyn.write("fleig='{}',".format(eig) + "\n")
        matdyn.write("la2F=.true.," + "\n")
        matdyn.write("dos=.false." + "\n")
        matdyn.write("q_in_cryst_coord=.true." + "\n")
        matdyn.write("/" + "\n")
        with open("scf_dir/kpathlines.dat", "r") as kpathline:
            lines=kpathline.readlines()
        nkpt=int(lines[1])
        matdyn.write(str(nkpt) + "\n")
        for i in range(nkpt):
            matdyn.write(lines[2+i].split()[0] + " " + lines[2+i].split()[1])
            matdyn.write(" " + lines[2+i].split()[2] + " " + str(0) + "\n")

main = matdyn_in
if __name__ == "__main__":
    matdyn_in()
