#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D."""
import sys
import numpy as np
def matdyn_in():
    """
    Prepares input file for phonon band structure and eigenvectors.

    Parameters:
    - mpid (str): Materials ID.
    - compound (str): Chemical compound.
    - prefix2 (str): Prefix for file names.

    Returns:
    Creates input file named matdyn-mpid-compound.in.
    """
    mpid = sys.argv[1]
    compound = sys.argv[2]
    prefix2=sys.argv[3]
    freq = prefix2.replace("'", "") + ".freq"
    frc = prefix2.replace("'", "") + ".fc"
    eig = prefix2.replace("'", "") + ".eig"
    mass=np.loadtxt('mass.dat')
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

if __name__ == "__main__":
    matdyn_in()
