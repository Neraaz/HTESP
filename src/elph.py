#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to prepare QE el-ph calculations"""
import sys
import numpy as np

def elph_in():
    """

    Default filename: 'elph-mpid-compound.in'
    The file is saved inside the 'elph_dir' folder.
    These scripts are executed within the 'create-inputs' bash script.

    The function then generates the input file for electron-phonon coupling calculation.

    Parameters:
    -----------------
    None

    Returns:
    -----------------
    None
    """
    mass=np.loadtxt('mass.dat')
    qpt=np.loadtxt('qpoint.dat')
    if mass.ndim > 0:
        nion = mass.shape[0]
    else:
        nion = 1
        mass = [mass]
    mpid = sys.argv[1]
    comp = sys.argv[2]
    prefix2=sys.argv[3]
    dynmat = prefix2.replace("'", "") + ".dyn"
    print("q-mesh: {}".format(qpt))
    with open("elph-{}-{}.in".format(mpid,comp), 'w') as elph_write:
        elph_write.write("electron phonon coupling \n")
        elph_write.write("&inputph" + "\n")
        #elph_write.write("alpha_mix=0.3, nmix_ph=8," + "\n")
        elph_write.write("tr2_ph=1.0d-16," + "\n")
        elph_write.write("prefix={},".format(prefix2) + "\n")
        elph_write.write("fildvscf='aldv'," + "\n")
        for i in range(1,nion+1):
            elph_write.write("amass({})={},".format(i,mass[i-1]) + "\n")
        elph_write.write("outdir='./'," + "\n")
        elph_write.write("fildyn='{}',".format(dynmat) + "\n")
        elph_write.write("electron_phonon='interpolated'," + "\n")
        elph_write.write("el_ph_sigma=0.005," + "\n")
        elph_write.write("el_ph_nsigma=10," + "\n")
        elph_write.write("trans=.true.," + "\n")
        elph_write.write("ldisp=.true." + "\n")
        elph_write.write("nq1={},nq2={},nq3={}".format(int(qpt[0]),int(qpt[1]),int(qpt[2])) + "\n")
        elph_write.write("/" + "\n")

if __name__ == "__main__":
    elph_in()
