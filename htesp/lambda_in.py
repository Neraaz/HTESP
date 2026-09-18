#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to write lambda.in files for calculating superconducting properties"""
import os
import sys
import glob

def lambda_in(compound,maxfreq,qgauss,smearing,mustar):
    """
    Create lambda.in file for lambda.x.

    Parameters:
    - compound (str): Name of the compound.
    - maxfreq (str): Maximum phonon frequency + 5 THz.
    - qgauss (str): Smearing for q-mesh integration.
    - smearing (str): Smearing type.
    - mustar (str): Coulomb potential.

    Returns:
    None
    Example:
    >>> lambda_in('SiO2', '30', '0.01', 'gauss', '0.13')
    Suppose elph.out file available. If calculation is performed in multiple
    times, elph.out (latest), elph.out1 (first), elph.out2 (second),.... files available
    """
    # Get a list of elph.out files and sort them to ensure consistent order
    elph_files = sorted(glob.glob("elph.out*"))
    if not elph_files:
        raise FileNotFoundError("no elph.out* files in {}".format(os.getcwd()))
    first_elph = elph_files.pop(0)
    elph_files.append(first_elph)
    # FIX(all): the three os.system() calls here -- a guarded "rm
    # touch_list.txt" (which never matched the file it meant to delete, the
    # scratch file is elph_list.txt), a "touch", and a "grep ... >>" per
    # elph.out file, with the file name interpolated into the shell string --
    # are replaced by an in-memory grep.  The scratch file is no longer
    # written at all; nothing else in the package reads it.
    elph_list = []
    for elph_file in elph_files:
        with open(elph_file, "r", errors="replace") as read_elph:
            for line in read_elph:
                if 'Number of q in the star' in line:
                    elph_list.append(line)
    # Read compound dynamical matrix and elph list
    with open("{}.dyn0".format(compound), "r") as read_dyn:
        dyn0 = read_dyn.readlines()
    elph_list = elph_list[::2]
    dyn0 = dyn0[2:]
    # Write lambda.in file
    with open("lambda.in", "w") as write_lambda:
        write_lambda.write(str(maxfreq) + " " + str(qgauss) + " " + smearing + "\n")
        nqsym = len(elph_list)
        write_lambda.write(str(nqsym) + "\n")
        for i in range(nqsym):
            qvec = dyn0[i].split("\n")[0]
            qweight = int(elph_list[i].split("\n")[0].split("=")[-1].split(" ")[-1])
            write_lambda.write(qvec + " " + str(qweight) + "\n")
        for i in range(nqsym):
            write_lambda.write("elph_dir/elph.inp_lambda.{}".format(i+1) + "\n")
        write_lambda.write(str(mustar)+"\n")
def main(argv=None):
    """
    Execute the main functionality.

    This function is the main entry point for executing the functionality to create
    the lambda.in file for lambda.x based on the provided command-line arguments.

    Parameters:
    None

    Returns:
    None
    """
    # FIX(21): ``main`` is the entry point htesp.workflow.run_helper() calls.
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) < 5:
        raise SystemExit("usage: lambda_in.py <compound> <maxfreq> <qgauss> "
                         "<smearing> <mustar>")
    compound, maxfreq, qgauss, smearing, mustar = argv[:5]
    lambda_in(compound,maxfreq,qgauss,smearing,mustar)
    # FIX(all): the temporary elph_list.txt is no longer written, so there is
    # nothing left to os.system("rm").
if __name__ == "__main__":
    main()
