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
    first_elph = elph_files.pop(0)
    elph_files.append(first_elph)
    # Create elph_list.txt containing the number of q points for each elph.out file
    if os.path.isfile("touch_list.txt"):
        os.system("rm touch_list.txt")
    os.system("touch elph_list.txt")
    for i,_ in enumerate(elph_files):
        os.system("grep 'Number of q in the star' {}".format(elph_files[i]) + ">> elph_list.txt")
    # Read compound dynamical matrix and elph list
    with open("{}.dyn0".format(compound), "r") as read_dyn:
        dyn0 = read_dyn.readlines()
    with open("elph_list.txt", "r") as read_elph_list:
        elph_list = read_elph_list.readlines()
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
def main():
    """
    Execute the main functionality.

    This function is the main entry point for executing the functionality to create
    the lambda.in file for lambda.x based on the provided command-line arguments.

    Parameters:
    None

    Returns:
    None
    """
    compound = sys.argv[1]
    maxfreq = sys.argv[2]
    qgauss = sys.argv[3]
    smearing = sys.argv[4]
    mustar = sys.argv[5]
    lambda_in(compound,maxfreq,qgauss,smearing,mustar)
    os.system("rm elph_list.txt")
if __name__ == "__main__":
    main()
