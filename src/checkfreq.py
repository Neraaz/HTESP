#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to check if frequency is less than cutoff value"""
import sys
import numpy as np
def check_freq(filename):
    """
    If frequencies less than -1 THz are found, a file 'freq.dat' is created.
    It runs within the 'checkfreq-scan' bash script.

    Parameters:
    - filename (str): Phonon frequency filename in 'name.freq.gp' format created in QE calculations.

    Returns:
    None
    """
    # Load the frequency data from the provided filename
    data = np.loadtxt(filename)
    data = data[:, 1:]
    # Check if any frequency is less than -1 THz
    if np.any(data < -33.356): # -33.356 corresponds to -1 THz
        print("{}:".format(filename) + " Negative frequency smaller than -1.0 THz present\n")
        with open('freq.dat', 'w') as write_freq:
            write_freq.write("Negative frequencies\n")
def main():
    """
    main function
    """
    file_name = sys.argv[1]
    check_freq(file_name)
if __name__ == "__main__":
    main()
