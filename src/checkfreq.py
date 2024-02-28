#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D."""
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
    data = np.loadtxt(filename)
    data = data[:, 1:]
    if np.any(data < -33.356):
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
