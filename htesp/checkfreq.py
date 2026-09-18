#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to check if frequency is less than cutoff value"""
import sys
import numpy as np
#: default name of the marker file written when a soft mode is found
DEFAULT_FLAG_FILE = "freq.dat"
def check_freq(filename, flag_file=DEFAULT_FLAG_FILE):
    """
    If frequencies less than -1 THz are found, a file 'freq.dat' is created.
    It runs within the 'checkfreq-scan' bash script.

    Parameters:
    - filename (str): Phonon frequency filename in 'name.freq.gp' format created in QE calculations.
    - flag_file (str): FIX(22) -- name of the marker file to write.  It used to
      be the hard-coded 'freq.dat' in the current directory, so two materials
      checked from the same directory overwrote each other's marker.  The
      default is unchanged.

    Returns:
    None
    """
    # Load the frequency data from the provided filename
    data = np.loadtxt(filename)
    data = data[:, 1:]
    # Check if any frequency is less than -1 THz
    if np.any(data < -33.356): # -33.356 corresponds to -1 THz
        print("{}:".format(filename) + " Negative frequency smaller than -1.0 THz present\n")
        with open(flag_file, 'w') as write_freq:
            write_freq.write("Negative frequencies\n")
def main(file_name=None, flag_file=None, argv=None):
    """
    main function.

    FIX(21)/FIX(22): ``main`` is the entry point htesp.workflow.run_helper()
    calls; an optional second argument names the marker file.
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    if file_name is None:
        if not argv:
            raise SystemExit("usage: checkfreq.py <name.freq.gp> [flag file]")
        file_name = argv[0]
    if flag_file is None:
        flag_file = argv[1] if len(argv) > 1 else DEFAULT_FLAG_FILE
    check_freq(file_name, flag_file)
if __name__ == "__main__":
    main()
