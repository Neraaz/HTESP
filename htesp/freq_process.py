#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to process phonon frequency"""
import sys
import numpy as np

def freq_process():
    """
    Function to process phonon frequencies from compound.freq.gp (gnuplot)
    and convert them into freq.plot (Xmgrace) format.

    This function reads phonon frequencies from a file named compound.freq.gp,
    which is in gnuplot format, and writes the frequencies to a file named freq.plot
    in Xmgrace format.

    Parameters:
    -----------------
    compound : str
        The name of the compound. It is provided as a command-line argument.

    Returns:
    -----------------
    None
    """
    compound = sys.argv[1]
    filename = compound+".freq.gp"
    # FIX(24): the output file was hard-coded; an optional second argument now
    # redirects it, so a caller can write per-material output.
    outfile = sys.argv[2] if len(sys.argv) > 2 else "freq.plot"
    data = np.loadtxt(filename)
    row = data[:,0]
    col = data[:,1:]
    with open(outfile, "w") as file1:
        nrow = col.shape[0]
        ncol = col.shape[1]
        for i in range(ncol):
            for j in range(nrow):
                # FIX(24): was `row[i]` - the *band* index was used to index the
                # q-axis, so every band got the same wrong abscissa and the
                # loop raised IndexError as soon as nband > nqpoint (which is
                # the normal case for a dense q-mesh).
                file1.write(str(row[j]) + " " + str(col[j,i]) + "\n")
            file1.write("\n")
if __name__ == "__main__":
    freq_process()
