#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D."""
import sys

def dos_in():
    """
    Function to write processing file for density of states (DOS).
    Default: dos-mpid-compound.in
    Also writes input file for partial DOS (PDOS) calculation.
    Default: pdos-mpid-compound.in
    Can find these files inside scf_dir folder.
    This script is run within 'create-inputs' bash script.

    Usage:
    ------
    python script.py mpid comp prefix2

    mpid : str
        Material ID.

    comp : str
        Compound name.

    prefix2 : str
        Prefix for DOS and PDOS files.

    Returns:
    --------
    None
    """
    mpid = sys.argv[1]
    comp = sys.argv[2]
    prefix2=sys.argv[3]
    dynmat = prefix2.replace("'", "") + ".dos"
    dynmat1 = prefix2.replace("'", "") + ".pdos"
    with open("scf_dir/dos-{}-{}.in".format(mpid,comp), 'w') as file1:
        file1.write("&dos" + "\n")
        file1.write("prefix={},".format(prefix2) + "\n")
        file1.write("outdir='./'," + "\n")
        file1.write("fildos='{}',".format(dynmat) + "\n")
        file1.write("DeltaE=0.01" + "\n")
        file1.write("/" + "\n")
    with open("scf_dir/pdos-{}-{}.in".format(mpid,comp), 'w') as file1:
        file1.write("&projwfc" + "\n")
        file1.write("prefix={},".format(prefix2) + "\n")
        file1.write("outdir='./'," + "\n")
        file1.write("filpdos='{}',".format(dynmat1) + "\n")
        file1.write("DeltaE=0.01" + "\n")
        file1.write("/" + "\n")
if __name__ == "__main__":
    dos_in()
