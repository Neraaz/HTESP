#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D."""
import sys

def main():
    """
    Main function to convert AXSF file "dynmat.axsf" to cell positions.

    This function is the main entry point for the script `qe_axsf2cellpos.py`,
    which converts an AXSF file to cell positions.
    It takes command-line arguments specifying the AXSF file name, the mode
    to be considered, and a scale factor.

    Parameters:
    None (relies on command-line arguments provided through sys.argv)

    Returns:
    None
    """
    if len(sys.argv) != 4:
        print("Usage: python qe_axsf2cellpos.py <axsf> <thismode> <scale>")
        sys.exit(1)
    input_filename = sys.argv[1]
    this_mode = int(sys.argv[2]) - 1
    scale = float(sys.argv[3])
    with open(input_filename, 'r') as infile:
        lines = infile.readlines()
    #nstep = int(lines[0].split()[1])
    nion = int(lines[7].split()[0])
    print("CELL_PARAMETERS angstrom")
    for i in range(3):
        print(lines[3 + i].rstrip())
    print("\nATOMIC_POSITIONS angstrom")
    for i in range(nion):
        line = lines[8 + this_mode * (nion + 2) + i].split()
        atom_type = line[0]
        x = float(line[1]) + float(line[4]) * scale
        y = float(line[2]) + float(line[5]) * scale
        z = float(line[3]) + float(line[6]) * scale
        print(f"{atom_type} {x:.6f} {y:.6f} {z:.6f}")
if __name__ == "__main__":
    main()
