#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Converts .axsf file to cell positions"""
import sys

def main(argv=None):
    """
    Main function to convert AXSF file "dynmat.axsf" to cell positions.

    FIX(21): ``main`` is the entry point htesp.workflow.run_helper() calls
    (with ``capture=True``); it accepts an explicit argument list.

    This function is the main entry point for the script `qe_axsf2cellpos.py`,
    which converts an AXSF file to cell positions.
    It takes command-line arguments specifying the AXSF file name, the mode
    to be considered, and a scale factor.

    Parameters:
    None (relies on command-line arguments provided through sys.argv)

    Returns:
    None
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 3:
        raise SystemExit("Usage: qe_axsf2cellpos.py <axsf> <thismode> <scale>")
    input_filename = argv[0]
    this_mode = int(argv[1]) - 1
    scale = float(argv[2])
    with open(input_filename, 'r') as infile:
        lines = infile.readlines()
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
