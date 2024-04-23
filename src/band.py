#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to write input files for QE band calculations"""
import sys

def write_band_input(material_id, compound_name, prefix, dynamic_matrix):
    """
    Function to write band.in file to process a bandstructure
    .freq.gp file in Quantum ESPRESSO calculations.

    Parameters:
    - material_id (str): Materials ID of the compound..
    - compound_name (str): Compound name.
    - prefix (str): Prefix used in the QE calculations.
    - dynamic_matrix (str): Bandstructure file to process.

    Returns:
    None
    """
    # Write input parameters to band.in file
    with open("scf_dir/band-{}-{}.in".format(material_id,compound_name), 'w') as band_in:
        band_in.write("&BANDS\n")
        band_in.write(f"prefix={prefix},\n")
        band_in.write("outdir='./',\n")
        band_in.write(f"filband='{dynamic_matrix}',\n")
        band_in.write("lsym=.true.\n")
        band_in.write("/\n")
    with open("scf_dir/bandproj-{}-{}.in".format(material_id,compound_name), 'w') as file1:
        file1.write("&projwfc" + "\n")
        file1.write(f"prefix={prefix}," + "\n")
        file1.write("outdir='./'," + "\n")
        file1.write("filproj='proj.out'," + "\n")
        file1.write("/" + "\n")
def main():
    """
    main function
    """
    material_id = sys.argv[1]
    compound_name = sys.argv[2]
    prefix = sys.argv[3]
    dynamic_matrix = prefix.replace("'", "") + ".dat"
    write_band_input(material_id, compound_name, prefix, dynamic_matrix)
if __name__ == "__main__":
    main()
