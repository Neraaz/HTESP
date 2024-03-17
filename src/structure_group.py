#!/usr/bin/env python
#"""Written by Niraj K. Nepal, Ph.D."""
"""Module to group structures from different database"""
# coding: utf-8
import os
from pymatgen.core import Structure
from pymatgen.io.pwscf import PWInput
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher

def read_structure_file(filepath):
    """
    Read structure from a file.

    Args:
        filepath (str): Path to the structure file.

    Returns:
        Structure: Pymatgen Structure object.
    """
    try:
        return Structure.from_file(filepath)
    except:
        return PWInput(filepath).structure

def create_structure_list(input_file):
    """
    Create a list of structures from the input file.

    Args:
        input_file (str): Path to the input file.

    Returns:
        list: List of Pymatgen Structure objects.
    """
    struc_list = []
    with open(input_file, "r") as read_input:
        lines = read_input.readlines()
    for line in lines:
        line = line.split(" ")
        mpid = line[1]
        compound = line[2].strip()
        try:
            struc_path = f"R{mpid}-{compound}/relax/POSCAR"
            struc = read_structure_file(struc_path)
        except:
            struc_path = f"scf_dir/scf-{mpid}.in"
            struc = read_structure_file(struc_path)
        struc = SpacegroupAnalyzer(structure=struc, symprec=0.1).get_primitive_standard_structure()
        site_prop = [mpid] * len(struc.sites)
        struc.add_site_property("id", site_prop)
        struc_list.append(struc)
    return struc_list

def write_unique_structures(unique_structures, output_file):
    """
    Write unique structures to an output file and copy files to a directory.

    Args:
        unique_structures (list): List of unique Pymatgen Structure objects.
        output_file (str): Path to the output file.
    """
    if not os.path.isdir("filtered_inputs"):
        os.mkdir("filtered_inputs")
    with open(output_file, "w") as write_mpid:
        for i, group in enumerate(unique_structures):
            extract_struc = group[0]
            mpid = extract_struc.site_properties['id'][0]
            comp = str(extract_struc.composition).replace(" ", "")
            write_mpid.write(f"v{i+1} {mpid} {comp}\n")
            try:
                os.system(f"cp -r R{mpid}* filtered_inputs/")
            except FileNotFoundError:
                os.system(f"cp scf_dir/scf-{mpid}.in filtered_inputs/")

def main():
    """
    Main function to process structures from an input file, group them into unique structures,
    and write the results to an output file.

    Reads structure information from the 'mpid.in' file, creates a list of structures,
    groups them into unique structures using StructureMatcher, and writes the unique structures
    along with associated metadata to the 'mpid-new.in' file.

    Returns:
        None
    """
    input_file = "mpid.in"
    struc_list = create_structure_list(input_file)
    unique_structures = StructureMatcher().group_structures(struc_list)
    output_file = "mpid-new.in"
    write_unique_structures(unique_structures, output_file)
    print("Check filtered_inputs folder and mpid-new.in file\n")

if __name__ == "__main__":
    main()
