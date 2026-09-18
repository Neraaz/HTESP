#!/usr/bin/env python
#"""Written by Niraj K. Nepal, Ph.D."""
"""Module to group structures from different database"""
# coding: utf-8
import glob
import os
import shutil
import warnings
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
    # FIX(35): read a Quantum ESPRESSO input directly rather than letting
    # pymatgen guess from the extension.  Its handler registry ends with a
    # deliberately generic ``("fleur-inpgen", "*.in*")`` catch-all, so every
    # scf-<mpid>.in matches Fleur and Structure.from_file() raises
    # ``ModuleNotFoundError: No module named 'pymatgen.io.fleur'`` -- which is
    # an ImportError, not one of the exceptions caught below, so the QE
    # fallback never ran and `mainprogram data-combine` died.
    if str(filepath).endswith(".in"):
        try:
            return PWInput.from_file(filepath).structure
        except (OSError, ValueError, IndexError, KeyError):
            pass          # not a QE input after all; let the generic path try
    try:
        return Structure.from_file(filepath)
    except (OSError, ValueError, ImportError):
        # FIX(17): PWInput's constructor takes a Structure, not a path -- the
        # old call built a PWInput whose ``structure`` attribute was the file
        # name string, and the ``.structure`` access below then returned that
        # string instead of a Structure.
        # FIX(35): ImportError is caught too -- see above.
        return PWInput.from_file(filepath).structure

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
        fields = line.split()
        if len(fields) < 3:
            warnings.warn("ignoring malformed line {!r} in {}".format(
                line.strip(), input_file), RuntimeWarning)
            continue
        mpid = fields[1]
        compound = fields[2].strip()
        struc_path = f"R{mpid}-{compound}/relax/POSCAR"
        try:
            struc = read_structure_file(struc_path)
        except (OSError, ValueError):
            # FIX(17): the fallback used to be guarded by a bare ``except:``,
            # which also swallowed KeyboardInterrupt.
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
    os.makedirs("filtered_inputs", exist_ok=True)
    with open(output_file, "w") as write_mpid:
        for i, group in enumerate(unique_structures):
            extract_struc = group[0]
            mpid = extract_struc.site_properties['id'][0]
            comp = str(extract_struc.composition).replace(" ", "")
            write_mpid.write(f"v{i+1} {mpid} {comp}\n")
            # ``cp -r R{mpid}*`` went through a shell with the id interpolated
            # into it, and its ``except FileNotFoundError`` could never fire
            # because os.system reports failure through a return code.
            matches = sorted(glob.glob("R{}*".format(mpid)))
            if matches:
                for source in matches:
                    dest = os.path.join("filtered_inputs", os.path.basename(source))
                    if os.path.isdir(source):
                        shutil.copytree(source, dest, dirs_exist_ok=True)
                    else:
                        shutil.copy(source, dest)
                continue
            scf_input = f"scf_dir/scf-{mpid}.in"
            if os.path.isfile(scf_input):
                shutil.copy(scf_input, "filtered_inputs/")
            else:
                warnings.warn(
                    "nothing to copy for {}: neither R{}* nor {}".format(
                        mpid, mpid, scf_input), RuntimeWarning)

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
