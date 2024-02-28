#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D."""
import sys
import spglib as spg
from pymatgen.core import structure
from pymatgen.io.pwscf import PWInput
from ase.io import espresso, vasp
from ase.cell import Cell
from htepc import MpConnect as mpc
def crystal_extract(filename):
    """
    Function to extract lattice parameters from Quantum Espresso (QE) input files using ASE package.

    Parameters:
    -----------
    filename : str
        QE input file or VASP POSCAR.

    Returns:
    --------
    None

    Prints lattice parameters including cell parameters, cell volume, crystal system, and spacegroup.
    """
    if 'scf' in filename:
        cell_inp = PWInput.from_file(filename).structure
    elif 'POSCAR' in filename or 'CONTCAR' in filename:
        cell_inp = structure.Structure.from_file(filename)
    else:
        print("Structure file not found\n")
    vol = cell_inp.lattice.volume
    cellpar = cell_inp.lattice.parameters
    sg_info = cell_inp.get_space_group_info()
    print("************ Printing structural parameters *****************\n")
    print("*         Cell par: {}".format(cellpar))
    print('*         Cell volume: {} '.format(round(vol,4)) + r'$\AA^3$')
    print("*         Spacegroup: {}              ".format(sg_info))
    print("************************************************************\n")

def convex_hull(comp):
    """
    Function to extract energy above hull from the Materials Project database.

    Parameters:
    -----------
    comp : str
        Compound name.

    Returns:
    --------
    tuple
        A tuple containing the energy above hull, formation energy per atom, and magnetic ordering information.
    """
    obj = mpc()
    print(comp)
    obj.setting(comp)
    order = obj.data['ordering']
    return obj.data['energy_above_hull'],obj.data['formation_energy_per_atom'],order
def main():
    """
    Main function to extract properties from Materials Project and relaxed outputs.

    Usage:
    ------
    python script.py filename comp cond

    filename : str
        Name of the QE input file.

    comp : str
        Compound name.

    cond : str
        Condition indicator.

    Returns:
    --------
    None
    """
    filename = sys.argv[1]
    #comp = sys.argv[2]
    #ehull,form_e,order=convex_hull(comp)
    #cond = sys.argv[3]
    #if cond == "initial":
    #    print("\n")
    #    print("******* Extracting properties from materials project and relaxed outputs *********\n")
    #    print("*         FE (eV/atom): {}                 ".format(round(form_e,4)))
    #    print("*         e_above_hull (eV/atom): {}                 ".format(ehull))
    #    print("*         Ordering: {}                   ".format(order) + "\n")
    crystal_extract(filename)
if __name__ == "__main__":
    main()
