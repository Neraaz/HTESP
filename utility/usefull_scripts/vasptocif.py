"""Writen by Niraj K. Nepal, Ph.D. (tug11655@temple.edu)"""
import sys
from pymatgen.io import cif
from pymatgen.io.vasp import inputs
filename = sys.argv[1]
data = inputs.Poscar.from_file(filename)
structure = data.structure
file = cif.CifWriter(structure)
file.write_file('relax.cif')
