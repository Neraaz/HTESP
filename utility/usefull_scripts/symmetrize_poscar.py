# coding: utf-8
import os
import sys
from pymatgen.core import structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
filename=sys.argv[1]

struc = structure.Structure.from_file(filename)
os.system("mv POSCAR POSCAR_old")
struc = SpacegroupAnalyzer(struc,symprec=0.1).get_primitive_standard_structure()
struc.to("POSCAR")
