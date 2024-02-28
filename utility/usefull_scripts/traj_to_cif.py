"""
extracting cif files from optimal population
"""
# coding: utf-8
import sys
from ase.io.trajectory import Trajectory
from ase.io import cif
filename = sys.argv[1]
data = Trajectory(filename)
for i,d in enumerate(data):
    cif.write_cif("mp-{}.cif".format(i+1),d)
