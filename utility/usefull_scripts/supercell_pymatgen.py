# coding: utf-8
#computing different supercell
import numpy as np
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import structure
from pymatgen.transformations.advanced_transformations import SupercellTransformation
data = structure.Structure.from_file("POSCAR")
sd = SupercellTransformation([[2,0,0],[0,2,0],[0,0,2]])
## supercell of (sqrt(2),sqrt(2),2) size
min_dist = np.sqrt(2)*data.lattice.a
min_dist = np.round(min_dist,2)-0.01
df = sd.from_boundary_distance(data,min_boundary_dist=min_dist,allow_rotation=True,max_atoms=28).apply_transformation(data)
df.make_supercell([1,1,2])
df.to("poscar.vasp",fmt='poscar')
