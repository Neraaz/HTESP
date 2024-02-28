# coding: utf-8
import numpy as np
from ase.io import vasp
data = vasp.read_vasp('POSCAR')
bandpath = data.cell.bandpath(npoints=300)
data1 = np.loadtxt("band.dat")
kpt = bandpath.get_linear_kpoint_axis()[0]
data2 = data1
for i in range(data2.shape[0]):
    data2[i][0] = np.round(kpt[i],4)
    
with open("data-band.txt", "w") as f:
    for i in range(data2.shape[0]):
        d = data2[i]
        for j in range(d.shape[0]):
            if j < d.shape[0] - 1:
                f.write(str(d[j]) + " ")
            else:
                f.write(str(d[j]) + "\n")
                
kspecial_point = bandpath.get_linear_kpoint_axis()[1]
kspecial_name = bandpath.get_linear_kpoint_axis()[2]
with open("kspecial.txt", "w") as g:
    for i in range(len(kspecial_point)):
        g.write(kspecial_name[i] + " " + str(round(kspecial_point[i],4)) + "\n")
