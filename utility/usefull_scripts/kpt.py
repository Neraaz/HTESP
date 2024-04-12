"""Writen by Niraj K. Nepal, Ph.D. (tug11655@temple.edu)
This script prints K-mesh over high-symmetry BZ path,
and performs cartesian to crystal coordinate conversion and vice-
versa.
"""
import sys
from ase.io import espresso, vasp

filename = sys.argv[1]
option = sys.argv[2]
dft = sys.argv[3]
if dft == 'QE':
    data = espresso.read_espresso_in(filename)
else:
    data = vasp.read_vasp(filename)

if option == "bandpath":
    path = sys.argv[3]
    band = data.cell.bandpath(path='GXSYGZURTZ', npoints=400)
    kpt = band.kpts
    with open("kpt.dat", "w") as f:
        for i in range(kpt.shape[0]):
            f.write(str(kpt[i][0]) + " ")
            f.write(str(kpt[i][1]) + " ")
            f.write(str(kpt[i][2]) + " " + str(0) +  "\n")

elif option == "cart2crys":
    pos = data.get_scaled_positions()
    sym = data.symbols
    with open("pos-cart-2-crys.dat", "w") as f:
        for i in range(pos.shape[0]):
            f.write(sym[i] + " " + str(pos[i][0]))
            f.write(" " + str(pos[i][1]) + " " + str(pos[i][2]) +  "\n")
elif option == "crys2cart":
    pos = data.positions
    sym = data.symbols
    with open("pos-crys-2-cart.dat", "w") as f:
        for i in range(pos.shape[0]):
            if dft == 'QE':
                f.write(sym[i] + " " + str(pos[i][0]))
                f.write(" " + str(pos[i][1]) + " " + str(pos[i][2]) +  "\n")
            else:
                f.write(str(pos[i][0]) + " " + str(pos[i][1]) + " " + str(pos[i][2]) + " " + sym[i] +  "\n")
else:
    print("option is either 'bandpath', crys2cart, or 'cart2crys'\n")
