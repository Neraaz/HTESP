#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
""" """
import sys
import os
from fractions import Fraction
from ase.io import espresso
from ase.cell import Cell
def main():
    """
    main function
    """
    comp = sys.argv[1]
    mpid = sys.argv[2]
    if os.path.isfile(comp):
        file_read = espresso.read_espresso_in(comp)
    elif os.path.isfile('inputs/qeinputs/new/{}.cif'.format(mpid)):
        print("cif file present")
        file_read = espresso.read_espresso_in('inputs/qeinputs/new/scf-{}.in'.format(mpid))
    else:
        print("Neither QE input scf.in file nor structure cif file present\n")
    #cell_orig = file_read.cell
    brav_lat = Cell.get_bravais_lattice(file_read)
    print(brav_lat.lattice_system)
    if not brav_lat.lattice_system == 'cubic':
        #y = x.cell.standard_form()
        #Q = y[1]
        #pos = x.get_positions()
        #finalpos = np.matmul(pos,Q)
        #mattransfer = np.array([[1,0,0],[1,1,0],[0,0,1]])
        #pos = np.transpose(x.get_positions())
        pos = file_read.get_scaled_positions()
        #print(pos)
        finalpos = pos
        #finalpos = np.transpose(np.matmul(Q,pos))
        #finalcell=y[0]
        #finalcell = l.tocell()
        finalcell = brav_lat
        #scaled position
        #finalpos = finalcell.scaled_positions(pos)
        #print(finalpos)
        #finalcell=np.matmul(mattransfer,finalcell)
        specieslist = file_read.symbols
        with open("scf_dir/kpoint-{}.dat".format(mpid), "r") as kpt_file:
            kplines = kpt_file.readlines()
        with open("scf_dir/temp.dat", "w") as temp:
            #f.write("ATOMIC_POSITIONS angstrom\n")
            temp.write("ATOMIC_POSITIONS crystal\n")
            for i,_ in enumerate(specieslist):
                temp.write(specieslist[i] + " ")
                temp.write(str(Fraction(str(finalpos[i][0])).limit_denominator(10))+ " ")
                temp.write(str(Fraction(str(finalpos[i][1])).limit_denominator(10)) + " ")
                temp.write(str(Fraction(str(finalpos[i][2])).limit_denominator(10)) + "\n")
            for i in range(2):
                temp.write(kplines[i])
            temp.write("CELL_PARAMETERS angstrom\n")
            for i in range(3):
                temp.write(str(finalcell[i][0]) + " ")
                temp.write(str(finalcell[i][1]) + " " + str(finalcell[i][2]) + "\n")
        with open('hexagon.dat', 'w') as lat_file:
            lat_file.write('x\n')
            lat_file.write('y\n')
    else:
        with open('hexagon.dat', 'w') as lat_file:
            lat_file.write('x\n')
if __name__ == "__main__":
    main()
