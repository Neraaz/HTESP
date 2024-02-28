"""Writen by Niraj K. Nepal, Ph.D. (tug11655@temple.edu)"""
"""
To manually obtain CDW, find q for which the phonon instability occurs for acoustic mode
Create a supercell according to that q.
Distort the unit cell part of the supercell according to eigenmode displacement and find the
ground-state by using all the eigenmode for the distortion.
Try to find symmetric structure out of it after relaxation. Calculate phonon dispersion.
This script create a supercell and apply distortion to its unit cell.
parameters
--------------------
file1 ==> undistorted structure
file2 ==> distorted structure (obtained either from QE .axsf file from DFPT or VASP+Phonopy band.yaml)
For QE, use src/qe_axsf2cellpos.py scripts, while for phonopy band.yaml use apply_displacement.py script.
Also provide dimension of supercell along each direction.
Usage: python cdw.py file1 file2 Nx Ny Nz output
This script simply computes distances of each sites in distorted structure and replace them in supercell
with minimum distances. Please check the final structure before proceeding further.
"""
import sys
from ase.io import espresso
from pymatgen.io import pwscf
from pymatgen.core import structure
import numpy as np
from scipy.linalg import norm

if __name__ == "__main__":
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    N1 = int(sys.argv[3])
    N2 = int(sys.argv[4])
    N3 = int(sys.argv[5])
    outfile = sys.argv[6]
    try:
        STRUC = structure.Structure.from_file(file1)
    except:
        STRUC = pwscf.PWInput.from_file(file1).structure
    STRUC.make_supercell([N1,N2,N3])
    STRUC.to("relax-{}{}{}.cif".format(N1,N2,N3),fmt="cif")
    try:
        STRUCB = structure.Structure.from_file(file2)
    except:
        STRUCB = pwscf.PWInput.from_file(file2).structure
    N = len(STRUC.sites)
    M = len(STRUCB.sites)
    super_i = []
    change_i = []
    for i in range(M):
        dist = np.zeros(N)
        for j in range(N):
            if STRUCB.sites[i].specie == STRUC.sites[j].specie and j not in super_i and norm(STRUC.sites[j].coords) < 5.42:
                dist[j] += norm(STRUCB.sites[i].coords-STRUC.sites[j].coords)
            else:
                dist[j] += 100
        try:
            IDX = int(np.where(dist == dist.min())[0])
            change_i.append(i)
            super_i.append(IDX)
            print(i,IDX)
        except IndexError:
            print("{} ion not found\n".format(M-len(super_i)))
            print("{} ion not found\n".format(i))
    ind = np.unique(super_i)
    if not  ind.shape[0] == M:
        print("Duplication detected\n")
    for i,_ in enumerate(change_i):
        STRUC.sites[super_i[i]].coords = STRUCB.sites[change_i[i]].coords
    cell = STRUC.lattice.matrix
    pos1 = STRUC.frac_coords
    pos2 = STRUC.cart_coords
    pos1 = pos2
    with open("cdw-supercell-{}{}{}.in".format(N1,N2,N3), "w") as g:
        g.write("ATOMIC_POSITIONS crystal\n")
        for j in range(pos1.shape[0]):
            g.write(STRUC.species[j].name + " ")
            g.write(str(pos1[j][0]) + " " + str(pos1[j][1]))
            g.write(" " + str(pos1[j][2]) + "\n")
        g.write("CELL_PARAMETERS angstrom\n")
        for i in range(3):
            g.write(str(cell[i][0]) + " " + str(cell[i][1]) + " " +  str(cell[i][2]) + "\n")
    STRUC.to(outfile,fmt="poscar")
