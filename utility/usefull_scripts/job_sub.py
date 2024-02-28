# coding: utf-8
from ase.io import espresso
import os
import sys
"""
For volume relaxation. provide scf-header-vcrelax.in (either change shape or keep it fix with cell_dofree param)
 and scf.in (fully relaxed ground-state scf.in file). Also turn on fac for scaling cell parameters.

For pressure, use scf.in file. Run with option = "press". Here pressure is kept fixed and volume is relaxed
with all cell params and angles. Turn on fac for pressure in kbar.

"""
cell = espresso.read_espresso_in('scf.in')
c = cell.cell
#fac = [0.90, 0.92, 0.94, 0.96, 0.98, 0.985, 0.99, 0.995, 1.0, 1.005, 1.01, 1.015, 1.02, 1.04, 1.06, 1.08,1.10]
#fac = [0.91, 0.93, 0.95, 0.97]
fac = [-150,-100,-50,50,100,150,200,250,300,350,400]
#fac=[200]
pos = cell.get_scaled_positions()
species = cell.symbols
option=sys.argv[1]
    
for j,fc in enumerate(fac):
    if option == "relax":
        cell = fc*c
        with open("cell-{}".format(fc), "w") as f:
            f.write("ATOMIC_POSITIONS crystal\n")
            for i in range(len(species)):
                f.write(species[i] + " " + str(pos[i][0])+ " " + str(pos[i][1]) + " " + str(pos[i][2]) + "\n")
            f.write("CELL_PARAMETERS angstrom\n")
            for i in range(3):
                f.write(str(cell[i][0]) + " " + str(cell[i][1]) + " " + str(cell[i][2]) + "\n")
        if not os.path.isdir("R{}/".format(j)):
            os.mkdir("R{}".format(j))
    elif option == "rerelax":
        os.system("""sed -n '/Begin final coordinates/,/End final coordinates/p' R{}/scf.out | sed '$d' | sed '1,4d'| sed '5d' > cell-{}""".format(j,fc))     
    elif option == "press":
        if not os.path.isdir("R{}/".format(fc)):
            os.mkdir("R{}".format(fc))
        if os.path.isfile("R{}/scf.out".format(fc)):
            with open("R{}/scf.out".format(fc), "r") as read_scfout:
                lines = read_scfout.readlines()
            for l in lines:
                if "bfgs converged" in l:
                    converged = True
        else:
            converged = False
        if converged:            
            os.system("""sed -n '/Begin final coordinates/,/End final coordinates/p' R{}/scf.out | sed '$d' | sed '1,4d'| sed '5d' > R{}/cell.dat""".format(fc,fc))     
            os.system("""sed -n '/&CONTROL/,/CELL_PARAMETERS/p' R{}/scf.in | sed '$d' > R{}/scf-header.in""".format(fc,fc))
            os.system("""cp R{}/scf.in R{}/scf-1.in""".format(fc,fc))
            os.system("""cp R{}/scf.out R{}/scf-1.out""".format(fc,fc))
            os.system("""cat R{}/scf-header.in R{}/cell.dat > R{}/scf.in""".format(fc,fc,fc))
        else:
            os.system("""sed "/&CELL/a   press = {}," scf.in > R{}/scf.in""".format(fc,fc))
    else:
        print("option = relax/rerelax\n")
    if option == "press":
        os.system("cp run.sh R{}/run-{}.sh".format(fc,j))
        os.chdir("R{}".format(fc))
        os.system("sbatch run-{}.sh".format(j))
        os.system("sleep 2s")
        os.chdir("../")
    else:
        os.system("""cat scf-header-vcrelax.in cell-{} > R{}/scf.in""".format(fc,j))
        os.system("cp run.sh R{}/run-{}.sh".format(j,fc))
        os.chdir("R{}".format(j))
        os.system("sbatch run-{}.sh".format(fc))
        os.system("sleep 2s")
        os.chdir("../")
