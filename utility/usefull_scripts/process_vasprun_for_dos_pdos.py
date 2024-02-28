# coding: utf-8
from pymatgen.io.vasp.outputs import Vasprun
data = Vasprun("vasprun.xml",parse_potcar_file=False)
fermi = data.efermi
element_list = data.complete_dos.get_element_dos().keys()
element_list = list(element_list)
with open("atom_dos.txt", "w") as f:
    for i,atom in enumerate(element_list):
        f.write("#" + str(atom) + ", E-fermi: " + str(fermi) + "\n")
        energies = data.complete_dos.get_element_dos()[element_list[i]].energies
        dos = data.complete_dos.get_element_dos()[element_list[i]].get_densities()
        for en,d in zip(energies,dos):
            f.write(str(en) + " " + str(d) + "\n")
            
spd_list = list(data.complete_dos.get_element_spd_dos(element_list[0]).keys())
with open("atom_dos_orb.txt", "w") as f:
    for i,atom in enumerate(element_list):
        f.write("#" + str(atom) + ", E-fermi: " + str(fermi) + " ")
        spd_list = list(data.complete_dos.get_element_spd_dos(element_list[i]).keys())
        dos = []
        energies = data.complete_dos.get_element_spd_dos(element_list[i])[spd_list[0]].energies
        for j in range(len(spd_list)):
            f.write(spd_list[j].name + " ")
            dos.append(data.complete_dos.get_element_spd_dos(element_list[i])[spd_list[j]].get_densities())
        f.write("\n")
        print(dos[0][0],dos[1][0],dos[2][0])
        for i in range(energies.shape[0]):
            f.write(str(energies[i]) + " ")
            for d in dos:
                f.write(str(d[i]) + " ")
            f.write("\n")
