# coding: utf-8
"""
NKN
"""
from matplotlib import pyplot as plt
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import DosPlotter
plotter = DosPlotter()
data = Vasprun("vasprun.xml",parse_potcar_file=False)
totaldos = data.complete_dos
atomdos = totaldos.get_element_dos()
ELEMENTLIST = list(totaldos.get_element_dos().keys())
for elm in ELEMENTLIST:
    ELEMENT = str(elm)
    plotter.add_dos(ELEMENT,atomdos[elm])
plotter.add_dos("total DOS",totaldos)
ax=plotter.get_plot(xlim=(-4,4),ylim=(0,100))
ax.ylabel("Density of states (states/eV)",fontweight='bold',fontsize=25)
ax.xlabel(r"Energy - E$_F$ (eV)",fontweight='bold',fontsize=25)
plt.savefig("new.pdf")
