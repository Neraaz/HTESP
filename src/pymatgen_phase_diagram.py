#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to compute convex hull"""
import os
import pandas as pand
from pymatgen.io import pwscf
from pymatgen.core import structure
from pymatgen.entries.computed_entries import ComputedStructureEntry as PDEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.phase_diagram import PDPlotter
from htepc import MpConnect
from check_json import config

def plot_phase(data):
    """
    Function to plot convex hull phase diagram using pymatgen.

    This function generates a convex hull phase diagram based on the provided data, which contains
    energies and compositions. It creates a 'convexhull.csv' file with thermodynamic stability
    information.

    Parameters:
    - data : pandas DataFrame
        Dataframe containing energies and compositions. Columns should include 'ID', 'comp', and 'energy'.

    Returns:
    None

    Generates:
    - convexhull.csv : CSV file with thermodynamic stability information.
    - convexhull.pdf : PDF file containing the convex hull phase diagram.

    """
    input_data = config()
    if 'chull_cutoff' not in input_data.keys():
        chull_cutoff = False
    else:
        chull_cutoff = input_data['chull_cutoff']
    entries = []
    for i in range(data.shape[0]):
        if os.path.isfile("R{}-{}/relax/POSCAR".format(data.ID[i],data.comp[i])):
            struc = structure.Structure.from_file("R{}-{}/relax/POSCAR".format(data.ID[i],data.comp[i]))
        elif os.path.isfile("R{}-{}/relax/scf.in".format(data.ID[i],data.comp[i])):
            struc = pwscf.PWInput.from_file("R{}-{}/relax/scf.in".format(data.ID[i],data.comp[i])).structure
        else:
            print("Neither scf.in nor POSCAR present inside R{}-{}/relax/".format(data.ID[i],data.comp[i]) + "\n")
        num_atoms = len(struc)
        pd_entry = PDEntry(struc,data.energy[i]*int(num_atoms))
        entries.append(pd_entry)
    phase_diag = PhaseDiagram(entries)
    with open("convexhull.csv", "w") as write_convexhull:
        write_convexhull.write("ID,total energy (eV),comp,FE(eV/atom),e_above_hull_calc,decomposition\n")
        for i,entry in enumerate(entries):
            write_convexhull.write(data.ID[i] + "," + str(round(entry.energy,4)) + ",")
            write_convexhull.write(entry.composition.formula.replace(' ','') + ",")
            write_convexhull.write(str(round(phase_diag.get_form_energy_per_atom(entry),4)))
            write_convexhull.write("," + str(round(phase_diag.get_e_above_hull(entry),4)) + "\n")
            #write_convexhull.write("," + str(phase_diag.get_phase_separation_energy(entry, stable_only=True)) + "\n")
            #write_convexhull.write("," + str(phase_diag.get_decomposition(entry.composition)) + "\n")
    plotter = PDPlotter(phase_diag,show_unstable=chull_cutoff,backend="matplotlib")
    #plotter = PDPlotter(phase_diag,show_unstable=0.25,backend="matplotlib")
    plt1 = plotter.get_plot()
    plt1.figure.savefig("convexhull.pdf")
    # Uncomment below if you want to extract data from materials project
    # If it exists

    #print("Adding formation energy and energy above hull from MP database for comparison ....\n")
    #data = pand.read_csv('convexhull.csv')
    #data = pand.DataFrame(data)
    #ebh = []
    #form_e = []
    #obj = MpConnect()
    #for i in range(data.shape[0]):
    #    mpid = data.ID[i]
    #    obj.setting(mpid)
    #    ebh.append(round(obj.data['energy_above_hull'],4))
    #    form_e.append(round(obj.data['formation_energy_per_atom'],4))
    #data['e_above_hull_MP'] = ebh
    #data['FE_MP (eV/atom)'] = form_e
    #data.to_csv('convexhull.csv',index=False)
def main():
    """
    Main function to generate and plot convex hull phase diagram.

    This function serves as the entry point to generate a convex hull phase diagram based on
    energy data provided in a CSV file. It reads the energy data from the 'econv.csv' file,
    generates the convex hull phase diagram using the plot_phase function, and prints a message
    indicating the results are stored in 'convexhull.csv' file.

    Parameters:
    None

    Returns:
    None

    """
    energyfile = "econv.csv"
    data = pand.read_csv(energyfile)
    data = pand.DataFrame(data)
    plot_phase(data)
    print("------------------------------------------------------------------------\n")
    print("Results are stored in 'convexhull.csv' file\n")
    print("------------------------------------------------------------------------\n")
if __name__ == "__main__":
    main()
