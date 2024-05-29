#!/usr/bin/env python
"""
This script generates band structure and density of states (DOS) plots
for a given material using Quantum ESPRESSO output files.

Usage: python plot_band_dos_QE.py <material_id> <compound_name> <nkpoint>
"""

import sys
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pylab
import os
from ase.io import espresso
from kpath import kpath

# Set font properties for the plots
font = {'weight': 'bold', 'size': 13}
matplotlib.rc('font', **font)
plt.rcParams["figure.autolayout"] = True

def plot(file, comp):
    """
    Plots the band structure and DOS for the given compound.

    Parameters:
    file (str): Path to the input file for the kpath calculation.
    comp (str): Compound name.
    projection (bool): Whether to plot the projection data. Default is False.
    """
    n = int(sys.argv[3])
    kcut = 0
    pt1, sympt, symlb, pt, _, _ = kpath(file, n, kcut)
    pt = np.array(pt)
    sympt = np.array(sympt)
    
    # Read Fermi energy from the DOS file
    with open(f"{comp}.dos", "r") as f:
        fermi = float(f.readlines()[0].split()[8])
    
    # Load DOS data
    datados = np.loadtxt(f"{comp}.dos", skiprows=1)
    datadosx = datados[:, 1]
    datadosy = datados[:, 0] - fermi

    # Load band structure data
    data = np.loadtxt(f'{comp}.dat.gnu')
    fig, ax = plt.subplots(nrows=1, ncols=2, gridspec_kw={'width_ratios': [3, 1], 'wspace': 0.05})
    nband = int(data.shape[0] / n)
    print(f"nkpoints: {n} - nband: {nband}")

    for i in range(nband):
        ax[0].plot(pt, data[0 + n * i:n * (i + 1)][:, 1] - fermi, 'b', lw=1)
    
    # Plot the total DOS
    ax[1].plot(datadosx, datadosy, 'orange', lw=1, label='Total dos')
    ax[1].set_xlabel('DOS', fontsize=15)
    ax[0].set_xticks(sympt)
    ax[0].plot([sympt[0], sympt[-1]], [0, 0], 'k--', lw=0.75)
    ax[1].plot([datadosx.min(), datadosx.max()], [0, 0], 'k--', lw=0.75)
    ax[0].set_xticklabels(symlb)
    
    # Set up colors for the projection plot
    color = ['k', 'r', 'b', 'g', 'cyan', 'lightgreen', 'orange', 'yellow', 'lightblue']

    # Read projected DOS data
    try:
        with open("filedos.in", "r") as p_dos:
            lines = p_dos.readlines()
    except FileNotFoundError:
        print("filedos.in not found\n")
        sys.exit()

    # Parse the projected DOS data
    dict_ = {}
    for line in lines:
        key, *values = line.strip().split()
        dict_[key] = values

    datalist = []
    datalist_name = []
    for key in dict_.keys():
        for value in dict_[key]:
            data = np.loadtxt(f"{key}-{value}.dat")
            print(data.shape)
            datalist.append(data)
            datalist_name.append(f'{key}-{value}')

    en_l = datalist[0][:, 0]
    for i, data in enumerate(datalist):
        ax[1].plot(data[:, 1], en_l - fermi, color=color[i], label=datalist_name[i], lw=1)

    if len(datalist) < 5:
        ax[1].legend(loc="best", frameon=False, fontsize=10, handlelength=1.0, bbox_to_anchor=(0.1, 0.65))
    else:
        ax[1].legend(ncol=2, loc="best", frameon=False, fontsize=10, handlelength=1.0, bbox_to_anchor=(0.1, 0.65))
    
    # Add vertical lines at high-symmetry points
    for i in range(sympt.shape[0]):
        ax[0].vlines(sympt[i], -4.0, 12.0, color='k', linestyles='dashed')
    
    # Set y-axis labels and limits
    ax[0].set_ylabel("E - Ef (eV)")
    ax[0].set_ylim((-5, 5)) # Adjust this depending on the energy levels
    ax[1].set_ylim((-5, 5)) # Adjust this depending on the energy levels

    # Set x-axis limits for DOS plot
    indices1 = np.where((datadosy >= -3.5) & (datadosy <= 8.5))
    dosmax = (np.max(datadosx[indices1])) + 1
    ax[1].set_xlim((0, dosmax)) # Adjust this depending on the energy levels

    # Remove tick labels from the DOS plot
    ax[1].set_xticklabels([])
    ax[1].set_yticklabels([])

    # Save the plot
    plt.savefig(comp + "-band.png")
    pylab.savefig(comp + '-band.pdf', format='pdf', bbox_inches='tight')

if __name__ == "__main__":
    print("Usage: python plot_band_dos_QE.py <material_id> <compound_name> <nkpoint>\n")
    mpid = sys.argv[1]
    comp = sys.argv[2]
    filename = 'scf.in'
    plot(filename, comp)
