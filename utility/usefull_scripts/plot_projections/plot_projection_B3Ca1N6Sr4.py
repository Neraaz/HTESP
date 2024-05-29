# coding: utf-8
"""
Python script to plot atomic projection on phonon bands in Quantum Espresso.
Works for 2 atoms solid.
requirements: 
a. scf.in (ground-state input file)
b. phonon dispersion in gnuplot format.
USAGE: python plot_projection.py phonon-name.eig.gp A B

"""
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from ase.io import espresso
filename_scf = 'scf.in'
filename_projection = sys.argv[1] #This file can be obtained using "python projection_phband.py"
filename_phonon_freq = sys.argv[2] #phonon dispersion data in gnu format obtained from QE phonon calculations.
outfile = sys.argv[3]
#projection_phband.py requires a file "phonproj.in", whose format is provided inside the script, also need "name.eig" file.

bandpath = 'GHNGPH,PN' #Provide bandpath and total number of kpoints along symmetric-path.
nkpoint = 200


file1 = espresso.read_espresso_in(filename_scf)
band = file1.cell.bandpath(npoints=nkpoint)
en = band.get_linear_kpoint_axis()[0]
ks = band.get_linear_kpoint_axis()[1]


data = np.loadtxt(filename_phonon_freq)
nband = data.shape[1]
nkpt = data.shape[0]
#fig,ax = plt.subplots(1, 1, sharex=True, sharey=True)
proj = np.loadtxt(filename_projection)
#provide ion names. Mg and O in MgO crystal
#ion1 = sys.argv[2]
#ion2 = sys.argv[3]
#mode = int(sys.argv[4]) 
proj1 = proj[:nkpt,:]
proj2 = proj[nkpt:2*nkpt,: ]
proj3 = proj[2*nkpt:3*nkpt,:]
proj4 = proj[3*nkpt:,:]
percmtoTHZ = 33.356    
print(proj1.shape, proj2.shape, proj3.shape, proj4.shape, nband)
proj5 = np.zeros_like(proj1)
proj5[np.where(proj1 > 0.5)] = 0.25
proj5[np.where(proj2 > 0.5)] = 0.50
proj5[np.where(proj3 > 0.5)] = 0.75
proj5[np.where(proj4 > 0.5)] = 1.0
cmap = matplotlib.colors.ListedColormap(['red', 'blue', 'green', 'cyan', 'k'])
#Al blue, Mo: green, 'red' C

# We simply plot projection (p) of one ion, and 1-p will corresponds to another ion.
def plot_colourline(x,y,c):
    #c = cm.jet((c-np.min(c))/(np.max(c) - np.min(c)))
    ax = plt.gca()
    for i in np.arange(len(x) - 1):
        if abs(c[i] - 0.25) < 0.0000001:
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], lw=1.5, color='red')
        elif abs(c[i] - 0.50) < 0.0000001:
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], lw=1.5, color='blue')
        elif abs(c[i] - 0.75) < 0.0000001:
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], lw=1.5, color='green')
        elif abs(c[i] - 1.0) < 0.0000001:
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], lw=1.5, color='cyan')
        else:
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], lw=1.5, color='k')
           
    return

for i in range(1,nband):
    print(i)
    x = en
    y = data[:,i]/percmtoTHZ
    c = proj5[:,i-1]
    plot_colourline(x,y,c)

#plt.colorbar(f)
norm = matplotlib.colors.Normalize(vmin=0, vmax=1.0)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
maxy = np.max(data)/percmtoTHZ + 2.0
miny = np.min(data)/percmtoTHZ - 1.0
#print(sm)
# For vertical dashed line. Change accordingly.
plt.plot([ks[1], ks[1]], [miny, maxy], 'k--', lw=2) 
plt.plot([ks[2], ks[2]], [miny, maxy], 'k--', lw=2) 
plt.plot([ks[3], ks[3]], [miny, maxy], 'k--', lw=2) 
plt.plot([ks[4], ks[4]], [miny, maxy], 'k--', lw=2) 
plt.plot([ks[5], ks[5]], [miny, maxy], 'k--', lw=2) 
#plt.title("Projection on bands. (0,1) represents ({},{})".format(ion2,ion1), fontsize=10)
#plt.colorbar(sm, ticks=np.linspace(0,1.0,5))
#plt.xlabel('K-point',fontsize=15)
plt.ylabel(r'$\omega$ (THz)',fontsize=20) 
#Change tickline accordingly
print(ks.shape,en.shape)
plt.xticks(ks, [r'$\Gamma$', 'H', 'N', r'$\Gamma$', 'P','', 'H|P', 'N'])    
plt.ylim(miny,maxy)
plt.xlim(ks[0],ks[-1])
plt.savefig(outfile)
