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
filename_phonon_freq = 'B2Ta4.freq.gp' #phonon dispersion data in gnu format obtained from QE phonon calculations.
filename_scf = 'scf.in'
filename_projection = sys.argv[1] #This file can be obtained using "python projection_phband.py"
#projection_phband.py requires a file "phonproj.in", whose format is provided inside the script, also need "name.eig" file.

bandpath = 'GMKGALHA,LM,KH' #Provide bandpath and total number of kpoints along symmetric-path.
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
ion1 = sys.argv[2]
ion2 = sys.argv[3]
#mode = int(sys.argv[4]) 
proj1 = proj[:nkpt,:]
proj2 = proj[nkpt:,: ]
percmtoTHZ = 33.356    
print(proj1.shape, proj2.shape, nband)
cmap = matplotlib.colors.ListedColormap(['blue', 'm', 'red'])

# We simply plot projection (p) of one ion, and 1-p will corresponds to another ion.
def plot_colourline(x,y,c):
    #c = cm.jet((c-np.min(c))/(np.max(c) - np.min(c)))
    ax = plt.gca()
    for i in np.arange(len(x) - 1):
        if c[i] <= 0.33:
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], lw=1.5, color='blue')
        elif c[i] > 0.33 and c[i] < 0.67:
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], lw=1.5, color='m')
        else:
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], lw=1.5, color='red')
           
    return

for i in range(1,nband):
    print(i)
    x = en
    y = data[:,i]/percmtoTHZ
    c = proj1[:,i-1]
    plot_colourline(x,y,c)

#plt.colorbar(f)
norm = matplotlib.colors.Normalize(vmin=0, vmax=1.0)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
#print(sm)
# For vertical dashed line. Change accordingly.
plt.plot([ks[1], ks[1]], [-1, 25], 'k--', lw=2) 
plt.plot([ks[2], ks[2]], [-1, 25], 'k--', lw=2) 
plt.plot([ks[3], ks[3]], [-1, 25], 'k--', lw=2) 
plt.plot([ks[4], ks[4]], [-1, 25], 'k--', lw=2) 
plt.plot([ks[5], ks[5]], [-1, 25], 'k--', lw=2) 
plt.plot([ks[6], ks[6]], [-1, 25], 'k--', lw=2) 
plt.plot([ks[7], ks[7]], [-1, 25], 'k--', lw=2) 
plt.plot([ks[8], ks[8]], [-1, 25], 'k--', lw=2) 
plt.plot([ks[9], ks[9]], [-1, 25], 'k--', lw=2) 
plt.title("Projection on bands. (0,1) represents ({},{})".format(ion2,ion1), fontsize=10)
plt.colorbar(sm, ticks=np.linspace(0,1.0,3))
plt.xlabel('K-point',fontsize=15)
plt.ylabel(r'$\omega$ (THz)',fontsize=20) 
#Change tickline accordingly
print(ks.shape,en.shape)
plt.xticks(ks, [r'$\Gamma$', 'X', 'M', r'$\Gamma$', 'Z', 'P', 'N', 'Z1', '', 'M|X', 'P'])    
plt.ylim(-6,20)
plt.xlim(ks[0],ks[-1])
plt.savefig('plot.png')
