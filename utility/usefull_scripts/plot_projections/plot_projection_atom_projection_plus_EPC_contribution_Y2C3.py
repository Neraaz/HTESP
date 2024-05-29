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
filename_phonon_freq = 'C12Y8.freq.gp' #phonon dispersion data in gnu format obtained from QE phonon calculations.
filename_scf = 'scf.in'
filename_projection = sys.argv[1] #This file can be obtained using "python projection_phband.py"
#projection_phband.py requires a file "phonproj.in", whose format is provided inside the script, also need "name.eig" file.

bandpath = 'GHNGPH' #Provide bandpath and total number of kpoints along symmetric-path.
nkpoint = 200


file1 = espresso.read_espresso_in(filename_scf)
band = file1.cell.bandpath(bandpath, npoints=nkpoint)
en = band.get_linear_kpoint_axis()[0]
ks = band.get_linear_kpoint_axis()[1]

with open('lambda.out') as f:
   lines = f.readlines()
x=float(lines[1].split()[11])
print(x)
dosef = float(lines[1].split()[11])/3289.9146
data_gam=np.loadtxt('elgam.dat')/1000
data_gam = data_gam.reshape(200,60)
print(dosef,data_gam.min(),data_gam.max())


data = np.loadtxt(filename_phonon_freq)
#data2 = np.loadtxt('freq.plot')[:,1]/33.356
data3=(data[:,1:]/33.356)**2.0
data3[np.where(data < 0)] *= -1
#np.seterr(invalid='ignore')
lqu=data_gam/(np.pi*dosef*data3)
#lqu[np.where(lqu > 0.1)] = 0
#lqu_min = lqu.max()*0.4
#lqu[np.where(data < 0.0)] = 0
#lqu[np.where(lqu < lqu_min)] = 0
#lqu[np.where(lqu > lqu_min)] = 2
lqu = lqu/lqu.max()
nband = data.shape[1]
nkpt = data.shape[0]
lqu = lqu.reshape(nkpt,nband-1)
#print(lqu)
#fig,ax = plt.subplots(1, 1, sharex=True, sharey=True)
proj = np.loadtxt(filename_projection)
#provide ion names. Mg and O in MgO crystal
ion1 = sys.argv[2]
ion2 = sys.argv[3]
#mode = int(sys.argv[4]) 
proj1 = proj[:nkpt,:]
proj2 = proj[nkpt:,: ]
percmtoTHZ = 33.356    
print(proj1.shape, proj2.shape, nband,lqu.shape)
cmap = matplotlib.colors.ListedColormap(['blue', 'm', 'red'])

# We simply plot projection (p) of one ion, and 1-p will corresponds to another ion.
def plot_colourline(x,y,c,d):
    #c = cm.jet((c-np.min(c))/(np.max(c) - np.min(c)))
    ax = plt.gca()
    for i in np.arange(len(x) - 1):
        if c[i] <= 0.33:
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], lw=0.5, color='blue')
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], linestyle="None", color='springgreen',marker='o', alpha=0.15, markersize=d[i],markeredgewidth=1,markerfacecolor=None,markeredgecolor='springgreen')
        elif c[i] > 0.33 and c[i] < 0.67:
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], lw=0.5, color='m')
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], linestyle="None", color='springgreen', alpha=0.15, marker='o', markersize=d[i],markeredgewidth=1,markerfacecolor=None,markeredgecolor='springgreen')
        else:
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], lw=0.5, color='red')
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], linestyle="None", color='springgreen', marker='o', alpha=0.15, markersize=d[i],markeredgewidth=1,markerfacecolor=None,markeredgecolor='springgreen')
           
    return

for i in range(1,nband):
    print(i)
    x = en
    y = data[:,i]/percmtoTHZ
    c = proj1[:,i-1]
    d = lqu[:,i-1]*15
    plot_colourline(x,y,c,d)

#plt.colorbar(f)
norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
# For vertical dashed line. Change accordingly.
plt.plot([ks[1], ks[1]], [-10, 50], 'k--', lw=0.5) 
plt.plot([ks[2], ks[2]], [-10, 50], 'k--', lw=0.5) 
plt.plot([ks[3], ks[3]], [-10, 50], 'k--', lw=0.5) 
plt.plot([ks[4], ks[4]], [-10, 50], 'k--', lw=0.5) 
plt.title("BCC", fontsize=25)
cbar=plt.colorbar(sm,label="Atomic projection on bands",aspect=50,orientation='horizontal')
cbar.set_ticks([0,0.5,1])
cbar.set_ticklabels(['C','C+Y','Y'], weight='bold', fontsize=10)
#plt.xlabel('K-point',fontsize=15)
plt.ylabel(r'$\omega$ (THz)',fontsize=25) 
#Change tickline accordingly
plt.xticks(ks, [r'$\Gamma$', 'H', 'N', r'$\Gamma$', 'P', 'H'])    
plt.ylim(-10,50)
plt.xlim(ks[0],ks[-1])
plt.savefig('plot-BCC.png')
plt.savefig('plot-BCC.pdf', dpi=600, format='pdf', bbox_inches='tight')
