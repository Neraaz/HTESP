#!/usr/bin/env python
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

font = {'weight' : 'bold','size'   : 13}
matplotlib.rc('font', **font)
plt.rcParams["figure.autolayout"] = True
def plot(file,comp,projection=False): 
    n = int(sys.argv[3])
    kcut = int(sys.argv[4])
    pt1,sympt,symlb,pt,_,_ = kpath(file,n,kcut)
    pt = np.array(pt)
    sympt = np.array(sympt)
    with open("{}.dos".format(comp), "r") as f:
        fermi = float(f.readlines()[0].split()[8])
    datados = np.loadtxt("{}.dos".format(comp), skiprows=1)
    datadosx = datados[:,1]
    datadosy = datados[:,0] - fermi
    data=np.loadtxt('{}.dat.gnu'.format(comp))
    fig,ax = plt.subplots(nrows=1,ncols=2,gridspec_kw={'width_ratios': [3, 1],'wspace': 0.05})
    nband=int(data.shape[0]/n)
    print("nkpoints:{}-nband:{}".format(n,nband))
    if not projection:
        for i in range(nband):
            ax[0].plot(pt, data[0+n*i:n*(i+1)][:,1]-fermi,'b', lw = 1)
    ax[1].plot(datadosx,datadosy, 'orange', lw=1, label='Total dos')
    ax[1].set_xlabel('DOS', fontsize=15)
    ax[0].set_xticks(sympt)
    ax[0].plot([sympt[0], sympt[-1]], [0, 0], 'k--', lw=0.75)
    ax[1].plot([datadosx.min(),datadosx.max()], [0, 0], 'k--', lw=0.75)
    
    ax[0].set_xticklabels(symlb)
    color = ['k', 'r', 'b', 'g','cyan','lightgreen','orange','yellow','lightblue']
    try:
        with open("filedos.in", "r") as p_dos:
            lines = p_dos.readlines()
    except FileNotFoundError:
        print("filedos.in not found\n")
        sys.exit()
    dict_ = {}
    for line in lines:
        dict_[line.split('\n')[0].split(' ')[0]] = line.split('\n')[0].split(' ')[1:]
    datalist = []
    datalist_name = []
    for key in dict_.keys():
        len_dict = len(dict_[key])
        for i in range(len_dict):
            #os.system("sumpdos.sh {} {}".format(key,dict_[key][i]))
            data = np.loadtxt("{}-{}.dat".format(key,dict_[key][i]))
            print(data.shape)
            datalist.append(data)
            datalist_name.append('{}-{}'.format(key,dict_[key][i]))
    en_l = datalist[0][:,0]
    for i,_ in enumerate(datalist):
        ax[1].plot(datalist[i][:,1],en_l - fermi,color=color[i],label=datalist_name[i],lw=1)
    if len(datalist) < 5:
        ax[1].legend(loc="best",frameon=False,fontsize=10,handlelength=1.0,bbox_to_anchor=(0.1, 0.65))
    else:
        ax[1].legend(ncol=2,loc="best",frameon=False,fontsize=10,handlelength=1.0,bbox_to_anchor=(0.1, 0.65))
    for i in range(sympt.shape[0]):
        ax[0].vlines(sympt[i],-4.0,12.0,color='k',linestyles='dashed')
    ax[0].set_ylabel("E - Ef (eV)")
    ax[0].set_ylim((-5, 5)) #change this depending on the energy levels.
    ax[1].set_ylim((-5, 5)) #change this depending on the energy levels.
    indices1 = np.where(np.logical_and(datadosy>=-3.5, datadosy<=8.5))
    dosmax= (np.max(datadosx[indices1])) + 1
    ax[1].set_xlim((0, dosmax)) #change this depending on the energy levels.
    ax[1].set_xticklabels([])
    ax[1].set_yticklabels([])
    #ax[0].set_title(comp)
    #dos_plot(fermi,ax[1])
    plt.savefig(comp + "-band.png")
    pylab.savefig(comp + '-band.pdf', format='pdf',bbox_inches='tight')
if __name__ == "__main__":
    mpid = sys.argv[1]
    comp = sys.argv[2]
    filename = 'scf.in'
    plot(filename,comp,projection=False)
