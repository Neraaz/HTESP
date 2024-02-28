"""Writen by Niraj K. Nepal, Ph.D. (tug11655@temple.edu)"""
import sys
import glob
import pylab
import numpy as np
from matplotlib import pyplot as plt

def plot_epw(filename,typ):
    """
    Function to plot EPW results.
    
    Parameters
    ----------
    filename : str
        File to plot. For gap-imag, Del-T, gap-pade, filename is prefix of the calculation.
    typ : str
        Type of plot.
        
    typ options
    -----------
    decay : Plotting 'decay.H', 'decay.v', etc.
    l_k_pair : Plotting .lambda_k_pairs.
    l_pair : Plotting .lambda_pair.
    a2f : Plotting .a2f file.
    imdelta : Plotting .imag_aniso_*.00 files.
    redelta : Plotting .pade_aniso_*.00 files.
    qdos : Plotting .qdos_*.00 files.
    gap-imag : Plotting .imag_aniso_gap0*.00 files. Distribution of superconducting gap w.r.t temperature T.
    Del-T : Plotting .imag_aniso_*.00 (only corresponding to zero frequency) w.r.t temperature T.
    gap-pade : Plotting .pade_aniso_gap0*.00 files similar to gap-imag.
    """
    if typ == 'decay':
        data = np.loadtxt(filename, skiprows=2)
        plt.plot(data[:,0], data[:,1], 'ro', linestyle='None')
        param = sys.argv[3]
        plt.text(5, data[:,1].max()/2,param,fontsize = 14)
        plt.xlabel(r"R($\AA$)",fontsize=12)
        plt.savefig('{}.png'.format(filename))
        pylab.savefig("{}.pdf".format(filename), format='pdf', bbox_inches='tight')
    elif typ == 'l_k_pair':
        data = np.loadtxt(filename, skiprows=1)
        plt.plot(data[:,0], data[:,1], 'r-')
        plt.xlabel(r"$\lambda_{nk}$",fontsize=25)
        #plt.xlim(-0.5,15)
        plt.yticks([0.0,0.5,1.0],["0", "0.5", "1.0"],fontsize=20)
        plt.xticks([0.0, 1.0, 2.0, 3.0, 4.0], ["0", "1", "2", "3", "4"],fontsize=20)
        plt.ylabel(r"$\rho(\lambda_{nk}$)",fontsize=25)
        plt.savefig('{}.png'.format(filename))
        pylab.savefig("{}.pdf".format(filename), format='pdf', bbox_inches='tight')
    elif typ == 'l_pair':
        data = np.loadtxt(filename, skiprows=1)
        plt.plot(data[:,0], data[:,1], 'r-')
        plt.xticks(fontsize=12)
        plt.xlabel(r"$\lambda_{nk,mk+q}$",fontsize=20)
        plt.ylabel(r"$\rho(\lambda_{nk,mk+q}$)",fontsize=20)
        plt.savefig('{}.png'.format(filename))
        pylab.savefig("{}.pdf".format(filename), format='pdf', bbox_inches='tight')
    elif typ == 'a2f':
        data = np.loadtxt(filename,skiprows=1)
        plt.plot(data[:,0], data[:,1], 'r-')
        plt.plot(data[:,0], data[:,11], 'k--')
        plt.ylim(-0.1,data[:,1].max()*1.2)
        plt.xlabel(r"$\omega (meV)$",fontsize=14)
        plt.ylabel(r"$\alpha^2F(\omega)$",fontsize=16)
        plt.savefig('{}.png'.format(filename))
        pylab.savefig("{}.pdf".format(filename), format='pdf', bbox_inches='tight')
    elif typ == 'imdelta':
        data = np.loadtxt(filename)
        _,axp = plt.subplots()
        plt.plot(data[:,0]*1000, data[:,3]*1000, 'ro',linestyle='None',markersize=2)
        plt.xlim(0,200)
        axp.margins(x=0.60,y=0.60)
        #plt.text(100,1,"T = {} K".format(Tp),fontsize=20)
        plt.xlabel(r"$i\omega (meV)$",fontsize=20)
        plt.ylabel(r"$\Delta_{nk}$ (meV)",fontsize=20)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.savefig('{}.png'.format(filename))
        pylab.savefig("{}.pdf".format(filename), format='pdf', bbox_inches='tight')
    elif typ == 'redelta':
        data = np.loadtxt(filename)
        plt.plot(data[:,0]*1000, data[:,4]*1000, 'r',linestyle='None',marker='.',markersize=1)
        plt.xlim(0,200)
        #plt.title("T = {} K".format(Tp), fontsize = 14)
        plt.ylim(-20,40)
        plt.xlabel(r"$\omega (meV)$",fontsize=14)
        plt.ylabel(r"$\Delta_{nk}$ (meV)",fontsize=16)
        plt.savefig('{}.png'.format(filename))
        pylab.savefig("{}.pdf".format(filename), format='pdf', bbox_inches='tight')
    elif typ == 'qdos':
        data = np.loadtxt(filename,skiprows=1)
        dos_val = 3.5263592194
        plt.plot(data[:,0]*1000, data[:,1]/dos_val, 'r-')
        plt.xlim(0,15)
        plt.ylim(-0.1,data[:,1].max()*1.20/dos_val)
        #plt.title("T = {} K".format(Tp), fontsize = 14)
        plt.xlabel(r"$\omega$ (meV)",fontsize=14)
        plt.ylabel(r"$\frac{N_s (\omega)}{N_F}$",fontsize=16)
        plt.savefig('{}.png'.format(filename))
        pylab.savefig("{}.pdf".format(filename), format='pdf', bbox_inches='tight')
    elif typ == 'gap-imag':
        files = glob.glob('{}.imag_aniso_gap0*.00'.format(filename))
        lenfile = len(files)
        print(lenfile)
        temp = []
        for i in range(lenfile):
            temp.append(float(files[i].split('_')[-1].split('.')[0]))
        temp = sorted(temp)
        #delta = [1.2180840690,0.82932518177,0.44159748048,0.21778346385,0.032814828592]
        for i in range(lenfile):
            if temp[i] < 10:
                data = np.loadtxt('{}.imag_aniso_gap0_00{}.00'.format(filename,int(temp[i])))
            else:
                data = np.loadtxt('{}.imag_aniso_gap0_0{}.00'.format(filename,int(temp[i])))
            #try:
            #    delta.append(data[:,1].mean())
            #except:
            #    delta.append(data[1])
            for j in range(data.shape[0]):
                try:
                    plt.plot([temp[i],data[j,0]], [data[j,1],data[j,1]], 'b-',lw=3)
                except:
                    plt.plot([temp[i],data[0]], [data[1],data[1]], 'b-',lw=3)
                plt.vlines(temp[i],0,3,'k',linestyle='solid')
        plt.xlabel('T (K)',fontsize=25)
        templabel = []
        for tempi in temp:
            templabel.append("{}".format(int(tempi)))
        print(templabel)
        #temp.pop(-2)
        #templabel.pop(-2)
        plt.xticks(temp,templabel,fontsize=20)
        plt.yticks([0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0], ["0", "0.5", "1", "1.5", "2", "2.5", "3"],fontsize=20)
        #plt.plot(temp,delta,'r--',label=r"$\Delta_{nk} (\omega = 0)@T(K)$")
        #plt.legend(fontsize=12)
        plt.ylabel(r"$\Delta_{nk}$ (meV)",fontsize=25)
        plt.savefig('{}-gap.png'.format(filename))
        pylab.savefig("{}-gap.pdf".format(filename), format='pdf', bbox_inches='tight')
    elif typ == 'Del-T':
        files = glob.glob('{}.imag_aniso_*.00'.format(filename))
        file_list = []
        for filei in files:
            if 'gap' not in filei:
                file_list.append(filei)
        files = file_list
        lenfile = len(files)
        print(lenfile)
        temp = []
        for i in range(lenfile):
            temp.append(float(files[i].split('_')[-1].split('.')[0]))
        print(temp)
        temp = sorted(temp)
        #delta = [1.2180840690,0.82932518177,0.44159748048,0.21778346385,0.032814828592]
        data = []
        for i in range(lenfile):
            if temp[i] < 10:
                data.append(np.loadtxt('{}.imag_aniso_00{}.00'.format(filename,int(temp[i])))[0,3])
            else:
                data.append(np.loadtxt('{}.imag_aniso_0{}.00'.format(filename,int(temp[i])))[0,3])
            #try:
            #    delta.append(data[:,1].mean())
            #except:
            #    delta.append(data[1])
        plt.plot(temp, data*1000, color='k',marker='o',markersize=5,lw=0.4)
        plt.xlabel('T (K)',fontsize=14)
        #plt.plot(temp,delta,'r--',label=r"$\Delta_{nk} (\omega = 0)@T(K)$")
        #plt.legend(fontsize=12)
        plt.ylabel(r"$\Delta_{nk}$ (meV)",fontsize=16)
        plt.savefig('{}-delta-T.png'.format(filename))
        pylab.savefig("{}-delta-T.pdf".format(filename), format='pdf', bbox_inches='tight')
    elif typ == 'gap-pade':
        files = glob.glob('{}.pade_aniso_gap0*.00'.format(filename))
        lenfile = len(files)
        temp = []
        for i in range(lenfile):
            temp.append(float(files[i].split('_')[-1].split('.')[0]))
        temp = sorted(temp)
        for i in range(lenfile):
            data = np.loadtxt('{}.pade_aniso_gap0_0{}.00'.format(filename,int(temp[i])))
            for j in range(data.shape[0]):
                try:
                    plt.plot([temp[i],data[j,0]], [data[j,1],data[j,1]], 'b-',lw=3)
                except:
                    plt.plot([temp[i],data[0]], [data[1],data[1]], 'b-',lw=3)
                plt.vlines(temp[i],0,1.6,'k',linestyle='solid')
        plt.xlabel('T (K)',fontsize=14)
        plt.ylabel(r"$\Delta_{nk}$ (meV)",fontsize=16)
        plt.savefig('{}-pade-gap.png'.format(filename))
        pylab.savefig("{}-pade-gap.pdf".format(filename), format='pdf', bbox_inches='tight')
    else:
        print("no option")
if __name__ == "__main__":
    FILENAME = sys.argv[1]
    TYP = sys.argv[2]
    plot_epw(FILENAME,TYP)
