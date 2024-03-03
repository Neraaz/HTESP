#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D."""
import sys
import os
import warnings
import json
import matplotlib
import matplotlib.pyplot as plt
import pylab
import numpy as np
from pymatgen.io.vasp.inputs import Kpoints
from pymatgen.io.vasp.sets import Vasprun
from pymatgen.core import Composition
from pymatgen.electronic_structure.core import OrbitalType
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.electronic_structure.plotter import BSPlotter
from kpath import kpath
warnings.filterwarnings("ignore")
matplotlib.use('Agg')

font = {'weight' : 'bold','size'   : 20}
matplotlib.rc('font', **font)
plt.rcParams["figure.autolayout"] = True
try:
    PWD = os.getcwd()
    if os.path.isfile(PWD+"/htepc.json"):
        JSONFILE = PWD+"/htepc.json"
    else:
        JSONFILE = "../../htepc.json"
    with open(JSONFILE, "r") as readjson:
        input_data = json.load(readjson)
except FileNotFoundError:
    print("htepc.json file not found\n")

def kptline():
    """
    Generate k-path for VASP-line format.

    Reads KPOINTS file and vasprun.xml file to generate the k-path for the VASP-line format.
    Returns:
    - nkpt (numpy.ndarray): Array containing distances between k-points.
    - sympoint (numpy.ndarray): Array containing the special k-points along the path.
    - symname (list): List of special k-point names formatted for LaTeX.

    This function computes the k-path for the VASP-line format based on the KPOINTS file and vasprun.xml file.
    It calculates the distances between k-points and identifies special k-points along the path.

    Example usage:
    nkpt, sympoint, symname = kptline()
    """
    kpt = Kpoints.from_file("KPOINTS")
    data = Vasprun("vasprun.xml")
    nkpt = np.array(data.get_band_structure("KPOINTS").distance)
    n = []
    npair = int(len(kpt.labels)/2)
    nline = int(kpt.num_kpts)
    for i in range(npair):
        n.append(nline*i)
        n.append(nline*(i+1)-1)
    sympoint = nkpt[n]
    symname = kpt.labels
    for i in range(sympoint.shape[0]):
        if sympoint[i-1] == sympoint[i]:
            if symname[i] != symname[i-1]:
                symadd = symname[i-1] + "|" + symname[i]
                symname[i-1] = ''
            else:
                symadd = symname[i]
            symname[i] = r'${}$'.format(symadd)
    symname[0] = r'${}$'.format("\\Gamma")
    return nkpt,sympoint,symname
def plot(plottype,file,comp):
    """
    Function to create quick plots.

    Parameters:
    -----------
    plottype : str
        Type of plot. 'band' for electronic bandstructure and DOS side by side,
        'phonband' for phonon bandstructure, 'gammaband' for phonon band with lambda projection,
        'a2f' for Eliashberg spectral function.
    file : str
        QE scf input file to generate k-point mesh.
    comp : str
        Compound name.

    Returns:
    ---------
    None

    This function creates various types of plots based on the input parameters.
    """
    if plottype == 'band':
        nkpoint = int(sys.argv[4])
        kcut = int(sys.argv[5])
        with open("../../input.in","r") as read_inputin:
            inputline = read_inputin.readlines()
        vasp_line = False
        for line in inputline:
            if "vasp-line" in line:
                vasp_line = True
        if vasp_line:
            pt_l,sympt,symlb = kptline()
            nkpoint = pt_l.shape[0]
        else:
            _,sympt,symlb,pt_l,_,_ = kpath(file,nkpoint,kcut)
        pt_l = np.array(pt_l)
        if os.path.isfile('scf.out'):
            os.system("grep Fermi scf.out | awk '{ print $5 }' > fermi.dat")
            os.system("""grep "number of electrons       =     " scf.out | tail -n 1 | awk '{print $5}' > band_fermi.dat""")
            with open("fermi.dat", "r") as read_fermi:
                fermi = float(read_fermi.readlines()[0].split("\n")[0])
        elif os.path.isfile('../relax/OUTCAR'):
            os.system("""grep NELECT ../relax/OUTCAR | awk '{print $3}' > band_fermi.dat""")
            os.system("grep 'E-fermi' ../relax/OUTCAR | awk '{print $3}' > fermi.dat")
            fermi = np.loadtxt('fermi.dat').item()
        else:
            print("No scf.out and OUTCAR files present\n")
        band_fermi = int(np.loadtxt("band_fermi.dat").item())
        os.system("""grep LSORBIT INCAR | wc -l > lsorbit""")
        lsorbit = float(np.loadtxt("lsorbit"))
        #lsorbit = 0
        if lsorbit > 0:
            band_fermi = band_fermi
        else:
            if not band_fermi % 2 == 0:
                band_fermi = int((band_fermi + 1)/2)
            else:
                band_fermi = int(band_fermi/2)
        if os.path.isfile("band_fermi.dat"):
            os.system("rm band_fermi.dat")
        sympt = np.array(sympt)
        data=np.loadtxt('{}.dat.gnu'.format(comp))
        fig,ax_p = plt.subplots()
        nband=int(data.shape[0]/nkpoint)
        print("nkpoints:{}-nband:{}".format(nkpoint,nband))
        with open("band_stat.csv", "w") as band_stat:
            band_stat.write("Band,Emin-Ef,Emax-Ef,Emin,Emax\n")
            for i in range(nband):
                if i == band_fermi - 1:
                #ax_p.plot(pt_l, data[0+nkpoint*i:nkpoint*(i+1)][:,1]-fermi,'k', lw = 1)
                    ax_p.plot(pt_l, data[0+nkpoint*i:nkpoint*(i+1)][:,1]-fermi,'r', lw = 1.0)
                else:
                    ax_p.plot(pt_l, data[0+nkpoint*i:nkpoint*(i+1)][:,1]-fermi,'k', lw = 0.5)
                band_stat.write(str(i) + ",")
                band_stat.write(str(np.min(data[0+nkpoint*i:nkpoint*(i+1)][:,1])-fermi) + ",")
                band_stat.write(str(np.max(data[0+nkpoint*i:nkpoint*(i+1)][:,1])-fermi) + ",")
                band_stat.write(str(np.min(data[0+nkpoint*i:nkpoint*(i+1)][:,1])) + ",")
                band_stat.write(str(np.max(data[0+nkpoint*i:nkpoint*(i+1)][:,1])) + "\n")
        ax_p.set_xticks(sympt)
        #ax_p.plot([sympt[0], sympt[-1]], [0, 0], 'k--', lw=0.75)
        ax_p.plot([sympt[0], sympt[-1]], [0, 0], 'k--', lw=0.1)
        ax_p.set_xticklabels(symlb)
        min_band = np.min(data[:,1])-fermi
        max_band = np.max(data[:,1])-fermi
        for i in range(sympt.shape[0]):
            ax_p.vlines(sympt[i],min_band-2,max_band+2,color='k',linestyles='dashed')
        ax_p.set_ylabel("Energy - Ef (eV)")
        #ax_p.set_ylim((min_band-2, max_band+2)) #change this depending on the energy levels.
        ylim = input_data['plot']['ylim']
        if ylim is not None:
            ax_p.set_ylim((ylim[0], ylim[1])) #change this depending on the energy levels.
        else:
            ax_p.set_ylim((-1, 1)) #change this depending on the energy levels.
        fig.set_figheight(8)
        fig.set_figwidth(12)
        #ax_p.set_ylim((-6, 4)) #change this depending on the energy levels.
        ax_p.set_title(comp)
        plt.savefig(comp + "-band.png")
        pylab.savefig(comp + '-band.pdf', format='pdf',bbox_inches='tight')

    if plottype == 'phonband':
        nkpoint = int(sys.argv[4])
        kcut = int(sys.argv[5])
        _,sympt,symlb,pt_l,_,_ = kpath(file,nkpoint,kcut)
        pt_l = np.array(pt_l)
        sympt = np.array(sympt)
        data=np.loadtxt('freq.plot')
        #dataph=np.loadtxt('phonon.dos')
        nband=int(data.shape[0]/nkpoint)
        print("nkpoints:{}-nband:{}".format(nkpoint,nband))
        percmtothz = 33.356
        #dataphx=dataph[:,0]/percmtothz
        fig,ax_p = plt.subplots()
        #color = ['r', 'b', 'g', 'k', 'm', 'y', 'c']
        #markers = ['o', 's', '+', '^', '*', 'D', 'p']
        for i in range(nband):
            ax_p.plot(pt_l, data[0+nkpoint*i:nkpoint*(i+1)][:,1]/percmtothz,'k', lw = 1)
        ax_p.set_ylabel(r'$\omega$'+'(THz)',fontsize=25)
        ax_p.set_xticks(sympt)
        ax_p.set_xticklabels(symlb)
        maxy = np.max(data[:,1])/percmtothz + 2.0
        miny = np.min(data[:,1])/percmtothz - 2.0
        for i in range(sympt.shape[0]):
            ax_p.vlines(sympt[i],0,maxy,color='k',linestyles='dashed')
        maxy = np.max(data[:,1])/percmtothz + 2.0
        miny = np.min(data[:,1])/percmtothz - 1.0
        ax_p.set_ylim((miny,maxy))
        ax_p.set_title(comp)
        #plt.savefig(comp+ "-phonon.png")
        pylab.savefig(comp + '-phonon.pdf', format='pdf',bbox_inches='tight')

    if plottype == 'gammaband':
        with open('lambda.dat') as lambd:
            lines = lambd.readlines()
        dosef=float(lines[2].split('\n')[0].split(' ')[-1])/3289.9146
        #dosef = float(lines[1].split()[11])/3289.9146
        nkpoint = int(sys.argv[4])
        kcut = int(sys.argv[5])
        _,sympt,symlb,pt_l,_,_ = kpath(file,nkpoint,kcut)
        pt_l = np.array(pt_l)
        sympt = np.array(sympt)
        data=np.loadtxt('gamma.plot')[:,1]/1000
        data2=np.loadtxt('freq.plot')[:,1]/33.356
        data3=data2**2.0
        #data3[np.where(data2 < -1)] *= -1
        np.seterr(invalid='ignore')
        lqu=data/(np.pi*dosef*data3)
        #lqu=np.nan_to_num(lqu)
        #lqu_norm = np.linalg.norm(lqu)
        lqu[np.where(data2 < 0.6)] = 0
        #lqu_min = lqu.max()*0.2
        #lqu = lqu/lqu.max()
        lqu_min = lqu.max()*0.1
        #print(dosef)
        #lqu = data2
        wph=data2
        nband=int(data.shape[0]/nkpoint)
        #if nband*n < data.shape[0]:
        #    n = 60
        #    nband=int(data.shape[0]/n)
        #    _,sympt,symlb,pt,_,_ = kpath(file,n,kcut)
        #    pt = np.array(pt)
        #    sympt = np.array(sympt)
        print("nkpoints:{}-nband:{}".format(nkpoint,nband))
        percmtothz = 33.356
        fig,ax_p = plt.subplots()
        for i in range(nband):
            x_data = pt_l
            y_data = wph[0+nkpoint*i:nkpoint*(i+1)]
            c_data = lqu[0+nkpoint*i:nkpoint*(i+1)]
            plot_gamma(x_data,y_data,c_data,lqu_min)
        #plt.xlabel("k-vector",fontsize=15)
        ax_p.set_ylabel(r'$\omega$'+'(THz)',fontsize=25)
        ax_p.set_xticks(sympt)
        ax_p.set_xticklabels(symlb)
        maxy = np.max(wph) + 1.0
        miny = np.min(wph) - 0.5
        for i in range(sympt.shape[0]):
            ax_p.vlines(sympt[i],0,maxy,color='k',linestyles='dashed')
        ax_p.set_ylim((miny,maxy))
        ax_p.set_title(comp)
        plt.savefig(comp+ "-gamma.png")
        pylab.savefig(comp + '-gamma.pdf', format='pdf',bbox_inches='tight')
    if plottype == 'a2f':
        fig,ax_p = plt.subplots()
        with open('lambda.out', 'r') as lambd:
            lines = lambd.readlines()
        values = lines[12]
        lam = float(values.split()[0])
        values = values.split()[1:]
        #omglog = float(values[0])
        #if "*" not in values[1]:
        #    Tc = float(values[1])
        #else:
        #    Tc = 0.0
        data = np.loadtxt('alpha2F.dat')
        en_data = data[:,0]
        a2f = data[:,2]
        delw = en_data[1]
        en_data[np.where(en_data == 0.0)] = 0.0000000001
        plt.plot(en_data,a2f,'g',lw=2)
        a2fint = 2.0*delw*a2f/en_data
        a2fint = a2fint.cumsum()
        #lmbd = round(a2fint[-1],4)
        plt.plot(en_data,a2fint,'k--', lw=2)
        if np.max(a2f) > 2.0:
            maxy = np.max(a2f) + 1.0
        else:
            maxy = np.max(a2f) + 0.5
        maxx = np.max(en_data) + 2
        miny = np.min(a2f) - 0.1
        #for tick in ax_p.yaxis.get_major_ticks():
        #    tick.label.set_fontsize(20)
        #for tick in ax_p.yaxis.get_major_ticks():
        #    tick.label.set_fontsize(20)
        ax_p.tick_params(axis='both',labelsize=20)
        fig.set_figheight(8)
        fig.set_figwidth(18)
        plt.ylim((miny,maxy))
        plt.xlim((0.0,maxx))
        plt.xlabel(r'$\omega$' + '(THz)', fontsize=15)
        plt.ylabel(r'$\alpha^2F$'+r'($\omega$)', fontsize=20)
        plt.text(maxx-20,lam+0.5, r'$\lambda$' + "=" + str(round(lam,2)), fontsize=20,fontweight='bold')
        #plt.text(maxx-20,maxy-0.1,r'$\omega_{log}$' + "=" + str(int(omglog)) + " K", fontsize=20,fontweight='bold')
        #plt.text(maxx-20,maxy-1.5, r'$T_c$' + "=" + str(round(Tc,2))+" K", fontsize=20,fontweight='bold')
        plt.savefig(comp+"-a2f.png")
        pylab.savefig(comp+ "-a2f.pdf", format='pdf', bbox_inches='tight')
    if plottype == '':
        print("plot all")
def plot_gamma(xdata,ydata,color,min_):
    """
    Function to plot different sections with respect to different colors.

    Parameters:
    -----------
    xdata : array-like
        x-axis data.
    ydata : array-like
        y-axis data.
    color : array-like
        Projection data.
    min_ : float
        Minimum threshold for binary classification.

    Returns:
    --------
    None

    This function plots different sections of the data with respect to different colors.
    If the color value is less than or equal to the minimum threshold, the section is plotted in black.
    Otherwise, the section is plotted in green with marker sizes proportional to the color values.

    """
    ax_p = plt.gca()
    for i in np.arange(len(xdata) - 1):
        marker_size = color[i]*3
        if color[i] <= min_:
            ax_p.plot([xdata[i],xdata[i+1]], [ydata[i], ydata[i+1]], lw=1.5, color='k')
        else:
            ax_p.plot([xdata[i],xdata[i+1]], [ydata[i], ydata[i+1]], lw=1.5, color='k')
            ax_p.plot([xdata[i],xdata[i+1]], [ydata[i], ydata[i+1]],linestyle='none',marker='o',color='green',markersize=marker_size,fillstyle='none')
def dos_plot(filedos,out='pdos.pdf'):
    """
    Function to plot density of states and partial density of states in different rows.

    Parameters:
    -----------
    filedos : str
        DOS file in .dos format obtained from QE calculations.
    out : str, optional
        Output plot in PDF format. Default is 'pdos.pdf'.

    Returns:
    --------
    None

    This function reads the DOS file and performs DOS and PDOS calculations.
    It requires a 'filedos.in' with ion and orbital contribution from different lines.
    It also uses 'sumpdos.sh' scripts as "sumpdos.sh element orbital" to create 'element-orbital.dat' files.
    For example, for boron and s orbital, element = B, orbital = s, B-s.dat file is created.
    The function then plots the density of states and partial density of states in different rows.

    """
    color = ['k', 'r', 'b', 'g','cyan','lightgreen','orange','yellow','lightblue']
    with open(filedos, "r") as dos:
        fermi = float(dos.readlines()[0].split()[8])
    #fermi = 12.3316
    try:
        with open("filedos.in", "r") as p_dos:
            lines = p_dos.readlines()
    except FileNotFoundError:
        print("filedos.in not found\n")
        sys.exit()
    dos = np.loadtxt(filedos)
    fig,ax_p = plt.subplots(2,1)
    dict_ = {}
    for line in lines:
        dict_[line.split('\n')[0].split(' ')[0]] = line.split('\n')[0].split(' ')[1:]
    datalist = []
    datalist_name = []
    for key in dict_.keys():
        len_dict = len(dict_[key])
        for i in range(len_dict):
            os.system("sumpdos.sh {} {}".format(key,dict_[key][i]))
            data = np.loadtxt("{}-{}.dat".format(key,dict_[key][i]))
            datalist.append(data)
            datalist_name.append('{} {}'.format(key,dict_[key][i]))
    en_l = datalist[0][:,0]
    for i,_ in enumerate(datalist):
        ax_p[0].plot(en_l-fermi,datalist[i][:,1],color=color[i],label=datalist_name[i],lw=3.0)
    ax_p[1].plot(dos[:,0]-fermi,dos[:,1],'k-',lw=3.0)
    ind = np.logical_and(dos[:,0]-fermi > -8,dos[:,0]-fermi < 4.1)
    dos_range = dos[ind][:,1]
    maxdos = dos_range.max()
    ax_p[0].plot([0,0], [0,maxdos], 'k-.', lw=0.75)
    ax_p[1].plot([0,0], [0,maxdos], 'k-.', lw=0.75)
    if len(datalist) < 5:
        ax_p[0].legend(loc="best",frameon=False,fontsize=20)
    else:
        ax_p[0].legend(ncol=2,loc="best",frameon=False,fontsize=20)
    ax_p[1].set_xlabel(r"E - E$_F$ (eV)", fontsize=30)
    ax_p[0].set_ylabel("PDOS (states/eV/cell)",fontsize=20)
    ax_p[1].set_ylabel("DOS (states/eV/cell)",fontsize=20)
    ax_p[0].set_xticklabels([])
    for tick in ax_p[0].yaxis.get_major_ticks():
        tick.label.set_fontsize(30)
    for tick in ax_p[1].yaxis.get_major_ticks():
        tick.label.set_fontsize(30)
    plt.xticks([-8,-6,-4,-2,0,2,4], ["-8", "-6", "-4", "-2", "0", "2", "4"], fontsize = 30)
    #plt.yticks([0, 10, 20, 30, 40], ["0", "10", "20", "30", "40"], fontsize=30)
    fig.set_figheight(8)
    fig.set_figwidth(18)
    plt.subplots_adjust(bottom=0.15)
    ylim = input_data['plot']['ylim']
    xlim = input_data['plot']['xlim']
    if ylim is not None:
        ax_p[0].set_ylim((ylim[0], ylim[1])) #change this depending on the energy levels.
        ax_p[1].set_ylim((ylim[0], ylim[1])) #change this depending on the energy levels.
    else:
        ax_p[0].set_ylim(0,maxdos)
        ax_p[1].set_ylim(0,maxdos)
    if xlim is not None:
        ax_p[0].set_xlim(xlim[0],xlim[1])
        ax_p[1].set_xlim(xlim[0],xlim[1])
    else:
        ax_p[0].set_xlim(-8,4)
        ax_p[1].set_xlim(-8,4)
    plt.savefig(out)
def band_wann_plot(fileband='ex_band.dat',fileout='plot.pdf'):
    """
    Function to plot Wannier interpolated bandstructure.

    Parameters:
    -----------
    fileband : str, optional
        Bandstructure file. Default is 'ex_band.dat'.
    fileout : str, optional
        Output file. Default is 'plot.pdf'.

    Returns:
    --------
    None

    This function reads the bandstructure file and plots the Wannier interpolated bandstructure.
    It extracts data from the provided files and plots the bandstructure accordingly.
    The function also checks for bands below the Fermi level and prints a message if found.

    """
    print("plotting wannier band\n")
    os.system("cat ex_band.labelinfo.dat | awk '{ print $2 }' | tail -n 1 > n.dat")
    os.system("cat *_band.labelinfo.dat | awk '{ print $1 }' > label")
    os.system("cat *_band.labelinfo.dat | awk '{ print $3 }' > lbpoint")
    sympt = np.loadtxt('lbpoint')
    nkpt = int(np.loadtxt("n.dat"))
    os.system("grep Fermi scf.out | awk '{ print $5 }' > fermi.dat")
    with open("fermi.dat", "r") as read_fermi:
        fermi = float(read_fermi.readlines()[0].split("\n")[0])
    with open("label", "r") as label:
        lines = label.readlines()
    symlb = []
    for line in lines:
        symlb.append(line.split()[0])
    color = ['r', 'b', 'r','g']
    linestyle = ['solid', 'dashed', 'dotted','-.']
    _,ax_p = plt.subplots()
    data=np.loadtxt(fileband)
    pt_k = data[:nkpt,0]
    nband=int(data.shape[0]/nkpt)
    print("nkpoints:{}-nband:{}".format(nkpt,nband))
    for i in range(nband):
        ax_p.plot(pt_k, data[0+nkpt*i:nkpt*(i+1)][:,1], lw = 0.5, linestyle=linestyle[0],color=color[0])
        if np.any(data[0+nkpt*i:nkpt*(i+1)][:,1] < fermi):
            print("Band below Fermi level: {} \n".format(i+1))
    ax_p.set_xticks(sympt)
    ax_p.set_xticklabels(symlb)
    #ax.set_ylim(fermi-1,fermi+1)
    #for i in range(sympt.shape[0]):
    #    ax.vlines(sympt[i],11.0,18.0,color='k',linestyles='dashed')
    ax_p.set_ylabel("Energy (eV)")
    ax_p.plot([sympt[0],sympt[-1]],[fermi,fermi],'k--', lw=0.5)
    pylab.savefig(fileout, format='pdf',bbox_inches='tight')

def plot_projection(scf_file,projection_file,phonon_freq,outfile,nkpt):
    """
    Function to plot atomic projection on phonon bandstructure.

    Parameters:
    -----------
    scf_file : str
        QE scf.in file to extract structure.
    projection_file : str
        File containing atomic projection data.
    phonon_freq : str
        File containing phonon dispersion data.
    outfile : str
        Output plot file.
    nkpt : int
        Number of q points.
    proj_cutoff : float, optional
        Cutoff to apply filter. Plotting only colors that have projection larger than the cutoff.
        If not provided, it will be calculated based on the number of atoms.

    Returns:
    --------
    None

    This function reads the necessary files and plots the atomic projection on the phonon bandstructure.
    It calculates the projection cutoff based on the number of atoms if not provided explicitly.
    The function generates a plot showing the atomic projection on the phonon dispersion.
    """
    _,sympt,symlb,pt_l,_,_ = kpath(scf_file,nkpt,0)
    pt_l = np.array(pt_l)
    sympt = np.array(sympt)
    #print(pt_l,sympt,symlb)
    data = np.loadtxt(phonon_freq)
    nband = data.shape[1]
    proj = np.loadtxt(projection_file)
    proj_list = []
    nat = int(proj.shape[0]/nkpt)
    if nat < 4:
        proj_cutoff = 0.6
    elif nat == 4:
        proj_cutoff = 0.5
    else:
        print("Compound has more than 4 ions. So not plotting\n")
        sys.exit()
    print("Cutoff filter for projection: {}\n".format(proj_cutoff))
    print("Change proj_cutoff variable inside plot.py, if you need otherwise\n")
    for i in range(nat):
        proj_list.append(proj[i*nkpt:(i+1)*nkpt,:])
    percmto_thz = 33.356
    proj_array = np.zeros_like(proj_list[0])
    color_list = ['red','blue','green','cyan','k']
    for i,proj_i in enumerate(proj_list):
        proj_array[np.where(proj_i > proj_cutoff)] = i
    _,ax_p = plt.subplots()
    for i in range(1,nband):
        xdata = pt_l
        ydata = data[:,i]/percmto_thz
        cdata = proj_array[:,i-1]
        for j in np.arange(len(xdata) - 1):
            ax_p.plot([xdata[j], xdata[j+1]], [ydata[j], ydata[j+1]], lw=1.5, color=color_list[int(cdata[j])])
    maxy = np.max(data)/percmto_thz + 2.0
    miny = np.min(data)/percmto_thz - 1.0
    plt.ylabel(r'$\omega$ (THz)',fontsize=20)
    for i in range(sympt.shape[0]):
        ax_p.vlines(sympt[i],0,maxy,color='k',linestyles='dashed')
    plt.xticks(sympt,symlb)
    plt.ylim(miny,maxy)
    plt.xlim(sympt[0],sympt[-1])
    plt.savefig(outfile)
def write_filedos(comp):
    """
    Function to write a 'filedos.in' file required for PDOS plot.

    Parameters:
    -----------
    comp : str
        Composition.

    Returns:
    --------
    None

    This function writes the 'filedos.in' file based on the composition provided.
    It extracts the elements and their electronic structure to determine the orbital contributions for PDOS.
    """
    if "-" in comp:
        comp = comp.split("-")[0]
    comp = Composition(comp)
    nelm = comp.elements
    with open('filedos.in', 'w') as write_pdos:
        for elm in nelm:
            els = elm.full_electronic_structure
            orb_list = []
            write_pdos.write(elm.symbol + " ")
            for orb in els:
                if 's' in orb and 's' not in orb_list:
                    orb_list.append('s')
                elif 'p' in orb and 'p' not in orb_list:
                    orb_list.append('p')
                elif 'd' in orb and 'd' not in orb_list:
                    orb_list.append('d')
                elif 'f' in orb and 'f' not in orb_list:
                    orb_list.append('f')
                else:
                    continue
            for i,orb in enumerate(orb_list):
                if i < len(orb_list) - 1:
                    write_pdos.write(orb + " ")
                else:
                    write_pdos.write(orb)
            write_pdos.write("\n")
def dos_plot_vasp(outfile="pdos.pdf"):
    """
    Function to plot DOS and partial DOS (pDOS) using the vasprun.xml file from VASP using the pymatgen package.

    Parameters:
    -----------
    outfile : str, optional
        Output plot file name (default is "pdos.pdf").

    Returns:
    --------
    None

    This function reads the 'filedos.in' file to determine which elements and orbitals to include in the pDOS plot.
    It uses the Vasprun object from the pymatgen package to parse the vasprun.xml file.
    The function plots the Total DOS and pDOS for the specified elements and orbitals.

    Note:
    -----
    'filedos.in' should be present in the current directory to specify the elements and orbitals for pDOS plotting.
    """
    try:
        with open("filedos.in", "r") as p_dos:
            lines = p_dos.readlines()
    except FileNotFoundError:
        print("filedos.in not found\n")
        sys.exit()
    try:
        result = Vasprun('vasprun.xml',parse_dos=True)
    except:
        result = Vasprun('vasprun.xml', parse_potcar_file=False)
    complete_dos = result.complete_dos
    nspin = len(complete_dos.densities.keys())
    plotter = DosPlotter()
    xlim1 = input_data['plot']['xlim']
    ylim1 = input_data['plot']['ylim']
    xmin = xlim1[0]
    xmax = xlim1[1]
    ymin = ylim1[0]
    ymax = ylim1[1]
    if nspin == 1:
        plotter.add_dos('Total DOS', result.tdos)
        dict_ = {}
        for line in lines:
            dict_[line.split('\n')[0].split(' ')[0]] = line.split('\n')[0].split(' ')[1:]
        for key in dict_.keys():
            value = dict_[key]
            len_value = len(value)
            pdos_ion = complete_dos.get_element_spd_dos(key)
            for i in range(len_value):
                plotter.add_dos("{}({})".format(key,value[i]),pdos_ion[OrbitalType[value[i]]])
        if xlim1 is not None:
            plot_axis = plotter.get_plot(xlim=(xmin,xmax))
        else:
            plot_axis = plotter.get_plot(xlim=(-4,4))
        plot_axis.legend(loc="upper right")
        plot_axis.set_ylabel("Density of States (states/eV)",fontweight='bold',fontsize=25)
        plot_axis.set_xlabel(r"Energy - E$_F$ (eV)",fontweight='bold',fontsize=25)
    elif nspin == 2:
        dos = complete_dos
        energies = dos.energies - dos.efermi
        spin = list(dos.densities.keys())
        dos_spin_up = dos.densities[spin[0]]
        dos_spin_down = -1*dos.densities[spin[1]]
        plt.plot(energies, dos_spin_up, label="Up", color='b')
        plt.plot(energies, dos_spin_down, label="Down", color='r')
        plt.axhline(0, color='black', linestyle='--', linewidth=0.5)
        plt.axvline(0, color='black', linestyle='--', linewidth=0.5)
        plt.legend(bbox_to_anchor=(1.05, 1),loc="upper right",frameon=False)
        plt.ylabel("DOS",fontweight='bold',fontsize=15)
        plt.xlabel(r"E - E$_F$ (eV)",fontweight='bold',fontsize=15)
        if xlim1 is not None:
            plt.xlim(xmin,xmax)
        else:
            plt.xlim(-4,4)
        if ylim1 is not None:
            plt.ylim(ymin,ymax)
        else:
            plt.xlim(-10,40)
        
    else:
        print("nspin should be 1 or 2\n")
    plt.savefig(outfile, dpi=500)
def band_plot_vasp_line(mpid,compound):
    """
    Function to plot the bandstructure in line mode.

    Parameters:
    -----------
    mpid : str
        Materials ID.
    compound : str
        Name of the compound.

    Returns:
    --------
    None

    This function generates a PDF plot of the bandstructure in line mode using data from the vasprun.xml file.
    The plot is saved with the filename format: "<mpid>-<compound>-band.pdf".
    The band structure is obtained using the Vasprun object from the pymatgen package.
    The BSPlotter object is used to generate the band structure plot with specified settings.

    Note:
    -----
    Ensure that the vasprun.xml file containing band structure data is present in the current directory.
    """
    outfile = mpid + "-" + compound + "-" + "band.pdf"
    vaspout = Vasprun("vasprun.xml")
    bandstr = vaspout.get_band_structure(line_mode=True)
    plt1 = BSPlotter(bandstr)
    ylim1 = input_data['plot']['ylim']
    if ylim1 is not None:
        plt2 = plt1.get_plot(ylim=[ylim1[0],ylim1[1]],vbm_cbm_marker=True,zero_to_efermi=True)
    else:
        plt2 = plt1.get_plot(ylim=[-1.8,1.8],vbm_cbm_marker=True,zero_to_efermi=True)
    plt2.legend('',frameon=False)
    plt2.figure.savefig(outfile)
def main():
    """
    Main function to execute plotting tasks based on command-line arguments.

    This function reads command-line arguments to determine the type of plot and other required parameters.
    It performs different plotting tasks based on the specified plot type.

    Command-line Arguments:
    -----------------------
    plottype : str
        Type of plot to generate.
    mpid : str
        Materials ID.
    comp : str
        Name of the compound.

    Returns:
    --------
    None

    Plotting Tasks:
    ---------------
    - For 'pdos' plot type:
        - Checks if 'filedos.in' exists, creates one if not found.
        - Determines the type of calculation (QE or VASP).
        - Calls 'dos_plot_vasp' function if VASP calculation is detected, otherwise 'dos_plot' function.

    - For 'wann_band' plot type:
        - Calls 'band_wann_plot' function.

    - For 'phonproj' plot type:
        - Reads the number of k-points from command-line argument.
        - Calls 'plot_projection' function with appropriate parameters.

    - For other plot types:
        - Calls 'plot' function with specified plot type, file name, and compound name.

    Note:
    -----
    Ensure proper command-line arguments are provided.
    Only QE and VASP outputs are allowed for plotting tasks.
    """
    plottype = sys.argv[1]
    mpid = sys.argv[2]
    comp = sys.argv[3]
    if os.path.isfile('scf.in'):
        filename = 'scf.in'
    elif os.path.isfile('POSCAR'):
        filename = 'POSCAR'
    else:
        print("No scf.in and POSCAR exists\n")
    #input_file = espresso.read_espresso_in('scf.in')
    if plottype == 'pdos':
        if not os.path.isfile("filedos.in"):
            print("filedos.in file not found, creating one\n")
            write_filedos(comp)
        if os.path.isfile('POSCAR'):
            print("POSCAR present. searching vasprun.xml..\n")
            dos_plot_vasp("pdos-{}.pdf".format(comp))
        elif os.path.isfile('scf.in'):
            print("scf.in present. searching {}.dos and *pdos.pdos* files..\n".format(comp))
            dos_plot("{}.dos".format(comp),"pdos-{}.pdf".format(comp))
        else:
            print("Only QE and VASP outputs are allowed\n")
    elif plottype == 'wann_band':
        band_wann_plot()
    elif plottype == 'phonproj':
        nkpt = int(sys.argv[4])
        plot_projection(filename,"phonon-{}.proj.gp".format(comp),"{}.freq.gp".format(comp),"plot-proj-{}-{}.pdf".format(mpid,comp),nkpt)
    #elif plottype == 'line':
    #    band_plot_vasp_line(mpid,comp)
    else:
        plot(plottype,filename,comp)
if __name__ == "__main__":
    main()
