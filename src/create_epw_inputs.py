#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D."""
import sys
import os
import glob
import json
import math as m
import numpy as np
from lmfit import Model
from scipy.special import erfc
import matplotlib.pyplot as plt
from kpoint_path import kpoint_path
from wannier90 import epw_bandcheck

def phonon_input(mass_file, qpoint_file, mpid, compound, prefix):
    """
    function to write ph.in file for phonon calculation
    parameters
    --------------
    mass_file : mass.in file automatically created in 'create-inputs' bash script
    qpoint_file : qpoint.in file provided
    mpid : MP id
    compound : MP compound name
    prefix : prefix used in QE calculations
    Returns: None
    --------------
    write a file in 'ph-mpid-compound.in' format
    """
    mass=np.loadtxt(mass_file)
    if not os.path.isfile("ph-q.in"):
        qpt=np.loadtxt(qpoint_file)
        print("q-mesh: {}".format(qpt))
    if mass.ndim > 0:
        nat = mass.shape[0]
    else:
        nat = 1
        mass = [mass]
    dynmat = prefix.replace("'", "") + ".dyn"
    with open("ph-{}-{}.in".format(mpid,compound), 'w') as ph_write:
        ph_write.write("phonon calculation \n")
        ph_write.write("&inputph" + "\n")
        ph_write.write("tr2_ph=1.0d-16," + "\n")
        ph_write.write("prefix={},".format(prefix) + "\n")
        ph_write.write("fildvscf='dvscf'," + "\n")
        for i in range(1,nat+1):
            ph_write.write("amass({})={},".format(i,mass[i-1]) + "\n")
        ph_write.write("outdir='./'," + "\n")
        ph_write.write("fildyn='{}',".format(dynmat) + "\n")
        if not os.path.isfile("ph-q.in"):
            print("No ph-q.in file provided\n")
            print("calculations will be performed for metal and qpoint obtained from qpoint.in\n")
            ph_write.write("ldisp=.true." + "\n")
            ph_write.write("nq1={},nq2={},nq3={}".format(int(qpt[0]),int(qpt[1]),int(qpt[2])) + "\n")
            ph_write.write("/" + "\n")
        else:
            print("Reading ph-q.in file\n")
            with open("ph-q.in", "r") as ph_read:
                lines = ph_read.readlines()
            nqpoint = lines[0].split("\n")[0].split("#")[0].split()
            is_metal = lines[1].split()[0]
            print("nq = " + nqpoint[0] + " " + nqpoint[1] + " ")
            print(nqpoint[2] + "," + "metal = " + is_metal + "\n")
            if is_metal in ('T','t'):
                ph_write.write("reduce_io = .true.\n")
                ph_write.write("/" + "\n")
                ph_write.write(nqpoint[0] + " " + nqpoint[1] + " " + nqpoint[2] + "\n")
            else:
                ph_write.write("reduce_io = .true.\n")
                ph_write.write("epsil = .true.\n")
                ph_write.write("/" + "\n")
                ph_write.write(nqpoint[0] + " " + nqpoint[1] + " " + nqpoint[2] + "\n")
def prepare_nscf(grid_file,what='nscf'):
    """
    function to write file for non-scf calculations
    parameters
    ---------------
    grid_file: kmesh.grid file automatically created from bash scripts
    what : type of nscf grid. 'wan' is nscf grid for wannierization
    Returns : None
    ---------------
    write a file with kmesh grid for non-scf calculation
    """
    print("******************************************************")
    print("******************************************************\n")
    print("_______________________________________________________")
    print("_____________________________________________________\n")
    grid = np.loadtxt(grid_file)
    if what == 'wan':
        os.system("kmesh_nscf.py {} {} {} wan".format(int(grid[0]), int(grid[1]), int(grid[2])) + " " + ">" + " " + "wann_grid.out")
    else:
        os.system("kmesh_nscf.py {} {} {}".format(int(grid[0]), int(grid[1]), int(grid[2])) + " " + ">" + " " + "nscf_grid.out")

def projwfc(prefix,mpid,compound):
    """
    function to write projwfc.in file for orbital projection calculation
    parameters
    --------------
    mpid : MP id
    compound : MP compound name
    prefix : prefix used in QE calculations
    Returns : None
    --------------
    write a file in 'projwfc-mpid-compound.in' format
    """
    dynmat = prefix.replace("'", "") + ".proj"
    with open("projwfc-{}-{}.in".format(mpid,compound), 'w') as projwfc_write:
        projwfc_write.write("&projwfc\n")
        projwfc_write.write("prefix={},".format(prefix) + "\n")
        projwfc_write.write("outdir='./'," + "\n")
        projwfc_write.write("lsym = .False.,\n")
        projwfc_write.write("filproj = '{}',".format(dynmat) + "\n")
        projwfc_write.write("/" + "\n")

def ciftoxsf(file_name,outfile):
    """
    function to change .cif file to .xsf format, suitable for xcrysden.
    parameters
    -------------
    file_name : .cif filename
    outfile : .xsf filename
    Returns : None
    """
    data=read(file_name,format='cif')
    data.write(outfile,format='xsf')

def func_erfc(xdata,aparam,bparam,func_type='erfc'):
    """
    erfc function erfc((x-a)/b)
    parameters
    --------------
    xdata : x data
    aparam and bparam are parameters
    func_type : Type of function
    Returns
    --------------
    ydata = 0.5 * erfc((xdata-aparam)/bparam)
    """
    if func_type == 'erfc':
        ydata = 0.5 * erfc((xdata - aparam)/bparam)
    elif func_type == 'gauss':
        ydata = m.exp(-(xdata - aparam)**2.0/bparam**2.0)
    else:
        print('Either erfc or gauss option available for entanglement\n')
    return ydata

def scdmfit(func,file_name,mpid,compound,fermi,entang='erfc'):
    """
    function to fit erfc function and calculate SCDM parameters
    to generate automatic projections for  wannierization
    parameters
    --------------
    func : function to fit
    file_name : Data in 'energy, projection' format
    mpid : MP id
    compound : compound
    fermi : Fermi level
    entang : type of entanglement. Default: 'erfc'
    Returns
    --------------
    mu_scdm, sigma_scdm parameters
    """
    gmodel = Model(func)
    data = np.loadtxt(file_name)
    if entang == 'erfc':
        idx = (np.abs(data[:,1]-0.5)).argmin()
    elif entang == 'gauss':
        idx = np.where(data[:,1] == np.max(data[:,1]))
    else:
        print('Either erfc or gauss option available for entanglement\n')
    params = gmodel.make_params(aparam=data[:,0][idx], bparam=1.0)
    result = gmodel.fit(data[:,1], params, xdata=data[:,0])
    afit = result.best_values['aparam']
    bfit = result.best_values['bparam']
    mu_scdm = afit - 3*bfit
    sigma_scdm = bfit
    with open("p_vs_e_fit-{}-{}.csv".format(mpid,compound), "w") as pe_fit:
        for i in range(data.shape[0]):
            pe_fit.write("{},{},{}".format(data[:,0][i],data[:,1][i],result.best_fit[i]) + "\n")
    plt.plot(data[:,0], data[:,1], 'bo', label='projection')
    plt.plot(data[:,0], result.best_fit, color="red", linestyle="solid", linewidth=2.0, label='fit')
    plt.plot(np.array([fermi,fermi]), np.array([0,1]), color="blue", linestyle="dashed", linewidth=2.0, label="E-fermi")
    plt.plot(np.array([mu_scdm,mu_scdm]), np.array([0,1]), color="black", linestyle="solid", linewidth=2.0, label=r"$\mu_c$")
    plt.xlabel("Energy (eV)", fontsize=20)
    plt.ylabel("Projection", fontsize=20)
    plt.legend(loc="lower left", fontsize=15)
    plt.title(r"$\mu_c$ = {}, $\sigma$ = {}".format(round(afit,4),round(bfit,4)))
    plt.savefig("scdm-proj-{}-{}.png".format(mpid,compound))
    return mu_scdm,sigma_scdm

#Now, let's test wannierization by producing wannier interpolated band structures.
# After that, we can perform epw calculations for superconductivity.
# Also, lambda_qnu projection or atomic projection to the phonon band.
def band_file_epw(mpid,compound):
    """
    function to include 'nosym = .true.' inside QE scf.in file
    for bandstructure calculation
    parameters
    -------------
    mpid : Materials id
    compound : compound name
    """
    if os.path.isfile("scf_dir/scf-{}-{}-band.in".format(mpid,compound)):
        os.system("sed '/&SYSTEM/a   nosym = .true.,'" + " scf_dir/scf-{}-{}-band.in".format(mpid,compound)  + " " + ">" + "scf_dir/{}-{}-band.in".format(mpid,compound))
    else:
        print("scf-{}-{}-band.in not found inside scf_dir.\n")
        print("Run the script 'mainprogram.py process', with process from 1 to 4".format(mpid,compound) + "\n")
        print("Check 'mainprogram.py qe-info' to find what processes from 1 to 4 do.\n")

def pw2wan_input(mpid,compound,prefix,proj=" "):
    """
    function to write pw2wan input file
    parameters
    -----------
    mpid : Materials id
    compound : compound name
    prefix : prefix used in QE calculations
    proj : 'scdm' if used SCDM projection, '' otherwise
    Returns
    -------------
    write 'pw2wan-mpid-compound.in' file
    """
    if os.path.isfile("scdm_dir/scdm-{}-{}".format(mpid,compound)):
        data=np.loadtxt("scdm_dir/scdm-{}-{}".format(mpid,compound))
        mu_scdm = data[0]
        sigma_scdm = data[1]
    if not os.path.isdir('epw_dir'):
        os.system("mkdir epw_dir")
    with open("epw_dir/pw2wan-{}-{}.in".format(mpid,compound), "w") as pw2wan_write:
        pw2wan_write.write("&inputpp" + "\n")
        pw2wan_write.write("prefix={}".format(prefix) + "\n")
        pw2wan_write.write("outdir='./'" + "\n")
        pw2wan_write.write("seedname = 'ex'" + "\n")
        pw2wan_write.write("write_amn = .true.\n")
        pw2wan_write.write("write_mmn = .true.\n")
        pw2wan_write.write("write_unk = .true.\n")
        if proj == "scdm":
            pw2wan_write.write("scdm_entanglement = 'erfc'\n")
            pw2wan_write.write("scdm_proj = .true.\n")
            pw2wan_write.write("scdm_mu = {}".format(mu_scdm) + "\n")
            pw2wan_write.write("scdm_sigma = {}".format(sigma_scdm) + "\n")
        pw2wan_write.write("/" + "\n")

#def epw_bandcheck(infile='scf.in',out="ex.win",proj=" ",num_iter=300,dis_num_iter=200):
#    """
#    function to write input file for wannier90 calculations. Default: 'ex.win'
#
#    parameters
#    -------------
#    infile : input file for scf calculation
#    out : output file. Default: 'ex.win'
#    proj : 'scdm' if used SCDM projection, '' otherwise
#    num_iter : number of iteration for wannierization
#    dis_num_iter : number of iteration for disentanglement process
#
#
#    """
#    if os.path.isfile("scf_dir/{}".format(infile)):
#        data = espresso.read_espresso_in("scf_dir/{}".format(infile))
#    else:
#        print("No {} inside scf_dir".format(infile) + "\n")
#    #l=data.cell.get_bravais_lattice()
#    cell=data.cell
#    symbol=data.get_chemical_symbols()
#    pos=data.get_scaled_positions()
#    kmesh = np.loadtxt('kmesh.grid')
#
#    with open(out, 'w') as epw_write:
#        epw_write.write("use_ws_distance = .true.\n")
#        epw_write.write("dis_num_iter = {}".format(dis_num_iter) + "\n")
#        epw_write.write("write_hr = .true.\n")
#        epw_write.write("iprint = 2\n")
#        epw_write.write("spinors = .false.\n")
#        if proj == 'scdm':
#            epw_write.write("auto_projections = .true.\n")
#        elif proj == "fromfile":
#            if os.path.isfile("projection.in"):
#                with open("projection.in", "r") as gfile:
#                    lines = gfile.readlines()
#                len_l = len(lines)
#                epw_write.write("begin projections\n")
#                for i in range(len_l):
#                    #f.write("proj({})".format(i+1) + "=" + "'{}'".format(lines[i].split("\n")[0]) + "\n")
#                    epw_write.write("{}".format(lines[i].split("\n")[0]) + "\n")
#                epw_write.write("end projections\n")
#            else:
#                print("projection.in file not found\n")
#                print("write projections in different line, 'X:s', 'Y:pz', ... so on\n")
#        else:
#            epw_write.write("begin projections\n")
#            epw_write.write("random\n")
#            epw_write.write("end projections\n")
#        epw_write.write("\n")
#        epw_write.write("num_bands = XXX\n")
#        epw_write.write("num_wann = XXX\n")
#        epw_write.write("num_iter = {}".format(num_iter) + "\n")
#        epw_write.write("dis_froz_min = XXX ! obtain from electronic bandstructure\n")
#        epw_write.write("dis_froz_max = XXX\n")
#        epw_write.write("dis_win_min = XXX\n")
#        epw_write.write("dis_win_max = XXX\n")
#        epw_write.write("\n")
#        epw_write.write("Begin Kpoint_Path\n")
#        with open("wannier_kpath.in", "r") as gfile:
#            lines = gfile.readlines()
#        for i,line in enumerate(lines):
#            epw_write.write(line)
#        epw_write.write("End Kpoint_Path\n")
#        epw_write.write("\n")
#        epw_write.write("fermi_surface_plot = .true.\n")
#        epw_write.write("bands_plot = .true.\n")
#        epw_write.write("wannier_plot = .true.\n")
#        epw_write.write("wannier_plot_supercell = 3\n")
#        epw_write.write("\n")
#        epw_write.write("begin unit_cell_cart\n")
#        epw_write.write("Ang\n")
#        for lat in cell:
#            epw_write.write(str(lat[0]) + " " + str(lat[1]) + " " + str(lat[2]) + "\n")
#        epw_write.write("end unit_cell_cart\n")
#        epw_write.write("\n")
#        epw_write.write("begin atoms_frac\n")
#        for sym,pos in zip(symbol,pos):
#            epw_write.write(sym + " " + str(pos[0]) + " " + str(pos[1]) + " " + str(pos[2]) + "\n")
#        epw_write.write("end atoms_frac\n")
#        epw_write.write("\n")
#        epw_write.write("mp_grid : {} {} {}".format(int(kmesh[0]),int(kmesh[1]),int(kmesh[2])) + "\n")
#        epw_write.write("begin kpoints\n")
#        with open('wann_grid.out', 'r') as gfile:
#            lines = gfile.readlines()
#        for line in lines:
#            epw_write.write(line)
#        epw_write.write("end kpoints\n")


def epw_sc_from_json(json_file, out="epw.in"):
    """
    Write input file for EPW calculations based on JSON configuration.

    Parameters:
    ----------------
    json_file : str
        Path to the JSON configuration file.
    out : str, optional
        EPW input file name. Default is 'epw.in'.
    """
    with open(json_file, 'r') as f:
        config = json.load(f)

    with open(out, "a") as epw_write:
        #epw_write.write("--\n")
        #epw_write.write("&inputepw\n")
        for key, value in config.items():
            if isinstance(value, bool):
                epw_write.write("{} = .{}.\n".format(key, str(value).lower()))
            elif isinstance(value, int) or isinstance(value, float):
                epw_write.write("{} = {}\n".format(key, value))
            elif isinstance(value, str):
                epw_write.write("{} = '{}'\n".format(key, value))
            elif key == 'wdata':
                for i,wdt in enumerate(config['wdata']):
                    epw_write.write("wdata({}) = '{}'".format(i+1,value[i]) + "\n")
            elif isinstance(value, list) and key != 'wdata':
                epw_write.write("{} = {}\n".format(key, ' '.join(map(str, value))))


def epw_sc(mpid,compound,prefix,qpoint_file,kpoint_file,proj=" ",out="epw.in"):
    """
    function to write input file for EPW calculations. Default: 'epw.in'

    parameters
    ----------------
    mpid : Materials id
    compound : compound name
    prefix : prefix used in QE calculations
    qpoint_file : 'qpoint.in' file for q-mesh
    kpoint_file : 'kpoint.in' file for k-mesh
    proj: projection types
    out : epw input file

    """
    qpt=np.loadtxt(qpoint_file)
    kpt=np.loadtxt(kpoint_file)
    with open(out, "w") as epw_write:
        epw_write.write("--\n")
        epw_write.write("&inputepw\n")
        epw_write.write("prefix={}".format(prefix) + "\n")
        epw_write.write("outdir = './'" + "\n")
        #epw_write.write("dvscf_dir = '../phonon/save'" + "\n")
        #epw_write.write("ep_coupling = .true.\n")
        #epw_write.write("elph = .true.\n")
        #epw_write.write("epwwrite = .true.\n")
        #epw_write.write("epwread = .false.\n")
        #epw_write.write("etf_mem = 1\n")
        #epw_write.write("\n")
        #epw_write.write("wannierize = .true.\n")
        #epw_write.write("nbndsub = {} ! change according to your need".format(nbndsub) + "\n")
        #epw_write.write("bands_skipped = 'exclude_bands = 1:{}' ! If possible exclude bands based on bandstructure".format(band_skip) + "\n")
        #epw_write.write("wdata(1) = 'fermi_surface_plot = .true.'\n")
        #epw_write.write("wdata(2) = 'dis_num_iter = {}'".format(dis_num_iter) + "\n")
        #epw_write.write("max_memlt = XXX !determine the memory requirement\n")
        if proj == "scdm":
            print("perform epw1, epw2, epw3, epw4, and proj calculations\n")
            data=np.loadtxt("scdm_dir/scdm-{}-{}".format(mpid,compound))
            mu_scdm = round(float(data[0]),5)
            sigma_scdm = round(float(data[1]),5)
            epw_write.write("auto_projections = .true.\n")
            epw_write.write("scdm_entanglement = 'erfc'\n")
            epw_write.write("scdm_proj = .true.\n")
            epw_write.write("scdm_mu = {}".format(mu_scdm) + "\n")
            epw_write.write("scdm_sigma = {}".format(sigma_scdm) + "\n")
        if proj == "fromfile":
            if os.path.isfile("projection.in"):
                with open("projection.in", "r") as gfile:
                    lines = gfile.readlines()
                len_l = len(lines)
                for i in range(len_l):
                    epw_write.write("proj({})".format(i+1) + "=" + "'{}'".format(lines[i].split("\n")[0]) + "\n")
            else:
                print("projection.in file not found\n")
                print("write projections in different line, 'X:s', 'Y:pz', ... so on\n")

    epw_write.write("\n")
    epw_write.write("#From epw.json file\n")
    epw_sc_from_json("epw.json",out='epw.in')
        #epw_write.write("num_iter = {}".format(num_iter) + "\n")
        #epw_write.write("dis_froz_min = XXX ! obtain from electronic bandstructure\n")
        #epw_write.write("dis_froz_max = XXX\n")
        #epw_write.write("dis_win_min = XXX ! not necessary after defining exclude_bands\n")
        #epw_write.write("dis_win_max = XXX\n")
        #epw_write.write("\n")
        #epw_write.write("iverbosity = 2\n")
        #epw_write.write("fsthick = 0.2\n")
        #epw_write.write("degaussw = 0.05\n")
        #epw_write.write("ephwrite = .true.\n")
        #epw_write.write("eliashberg = .true.\n")
        #epw_write.write("\n")
        #epw_write.write("{} = .true.".format(elphtype) + "\n")
        #epw_write.write("limag = .true.\n")
        #epw_write.write("lpade = .true.\n")
        #epw_write.write("lacon = .false.\n")
        #epw_write.write("\n")
        #epw_write.write("nsiter = 300 \n")
        #epw_write.write("conv_thr_iaxis = 1.0d-4\n")
        #epw_write.write("wscut = 1.0 ! 10 times of maximum phonon frequency.\n")
        #epw_write.write("muc = 0.16\n")
        #epw_write.write("nstemp = 10\n")
        #epw_write.write("temps = 10 55\n")
        #epw_write.write("\n")
        #epw_write.write("! Other properties\n")
        #epw_write.write("!#Calculate the electron spectral function from the e-ph interaction\n")
        #epw_write.write("specfun_el = .false.\n")
        #epw_write.write("scattering = .false. !computes scattering rates\n")
        #epw_write.write("!scattering rates are calculated using 0th order relaxation time approximation.\n")
        #epw_write.write("scattering_0rta = .false.\n")
        #epw_write.write("!scattering rates in the self-energy relaxation time approximation.\n")
        #epw_write.write("scattering_serta = .false.\n")
        #epw_write.write("!el-ph matrix elements are screened by the RPA or TF dielectric function\n")
        #epw_write.write("lscreen = .false.\n")
        #epw_write.write("!Lindhard screening, if 1 the Thomas-Fermi screening. Only relevant if lscreen = .true.\n")
        #epw_write.write("scr_typ = 0 \n")
        #epw_write.write("!phonon self-energy from the el-ph interaction.\n")
        #epw_write.write("phonselfen = .false.\n")
        #epw_write.write("!Eliashberg spectral function, 𝛼2𝐹(𝜔), transport Eliashberg spectral function 𝛼2𝐹tr(𝜔)\n")
        #epw_write.write("!phonon density of states 𝐹(𝜔). Only allowed in the case of phonselfen = .true.\n")
        #epw_write.write("a2f = .false.\n")
        #epw_write.write("!Calculate the electronic nesting function.\n")
        #epw_write.write("nest_fn = .false.\n")
        #epw_write.write("!enable the correct Wannier interpolation in the case of polar material.\n")
        #epw_write.write("lpolar = .false. \n")
        #epw_write.write("!only the long-range part of the electron-phonon matrix elements are calculated\n")
        #epw_write.write("!works only with lpolar=.true.\n")
        #epw_write.write("longrange = .false.\n")
        #epw_write.write("!Calculate the electron self-energy from the el-ph interaction\n")
        #epw_write.write("elecselfen = .false.\n")
    epw_write.write("\n")
    with open(out, "a") as epw_write:
        epw_write.write("nk1 = {}".format(int(kpt[0])) + "\n")
        epw_write.write("nk2 = {}".format(int(kpt[1])) + "\n")
        epw_write.write("nk3 = {}".format(int(kpt[2])) + "\n")
        epw_write.write("\n")
        epw_write.write("nq1 = {}".format(int(qpt[0])) + "\n")
        epw_write.write("nq2 = {}".format(int(qpt[1])) + "\n")
        epw_write.write("nq3 = {}".format(int(qpt[2])) + "\n")
        epw_write.write("\n")
        epw_write.write("mp_mesh_k = .true.\n")
        epw_write.write("fermi_plot = .true.\n")
        epw_write.write("nkf1 = {}".format(int(10*kpt[0])) + "\n")
        epw_write.write("nkf2 = {}".format(int(10*kpt[1])) + "\n")
        epw_write.write("nkf3 = {}".format(int(10*kpt[2])) + "\n")
        epw_write.write("\n")
        epw_write.write("nqf1 = {}".format(int(10*qpt[0])) + "\n")
        epw_write.write("nqf2 = {}".format(int(10*qpt[1])) + "\n")
        epw_write.write("nqf3 = {}".format(int(10*qpt[2])) + "\n")
        epw_write.write("/")
def main():
    """
    main function
    """
    try:
        pwd = os.getcwd()
        if os.path.isfile(pwd+"/htepc.json"):
            jsonfile = pwd+"/htepc.json"
        else:
            jsonfile = "../../htepc.json"
        with open(jsonfile, "r") as readjson:
            input_data = json.load(readjson)
    except FileNotFoundError:
        print("htepc.json file not found\n")
    dft = input_data['download']['inp']['calc']
    mpid = sys.argv[1]
    compound = sys.argv[2]
    prefix = sys.argv[3]
    if dft in ("QE","qe"):
        projwfc(prefix,mpid,compound)
    condition = sys.argv[4]
    if condition == "ph":
        phonon_input('mass.dat', 'qpoint.dat', mpid, compound, prefix)
        print("Generating phonon input file\n")
    elif condition == "nscf":
        if not os.path.isfile("kmesh.grid"):
            os.system("echo 3 3 3 > kmesh.grid")
        prepare_nscf('kmesh.grid')
        print("Generating file for non self consistent calculations\n")
    elif condition == "ciftoxsf":
        pwd = os.getcwd()
        file_path = glob.glob(pwd+"/**/cif/{}.cif".format(mpid), recursive=True)[0]
        print(file_path)
        outfile = pwd + "/scf_dir/" + str(mpid) + ".xsf"
        ciftoxsf(file_path,outfile)
    elif condition == "projection":
        print("inside proj")
        projwfc(prefix,mpid,compound)
        print("Generating input file for projwfc.x\n")
    elif condition == "scdm":
        fermi = float(sys.argv[5])
        file_name = "R{}-{}/epw/p_vs_e.dat".format(mpid,compound)
        print("Obtaining SCDM parameters from {} projectibility".format(compound) + "\n")
        print("Please cite ==> https://doi.org/10.1137/17M1129696, if you use SCDM procedure\n")
        print("Vitale, V., Pizzi, G., Marrazzo, A. et al.\n")
        print("Automated high-throughput Wannierisation. npj Comput Mater 6, 66 (2020)\n")
        aconst,bconst = scdmfit(func_erfc,file_name,mpid,compound,fermi,'erfc')
        with open("scdm-{}-{}".format(mpid,compound), "w") as file1:
            file1.write("{} {}".format(aconst,bconst) + "\n")
    elif condition == "band":
        band_file_epw(mpid,compound)
    elif condition == "pw2wan":
        projection = sys.argv[5]
        pw2wan_input(mpid,compound,prefix,proj=projection)
    elif condition == "kpathwan":
        if os.path.isfile("scf_dir/scf-{}-{}.in".format(mpid,compound)):
            file_name="scf_dir/scf-{}-{}.in".format(mpid,compound)
        else:
            print("scf-{}-{}.in not found in scf_dir/".format(mpid,compound) + "\n")
            print("Complete mainprogram.py from 0 to 6\n")
        kpoint_path(file_name)
    elif condition == "wankmesh":
        prepare_nscf('kmesh.grid',what="wan")
    elif condition == "epw_band":
        #with open("band.dat", "r") as file:
        #    lines = file.readlines()
        #num_bands = int(lines[0].split("\n")[0].split("=")[1].split(",")[0])
        projection = sys.argv[5]
        if dft in ('QE', 'qe'):
            epw_bandcheck("scf-{}-{}.in".format(mpid,compound),proj=projection,out="ex.win",dft=dft)
            os.system("mv ex.win epw_dir/ex-{}-{}.win".format(mpid,compound))
        else:
            epw_bandcheck("POSCAR",proj=projection,out="ex.win",dft=dft)
    elif condition == "epw":
        print("creating epw inputs for superconductivity\n")
        #with open("band.dat", "r") as file1:
        #    lines = file1.readlines()
        projection = sys.argv[5]
        epw_sc(mpid,compound,prefix,'qpoint.dat','kpoint.dat',proj=projection,out="epw.in")
    else:
        print("bad inputs")
if __name__ == "__main__":
    main()
