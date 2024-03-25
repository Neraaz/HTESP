#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to prepare input files for wannier90 and EPW calculations"""
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
                It has mass in a.m.u in each line as in QE scf.in file
    qpoint_file : qpoint.in file provided
    mpid : MP id
    compound : MP compound name
    prefix : prefix used in QE calculations
    Returns: None
    --------------
    write a file in 'ph-mpid-compound.in' format
    Example:
    >>> phonon_input('mass.in', 'qpoint.in', 'mp-123', 'SiO2', 'si2o4')
    """
    # Load mass data from mass_file
    mass=np.loadtxt(mass_file)
    if not os.path.isfile("ph-q.in"):
        qpt=np.loadtxt(qpoint_file)
        print("q-mesh: {}".format(qpt))
    # Determine unqiue species in the system
    if mass.ndim > 0:
        nat = mass.shape[0]
    else:
        nat = 1
        mass = [mass]
    # Set the name of the dynamical matrix file
    dynmat = prefix.replace("'", "") + ".dyn"
    # Write the input file for phonon calculation
    with open("ph-{}-{}.in".format(mpid,compound), 'w') as ph_write:
        ph_write.write("phonon calculation \n")
        ph_write.write("&inputph" + "\n")
        ph_write.write("tr2_ph=1.0d-16," + "\n")
        ph_write.write("prefix={},".format(prefix) + "\n")
        ph_write.write("fildvscf='dvscf'," + "\n")
        # Write atomic masses for each atom
        for i in range(1,nat+1):
            ph_write.write("amass({})={},".format(i,mass[i-1]) + "\n")
        ph_write.write("outdir='./'," + "\n")
        ph_write.write("fildyn='{}',".format(dynmat) + "\n")
        # Extract a generic q-point from ph-q.in file if provided
        # Otherwise, use q-point mesh from qpoint.in file
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
                # For non metal at q = 0 0 0, turn epsil flag .true.
                ph_write.write("epsil = .true.\n")
                ph_write.write("/" + "\n")
                ph_write.write(nqpoint[0] + " " + nqpoint[1] + " " + nqpoint[2] + "\n")
def prepare_nscf(grid_file,what='nscf'):
    """
    function to write file for non-scf calculations
    parameters
    ---------------
    grid_file: kmesh.grid file automatically created from bash scripts
               It just has 3 numbers, representing the k-mesh, nkx nky nkz (no offset)
    what : type of nscf grid. 'wan' is nscf grid for wannierization
    Returns : None
    ---------------
    write a file with kmesh grid for non-scf calculation

    Example:
    >>> prepare_nscf('kmesh.grid', what='wan')
    """
    print("******************************************************")
    print("******************************************************\n")
    print("_______________________________________________________")
    print("_____________________________________________________\n")
    grid = np.loadtxt(grid_file)
    if what == 'wan':
        # Generate the nscf grid for wannierization
        os.system("kmesh_nscf.py {} {} {} wan".format(int(grid[0]), int(grid[1]), int(grid[2])) + " " + ">" + " " + "wann_grid.out")
    else:
        # Generate the standard nscf grid
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

    Example:
    --------
    >>> projwfc('prefix', 'mp-123', 'FeS2')
    This will create a file named 'projwfc-mp-123-FeS2.in' with the necessary parameters.
    """
    # Define the name for the projected wavefunction file
    dynmat = prefix.replace("'", "") + ".proj"
    # Write projwfc.in file to be executed by projwfc.x
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

    Example:
    --------
    >>> scdmfit(erfc, 'data.csv', 'mp-123', 'FeS2', 5.0, entang='erfc')
    This will fit the erfc function to the data in 'data.csv',
    calculate SCDM parameters, and return mu and sigma.
    """
    # Create a model for fitting
    gmodel = Model(func)
    # Load data
    data = np.loadtxt(file_name)
    # Find the index close to mean of the distribution
    if entang == 'erfc':
        idx = (np.abs(data[:,1]-0.5)).argmin()
    elif entang == 'gauss':
        idx = np.where(data[:,1] == np.max(data[:,1]))
    else:
        print('Either erfc or gauss option available for entanglement\n')
    # Set initial parameters for fitting
    params = gmodel.make_params(aparam=data[:,0][idx], bparam=1.0)
    # Fit the model to the data
    result = gmodel.fit(data[:,1], params, xdata=data[:,0])
    # Calculate SCDM parameters
    afit = result.best_values['aparam']
    bfit = result.best_values['bparam']
    mu_scdm = afit - 3*bfit
    sigma_scdm = bfit
    # Write fitted data to a file
    with open("p_vs_e_fit-{}-{}.csv".format(mpid,compound), "w") as pe_fit:
        for i in range(data.shape[0]):
            pe_fit.write("{},{},{}".format(data[:,0][i],data[:,1][i],result.best_fit[i]) + "\n")
    # Plot fitted data
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

    Example:
    >>> band_file_epw("mp-1234", "Si")
    """
    # Check if the scf.in file for bandstructure calculation exists
    if os.path.isfile("scf_dir/scf-{}-{}-band.in".format(mpid,compound)):
        # Add 'nosym = .true.' to the scf.in file
        os.system("sed '/&SYSTEM/a   nosym = .true.,'" + " scf_dir/scf-{}-{}-band.in".format(mpid,compound)  + " " + ">" + "scf_dir/{}-{}-band.in".format(mpid,compound))
    else:
        # Print an error message if the file is not found
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

    Example:
    >>> pw2wan_input("mp-1234", "Si", "Si", proj="scdm")
    """
    # Check if the SCDM projection file exists
    if proj == "scdm":
        # Load SCDM parameters if the file exists
        if os.path.isfile("scdm_dir/scdm-{}-{}".format(mpid,compound)):
            data=np.loadtxt("scdm_dir/scdm-{}-{}".format(mpid,compound))
            mu_scdm = data[0]
            sigma_scdm = data[1]
    # Create a directory for epw_dir if it doesn't exist
    if not os.path.isdir('epw_dir'):
        os.system("mkdir epw_dir")
    # Write pw2wan input file
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

def epw_sc_from_json(json_file, out="epw.in"):
    """
    Write input file for EPW calculations based on JSON configuration.

    Parameters:
    ----------------
    json_file : str
        Path to the JSON configuration file.
    out : str, optional
        EPW input file name. Default is 'epw.in'.
    Example:
    >>> epw_sc_from_json("epw.json", out="epw.in")
    """
    # Load JSON data
    with open(json_file, 'r') as f:
        config = json.load(f)
    # Write EPW input file
    with open(out, "a") as epw_write:
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

    Example:
    >>> epw_sc("mp-123", "SiO2", "Si2O4", "qpoint.in", "kpoint.in", proj="", out="epw.in")
    """
    # Load q-point and k-point files
    qpt=np.loadtxt(qpoint_file)
    kpt=np.loadtxt(kpoint_file)
    # Write EPW input file
    with open(out, "w") as epw_write:
        epw_write.write("--\n")
        epw_write.write("&inputepw\n")
        epw_write.write("prefix={}".format(prefix) + "\n")
        epw_write.write("outdir = './'" + "\n")
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
            # Use projections from a file
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
    # Write input from epw.json file
    epw_sc_from_json("epw.json",out='epw.in')
    epw_write.write("\n")
    # Write k-point and q-point information
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
    # Load input data from JSON file
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
    # Extract relevant data from input
    dft = input_data['download']['inp']['calc']
    mpid = sys.argv[1]
    compound = sys.argv[2]
    prefix = sys.argv[3]
    # Perform tasks based on input conditions
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
        # Convert CIF file to XSF format
        pwd = os.getcwd()
        file_path = glob.glob(pwd+"/**/cif/{}.cif".format(mpid), recursive=True)[0]
        print(file_path)
        outfile = pwd + "/scf_dir/" + str(mpid) + ".xsf"
        ciftoxsf(file_path,outfile)
    elif condition == "projection":
        # Generate input file for projwfc.x
        print("inside proj")
        projwfc(prefix,mpid,compound)
        print("Generating input file for projwfc.x\n")
    elif condition == "scdm":
        # Perform SCDM parameter fitting
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
        # Modify band file
        band_file_epw(mpid,compound)
    elif condition == "pw2wan":
        # Generate input file for pw2wannier90.x
        projection = sys.argv[5]
        pw2wan_input(mpid,compound,prefix,proj=projection)
    elif condition == "kpathwan":
        # Generate k-path for plotting Wannierized bands
        if os.path.isfile("scf_dir/scf-{}-{}.in".format(mpid,compound)):
            file_name="scf_dir/scf-{}-{}.in".format(mpid,compound)
        else:
            print("scf-{}-{}.in not found in scf_dir/".format(mpid,compound) + "\n")
            print("Complete mainprogram.py from 0 to 6\n")
        kpoint_path(file_name)
    elif condition == "wankmesh":
        # Generate k-mesh for Wannierization
        prepare_nscf('kmesh.grid',what="wan")
    elif condition == "epw_band":
        # Prepare wannier90 input file
        projection = sys.argv[5]
        if dft in ('QE', 'qe'):
            epw_bandcheck("scf-{}-{}.in".format(mpid,compound),proj=projection,out="ex.win",dft=dft)
            os.system("mv ex.win epw_dir/ex-{}-{}.win".format(mpid,compound))
        else:
            epw_bandcheck("POSCAR",proj=projection,out="ex.win",dft=dft)
    elif condition == "epw":
        print("creating epw inputs for superconductivity\n")
        # Create EPW inputs for superconductivity
        projection = sys.argv[5]
        epw_sc(mpid,compound,prefix,'qpoint.dat','kpoint.dat',proj=projection,out="epw.in")
    else:
        print("bad inputs")
if __name__ == "__main__":
    main()
