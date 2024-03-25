#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Main command line scripts"""
import os
import sys
#import glob
import json
from cif_to_gsinput import pos_to_kpt
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
def main():
    """
    main function
    """
    with open("log", "w") as flog:
        flog.write("_____________________________________________________________________________________________________________________________________\n")
        flog.write("\n")
        flog.write("\n")
        flog.write("                        ****    ****   ***********   *********   ********       ********                 \n")
        flog.write("                        |  |    |  |       | |       | |____     | |            |  |  \\ \\                 \n")
        flog.write("                        |  |____|  |       | |       |  ____|    | |*****       |  |__| |                       \n")
        flog.write("                        |   ____   |       | |       | |_____          | |      |  _____/                       \n")
        flog.write("                        |  |    |  |       | |       |_______|   ******| |      |_ |                       \n")
        flog.write("                        |__|    |__|       |_| *********************************************                                          \n")
        flog.write("                        **********************                                                               \n")
        flog.write("                                      High Throughput Electron-Structure Package                              \n")
        flog.write("\n")
        flog.write("                                              Program written by\n")
        flog.write("\n")
        flog.write("                                Niraj K Nepal, PhD       &       Lin-Lin Wang, PhD                        \n")
        flog.write("                          Email: nnepal@ameslab.gov                              \n")
        flog.write("                                 tug11655@temple.edu                              \n")
        flog.write("\n")
        flog.write("\n")
        flog.write("_____________________________________________________________________________________________________________________________________\n")
        flog.write("\n")
        flog.write("\n")
        if os.path.isfile('input.in'):
            with open('input.in', 'r') as read_in:
                lines = read_in.readlines()
            if len(lines) >= 4:
                start = int(lines[0])
                end = int(lines[1])
                kpt = lines[2].split("\n")[0].split(" ")
                if len(kpt) > 1:
                    nkpt = int(kpt[0])
                else:
                    nkpt = int(lines[2])
                element = lines[3].split()[0]
                plot_type = lines[4].split("\n")[0].split(" ")
            else:
                print("Input file has less than 4 lines\n")
                print("Please type 'mainprogram inputinfo' for detail\n")
        else:
            print("No input.in file detected, creating one\n")
            with open('input.in', 'w') as input_file:
                start=1
                end=2
                nkpt=200
                element='mpid-list.in'
                plot_type='phband'
                input_file.write(str(start)+"\n")
                input_file.write(str(end) + "\n")
                input_file.write(str(nkpt) + " " + str(0) + "\n")
                input_file.write(element + "\n")
                input_file.write(plot_type + "\n")
                input_file.write(f"DFT = {input_data['download']['inp']['calc']}\n")
        flog.write("#########################################################################")
        flog.write("#######################################################################\n")
        flog.write("#start: {}, end: {}, tracking_file: {}, nkpoint: {}".format(start,end,element,nkpt) + "#\n")
    process = sys.argv[1]
    if not process.isdigit():
        if process == 'compound':
            os.system("info-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 'jobscript':
            os.system("generate_submission.py")
        elif process == 'download':
            download = input_data['download']
            if download['info']['mode'] in ('element', 'chemsys'):
                os.system("download-input" + " " + str(start) + " " + str(end) + " " + element)
            elif download['info']['mode'] == 'fromcif':
                cif2cell=download['inp']['use_cif2cell']
                print("\n")
                print("--------------------------------------------------------------\n")
                print("processing .cif files for {} calculations".format(download['inp']['calc']) + "\n")
                if cif2cell:
                    print("Install cif2cell package using 'pip install cif2cell'\n")
                # Check CIF2CELL parameters in cif_to_gsinput.py
                # CIF2CELL if True, uses cif2cell package to create POSCAR from given .cif files.
                # if False, It uses pymatgen cifparser to read and produce cif output, which then explicitely
                # read and uses to write POSCAR.
                os.system("cif_to_gsinput.py {}".format(download['inp']['calc']))
                if download['inp']['calc'] in ('VASP','vasp'):
                    print("Check Rmp-mpid-compound folders for inputs and mpid.in tracking file\n")
                if download['inp']['calc'] in ('QE','qe'):
                    print("Check scf_dir folder for inputs and mpid.in tracking file\n")
                print("--------------------------------------------------------------\n")
            else:
                print("\n")
                print("--------------------------------------------------------------\n")
                print("processing .vasp files for {} calculations".format(download['inp']['calc']) + "\n")
                os.system("poscar_to_vasp.py")
        elif process == 'convtest':
            os.system("convergence_test.py calculate")
        elif process == 'pressure-input':
            os.system("pressure-input" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 'charge-input':
            os.system("charge-input" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 'singlemode':
            os.system("distortion.sh")
        elif process == 'fermisurface':
            os.system("ifermi-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 'oqmd-search':
            os.system("oqmd_extract.py search")
        elif process == 'oqmd-download':
            os.system("oqmd_extract.py download")
        elif process == 'aflow-search':
            os.system("aflow_extract.py search")
        elif process == 'aflow-download':
            os.system("aflow_extract.py download")
        elif process == 'data-combine':
            os.system("structure_group.py")
        elif process == 'elastic-input':
            os.system("elastic.py input")
        elif process == 'compute-elastic':
            os.system("elastic.py compute_elastic")
        elif process == 'magenum':
            os.system("magnetic.py")
        elif process == 'search':
            if os.path.isfile('htepc.json') or os.path.isfile("../../htepc.json"):
                download = input_data['download']
                mode = download['info']['mode']
                print("--------------------------------------------------------------\n")
                print("'{}' mode found".format(download['info']['mode']) + "\n")
            else:
                mode = 'element'
            if mode not in ('fromcif','fromvasp'):
                print("--------------------------------------------------------------\n")
                print("   Searching compounds in materials project database          \n")
                print("--------------------------------------------------------------\n")
            if os.path.isfile("mpid-list.in"):
                print("A mpid-list.in file found, renaming mpid-list-2.in\n")
                os.system("mv mpid-list.in mpid-list-2.in")
            os.system("element_extract.py")
            if mode not in ('fromcif','fromvasp'):
                print("--------------------------------------------------------------\n")
                print("Change indices in input.in and execute 'mainprogram download'\n")
                print("to download input files\n")
                print("--------------------------------------------------------------\n")
            else:
                print("Execute 'mainprogram download' to prepare input files from .cif or .vasp files\n")
        elif process == 'basicinfo':
            print("___________________________________________________________\n")
            print("**************Basic instructions**********************************\n")
            print("\n")
            print("# run 'mainprogram process'\n")
            print("Look for htepc.json file in utility/input_files/ and copy that to the working directory\n")
            print("process = jobscript generates the job scripts for the calculations\n")
            print("process = search, search for data in materials project database\n")
            print("For process = download, download QE and VASP input files")
            print("process = oqmd-search, search for data in oqmd database\n")
            print("For process = oqmd-download, download QE and VASP input files\n")
            print("process = aflow-search, search for data in aflow database\n")
            print("process = aflow-download, download QE and VASP input files\n")
            print("process = data-combine, combining and eliminating duplicate inputs for different database\n")
            print("For information about QE+VASP calculations, process = process-info\n")
            print("process = epw-info for EPW calculations \n")
            print("process = wt-info for wanniertools calculations \n")
            print("process = vasp-info for vasp+phonon calculations \n")
            print("process = elastic-input, to create vasp input files with deformations\n")
            print("process = compute-elastic, to compute elastic properties\n")
            print("process = magenum, to create vasp input files for different magnetic state\n")
            print("process = magmom_extract, extract magnetic moment\n")
            print("process = fermisurface to plot fermi surface from vasprun.xml\n")
            print("process = charge-input for creating input files for system with non-zero net charge\n")
            print("process = pressure-input for creating input files for different pressure")
            print("For QE, 'pressure.in' file is provided")
            print(" with v1 pressure1, v2 pressure2 on different lines\n")
            print("For pressure calculations, mpid-pressure.in file is created and")
            print(" pressure value is inserted to scf_dir/scf-mpid.in files to get scf_dir/scf-mpid-pressure.in\n")
            print("For vasp, 'pressure.in' file has scaling factor for isotropic volume change")
            print(" with v1 scale1, v2 scale 2, on different lines, where scale can be 0.94, 0.96, .. ")
            print("This change scaling factor for cell in POSCAR\n")
            print("For information about compound, Run after relaxation. process = compound\n")
            print("To check the status of el-ph calculation, use process = checkph\n")
            print("To check negative frequency in the phonon band, use process = checkfreq\n")
            print("process = singlemode provides info for single-mode phonon calculations.")
            print(" Also do 'mainprogram process-info' and look for process 23-25 for automated calculations.")
            print(" Requires run-dynmat.sh, run-scf.sh, and Vasp.pm files in working directory\n")
            vasp_potcar="""#For POTCAR. Suppose we have POTCARS as  POT_GGA_PAW_PBE/Mg_p/POTCAR
                          #Mg_p can be found in POTCAR.spec
                          #pmg config -p /path_to/POT_GGA_PAW_PBE PBE52
                          # After that add path to .pmgrc.yaml
                          #pmg config --add PMG_VASP_PSP_DIR PBE52"""
            print(vasp_potcar + "\n")
            print("kmesh can be changed using process=change_k\n")
            print("Use kpoint.in similar to qpoint.in file.")
            print(" Either provide new kmesh or fractional number to scale old k-mesh\n")
            print("process = history to print latest 10 mainprogram command executed\n")
            print("First execute 'history -a' in command line before process = history.\n")
            print("For MacOs, replace it by ~/.zsh_history in the mainprogram file\n")
            print("####################################################################")
            print("Please refer to the Online Documentation for further information and guidance.\n") 
            print("##################################################################\n")
        elif process == 'vasp-info':
            print("Strucutural relaxation, substitution, pressure, magnetic orderings")
            print(" (isotropic, changing scaling factor of lattice) can be performed")
            print(" similary as of QE, but using DFT = vasp (or VASP) in input.in file\n")
            print("For process = 1 - 3, Structural relaxation similar to QE\n")
            print("process = 13 ==> Submit bandstructure calculations\n")
            print("process = 15 ==> processing bandstructure\n")
            print("process = 16 ==> Submit DOS and pDOS calculations\n")
            print("process = 19 ==> plotting bandstructure or DOS/pDOS")
            print(" Use 'eband', or  'pdos' in input.in\n")
            print("Thermodynamic quantities can be calculated using VASP and phonopy (vp-ph)\n")
            print("process = primtoconv, to change structure into conventional unit cell")
            print(",useful for vasp+phonopy calculations\n")
            print("process = vp-pd, computing thermodynamic stability")
            print(" using pymatgen with 'econv_vasp.csv' file\n")
            print("process = phono1, to make supercell")
            print(" and submit scf calculations for different displacement\n")
            print("process = phono2, computing force constant\n")
            print("process = phono3, computing and plotting thermodynamic properties\n")
            print("process = phono4, computing and plotting phonon band\n")
            print("process = phono5, printing symmetry analysis\n")
            print("process = phono-qha,")
            print(" Computing temperature and pressure dependent thermal properties\n")
            print("process = ev-collect, extracting the total energies")
            print(" for different isotropic volumes from VASP calculations\n")
            print("Do 'mainprogram 26' calculation before process = ev-collect\n")
            print("process = phono1-pressure,")
            print(" submit phono1 calculations for different isotropic volumes\n")
            print("process = phono2-pressure,")
            print(" submit phono2 calculations for different isotropic volumes\n")
            print("process = phono3-pressure,")
            print(" submit phono3 calculations for different isotropic volumes\n")
            print("process = phono4-pressure,")
            print(" submit phono4 calculations for different isotropic volumes\n")
            print("process = eos-bm, equation of state fitting using Birch-Murnaghan fit\n")
            print("process = eos-vinet, equation of state fitting using vinet fit\n")
            print("process = vp-ph-help, printing help page\n")
        elif process == 'phono-qha':
            os.system("phonopy-scan" + " " + str(start) + " " + str(end) + " " + element + " " + "vp-ph-qha")
        elif process == 'eos-bm':
            os.system("phonopy-scan" + " " + str(start) + " " + str(end) + " " + element + " " + "eos-bm")
        elif process == 'eos-vinet':
            os.system("phonopy-scan" + " " + str(start) + " " + str(end) + " " + element + " " + "eos-vinet")
        elif process == 'ev-collect':
            os.system("phonopy-scan" + " " + str(start) + " " + str(end) + " " + element + " " + "ev-collect")
        elif process == 'vp-pd':
            os.system("pymatgen_phase_diagram.py")
        elif process == 'primtoconv':
            with open(element, "r") as read_elm:
                lines = read_elm.readlines()
            folders = lines[start-1:end-1]
            currdir = os.getcwd()
            for fold in folders:
                fold_list = fold.split("\n")[0].split(" ")
                fld = "R{}-{}".format(fold_list[1],fold_list[2])
                os.chdir(currdir + "/" + fld)
                os.system("cp -r relax/ relax_prim/")
                os.chdir(currdir + "/" + fld + "/relax")
                os.system("vasp_process.py conventional")
                os.chdir("../../")
        elif process == 'e0':
            os.system("phonopy-scan" + " " + str(start) + " " + str(end) + " " + element + " " + str(0))
        elif process == 'phono1':
            os.system("phonopy-scan" + " " + str(start) + " " + str(end) + " " + element + " " + str(1))
        elif process == 'phono2':
            os.system("phonopy-scan" + " " + str(start) + " " + str(end) + " " + element + " " + str(2))
        elif process == 'phono3':
            os.system("phonopy-scan" + " " + str(start) + " " + str(end) + " " + element + " " + str(3))
        elif process == 'phono4':
            os.system("phonopy-scan" + " " + str(start) + " " + str(end) + " " + element + " " + str(4))
        elif process == 'phono5':
            os.system("phonopy-scan" + " " + str(start) + " " + str(end) + " " + element + " " + str(5))
        elif process == 'phono1-pressure':
            os.system("phonopy-scan" + " " + str(start) + " " + str(end) + " " + element + " " + "vp-ph2-pressure")
        elif process == 'phono2-pressure':
            os.system("phonopy-scan" + " " + str(start) + " " + str(end) + " " + element + " " + "vp-ph3-pressure")
        elif process == 'phono3-pressure':
            os.system("phonopy-scan" + " " + str(start) + " " + str(end) + " " + element + " " + "vp-ph4-pressure")
        elif process == 'phono4-pressure':
            os.system("phonopy-scan" + " " + str(start) + " " + str(end) + " " + element + " " + "vp-ph5-pressure")
        elif process == "history":
            os.system("history.sh")
        elif process == 'checkph':
            os.system("phcheck-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 'checkfreq':
            os.system("checkfreq-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 'change_k':
            os.system("double-kmesh" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 'magmom_extract':
            os.system("magmom-extract" + " " + str(start) + " " + str(end) + " " + element)
        elif process == "epw1":
            os.system("epw-bash-scripts" + " " + str(start) + " " + str(end) + " " + element + " " + "epw1")
        elif process == "qe-ph":
            os.system("epw-bash-scripts" + " " + str(start) + " " + str(end) + " " + element + " " + "epw2")
        elif process == "epw2":
            os.system("epw-bash-scripts" + " " + str(start) + " " + str(end) + " " + element + " " + "epw3")
        elif process == "epw3":
            os.system("epw-bash-scripts" + " " + str(start) + " " + str(end) + " " + element + " " + "epw4")
        elif process == "epw4":
            os.system("epw-bash-scripts" + " " + str(start) + " " + str(end) + " " + element + " " + "proj")
        elif process == "wann-scdm":
            os.system("epw-bash-scripts" + " " + str(start) + " " + str(end) + " " + element + " " + "band_wann" + " " + "scdm")
        elif process == "wann-file":
            os.system("epw-bash-scripts" + " " + str(start) + " " + str(end) + " " + element + " " + "band_wann" + " " + "fromfile")
        elif process == "wann-random":
            os.system("epw-bash-scripts" + " " + str(start) + " " + str(end) + " " + element + " " + "band_wann" + " " + "random")
        elif process == "epw5":
            os.system("epw-bash-scripts" + " " + str(start) + " " + str(end) + " " + element + " " + "band_wann2")
        elif process == "epw-scdm":
            os.system("epw-bash-scripts" + " " + str(start) + " " + str(end) + " " + element + " " + "epw" + " " + "scdm")
        elif process == "epw-file":
            os.system("epw-bash-scripts" + " " + str(start) + " " + str(end) + " " + element + " " + "epw" + " " + "fromfile")
        elif process == "epw-random":
            os.system("epw-bash-scripts" + " " + str(start) + " " + str(end) + " " + element + " " + "epw" + " " + "random")
        elif process == "epw-info":
            print("perform relaxation and ground-state calculations with process from 1 - 4\n")
            print("Locate wannier90.json file in utility/input_files, copy to working directory, and adjust parameters\n")
            print("mainprogram epw1 ==> preparing input files for scf, non-scf, phonon calculations \n")
            print("mainprogram qe-ph ==> scf and phonon calculations\n")
            print("mainprogram epw2 ==> copy phonon files in save directory\n")
            print("mainprogram epw3 ==> projection calculations for scdm projection\n")
            print("mainprogram epw4 ==> fitting procedure to obtain scdm parameters\n")
            print("mainprogram wann-scdm ==> preparing input files")
            print(" for wannierization using scdm projections\n")
            print("mainprogram wann-file ==> preparing input files")
            print(" for wannierization taking projections from projection.in file\n")
            print("mainprogram wann-random ==> preparing input files")
            print(" for wannierization using random projections\n")
            print("mainprogram 13-18 for bandstructure and DOS calculations")
            print(" for analyzing and determining different windows\n")
            print("mainprogram epw5 ==> preparing inputfiles for QE bandstructure")
            print(" calculation using kpoints from wannier calculation (to obtain bands on same k-points)\n")
            print("mainprogram epw-scdm, epw-file, epw-random ==> preparing input")
            print(" files for epw calculations (anisotropic Eliashberg-Migdel approximations)")
            print(" with different projection schemes \n")
        elif process == "wt-info":
            print("************************************************************************************\n")
            print("First repeat all the calculations as described in 'mainprogram epw-info'")
            print(" command upto wannierization\n")
            print("process = wt1, prepare input file wt.in required for initial bulk bandgap calculation\n")
            print(" if not found, it will create a default one\n")
            print("copy 'wt-mpid-compound.in' file from 'WT_dir' to")
            print("Rmpid-compound/epw/ folder where wannierization process was done\n")
            print("Please include slab dimension even in bulk calculation,")
            print(" so that it produces 'POSCAR-slab' file")
            print(" which is used by ASE package to create 'KPATH_SLAB' for slab system\n")
            print("process = wt2, prepare input file for other calculations including surfaces\n")
            print("Edit 'wanniertool_input' key in htepc.json according to properties of interest\n")
            print("************************************************************************************\n")
        elif process == "wt1":
            os.system("wt-bash-scripts" + " " + str(start) + " " + str(end) + " " + element + " " + "wt1-b")
        elif process == "wt2":
            os.system("wt-bash-scripts" + " " + str(start) + " " + str(end) + " " + element + " " + "wt1-s")
        elif process == "process-info":
            print("******************************************************************************************************************************************************************\n")
            print("Follow these instruction, start calculations with process number\n")
            print("**********************************************************************\n")
            print("Adjust start and end in input.in according to mpid-list.in\n")
            print("*************************************************************************\n")
            print("process = e0, to extract the total energies per atom and store in econv_vasp.csv file (QE+VASP)\n")
            print("For process = 1, relax-scan. This will relax the structure for the first time (QE+VASP).\n")
            print("Now a Rmpid-compound and Rmpid-compound/relax folders are created\n")
            print("process = 2 updates the input file with new structure (QE+VASP)\n")
            print("For process = 2, further-relax-input.")
            print("process = 3 resubmit the relaxation with updated input files\n")
            print("For process = 3, further-relax-scan \n")
            print("Repeat process = 2 and 3 for more relaxation (QE+VASP)\n")
            print("For process = 4, create-inputs, (QE).")
            print("Default: qmesh = kmesh/2 along each direction for el-ph calculations\n")
            print("provide qpoint.in file to provide qpoint mesh for phonon calculation\n")
            print("All the necessary inputs are created")
            print(" inside folders scf_dir,matdyn_dir,elph_dir,q2r_dir,kpath\n")
            print("For process = 5, fine-scan, This perform scf calculations with fine k grid\n")
            print("Now a Rmpid-compound/calc folder is created (QE)\n")
            print("For process = 6, coarse-scan, performs scf calculations with a coarse k grid (QE)\n")
            print("For process = 7, ph-scan. Performs ELECTRO-PHONON coupling (EPC) calculations (QE)\n")
            print("For process = 8, q2r-scan (QE)\n")
            print("For process = 9, matdyn-scan (QE)\n")
            print("For process = 10, matdyn-dos-scan. Phonon DOS calculation (QE)\n")
            print("For process = 11, lambda-scan (QE)\n")
            print("For process = 12, phonband-scan. Processing phonon dos (QE)\n")
            print("For process = 13, bandscf-scan (QE+VASP)\n")
            print("Rmpid-compound/bands folder is created")
            print(" for electronic bandstructure and density of states calculations\n")
            print("For process = 14, band-scan. NonSCF band Structure calculation (QE)\n")
            print("For process = 15, bandp-scan. Processing Bandstructure data (QE+VASP)\n")
            print("For process = 16, dos-scan. eDOS calculations (QE+VASP)\n")
            print("For process = 17, dosp-scan. Processing totalDOS (QE)\n")
            print("For process = 18, pdos-scan Processing partial DOS (QE)\n")
            print("For process = 19, plot-scan (QE+VASP)\n")
            print("For process = 20, clean-scan, Removing wavefunctions and bulky folders (QE)\n")
            print("For process = 21, extract-scan, Extracting EPC results and store in result.csv file (QE)\n")
            print("For process = 22,  Extracting total energy (QE+VASP)")
            print("for different plane wave cutoff and kpoint and store in {param}-{id}-{name}.txt file.")
            print(" run after process=convtest\n")
            print("For process = 23, dynmat-scan, Obtain atomic displacement files (QE)")
            print(" for vibrational mode at Gamma point\n")
            print("For process = 24, distortion-relax-scan, relaxing distorted structure (QE)\n")
            print("For process = 25, distortion-energy-scan,")
            print(" collecting distorted structure relaxation results (QE)\n")
            print("For process = 26, pressure-relax-scan, SCF calculations for different pressure (QE+VASP).")
            print("Use 'pressure.in' file with v1 pressure1, v2 pressure2, .... in different line\n")
            print("For process = 27, pressure-ph-scan, phonon calculation for different pressure (QE).")
            print("Create input file ph-mpid-compound.in with 'mainprogram epw1'.")
            print("ph-q.in file is provided for phonon calculation at particular q point,")
            print("otherwise, provide qpoint.in file for direct generic phonon calculation.")
            print("File ph-q.in file has nq1 nq2 nq3 and metal info on different line.")
            print("if T or t are used, calculation is performed for metal.\n")
            print("For process = 28, delete pressure folder (QE)\n")
            print("For process = 29, element substitution. check 'site_subs.py h' (QE+VASP)\n")
            print("For process = 30, creating inputs for different magnetic orderings (VASP)\n")
            print("For process = convtest, perform convergence tests for Ecut and kpoint mesh (QE+VASP)\n")
            print("For process = compound, Get info about compounds (QE+VASP)\n")
            print("**********************************************************************************")
            print("*********************************************************************************\n")
        else:
            print("Bad input, Do 'mainprogram basicinfo'\n")
    else:
        process = int(process)
        if process == 0:
            if not os.path.isdir('scf_dir') or not os.path.isdir("elph_dir") or not os.path.isdir("matdyn_dir") or not os.path.isdir("q2r_dir"):
                os.system("mkdir scf_dir elph_dir matdyn_dir q2r_dir")
        elif process == 1:
            os.system("relax-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 2:
            os.system("further-relax-input" + " " + str(start) + " " + str(end) + " " + element + " first")
        elif process == 3:
            os.system("further-relax-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 4:
            os.system("create-inputs" + " " + str(start) + " " + str(end) + " " + element + " " + str(nkpt))
        elif process == 5:
            os.system("fine-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 6:
            os.system("coarse-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 7:
            os.system("ph-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 8:
            os.system("q2r-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 9:
            os.system("matdyn-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 10:
            os.system("matdyn-dos-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 11:
            os.system("lambda-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 12:
            os.system("phonband-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 13:
            os.system("bandscf-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 14:
            os.system("band-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 15:
            os.system("bandp-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 16:
            os.system("dos-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 17:
            os.system("dosp-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 18:
            os.system("pdos-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 19:
            for plot in plot_type:
                print("plotting: {}".format(plot) + "\n")
                os.system("plot-scan" + " " + str(start) + " " + str(end) + " " + element + " " + str(nkpt) + " " + plot)
        elif process == 20:
            os.system("clean-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 21:
            os.system("extract-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 22:
            os.system("convergence_test.py extract")
        elif process == 23:
            os.system("dynmat-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 24:
            os.system("distortion-relax-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 25:
            os.system("distortion-energy-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 26:
            os.system("pressure-relax-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 27:
            os.system("pressure-ph-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 28:
            os.system("pressure-reset" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 29:
            os.system("sitesub-scan" + " " + str(start) + " " + str(end) + " " + element)
        elif process == 30:
            os.system("magnetic-scan" + " " + str(start) + " " + str(end) + " " + element)
        else:
            print("Use python mainprogram process, 0 <= process <= 29\n")
if __name__ == "__main__":
    main()
