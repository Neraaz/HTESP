#!/usr/bin/env python
"""The four long help blocks printed by ``mainprogram``.

The original dispatcher spent 165 ``print()`` calls over 170 lines on this
text, which made ``mainprogram.py`` hard to read and let ``docs/command.rst``
drift away from it (the documentation still described ``mainprogram 20`` --
``clean-scan`` -- as the partial-DOS command).  Keeping the text here as data
means

* the dispatcher is a dict lookup,
* the prose can be edited as prose, and
* ``docs/command.rst`` is generated from these constants by
  ``tools/gen_command_rst.py``, so the two can no longer disagree.

Each value is printed verbatim by ``mainprogram <name>``.
"""
from __future__ import annotations

HELP: dict[str, str] = {}

HELP['basicinfo'] = r"""
___________________________________________________________

**************Basic instructions**********************************

# run 'mainprogram process'

A fully-populated config.json ships with the package; 'mainprogram config-init'
writes a copy of it here to edit, and 'mainprogram config-validate'
shows which file is in use and what is wrong with it.  To start from a copy, take
utility/input_files/config.json into the working directory and edit it.  Any key you
leave out falls back to the packaged default, so an older config.json still works.
Set the Materials Project API key in the environment: export MP_API_KEY=<your key>

process = jobscript generates the job scripts for the calculations

process = search, search for data in materials project database

For process = download, download QE and VASP input files
process = oqmd-search, search for data in oqmd database

For process = oqmd-download, download QE and VASP input files

process = aflow-search, search for data in aflow database

process = aflow-download, download QE and VASP input files

process = data-combine, combining and eliminating duplicate inputs for different database

For information about QE+VASP calculations, process = process-info

process = epw-info for EPW calculations

process = wt-info for wanniertools calculations

process = elastic-input, to create input files with deformations

process = compute-elastic, to compute elastic properties

process = magenum, to create vasp input files for different magnetic state

process = magmom_extract, extract magnetic moment

process = fermisurface to plot fermi surface from vasprun.xml

process = charge-input for creating input files for system with non-zero net charge

process = pressure-input for creating input files for different pressure
For QE, 'pressure.in' file is provided with v1 pressure1, v2 pressure2 on different lines

For pressure calculations, mpid-pressure.in file is created and pressure value is inserted to scf_dir/scf-mpid.in files to get scf_dir/scf-mpid-pressure.in

For vasp, 'pressure.in' file has scaling factor for isotropic volume change with v1 scale1, v2 scale 2, on different lines, where scale can be 0.94, 0.96, ..
This change scaling factor for cell in POSCAR

For information about compound, Run after relaxation. process = compound

To check the status of el-ph calculation, use process = checkph

To check negative frequency in the phonon band, use process = checkfreq

process = singlemode provides info for single-mode phonon calculations. Also do 'mainprogram process-info' and look for process 23-25 for automated calculations. Requires run-dynmat.sh, and run-scf.sh files in working directory

kmesh can be changed using process=change_k

Use kpoint.in similar to qpoint.in file. Either provide new kmesh or fractional number to scale old k-mesh

process = history to print latest 10 mainprogram command executed

First execute 'history -a' in command line before process = history.

For MacOs, replace it by ~/.zsh_history in the mainprogram file

####################################################################
Please refer to the Online Documentation for further information and guidance.

##################################################################
"""

HELP['process-info'] = r"""
******************************************************************************************************************************************************************

Follow these instruction, start calculations with process number

**********************************************************************

Adjust start and end in input.in according to mpid-list.in

*************************************************************************

process = e0, to extract the total energies per atom and store in econv.csv file (QE+VASP)

For process = 1, relax-scan. This will relax the structure for the first time (QE+VASP).

Now a Rmpid-compound and Rmpid-compound/relax folders are created

process = 2 updates the input file with new structure (QE+VASP)

For process = 2, further-relax-input.
process = 3 resubmit the relaxation with updated input files

For process = 3, further-relax-scan

Repeat process = 2 and 3 for more relaxation (QE+VASP)

For process = 4, create-inputs, (QE).
Default: qmesh = kmesh/2 along each direction for el-ph calculations

provide qpoint.in file to provide qpoint mesh for phonon calculation

All the necessary inputs are created inside folders scf_dir,matdyn_dir,elph_dir,q2r_dir,kpath

For process = 5, fine-scan, This perform scf calculations with fine k grid

Now a Rmpid-compound/calc folder is created (QE)

For process = 6, coarse-scan, performs scf calculations with a coarse k grid (QE)

For process = 7, ph-scan. Performs ELECTRO-PHONON coupling (EPC) calculations (QE)

For process = 8, q2r-scan (QE)

For process = 9, matdyn-scan (QE)

For process = 10, matdyn-dos-scan. Phonon DOS calculation (QE)

For process = 11, lambda-scan (QE)

For process = 12, phonband-scan. Processing phonon dos (QE)

For process = 13, bandscf-scan (QE+VASP)

Rmpid-compound/bands folder is created for electronic bandstructure and density of states calculations

For process = 14, band-scan. NonSCF band Structure calculation (QE)

For process = 15, bandp-scan. Processing Bandstructure data (QE+VASP)

For process = 16, dos-scan. eDOS calculations (QE+VASP)

For process = 17, dosp-scan. Processing totalDOS (QE)

For process = 18, pdos-scan Processing partial DOS (QE)

For process = 19, plot-scan (QE+VASP)

For process = 20, clean-scan, Removing wavefunctions and bulky folders (QE)

For process = 21, extract-scan, Extracting EPC results and store in result.csv file (QE)

For process = 22,  Extracting total energy (QE+VASP)
for different plane wave cutoff and kpoint and store in {param}-{id}-{name}.txt file. run after process=convtest

For process = 23, dynmat-scan, Obtain atomic displacement files (QE) for vibrational mode at Gamma point

For process = 24, distortion-relax-scan, relaxing distorted structure (QE)

For process = 25, distortion-energy-scan, collecting distorted structure relaxation results (QE)

For process = 26, pressure-relax-scan, SCF calculations for different pressure (QE+VASP).
Use 'pressure.in' file with v1 pressure1, v2 pressure2, .... in different line

For process = 27, pressure-ph-scan, phonon calculation for different pressure (QE).
Create input file ph-mpid-compound.in with 'mainprogram epw1'.
ph-q.in file is provided for phonon calculation at particular q point,
otherwise, provide qpoint.in file for direct generic phonon calculation.
File ph-q.in file has nq1 nq2 nq3 and metal info on different line.
if T or t are used, calculation is performed for metal.

For process = 28, delete pressure folder (QE)

For process = 29, element substitution. (QE+VASP)

For process = convtest, perform convergence tests for Ecut and kpoint mesh (QE+VASP)

For process = compound, Get info about compounds (QE+VASP)

process = primtoconv, to change structure into conventional unit cell
,useful for phonopy calculations

process = pd, computing thermodynamic stability (QE+VASP) using pymatgen with 'econv.csv' file

process = phono1, to make supercell (QE+VASP) and submit scf calculations for different displacement

process = phono2, computing force constant (QE+VASP)

process = phono3, computing and plotting thermodynamic properties (QE+VASP)

process = phono4, computing and plotting phonon band (QE+VASP)

process = phono5, printing symmetry analysis (QE+VASP)

process = phono-qha, Computing temperature and pressure dependent thermal properties

process = ev-collect, extracting the total energies (QE+VASP) for different isotropic volumes from VASP calculations

Do 'mainprogram 26' calculation before process = ev-collect (QE+VASP)

process = phono1-pressure, (QE+VASP) submit phono1 calculations for different isotropic volumes

process = phono2-pressure, (QE+VASP) submit phono2 calculations for different isotropic volumes

process = phono3-pressure, (QE+VASP) submit phono3 calculations for different isotropic volumes

process = phono4-pressure, (QE+VASP) submit phono4 calculations for different isotropic volumes

process = eos-bm, equation of state fitting using Birch-Murnaghan fit (QE+VASP)

process = eos-vinet, equation of state fitting using vinet fit (QE+VASP)

**********************************************************************************
*********************************************************************************
"""

HELP['epw-info'] = r"""
perform relaxation and ground-state calculations with process from 1 - 4

Locate wannier90.json file in utility/input_files, copy to working directory, and adjust parameters

mainprogram epw1 ==> preparing input files for scf, non-scf, phonon calculations

mainprogram qe-ph ==> scf and phonon calculations

mainprogram epw2 ==> copy phonon files in save directory

mainprogram epw3 ==> projection calculations for scdm projection

mainprogram epw4 ==> fitting procedure to obtain scdm parameters

mainprogram wann-scdm ==> preparing input files for wannierization using scdm projections

mainprogram wann-file ==> preparing input files for wannierization taking projections from projection.in file

mainprogram wann-random ==> preparing input files for wannierization using random projections

mainprogram 13-18 for bandstructure and DOS calculations for analyzing and determining different windows

mainprogram epw5 ==> preparing inputfiles for QE bandstructure calculation using kpoints from wannier calculation (to obtain bands on same k-points)

mainprogram epw-scdm, epw-file, epw-random ==> preparing input files for epw calculations (anisotropic Eliashberg-Migdel approximations) with different projection schemes
"""

HELP['wt-info'] = r"""
************************************************************************************

First repeat all the calculations as described in 'mainprogram epw-info' command upto wannierization

process = wt1, prepare input file wt.in required for initial bulk bandgap calculation

 if not found, it will create a default one

copy 'wt-mpid-compound.in' file from 'WT_dir' to
Rmpid-compound/epw/ folder where wannierization process was done

Please include slab dimension even in bulk calculation, so that it produces 'POSCAR-slab' file which is used by ASE package to create 'KPATH_SLAB' for slab system

process = wt2, prepare input file for other calculations including surfaces

Edit 'wanniertools_input' key in config.json according to properties of interest

************************************************************************************
"""


HELP = {name: text.lstrip("\n") for name, text in HELP.items()}


#: short one-line description per command, used by ``mainprogram --list`` and
#: by ``tools/gen_command_rst.py``.
SUMMARY: dict[str, str] = {
    "basicinfo": "print the basic instructions",
    "process-info": "print the numbered-process reference",
    "epw-info": "print the EPW / Wannier90 pipeline reference",
    "wt-info": "print the WannierTools pipeline reference",
    "jobscript": "generate the batch submission scripts",
    "search": "search the Materials Project database",
    "download": "create input files for the search results",
    "oqmd-search": "search the OQMD database",
    "oqmd-download": "create input files from OQMD results",
    "aflow-search": "search the AFLOW database",
    "aflow-download": "create input files from AFLOW results",
    "data-combine": "merge the three databases, dropping duplicates",
    "elastic-input": "create deformed structures for elastic constants",
    "compute-elastic": "compute the elastic tensor and derived moduli",
    "magenum": "enumerate magnetic configurations",
    "magmom_extract": "extract magnetic moments from the finished runs",
    "fermisurface": "plot the 3D Fermi surface with IFermi",
    "charge-input": "create inputs for a non-zero net charge",
    "pressure-input": "create inputs for a pressure or volume series",
    "compound": "print structural information for the tracked compounds",
    "checkph": "classify the state of every electron-phonon run",
    "checkfreq": "flag imaginary frequencies in the phonon bands",
    "singlemode": "help for single-mode phonon calculations",
    "change_k": "rewrite the k-mesh from kpoint.in",
    "primtoconv": "convert the relaxed cells to conventional settings",
    "history": "print the last 10 mainprogram commands",
    "convtest": "prepare the cutoff and k-mesh convergence tests",
    "pd": "compute thermodynamic stability from econv.csv",
    "e0": "collect total energies per atom into econv.csv",
    "config-init": "write the packaged default config.json here",
    "config-validate": "check config.json and report problems",
    "phono1": "build the phonopy supercells and submit them",
    "phono2": "compute the force constants",
    "phono3": "compute and plot the thermodynamic properties",
    "phono4": "compute and plot the phonon band structure",
    "phono5": "print the phonopy symmetry analysis",
    "phono-qha": "temperature and pressure dependent thermal properties",
    "ev-collect": "collect energies versus volume",
    "eos-bm": "Birch-Murnaghan equation-of-state fit",
    "eos-vinet": "Vinet equation-of-state fit",
    "phono1-pressure": "phono1 for every volume of the series",
    "phono2-pressure": "phono2 for every volume of the series",
    "phono3-pressure": "phono3 for every volume of the series",
    "phono4-pressure": "phono4 for every volume of the series",
    "epw1": "write the scf, nscf and phonon inputs for EPW",
    "qe-ph": "run the scf and phonon steps",
    "epw2": "copy the phonon files into the save directory",
    "epw3": "run the projection calculation for SCDM",
    "epw4": "fit the SCDM parameters",
    "epw5": "QE band structure on the Wannier k-points",
    "wann-scdm": "wannierisation with SCDM projections",
    "wann-file": "wannierisation with projections from projection.in",
    "wann-random": "wannierisation with random projections",
    "epw-scdm": "EPW (anisotropic Eliashberg) with SCDM projections",
    "epw-file": "EPW with projections from projection.in",
    "epw-random": "EPW with random projections",
    "wt1": "write wt.in for the bulk WannierTools calculation",
    "wt2": "write wt.in for the surface WannierTools calculations",
}

#: numbered processes -> one-line description
PROCESS_SUMMARY: dict[str, str] = {
    "0": "create scf_dir, elph_dir, matdyn_dir and q2r_dir",
    "1": "relax-scan: first structure relaxation (QE+VASP)",
    "2": "further-relax-input: rebuild the input from the relaxed structure",
    "3": "further-relax-scan: resubmit the relaxation",
    "4": "create-inputs: all downstream QE inputs from the relaxed structure",
    "5": "fine-scan: scf on the fine k-grid (QE)",
    "6": "coarse-scan: scf on the coarse k-grid (QE)",
    "7": "ph-scan: electron-phonon coupling with ph.x (QE)",
    "8": "q2r-scan: force constants in real space (QE)",
    "9": "matdyn-scan: phonon dispersion (QE)",
    "10": "matdyn-dos-scan: phonon density of states (QE)",
    "11": "lambda-scan: lambda.x (QE)",
    "12": "phonband-scan: post-process the phonon bands (QE)",
    "13": "bandscf-scan: stage the band-structure calculation (QE+VASP)",
    "14": "band-scan: non-scf band structure (QE)",
    "15": "bandp-scan: post-process the band structure (QE+VASP)",
    "16": "dos-scan: electronic density of states (QE+VASP)",
    "17": "dosp-scan: post-process the total DOS (QE)",
    "18": "pdos-scan: post-process the partial DOS (QE)",
    "19": "plot-scan: write the PDFs into plots/ (QE+VASP)",
    "20": "clean-scan: DELETE wavefunctions and move finished runs to completed/",
    "21": "extract-scan: collect Tc and lambda into result.csv (QE)",
    "22": "extract the convergence-test energies (run after convtest)",
    "23": "dynmat-scan: eigenmode displacements at Gamma (QE)",
    "24": "distortion-relax-scan: relax the distorted structures (QE)",
    "25": "distortion-energy-scan: collect the distortion energies (QE)",
    "26": "pressure-relax-scan: scf for every pressure or volume (QE+VASP)",
    "27": "pressure-ph-scan: phonons for every pressure (QE)",
    "28": "pressure-reset: remove the pressure folders (QE)",
    "29": "sitesub-scan: element substitution (QE+VASP)",
}
