#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D."""
MSG="""***********************************************************
*****************************************************************
Please follow following instructions.

We need "mpi_key.py" file with materials project API key.

Create a working directory (suppose work_dir)

cd work_dir
***************************************
Execute "mainprogram.py search" command
***************************************

It will create files download.py, input.in, and mpid-list.in files and a folder "download"

------------------------------------------------------------------------------------------------------------------------------
          process = search, files and  parameters
------------------------------------------------------------------------------------------------------------------------------
A. download.py ==> control file for inputs

    It has structure as follows.
  
    info={'mode': 'element', 'metal': True, 'FE': True, 'thermo_stable': False, 'exclude': ['Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Th','Pa','U','Np','Pu','O'], 'ntype': (2, 8), 'elm': ['B','C'], 'prop': ['material_id', 'formula_pretty', 'structure', 'formation_energy_per_atom', 'band_gap', 'energy_above_hull', 'total_magnetization', 'ordering', 'total_magnetization_normalized_formula_units', 'num_magnetic_sites', 'theoretical', 'composition_reduced'],'ordering':'NM','nsites':10}
inp={'start': 1, 'end': 2, 'nband': 200, 'plot': 'phband', 'calc': 'QE'}
chemsys={'entries': ['B', 'Mg'], 'size_constraint': 20, 'ntype_constraint': 5, 'must_include': ['Mg', 'B'],'form_en':False,'metal':False,'magnetic':False}
    
    
    We have 3 dictionaries "info", "inp", and chemsys. "info" controls the nature of the systems while "inp" controls the nature of calculations. If 'mode' is 'chemsys', then 3rd dictionary 'chemsys' is utilized to create inputs.
    
    Dictionary "info" parameters. 
     a. metal ==> True for zero band gap system
     b. FE ==> True for negative formation energy system
     c. Thermo_stable ==> True for system at the convex hull
     d. exclude ==> list of elements to exclude
     e. ntype ==> integer or tuple for different type of elements in the system
        (1,2) represents 2 different type of elements in the system
     f. prop ==> various properties to extract
     g. ordering ==> magnetic ordering
     h. mode ==> How to extract info from MP database. 
                 Available options
                 'element' : query from MP database with mpr.materials.search module. Suitable for searching with element
                 'chemsys' : query from MP database with mpr.get_entries_in_chemsys() function. Suitable for searching with elemental combinations
                 'fromcif' : check and convert provided .cif files to input files. Need .cif files inside working directory
     i. nsites ==> Total number of atoms in a system.
    
    Depending on these paramters a folder "download" is created with download.csv files within it. That has all the properties extracted according to "download.py" file. Also mpid-list.in with material id and compound name is created to track the system and calculations
    
    Dictionary "inp" parameters
    
     a. start ==> first index to chose system to download from the file mpid-list.in or download/download.csv
     b. end ==> last index to chose the system 
     c. nband ==> the total number of k-points along the high-symmetry path for band calculations
     d. plot ==> type of plots [options: phband (phonon band structure), eband (electronic band structure), gammaband (electron-phonon coupling strength projected phonon band), pdos (DOS and partial DOS), wann_band (wannier interpolated bandstructure), phonproj (compute atom-projected phonon bandstructure)]
     e. calc ==> type of DFT calculations (options: QE or VASP)

    Dictionary "chemsys" parameters
    
     a. entries ==> List of elements to search
     b. size_constraint ==> Total number of ions in the system
     c. ntype_constraint ==> number of distinct elements in the system
     d. must_include ==> Elements to include besides "entries"
     e. form_en: This parameter is a logical flag that indicates whether the compound should possess a negative formation energy.
     f. metal: This parameter is a logical flag used to determine if the compound is classified as a metal.
     g. magnetic: This parameter is a logical flag used to determine if the compound is classified as magnetic.

-----------------------------------------------------------------------------------------------------------------------------
          process = download, INPUT parameters for download
------------------------------------------------------------------------------------------------------------------------------
A. input.in is the file to control the calculation process.
   same paramters of dictionary "inp" is printed in different lines.
   ************************************
   format of input.in
   ************************************
   1
   2
   200 0
   mpid-list.in
   phband
   DFT = QE
   ************************************
   First and second line respectively start and end index to look in mpid-list.in, 3rd line (1st column) represents number of k-point meshes in high-symmetry path (from atomic simulation environment (ASE)) and 2nd column represents cutoff for full Brillouin zone path. 0 means using full path and 1 means ignoring last point.

B. mpid-list.in ==> tracking files in the format of "v1 materials_id compound_prefix" in different lines

C. We can tune download.py file according to our own need, and download input files with
   ******************************************
   Execute "mainprogram.py download" command
   ******************************************

------------------------------------------------------------------------------------------------------------------------------
          Quantum Espresso parameters
------------------------------------------------------------------------------------------------------------------------------

A. Kinetic energy cutoff can be employed using "pseudo.py" file. If not provided, default one is created

B. "pwscf_in.py" file contains the parameters for QE ground-state calculations. If not provided, default one is created.

------------------------------------------------------------------------------------------------------------------------------
          VASP parameters
------------------------------------------------------------------------------------------------------------------------------

For vasp, we can configure to look for potcars using pymatgen as,
provide pseudo.py file similar to htepc/utility/some_inputs_required/vasp_potcar.py
copy that file to pseudo.py

Suppose we have POTCARS as  POT_GGA_PAW_PBE/Mg/POTCAR 

pmg config -p ~/path_along_directory_POT_GGA_PAW_PBE PBE52

After that add path to .pmgrc.yaml

pmg config --add PMG_VASP_PSP_DIR PBE52

Provide "vasp.in" file to control the INCAR input file.

Within the "vasp.in" file, keywords that require replacement are listed alongside their corresponding values. 
Additionally, there are keywords without associated values, and they are included to indicate the removal of 
those specific keywords from the INCAR file.
------------------------------------------------------------------------------------------------------------------------------
          After downloading files
------------------------------------------------------------------------------------------------------------------------------

With "mainprogram.py download", a tracking file "mpid.in" file is created with materials IDs and compound name.
Use mpid.in file for further process. 

For QE, put a pseudopotential (PP) folder naming "pp" in the working directory (work_dir). Name PP files in the element.upf format 

Do "mainprogram.py inputinfo" to learn about input file formats and generation of batch scripts to submit calculations in clusters.

F. Finally use generate_submission_file.sh scripts to generate batch script files required for various calculations
   type "generate_submission_file.sh help"
*********************************************************
**********************************************************
******************************"""

print(MSG)
