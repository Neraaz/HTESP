Do "mainprogram basicinfo" for various information. 

Key processes are
    
    a. install ==> Information about required packages and installation
    b. extract ==> Information about extracting input files from materials project database
    c. inputinfo ==> Some info about input files
    d. gs-info ==> Information about ground-state QE and VASP calculations
    e. elph-info ==> electron-phonon coupling calculations with QE
    f. epw-info ==> Informationa about EPW calculations
    g. vasp-info ==> Information about various vasp calculations
    h. wt-info ==> Information about wanniertools input generation
    i. compound ==> Extract various structural and energetic information from calculations, works only for QE
    j. convtest ==> convergence tests for kinetic energy cutoff
    k. change_k ==> change k-point sampling according to kpoint.in file. Default: 2 times original K-point mesh
       kpoint.in has the same structure of qpoint.in (Do mainprogram inputinfo for more)
    l. pressure-input ==> insert pressure in QE input files and create input files with different pressure and new mpid-pressure.in file to track.
    m. fermisurface ==> plot Fermi surface using IFERMI package and vasprun.xml file
    n. search ==> search .cif files in working directory or data in materials project (MP)
    o. download ==> convert structure (.cif) files to inputs or download files from MP

--------------------------------------
  Preparing folder for calculations
--------------------------------------

    1. Create a folder to work on: mkdir work_dir ==> cd work_dir

    2. Type "mainprogram extract" for more details.

    3. Find htepc.json file from utility/parameters directory, which serves as the input file for the code.
       json file consists of dictionaries with following keywords,

    4. Apart from htepc.json file, there are several input files in ".in" format,
        to be directly processed by bash scripts. some of these files are either produced internally (IP,internally produced) or
        some of them should be provided by the users (UP, user provided). Here are some of the major files,

         input.in (IP),
         mpid-list.in (IP),
         mpid.in (IP),
         mpid-substitute.in (IP),
         mpid-magnetic.in (IP),
         mpid-pressure.in (IP),
         kpoint.in (UP),
         vasp.in (UP),
         doping.in (UP),
         pressure.in (UP),
         filedos.in (UP or IP),
         vasp-band.in (UP or IP),
         vasp-phonopy.in (UP or IP),
         phonproj.in (IP),
         proj-wt.in (UP),
         qpoint.in (UP),
         ph-q.in (UP),
       
    5. "mainprogram search" search the data either in database or .cif or .vasp files in working directory,
       input.in and mpid-list.in files are created. Now, change first and second indices in input.in to create
       input files according to mpid-list.in file.

    6. "mainprogram download" to download input files.

--------------------------------------
 Setting QE calculations
--------------------------------------

Once we have "input.in" and "mpid-list.in" files, then "mainprogram download" with "DFT = QE" will download QE input files inside "scf_dir" folder in "scf-mpid.in" format, where "mpid" may be materials id of the system from Materials Project Database. Once downloaded, do "mv mpid.in mpid-list.in", before starting further process. Prefix inside QE 'scf.in' inputs has 'AmBn' form where the compound has m number of 'A' ions and n number of 'B' ions. Rmpid-compound folder in QE case should be Rmpid-prefix.

Type "mainprogram qe-info" to get more detail about QE calculations

Put pseudopotentials (PP) files inside "pp" folder and named PPs in "element.upf" format. For Hydrogen, "H.upf". "pp" folder should be inside working directory (work_dir). Available PP with this package only includes ultrasoft and/or norm
conserving ones, which is utilized for electron-phonon coupling calculations. Please, check PPs yourself !

--------------------------------------
Submission scripts generation
--------------------------------------

Generate submission files with "generate_submission_file.sh" scripts.

generate_submission.sh help ==> help with the scripts

Make ready "batch.header" file

Type: generate_submission_file.sh qe-elph mpirun 48 to create submission files for QE el-ph coupling calculations

with mpirun (parallel command) and 48 cores. Adjust "batch.header" accordingly for number of cores.

Default parameters in the code are used in following work


Input of QE scf input file can be controlled using "pwscf_in" keyword with htepc.json file. file in working directory, which has dictionary objects for different section as readable by pymatgen package. If it is not provided, it will use default value and create this file. One can edit this file to use desired parameters for the calculations.

--------------------------------------
 Electron-phonon coupling from QE
--------------------------------------

Type: "mainprogram 1" to perform the first relaxation. It creates Rmpid-compound/relax folder and perform scf relaxation. For MgB2, mpid = mp-763 and compound = 'B2Mg1', Rmp-763-B2Mg1/relax folders are created.

"mainprogram 1 (or 2 or 3 )" perform submit/resubmit relaxation. 

"mainprogram 4-12" completes electron-phonon calculations after creating calc folder inside Rmpid-compound/
Process = 7 submit/resubmit the el-ph calculations

"mainprogram 19" performs plotting (For eg. use gammaband in input.in). Other options are eband (electronic band), phband (phonon band), wann_band (wannier-interpolated band), pdos (DOS,pDOS), phonproj (atom projected phonon dispersion)

"mainprogram 21" extract the results (result.csv file is created and relaxed structures are stored in "cif" folder).

"mainprogram 20" clean the heavy files and copy to a folder "completed".

--------------------------------------
 Atom-projected phonon dispersion
--------------------------------------

1. Make sure .eig file is present inside Rmpid-compound/calc folder

2. Use 'phonproj' keyword in input.in file in the plot section, and perform 'mainprogram 19'.
   This step creates following files inside Rmpid-compound/calc/ folder.

phonon-name.proj.gp ==> Has atomic projection for one atoms followed by others. name would be 'MgB2'

phonon-Mg.proj ==> separate file for Mg projections

phonon-B.proj ==> separate file for B projection


--------------------------------------
 Band structure and DOS calculation
--------------------------------------

Use process 13-15 perform the bandstructure, 16-17 perform the density of states (DOS).

20 perform partial DOS (PDOS)

Use eband,pdos in input.in to produce quick bandstructure or DOS/pDOS plot. 

For VASP, use process 13 and 15 for bandstructure, 16 for DOS and pDOS. Finally use mainprogram 19 for plotting.
 Convergence tests for kinetic energy cutoff

--------------------------------------
Some other processes
--------------------------------------

Use process = convtest to run convergence tests for ecut = 30 Ry to 110 Ry with 10 Ry interval. 

Use process = 22 to collect the total energies in Energy-system.csv file

 Distorted ground-states for any mode

Make sure you have "compound.dyn" file inside Rmpid-compound/calc/ folder, mpid and compound can be seen in "mpid-list.in" file.

For process = 23, dynmat-scan, Obtain atomic displacement files for a phonon mode

For process = 24, distortion-relax-scan, relaxing distorted structure. Perform relaxation of 3N modes of system having
N number of ions.

For process = 25, distortion-energy-scan, collecting distorted structure relaxation results

Type "mainprogram 26" ==> perform scf relaxation

Type "mainprogram 27" ==> phonon calculation

Apart from this, one can also use "mpid-pressure.in" in input.in and perform other calculations.

Type "mainprogram 28" to clean pressure folders.

Copy extract_single_distort and distort-extract.py from utility folder. Type "./extract_single_distort start end mpid-list.in"

This will extract unique ground-state energies for any compound and store in "distorted-energy.csv" file.

--------------------------------------
  Pressure calculation
--------------------------------------

Use "pressure.in" file with following format

   all
   
   v1 50
   
   v2 100
   
   v3 150
   
First line represents the cell_dofree parameter. 'all' means all angle and axis are moved.

Here v1, v2,.. run the indices while 50, 100, and 150 represent pressure in kbar.

Type "mainprogram press-qe" to create separate files for pressure inside scf_dir/scf-mpid-pv.in, where pv means pressure value substituted. Also mpid-pressure.in is created.

For vasp, use pressure.in in following format.
v1 0.92
v2 0.94
...
...
Here 0.92, 0.94, ... are scaling factor to isotropically scale the lattice parameters

 For phonon calculation with pressure
 With "ph-q.in" file. For Gamma point. Can use other generic q-points
 
     0 0 0
     T

'T' is used for metal, otherwise nonmetallic. Do "mainprogram epw1" to create input files for phonon calculation.


----------------------------------------------------
 EPW for anisotropic superconductivity calculations
----------------------------------------------------

Compile epw package with "make epw" inside Quantum Espresso folder.

Create submission scripts from generate_submission_file.sh using "epw-elph" and "wannier_band".

Type "mainprogram epw-info" to look for more information about EPW calculations.

1. Before performing EPW calculations, perform relaxation from process = 1 to 3, and prepare other input files with process = 4.

2. perform bandstructure and DOS calculations (process = 13-18) to figure out frozen window, disentanglement window, and projectors (from PDOS) for wannierization.

3. Type "mainprogram epw1" ==> for preparing scf, non-scf, and phonon calculations required for anisotropic superconductivity calculations.

4. "mainprogram epw2" ==> for submitting these calculations. Creates a folder "phonon" inside Rmpid-compound.

5. "mainprogram epw3" ==> copy phonon files in save directory 

6. "mainprogram epw4" ==> projection calculations for SCDM erfc projections for wannierization (https://doi.org/10.1038/s41524-020-0312-y). Creates a folder "epw" inside Rmpid-compound

7. "mainprogram epw5" ==> fitting process for SCDM parameters. 
    Install "lmfit" package: https://lmfit.github.io/lmfit-py/

8. "mainprogram epw6-scdm" ==> prepare input files for wannierization using automatic scdm projection

9. "mainprogram epw6-random" ==> prepare input files for wannierization with random projectors. Change various parameters including projectors, energy windows, number of bands, number of wannier function projectors, etc. Go to Rmpid-compound/epw and submit the submission scripts.

With calculations by default, we obtain wannier interpolated band structure, and fermi surface file.
Use plot_band_wannier.py from utility to plot bandstructure. We can use Fermisurfer to plot Fermi surface. Using these, we can check the quality of our wannier functions.

"python plot_band_wannier.py bandfile fermi_energy output_file"

10. "mainprogram epw8-scdm" ==> prepare input files for EPW calculations with scdm projection. Creates a folder "EPW" inside Rmpid-compound

11. "mainprogram epw8-random" ==> prepare input files for EPW calculations with random projection. Creates a folder "EPW" inside Rmpid-compound. Now update similar parameters for wannierization after reproducing similar bands and Fermi surface as of Kohn-Sham calculations. Update "epw.in" file according to the need of the calculations.

Go to Rmpid-compound/EPW and submit job script.

EPW calculations can't be perform automatically due to parameters choice for wannierization. However, creating input files automatically can save a lot of time.

----------------------------------------------------
 Substitution calculations
----------------------------------------------------

Requires: bsym package

Type "site_subs.py h" to get help 

require "substitute" keyword in htepc.json

Suppose for MgB2, we need substitution for 'B'. Then substitute key will be as follows.

elm = 'B'

sub={'B':1, 'C':1},{'B':0,'C':2}

This will create 2 other files with "MgBC" and 'MgC2". "mpid" and "compound name" will be added "mpid.in" file. Also
2 input files will be created inside scf_dir as scf-mpid-1.in and scf-mpid-2.in

For VASP, Rmp-763-1-MgB2 and Rmp-763-2-MgB2 folders will be created with INCAR, KPOINTS, and POSCAR. Create POTCAR separately. Again "mpid" and "compound name" will be added "mpid.in" file. Now use "mpid.in" file name in input.in to perform other calculations.

----------------------------------------------------
 Setting VASP calculations
----------------------------------------------------

Once we have "input.in" and "mpid-list.in" files, then "mainprogram download" with "DFT = vasp" will download vasp input files (INCAR,POSCAR,POTCAR,KPOINTS) inside Rmpid-compound/relax/ folder, where mpid and compound are the materials id and name respectively.

For MgB2, mpid = mp-763 and compound = 'B2Mg1', Rmp-763-B2Mg1/relax folders are created.

We can also provide vasp.in file to override input parameters of INCAR downloading from Materials Project. A simple format of vasp.in can be found in HTEPC/utility/vasp.in.

Configure POTCARS before performing "mainprogram download" according to pymatgen instructions.

Type: "mainprogram 1-3" to perform the relaxations

 Fermi surface calculations using IFERMI

Install IFERMI package: https://fermisurfaces.github.io/IFermi/introduction.htmlinstallation

Once vasprun.xml is created inside Rmpid-compound/relax/, then type "mainprogram fermisurface" to create fermi surfaces in different formats using IFERMI package. 

----------------------------------------------------
 VASP and Phonopy Calculations
----------------------------------------------------

Install Phonopy package: https://phonopy.github.io/phonopy/

After performing relaxation, one can extract the total energy per atom and compound information using "mainprogram vp-ph1" command.

Usage: "mainprogram arg"

if arg = vp-ph1 ==> extract the total energy per atom in econv_vasp.csv

if arg = vp-pd, computing thermodynamic stability (Phase diagram)
 using pymatgen with 'econv_vasp.csv' file

if arg = vp-ph2 ==> creates 2x2x2 supercell and submit scf relaxations. It requires vasp-phonopy.in file. Look at utility folder for a demo.

if arg = vp-ph3 ==> create force constant

if arg = vp-ph4 ==> thermal properties with 48x48x48 mesh

if arg = vp-ph5 ==> phonon bandstructure

if arg = vp-ph6 ==> print symmetries in symmetry_analysis.in file

if arg = vp-ph-qha,
 Computing temperature and pressure dependent thermal properties

if arg = ev-collect, extracting the total energies
 for different isotropic volumes from VASP calculations

Do 'mainprogram 26' calculation before process = ev-collect

if arg = vp-ph2-pressure,
 submit vp-ph2 calculations for different isotropic volumes

if arg = vp-ph3-pressure,
 submit vp-ph3 calculations for different isotropic volumes

if arg = vp-ph4-pressure,
 submit vp-ph4 calculations for different isotropic volumes

if arg = vp-ph5-pressure,
 submit vp-ph5 calculations for different isotropic volumes

if arg = eos-bm, equation of state fitting using Birch-Murnaghan fit

if arg = eos-vinet, equation of state fitting using vinet fit

Do 'mainprogram vasp-info' for more details

----------------------------------------------------
 Convex Hull using Pymatgen
----------------------------------------------------

After getting "econv_vasp.csv" file, Do 'mainprogram vp-pd' to plot convex hull using pymatgen package.

This also output "convexhull.csv" with last 2 columns are formation energies and energy above hull data extracted from

Materials Projects Database.

 Utility folder contains additional scripts for postprocessing and analysis. README file contains more info.

