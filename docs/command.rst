.. _command:

-------------------------------
Commands (mainprogram process)
-------------------------------

Get the basic information about the code

.. code-block:: bash

    mainprogram basicinfo

It provides information about these subprocesses.

.. code-block:: bash

    - process = jobscript generates the job scripts for the calculations
    
    - process = search, search for data in materials project database
    
    - process = download, download QE and VASP input files
    
    - process = oqmd-search, search for data in oqmd database
    
    - process = oqmd-download, download QE and VASP input files
    
    - process = aflow-search, search for data in aflow database

    - process = aflow-download, download QE and VASP input files

    - process = data-combine, combining and eliminating duplicate inputs from different database
    
    - process = process-info, information about various QE+VASP calculations,
    
    - process = epw-info for EPW calculations 
    
    - process = wt-info for wanniertools calculations 
    
    - process = vasp-info for vasp+phonon calculations 
    
    - process = elastic-input, to create vasp input files with deformations
    
    - process = compute-elastic, to compute elastic properties
    
    - process = magenum, to create vasp input files for different magnetic state
    
    - process = fermisurface to plot fermi surface from vasprun.xml
    
    - process = charge-input for creating input files for systems with non-zero net charge
    
    - process = pressure-input for creating input files for different pressure

    - process = change_k for updating kpoint mesh according "kpoint.in" file.

    - process = magmom_extract for exacting final magnetic moment for magnetic calculations, available for VASP.
      Turn on LORBIT to tag to print final magnetization.

    - process = history to print latest 10 mainprogram command executed
      First execute 'history -a' in command line before process = history.
      For MacOs, replace it by ~/.zsh_history in the mainprogram file


Let's look at the "process-info" command,

.. code-block:: bash

    mainprogram process-info

.. code-block:: bash

    - process = e0, to extract the total energies per atom and store in econv_vasp.csv file (QE+VASP)
    
    - process = 1, relax-scan. This will relax the structure for the first time (QE+VASP).
      Now a R{id}-{name} and R{id}-{name}/relax folders are created
    
    - process = 2 updates the input file with new structure (QE+VASP)
      For process = 2, further-relax-input.
      process = 3 resubmit the relaxation with updated input files
      For process = 3, further-relax-scan 
      Repeat process = 2 and 3 for more relaxation (QE+VASP)
    
    - process = 4, create-inputs, (QE).
      Default: qmesh = kmesh/2 along each direction for el-ph calculations
      provide qpoint.in file to provide qpoint mesh for phonon calculation
      All the necessary inputs are created
      inside folders scf_dir,matdyn_dir,elph_dir,q2r_dir,kpath
    
    - process = 5, fine-scan, This perform scf calculations with fine k grid
      Now a R{id}-{name}/calc folder is created (QE)
    
    - process = 6, coarse-scan, performs scf calculations with a coarse k grid (QE)
    
    - process = 7, ph-scan. Performs ELECTRO-PHONON coupling (EPC) calculations (QE)

    - process = checkph, to check the status of the EPC calculations.
    
    - process = 8, q2r-scan (QE)
    
    - process = 9, matdyn-scan (QE)
    
    - process = 10, matdyn-dos-scan. Phonon DOS calculation (QE)
    
    - process = 11, lambda-scan (QE)
    
    - process = 12, phonband-scan. Processing phonon dos (QE)
    
    - process = 13, bandscf-scan (QE+VASP)
      R{id}-{name}/bands folder is created
      for electronic bandstructure and density of states calculations
    
    - process = 14, band-scan. NonSCF band Structure calculation (QE)
    
    - process = 15, bandp-scan. Processing Bandstructure data (QE+VASP)
    
    - process = 16, dos-scan. eDOS calculations (QE+VASP)
    
    - process = 17, dosp-scan. Processing totalDOS (QE)
    
    - process = 18, pdos-scan Processing partial DOS (QE)
    
    - process = 19, plot-scan (QE+VASP)
    
    - process = 20, clean-scan, Removing wavefunctions and bulky folders (QE)
    
    - process = 21, extract-scan, Extracting EPC results and store in result.csv file (QE)
    
    - process = 22, Extracting total energy (QE+VASP) of convergence tests
      for different plane wave cutoff and kpoint mesh.
      run after process=convtest.
    
    - process = 23, dynmat-scan, Obtain atomic displacement files (QE)
      for vibrational mode at Gamma point
    
    - process = 24, distortion-relax-scan, relaxing distorted structure (QE)
    
    - process = 25, distortion-energy-scan,
      collecting distorted structure relaxation results (QE)
    
    - process = 26, pressure-relax-scan, SCF calculations for different pressure (QE+VASP).
      Use 'pressure.in' file with v1 pressure1, v2 pressure2, .... in different line
    
    - process = 27, pressure-ph-scan, phonon calculation for different pressure (QE).
      Create input file ph-{id}-{name}.in with 'mainprogram epw1'.
      ph-q.in file is provided for phonon calculation at particular q point,
      otherwise, provide qpoint.in file for direct generic phonon calculation.
      File ph-q.in file has nq1 nq2 nq3 and metal info on different line.
      if T or t are used, calculation is performed for metal.
    
    - process = 28, delete pressure folder (QE)
    
    - process = 29, element substitution. check 'site_subs.py h' (QE+VASP)
    
    - process = convtest, perform convergence tests for Ecut and kpoint mesh (QE+VASP)
    
    - process = compound, obtain details about compounds, such as structural information before and after relaxation, electron count, Fermi level, kinetic energy cutoff, k-point mesh, etc. (QE+VASP).

.. code-block:: bash

    mainprogram epw-info

.. code-block:: bash

    - perform relaxation and ground-state calculations with process from 1 - 4
    
    - process = epw1 , preparing input files for scf, non-scf, phonon calculations (QE (all), VASP (scf))
    
    - process = qe-ph , scf and phonon calculations (QE)
    
    - process = epw2 , copy phonon files in save directory (QE)
    
    - process = epw3 , projection calculations for scdm projection (QE)
    
    - process = epw4 , fitting procedure to obtain scdm parameters (QE)
    
    - process = wann-scdm , preparing input files (QE)
      for wannierization using scdm projections
    
    - process = wann-file , preparing input files (QE+VASP)
      for wannierization taking projections from projection.in file
    
    - process = wann-random , preparing input files (QE+VASP)
      for wannierization using random projections
    
    - process = 13,14,15,16,17,18 for bandstructure and DOS calculations
      for analyzing and determining different windows (QE (all) + VASP (check with vasp-info))
    
    - process = epw5 , preparing inputfiles for QE bandstructure (QE)
      calculation using kpoints from wannier calculation (to obtain bands on same k-points)
    
    - process = epw-scdm, epw-file, epw-random , preparing input (QE)
      files for epw calculations (anisotropic Eliashberg-Migdel approximations)
      with different projection schemes 

.. code-block:: bash

    # Note1: Utilize atom- and orbital-resolved density of states to provide initial Wannier projections.

    # Note2: To assess the quality of the Wannier orbitals, compare the Wannier interpolated band structure and Fermi surface with those obtained from KS orbitals.

    # Note3: Consider symmetrization of the Wannier interpolated Hamiltonian, if necessary.


.. code-block:: bash

    mainprogram wt-info

.. code-block:: bash

    - First repeat all the calculations as described in 'mainprogram epw-info' command upto wannierization
    
    - process = wt1, prepare input file wt.in required for initial bulk bandgap calculation
      if not found, it will create a default one
      copy 'wt-{id}-{name}.in' file from 'WT_dir' to
      R{id}-{name}/epw/ folder where wannierization process was done
      Please include slab dimension even in bulk calculation,
      so that it produces 'POSCAR-slab' file
      which is used by ASE package to create 'KPATH_SLAB' for slab system
    
    - process = wt2, prepare input file for other calculations including surfaces
  Edit 'wanniertool_input' key in config.json according to properties of interest

.. code-block:: bash

    mainprogram vasp-info

.. code-block:: bash

    - Strucutural relaxation, substitution, pressure, magnetic orderings
      (isotropic, changing scaling factor of lattice) can be performed
      similary as of QE, but using DFT = vasp (or VASP) in input.in file
    
    - For process = 1 - 3, Structural relaxation similar to QE
    
    - process = 13 , Submit bandstructure calculations
    
    - process = 15 , processing bandstructure
    
    - process = 16 , Submit DOS and pDOS calculations
    
    - process = 19 , plotting bandstructure or DOS/pDOS
      Use 'eband', or  'pdos' in input.in
    
    - Thermodynamic quantities can be calculated using VASP and phonopy (vp-ph)
    
    - process = primtoconv, to change structure into conventional unit cell
      ,useful for vasp+phonopy calculations
    
    - process = vp-pd, computing thermodynamic stability
      using pymatgen with 'econv_vasp.csv' file
    
    - process = phono1, to make supercell and submit scf calculations for different displacement.
      To specify the dimensions of the supercell, utilize the "setting.conf" file. Without this specification, the code will default to creating a 2 x 2 x 2 supercell. Note: Consider starting from primitive cell to create supercell.
    
    - process = phono2, computing force constant
    
    - process = phono3, computing and plotting thermodynamic properties
    
    - process = phono4, computing and plotting phonon band
    
    - process = phono5, printing symmetry analysis
    
    - process = phono-qha,
      Computing temperature and pressure dependent thermal properties
    
    - process = ev-collect, extracting the total energies
      for different isotropic volumes from VASP calculations
      Do 'mainprogram 26' calculation before process = ev-collect
    
    - process = phono1-pressure,
      submit phono1 calculations for different isotropic volumes
    
    - process = phono2-pressure,
      submit phono2 calculations for different isotropic volumes
    
    - process = phono3-pressure,
      submit phono3 calculations for different isotropic volumes
    
    - process = phono4-pressure,
      submit phono4 calculations for different isotropic volumes
    
    - process = eos-bm, equation of state fitting using Birch-Murnaghan fit
    
    - process = eos-vinet, equation of state fitting using vinet fit
    

