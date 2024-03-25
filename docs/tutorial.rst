--------------------------------------
Submission scripts generation
--------------------------------------

**Bash script**:

Generate submission files with "generate_submission_file.sh" scripts.

generate_submission.sh help ==> help with the scripts

Make ready "batch.header" file

.. code-block:: bash

    #!/bin/bash
    
    ## SBATCH commands
    
    #SBATCH ....
    
    #module load ....

    submission here
    
Add 'submission here' in the last line to :ref:`batch <batch-label>`. The bash script will replace this placeholder with the appropriate commands. Execute the following command to generate submission files for Quantum Espresso electron-phonon coupling calculations using mpirun (parallel command) and 48 cores. Make sure to customize the "batch.header" file according to the desired number of cores.

.. code-block:: bash

    generate_submission_file.sh qe-elph mpirun 48 

It generates bunch of submission scripts for electron-phonon coupling (EPC) with QE package.

Other available options are epw-elph, wannier_band, and vasp.     

**command line**:

This requires adjusting appropriate keys and values in :ref:`jobscript <job-label>` dictionary in htepc.json file.

.. code-block:: bash

    mainprogram jobscript


--------------------------------------
  Preparing folder for calculations
--------------------------------------

1.Begin by creating a working directory:

.. code-block:: bash

    mkdir work_dir

    cd work_dir

2.Locate the ``htepc.json`` file from the ``utility/input_files`` directory, referenced as :ref:`htepc.json <json-label>`, which functions as the code's input file.

3.In addition to the ``htepc.json`` file, there exist several input files in the ".in" format, intended for direct processing by bash scripts. For VASP calculations, utilize the :ref:`vasp.in <vasp-label>` file.

4.Execute the 

.. code-block:: bash

    mainprogram search 

command to search for data either in the database or in .cif or .vasp files within the working directory. This action generates ``input.in`` and ``mpid-list.in`` files, excluding ``fromcif`` and ``fromvasp`` mode. Proceed to modify the first and second indices in ``input.in`` to generate input files in accordance with the ``mpid-list.in`` file.


5.Utilize 

.. code-block:: bash

    mainprogram download

to generate :ref:`input files <pwd-label>`. This will generate VASP (or QE) input files, if ``DFT = VASP (or QE)`` is set in :ref:`input.in <inputin-label>`. Now the :ref:`mpid.in <mpid-label>` is created. Always set ``DFT = VASP (or QE)`` in :ref:`input.in <inputin-label>` to proceed with calculations.

For QE, it creates ``scf_dir`` directory and put ``scf-{id}.in`` files inside it. For VASP, it creates the ``R{id}-{name}/relax`` folder and put ``INCAR POSCAR POTCAR KPOINTS`` inside it. For example, for MgB2 where ``mpid = mp-763`` and ``name = 'B2Mg1'``, folders named ``Rmp-763-B2Mg1/relax`` will be created to store the downloaded files for VASP, whereas ``scf-mp-763.in`` will be present inside ``scf_dir``. See :ref:`working folder <pwd-label>`.

In addition to generating input files, this process exclusively updates the ``INCAR`` file under the condition that a :ref:`vasp.in <vasp-label>` file is provided, without generating it from scratch. This update only occurs if an ``INCAR`` file is already present within the ``R{id}-{name}/relax/`` directory.

Note: Before executing ``mainprogram download``, configure the ``POTCARs`` according to the instructions provided by ``pymatgen``.

--------------------------------------
Input generation from structure files
--------------------------------------

There are 2 modes to generate input files from structure files.

- **cif Mode**:

This mode is activated when ``'mode':'fromcif'`` is set in the :ref:`download <download-label>` section.

.. code-block:: json

    "download": {
      "info": {
        "mode": "fromcif",
        "metal": false,
        "FE": false,
        "thermo_stable": false,
        "exclude": ["Lu"],
        "...": "..."
      }
    }

By default, it utilizes Pymatgen to explore CIF files.

If ``"use_cif2cell"`` is set to ``true``, it employs the cif2cell package to convert ``.cif`` files to input files.

.. code-block:: json

    "inp": {
      "start": 1,
      "end": 65,
      "nkpt": 200,
      "evenkpt": false,
      "plot": "phband",
      "calc": "VASP",
      "use_cif2cell": true
    }

 
- **VASP Mode**:

In this mode, the code utilizes structure files in ``.vasp`` format to generate input files.

- ``"mode": "fromvasp"`` is turn on.

- We employ the ``.vasp`` format instead of ``POSCAR`` to prevent interference with other activities that specifically require a ``POSCAR`` file format.


Now, execute the following commands to generate input files.

.. code-block:: bash

    mainprogram search

    mainprogram download

-------------------------------------
Inputs with magnetic ordering
-------------------------------------

For QE, one can generate input file with FM ordering by setting ``magnetic`` flag to ``true``.

.. code-block:: json

    "pwscf_in": {
      "magnetic": true,
      "control": {},
      "systems": {},
      "electrons": {}}

Now, any command used to generate Quantum Espresso (QE) input files will create inputs with "FM" (ferromagnetic) ordering.

.. code-block:: bash

    mainprogram download # aflow-download, oqmd-download, substitutiton, ....

For VASP, follow :ref:`this <magenum-label>`, or use ``vasp.in`` to update ``INCAR`` file.


---------------------------------------
Combining data from different database
---------------------------------------

In this section, we explored techniques for extracting data and generating input files from three distinct databases, subsequently amalgamating them. To facilitate this process, we employed the following ``download`` keyword in the :ref:`htepc.json <json-label>` input configuration. In this context, we will delve into the phase diagram of MgB2 using the convex hull method. To achieve this, we require the ground-state configurations of Mg, B, and Mg-B compounds, which we can extract from databases that offer extensive resources.


.. code-block:: json

  "download": {
    "info": {
      "mode": "chemsys",
      "metal": false,
      "FE": false,
      "thermo_stable": false,
      "exclude": ["Lu"],
      "ntype": [1, 2],
      "elm": ["B"],
      "prop": ["material_id", "formula_pretty", "structure", "formation_energy_per_atom", "band_gap", "energy_above_hull", "total_magnetization", "ordering", "total_magnetization_normalized_formula_units", "num_magnetic_sites", "theoretical", "nsites"],
      "ordering": "NM",
      "nsites": 10,
      "spacegroup": null
    },
    "inp": {
      "start": 1,
      "end": 65,
      "nkpt": 200,
      "evenkpt": false,
      "plot": "phband",
      "calc": "VASP",
      "use_cif2cell": false
    },
    "chemsys": {
      "entries": ["Mg", "B"],
      "size_constraint": 60,
      "ntype_constraint": 3,
      "must_include": ["Mg","B"],
      "form_en": false,
      "metal": false,
      "magnetic": true,
      "spacegroup": null
    },
    "oqmd": {
      "limit": 1000,
      "entries": ["Mg", "B"],
      "size_constraint": 60,
      "ntype_constraint": 3,
      "must_include": [],
      "form_en": true,
      "metal": false,
      "magnetic": true,
      "spacegroup": null,
      "thermo_stable": true,
      "FE": true,
      "prop": ["composition", "spacegroup", "volume", "band_gap", "stability"]
      },
    "aflow": {
        "elm": ["Mg","B"],
        "nelm": 2,
        "nsites": 60,
        "metal": false,
        "FE": false,
        "spacegroup": null,
        "filter": false,
        "limit": 5000,
        "prop": [
            "spacegroup_relax", "Pearson_symbol_relax"
        ]
    }
  },


. **Data from Materials Project**:

- Here, we use ``chemsys`` mode.

- In ``chemsys`` dictionary, we set necessary :ref:`parameters <download-label>`.

- Here, we are preparing ``VASP`` input files.

- Now execute:

.. code-block:: bash

    #To perform search.

    mainprogram search

    # This creates a list of compounds and stored in ``mpid-list.in`` file. Now edit ``input.in`` file to include all the compounds in ``mpid-list.in``. To download, execute:

    mainprogram download

    # Now VASP input files are stored within R{id}-{name} folder and the {id} and {name} are stored in ``mpid.in`` file.


. **Data from OQMD database**:

Similarly, once adjusting parameters in ``oqmd`` dictionary, we execute following two commands:

.. code-block:: bash

    mainprogram oqmd-search

    # This creates a list of compounds and stored in ``mpid-list.in`` file. Here, instead of editing input.in, here we directly change ``start`` and ``end`` keyword in ``inp`` dictionary. To download, execute:

    mainprogram oqmd-download

    # Likewise, {id} and {name} are appended to the mpid.in file, with inputs generated within the R{id}-{name} directory.


. **Data from AFLOW database**:

Now, we adjust ``aflow`` dictionary as in example. Similarly, we execute following two commands:

.. code-block:: bash

    mainprogram aflow-search

    mainprogram aflow-download

This will create inputs for ``Mg-B`` binaries. Subsequently, we can utilize ``'elm': ["Mg"]`` or ``'elm': ["B"]`` with ``'nelm': 1`` to explore and download elemental solids.

Each process now generates input files within the ``R{id}-{name}`` directory and updates the ``mpid.in`` file accordingly. The ``mpid.in`` file after executing all of the above process looks like as:

.. code-block:: bash

    v1 mp-110 Mg1
    v2 mp-1056351 Mg1
    v3 mp-1094122 Mg9
    v4 mp-1247180 Mg10
    v5 mp-1055956 Mg1
    v6 mp-973364 Mg4
    v7 mp-1056702 Mg1
    v8 mp-153 Mg2
    v9 mp-978275 Mg4B28
    v10 mp-1016262 Mg7B1
    v11 mp-1016250 Mg3B1
    v12 mp-763 Mg1B2
    v13 mp-1222002 Mg2B6
    v14 mp-365 Mg4B16
    v15 mp-1023515 Mg15B1
    v16 mp-632401 B12
    v17 mp-22046 B50
    v18 mp-1202723 B48
    v19 mp-729184 B40
    v20 mp-570316 B48
    v21 mp-1055985 B1
    v22 mp-1193675 B28
    v23 mp-1196985 B48
    v24 mp-1182425 B12
    v25 mp-1104251 B15
    v26 mp-160 B12
    v27 mp-570602 B50
    v28 oqmd-5087 MgB2
    v29 oqmd-598035 B
    v30 oqmd-752496 Mg
    v31 oqmd-752493 Mg
    v32 oqmd-752498 Mg
    v33 oqmd-752503 Mg
    v34 oqmd-752513 Mg
    v35 oqmd-752518 Mg
    v36 aflow-0ff128f86aa02a78 Mg4B28
    v37 aflow-00cc032d65b23a0e Mg1B3
    v38 aflow-06e000e3dfa7f1a6 Mg1B3
    v39 aflow-1557ff34d66462bd Mg1B3
    v40 aflow-240ad1d192e37bd1 Mg2B6
    v41 aflow-30aabb3b48b8fe87 Mg1B3
    v42 aflow-372f21c264152ba1 Mg1B3
    v43 aflow-4c199fa89b1b5b52 Mg1B3
    v44 aflow-02cc1bd8affd16ae Mg2B4
    v45 aflow-0d34f75d8bfe3ae7 Mg2B4
    v46 aflow-0ec6d7c681e7b2ec Mg2B4
    v47 aflow-1df296ebd052995b Mg2B4
    v48 aflow-21ce6b6bac34f308 Mg1B2
    v49 aflow-258fa48850d77a27 Mg4B8
    v50 aflow-28a586560c918d8b Mg2B4
    v51 aflow-2ac2cc83dd458e57 Mg2B4
    v52 aflow-2ca5e61d6888c369 Mg1B2
    v53 aflow-3fb674873248b3f2 Mg2B4
    v54 aflow-4927ade1c3ed0756 Mg2B4
    v55 aflow-43c18edea03be7cb Mg3B5
    v56 aflow-0cc480408f59af4c Mg6B7
    v57 aflow-1465604b7bf98c6f Mg2B2
    v58 aflow-19d6ec62450aaafc Mg1B1
    v59 aflow-21114a3635ea4f5a Mg1B1
    v60 aflow-25b6828350dc1b7c Mg2B2
    v61 aflow-3149e8cc44907844 Mg6B6
    v62 aflow-3e0f6494b0887580 Mg2B2
    v63 aflow-3e7e4bdf16a54268 Mg2B2
    v64 aflow-3f62c1822fd74775 Mg4B4
    v65 aflow-41b3f291ff574622 Mg1B1
    v66 aflow-0522c21699f6160f Mg2B1
    v67 aflow-0dc1bf12a4577363 Mg4B2
    v68 aflow-1478d09a2eadbf28 Mg4B2
    v69 aflow-35ef05a738f98a62 Mg2B1
    v70 aflow-415ebc0814f84af4 Mg4B2
    v71 aflow-0dead30f9e4daf3a Mg3B1
    v72 aflow-0e54cf897161c5a2 Mg3B1
    v73 aflow-137c05f64299cd85 Mg3B1
    v74 aflow-28698c7a4ef12281 Mg3B1
    v75 aflow-13d6d9dea258ceeb Mg4B1
    v76 aflow-1ea75f4a6b073ff2 Mg4B1
    v77 aflow-04993233df287621 Mg5B1
    v78 aflow-429a554a5809b1ab Mg5B1
    v79 aflow-4f396289df6f3900 Mg7B1


Now, we amalgamate these datasets from various databases and generate a new set of data having distinct spacegroup by executing:

.. code-block:: bash

    mainprogram data-combine

This creates a new file ``mpid-new.in`` and input files within ``filtered_inputs`` directory. This command is not only useful for combining data from different databases but will also help filter out duplicate entry of a compound. Please delete all existing inputs within the working directory and transfer input files from the ``filtered_inputs`` folder. Additionally, replace ``mpid.in`` with ``mpid-new.in``.

A new ``mpid-new.in`` has following data:

.. code-block:: bash

    v1 mp-632401 B12
    v2 mp-22046 B50
    v3 mp-1202723 B12
    v4 mp-729184 B40
    v5 mp-570316 B48
    v6 mp-1055985 B1
    v7 mp-1193675 B28
    v8 mp-1196985 B48
    v9 mp-1182425 B12
    v10 mp-1104251 B15
    v11 mp-160 B12
    v12 mp-570602 B50
    v13 mp-978275 Mg4B28
    v14 mp-365 Mg4B16
    v15 mp-1222002 Mg2B6
    v16 aflow-00cc032d65b23a0e Mg1B3
    v17 aflow-06e000e3dfa7f1a6 Mg1B3
    v18 aflow-1557ff34d66462bd Mg1B3
    v19 aflow-240ad1d192e37bd1 Mg2B6
    v20 aflow-30aabb3b48b8fe87 Mg1B3
    v21 aflow-372f21c264152ba1 Mg1B3
    v22 aflow-4c199fa89b1b5b52 Mg1B3
    v23 mp-763 Mg1B2
    v24 aflow-02cc1bd8affd16ae Mg2B4
    v25 aflow-0ec6d7c681e7b2ec Mg2B4
    v26 aflow-1df296ebd052995b Mg2B4
    v27 aflow-21ce6b6bac34f308 Mg1B2
    v28 aflow-258fa48850d77a27 Mg4B8
    v29 aflow-28a586560c918d8b Mg2B4
    v30 aflow-2ac2cc83dd458e57 Mg2B4
    v31 aflow-2ca5e61d6888c369 Mg1B2
    v32 aflow-3fb674873248b3f2 Mg1B2
    v33 aflow-4927ade1c3ed0756 Mg1B2
    v34 aflow-43c18edea03be7cb Mg3B5
    v35 aflow-0cc480408f59af4c Mg6B7
    v36 aflow-1465604b7bf98c6f Mg2B2
    v37 aflow-19d6ec62450aaafc Mg1B1
    v38 aflow-21114a3635ea4f5a Mg1B1
    v39 aflow-25b6828350dc1b7c Mg2B2
    v40 aflow-3149e8cc44907844 Mg2B2
    v41 aflow-3e0f6494b0887580 Mg2B2
    v42 aflow-3e7e4bdf16a54268 Mg2B2
    v43 aflow-3f62c1822fd74775 Mg4B4
    v44 aflow-41b3f291ff574622 Mg1B1
    v45 aflow-0522c21699f6160f Mg2B1
    v46 aflow-0dc1bf12a4577363 Mg4B2
    v47 aflow-1478d09a2eadbf28 Mg2B1
    v48 aflow-35ef05a738f98a62 Mg2B1
    v49 aflow-415ebc0814f84af4 Mg4B2
    v50 mp-1016250 Mg3B1
    v51 aflow-0dead30f9e4daf3a Mg3B1
    v52 aflow-0e54cf897161c5a2 Mg3B1
    v53 aflow-137c05f64299cd85 Mg3B1
    v54 aflow-28698c7a4ef12281 Mg3B1
    v55 aflow-13d6d9dea258ceeb Mg4B1
    v56 aflow-1ea75f4a6b073ff2 Mg4B1
    v57 aflow-04993233df287621 Mg5B1
    v58 aflow-429a554a5809b1ab Mg5B1
    v59 mp-1016262 Mg7B1
    v60 aflow-4f396289df6f3900 Mg7B1
    v61 mp-1023515 Mg15B1
    v62 mp-110 Mg1
    v63 mp-1094122 Mg3
    v64 mp-1247180 Mg10
    v65 mp-1055956 Mg1
    v66 mp-973364 Mg4
    v67 mp-1056702 Mg1
    v68 mp-153 Mg2

-------------------------------------
Convergence tests
-------------------------------------

Please check :ref:`conv_test <convtest-label>` dictionary:

.. code-block:: bash

    # Run Convergence Tests
    mainprogram convtest
    
    # Collect Total Energies and Store Results in convergence_result folder.
    mainprogram 22


.. _relax-label:

------------------------------
Structure relaxation
------------------------------

To perform the structure relaxation, execute:

.. code-block:: bash

    mainprogram 1

This command creates the R{id}-{name}/relax folder and conducts scf relaxation. For MgB2, where mpid = mp-763 and compound = 'B2Mg1', the Rmp-763-B2Mg1/relax folder is generated.

To update the input files with the relaxed structure, use:

.. code-block:: bash

    mainprogram 2

For subsequent scf relaxations without folder creation, execute:

.. code-block:: bash

    mainprogram 3

Repeat processes 2 and 3 multiple times until the system is fully relaxed.

.. _EPC:

--------------------------------------
 Electron-phonon coupling from QE
--------------------------------------

Perform :ref:`structure relaxation <relax-label>`.

For extracting the structure and generating necessary input files for various (scf,bands,dos,phonon,epw,wannier,wanniertool,pdos,electron-phonon) calculations with Quantum Espresso (QE), use:

.. code-block:: bash

    mainprogram 4

``This step is crucial, don't forget to execute after structure relaxation.``

To conduct electron-phonon calculations, follow these steps:

.. code-block:: bash

    # For SCF calculations using a fine k-grid (twice that of coarse grid) for interpolating EPC quantities, create the "calc" folder inside "R{id}-{name}/".

    mainprogram 5

    # Now SCF calculations using a coarse k-grid, in which EPC calculations is performed.

    mainprogram 6


Here, the el-ph calculations (EPC) is started by executing following command:

.. code-block:: bash

    mainprogram 7

The behavior of this command is described below.

    **(a)** For fresh calculations, it run EPC calculations with `alpha_mix <https://www.quantum-espresso.org/Doc/INPUT_PH.html#idm36>`_ as ``alpha_mix = 0.7``.

    **(b)** If the calculations do not convergence, consider increasing `niter_ph <https://www.quantum-espresso.org/Doc/INPUT_PH.html#idm30>`_ value.

    **(c)** If the calculations still fail to converge, they are resubmitted with ``alpha_mix = 0.3``.

    **(d)** If the calculations fail again with lower ``alpha_mix``, then the script adjusts `nmix_ph <https://www.quantum-espresso.org/Doc/INPUT_PH.html#idm39>`_ as ``alpha_mix = 0.3, nmix_ph = 8``.

Check the status of the EPC calculations with

.. code-block:: bash

    mainprogram checkph

After a converged EPC calculations, postprocessing is performed:

.. code-block:: bash

    mainprogram 8-12

Perform plotting with:

.. code-block:: bash

    mainprogram 19

Options include gammaband in input.in, as well as eband (electronic band), phband (phonon band), wann_band (wannier-interpolated band), pdos (DOS, pDOS), and phonproj (atom projected phonon dispersion).

To extract the results, creating the result.csv file and storing relaxed structures in the "cif" folder, use:

.. code-block:: bash

    mainprogram 21

Finally, clean heavy files and copy to a folder named "completed" with:

.. code-block:: bash

    mainprogram 20

With the ``substitution``, ``pressure``, and ``charge calculations`` modes, we can progress towards identifying those near the convex hull and exploring phonon-mediated superconductivity.

--------------------------------------
 Atom-projected phonon dispersion
--------------------------------------

1.Make sure .eig file is present inside R{id}-{name}/calc folder


2.Use ``phonproj`` keyword in :ref:`input.in <inputin-label>` file in the plot section, and perform ``mainprogram 19``.
   This step creates following files inside R{id}-{name}/calc/ folder.

   phonon-name.proj.gp ==> Has atomic projection for one atoms followed by others. name would be 'MgB2'
   
   phonon-Mg.proj ==> separate file for Mg projections
   
   phonon-B.proj ==> separate file for B projection

3.Check ``plot-proj-{id}-{name}.pdf`` insides plots directory. Here is an example of Y2C3: Y (red) and C (blue).

.. image:: _static/atom_proj.jpg
   :align: center
   :width: 500px
   :height: 300px


.. _band-dos-label:

--------------------------------------
 Band structure and DOS calculation
--------------------------------------

Perform :ref:`structure relaxation <relax-label>`.

For QE, after structure relaxation, execute:

.. code-block:: bash

    mainprogram 4

.. code-block:: bash

    # Perform Bandstructure Calculation (QE). R{id}-{name}/bands folder created.
    mainprogram 13-15
    
    # Perform Density of States (DOS) Calculation (QE). R{id}-{name}/dos folder created.
    mainprogram 16-17
    
    # Perform Partial Density of States (PDOS) Calculation (QE)
    mainprogram 20
    
    # Perform Bandstructure Calculation (VASP)
    mainprogram 13 15
    
    # Perform Density of States (DOS) Calculation (VASP)
    mainprogram 16
    
    # Perform Partial Density of States (PDOS) Calculation (VASP)
    mainprogram 16
    
    # Plotting (QE and VASP)
    mainprogram 19 # Check band_stat.csv inside R{id}-{name}/bands/, which store minimum and maximum eigenvalues of different bands, useful to locate energy windows for wannierization process.

Here, each process needs to execute one at a time. For example, ``13-15`` means executing ``13``, ``14``, and ``15`` in serial mode, while ``13 15`` means executing ``13`` and ``15`` individually.


----------------------------------------
Applying distortion following eigenmode
----------------------------------------

This functionality is only available for Quantum ESPRESSO (QE).

Ensure that the ``{name}.dyn`` file exists inside the R{id}-{name}/calc/ folder, where {id} and {name} correspond to identifiers and names found in the ``mpid-list.in`` file.

.. code-block:: bash

    # Execute dynmat.x and create dynmat.axsf file with eigenmodes.

    mainprogram 23

    # Obtain atomic displacement files for a phonon mode. Relax distorted structures. Perform relaxation of 3N modes for systems with N ions. It creates ``R{i}`` folders, ``i = 1 to N`` inside R{id}-{name}/

    mainprogram 24

    #Collect results "Energy-mode.csv" from the relaxation of distorted structures, can be found in R{id}-{name}/

    mainprogram 25

To obtain unique relaxed energies, follow these steps:

1. Copy ``extract_single_distort`` and ``distort-extract.py`` from the ``utility/distortion`` folder.

2. Execute ``./extract_single_distort start end mpid-list.in`` in your terminal. Replace ``start`` and ``end`` with appropriate indices, and ``mpid-list.in`` with the relevant file containing compound information.

3. This process extracts unique ground-state energies for any compound and stores the results in the ``distorted-energy.csv`` file.    

.. _pressure-label:

--------------------------------------
  Pressure calculation
--------------------------------------

Perform :ref:`structure relaxation <relax-label>`.


Use the :ref:`pressure.in <pressure-label>` file with the following format:

.. code-block:: bash

    all

    v1 50

    v2 100

    v3 150

- The first line represents the ``cell_dofree`` parameter. ``all`` indicates that all angles and axes are moved.

- The subsequent lines, e.g., ``v1 50``, denote the indices (v1, v2, etc.) and the corresponding pressure values in kbar.

For Quantum ESPRESSO (QE):

- Execute ``mainprogram pressure-input`` to create separate files for pressure inside ``scf_dir/scf-mpid-pv.in``, where ``pv`` signifies the pressure value substituted. Additionally, ``mpid-pressure.in`` is created.

For VASP:

- Utilize the "pressure.in" file with the format:

.. code-block:: bash

    v1 0.92 

    v2 0.94

    to isotropically scale the lattice parameters.

For phonon calculations with pressure:

- Use the "ph-q.in" file for the Gamma point calculation. Other generic q-points can also be used.
  
.. code-block:: bash

    0 0 0 
    T

  Here, ``T`` denotes a metal; otherwise, it is considered nonmetallic.

- Execute ``mainprogram epw1`` to generate input files for phonon calculations.

Perform the following operations:

.. code-block:: bash

    mainprogram 26 : Perform SCF relaxation.

    mainprogram 27: Perform phonon calculation.

Additionally, you can incorporate ``mpid-pressure.in`` in ``input.in`` and conduct other calculations.

Execute ``mainprogram 28`` to clean pressure folders.

----------------------------------------------------
 Substitution calculations
----------------------------------------------------


To access help, type ``site_subs.py h``.

In the ``htepc.json`` file, ensure the existence of the :ref:`substitute <substitute-label>`  keyword. To execute the substitution, run:

.. code-block:: bash

    mainprogram 29

There are two modes of substitution:

1.Element Replacement Mode:
   
Replace an element in a parent compound with a dictionary of elements and the number of sites as key-value pairs.

For example, for the compound MgB2 requiring substitution for ``B``, the substitution key will be as follows:

.. code-block:: bash

    "mode": 1

    "elm": 'B'

    "sub": {'B': 1, 'C': 1}, {'B': 0, 'C': 2}

This action will generate two additional files named ``MgBC`` and ``MgC2``. The ``mpid`` and the compound name will be added to the ``mpid.in`` file. Additionally, two input files will be created inside the ``scf_dir`` directory as ``scf-mpid-1.in`` and ``scf-mpid-2.in``.

For VASP, ``Rmp-763-1-MgB2`` and ``Rmp-763-2-MgB2`` folders will be created with INCAR, KPOINTS, and POSCAR files. Create POTCAR separately. Once again, the ``mpid`` and the compound name will be added to the ``mpid.in`` file. Now, use the name of the ``mpid.in`` file in the ``input.in`` file to perform other calculations.

To use the functionality, ensure the presence of the ``bsym`` package.

2.Dictionary Replacement Mode:
   
Utilize a dictionary in which all keys are replaced by their corresponding value pairs.

.. code-block:: bash

    "mode": 2

    "new_sub": {'B': 'C'}


--------------
Fermi Surface
--------------

Perform :ref:`structure relaxation <relax-label>`.

First, generate a :ref:`job script <job-label>` ``run-ifermi.sh``. We need :ref:`ifermi.json <ifermi>` file. Once, we edit based on our need,
execute:

.. code-block:: bash

    mainprogram jobscript

For Fermi surface calculations using IFERMI:

1. Install the IFERMI package from [https://fermisurfaces.github.io/IFermi/introduction.html#installation](https://fermisurfaces.github.io/IFermi/introduction.html#installation).

2. Once the ``vasprun.xml`` file is created inside ``R{id}-{name}/relax/``, execute 

.. code-block:: bash

    mainprogram fermisurface

to generate Fermi surfaces in different formats using the IFERMI package.


---------------------------------------
Thermodynamic stability (Convex Hull)
---------------------------------------

A.Prepare input files

Suppose we are computing phase diagram of ``MgB2``, then we need to download all the elemental solids, binary solids corresponding to composition ``Mg-B``. We do that by switching on the ``chemsys`` mode in :ref:`download <download-label>` keyword in :ref:`htepc.json <json-label>`.

.. code-block:: bash

    "download": {
    "info": {
      "mode": "element", ==> change this to "chemsys"
      "metal": false,
      "FE": false,
      "thermo_stable": false,
      ......
      ......}

Now edit following portion of the ``htepc.json`` file:

.. code-block:: bash

    "chemsys": {
      "entries": ["Mg", "B"],
      "size_constraint": 60, ==> Optimal size of the compounds.
      "ntype_constraint": 3,  ==> Compounds containing fewer than 3 different species, specifically 2.
      "must_include": ["Mg", "B"], ==> If "B" is not included, only compounds containing "Mg" and the binary compound "Mg-B" are extracted since "Mg" must be included. 
      "form_en": false,
      "metal": false,
      "magnetic": true,
      "spacegroup": null},

By switching all other keyword to ``false``, we don't apply any filter on these properties. Now, we will search these
queries on Materials Project database by executing:

.. code-block:: bash

    mainprogram search

Now, we will adjust ``start``, ``end``, and ``mpid-list.in`` in :ref:`input.in <inputin-label>`. We now execute ``download`` command to generate input files. If you want to use ``VASP``, then set ``DFT = VASP``.

.. code-block:: bash

    mainprogram download

Now, the folders appeared in ``R{id}-{name}`` format, with ``id`` and ``name`` recorded in ``mpid.in`` file. Now change ``mpid-list.in`` to ``mpid.in`` in ``input.in`` file and perform structure relaxation by executing:

Perform :ref:`structure relaxation <relax-label>`.

Repeat the processes ``process = 2 and 3`` several times to ensure full relaxation. If the relaxation completes in 1 ionic step, using ``mainprogram 2`` will change the ``NSW`` keyword inside ``INCAR`` to ``0``. Finally, execute ``mainprogram 3`` for VASP to obtain accurate total energies. In QE, the code always performs one more electronic self-consistent field (SCF) after the convergence of each ionic relaxation, therefore, there is no need to worry about it.

Execute:

.. code-block:: bash

    mainprogram e0

to collect the total energy per atom. It will create ``econv_vasp.csv`` file. Finally, execute:

.. code-block:: bash

    mainprogram vp-pd

to compute phase diagram. Data are stored in ``convexhull.csv`` and ``convexhull.pdf`` plot is created. 

.. _magenum-label:

-----------------------
Magnetic enumeration
-----------------------

Perform :ref:`structure relaxation <relax-label>`.

There are 2 types of magnetic enumeration process with ``VASP``, ``ordering`` and ``magnetic anisotropy``. This requires update in "magmom":

- **ordering**:

.. code-block:: bash

    "magmom": {
        "magmom": {
          "Mn": 5,
          "Cr": 5,
          "Fe": 5,
          "B":0},
        "type":"ordering",
        "saxis":[[0,0,1],[1,0,0],[1,1,0],[1,1,1]],
        "order": ["ferromagnetic", "antiferromagnetic", "ferrimagnetic_by_motif"]}


- **Magnetic Anisotropy**:

.. code-block:: bash

    "magmom": {
        "magmom": {
          "Mn": 5,
          "Cr": 5,
          "Fe": 5,
          "B":0},
        "type":"anisotropy",
        "saxis":[[0,0,1],[1,0,0],[1,1,0],[1,1,1]],
        "order": ["ferromagnetic", "antiferromagnetic", "ferrimagnetic_by_motif"]}

To obtain input files for magnetic anisotropy calculations (MAEs), we can achieve in two step. First, add ``LSORBIT .TRUE.`` key in :ref:`vasp.in <vasp-label>` file (before lines with single column), and update the ``INCAR`` by executing:

.. code-block:: bash

    mainprogram download

This will write MAGMOM in ``mx my mz`` format. Now execute:

.. code-block:: bash

    mainprogram magenum 

This creates input files with ``mpid-magnetic.in`` file. Update ``mpid-magnetic.in`` in ``input.in`` execute relaxation command:

.. code-block:: bash

    mainprogram 1

Similarly repeat ``process 2 and 3`` for complete relaxation.

---------------------------
Magnetic force theorem
---------------------------

----------------------------
Computing elastic constants
----------------------------

Perform :ref:`structure relaxation <relax-label>`.

To execute strain calculations, modify the ``strain`` keyword in the :ref:`htepc.json <json-label>` file. Then, using the ``input.in`` and ``mpid.in`` files, execute the following commands:

.. code-block:: bash

    mainprogram elastic-input

This command applies strain to the conventional unit cell, generates deformed structures, and submits calculations.

Once the calculations are completed, execute the following command:

.. code-block:: bash

    mainprogram compute-elastic

This command computes the elastic constants and stores the results in the ``elastic.csv`` file. Make sure, your results are converged with respect to ``plane-wave energy cutoff`` and ``k-point mesh``.


----------------------------------
Printing compound information
----------------------------------

Now, basic information about systems under calculations can be printed by executing:

.. code-block:: bash

    mainprogram compound

It prints data as follows:

.. code-block:: bash

    Printing info about compounds. Do 'mainprogram compound > compound.txt' to save in 'compound.txt' file
    
    *               Compound: Mg1B2                         
    
    Looking for scf-mp-763.in file. Structural parameters before relaxation
    
    ************ Printing structural parameters *****************
    
             Cell par: (3.0627622617599624, 3.0627622617599624, 3.52087, 90.0, 90.0, 120.00000000000001)
    
             Cell volume: 28.6027 $\AA^3$
    
             Spacegroup: ('P6/mmm', 191)              
    
    
    Looking for scf-relax-mp-763-Mg1B2.in file. Structural parameters after relaxation
    
    ************ Printing structural parameters *****************
    
             Cell par: (3.0691304614624833, 3.0691304614624833, 3.515658311, 90.0, 90.0, 119.99999998841307)
    
             Cell volume: 28.6793 $\AA^3$
    
             Spacegroup: ('P6/mmm', 191)              
    
             Valence Electrons: 16.00                    
    
             Fermi Energy: 9.1373 eV                    
    
             KEcutoff: 35, Ry                         
    
             K-mesh info: (automatic)  16 16 16 0 0 0                          
    
             q-mesh info: [nq1=4,nq2=4,nq3=3]                                  
                                                        
    ******************************************************
    
    all done


-----------------------------------
Wannier interpolated bandstructure
-----------------------------------

Perform :ref:`structure relaxation <relax-label>`.

Additionally, compute :ref:`DFT band structures and atom and orbital projected density of states (DOS) <band-dos-label>`. The DOS information can be valuable for selecting initial projections, while band structures stored in the band_stat.csv file can aid in choosing energy windows.


-------------------------
Phonon calculations
-------------------------

- **DFPT**: 

Perform :ref:`structure relaxation <relax-label>`.

Prepare, QE input files for scf calculations to obtain the charge density by executing:

.. code-block:: bash

    mainprogram 4

One can extract phonons from :ref:`electron-phonon coupling (EPC) calculations <EPC>`.

Similarly, one can perform phonon calculations by executing:

.. code-block:: bash

    mainprogram epw1

    mainprogram qe-ph

Post-processing processes to obtain phonon dispersion plots are:

.. code-block:: bash

    mainprogram 8

    mainprogram 9

    mainprogram 12

    mainprogram 19

Please, refer to :ref:`command line info <command>` for description of these commands.


- **Supercell method**:

For this method, we utilize Phonopy with either Quantum ESPRESSO (QE) or VASP. Please install `Phonopy <https://phonopy.github.io/phonopy/install.html>`_ and ensure that the ``phonopy`` command is accessible.


.. code-block:: bash

    #Prepare input files, create supercells, generate structures with displaced ions, and submit SCF calculations for each displacement.
    mainprogram phono1

    #Compute force constants
    mainprogram phono2

    #Generate band.conf file and plot phonon dispersion.
    mainprogram phono4

--------------------------
Equation of states
--------------------------

To collect energy-volume data for different pressures and perform relaxation, use the ``ev-collect`` command. This command generates an ``e-v.dat`` file in each directory corresponding to a specific pressure. Follow the steps below:

1. First, create input files for different pressures and perform relaxation (:ref:`repeat this <pressure-label>`).

2. After relaxation, execute the main program ``ev-collect`` to extract energy-volume data. This program automatically generates an ``e-v.dat`` file in each directory corresponding to a specific pressure.

.. code-block:: bash

    mainprogram ev-collect

The ``eos-bm`` command is used to obtain the equation of state (EOS) using the Birch-Murnaghan EOS. To use this command, follow these steps:

**Generate Energy-Volume Data:** Before using ``eos-bm``, ensure you have collected energy-volume data using the ``ev-collect`` command.

Execute the main program ``eos-bm`` to analyze the energy-volume data and obtain the equation of state using the Birch-Murnaghan EOS.

.. code-block:: bash

    mainprogram eos-bm

To plot various volume-energy, pressure-volume, and pressure-enthalpy curves, follow these steps:

1. **Copy Plotting Script:** Copy the script ``birch_murnaghan_enthalpy.py`` from the utility ``useful_scripts`` directory.

2. **Execute the Script:** Execute the copied script using the command ``python birch_murnaghan_enthalpy.py``. Use the "help" option to get started and understand the plotting options available.

3. **Rename ``e-v.dat`` Files:** Ensure that each ``e-v.dat`` file from different directories is renamed to follow the format ``e-v-1.dat``, ``e-v-2.dat``, and so on. This ensures that the script can process multiple energy-volume datasets.


--------------------------
Charged calculations
--------------------------

First prepare :ref:`charge.in <charge-input>`.

Execute:

.. code-block:: bash

    mainprogram charge-input

To create input files with different net charges:

For QE:
- Input files named ``scf-{mpid}-{icharge}.in`` are created within the ``scf_dir`` folder, where ``{icharge}`` takes values from 1 onwards.

For VASP:
- Folders named ``R{mpid}-{icharge}`` are created, each containing the necessary input files for a specific net charge.

In both cases, a file named ``mpid-charge.in`` is generated to list the material ID and compound name of these input files.


------------------------------
Checking running calculations
------------------------------

.. code-block:: bash

    check_calc <your_queue_command> <your_account_id>

    #Here, queue command could be ``squeue``.

--------------------------------
Printing history of mainprogram
--------------------------------

To view the ten most recent occurrences of the "mainprogram" command within the current session, execute the following two commands.

.. code-block:: bash

    history -a
   
    mainprogram history

--------------------------------
Aborting jobs with job_id
--------------------------------

.. code-block:: bash

    cancel_job <job_id_start> <number_of_jobs> <Your_job_cancellation_command>
