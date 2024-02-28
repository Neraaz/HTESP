This section has descriptions of files and parameters, required for the package.

-------------------
htepc.json
-------------------

This json file serves as the main input file of the package. It has dictionary with following keys,

 
       'mpi_key' ==> materials project API key 

       'download' ==> information required for downloading and preparing inputs

       'magmom'   ==> magnetic moment for VASP calculations

       'pseudo'   ==> pseudopotential file, has 2 subkeys. 'pot' for VASP, and 'PSEUDO' for QE.

       'substitute' ==> parameters for preparing input files with substitution

       'pwscf_in' ==> input parameters for QE

       'bgw_input' ==> input parameters for QE+BerkeleyGW calculations

       'strain_inp' ==> strain inputs

       'wanniertools_input' ==> input parameters for wanniertools

       'kptden ==> kpoint density
    

--------------------
mpi_key
--------------------

Find your materials project key here, under API key section.

https://next-gen.materialsproject.org/api#api-key

--------------------
download
--------------------

It has a dictionary of the form.

"download": {
    "info": {
      "mode": "element",
      "metal": false,
      "FE": false,
      "thermo_stable": false,
      "exclude": ["Lu"],
      "ntype": [1, 2],
      "elm": ["B"],
      "prop": ["material_id", "formula_pretty", "structure", "formation_energy_per_atom", "band_gap", "energy_above_hull", "total_magnetization", "ordering", "total_magnetization_normalized_formula_units", "num_magnetic_sites", "theoretical", "nsites"],
      "ordering": "NM",
      "nsites": 10,
      "spacegroup": None},
    "inp": {
      "start": 1,
      "end": 2,
      "nkpt": 200,
      "evenkpt": false,
      "plot": "phband",
      "calc": "QE"},
    "chemsys": {
      "entries": ["Fe", "Pd", "I"],
      "size_constraint": 60,
      "ntype_constraint": 4,
      "must_include": ["Fe","Pd","I"],
      "form_en": false,
      "metal": false,
      "magnetic": true,
      "spacegroup": null},
    "oqmd": {
      "limit": 400,
      "entries": ["Mg", "B"],
      "size_constraint": 60,
      "ntype_constraint": 4,
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
      "elm": ["B"],
      "nelm": 5,
      "nsites": 10,
      "metal": true,
      "FE": true,
      "spacegroup": null,
      "limit": 5000,
      "prop": ["Bravais_lattice_orig", "Bravais_lattice_relax", "composition", "density", "dft_type", "eentropy_cell", "enthalpy_cell", "enthalpy_atom","enthalpy_formation_cell", "enthalpy_formation_atom", "entropic_temperature","kpoints", "volume_cell", "volume_atom"]}}

Here, the keyword 'info' is to control the preparing input files.

    a. There are 4 'modes' for preparing input files.
       
       element ==> Extract data from MP database using element based search method

       chemsys ==> Extract data from MP database using chemsys (eg. Mg-O) based search method

       fromcif ==> Convert .cif files into QE and VASP input files

       fromvasp ==> Convert poscar files in .vasp format into VASP input files
    
    b. metal:True selects zero bandgap systems.

    c. FE;True selects compounds with negative formation energy.

    d. thermo_stable:True selects the compounds on convex hull.

    e. exclude list of the elements to exclude from the search.

    f.  ntype number of different types of elements in compounds.

    g.  elm is the element to search in compounds. Can use upto 2 elements. For example, ['B','C'] to search boron and/or carbon containing compounds.

    h.  'prop' is list of the properties to extract.

    i.  'ordering' magnetic ordering to search. "NM" for nonmagnetic, "FM" for ferro, and "AFM" for antiferromagnetic, and so on.

    j.  'nsites' total number of ions in the compound.

    k.  'spacegroup', spacegroup symbol to chose. If None, then the code doesn't select based on spacegroup.

The keyword 'inp' controls the calculations, and also writes input.in file. It's keywords are

    a. "start" 

    b. end

    c. nkpt

    d. evenkpt

    e. plot

    f. calc

--------------------
magmom
--------------------

--------------------
pseudo
--------------------

--------------------
substitute
--------------------

--------------------
pwscf_in
--------------------

--------------------
strain
--------------------

--------------------
wanniertools_input
--------------------

--------------------
kptden
--------------------

--------------------
bgw_input
--------------------

--------------------
input.in
--------------------

  1    #Starting index

  2    #Ending index (not included)

  400 0  #number of kpoints for band calculations and cutoff of Brillouin-zone high-symmetry path, counting from last point.

  mpid-list.in  #filename with materials id and compound name to work with.

  phband   #type of plot 

  DFT = QE   #DFT calculation

---------------------
mpid-list.in
---------------------
 
 v1 mp-763 B2Mg1

 v2 mp-944 B2Al1

First column is just the identifier, second column is material id, and the third column is compound name.

---------------------
mpid.in
---------------------

Similar file as of mpid-list.in. But only written after executing "mainprogram download" command. This simply
checks duplication and only updated for new materials id, and acts as the tracking file.



