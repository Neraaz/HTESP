copy batch.header and config.json file to tutorial folders.

suppose config.json file turns into input_data dictionary.

Always set input_data['download']['inp']['calc'] to "VASP"
to generate VASP input files

POTCAR hasn't been provided. If you have POTCAR files, one can configure using pymatgen
and code automatically create necessary POTCAR based on POSCAR file.

tutorial1: Generate submission file.

tutorial2: Extract data from materials project in element mode

tutorial3: Extract data from materials project in chemsys mode

tutorial4: Extract data from OQMD database

tutorial5: Extract data from AFLOW database

tutorial6: Extract data in magnetic configuration

tutorial7: Combine database from 3 database

tutorial8: Generate input files from .cif files

tutorial9: Perform structure relaxation

tutorial10: Perform convergence tests

tutorial11: Bandstructure and DOS calculations

tutorial12: Preparing input files for different pressures or volumes

tutorial13: Generating input files with substitutions

tutorial14: Computing elastic constants

tutorial15: Thermodynamic phase stability

tutorial16: Automated phonopy+VASP calculations to compute phonon bandstructure

tutorial17: Equation of State

tutorial18: Wannier interpolated bandstructure

tutorial19: Preparing input files for systems with non-zero net charge

tutorial20: Create input files for different magnetic ordering

tutorial21: 3D Fermi Surface using IFermi package.
