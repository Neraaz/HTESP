______________


                        ****    ****   ***********   *********   ********       ********                 
                        |  |    |  |       | |       | |____     | |            |  |  \ \                 
                        |  |____|  |       | |       |  ____|    | |*****       |  |__| |                       
                        |   ____   |       | |       | |_____          | |      |  _____/                       
                        |  |    |  |       | |       |_______|   ******| |      |_ |                       
                        |__|    |__|       |_| *********************************************                                          
                        **********************                                                               
                                      High Throughput Electron-Structure Package                              

                                              Program written by

                                Niraj K Nepal, PhD       &       Lin-Lin Wang, PhD                        
                          Email: nnepal@ameslab.gov                              
                                 tug11655@temple.edu                              


_____________________________________________________________________________________________________________________________________
######################################################################################################################################

   # HTESP: DOCUMENTATION

##                                                                  Utilities:
                                
## 1. Downloading and preparing input files from Materials Project database for QE and VASP calculations.

## 2. Supports Quantum Espresso (QE) 
      a. Ground-state calculations (structure relaxation, band, and dos calculations), convergence tests for kinetic energy cutoff
      b. electron-phonon calculations and superconductivity using isotropic Eliashberg approximation including plotting of spectral function (alpha^2F), phonon dispersion (with or without atomic projections), bandstructure, density of states
      c. Prepares input files for wannier90, EPW (anisotropic superconductivity) calculations, and wanniertools calculations. Need to determine energy windows for wannierization processes
      d. Performing substitution, single phonon mode, distortion calculations for all modes at a single q-point with dynmat.axsf file, system under pressure calculations  
      
## 3. Supports VASP
      a. phonon and thermodynamic calculations with phonopy package and VASP
      b. Ground-state calculations to compute thermodynamic phase diagram (convex hull) using pymatgen
      c. Fermi surface calculations using IFERMI package
      d. Various automated calculations with vasp and pymatgen such as computing elastic properties, magnetic ordering,etc.
      

## Requirements

#### Current package is tested only for Linux Distribution, with Python and Bash languages.
Basic requirements

numpy, scipy, pandas, matplotlib

Pymatgen: https://pymatgen.org/

ASE: https://wiki.fysik.dtu.dk/ase/

mp_api: https://next-gen.materialsproject.org/api

Extra packages

lmfit: https://lmfit.github.io/lmfit-py/. conda install -c conda-forge lmfit. This is needed for SCDM fit to calculate initial projections for wannierization.

IFERMI: https://fermisurfaces.github.io/IFermi/introduction.html#installation. Required for Fermi surface generation.

bsym: https://bsym.readthedocs.io/en/latest/index.html. Required for substitutions of elements in crystal.

Phonopy: https://phonopy.github.io/phonopy/

## Package structure
    Parent folder: HTESP
    Sub folders: docs (Documentation), examples, src, utility (various additional scripts)
    installation script: setup.py
    License file: LICENSE
    Readme file: README.md
    Inside src, there is a "bash" folder that has bash scripts for running calculations

## Installation:
## Conda environment
conda create --name myenv python==3.9.12

source activate myenv

pip install -r requirements.txt

# Also install phonopy in the conda environment
# Download software 
git clone https://github.com/Neraaz/HTESP.git

cd HTESP
# Install HTESP package

pip install .

# check executable

which mainprogram

# Do "mainprogram basicinfo" to begin.
#Install phonopy to perform phonopy calculations

#Alternatively,

python setup.py develop

Look for executable with

which mainprogram.py

# After installation,

Provide path to ~/src/bash folder in ~/.bashrc

export PATH="path_to_HTESP/src/bash:$PATH"

Provide path to src file

export PYTHONPATH="path_to_HTESP/src:$PYTHONPATH"

### License

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


### Contributors

Written and maintained by

#### Niraj K. Nepal (nnepal@ameslab.gov, niraj.nepal@temple.edu)

Postdoctoral Researcher, Ames National Laboratory

#### Lin-Lin Wang

Staff Scientist, Ames National Laboratory 

## Citing HTESP

### Online Documentation in Progress .....
