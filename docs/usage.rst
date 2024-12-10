----------------------------
 Package structure
----------------------------

.. code-block:: bash

    Parent folder: HTESP

    Sub folders: docs (Documentation), examples, src, utility (various additional scripts)

    Dependencies: requirements.txt

    installation script: setup.py

    License file: LICENSE

    Readme file: README.md

    Inside src, there is a "bash" folder that has bash scripts for running calculations


----------------------------
Requirements
----------------------------


Current package is tested only for Linux Distribution, with Python and Bash languages.

Basic requirements

Numpy, Scipy, Pandas, Matplotlib, etc.

ASE: https://wiki.fysik.dtu.dk/ase/.

Pymatgen: https://pymatgen.org/.

mp_api: https://next-gen.materialsproject.org/api

qmpy-rester: https://github.com/mohanliu/qmpy_rester

----------------------------
Extra packages
----------------------------

lmfit: https://lmfit.github.io/lmfit-py/. This is needed for SCDM fit to calculate initial projections for wannierization.

IFERMI: https://fermisurfaces.github.io/IFermi/introduction.html#installation. Required for Fermi surface generation.

bsym: https://bsym.readthedocs.io/en/latest/index.html. Required for substitutions of elements in crystal.

Phonopy: https://phonopy.github.io/phonopy/

----------------------------
Download software 
----------------------------

.. code-block:: bash

    git clone https://github.com/Neraaz/HTESP.git

or go to `GitHub <https://github.com/Neraaz/HTESP>`_, download a zip file under code section and unzip it.

Go to the directory,

.. code-block:: bash

    cd HTESP

----------------------------
Conda environment
----------------------------

Make sure the conda is available either via `miniconda <https://docs.anaconda.com/free/miniconda/>`_ or `anaconda <https://www.anaconda.com/download/success>`_ installation

.. code-block:: bash

    conda create --name myenv python==3.10.0 (Please use python version newer than 3.10)
    
    source activate myenv

----------------------------
Install requirements
----------------------------
    
    pip install -r requirements.txt

Also install phonopy in the conda environment


----------------------------
Install HTESP package
----------------------------

.. code-block:: bash

    pip install .

----------------------------
check executable
----------------------------

.. code-block:: bash

    which mainprogram

Do "mainprogram basicinfo" to begin.

Install phonopy to perform phonopy calculations

#Alternatively, install in development version using following command:

.. code-block:: python

    python setup.py develop

.. code-block:: bash

    which mainprogram

----------------------------
After installation,
----------------------------

.. code-block:: bash

    Provide path to ~/src/bash folder in ~/.bashrc
    
    export PATH="path_to_HTESP/src/bash:$PATH"
    
    Provide path to src file
    
    export PYTHONPATH="path_to_HTESP/src:$PYTHONPATH"
    
