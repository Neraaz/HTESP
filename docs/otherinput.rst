

.. _inputin-label:

--------------------
input.in
--------------------

.. code-block:: bash

  1

  2 

  400 0

  mpid-list.in

  phband  

  DFT = QE 

In the file input.in, the starting index (1) and ending index (2) are indicated by the first and second rows, respectively. The third row specifies the number of kpoints for band calculations and the cutoff of Brillouin-zone high-symmetry path, counted from the last point. The fourth row provides the filename containing the materials id and compound name for processing. Following this information, the type of plot and the type of calculations are specified.


---------------------
mpid-list.in
---------------------
 
.. code-block:: bash

  v1 mp-763 B2Mg1

  v2 mp-944 B2Al1

First column is just the identifier, second column is material id, and the third column is compound name.

.. _mpid-label:

---------------------
mpid.in
---------------------

Similar file as of :ref:`mpid-list.in <pwd-label>`. But only written after executing "mainprogram download" command. This simply
checks duplication and only updated for new materials id, and acts as the tracking file.

.. _vasp-label:

---------------------
vasp.in
---------------------

Here is an example of :ref:`vasp.in <pwd-label>` file.

.. code-block:: bash

  EDIFF 1E-06

  EDIFFG   -0.01

  NPAR 2

  NSIM 2

  ISMEAR 0

  SIGMA 0.05

  IVDW 11

  NELM 300

  NSW 200

  IALGO 38

  ENCUT 400

  ENAUG 800

  PREC Accurate

  ISPIN 1

  LREAL .False.

  ISIF 3

  IBRION 2

  ALGO

  LASPH

  LMAXMIX

  LORBIT

  LWAVE

  MAGMOM

The keyword with a value is utilized for updating the INCAR file, whereas keywords without values are used to remove them from the INCAR if they exist. Keys with values must be followed by keys without values. For VASP+PHONOPY, a similar file named "vasp-phonopy.in" is utilized, while "vasp-band.in" is employed to generate input files for bandstructure calculations.


.. _kpoint-label:

---------------------
kpoint.in
---------------------

    It is used to change the k-mesh of the input file, using ``mainprogram change_k``. Four options are available for the content of this input.
    
.. code-block:: bash

   2

Here, the existing k-mesh is scaled by ``1/2``. To double it, use ``0.5``.

.. code-block:: bash

   2 0 0 0

In addition to the first case, here we also define k-mesh offset.

.. code-block:: bash

   5 5 5 

Here, the k-mesh is updated with ``5 5 5``.

.. code-block:: bash

   5 5 5 1 1 1

Here the k-mesh is updated with an offset of ``1 1 1``.

The offset for the first and the third case will be ``0 0 0``. The last format is useful for BerkeleyGW calculations.
    

---------------------
qpoint.in
---------------------

    It is used to compute electron-phonon coupling using Density Functional Perturbation Theory (DFPT). Its structure resembles that of :ref:`kpoint.in <kpoint-label>` in the first and third cases. By default, the code uses a q-mesh that is half the size of the k-mesh.


---------------------
ph-q.in
---------------------

The ``ph-q.in`` file is provided for phonon calculation (not electron-phonon coupling) at a particular q point. 
If not specified, provide a ``qpoint.in`` file for direct generic phonon calculation.

The ``ph-q.in`` file consists of the following information on different lines:

.. code-block:: bash

    xq1 xq2 xq3
    metal_info

Where:
- ``xq1 xq2 xq3`` represent the cartesian wavevector coordinates in the unit of 2*pi/(lattice parameters).
- ``metal_info`` indicates the material type. 'T' or 't' (true) for metals.

For example, to perform a phonon calculation for a metal at a specific Gamma point:

.. code-block:: bash

    0 0 0
    T

For non-metallic materials, replace ``T`` with any other character, this will add epsil=.true. for
non magnetic and at q = 0. Don't set this if q != 0 or metallic system. `epsil <https://www.quantum-espresso.org/Doc/INPUT_PH.html#idm72>`_

---------------------
projection.in
---------------------

This file is necessary for specifying projections when running process = epw6-file or epw8-file.

.. code-block:: bash

    X:s;a1;a1;..

    Y:pz;b1;b2;..

    f=0.25,0.25,0.25:s

    f=-0.25,-0.25,-0.25:s

In the above example, X and Y represent different species with various orbital projections separated by ";" 
Additionally, one can utilize different orbitals with coordinates as shown in the third and fourth lines.

---------------------
proj-wt.in
---------------------

This file is necessary for specifying projections for `wanniertools input <http://www.wanniertools.com/input.html#projectors>`_. 

.. code-block:: bash

  X pz px py
  Y pz px py
  Z s

The default order in Wannier90 is s, pz, px, py, dz2, dxz, dyz, dx2-y2, dxy.


.. _pressure-label:

---------------------
pressure.in
---------------------

For QE

.. code-block:: bash

   all

   v1 50 (or .92)

   v2 100 (or .98)

   v3 150 (or 1.06)

Here, the first row represents the `cell_dofree <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm1158>`_ parameter. 'x1' can represent either the scaling factor for the lattice, where the volume is scaled by (x1)^3, or it can denote pressure in kilobars (kbar).


For VASP

.. code-block:: bash

   v1 50 (or .92)

   v2 100 (or .98)

   v3 150 (or 1.06)

Here, we don't have cell_dofree parameters. Default is using pressure in kbar.

.. _charge-input:

---------------------
charge.in
---------------------

This prepares the input files for systems having non-zero net charge.

For QE

.. code-block:: bash

   v1 1

   v2 -1

In QE, the tot_charge <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm289>_ is +1 when one electron is missing and -1 when one electron is added.

For VASP

.. code-block:: bash

   v1 1

   v2 -1

In VASP, the NELECT represents the number of valence electrons, which behaves oppositely to the tot_charge.
