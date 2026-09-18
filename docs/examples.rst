.. _examples-label:

=====================
Worked examples
=====================

``examples/`` ships 42 complete tutorials: ``examples/QE/tutorial1`` through
``tutorial21`` and ``examples/VASP/tutorial1`` through ``tutorial21``.  Each
directory holds the ``config.json``, ``input.in``, ``batch.header`` and any other
input the topic needs, a ``README`` naming the commands to run in order, and a
``reference*.tar.gz`` with the expected output.

Copy one into a scratch directory and follow its ``README``, or run the whole
tree with the :doc:`tutorial runner <tutorial_runner>`:

.. code-block:: bash

    cp -r examples/QE/tutorial9 ~/scratch/relax && cd ~/scratch/relax
    # or
    htesp-tutorials --dry-run --only QE/9

---------------------
The numbering
---------------------

Both trees cover the same topics in the same order up to tutorial 10.  From 11
onward the numbers are offset by one: QE/11 is the DFPT electron-phonon
tutorial, which has no VASP counterpart, so **VASP** *n* covers the same topic
as **QE** *n+1* for *n* ≥ 11.  VASP then adds one topic of its own, the IFermi
Fermi surface, as its 21st.

.. list-table::
   :header-rows: 1
   :widths: 10 10 80

   * - QE
     - VASP
     - Topic
   * - 1
     - 1
     - Generate submission scripts from ``batch.header``
   * - 2
     - 2
     - Materials Project search, element mode
   * - 3
     - 3
     - Materials Project search, chemsys mode
   * - 4
     - 4
     - OQMD search and input generation
   * - 5
     - 5
     - AFLOW search and input generation
   * - 6
     - 6
     - Input generation in magnetic configuration
   * - 7
     - 7
     - Combine data from several databases
   * - 8
     - 8
     - Input generation from ``.cif`` files
   * - 9
     - 9
     - Structural relaxation -- the hub every later tutorial starts from
   * - 10
     - 10
     - Cutoff and k-point convergence tests
   * - 11
     - --
     - Electron-phonon coupling and superconducting Tc (QE only)
   * - 12
     - 11
     - Band structure and density of states
   * - 13
     - 12
     - Input files for different pressures or volumes
   * - 14
     - 13
     - Input files with site substitutions
   * - 15
     - 14
     - Elastic constants
   * - 16
     - 15
     - Thermodynamic stability (convex hull)
   * - 17
     - 16
     - Phonon band structure with phonopy
   * - 18
     - 17
     - Equation of state
   * - 19
     - 18
     - Wannier-interpolated band structure
   * - 20
     - 19
     - Input files for a non-zero net charge
   * - 21
     - 20
     - Enumeration of magnetic orderings
   * - --
     - 21
     - 3D Fermi surface with IFermi (VASP only)

This table is the same mapping the tutorial runner's catalogue uses
(``tutorials/catalog.py``), and its self-tests check that the two agree.

---------------------
Before you start
---------------------

* Export ``MP_API_KEY``.  The shipped ``config.json`` files carry the
  placeholder ``use_your_API_KEY``, so the database tutorials -- QE/2, QE/3,
  QE/6, QE/16 and the VASP equivalents -- cannot run without a key of your own.
* Tutorial 9 is the one to do first.  Its relaxed structures are the starting
  point for tutorials 11 onward in both trees.
* For VASP, configure ``PMG_VASP_PSP_DIR`` with ``pmg config``.  POTCAR files
  are licensed and are not shipped.
* ``examples/VASP/tutorial21`` ships only its reference archive; copy a relaxed
  ``vasprun.xml`` into ``R{id}-{name}/relax/`` before running it.
