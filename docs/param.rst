This section has descriptions of files and parameters, required for the package.

.. _pwd-label:

--------------------
Working folder
--------------------

When initializing the calculations, the contents of the working folder would resemble the following:

.. code-block:: bash

    .
    ├── batch.header (content of batch jobmission script except commands)
    ├── config.json (Main input file)
    ├── input.in (Main control file for calculations)
    ├── mpid-list.in (Files to stored {id} and {name} of compounds before download process)
    ├── mpid.in (after download process)
    ├── run-vasp.sh (VASP submission script)
    ├── run-scf.sh, run-elph.sh, run-q2r.sh, ... etc (QE submission scripts)
    ├── pp (QE pseudopotentials)
    │   ├── B.upf
    │   └── Mg.upf
    ├── Rmp-763-B2Mg1 (VASP input files stored here)
    │   └── relax
    │       ├── INCAR
    │       ├── KPOINTS
    │       ├── POSCAR
    │       └── POTCAR
    ├── scf_dir (QE input files stored here)
    │   └── scf-mp-763.in
    ├── Rmp-763-B2Mg1 (unlike VASP, QE folder created only after process = 1)
    │   └── relax
    │       └── scf.in
    └── vasp.in (File to update VASP INCAR file)

.. _batch-label:

-------------------
batch.header
-------------------

This file is used to generate job submission scripts for different calculations.
It holds the scheduler directives and the module loads, and nothing else:

.. code-block:: bash

    #!/bin/bash

    #SBATCH --partition=dense --ntasks=1 --cpus-per-task=48 --time=1-0

    #module load ....

``mainprogram jobscript`` copies this file and appends the run command of each
entry of ``job_script.command_list`` after its last line, writing one
``run-{command}.sh`` per entry.  Nothing has to be reserved in the header for
that.

The older ``generate_submission_file.sh`` script works differently: it
*substitutes* a line reading exactly ``submission here``.  If you use that
script, add such a line at the point where the run command belongs.  The headers
shipped in ``examples/`` do not have one, because they are written for
``mainprogram jobscript``.

Generating a first header
^^^^^^^^^^^^^^^^^^^^^^^^^

The headers in ``examples/`` name a partition (``dense``) and a module that
exist on one cluster and nowhere else, so a copied header is usually rejected
by the scheduler before any calculation starts.  ``--init-header`` writes one
from what *this* machine reports instead:

.. code-block:: bash

    mainprogram jobscript --init-header qe      # or: vasp

``qe``, ``QE``, ``QuantumEspresso``, ``quantum-espresso`` and ``espresso`` all
name the same code; the module search behind it looks for each of those
spellings too, because a site that installs ``QuantumESPRESSO/7.1`` would
otherwise get a ``TODO`` where the ``module load`` belongs.

Module dependencies
^^^^^^^^^^^^^^^^^^^

Most HPC centres run a **hierarchical** Lmod, where an application module is
not visible at all until the compiler and MPI it was built against are loaded.
On TACC Vista ``qe/7.3`` lives under
``/opt/apps/nvidia24/openmpi5/modulefiles``, so a job script containing only
``module load qe`` dies on its first line with *"these module(s) exist but
cannot be loaded as requested"*.  Typing the same command interactively works,
because the login shell already has the toolchain -- which is exactly the kind
of difference that only shows up once the job has been queued.

``--init-header`` therefore asks ``module spider <name>/<version>`` which
combinations make the module available, picks the one sharing most modules
with the current environment (that combination is demonstrably working here),
and writes it out first:

.. code-block:: bash

    # qe/7.3 is built against a specific compiler/MPI; 'module spider qe/7.3'
    # reports this combination, so load it first:
    #   nvidia/24.7  cuda/12.5  openmpi/5.0.5
    module load nvidia cuda openmpi
    module load qe

The prerequisites are written unversioned for the same reason the code module
is: Lmod resolves each to the site default, and the exact combination found is
kept in the comment above them for pinning by hand.

Two sources are read, in this order:

``module help`` -- asked first, and believed
    This is whoever built the module speaking.  On Bridges-2 it says::

        > module load intel-oneapi QuantumEspresso/7.5-intel

    ``intel-oneapi`` carries the Intel MPI and MKL runtimes that build is
    linked against; Lmod's hierarchy does not model it at all, and that site
    reports ``This module can be loaded directly``.  Loading Quantum ESPRESSO
    alone there puts ``pw.x`` on ``PATH`` and then fails at run time on a
    missing shared library -- a much more confusing failure than a module that
    refuses to load.  The exact build is asked first and the bare name second,
    since ``module help qe`` resolves to the site default.

``module spider`` -- only when help names no prerequisite
    The hierarchy block lists every combination that makes the module
    visible; the one sharing most modules with the current environment is
    chosen, because that combination is demonstrably working here.  Its own
    copy of the Help block is read too, and merged.

"Help first" is a preference, not an exclusion.  Most help text says no more
than ``module load qe/7.3`` -- Vista's does -- and there ``spider`` is the only
thing that reveals ``nvidia cuda openmpi``.  Skipping ``spider`` whenever help
merely *exists* would break every hierarchical site.

Only the tokens *before* the module itself are taken, and only from a line
that names it -- so Vista's help, which says plainly ``module load qe/7.3``,
correctly yields nothing.  Prose is rejected: every token on a candidate line
has to look like ``name`` or ``name/version``, or the ``To run codes in
quantum espresso include the following lines`` that follows the hierarchy
block would parse as eight modules.

The version probed is the **newest** build, since its help describes the
toolchain the site currently expects; what the header actually loads is still
the bare name.

.. warning::

   A generated header is a starting point, not a working job script, and
   ``--init-header`` says so every time it writes one.  Read the file and make
   sure every module the build needs is loaded, **dependencies included** --
   a chain that is one module short is accepted by the scheduler and then
   fails inside the job, minutes later, with an error naming a shared library
   rather than a module.  The quickest check is a login shell::

       source batch.header && which pw.x      # or vasp_std

   If that prints nothing, the module chain is incomplete.

It asks SLURM for the partitions (``sinfo``), the accounts this user may
charge (``sacctmgr``), the cores per node of the chosen partition, and whether
any generic resource is configured at all (``scontrol show config``); it asks
Lmod for the ``qe``/``quantum-espresso`` or ``vasp`` modules, preferring the
one Lmod marks ``(D)``; and it looks for ``ibrun``, ``srun`` or ``mpirun`` on
``$PATH``.  Anything it cannot determine is written as a ``# TODO`` comment
rather than guessed, and the other partitions and accounts it found are listed
in comments so the choice can be changed by uncommenting.

No run command is written into the header.  ``mainprogram jobscript`` appends
that line itself, from ``job_script.parallel_command`` and ``job_script.nproc``
in ``config.json``; a command in the header would be run *as well*, before the
one that matters.  The launcher the generator found is therefore reported as a
comment naming the two configuration values to set.

Two limits are worth knowing.  ``--gres=gpu:N`` is emitted only on a cluster
whose SLURM configuration actually defines GRES types; a GPU machine that
schedules whole nodes reports ``GresTypes = (null)`` and gets no ``--gres``
line, because emitting one there produces a job the scheduler refuses.  And
the wall time and node count are placeholders -- the generator knows nothing
about the size of the study.  Read the file before submitting anything with
it.

An existing ``batch.header`` is never overwritten; pass ``--force`` to replace
one.

   
.. _json-label:

-------------------
config.json
-------------------

This `JSON <https://docs.python.org/3/library/json.html>`_ file serves as the main input file of the package. It contains a dictionary with various keys to configure different aspects of the package's functionality.

The same mapping in JSON and in Python:

.. code-block:: json

    {
      "name": "John",
      "age": 30,
      "is_student": false,
      "favorite_fruits": ["apple", "banana", "orange"],
      "address": null
    }

.. code-block:: python

    {
        "name": "John",
        "age": 30,
        "is_student": False,
        "favorite_fruits": ["apple", "banana", "orange"],
        "address": None,
    }

In JSON, boolean values are written ``true`` and ``false``, where Python writes
``True`` and ``False``.  JSON writes ``null`` where Python writes ``None``.
Lists use square brackets ``[]`` in both, with comma-separated values.

Where the file is read from
---------------------------

``config.json`` is looked for in the working directory and in up to five parent
directories, so a helper running in ``R<id>-<name>/pressure/R.../relax`` finds
the campaign's file without a copy in every folder.  ``$HTESP_CONFIG`` (or
``--config FILE``) overrides the search and may name either a file or a
directory.

Whatever is found is deep-merged over the packaged default,
``htesp/data/config.json``, which carries every key the code reads.  A file
written for an older version of HTESP therefore keeps working: keys it does not
mention come from the default rather than raising ``KeyError``.  The merged
result is cached, so a long campaign parses the file once.

Check a file before a campaign with:

.. code-block:: bash

    mainprogram config-validate

It prints the path actually used, then every problem it can find: missing
sections, values of the wrong type, an ``elph_mode`` that is not one of the four
accepted values, and whether a Materials Project API key is available.  It exits
non-zero when something is wrong.

The full shipped configuration
------------------------------

This is ``htesp/data/config.json`` verbatim -- the default every user file is
merged over, and a working starting point to copy into the working directory.

.. config-json-start
   Generated from htesp/data/config.json by docs/gen_param_block.py.
   Run that script after changing the shipped default; do not edit by hand.

.. code-block:: json

    {
      "job_script": {
        "batch": "batch.header",
        "which_calc": "qe",
        "parallel_command": "mpirun",
        "nproc": "4",
        "command_list": [
          "scf",
          "elph"
        ],
        "command_combine": false,
        "calc_visible_with": "id"
      },
      "mpi_key": {
        "API_KEY": {
          "key": "use_your_API_KEY"
        }
      },
      "download": {
        "mode": "element",
        "element": {
          "metal": false,
          "FE": false,
          "thermo_stable": false,
          "exclude": [
            "Lu"
          ],
          "ntype": [
            1,
            2
          ],
          "elm": [
            "B"
          ],
          "prop": [
            "material_id",
            "formula_pretty",
            "structure",
            "formation_energy_per_atom",
            "band_gap",
            "energy_above_hull",
            "total_magnetization",
            "ordering",
            "total_magnetization_normalized_formula_units",
            "num_magnetic_sites",
            "theoretical",
            "nsites"
          ],
          "ordering": [
            "NM",
            "Unknown"
          ],
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
          "use_cif2cell": false,
          "kpath_pbc": null
        },
        "chemsys": {
          "entries": [
            "Mg",
            "B"
          ],
          "size_constraint": 60,
          "ntype_constraint": 3,
          "must_include": [
            "Mg",
            "B"
          ],
          "FE": false,
          "thermo_stable": 0.08,
          "metal": false,
          "magnetic": true,
          "spacegroup": null
        },
        "oqmd": {
          "limit": 1000,
          "entries": [
            "Mg",
            "B"
          ],
          "size_constraint": 60,
          "ntype_constraint": 3,
          "must_include": [],
          "metal": false,
          "magnetic": true,
          "spacegroup": null,
          "thermo_stable": true,
          "FE": true,
          "prop": [
            "composition",
            "spacegroup",
            "volume",
            "band_gap",
            "stability"
          ]
        },
        "aflow": {
          "elm": [
            "Mg",
            "B"
          ],
          "nelm": 2,
          "nsites": 60,
          "metal": false,
          "FE": false,
          "spacegroup": null,
          "filter": false,
          "limit": 5000,
          "prop": [
            "spacegroup_relax",
            "Pearson_symbol_relax"
          ]
        }
      },
      "conv_test": {
        "param": "ecut",
        "ecut": [
          400,
          500,
          600
        ],
        "kpoint": [
          [
            6,
            6,
            6
          ],
          [
            12,
            12,
            12
          ],
          [
            18,
            18,
            18
          ]
        ]
      },
      "magmom": {
        "magmom": {
          "Cr": 5,
          "Pd": 0,
          "I": 0
        },
        "type": "anisotropy",
        "saxis": [
          [
            0,
            0,
            1
          ],
          [
            1,
            0,
            0
          ],
          [
            1,
            1,
            0
          ],
          [
            1,
            1,
            1
          ]
        ],
        "order": [
          "ferromagnetic",
          "antiferromagnetic",
          "ferrimagnetic_by_motif"
        ],
        "force_theorem": false
      },
      "pseudo": {
        "pot": {
          "H": "H",
          "He": "He",
          "Li": "Li_sv",
          "Be": "Be",
          "B": "B",
          "C": "C",
          "N": "N",
          "O": "O",
          "F": "F",
          "Ne": "Ne",
          "Na": "Na_pv",
          "Mg": "Mg",
          "Al": "Al",
          "Si": "Si",
          "P": "P",
          "S": "S",
          "Cl": "Cl",
          "Ar": "Ar",
          "K": "K_sv",
          "Ca": "Ca_sv",
          "Sc": "Sc_sv",
          "Ti": "Ti_sv",
          "V": "V_sv",
          "Cr": "Cr_pv",
          "Mn": "Mn_pv",
          "Fe": "Fe",
          "Co": "Co",
          "Ni": "Ni",
          "Cu": "Cu",
          "Zn": "Zn",
          "Ga": "Ga_d",
          "Ge": "Ge_d",
          "As": "As",
          "Se": "Se",
          "Br": "Br",
          "Kr": "Kr",
          "Rb": "Rb_sv",
          "Sr": "Sr_sv",
          "Y": "Y_sv",
          "Zr": "Zr_sv",
          "Nb": "Nb_sv",
          "Mo": "Mo_sv",
          "Tc": "Tc_pv",
          "Ru": "Ru_pv",
          "Rh": "Rh_pv",
          "Pd": "Pd",
          "Ag": "Ag",
          "Cd": "Cd",
          "In": "In_d",
          "Sn": "Sn_d",
          "Sb": "Sb",
          "Te": "Te",
          "I": "I",
          "Xe": "Xe",
          "Cs": "Cs_sv",
          "Ba": "Ba_sv",
          "La": "La",
          "Ce": "Ce",
          "Pr": "Pr_3",
          "Nd": "Nd_3",
          "Pm": "Pm_3",
          "Sm": "Sm_3",
          "Eu": "Eu_2",
          "Gd": "Gd_3",
          "Tb": "Tb_3",
          "Dy": "Dy_3",
          "Ho": "Ho_3",
          "Er": "Er_3",
          "Tm": "Tm_3",
          "Yb": "Yb_2",
          "Lu": "Lu_3",
          "Hf": "Hf_pv",
          "Ta": "Ta_pv",
          "W": "W_sv",
          "Re": "Re",
          "Os": "Os",
          "Ir": "Ir",
          "Pt": "Pt",
          "Au": "Au",
          "Hg": "Hg",
          "Tl": "Tl_d",
          "Pb": "Pb_d",
          "Bi": "Bi_d",
          "Po": "Po_d",
          "At": "At",
          "Rn": "Rn",
          "Fr": "Fr_sv",
          "Ra": "Ra_sv",
          "Ac": "Ac",
          "Th": "Th",
          "Pa": "Pa",
          "U": "U",
          "Np": "Np",
          "Pu": "Pu",
          "Am": "Am",
          "Cm": "Cm"
        },
        "PSEUDO": {
          "H": 60,
          "He": 50,
          "Li": 40,
          "Be": 40,
          "B": 35,
          "C": 45,
          "N": 60,
          "O": 60,
          "F": 45,
          "Ne": 50,
          "Na": 40,
          "Mg": 30,
          "Al": 30,
          "Si": 30,
          "P": 30,
          "S": 35,
          "Cl": 40,
          "Ar": 60,
          "K": 60,
          "Ca": 30,
          "Sc": 40,
          "Ti": 35,
          "V": 35,
          "Cr": 40,
          "Mn": 65,
          "Fe": 90,
          "Co": 45,
          "Ni": 45,
          "Cu": 55,
          "Zn": 40,
          "Ga": 70,
          "Ge": 40,
          "As": 35,
          "Se": 30,
          "Br": 30,
          "Kr": 45,
          "Rb": 30,
          "Sr": 30,
          "Y": 35,
          "Zr": 30,
          "Nb": 40,
          "Mo": 35,
          "Tc": 30,
          "Ru": 35,
          "Rh": 35,
          "Pd": 45,
          "Ag": 50,
          "Cd": 60,
          "In": 50,
          "Sn": 60,
          "Sb": 40,
          "Te": 30,
          "I": 35,
          "Xe": 60,
          "Cs": 30,
          "Ba": 30,
          "La": 40,
          "Ce": 50,
          "Pr": 40,
          "Nd": 40,
          "Pm": 40,
          "Sm": 40,
          "Eu": 40,
          "Gd": 40,
          "Tb": 40,
          "Dy": 40,
          "Ho": 40,
          "Er": 40,
          "Tm": 40,
          "Yb": 40,
          "Lu": 45,
          "Hf": 50,
          "Ta": 45,
          "W": 30,
          "Re": 30,
          "Os": 40,
          "Ir": 55,
          "Pt": 35,
          "Au": 45,
          "Hg": 50,
          "Tl": 50,
          "Pb": 40,
          "Bi": 45,
          "Po": 30,
          "At": 30,
          "Rn": 30,
          "Fr": 30,
          "Ra": 30,
          "Ac": 30,
          "Th": 30,
          "Pa": 30,
          "U": 30,
          "Np": 30,
          "Pu": 30,
          "Am": 30,
          "Cm": 30,
          "Bk": 60,
          "Cf": 60,
          "Es": 60,
          "Fm": 60,
          "Md": 60,
          "No": 60,
          "Lr": 60
        }
      },
      "substitute": {
        "mode": 2,
        "elm": "Cr",
        "sub": {
          "Cr": 0,
          "Al": 1
        },
        "new_sub": {
          "Cr": "Fe",
          "Pd": "Pd",
          "I": "I"
        }
      },
      "pwscf_in": {
        "magnetic": false,
        "control": {
          "calculation": "vc-relax",
          "nstep": 300,
          "restart_mode": "from_scratch",
          "pseudo_dir": "../../pp/",
          "outdir": "./",
          "tprnfor": true,
          "tstress": true,
          "etot_conv_thr": 1e-05,
          "forc_conv_thr": 0.0001
        },
        "system": {
          "smearing": "gauss",
          "occupations": "smearing",
          "degauss": 0.02
        },
        "electrons": {
          "diagonalization": "david",
          "mixing_mode": "plain",
          "mixing_beta": 0.7,
          "conv_thr": 1e-16,
          "electron_maxstep": 300
        }
      },
      "strain": [
        -0.01,
        -0.005,
        0.005,
        0.01
      ],
      "wanniertools_input": {
        "tb_file": {
          "Hrfile": "'ex_hr.dat'",
          "Package": "'QE'"
        },
        "control": {
          "BulkBand_calc": "T",
          "BulkFS_calc": "F",
          "BulkGap_cube_calc": "F",
          "BulkGap_plane_calc": "F",
          "FindNodes_calc": "F",
          "SlabBand_calc": "F",
          "WireBand_calc": "F",
          "Dos_calc": "F",
          "JDos_calc": "F",
          "SlabSS_calc": "F",
          "SlabArc_calc": "F",
          "SlabQPI_calc": "F",
          "SlabSpintexture_calc": "F",
          "wanniercenter_calc": "F",
          "Z2_3D_calc": "F",
          "Chern_3D_calc": "F",
          "WeylChirality_calc": "F",
          "BerryPhase_calc": "F",
          "BerryCurvature_calc": "F",
          "AHC_calc": "F"
        },
        "parameters": {
          "E_arc": "0.0",
          "Eta_Arc": "0.001",
          "OmegaMin": "0.0",
          "OmegaMax": "0.0",
          "OmegaNum": "100",
          "Nk1": "10",
          "Nk2": "10",
          "Nk3": "10",
          "NP": "2",
          "Gap_threshold": "0.1"
        },
        "system": {
          "NSlab": "10",
          "NSlab1": "1",
          "NSlab2": "1",
          "NumOccupied": "1",
          "SOC": "1",
          "E_FERMI": "0.0",
          "Bx": "0.0",
          "By": "0.0",
          "Bz": "0.0",
          "surf_onsite": "0.0"
        },
        "surface": {
          "surface": [
            [
              1,
              0,
              0
            ],
            [
              0,
              1,
              0
            ],
            [
              0,
              0,
              1
            ]
          ],
          "KPATH_SLAB": {},
          "KPLANE_SLAB": {},
          "EFFECTIVE_MASS": "0.0",
          "SELECTED_ATOMS": {}
        }
      },
      "kptden": 0.025,
      "chull_cutoff": 0.04,
      "kpt_opt": true,
      "elph_mode": "serial",
      "plot": {
        "xlim": [
          -8,
          10
        ],
        "ylim": [
          -5,
          5
        ],
        "atomproj": 0.6,
        "bandproj": {
          "proj_type": "element-orbital",
          "proj": {
            "Zr": "dxz",
            "Pd": "dxz"
          },
          "colormap": "Reds"
        },
        "a2f_smearing": null,
        "dos_ylim": null
      }
    }

.. config-json-end


- **job_script**: Information about creating job submission scripts. :ref:`here <job-label>`

- **mpi_key**: Legacy home of the Materials Project API key. Set ``$MP_API_KEY`` instead; leave the placeholder in the file. Without a key, extraction from the Materials Project is not possible, but the OQMD and AFLOW databases are still accessible. :ref:`here <mpikey-label>`

- **download**: Information required for downloading and preparing inputs. :ref:`here <download-label>`

- **conv_test**: Parameters for performing convergence tests.

- **magmom**: Specifies magnetic moment configurations. This key is essential for magnetic calculations. :ref:`here <magmom-label>`

- **pseudo**: Pseudopotential file, with two subkeys: ``pot`` for VASP and ``PSEUDO`` for QE. :ref:`here <pseudo-label>`

- **substitute**: Parameters for preparing input files with substitution. :ref:`here <substitute-label>`

- **pwscf_in**: Input parameters for Quantum Espresso (QE). :ref:`here <pwscfin-label>`

- **strain**: Strain inputs, useful for elastic properties calculations.

- **wanniertools_input**: Input parameters for WannierTools. :ref:`here <wanniertoolsinput-label>`

- **kptden**: Kpoint density.

- **chull_cutoff**: It controls the display of thermodynamically unstable systems in convexhull plot. It can be logical or value in eV/atom. Default is ``false``, doesn't show systems with energy above hull (EBH) larger than zero. If ``true``, it shows all the systems in consideration. If it has the value ``0.04``, it shows the systems with EBH less than ``0.04 eV/atom``. 

- **kpt_opt**: Flag for generating KPOINTS_OPT file to perform band structure in VASP

- **elph_mode**: Parameter for performing DFPT calculations in various modes in QE.

  - **serial**: Serial mode (Default)

  - **parallel_q**: Parallel over q-points

  - **parallel_irr**: Parallel over q-points and irreducible representations

  - **only_init**: Run only the initialisation that writes the irreducible
    representations.  This is a fourth value of ``elph_mode``, not a sub-option
    of ``parallel_irr``: set ``"elph_mode": "only_init"`` and run process ``7``
    first, then set ``"elph_mode": "parallel_irr"`` and run process ``7`` again.

- **plot**: Plot variables: the energy window, the projections and the two
  optional keys described under :ref:`plot <plot-label>`.


.. _job-label:

--------------------
job_script
--------------------

.. code-block:: text

    "job_script": {
      "batch": "batch.header",
      "which_calc": "qe",
      "parallel_command": "mpirun",
      "nproc": "1",
      "command_list": ["scf", "elph"],
      "command_combine":false,
      "calc_visible_with":"id"}

- **batch**: Usual batch script

- **which_calc**: Type of calculations. One of ``"qe"``, ``"vasp"``, ``"wannier"`` or ``"epw"``. The value is normalised before use, so ``"QE"``, ``"VASP"``, ``"WANNIER"`` and ``"EPW"`` work exactly as well as the lower-case spellings. Anything else is rejected with a message naming the four accepted values.

- **parallel_command**: The MPI launcher, e.g. ``"mpirun"``, ``"srun"`` or ``"ibrun"``. How ``nproc`` is passed to it follows the launcher: the ``mpirun`` family (``mpirun``, ``mpiexec``, ``mpiexec.hydra``, ``orterun``) gets ``-np N``, ``srun`` and ``aprun`` get ``-n N``, and ``ibrun`` is given no count at all because it runs the whole SLURM allocation and rejects one. A launcher this list does not know gets ``-np N``, which is right for the mpirun-family wrappers sites usually install. If the value already contains a flag of its own -- ``"srun --cpu-bind=cores -n 8"`` -- it is used exactly as written and no count is added. An empty value runs the executable directly, with no launcher.

- **nproc**: Number of processes to launch. With ``ibrun`` it does not reach the launcher, so there it only has to agree with the ``#SBATCH`` directives in :ref:`batch.header <batch-label>`.

- **command_list**: List of command to execute.

  - **QE commands**:

    - Work with ``"which_calc":"qe"`` or ``"QE"``

    - **scf**: scf calculations, ``pw.x < scf.in > scf.out``

    - **band**: band calculations, ``pw.x < scf-band.in > nscf-band.out``

    - **bandp**: band processing, ``bands.x < band.in > band.out``

    - **dos**: Density of States (DOS) calculations, ``pw.x < scf-dos.in > scf-dos.out``

    - **dosp**: DOS postprocessing, ``dos.x < dos.in > dos.out``

    - **elph**: DFPT electron-phonon coupling (EPC) calculations, ``ph.x < elph.in > elph.out``

    - **q2r**: Converting force constant in reciprocal to real space, ``q2r.x < q2r.in > q2r.out``

    - **matdyn**: Calculating phonon frequencies in generic q points, ``matdyn.x < matdyn.in > matdyn.out``

    - **matdyn-dos**: Calculating phonon DOS, ``matdyn.x < matdyn-dos.in > matdyn-dos.out``

    - **lambda**: Calculating EPC strength constant, lambda. ``lambda.x < lambda.in > lambda.out``

    - **pdos**: Calculating partial DOS, ``projwfc.x < pdos.in > pdos.out``

  - **VASP commands**:

    - Work with ``"which_calc":"vasp"`` or ``"VASP"``

    - **vasp**: Performing VASP calculations, ``vasp_std``. Modify the code in ``generate_submission.py`` as needed, and then proceed to reinstall the package if you require commands other than ``vasp_std``.

    - **wannier**: Performing WANNIER90 calculations with VASP. This add ``wannier90.x wannier90`` in ``batch.header`` file.

    - **ifermi**: Utilizing ifermi package to compute Fermi surface related properties, ``'ifermi'``. Please locate the ``ifermi.json`` file in the ``utility/input_files`` directory and move it to the current working directory.

    - **Note**: Please checkout original documentation of the `ifermi package <https://fermisurfaces.github.io/IFermi/cli.html>`_. Use ``true`` or ``false`` for keys that don't have values. When set to ``true,`` the key will be included as a flag in the command.

.. _ifermi:

ifermi.json
^^^^^^^^^^^

.. code-block:: json

            {
                "info": {
                    "--filename": "vasprun.xml",
                    "--mu": 0.0,
                    "--wigner": false,
                    "--interpolation-factor": 8.0,
                    "--property": "velocity",
                    "--projection-axis": "0 0 1",
                    "--decimate-factor": 0.8,
                    "--smooth": false,
                    "--norm": false,
                    "--precision": 4
                },
                "plot": {
                    "--filename": "vasprun.xml",
                    "--mu": 0.0,
                    "--wigner": false,
                    "--interpolation-factor": 8.0,
                    "--property": "velocity",
                    "--projection-axis": "0 0 1",
                    "--decimate-factor": 0.8,
                    "--smooth": false,
                    "--output": "output.png",
                    "--symprec": 0.001,
                    "--azimuth": 45.0,
                    "--elevation": 35.0,
                    "--type": "plotly",
                    "--color-property": false,
                    "--property-colormap": "viridis",
                    "--vector-property": false,
                    "--vector-colormap": "plasma",
                    "--vector-spacing": 0.2,
                    "--cmin": 0.0,
                    "--cmax": 1.0,
                    "--vnorm": 1.0,
                    "--scale-linewidth": false,
                    "--hide-surface": false,
                    "--plot-index": 1,
                    "--hide-labels": false,
                    "--hide-cell": false,
                    "--spin": "up",
                    "--slice": "0.5 0.5 0.5 0.0",
                    "--scale": 4
                }
            }

      


  - **WANNIER commands**:

    - Work with ``"which_calc":"wannier"`` or ``"WANNIER"``

    - **scf**: QE scf command, ``pw.x < scf.in > scf.out``

    - **nscf**: QE nonscf command, ``pw.x < scf.in > nscf.out``

    - **wannier_prepare**: Preparing for wannier calculations, ``wannier90.x -pp ex``

    - **pw2wannier90**: Generating input files for WANNIER90 from QE output, ``pw2wannier90.x -in pw2wan.in > pw2wan.out``

    - **wannier_band**: Performing wannierization, ``wannier90.x ex``

  - **EPW commands**:

    - Work with ``"which_calc":"epw"`` or ``"EPW"``

    - **scf**: QE scf command, ``pw.x < scf.in > scf.out``

    - **ph**: QE phonon calculations, ``ph.x < elph.in > elph.out``

    - **proj**: QE nonscf calculations to run with projwfc.x, ``pw.x < nscf-proj.in > nscf-proj.out``

    - **epw_nscf**: QE nonscf calculations for EPW calculations, ``pw.x < nscf_epw.in > nscf_epw.out``

    - **epw**: EPW calculations, ``epw.x -npools -nproc -i epw.in > epw.out``

- **command_combine**: When set to true, commands will be consolidated into one file named run-{last_command}.sh, following the order specified in the command_list. Otherwise, each command will be written to separate files named :ref:`run-{command}.sh <pwd-label>`.

- **calc_visible_with**: It determines the naming convention for the job submission scripts after submission. Available options include ``"id"``, ``"name"``, or ``"id-name"``. If not provided, the scripts will be named simply as run-{command}.sh files for QE, and run.sh for VASP calculations.

  - **id**: Materials_id, ``CALC_VISIBLE_WITH_ID`` file created.

  - **name**: Compound_name, ``CALC_VISIBLE_WITH_NAME`` file created.

  - **id-name**: Materials_id-Compound_name, ``CALC_VISIBLE_WITH_ID-NAME`` file created.

  - If you prefer not to display the compounds information, please provide an empty string ``""``.


.. _mpikey-label:

--------------------
mpi_key
--------------------

The Materials Project API key is **not** kept in ``config.json``.  Export it
instead:

.. code-block:: bash

    export MP_API_KEY=your_key_here

Get a key from https://next-gen.materialsproject.org/api#api-key, under the API
key section.

For a permanent setting without putting the key in your shell history, write it
to ``~/.config/htesp/credentials``:

.. code-block:: bash

    mkdir -p ~/.config/htesp
    printf 'MP_API_KEY=your_key_here\n' > ~/.config/htesp/credentials
    chmod 600 ~/.config/htesp/credentials

The key is looked up in this order:

1. ``$MP_API_KEY``;
2. a ``MP_API_KEY=...`` line in ``~/.config/htesp/credentials``;
3. ``mpi_key.API_KEY.key`` in ``config.json``.

The ``mpi_key`` section still exists so that old configurations load, and every
shipped example carries the placeholder ``"use_your_API_KEY"``, which is treated
as "no key".  Leave it as the placeholder: a key committed to a repository or
copied between working directories is a key you have to revoke.

Without a key, ``mainprogram search`` and ``mainprogram download`` stop with a
message telling you how to set one.  The OQMD and AFLOW searches need no key.

.. _download-label:

--------------------
download
--------------------

It has a dictionary of the form.

.. code-block:: text

    "download": {
        "mode": "element",
        "element": {
          "metal": false,
          "FE": false,
          "thermo_stable": false,
          "exclude": ["Lu"],
          "ntype": [1, 2],
          "elm": ["B"],
          "prop": ["material_id", "formula_pretty", "structure", "formation_energy_per_atom", "band_gap", "energy_above_hull", "total_magnetization", "ordering", "total_magnetization_normalized_formula_units", "num_magnetic_sites", "theoretical", "nsites"],
          "ordering": ["NM", "Unknown"],
          "nsites": 10,
          "spacegroup": null},
        "inp": {
          "start": 1,
          "end": 2,
          "nkpt": 200,
          "evenkpt": false,
          "plot": "phband",
          "calc": "QE",
          "use_cif2cell": false,
          "kpath_pbc": null},
        "chemsys": {
          "entries": ["Fe", "Pd", "I"],
          "size_constraint": 60,
          "ntype_constraint": 4,
          "must_include": ["Fe", "Pd", "I"],
          "FE": false,
          "thermo_stable": 0.08,
          "metal": false,
          "magnetic": true,
          "spacegroup": null},
        "oqmd": {
          "limit": 400,
          "entries": ["Mg", "B"],
          "size_constraint": 60,
          "ntype_constraint": 4,
          "must_include": [],
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
          "filter": false,
          "prop": ["spacegroup_relax", "Pearson_symbol_relax"]}}

- **(A)mode**:

  - The keyword ``'mode'`` controls the preparation of input files.

    - There are 2 modes for preparing input files from materials project database. In addition, input preparation utilizing ``.cif`` and ``.vasp`` structure files are also available:

      - **element**: Extracts data from the `Materials Project (MP) <https://next-gen.materialsproject.org/>`_ database using element-based search method. It uses the parameters within the ``element`` dictionary.

      - **chemsys**: Extracts data from the MP database using chemsys (e.g., Mg-O) based search method, utilizing parameters from ``chemsys`` dictionary.

      - **fromcif**: Converts generic .cif files into Quantum Espresso (QE) and VASP input files.

      - **fromvasp**: Converts POSCAR files in .vasp format into VASP input files.

      - **Note**: Once you have finished generating the input, set the flag to empty string ``""`` if you wish only to update the ``INCAR`` file according to ``vasp.in``.

- **(B)element**: 

  - **metal**: 
    - ``true`` selects zero bandgap systems.
  
  - **FE**: 
    - ``true`` selects compounds with negative formation energy.
  
  - **thermo_stable**: 
    - ``true`` selects compounds on the convex hull.
  
  - **exclude**: 
    - List of elements to exclude from the search.
  
  - **ntype**: 
    - Number of different types of elements in compounds.
  
  - **elm**: 
    - The element to search in compounds. Can use up to 2 elements. For example, ``['B', 'C']`` to search boron and/or carbon-containing compounds.
  
  - **prop**: 
    - List of properties to extract.
  
  - **ordering**: 
    - Magnetic ordering to keep. Accepts three forms:
    - a **string** -- ``"NM"`` nonmagnetic, ``"FM"`` ferromagnetic, ``"AFM"`` antiferromagnetic, ``"FiM"`` ferrimagnetic, ``"Unknown"`` (see below);
    - a **list** -- keep any of its values, e.g. ``["NM", "Unknown"]``. **This is the shipped default.**
    - ``null`` -- apply no ordering filter at all.
    - **About** ``"Unknown"``: Materials Project reports it for any material that has no magnetism calculation, and that is now a large fraction of the database. In a search over boron binaries it was 202 of 465 hits, and 158 of the 222 that survived the metal and formation-energy filters. ``"Unknown"`` therefore means "nobody computed the ordering", not "this material is magnetic" -- so excluding it silently discards most of the database. A bare ``"NM"`` is easy to mistake for a broken search: 465 compounds in, 2 rows out.
    - Use ``"NM"`` alone only when you want materials explicitly determined to be nonmagnetic, and are content to lose the uncomputed ones.
  
  - **nsites**: 
    - Total number of ions in the compound.
  
  - **spacegroup**: 
    - Spacegroup symbol to choose. If ``null``, then the code doesn’t select based on spacegroup.

- **(B)inp**:

  - **start**: Determines the starting index within various tracking files identified by the "mpid".

  - **end**: Determines the ending index (exclusive) within various tracking files identified by the "mpid". Tracking file may be  mpid-list.in, used as in :ref:`input.in <inputin-label>`.

  - **nkpt**: Specifies the total number of k-points in the Brillouin zone sample along the high symmetry path, crucial for electronic and phonon bandstructure calculations, particularly when VASP line mode is not utilized.

  - **evenkpt**: Indicates whether an even k-mesh should be used, particularly relevant for electron-phonon coupling calculations where selecting k-mesh integer multiples of q is beneficial.

  - **plot**: Defines the type of plot to be generated. Options include:

    - **eband**: Bandstructure excluding line-mode kpath, applicable to QE/VASP.

    - **vasp-line**: Bandstructure for VASP line mode k-path.

    - **phband**: Phonon band for QE.

    - **pdos**: Density of states and partial DOS, applicable to QE/VASP.

    - **gammaband**: Electron-phonon coupling strength projected phonon band for QE-DFPT.

    - **wann_band**: Bandstructure from wannier calculation.

    - **phonproj**: Atom-projected phonon for QE-DFPT. Provide value for ``atomproj`` cutoff for plotting atomic projections on phonon dispersion.

    - **bandproj**: Atom or orbital projected band structure. (Turn on only after processing/plotting band structure with eband/vasp-line)

  - **calc**: Type of DFT calculations. Options are QE/VASP.

  - **use_cif2cell**: If True, cif2cell package will be utilized to read .cif files. Install `cif2cell <https://pypi.org/project/cif2cell/>`_ to use this function.

  - **kpath_pbc**: Determine periodic boundary condition for kpath with ``default [1, 1, 1] for 3D``. Generating kpath for 2D systems with vacuum along z-direction requires [1, 1, 0]

- **(C)chemsys**:

  - The "chemsys" keyword mirrors the construction of the Materials Project database and is utilized to search for compounds.

  - **entries**:

    - It employs the "mp_api.client.MPRester.get_entries_in_chemsys" function to explore atoms, binary, ternary, and other combinations based on the "entries" keyword.
    - This functionality is valuable in studying thermodynamic stability using convex hull phase diagrams.
    - Besides "entries", other keys within chemsys include:

  - **size_constraint**: 
    - ``60``: Maximum size allowed. Search includes structures with less than 60 ions per cell.
  
  - **ntype_constraint**: 
    - ``4``: Limit of different species allowed. Search includes structures with less than 4 different species.
  
  - **must_include**: 
    - ``["Fe", "Pd", "I"]``: All of these elements are included in the search. If only ``["Fe", "Pd"]``, then "I" only elements are not searched.
  
  - **FE**: 
    - ``false``: If ``true``, compounds with negative formation energies are searched. Otherwise, include all.

  - **thermo_stable**:
    - If ``true``, compounds with energy above hull less than 0.001 eV/atom is searched.
    - If ``0.08``, compounds with energy above hull less than 0.08 eV/atom is searched.
    - If It holds the false value, no filter is applied.
  
  - **metal**: 
    - ``false``: If ``true``, compounds with zero bandgap are searched. Otherwise, include all.
  
  - **magnetic**: 
    - ``true``: If ``true``, all magnetic orderings are searched.
  
  - **spacegroup**: 
    - ``null``: Spacegroup symbol to search. Otherwise, include all.

- **(D)oqmd**:

  - Extract data from `OQMD <https://oqmd.org/>`_  database, utilizing parameters from ``oqmd`` dictionary.

  - It needs `qmpy-rester <https://pypi.org/project/qmpy-rester/0.1.8/>`_ API package.

  - **limit**: 400 - Specifies the initial limit for the search. If the search doesn't succeed, the number is decreased by a factor of 0.2 of this limit and the search is performed again

  - **entries**: ["Mg", "B"] - List of elements to include in the search

  - **size_constraint**: 60 - Maximum size allowed for compounds

  - **ntype_constraint**: 4 - Limit of different species allowed

  - **must_include**: [] - List of elements must be included in the search. Here, Here, its role differs from the ``must_include`` parameter in Materials Project's ``chemsys`` mode. For instance, in the configuration ``"element_set": "(Fe-Mn),O"`` (`here <https://github.com/mohanliu/qmpy_rester>`_), the ``entries`` consist of ``["Fe", "Mn"]``, while the ``must_include`` parameter contains ``["O"]``.

  - **metal**: false - If true, searches for compounds with zero bandgap

  - **magnetic**: true - If true, searches for compounds with magnetic properties

  - **spacegroup**: null - Spacegroup symbol to search; if null, the code doesn’t select based on spacegroup

  - **thermo_stable**: true - If true, selects compounds on the convex hull indicating thermodynamic stability

  - **FE**: true - If true, selects compounds with negative formation energy

  - **prop**: ["composition", "spacegroup", "volume", "band_gap", "stability"] - List of properties to extract for each compound

- **(E)aflow**:

  - Extract data from `AFLOW Database <http://www.aflowlib.org/>`_, utilizing parameters from ``aflow`` dictionary.

  - `The AFLUX search API <http://www.aflowlib.org/documentation/>`_

  - **elm**: ["B"] - Element to search for in compounds

  - **nelm**: 5 - Number of different types of elements in compounds

  - **nsites**: 10 - Total number of ions in the compound

  - **metal**: true - If true, selects compounds with metallic properties

  - **FE**: true - If true, selects compounds with negative formation energy

  - **spacegroup**: null - Spacegroup symbol to search; if null, the code doesn’t select based on spacegroup

  - **limit**: 5000 - Specifies the limit for the search

  - **filter**: 
    If set to True, applies filters such as ``metal``, ``FE``, ``nsites``, ``spacegroup``. Many entries in the database may not have these data, which can limit the search space. To broaden the search and include more structures, set ``filter`` to False. 
    Example: ``'filter': false``.
  
  - **prop**: 
    Specify properties to retrieve. By default, set to ``["spacegroup_relax", "Pearson_symbol_relax"]`` with ``filter`` set to False to search for additional structures. Other available properties include:
    ``["Bravais_lattice_orig", "Bravais_lattice_relax", "composition", "density", "dft_type", "eentropy_cell", "enthalpy_cell", "enthalpy_atom", "enthalpy_formation_cell", "enthalpy_formation_atom", "entropic_temperature", "kpoints", "volume_cell", "volume_atom"]``.
    Example: ``"prop": ["spacegroup_relax", "Pearson_symbol_relax"]`` with ``"filter": false``.
  
  - **Note**: 
    Setting ``filter`` to True or using more than 2 properties in ``prop``, or both, can significantly narrow down the search space. Many entries may not have these specified data, which can limit the results.
  

.. _convtest-label:

-------------------
conv_test
-------------------

.. code-block:: text

  "conv_test": {
     "param": "ecut",
     "ecut": [400, 500, 600],
     "kpoint": [[6, 6, 6], [12, 12, 12], [18, 18, 18]]
    },


- **conv_test**: This provides the parameters for the convergence tests.

  - **param**: Parameter to perform convergence tests. Available options are ``ecut`` and ``kpoint``.

  - **ecut**: List of kinetic energy cutoff (eV for VASP and Ry for QE).

  - **kpoint**: List of kpoint mesh.

  - **Note**: For ``ecut`` convergence, smallest ``kmesh`` is utilized and vice-versa.


.. _magmom-label:

--------------------
magmom
--------------------

.. code-block:: text

    "magmom": {
        "magmom": {
          "Sr": 5,
          "Al": 0,
          "I": 0},
        "type":"anisotropy",
        "saxis":[[0,0,1],[1,0,0],[1,1,0],[1,1,1]],
        "order": ["ferromagnetic", "antiferromagnetic", "ferrimagnetic_by_motif"]}

- **magmom**:

  - **magmom**: This holds a dictionary with elements as keys and initial magnetic moment in :math:`\mu_B` as values (available only for VASP).

  - **type**: This parameter is required when the magnetic enumeration (mainprogram magenum) process is run. Available options are:

    - **ordering**: This creates structures with different magnetic ordering using the ``order`` list. It uses `MagneticStructureEnumerator class <https://pymatgen.org/pymatgen.analysis.magnetism.html>`_.

    - **anisotropy**: This option changes the `SAXIS  <https://www.vasp.at/wiki/index.php/SAXIS>`_ keyword to calculate magnetic anisotropic calculations with values given by the ``saxis`` list.

    - **Note**: Don't include ``LSORBIT`` tag in :ref:`vasp.in <vasp-label>`, if you are not doing spin-orbit calculations. If you are not utilizing the ``magenum`` process, remember to configure the ``type`` setting to either ``null`` or an empty string ("").


.. _pseudo-label:

--------------------
pseudo
--------------------

.. code-block:: text

    "pseudo": {
        "pot": {
          "H":"H",
          "He":"He",
          "Li":"Li_sv",
          "Be":"Be"},
        "PSEUDO": {
          "H": 60,
          "Li": 40,
          "Be": 40,
          "N": 60,
          "F": 45}}

- **pseudo**: Pseudopotential information

  - **pot**: `Recommended pseudopotential for VASP <https://www.vasp.at/wiki/index.php/Available_PAW_potentials>`_. To use VASP POTCARs, please consult vasp developers `<https://www.vasp.at/>`_. 

    - Suppose we have POTCARS as  POT_GGA_PAW_PBE/Mg/POTCAR 
    - pmg config -p ~/path_along_directory_POT_GGA_PAW_PBE PBE52
    - After that add path to .pmgrc.yaml
    - pmg config --add PMG_VASP_PSP_DIR PBE52
    - Check your ~/.pmgrc.yaml for the path to pseudopotential, if you can find the POTCARs in ".gz" extension. 


  - **PSEUDO**: 

    - elements and energy cutoff (in Ry) as dictionary key-value pairs.
    - Checkout `QE pseudopotentials <https://www.quantum-espresso.org/pseudopotentials/>`_ and energy cutoff at `SSSP <https://www.materialscloud.org/discover/sssp/table/efficiency>`_.
    - Put your pseudopotentials in :ref:`pp <pwd-label>` folder in the working directory in "element.upf" format.


.. _substitute-label:

--------------------
substitute
--------------------

.. code-block:: text

    "substitute": {
        "mode": 2,
        "elm": "Cr",
        "sub": {"Cr": 0,"Al": 1},
        "new_sub": {"Cr": "Fe","Pd": "Pd","I": "I"}}

Parents compound needed for substitution.

- **keyword mode**: "Mode" of the substitution. Two options are available:

  - **Mode 1**: Changes the substitution on the "elm" keyword with compositions defined by the "sub" dictionary. If multiple substitutions need to be performed, simply use a list of dictionaries for different compositions. This performs substitution according to `<https://bsym.readthedocs.io/en/latest/api/interface/pymatgen.html?highlight=unique%20structures#>`_.

  - **Mode 2**: Simply replaces the elements using the "new_sub" dictionary. Here, each key is replaced by its corresponding value pair.

  - **elm**: Element to be replaced. Here "Cr" is being replaced by "Al", defined by "sub" dictionary.
  

.. _pwscfin-label:

--------------------
pwscf_in
--------------------

.. code-block:: text

    "pwscf_in": {
        "magnetic": false,
        "control": {
          "calculation": "vc-relax",
          "nstep": 300,
          "restart_mode": "from_scratch",
          "pseudo_dir": "../../pp/",
          "outdir": "./",
          "tprnfor": true,
          "tstress": true,
          "etot_conv_thr": 1e-05,
          "forc_conv_thr": 0.0001},
        "system": {
          "smearing": "gauss",
          "occupations": "smearing",
          "degauss": 0.02},
        "electrons": {
          "diagonalization": "david",
          "mixing_mode": "plain",
          "mixing_beta": 0.7,
          "conv_thr": 1e-16,
          "electron_maxstep": 300}}

- **magnetic**: if ``true``, it will assigns ``starting_magnetization`` in ``Ferromagnetic`` ordering.

- **control**: `control (QE) <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm36>`_

- **system**: `system (QE) <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm223>`_

- **electrons**: `electrons (QE) <https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm802>`_



--------------------
strain
--------------------

.. code-block:: text

    "strain": [-0.01, -0.005, 0.005, 0.01]

List of strain (both tensile (+ve) and compressive (-ve)) for deforming relaxed crystal structure to compute elastic properties.


.. _wanniertoolsinput-label:

--------------------
wanniertools_input
--------------------

.. code-block:: text

    "wanniertools_input": {
        "tb_file": {
          "Hrfile": "'ex_hr.dat'",
          "Package": "'QE'"},
        "control": {
          "BulkBand_calc": "T",
          "BulkFS_calc": "F",
          "BulkGap_cube_calc": "F",
          "BulkGap_plane_calc": "F",
          "FindNodes_calc": "F",
          "SlabBand_calc": "F",
          "WireBand_calc": "F",
          "Dos_calc": "F",
          "JDos_calc": "F",
          "SlabSS_calc": "F",
          "SlabArc_calc": "F",
          "SlabQPI_calc": "F",
          "SlabSpintexture_calc": "F",
          "wanniercenter_calc": "F",
          "Z2_3D_calc": "F",
          "Chern_3D_calc": "F",
          "WeylChirality_calc": "F",
          "BerryPhase_calc": "F",
          "BerryCurvature_calc": "F",
          "AHC_calc": "F"},
        "parameters": {
          "E_arc": "0.0",
          "Eta_Arc": "0.001",
          "OmegaMin": "0.0",
          "OmegaMax": "0.0",
          "OmegaNum": "100",
          "Nk1": "10",
          "Nk2": "10",
          "Nk3": "10",
          "NP": "2",
          "Gap_threshold": "0.1"},
        "system": {
          "NSlab": "10",
          "NSlab1": "1",
          "NSlab2": "1",
          "NumOccupied": "1",
          "SOC": "1",
          "E_FERMI": "0.0",
          "Bx": "0.0",
          "By": "0.0",
          "Bz": "0.0",
          "surf_onsite": "0.0"},
        "surface": {
          "surface": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
          "KPATH_SLAB": {},
          "KPLANE_SLAB": {},
          "EFFECTIVE_MASS": "0.0",
          "SELECTED_ATOMS": {}}}


This keyword provides input parameters for WannierTools calculations. While it has the capability to function with codes other than QE, only QE is currently implemented in our codebase. The definition of these input parameters can be found in the `WannierTools Documentation <http://www.wanniertools.com/input.html>`_.

.. _plot-label:

----------------------
plot
----------------------

It has following keys and values.

- **xlim/ylim**: List of two numbers giving the range of the plot along x or y.
  For band structures and DOS plots ``ylim`` is the *energy* window in eV
  relative to the Fermi level.

- **dos_ylim**: List of two numbers bounding the DOS axis -- the height of the
  density of states, not the energy.  Optional; when it is ``null`` the DOS axis
  is scaled to the data.  Set it when you want several compounds drawn on the
  same DOS scale.  ``ylim`` no longer has any effect on the DOS axis.

- **a2f_smearing**: One-based index of the electron-phonon smearing whose
  Eliashberg function :math:`\alpha^2F(\omega)` is plotted.  QE writes one
  ``a2F`` file per value of ``degauss`` in the el-ph calculation; ``1`` is the
  first (smallest) smearing, ``2`` the second, and so on.  Optional; when it is
  ``null`` the previous behaviour is kept and the last smearing found is used.
  Check the ``lambda.out`` file for the list of smearings and pick the one in
  the converged plateau.

- **atomproj**: Cutoff for plotting atomic projection on phonon dispersion.

- **bandproj**: Dictionary required for the projected band structure.  It has
  three keys, ``proj_type``, ``proj`` and ``colormap``.

  - **proj_type**: What ``proj`` holds.  One of:

    - **element**: a species, or a list of species -- ``"Mg"`` or
      ``["Mg", "B"]``.

    - **orbital**: an orbital, or a list of orbitals -- ``"px"``, ``"dxy"`` for
      VASP, and ``"2px"``, ``"2dxy"`` for QE.

    - **element-orbital**: a dictionary mapping species to orbitals.

      For QE: ``{"Mg": "3s", "B": "2p"}``,
      ``{"Mg": "3s", "B": ["2px", "2py", "2pz"]}``,
      ``{"Mg": "3s", "B": "2px"}``.

      For VASP: ``{"Mg": "s", "B": "p"}``,
      ``{"Mg": "s", "B": ["px", "py", "pz"]}``, ``{"Mg": "s", "B": "px"}``.

  - **proj**: the projection itself, in whichever of the three forms
    ``proj_type`` names.

  - **colormap**: colour map to use for the projection.  Sequential maps such as
    ``"Reds"`` or ``"Greens"`` read best; see the `matplotlib colormaps
    <https://matplotlib.org/stable/users/explain/colors/colormaps.html>`_ for
    the full list.

Turn ``bandproj`` on only after the band structure has been processed and
plotted with ``eband`` or ``vasp-line``.

.. _wannier90-label:

------------------------
wannier90.json
------------------------

.. code-block:: json

    {
      "config_settings": {
        "use_ws_distance": ".true.",
        "dis_num_iter": 200,
        "write_hr": ".true.",
        "iprint": 2,
        "spinors": ".false."
      },
      "plot_settings": {
        "fermi_surface_plot": ".true.",
        "bands_plot": ".true.",
        "wannier_plot": ".true.",
        "wannier_plot_supercell": 3
      }
    }

- **config_settings**: Configuration settings for `WANNIER90 <https://wannier90.readthedocs.io/en/latest/user_guide/wannier90/parameters/>`_ calculations. This dictionary contains most parameters required for ``WANNIER90`` calculations, excluding plotting.

- **plot_settings**: Parameters for plotting are presented in this dictionary.


------------------------
epw.json
------------------------

.. code-block:: json

    {
      "prefix": "MgB2",
      "outdir": "./",
      "dvscf_dir": "~/path_to_phonon/save",
      "ep_coupling": true,
      "elph": true,
      "epwwrite": true,
      "epwread": false,
      "etf_mem": 1,
      "wannierize": true,
      "nbndsub": 1,
      "bands_skipped": "exclude_bands = 1:1",
      "wdata": ["fermi_surface_plot = .true.", "dis_num_iter = 2000"],
      "max_memlt": "XXX",
      "auto_projections": true,
      "scdm_entanglement": "erfc",
      "scdm_proj": true,
      "scdm_mu": 0.123,
      "scdm_sigma": 0.456,
      "num_iter": 500,
      "dis_froz_min": "XXX",
      "dis_froz_max": "XXX",
      "dis_win_min": "XXX",
      "dis_win_max": "XXX",
      "iverbosity": 2,
      "fsthick": 0.2,
      "degaussw": 0.05,
      "ephwrite": true,
      "eliashberg": true,
      "laniso": true,
      "limag": true,
      "lpade": true,
      "lacon": false,
      "nsiter": 300,
      "conv_thr_iaxis": 1.0e-4,
      "wscut": 1.0,
      "muc": 0.16,
      "nstemp": 10,
      "temps": [10, 55],
      "specfun_el": false,
      "scattering": false,
      "scattering_0rta": false,
      "scattering_serta": false,
      "lscreen": false,
      "scr_typ": 0,
      "phonselfen": false,
      "a2f": false,
      "nest_fn": false,
      "lpolar": false,
      "longrange": false,
      "elecselfen": false
    }

This keyword provides input parameters for EPW calculations. The definition of these input parameters can be found in the `EPW Documentation <https://docs.epw-code.org/doc/Inputs.html>`_.
