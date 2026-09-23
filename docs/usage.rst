.. _usage-label:

----------------------------
Package structure
----------------------------

.. code-block:: text

    HTESP/
    ├── htesp/          the importable package (every module, plus data/config.json)
    ├── bin/            52 compatibility shims with the old scan-script names
    ├── legacy/bash/    the original bash scan scripts, unmodified, for reference
    ├── tutorials/      the tutorial runner (htesp-tutorials)
    ├── tests/          the test suite
    ├── tools/          check_names.py, gen_command_rst.py
    ├── examples/       42 worked tutorials, QE/tutorial1..21 and VASP/tutorial1..21
    ├── utility/        input-file templates and stand-alone analysis scripts
    ├── docs/           this documentation
    ├── pyproject.toml  package metadata and dependencies
    ├── LICENSE
    └── README.md

``htesp/workflow.py`` holds the workflow layer that replaced the bash scan
scripts; see :doc:`workflow`.

----------------------------
Requirements
----------------------------

Linux and macOS, Python 3.10 or newer.  Quantum ESPRESSO and/or VASP and a SLURM
scheduler are needed to run calculations, but not to prepare inputs or to read
results.

The required Python packages are declared in ``pyproject.toml`` and are
installed for you by ``pip install .``:

``numpy``, ``scipy``, ``pandas``, ``matplotlib``, ``pymatgen``,
``mp-api`` (imported as ``mp_api``),
`ASE <https://wiki.fysik.dtu.dk/ase/>`_, ``spglib``, ``PyYAML``,
`bsym <https://bsym.readthedocs.io/>`_ (site substitutions),
`lmfit <https://lmfit.github.io/lmfit-py/>`_ (the SCDM fit that produces initial
projections for wannierisation) and
`qmpy-rester <https://github.com/mohanliu/qmpy_rester>`_ (the OQMD searches;
imported as ``qmpy_rester``).

Optional features live in extras:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Extra
     - What it adds
   * - ``htesp[ml]``
     - ``matminer`` and ``scikit-learn``, for the machine-learning helpers
   * - ``htesp[fermisurface]``
     - `IFermi <https://fermisurfaces.github.io/IFermi/>`_ and ``plotly``, for
       ``mainprogram fermisurface``
   * - ``htesp[test]``
     - ``pytest`` and ``pytest-cov``
   * - ``htesp[docs]``
     - ``sphinx`` and the theme used to build this documentation
   * - ``htesp[all]``
     - ``ml`` and ``fermisurface`` together

`Phonopy <https://phonopy.github.io/phonopy/>`_ is not a Python dependency of
HTESP: install it separately (``conda install -c conda-forge phonopy``) and make
sure the ``phonopy`` command is on ``$PATH`` before using ``mainprogram phono1``
and its siblings.  `cif2cell <https://pypi.org/project/cif2cell/>`_ is likewise
optional and only needed when ``download.inp.use_cif2cell`` is ``true``.

``requirements.txt`` and ``INSTALL/requirements1.txt`` are kept for people who
expect them, and they are not the source of truth -- ``pyproject.toml`` is.
They also do different jobs, which is worth knowing before either one surprises
you:

* ``requirements.txt`` pins the **required** set.  The extras are deliberately
  absent, so ``pip install -r requirements.txt`` leaves ``matminer``,
  ``scikit-learn``, ``ifermi`` and ``plotly`` uninstalled.
* ``INSTALL/requirements1.txt`` is the same pinned set plus
  ``htesp[fermisurface]``, and carries the reasoning behind the versions.

Both also pin ``emmet-core``, which no extra declares: ``mp-api`` 0.46.0 asks
for ``emmet-core>=0.86.3`` with no upper bound but imports a name that 0.87
moved, so an unpinned resolve produces an installation whose
``import mp_api.client`` fails.  A test keeps the versions in the two files
identical.

If something is missing, ``htesp-check`` names the extra that provides it.

----------------------------
Download software
----------------------------

.. code-block:: bash

    git clone https://github.com/Neraaz/HTESP.git

or go to `GitHub <https://github.com/Neraaz/HTESP>`_, download a zip file under
the code section and unzip it.

Go to the directory:

.. code-block:: bash

    cd HTESP

----------------------------
Conda environment
----------------------------

Make sure conda is available, either via
`miniconda <https://docs.anaconda.com/free/miniconda/>`_ or
`anaconda <https://www.anaconda.com/download/success>`_.

.. code-block:: bash

    conda create --name htesp python=3.11
    conda activate htesp

Any Python from 3.10 onwards works; 3.11 is what the package is developed
against.

----------------------------
Install HTESP
----------------------------

.. code-block:: bash

    pip install .

That is the whole installation: it installs the ``htesp`` package, every
required dependency, the ``mainprogram`` entry point and the ``bin/`` shims.
Add the extras you need:

.. code-block:: bash

    pip install ".[ml,fermisurface]"

For an editable install with the test and documentation tooling:

.. code-block:: bash

    pip install -e ".[test,docs]"

There is no ``PYTHONPATH`` step.  ``mainprogram`` imports ``htesp``, which pip
has put on the path; setting ``PYTHONPATH`` to a source directory now only
creates a second, shadowing copy.

Check the executable:

.. code-block:: bash

    which mainprogram
    mainprogram --version
    mainprogram basicinfo

``mainprogram --list`` prints every process number and command name with a
one-line description; the same tables are in the :ref:`command reference
<command-label>`.

.. _cli-label:

----------------------------
Command-line reference
----------------------------

Five commands are installed.  ``mainprogram`` runs the science; the others
support it.

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Command
     - Purpose
   * - ``mainprogram``
     - the campaign driver (also installed as ``htesp``)
   * - ``htesp-check``
     - environment report, and the one-off configuration steps
   * - ``htesp-tutorials``
     - run the 42 worked examples end to end
   * - ``htesp-workflow``
     - invoke one former scan script directly
   * - ``bin/<name>``
     - 52 compatibility shims, one per 1.x bash script

mainprogram
===========

.. code-block:: text

    mainprogram <process> [options]

``<process>`` is a number (0-29) or one of 58 named commands.  The
program prints its own list, which is generated from the same text as
:doc:`command`, so the two cannot disagree:

.. code-block:: bash

    mainprogram --list          # every process and command, one line each
    mainprogram basicinfo       # the introduction
    mainprogram process-info    # the numbered processes
    mainprogram epw-info        # the EPW / Wannier90 pipeline
    mainprogram wt-info         # the WannierTools pipeline

Global options:

.. code-block:: text

     --workers WORKERS  per-material process pool size (default: min(cpu, 8))
     --dry-run          prepare every input file but do not submit anything
     --root ROOT        project directory (default: .)
     --config CONFIG    config.json to use instead of searching for one
     --force            overwrite a file the command would otherwise refuse to
                        replace (config-init)
     -v, --verbose      debug logging
     --list             list every command and exit
     --version
     -h, --help

**Exit status** is ``0`` on success, ``1`` when one or more materials failed,
``2`` for a bad command or a malformed ``input.in``, and ``130`` on interrupt.
Check it in batch scripts: in 1.x every exit code was discarded, so a stage in
which every material crashed still looked successful and the next stage ran on
nothing.

htesp-check
===========

Reports what this machine can actually run -- interpreter, architecture, memory
page size, and every dependency imported in a *child* interpreter, so a wheel
built for the wrong architecture is reported rather than taking the report down
with it.  It also carries the one-off configuration steps.

.. code-block:: text

     -h, --help            show this help message and exit
     --json                print the report as JSON
     --no-extras           probe only the required dependencies
     --executables         also look for pw.x, vasp_std, sbatch, ...
     --config_vasp_pot DIR
                           point pymatgen at a VASP POTCAR tree and exit; give
                           either POT_GGA_PAW_PBE or its parent
     --clean               remove build artifacts (build/, *.egg-info/,
                           __pycache__/, *.pyc) from the source checkout,
                           returning it to its pre-build state; leaves enumlib
                           and all configuration alone
     --root DIR            the checkout --clean acts on (default: the one this
                           htesp package lives in)
     --set_mp_api KEY      write the Materials Project API key to
                           ~/.config/htesp/credentials and verify it; survives
                           new shells and batch jobs, unlike an exported
                           MP_API_KEY
     --dry-run             with --clean, list what would be removed and remove
                           nothing
     --no-verify           with --set_mp_api, skip the live check
     --install-enumlib     build enumlib from source and install enum.x and
                           makestr.x (needed by 'mainprogram magenum'); it is not
                           on PyPI, so pip cannot install it
     --prefix DIR          where --install-enumlib puts the executables (default:
                           the bin/ of this Python environment)
     --fortran-compiler FC
                           compiler for --install-enumlib (default: gfortran)
     --install-phonopy     install phonopy, which the eleven 'mainprogram
                           phono*' commands shell out to; conda-forge by
                           default, as INSTALL/README specifies
     --installer {auto,conda,mamba,pip}
                           how --install-phonopy installs (default: auto --
                           mamba, then conda, then pip)
     --phonopy-version V   exact phonopy version for --install-phonopy (default:
                           newest the channel offers)

``--install-phonopy`` prefers conda-forge because that channel ships phonopy
prebuilt; ``-p`` pins the install to the interpreter running ``htesp-check``,
so it cannot land in the ``base`` environment where HTESP would not see it.  On
a non-conda interpreter it falls back to pip, where phonopy builds from source
and needs a C compiler.  An existing ``phonopy`` on ``$PATH`` is left alone.
Afterwards the executable is run once (``phonopy --help``; phonopy has no
``--version`` flag) so that an installer exiting ``0`` without putting anything
on ``$PATH`` is reported rather than believed.

Exit status is ``0`` only when every required dependency imports.

htesp-tutorials
===============

Runs the worked examples and, when it stops, says exactly where.  See
:doc:`tutorial_runner` for the full description.

.. code-block:: text

     -h, --help            show this help message and exit
     --workdir WORKDIR     root for the work directories, logs, checkpoint and
                           report (default: ./tutorial_runs_root)
     --resume              re-run only what is not already done (the default)
     --restart             forget the checkpoint and run everything again
     --keep_output {yes,no}
                           keep the generated tutorial_runs/ directories when the
                           run ends (default: yes). 'no' removes them, including
                           those of failed tutorials, so keep 'yes' while
                           debugging. Logs, the report and the checkpoint are
                           kept either way
     --only ONLY           comma-separated tutorial codes to run, e.g.
                           QE/9,VASP/14. A bare tree name means all of it: --only
                           QE
     --skip SKIP           comma-separated tutorial codes to leave out; a bare
                           tree name works here too
     --from FROM_STEP      start each selected tutorial at this step id (the
                           retry line in the report uses this)
     --workers WORKERS     passed through to mainprogram --workers (default: 1).
                           A tutorial works on one or two materials, so a bigger
                           pool buys nothing and multiplies with --jobs: four
                           tutorials at mainprogram's own default of 8 is 36
                           processes, and a login node allows 100
     --jobs N              run N tutorials at once (default: 1). Most of a sweep
                           is spent waiting on the Materials Project, OQMD and
                           AFLOW servers, and those tutorials are independent, so
                           the waiting overlaps. Keep it modest: the same APIs
                           rate-limit
     --force               run even when preflight reports errors
     --output              after the run, list every file it produced and what
                           that file is for, grouped by the step that wrote it. A
                           step that wrote nothing is named as such -- which is
                           how two wrong artefact declarations were found
     --list                print the catalogue and exit
     -v, --verbose         debug logging

Exit codes: ``0`` all finished, ``1`` something failed or was blocked, ``2`` a
preflight check failed before anything ran, ``130`` interrupted.

htesp-workflow and the bin/ shims
=================================

Every former ``src/bash`` scan script is a method of ``HTESPWorkflow``.  These
two call the same code:

.. code-block:: bash

    htesp-workflow relax-scan 1 5 mpid.in --workers 8
    relax-scan 1 5 mpid.in --workers 8      # the bin/ shim, with bin/ on PATH

The argument contract is unchanged from 1.x -- ``start end trackfile [extra]``,
with ``end`` **exclusive** -- and both accept ``--workers N``, ``--dry-run``,
``--keep-scratch``, ``--root DIR`` and ``-v``.  51 script names are
recognised; ``htesp-workflow --help`` lists them.

.. _environment-label:

Environment variables
=====================

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Variable
     - Effect
   * - ``MP_API_KEY``
     - Materials Project key.  Takes precedence over
       ``~/.config/htesp/credentials``; see `Credentials`_ below.
   * - ``HTESP_CONFIG``
     - A ``config.json`` to use instead of searching the working directory and
       its parents.
   * - ``HTESP_WORKERS``
     - Default size of the per-material process pool.  ``--workers`` overrides
       it; ``1`` disables multiprocessing.
   * - ``HTESP_EXAMPLES``
     - Where ``htesp-tutorials`` looks for the example tree, before
       ``./examples`` and the package-relative path.
   * - ``HTESP_PYTHON``
     - The interpreter the ``bin/`` shims exec (default ``python3``).
   * - ``PMG_VASP_PSP_DIR``
     - pymatgen's POTCAR location.  Set it with
       ``htesp-check --config_vasp_pot``; nothing that writes VASP inputs works
       without it.

.. _credentials-label:

Credentials
===========

The Materials Project key is read from, in order:

1. ``$MP_API_KEY``
2. ``~/.config/htesp/credentials`` (a ``MP_API_KEY=<key>`` line, mode ``0600``)
3. ``config.json`` -- **supported but discouraged**

Prefer the credentials file: an exported variable lives only in the shell that
exported it, so a batch job, a ``nohup``-ed sweep or a fresh terminal loses it
and every database command starts skipping or failing.  Write it once with

.. code-block:: bash

    htesp-check --set_mp_api <your key>

which verifies the key against the API *before* storing it, so a typo never
replaces a working key.

Do not put the key in ``config.json``.  In 1.x a real key was committed to 218
tracked files; every shipped configuration now carries the placeholder
``use_your_API_KEY`` instead.  If you used 1.x, rotate that key -- it is still
in the git history of the original repository.

----------------------------
Reinstalling from scratch
----------------------------

``pip uninstall htesp`` removes the package, the five console scripts
(``mainprogram``, ``htesp``, ``htesp-workflow``, ``htesp-tutorials``,
``htesp-check``) and the 52 shims in ``bin/``.  From the repository root:

.. code-block:: bash

    pip uninstall -y htesp
    pip install .                     # or: pip install -e ".[test,docs]"
    htesp-check                      # confirm what is now installed

**If HTESP 1.x was ever installed in the same environment**, uninstalling
``htesp`` is not enough.  1.x was a different distribution, named ``HTESP``, and
its ``find_packages()`` installed a top-level package literally called ``src``
along with 46 modules as bare commands on ``$PATH`` (``plot.py``, ``band.py``,
``dos.py``, ``crystal.py``, ...).  pip treats the two as unrelated, so they sit
side by side and whichever comes first on ``sys.path`` wins:

.. code-block:: bash

    pip uninstall -y htesp HTESP
    python -c "import src" 2>&1 | tail -1   # should say No module named 'src'
    pip install .

``htesp-check`` reports which distribution is installed, which directory the
running ``htesp`` package was imported from, and warns when both are present or
when a stray ``src`` package is importable.

Two shell settings from the 1.x instructions are also worth removing from
``~/.bashrc``, because they put the old tree back on the path:
``export PYTHONPATH="path_to_HTESP/src"`` and any ``PATH`` entry ending in
``/src/bash`` (the replacement is ``path_to_HTESP/bin``).

An editable install (``pip install -e .``) points at the source tree, so after
``git pull`` or an ``rsync`` there is nothing to reinstall -- but the console
scripts are only rewritten by a real install, so re-run ``pip install -e .``
after changing ``[project.scripts]`` or adding a shim to ``bin/``.

----------------------------
After installation
----------------------------

**Check the machine first.**  HTESP itself is pure Python and behaves
identically on x86_64 and on arm64/aarch64: there is no compiled code in
``htesp/``, no ``platform.machine()`` branch and every shim in ``bin/`` is
POSIX ``sh``.  What differs between architectures is the compiled wheels
underneath -- NumPy, SciPy, pymatgen, spglib, pyarrow and ``deltalake``, which
``mp_api`` imports.  A wheel built for the wrong architecture does not always
raise ``ImportError``; it can abort the interpreter with ``SIGABRT``, which no
``try``/``except`` can catch:

.. code-block:: bash

    htesp-check                 # required dependencies and the extras
    htesp-check --executables   # also look for pw.x, vasp_std, sbatch
    htesp-check --json          # the same report as data

The Materials Project key can be stored once instead of exported in every
shell:

.. code-block:: bash

    htesp-check --set_mp_api <your key>

It is verified against the API before anything is written -- a rejected key
leaves the existing credentials untouched -- and stored ``0600`` in
``~/.config/htesp/credentials``, which ``htesp/config.py`` reads.  This matters
for batch jobs and for ``htesp-tutorials``: an exported ``$MP_API_KEY`` is lost
by any new shell, and the database tutorials then skip themselves.  ``$MP_API_KEY``
still takes precedence when it is set, and the command says so if the two differ.

VASP POTCARs are licensed and are not shipped with anything, so pymatgen has to
be told where yours are or every command that writes a VASP input fails with
``PmgVaspPspDirError: PMG_VASP_PSP_DIR is not set``:

.. code-block:: bash

    htesp-check --config_vasp_pot /path/to/POT_GGA_PAW_PBE

Give it either the ``POT_GGA_PAW_PBE`` directory or the directory above it --
pymatgen's setting names the *parent*, which is the easy thing to get wrong.
The option writes ``PMG_VASP_PSP_DIR`` into pymatgen's configuration (the same
file ``pmg config --add`` uses) and then proves it by writing a POTCAR, so a
tree in the wrong layout is reported immediately rather than at the first
material.

Each dependency is imported in a child interpreter, so one that aborts is
reported as ``ABORTED  <name>  killed by SIGABRT`` instead of taking the check
down with it, and the report says which kind of failure it was: a wheel for the
wrong architecture, a wheel that uses instructions this CPU lacks, or -- the
case seen on a 64 KiB-page aarch64 node -- a wheel that is right for the
architecture but whose bundled allocator was built assuming 4 KiB memory pages.
Only the first of those is fixed by reinstalling; the report gives the command
for each. The exit status is ``0`` only when every required dependency imports.

**Put the scan-script shims on ``$PATH``.**  ``bin/`` holds one shim per former
bash scan script, with the same name and the same argument contract
(``start end trackfile [extra]``, ``end`` exclusive).  Add to ``~/.bashrc``:

.. code-block:: bash

    export PATH="path_to_HTESP/bin:$PATH"

This replaces the old ``path_to_HTESP/src/bash`` entry.  Every habit and every
script that called ``relax-scan 1 5 mpid.in`` keeps working;
``bin/jobscript.sh`` is still sourceable.  The originals are kept unmodified in
``legacy/bash/`` for reference and are not on ``$PATH``.

**Give the project its own configuration.**  With no ``config.json`` in the
working directory (or up to five parents above it, or ``$HTESP_CONFIG``), HTESP
runs on the packaged default: every key has a value, so nothing fails, but the
cutoffs and k-point density are whatever the package ships rather than what the
study needs.  ``mainprogram config-validate`` reports which file is in force,
and ``config-init`` writes one to edit:

.. code-block:: bash

    mainprogram config-init          # writes ./config.json, refuses to clobber
    mainprogram config-init --force  # ... unless told to
    mainprogram config-validate      # names the file and checks it

Every run also records the resolved path in its ``log``, so which configuration
a finished campaign used is answerable afterwards.

**Set the Materials Project API key.**  It is read from the environment, not
from ``config.json``:

.. code-block:: bash

    export MP_API_KEY=your_key_here

Get a key at https://next-gen.materialsproject.org/api#api-key.  To avoid
putting it in your shell history, write it to a credentials file instead:

.. code-block:: bash

    mkdir -p ~/.config/htesp
    printf 'MP_API_KEY=your_key_here\n' > ~/.config/htesp/credentials
    chmod 600 ~/.config/htesp/credentials

The OQMD and AFLOW searches need no key.

**Make the DFT executables reachable.**  The directory holding ``pw.x`` and the
other QE binaries, or ``vasp_std``, must be on ``$PATH`` inside the job script;
put the ``module load`` lines in ``batch.header``.

**Configure the VASP pseudopotentials** with pymatgen's ``pmg`` command; see the
:ref:`pseudo <pseudo-label>` section.  POTCAR files are licensed and are not
shipped with HTESP.

----------------------------
Global options
----------------------------

These may be added to any ``mainprogram`` command:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Option
     - Meaning
   * - ``--workers N``
     - Size of the per-material process pool.  The default is
       ``min(cpu_count, 8)``; ``$HTESP_WORKERS`` changes that default, and
       ``--workers 1`` disables multiprocessing.
   * - ``--dry-run``
     - Build every input file but never call the scheduler.  Nothing is
       submitted and nothing is deleted.
   * - ``--root DIR``
     - Run against ``DIR`` instead of the working directory.
   * - ``--config FILE``
     - Use ``FILE`` instead of searching for ``config.json``.  Same as
       ``$HTESP_CONFIG``.
   * - ``--init-header CODE``
     - Only with ``jobscript``: write a starting ``batch.header`` for that
       code, filled in from what SLURM and Lmod report on this machine
       (partition, account, cores per node, module, MPI launcher).  ``CODE`` is
       ``vasp``, or any spelling of Quantum ESPRESSO -- ``qe``, ``QE``,
       ``QuantumEspresso``, ``quantum-espresso``.  Refuses to overwrite an
       existing header unless ``--force`` is given.
   * - ``--force``
     - Let a command overwrite a file it would otherwise refuse to replace
       (``config-init``, ``jobscript --init-header``).
   * - ``-v``, ``--verbose``
     - Debug logging.

and, on their own:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Option
     - Meaning
   * - ``--list``
     - Print every process number and command name with a one-line description.
   * - ``--version``
     - Print the installed version.
   * - ``-h``, ``--help``
     - Print the usage.

----------------------------
Exit codes
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 12 88

   * - Code
     - Meaning
   * - ``0``
     - Every material in the range succeeded.
   * - ``1``
     - One or more materials failed.  The failures are named on stderr.
   * - ``2``
     - The command does not exist, or ``input.in`` could not be parsed.
   * - ``130``
     - Interrupted (Ctrl-C).

A stage in which every material failed exits non-zero, so ``&&`` chains and job
scripts stop where they should.

----------------------------
Check basic commands
----------------------------

.. code-block:: bash

    mainprogram basicinfo       # the introduction
    mainprogram process-info    # the numbered processes
    mainprogram epw-info        # the EPW / Wannier90 pipeline
    mainprogram wt-info         # the WannierTools pipeline
    mainprogram config-validate # check config.json before a campaign
