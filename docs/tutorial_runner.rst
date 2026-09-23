.. _tutorial-runner-label:

=====================
The tutorial runner
=====================

``examples/`` ships 42 worked tutorials -- ``examples/QE/`` and
``examples/VASP/``, ``tutorial1`` through ``tutorial21`` in each.  Each is a prose ``README`` naming the
``mainprogram`` commands to run, in order, and the files to copy in from an
earlier tutorial first.  ``tutorials/`` turns that prose into a catalogue and
runs it: one command walks the whole tree, waits for the cluster jobs each step
submits, checkpoints after every step, and when it stops says exactly where.

--------------------
Running it
--------------------

.. code-block:: bash

    htesp-tutorials                           # the installed entry point
    python -m tutorials.run_tutorials         # the same thing from a checkout
    sbatch tutorials/submit_tutorials.sh --only QE

--------------------
Where ``examples/`` is found
--------------------

``examples/`` is 185 MB of reference data and is **not** shipped inside the
wheel, so after an ordinary ``pip install .`` there is no ``examples/`` beside
the installed package -- which is what ``the example tree
.../site-packages/examples is missing`` means.  The runner looks for the tree
in this order:

1. ``$HTESP_EXAMPLES``
2. ``./examples`` under the current directory
3. beside the package -- a source checkout, or ``pip install -e .``

in that order:

.. code-block:: bash

    HTESP_EXAMPLES=~/HTESP_claude/examples htesp-tutorials --workdir ~/tutorial_runs
    export HTESP_EXAMPLES=~/HTESP_claude/examples      # or set it once

``examples/`` is read-only input: nothing is ever written there.  ``--workdir``
is where runs, logs and the checkpoint go, and a ``--workdir`` inside the
example tree is refused before anything is created.

--------------------
The three modes
--------------------

There are no modes.  Every step is invoked as ``mainprogram <cmd> --dry-run``:
all the input generation and file plumbing, nothing submitted, nothing
deleted.  A laptop is enough.  Driving a real campaign is ``mainprogram``'s
job, not this driver's -- see :doc:`usage`.

Steps that read the **output** of a real DFT run -- ``phono2``..``phono4``,
``ev-collect``, ``compute-elastic``, the plotting steps -- are recorded as
*skipped*, with that reason and with the tutorial's own instructions, rather
than failing for a reason that is not their own.  The steps that need only a
*relaxation* do run: its output is seeded from the reference (see
:data:`~tutorials.catalog.REFERENCE_OUTPUT`).

--------------------
Selecting tutorials
--------------------

.. code-block:: bash

    htesp-tutorials --list                      # the catalogue, every step
    htesp-tutorials --only QE                   # one whole example tree
    htesp-tutorials --only QE/9,QE/12           # just these
    htesp-tutorials --skip QE/4,QE/5            # all but these
    htesp-tutorials --from band-scf             # start each one at this step

Other options: ``--workdir DIR`` (root for work directories, logs, checkpoint
and report), ``--workers N`` (passed through to ``mainprogram --workers``),
``--timeout HOURS`` (how long to wait for a step's cluster jobs, default 24).
The poll interval (60 s) and the limit on one ``mainprogram`` call (6 h) are
fixed, because neither depends on the study,
``--force`` (run even when preflight found errors) and ``-v``.

Exit codes: ``0`` all finished, ``1`` something failed or was blocked, ``2`` a
preflight check failed before anything ran, ``130`` interrupted.

------------------------------------
examples/ is never written to
------------------------------------

Each tutorial runs in ``<workdir>/tutorial_runs/<CODE>/``, seeded in this order:

1. ``examples/<QE|VASP>/`` -- the shared ``batch.header``, ``config.json``,
   ``input.in`` / ``vasp.in``, and for QE a symlink to the pseudopotentials in
   ``examples/QE/pp/``, so 42 work directories do not each cost a copy of them;
2. whatever the tutorials this one **depends on** produced -- ``R<mpid>-<name>/``,
   ``scf_dir/``, ``mpid.in``, ``econv.csv``.  Tutorial 9 (relaxation) is the hub
   most later tutorials start from;
3. the tutorial's own shipped files, which therefore win.

``README``, ``log`` and every ``reference*.tar.gz`` are skipped: those are the
expected answer, not input.  The one exception is the ``.cif`` tutorial, whose
inputs genuinely live inside its reference archive and are unpacked from it.

--------------------
Resuming
--------------------

Every step writes ``<workdir>/state.json`` the moment it finishes: its status
(``pending``, ``running``, ``done``, ``failed``, ``skipped``, ``blocked``), exit
code, duration, log path, the artefacts expected and which were missing, and the
job ids it submitted.

``--resume`` is the default and re-runs only what is not ``done`` or
``skipped``.  ``--restart`` forgets the file and starts over.  An interrupted
run is saved too, so Ctrl-C is a pause, not a loss.

--------------------
The stop report
--------------------

On any stop -- a non-zero exit, a step that exited 0 and produced nothing, a
timeout, a dependency that blocked a tutorial, or Ctrl-C -- the driver prints a
summary and writes ``<workdir>/report.md`` and ``<workdir>/report.json``
containing:

* a one-screen table of all 42 tutorials: done / failed / blocked / skipped,
  with step counts;
* **where it stopped**: the tutorial code and title, the step number and id, the
  exact command line, the working directory, the exit code, the last 40 lines of
  that step's log, the artefacts expected and which were missing, and how long
  it ran;
* the exact retry line, for example
  ``htesp-tutorials --resume --only QE/12 --from band-scf --workdir ...``;
* everything that was **blocked**, and by which dependency;
* any step that finished but could not be **verified** -- no ``squeue`` on
  ``$PATH``, or no job ids recorded.  Those are called out, not counted as
  success.

Logs are at ``<workdir>/logs/<CODE>/<NN>-<step-id>.log``, one per step, with the
command and working directory in the header.

--------------------
How a step is judged
--------------------

A step is ``done`` when the command exits ``0`` **and** the artefacts the
catalogue declares for it exist.  "Exited 0 and wrote nothing" is this package's
most common silent failure, so it is treated as a failure and the report names
the glob that matched nothing.

After a step that submits, the runner reads the job ids the
:doc:`workflow layer <workflow>` recorded in ``<stage dir>/.htesp_job.json`` and
polls ``squeue -h -o %i -j <ids>`` until none remain.  The queue is never grepped
for a compound name.

--------------------
Prerequisites
--------------------

``preflight()`` reports every problem before the first step runs:

* ``python -m htesp`` must start;
* ``examples/<code>/batch.header`` must exist, and QE needs
  ``examples/QE/pp/*.upf``;
* a **Materials Project API key** is needed by the database tutorials (QE/2,
  QE/3, QE/6, QE/16 and the VASP equivalents), since it was removed from every
  shipped ``config.json``.  It is resolved the way HTESP resolves it
  everywhere -- ``$MP_API_KEY``, then ``~/.config/htesp/credentials``, then
  ``config.json`` -- so ``htesp-check --set_mp_api`` is enough and nothing has
  to be exported.  Its absence is a warning, not an error: the four tutorials
  that need it are skipped and the rest proceed.

Nothing else is required: no scheduler, no ``pw.x``, no ``vasp_std``.

Beyond that, the optional dependencies of the tutorials that use them must be
installed -- ``ase``, ``pymatgen``, ``bsym``, ``phonopy``, ``qmpy-rester``,
``ifermi``.  A missing one shows up in the report as the traceback of the step
that needed it.

--------------------
Known limitations
--------------------

* VASP needs ``POTCAR`` files, which are licensed and are not shipped; configure
  pymatgen's ``PMG_VASP_PSP_DIR`` before running the VASP tutorials for real.
* ``examples/VASP/tutorial21`` (3D Fermi surface) ships only ``ifermi.tar.gz`` --
  no ``config.json``, no ``input.in``, no ``vasprun.xml``.  It is flagged as a
  stub in the catalogue and cannot run as shipped; preflight warns about it
  and the report names it, so leave it out with ``--skip`` if it is in the
  way.
* The relaxation loop (``2`` → ``3`` → ``e0`` until ``niteration < 3``) runs at
  most four cycles by default.  A system that has not converged by then is
  reported, not looped forever.

--------------------
Self-tests
--------------------

.. code-block:: bash

    python -m unittest tutorials.selftest -v
    pytest tutorials/selftest.py

Thirty tests, standard library only, no ``pymatgen`` and no scheduler.  They
cover catalogue integrity, the checkpoint round trip, ``input.in`` patching,
artefact verification, job-id collection, the missing-``squeue`` path, the
contents of a stop report, and a full end-to-end run against a stub
``mainprogram``.
The batch header each tutorial runs with
-----------------------------------------

``examples/<code>/batch.header`` says ``--partition=dense`` and loads no
module.  It was written for one machine, so everywhere else ``sbatch`` rejects
it before any calculation starts.  When ``sinfo`` is present, seeding therefore
overwrites the copy in each work directory with one built by
:mod:`htesp.batch_header` -- this cluster's partition, account, cores per node
and ``qe``/``vasp`` module -- and points ``job_script.parallel_command`` at
whichever of ``ibrun``, ``srun``, ``mpirun`` or ``mpiexec`` is on ``$PATH``.
``ibrun`` exists only at TACC; elsewhere this lands on ``srun`` or ``mpirun``.

``job_script.nproc`` is never touched: how many ranks a study wants is not
something the machine can answer.

The ``module load`` line names the module **without a version** -- ``module
load qe``, not ``module load qe/7.3`` -- so Lmod resolves it to whatever the
site has marked default and the header keeps working after 7.3 is retired.
The versions found are written underneath as a comment, so a study that needs
one exact build can still pin it by writing the version out.

Building the submission scripts
--------------------------------

``HTESPWorkflow.stage_and_submit`` copies ``run-<stage>.sh`` from the work
directory into the stage directory and submits *that* -- ``run-scf.sh`` for
Quantum ESPRESSO, ``run-vasp.sh`` for VASP, both built from ``batch.header``
and ``job_script.command_list``.  When the script is missing it records the
material as *skipped* and carries on, so a tutorial that never ran
``mainprogram jobscript`` wrote no run-*.sh and the step still exited 0.

Every tutorial with a submitting step therefore gets a ``jobscript`` step
prepended (:func:`tutorials.catalog._with_job_scripts`), and a submitting step
that records no job id is now a failure rather than an "unverifiable": whatever
the cause, nothing was submitted.

On a machine with no ``sinfo`` nothing is generated -- the probes would have
nothing to say, and this runner submits nothing anyway.
What the run produced
---------------------

The status table answers "did it work?".  ``--output`` answers "what did it
give me?" -- every file the run wrote, grouped by the step that wrote it, with
a line saying what that file is for:

.. code-block:: text

    QE/9  Structural relaxation (the hub every later tutorial starts from) (QE)
        jobscript          run-scf.sh
                               -- submission script: batch.header with one
                                  stage's run command appended
        relax-submit       Rmp-763-Mg1B2/relax/scf.in
                               -- Quantum ESPRESSO relaxation input
        resubmit           (wrote nothing but the log)

Every step records the paths that appeared or changed while it ran, so this is
what the step *did*, not what its artefact globs merely matched.  The
difference is not academic: three defects were found exactly this way -- two
steps declaring artefacts they never write (``pressure-input`` claimed a
directory that ``mainprogram 26`` creates later, ``update-input`` one written
only when a structure is *not* relaxed), and a tutorial whose declared
artefact is shipped inside ``examples/`` itself, so the glob matched whether
or not the step ran.

``log`` is excluded from the judgement, though still listed: every
``mainprogram`` call appends to it, so counting it would mean no step ever
looks empty -- and "this step passed without writing anything" is the signal
worth having.  ``report.md`` carries the same table, plus a list of exactly
those steps.

Symlinked directories are not walked.  QE work directories link ``pp/`` at the
shared pseudopotential tree, which no step writes to.
Steps this runner cannot perform
--------------------------------

Forty-one steps read the output of a real Quantum ESPRESSO or VASP run --
FORCE_SETS, ``ph.x`` output, computed bands, the deformed cells of an elastic
calculation, relaxations across a whole material set.  Nothing here produces
those, so they are recorded as skipped **with the tutorial's own written
instructions**::

    3. energy  SKIPPED (reads the output of a real DFT run, which this
       runner never performs; to run it for real, follow
       examples/QE/tutorial9/README)

Twelve tutorials ship no ``README`` of their own -- QE/15 and eleven VASP
ones.  The two trees cover the same topics in the same order, so the
counterpart's instructions are used: ``VASP/9`` has none, and ``QE/9``
describes the same relaxation.  The numbering diverges from 11 onward, which
:func:`~tutorials.catalog.vasp_number_to_qe_number` already accounts for, so
``VASP/19`` resolves to ``QE/20``.  ``README.txt`` is accepted as well as
``README``.

Only VASP/21 has instructions nowhere -- its topic, IFermi Fermi surfaces, has
no QE counterpart -- and it gets no pointer rather than a path that does not
exist.

The pointer is added to **that skip reason only**.  A step skipped because the
machine has no POTCAR, no API key or no enumlib already says what to install,
and "read the README" would be a downgrade; a step skipped because the step
before it was skipped names *that* step, which is the actual cause.

``report.md`` closes with the tutorials nothing ran for at all -- the two
convex-hull tutorials, which need relaxations for a whole material set rather
than the single reference material, and VASP/21 -- each with the same
pointer.
