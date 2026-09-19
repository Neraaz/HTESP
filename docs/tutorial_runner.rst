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

    htesp-tutorials --dry-run                 # the installed entry point
    python -m tutorials.run_tutorials --dry-run   # the same thing from a checkout
    sbatch tutorials/submit_tutorials.sh --code QE

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

and ``--examples DIR`` overrides all three:

.. code-block:: bash

    htesp-tutorials --examples ~/HTESP_claude/examples --workdir ~/tutorial_runs
    export HTESP_EXAMPLES=~/HTESP_claude/examples      # or set it once

``examples/`` is read-only input: nothing is ever written there.  ``--workdir``
is where runs, logs and the checkpoint go, and a ``--workdir`` inside the
example tree is refused before anything is created.

--------------------
The three modes
--------------------

.. list-table::
   :header-rows: 1
   :widths: 16 54 30

   * - Mode
     - What it does
     - What it needs
   * - ``--dry-run``
     - Runs every step as ``mainprogram <cmd> --dry-run``: all input generation
       and file plumbing, nothing submitted.
     - A laptop.
   * - ``--no-dft``
     - Prepares everything and reports what *would* be submitted at each
       submission point, but submits nothing.
     - A laptop.
   * - *(default)*
     - The real campaign: submits jobs, polls ``squeue``, waits, verifies.
     - QE or VASP, and SLURM.

Steps that read the **output** of a real DFT run -- ``e0``, ``2``, ``4``,
``phono2``..``phono4``, ``ev-collect``, ``compute-elastic``, the plotting steps
-- are recorded as *skipped*, with that reason, under ``--dry-run``, rather than
failing for a reason that is not their own.

--------------------
Selecting tutorials
--------------------

.. code-block:: bash

    htesp-tutorials --list                      # the catalogue, every step
    htesp-tutorials --code QE                   # one example tree
    htesp-tutorials --only QE/9,QE/12           # just these
    htesp-tutorials --skip QE/4,QE/5            # all but these
    htesp-tutorials --from band-scf             # start each one at this step
    htesp-tutorials --skip-stubs                # leave out what cannot run as shipped

Other options: ``--workdir DIR`` (root for work directories, logs, checkpoint
and report), ``--workers N`` (passed through to ``mainprogram --workers``),
``--poll-interval S`` (seconds between ``squeue`` polls, default 60),
``--job-timeout S`` (how long to wait for a step's jobs, default 24 h),
``--step-timeout S`` (how long one ``mainprogram`` call may take, default 6 h),
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
* under ``--no-dft``, everything that *would* have been submitted;
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

``preflight()`` reports every problem before the first job is submitted:

* ``python -m htesp`` must start;
* ``examples/<code>/batch.header`` must exist, and QE needs
  ``examples/QE/pp/*.upf``;
* **``MP_API_KEY`` must be exported** for the database tutorials (QE/2, QE/3,
  QE/6, QE/16 and the VASP equivalents).  The API key was removed from every
  shipped ``config.json``, so those tutorials cannot run without it.  Under
  ``--dry-run`` this is a warning and the steps are skipped; in the real mode it
  is an error;
* in the real mode: ``sbatch`` and ``squeue`` on ``$PATH``, plus ``pw.x`` (QE) or
  ``vasp_std`` (VASP).

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
  stub in the catalogue and cannot run as shipped; ``--skip-stubs`` leaves it
  out.
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
Running a real campaign in waves
--------------------------------

In ``--dry-run`` and ``--no-dft`` nothing is ever queued, so the whole
selection runs in one pass and there is nothing to wait for.  A real run is
different: wave 0 ends with the relaxations sitting in the scheduler's queue,
and the tutorials that read the relaxed structure cannot honestly start until
those have finished.  Running everything in one pass either blocks for days
inside ``squeue`` polling or proceeds on a structure that is not relaxed yet --
and the second is silently wrong, which is worse than failing.

``--wave`` runs one dependency wave and stops:

.. code-block:: bash

    htesp-tutorials --list-waves            # the plan, and what is done so far
    htesp-tutorials --wave next             # run the first unfinished wave
    # ... wait for its jobs to drain ...
    htesp-tutorials --resume --wave next    # the next one
    htesp-tutorials --wave 1                # or name one exactly

The waves are derived from ``depends_on`` alone, so they stay correct as the
catalogue changes.  Wave 0 is every tutorial that depends on nothing: input
generation, the database front ends and the relaxation hubs.  Wave 1 is
everything that reads what wave 0 produced.

Each wave's status is written to ``state.json`` under ``waves``, so days later
the runner can still say which wave finished and which is next.  A wave counts
as done only when every tutorial in it did; an interrupted or failed wave is
offered again by ``--wave next`` rather than stepped over, because the wave
after it would run on structures that were never produced.  ``--restart --only
QE/9`` un-finishes any wave containing ``QE/9`` for the same reason.

A selection that leaves a dependency out makes the dependent a root -- wave 0
-- because :func:`~tutorials.catalog.topological_order` drops dependencies
outside the selection.  ``--list-waves`` names those cases rather than letting
the wave numbers imply an ordering that will not happen.
What "the step finished" means in a real run
--------------------------------------------

Two things that look like success are not, and both used to be recorded DONE.

**A job leaving the queue is not a job succeeding.**  ``COMPLETED``,
``FAILED``, ``TIMEOUT``, ``CANCELLED`` and ``OUT_OF_MEMORY`` all make a job id
disappear from ``squeue`` alike.  After the wait, the runner asks ``sacct -X -o
JobID,State`` how each one ended and fails the step unless every job is
``COMPLETED``, naming the state in words (``7 TIMEOUT (hit the wall time)``).
Where accounting is not enabled ``sacct`` answers nothing, and that is treated
as "cannot tell", not as failure.

**A relaxation finishing is not a relaxation converging.**  Quantum ESPRESSO
stops at ``nstep`` and VASP at ``NSW``, both exiting cleanly with a complete set
of output files -- so the artefact globs match and the job is ``COMPLETED``.
Steps marked ``check_converged`` therefore read the output for the same markers
``htesp/workflow.py`` uses:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - file
     - what must be in it
   * - ``OUTCAR``
     - ``reached required accuracy - stopping structural energy minimisation``
   * - ``relax.out``
     - ``End of BFGS Geometry Optimization``
   * - ``scf.out``
     - ``JOB DONE.``

``The maximum number of steps has been reached`` and ``convergence NOT
achieved`` are proof of the opposite.  Four steps carry the flag:
``relax-submit`` and ``resubmit`` in the relaxation hub, ``relax-deformed`` in
the elastic tutorial and ``relax-volumes`` in the equation-of-state series.

Both checks run only in real mode; ``--dry-run`` and ``--no-dft`` never produce
an ``OUTCAR``, so applying them there would fail every relaxation for the wrong
reason.  Together with the wave gate above, a relaxation that did not converge
now fails its step, which fails its wave, which leaves every dependent tutorial
BLOCKED -- instead of them quietly starting from an unrelaxed cell.

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
``mainprogram jobscript`` submitted nothing in real mode and still exited 0.

Every tutorial with a submitting step therefore gets a ``jobscript`` step
prepended (:func:`tutorials.catalog._with_job_scripts`), and a submitting step
that records no job id is now a failure rather than an "unverifiable": whatever
the cause, nothing was submitted.

On a machine with no ``sinfo`` nothing is generated -- the probes would have
nothing to say, and a laptop ``--dry-run`` submits nothing anyway.
