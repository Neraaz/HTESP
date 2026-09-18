.. _workflow-label:

====================================
Parallelism and the workflow layer
====================================

Every stage of an HTESP campaign is a loop over materials: read the tracking
file, take entries ``start`` through ``end`` (exclusive), and do the same thing
to each.  That loop lives in ``htesp/workflow.py``, in the class
``HTESPWorkflow``, and it runs in parallel.

--------------------
HTESPWorkflow
--------------------

``HTESPWorkflow`` is one method per stage.  Each method keeps the name and the
argument contract of the scan script it replaces -- ``start``, ``end``
(exclusive), the tracking file, and an optional extra argument -- so nothing
about how you drive a campaign has changed:

.. code-block:: python

    from htesp.workflow import HTESPWorkflow

    wf = HTESPWorkflow()                       # reads ./input.in
    wf.relax_scan(1, 5, "mpid.in")
    wf.create_inputs(1, 5, "mpid.in", nkpt=50)

or from the command line:

.. code-block:: bash

    python -m htesp.workflow relax-scan 1 5 mpid.in
    python -m htesp.workflow create-inputs 1 5 mpid.in 50 --workers 8

You do not normally call it directly.  ``mainprogram`` dispatches to it, and the
shims in ``bin/`` do too.

----------------------------
bin/ and legacy/bash/
----------------------------

``bin/`` holds 52 one-line shims, one per former bash scan script, with the same
names and the same argument contract.  Putting ``bin/`` on ``$PATH``:

.. code-block:: bash

    export PATH="path_to_HTESP/bin:$PATH"

reproduces exactly what ``export PATH=$PATH:path_to_HTESP/src/bash`` used to do.
``relax-scan 1 5 mpid.in`` still works and still means the same thing;
``bin/jobscript.sh`` is still sourceable, and now records the job ids it submits.

The original bash scripts are kept unmodified under ``legacy/bash/``.  They are
there to be read -- to check what the old behaviour was -- not to be run, and
they are not installed and not on ``$PATH``.

--------------------
--workers
--------------------

The per-material body runs in a process pool:

.. code-block:: bash

    mainprogram 1 --workers 16
    export HTESP_WORKERS=4      # the default for every later command

The default pool size is ``min(cpu_count, 8)``.  ``--workers 1`` runs the loop
serially, which is what you want when you are debugging a single material or
reading interleaved output.  A pool is not created at all when there is only one
material to process.

This is worth having on the login node too, where the work per material is
building input files rather than running DFT: a 300-compound ``mainprogram 4``
is bounded by ``pymatgen`` and file writing, and a pool cuts it by close to the
number of workers.

-------------------------------
Per-material scratch directories
-------------------------------

The old bash layer wrote fixed-name scratch files into the project root:
``mass.dat``, ``qpoint.dat``, ``kpoint.dat``, ``BZ.pdf``,
``scf_dir/kpathlines.dat``, ``temp*.in``.  With one material at a time that was
merely untidy.  It is why the loop could not simply be backgrounded.

Each material body now runs inside its own private scratch directory, carrying a
local ``scf_dir``.  The helper modules keep their existing "write next to me"
behaviour, and the artefacts they produce are moved to their final destination
when the body finishes.  Two materials processed at the same time cannot see
each other's temporary files.

The working directory is also restored after every iteration, and a failed
``cd`` raises instead of letting the rest of the body -- ``rm -r ...``,
``sbatch ...`` -- run in the project root.

--------------------
Deterministic output
--------------------

Results are collected and re-ordered by material index before anything is
written.  ``result.csv``, ``econv.csv``, ``elastic.csv`` and the ``mpid-*.in``
lists therefore come out in tracking-file order, whatever order the processes
happened to finish in.  Two runs over the same inputs produce byte-identical
files, so they can be diffed and committed.

--------------------
--dry-run
--------------------

.. code-block:: bash

    mainprogram 4 --dry-run
    mainprogram 20 --dry-run

``--dry-run`` builds every input file and does all the file plumbing, but never
calls the scheduler.  Use it to check a campaign before committing cluster time.

It is not a guard against the two destructive processes: ``20``
(``clean-scan``, which deletes wavefunctions and moves the run to
``completed/``) and ``28`` (``pressure-reset``, which removes the ``pressure/``
tree) delete whatever the ``start``/``end`` range in ``input.in`` selects,
``--dry-run`` or not.  Check the range before running either.

--------------------
Job ids
--------------------

Jobs are submitted with ``sbatch --parsable`` and the returned ids are written
to ``.htesp_job.json`` in the stage directory:

.. code-block:: json

    {"ph": [{"job": 4831201, "time": 1758067200}]}

``mainprogram checkph`` reads those ids and asks ``squeue -h -o %i -j <ids>``
whether they are still queued.  The old check grepped the whole queue for the
compound name, so a job named ``B`` matched every job in the queue containing a
``B``, and elemental compounds were reported as running forever.  A stage whose
job ids were not recorded, or a machine with no ``squeue``, is now reported as
*unverifiable* rather than silently passed.

The same file is what the :doc:`tutorial runner <tutorial_runner>` polls when it
waits for a step to finish.

--------------------
Failure accounting
--------------------

Each material's result carries its own status.  A stage in which every material
failed exits ``1``; a malformed ``input.in`` or an unknown command exits ``2``;
Ctrl-C exits ``130``.  In HTESP 1.x every stage exited ``0``, because the bash
layer discarded the exit code of each helper it called, so a campaign could fail
completely while every ``&&`` in a job script carried on.
