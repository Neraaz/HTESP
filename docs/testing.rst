.. _testing-label:

==========================
Testing and contributing
==========================

--------------------
Running the tests
--------------------

.. code-block:: bash

    pip install -e ".[test]"
    pytest tests/ -q

or, with nothing but the standard library:

.. code-block:: bash

    python -m unittest discover -s tests -t . -v

The tests are ``unittest`` classes so that both runners work.  That matters
because most of them need nothing installed: ``mainprogram --help``, the
configuration loader, ``input.in``, the tracking-file reader, the Quantum
ESPRESSO text handling and the whole command dispatcher are checked on a machine
with no ``pymatgen``, no ``ase`` and no cluster.  Anything that genuinely needs a
scientific package is skipped with a message naming it, rather than erroring.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - File
     - What it pins
   * - ``tests/test_config.py``
     - the search order, the deep merge over the packaged default,
       ``MP_API_KEY`` precedence, ``config-validate``
   * - ``tests/test_inputin.py``
     - ``input.in`` parsing, including short and ``DFT``-less files
   * - ``tests/test_qetext.py``
     - the QE reader: ``vc-relax``, ``relax``, timed-out outputs, cards,
       namelist edits
   * - ``tests/test_workflow_core.py``
     - tracking files (``end`` is exclusive), material paths, mesh files,
       ``pushd``, job-id capture, failure accounting
   * - ``tests/test_cli.py``
     - every command's dispatch, the plot-type loop, exit codes, the help blocks
   * - ``tests/test_imports.py``
     - every module imports, and nothing reads a file at import time
   * - ``tests/test_packaging.py``
     - the shims, the entry points, no committed API key, no ``os.system``, no
       bare ``except:``, no ``eval``
   * - ``tests/test_regressions.py``
     - one test per defect from the review, each naming the original symptom
   * - ``tests/test_docs.py``
     - the documentation cannot name a command that does not exist, cannot tell
       the reader to run the destructive process ``20`` for the partial DOS, its
       JSON blocks must parse, its labels must be unique and its references must
       resolve
   * - ``tests/test_tutorials.py``
     - the tutorial catalogue, tied to the dispatcher

--------------------
The other checks
--------------------

.. code-block:: bash

    python tools/check_names.py htesp tutorials tests   # undefined-name scan
    python tools/gen_command_rst.py --check             # docs/command.rst is current
    python docs/gen_param_block.py --check              # the config block in param.rst
    python docs/check_docs.py                           # the documentation checks
    python -m unittest tutorials.selftest               # the tutorial runner alone

``tools/check_names.py`` is a small standard-library substitute for ``pyflakes``
for machines where ``ruff`` cannot be installed.  It catches the one class of
defect that bit this package repeatedly: a name used at module or function scope
that is never bound anywhere.

--------------------
Generated pages
--------------------

Two pages of this documentation are generated and must not be hand-edited.

``docs/command.rst`` comes from ``htesp/help_text.py`` -- the same constants
``mainprogram basicinfo``, ``process-info``, ``epw-info`` and ``wt-info`` print.
Change the help text and regenerate:

.. code-block:: bash

    python tools/gen_command_rst.py

The full ``config.json`` listing in ``docs/param.rst``, between the
``.. config-json-start`` and ``.. config-json-end`` comments, comes from
``htesp/data/config.json``:

.. code-block:: bash

    python docs/gen_param_block.py

Both have a ``--check`` mode that exits non-zero when the file on disk is stale,
and ``tests/test_docs.py`` runs the first of them.

--------------------
Checking the docs
--------------------

``docs/check_docs.py`` is the documentation lint that does not need Sphinx
installed.  It checks that

* every ``.. code-block:: json`` body parses as JSON and every
  ``.. code-block:: python`` body compiles;
* every ``:ref:`` and ``:doc:`` target exists;
* no label is defined twice;
* every page named in a ``toctree`` exists, and every page is reachable from
  ``index.rst``;
* every ``mainprogram <name>`` in the prose is a real command or process, as
  ``htesp/help_text.py`` lists them;
* every repository path mentioned in a literal exists in the tree;
* no real-looking Materials Project API key has crept back in.

.. code-block:: bash

    python docs/check_docs.py          # report and exit non-zero on problems
    python docs/check_docs.py -v       # also list what it checked

Build the HTML with Sphinx when it is available:

.. code-block:: bash

    pip install -e ".[docs]"
    cd docs && make html

--------------------
Conventions
--------------------

* **A test that pins a fix names the original symptom**, in its docstring or in
  its name.  A test that only asserts current behaviour tells a later reader
  nothing about why it exists; a test that says "process 20 is ``clean-scan``,
  and the tutorial used to tell the reader to run it for the partial DOS" tells
  them everything.
* ``tests/helpers.py`` has ``TempProject``, which gives each test a throw-away
  project directory, makes it current, and clears ``$HTESP_CONFIG``,
  ``$MP_API_KEY`` and the configuration cache.  A stale value in a developer's
  shell would otherwise defeat the configuration tests.
* Negative source assertions go through ``code_only()``, which strips comments
  and strings.  Every fix carries a ``# FIX(n):`` comment naming the old
  pattern, so a naive "the old pattern is gone" check would match the comment
  itself.
* Documentation follows the code, not the other way round.  If a page is wrong
  because a command was renamed, fix the help text and regenerate; if it is
  wrong because the prose is wrong, fix the prose and add the wrong name to
  ``tests/test_docs.py`` so it cannot come back.
* Never commit an API key.  ``tests/test_packaging.py`` checks for one.

--------------------
Who to contact
--------------------

See :doc:`contrib`.
