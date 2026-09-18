"""Test suite for HTESP.

The tests are written as ``unittest`` classes on purpose: they run under
``pytest tests/`` and under ``python -m unittest discover -s tests`` alike, so
the core of the package can be checked on a machine that has neither pytest nor
the heavy scientific stack installed.

Anything that needs ``pymatgen``/``ase``/``spglib`` is guarded with
``skip_without`` from :mod:`tests.helpers`, so the suite is green wherever it
runs and simply reports what it could not exercise.
"""
