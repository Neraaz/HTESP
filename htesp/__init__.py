"""HTESP -- High Throughput Electronic Structure Package.

High-throughput driver for Quantum ESPRESSO and VASP: database search and
download, relaxation, electron-phonon coupling and superconducting Tc, band
structures, densities of states, phonons, elastic constants, EPW/Wannier90 and
WannierTools.

Layout
------
``htesp.cli``        the ``mainprogram`` command line
``htesp.workflow``   every former ``src/bash`` scan script, as a method
``htesp.config``     configuration loading and validation
``htesp.inputin``    the ``input.in`` control file
``htesp.help_text``  the long help blocks (``docs/command.rst`` is generated from these)
everything else      the science modules, unchanged in purpose
"""
from htesp.banner import __version__

__all__ = ["__version__"]
