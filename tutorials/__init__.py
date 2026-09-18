"""End-to-end tutorial runner for HTESP.

The modules here drive every worked example under ``examples/`` from one
command, checkpoint their progress and -- the point of the exercise -- say
precisely *where* a campaign stopped when it stops.

Nothing in this package imports :mod:`htesp` at module scope: the catalogue,
the checkpoint file, the report writer and the command line are usable (and
testable) on a laptop that has neither ``pymatgen`` nor Quantum ESPRESSO
installed.  ``mainprogram`` is invoked as a *subprocess*, so a tutorial that
dies cannot take the driver down with it.
"""
from __future__ import annotations

__all__ = ["catalog", "state", "runner", "report", "run_tutorials"]
