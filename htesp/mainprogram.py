#!/usr/bin/env python
"""Backwards-compatible entry point.

The dispatcher itself now lives in :mod:`htesp.cli`; this module is kept so
``python -m htesp.mainprogram`` and ``from htesp.mainprogram import main``
keep working.
"""
from htesp.cli import main

__all__ = ["main"]

if __name__ == "__main__":
    raise SystemExit(main())
