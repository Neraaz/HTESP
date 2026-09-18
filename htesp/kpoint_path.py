#!/usr/bin/env python
#"""Written by Niraj K. Nepal, Ph.D"""
"""Module to write the ``Kpoint_Path`` card used by wannier90 / WannierTools."""
import os
import re
from ase.io import espresso, vasp
from ase.cell import Cell
from htesp.kpath import kpath

# FIX(10): ASE reports a discontinuity in the band path as a single combined
# label ("K|U", sometimes "K,U").  The old code appended that combined string
# verbatim and then paired *every* consecutive entry, which both raised a
# KeyError when looking the combined label up in ``special_points`` and emitted
# a bogus "K -> U" segment straight across the break in the path.
_BREAK = re.compile(r"[|,]")


def split_path_segments(labels):
    """
    Split a list of high-symmetry labels into continuous segments.

    A label containing ``|`` or ``,`` (e.g. ``"K|U"``) marks a *break* in the
    band path: the part before the separator closes the running segment and the
    part after it opens a new one.  Segments shorter than two labels carry no
    line and are dropped.

    Parameters
    ----------
    labels : list of str
        High-symmetry point names as returned by ``kpath``.

    Returns
    -------
    list of list of str
        One list of labels per continuous stretch of the path.

    Example
    -------
    >>> split_path_segments(['G', 'X', 'W', 'K|U', 'X'])
    [['G', 'X', 'W', 'K'], ['U', 'X']]
    """
    segments = [[]]
    for raw in labels:
        parts = [part for part in _BREAK.split(str(raw)) if part]
        if not parts:
            continue
        segments[-1].append(parts[0])
        for part in parts[1:]:
            segments.append([part])
    return [seg for seg in segments if len(seg) > 1]


def _write_segments(handle, segments, coords, ndim=3):
    """Write ``label x y z <tab> label x y z`` pairs for each continuous segment."""
    for segment in segments:
        for i in range(len(segment) - 1):
            for which, label in enumerate((segment[i], segment[i + 1])):
                point = coords[label]
                handle.write(label + " ")
                handle.write(" ".join(str(round(point[j], 5)) for j in range(ndim)))
                handle.write("\t" if which == 0 else "\n")


def _missing_labels(segments, coords):
    """Return the labels used by ``segments`` that ``coords`` does not define."""
    return sorted({lbl for seg in segments for lbl in seg if lbl not in coords})


def kpoint_path(file_name, out="wannier_kpath.in", slab=False,
                slab_file="POSCAR-slab", slab_out="wannier_kpath_slab.in"):
    """
    Write the ``Kpoint_Path`` card for wannier90 / WannierTools.

    Parameters
    ----------
    file_name : str
        QE ``scf.in`` file (a VASP ``POSCAR`` is accepted as a fallback).
    out : str
        Output file. Default ``'wannier_kpath.in'``.
    slab : bool
        Also write the 2D path for ``slab_file`` when that file exists.
        Used by :mod:`htesp.create_wt_inputs`. Default ``False``.
    slab_file : str
        Slab structure to read for the 2D path. Default ``'POSCAR-slab'``.
    slab_out : str
        Output file for the 2D path. Default ``'wannier_kpath_slab.in'``.

    Returns
    -------
    None
    """
    _, _, _, _, sym, _ = kpath(file_name, 1, 0)
    try:
        data = espresso.read_espresso_in(file_name)
    # FIX: a bare ``except`` swallowed everything, KeyboardInterrupt included.
    except (OSError, ValueError, IndexError, KeyError):
        data = vasp.read_vasp(file_name)
    band = Cell.bandpath(data.cell)
    band_dict = band.special_points
    segments = split_path_segments(sym)
    missing = _missing_labels(segments, band_dict)
    if missing:
        raise KeyError(
            "high-symmetry labels {} are not defined in the band path of {}".format(
                missing, file_name))
    with open(out, 'w') as wan_kpath_write:
        _write_segments(wan_kpath_write, segments, band_dict, ndim=3)

    # FIX(7): the 2D slab path used to live in a duplicate copy of this
    # function inside create_wt_inputs.py; it is folded in here behind a flag.
    if slab and os.path.isfile(slab_file):
        slab_data = vasp.read_vasp(slab_file)
        band_s = Cell.bandpath(slab_data.cell, pbc=[1, 1, 0])
        slab_dict = band_s.special_points
        sym_s = band_s.get_linear_kpoint_axis()[2]
        segments_s = split_path_segments(sym_s)
        missing_s = _missing_labels(segments_s, slab_dict)
        if missing_s:
            raise KeyError(
                "high-symmetry labels {} are not defined in the slab band path "
                "of {}".format(missing_s, slab_file))
        with open(slab_out, 'w') as wan_kpath_write:
            _write_segments(wan_kpath_write, segments_s, slab_dict, ndim=2)
