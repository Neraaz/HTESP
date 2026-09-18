#!/usr/bin/env python
"""Does the symmetry tolerance change the k-mesh that `getkpt` produces?

1.x symmetrised in `MpConnect.getkpt` with symprec=0.1 and pymatgen's default
`international_monoclinic=True`, while `setting_qeinput` used symprec=0.01 with
`international_monoclinic=False` -- so the mesh was computed for a different
cell from the one written to scf-<mpid>.in.  2.0 unifies both on
(0.01, False) via htepc.SYMPREC / htepc.INTERNATIONAL_MONOCLINIC.

This prints the cell and the resulting mesh under both settings for a few
structures chosen to stress the difference.  Finding: the *cell* can differ
visibly (space group, monoclinic setting) while the *mesh* stays the same,
because pos_to_kpt depends only on reciprocal lattice vector lengths.  Run it
before claiming a k-mesh changed.

    python <repo>/tools/diagnostics/kmesh_symprec.py
"""
import os
import shutil
import tempfile
import warnings

warnings.filterwarnings("ignore")
import _path  # noqa: F401,E402

from pymatgen.core import Structure, Lattice  # noqa: E402
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer  # noqa: E402
from pymatgen.io.vasp.sets import MPRelaxSet  # noqa: E402
from htesp.htepc import pos_to_kpt  # noqa: E402

KPTDEN = 0.025          # the packaged default

SETTINGS = (("1.x getkpt (symprec=0.1 , international=True )", 0.1, True),
            ("2.0 getkpt (symprec=0.01, international=False)", 0.01, False))


def cases():
    a = 4.00
    yield ("near-cubic, c/a = 1.004",
           Structure(Lattice.from_parameters(a, a, a * 1.004, 90, 90, 90),
                     ["Sr", "Ti", "O", "O", "O"],
                     [[0, 0, 0], [.5, .5, .5], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]]))
    yield ("C-centred monoclinic, beta = 108.5",
           Structure(Lattice.from_parameters(9.2, 5.1, 6.4, 90, 108.5, 90),
                     ["Al", "Al", "O", "O", "O", "O"],
                     [[0, 0, 0], [.5, .5, 0], [.22, 0, .31], [.78, 0, .69],
                      [.72, .5, .31], [.28, .5, .69]]))


def main():
    previous = os.getcwd()
    workdir = tempfile.mkdtemp(prefix="htesp-kmesh-")
    try:
        os.chdir(workdir)
        for name, struct in cases():
            print("=" * 72)
            print(name)
            for label, symprec, intl in SETTINGS:
                s = SpacegroupAnalyzer(struct, symprec=symprec)\
                    .get_primitive_standard_structure(international_monoclinic=intl)
                MPRelaxSet(structure=s).poscar.write_file("POSCAR")
                mesh = pos_to_kpt("POSCAR", KPTDEN)
                sg = SpacegroupAnalyzer(s, symprec=symprec).get_space_group_symbol()
                abc, ang = s.lattice.abc, s.lattice.angles
                print("  %s" % label)
                print("      %-9s %d atoms  abc=(%.4f, %.4f, %.4f)  "
                      "angles=(%.2f, %.2f, %.2f)"
                      % (sg, len(s), abc[0], abc[1], abc[2], ang[0], ang[1], ang[2]))
                print("      K_POINTS automatic -> %s" % mesh)
    finally:
        os.chdir(previous)
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
