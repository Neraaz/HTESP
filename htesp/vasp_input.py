#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D.
"""This script will download vasp input files from materials project database.
Script is run within 'download-input' bash script."""
import os
import sys
import shutil
import warnings
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from htesp.cif_to_gsinput import pos_to_kpt, register_mpid
from htesp.write_potcar import poscar2potcar
from htesp.htepc import MpConnect
from htesp.check_json import config
# To ignore all warnings
warnings.filterwarnings("ignore")
def vasp_input(mpid,compound):
    """
    Download VASP input files from the Materials Project.
    Downloads INCAR, POSCAR, POTCAR, and KPOINTS inside Rmpid-compound/relax folder.
    Updates 'mpid.in' with entry number, mpid, and compound name.

    Parameters:
    - mpid (str): Materials ID.
    - compound (str): Compound name.
    """
    input_data = config()
    obj = MpConnect()
    obj.setting(mpid)
    obj.download()
    # FIX(all): makedirs(exist_ok=True) + shutil.move instead of isdir-then-mkdir
    # and os.system("mv ..."), so mpid is never interpolated into a shell string.
    os.makedirs("input_cif", exist_ok=True)
    cif_name = "{}.cif".format(mpid)
    if os.path.isfile(cif_name):
        shutil.move(cif_name, os.path.join("input_cif", cif_name))
    #obtain vasp inputs from MPRelaxSet
    structure = SpacegroupAnalyzer(obj.structure, symprec=0.1).get_primitive_standard_structure()
    relax_set = MPRelaxSet(structure=structure)
    #create a folder with structures
    relax_dir = os.path.join("R{}-{}".format(mpid, compound), "relax")
    os.makedirs(relax_dir, exist_ok=True)
    relax_set.write_input(output_dir=relax_dir)
    relax_set.poscar.write_file(os.path.join(relax_dir, "POSCAR"))
    if os.path.isfile("config.json"):
        shutil.copy("config.json", relax_dir)
    # FIX(all): poscar2potcar() reads POSCAR from the working directory, so the
    # chdir stays -- but it is now restored in a finally block, which the old
    # os.chdir("../../") was not (an exception left the process in relax/).
    pwd = os.getcwd()
    try:
        os.chdir(relax_dir)
        poscar2potcar()
    finally:
        os.chdir(pwd)
    # check_json.config() always returns a populated dict now, so the guard
    # that used to leave ``d`` unbound is gone.
    d = input_data['download']
    evenkpt = d['inp']['evenkpt']
    kptden = input_data['kptden']
    poscar_path = os.path.join(relax_dir, "POSCAR")
    if evenkpt:
        print("Even kpoint mesh is utilized\n")
        pos_to_kpt(poscar_path,kptden,True)
    else:
        pos_to_kpt(poscar_path,kptden)
    # pos_to_kpt() writes KPOINTS into the working directory.
    if os.path.isfile("KPOINTS"):
        shutil.move("KPOINTS", os.path.join(relax_dir, "KPOINTS"))
    # FIX(18): atomic, idempotent, densely renumbered registry update.
    register_mpid(obj.mpid, compound)
def main(mpid=None, compound=None, argv=None):
    """
    main function.

    ``mpid``/``compound`` may be passed directly; they otherwise come from the
    command line, keeping the ``vasp_input.py <mpid> <compound>`` contract used
    by ``htesp.workflow.run_helper``.
    """
    if mpid is None or compound is None:
        argv = list(sys.argv[1:] if argv is None else argv)
        if len(argv) < 2:
            raise SystemExit("usage: vasp_input.py <mpid> <compound>")
        mpid, compound = argv[0], argv[1]
    vasp_input(mpid,compound)
if __name__ == "__main__":
    main()
