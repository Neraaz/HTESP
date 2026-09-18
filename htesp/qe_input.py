#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D.
"""This script will download QE input files from materials project database.
Script is run within 'download-input' bash script."""
import sys
import os
import shutil
from htesp.htepc import MpConnect
from htesp.check_json import config
from htesp.cif_to_gsinput import find_mpid, register_mpid

def qe_input(mpid):
    """
    Extract QE input files for ground state calculations.

    This function retrieves QE input files for ground state calculations using the provided
    Materials ID (mpid). It saves the structures in CIF format inside the "input_cif"
    folder and the QE scf-mpid.in file inside the "scf_dir" folder. It updates the "mpid.in"
    file with the entry number, mpid, and compound name.

    Parameters:
    - mpid (str): Materials ID

    Returns:
    None
    """
    input_data = config()
    magnetic = input_data['pwscf_in']['magnetic']
    obj = MpConnect()
    obj.setting(mpid)
    obj.maxecut_sssp()
    obj.getkpt()
    if magnetic:
        default_magmoms = input_data['magmom']['magmom']
        obj.structure.add_spin_by_element(default_magmoms)
    # check_json.config() always returns a populated dict now, so the guard
    # that used to leave ``d`` unbound is gone.
    d = input_data['download']
    evenkpt = d['inp']['evenkpt']
    if evenkpt:
        print("Utilizing even kpoint mesh\n")
        obj.getevenkpt()
    # FIX(18): the registry update used to be an open('mpid.in').readlines()
    # followed by an append numbered v<len(lines)+1>, inside this otherwise
    # per-material function -- two concurrent materials produced two "v7"
    # lines and every consumer resolves a material with grep "v$ii ".  It is
    # now a read, an idempotency check, and one atomic, densely renumbered
    # rewrite.  A material already present is not downloaded again, which is
    # what the old ``if not any(mpid in line ...)`` test meant.
    if find_mpid(mpid) is None:
        obj.download()
        if magnetic:
            obj.setting_qeinput(magnetic=True,pseudo_dir='../../pp/')
        else:
            obj.setting_qeinput(pseudo_dir='../../pp/')
        register_mpid(obj.mpid, obj.prefix)
    # FIX(all): os.makedirs(..., exist_ok=True) and shutil.move instead of
    # isdir-then-mkdir and os.system("mv ...").
    os.makedirs("input_cif", exist_ok=True)
    os.makedirs("scf_dir", exist_ok=True)
    scf_name = "scf-{}.in".format(mpid)
    if os.path.isfile(scf_name):
        shutil.move(scf_name, os.path.join("scf_dir", scf_name))
    cif_name = "{}.cif".format(mpid)
    if os.path.isfile(cif_name):
        shutil.move(cif_name, os.path.join("input_cif", cif_name))
def main(mpid=None, argv=None):
    """
    main function.

    ``mpid`` may be passed directly; it otherwise comes from the command line,
    keeping the ``qe_input.py <mpid>`` contract.
    """
    if mpid is None:
        argv = list(sys.argv[1:] if argv is None else argv)
        if not argv:
            raise SystemExit("usage: qe_input.py <mpid>")
        mpid = argv[0]
    qe_input(mpid)
if __name__ == "__main__":
    main()
