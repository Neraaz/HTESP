#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D.
"""This script will download QE input files from materials project database.
Script is run within 'download-input' bash script."""
import sys
import os
from htepc import MpConnect
from check_json import config

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
    if os.path.isfile("config.json") or os.path.isfile("../../config.json"):
        d = input_data['download']
    evenkpt = d['inp']['evenkpt']
    if evenkpt:
        print("Utilizing even kpoint mesh\n")
        obj.getevenkpt()
    if not os.path.isfile('mpid.in'):
        kind = 0
        obj.download()
        if magnetic:
            obj.setting_qeinput(magnetic=True,pseudo_dir='../../pp/')
        else:
            obj.setting_qeinput(pseudo_dir='../../pp/')
        with open("mpid.in", "w") as write_mpid:
            write_mpid.write("v{}".format(kind+1) + " " + obj.mpid + " " + obj.prefix + "\n")
    else:
        with open('mpid.in', 'r') as read_mpid:
            lines = read_mpid.readlines()
        kind = len(lines)
        if not any(mpid in line for line in lines):
            obj.download()
            if magnetic:
                obj.setting_qeinput(magnetic=True,pseudo_dir='../../pp/')
            else:
                obj.setting_qeinput(pseudo_dir='../../pp/')
            with open("mpid.in", "a") as write_mpid:
                write_mpid.write("v{}".format(kind+1) + " " + obj.mpid + " " + obj.prefix + "\n")
    if not os.path.isdir("input_cif"):
        os.mkdir("input_cif")
    if not os.path.isdir('scf_dir'):
        os.mkdir("scf_dir")
    if os.path.isfile("scf-{}.in".format(mpid)):
        os.system("mv scf-{}.in scf_dir/".format(mpid))
    if os.path.isfile("{}.cif".format(mpid)):
        os.system("mv {}.cif input_cif/".format(mpid))
def main():
    """
    main function
    """
    mpid=sys.argv[1]
    qe_input(mpid)
if __name__ == "__main__":
    main()
