#!/usr/bin/env python
# coding: utf-8
# Written by Niraj K. Nepal, Ph.D.
"""
Creates different magnetic ordering according to MagneticStructureEnumerator functions.
"""
import logging
import os
import shutil
import subprocess
import sys
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core import structure
from pymatgen.io.pwscf import PWInput
from pymatgen.analysis.magnetism import MagneticStructureEnumerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from htesp.cif_to_gsinput import pos_to_kpt
from htesp.write_potcar import poscar2potcar
from htesp.htepc import MpConnect
from htesp.check_json import config

LOG = logging.getLogger("htesp")

# FIX(16): the k-point density was hard-coded to 0.025 here while every other
# module reads config['kptden'], so a campaign that changed kptden silently got
# a different mesh for its magnetic orderings.
DEFAULT_KPTDEN = 0.025


def run_vasp_process(what, cwd=None):
    """Run ``vasp_process.py <what>`` without a shell, optionally in ``cwd``."""
    # FIX(36): ``vasp_process.py`` was an executable only while 1.x put src/ on
    # $PATH.  In 2.0 it is a module, so the bare name is a FileNotFoundError --
    # "[Errno 2] No such file or directory: 'vasp_process.py'" -- and every
    # command routed through here died.  workflow.py and oqmd_extract.py were
    # converted; these four call sites were missed.
    subprocess.run([sys.executable, "-m", "htesp.vasp_process", what],
                   check=True, cwd=cwd)


def rewrite_incar(path, drop=(), add=()):
    """Drop ``drop`` keys from an INCAR and append the ``(key, value)`` ``add``.

    The key is matched at the start of the line, so dropping ``NSW`` does not
    also drop ``NSWRITE`` the way ``sed -i '/NSW/d'`` did.
    """
    drop = {key.upper() for key in drop}
    try:
        with open(path, "r") as handle:
            lines = handle.read().splitlines()
    except OSError:
        lines = []
    lines = [ln for ln in lines if ln.split("=")[0].strip().upper() not in drop]
    for key, value in add:
        lines.append("{} = {}".format(key, value))
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")


def magnetic_structure(obj,mpid,compound,magconfig,dft):
    """
    Functions to generate possible magnetic structures.

    Parameters:
    obj (MpConnect): An object representing an instance of the MpConnect class.
    mpid (str): Materials ID.
    compound (str): Compound name.
    magconfig (list): List of possible collinear magnetic configurations, e.g., ['ferromagnetic', 'antiferromagnetic'].
    dft (str): Type of DFT calculation, e.g., 'VASP' or 'QE'.

    Look for pymatgen.analysis.magnetism.MagneticStructureEnumerator class for more details.
    Install Enumlib library: https://github.com/msg-byu/enumlib to run this module.

    Example:
    >>> obj = MpConnect()
    >>> magnetic_structure(obj, 'mp-123', 'MnO', ['ferromagnetic', 'antiferromagnetic'], 'VASP')
    """
    input_data = config()
    # load the structure
    try:
        strucinit = structure.Structure.from_file("R{}-{}/relax/POSCAR".format(mpid,compound))
    except FileNotFoundError:
        strucinit = PWInput.from_file("R{}-{}/relax/scf.in".format(mpid,compound)).structure
    # Obtain magnetic structures based on given magnetic configurations
    if input_data.get('magmom', {}).get('magmom'):
        default_magmoms = input_data['magmom']['magmom']
        order = input_data['magmom']['order']
        newstructure = MagneticStructureEnumerator(strucinit,default_magmoms=default_magmoms,strategies=order,truncate_by_symmetry=True).ordered_structures
    else:
        newstructure = MagneticStructureEnumerator(strucinit,strategies=magconfig,truncate_by_symmetry=True).ordered_structures
    print(len(newstructure))
    # FIX(16): honour config['kptden'], falling back to the historical 0.025
    kptden = input_data.get('kptden') or DEFAULT_KPTDEN
    # Check the entry number
    if not os.path.isfile('mpid-magnetic.in'):
        entry = 0
    else:
        with open('mpid-magnetic.in', 'r') as mpid_read:
            lines = mpid_read.readlines()
        entry = len(lines)
    # Process each generated structure
    if len(newstructure) > 0:
        for i,struc in enumerate(newstructure):
            # Refine the structure
            struc = SpacegroupAnalyzer(struc,symprec=0.1).get_refined_structure()
            if dft in ('vasp', 'VASP'):
                obj.prefix = struc.composition.alphabetical_formula.replace(" ","")
                poscar = Poscar(structure=struc,comment=obj.prefix)
                relax_dir = "R{}-{}-{}/relax".format(mpid,i+1,obj.prefix)
                os.makedirs(relax_dir, exist_ok=True)
                orig_prefix = compound
                shutil.copy("R{}-{}/relax/INCAR".format(mpid,orig_prefix), relax_dir)
                poscar.write_file(filename="{}/POSCAR".format(relax_dir))
                evenkpt = input_data['download']['inp']['evenkpt']
                if evenkpt:
                    print("Even kpoint mesh is utilized\n")
                    pos_to_kpt("{}/POSCAR".format(relax_dir),kptden,True)
                else:
                    pos_to_kpt("{}/POSCAR".format(relax_dir),kptden)
                shutil.move("KPOINTS", "{}/KPOINTS".format(relax_dir))
                if os.path.isfile("vasp.in"):
                    shutil.copy("vasp.in", relax_dir)
                if os.path.isfile("config.json"):
                    shutil.copy("config.json", relax_dir)
                # poscar2potcar() and vasp_process work on the current
                # directory; the chdir is wrapped so that an exception cannot
                # leave the process inside the ordering directory.
                pwd = os.getcwd()
                try:
                    os.chdir(relax_dir)
                    # Process POTCAR and INCAR
                    poscar2potcar()
                    run_vasp_process("POSCAR")
                    # Modify INCAR for magnetic calculations.  The old code
                    # deleted MAGMOM by line number while also deleting the
                    # ISPIN lines, which shifts the numbering.
                    maglist = ""
                    for _,specie in enumerate(struc.species):
                        if specie.spin is not None:
                            maglist += str(specie.spin) + " "
                        else:
                            maglist += "0 "
                    rewrite_incar("INCAR", drop=("MAGMOM", "ISPIN"),
                                  add=(("ISPIN", 2), ("MAGMOM", maglist.strip())))
                finally:
                    os.chdir(pwd)
            else:
                # Convert IStructure to Structure
                struc = structure.Structure(struc.lattice, struc.species, struc.cart_coords, coords_are_cartesian=True)
                obj.structure = struc
                comp_list = []
                for composition in struc.composition.elements:
                    comp_list.append(composition.name)
                obj.comp_list = comp_list
                obj.getkpt(primitive=False)
                evenkpt = input_data['download']['inp']['evenkpt']
                if evenkpt:
                    print("Utilizing even kpoint mesh\n")
                    obj.getevenkpt()
                obj.maxecut_sssp_for_subs()
                obj.prefix = struc.composition.alphabetical_formula.replace(" ","")
                # FIX(15): setting_qeinput has no 'monoclinic' parameter -- the
                # call raised TypeError before it wrote anything.  The intent
                # (matching obj.getkpt(primitive=False) above) is to keep the
                # ordering exactly as MagneticStructureEnumerator produced it,
                # which is primitive=False; the monoclinic setting is now the
                # module constant htesp.htepc.INTERNATIONAL_MONOCLINIC.
                obj.setting_qeinput(magnetic=True,primitive=False,pseudo_dir='../../pp')
                os.makedirs("scf_dir", exist_ok=True)
                shutil.move("scf-{}.in".format(obj.mpid),
                            "scf_dir/scf-{}-{}.in".format(mpid,i+1))
            with open("mpid-magnetic.in", "a") as mpid_append:
                mpid_append.write("v{}".format(entry+1+i) + " " + mpid + "-{}".format(i+1) + " " + obj.prefix + "\n")
    else:
        print("Structures not created\n")


def changeaxis(mpid,comp,dft):
    """
    Function to change magnetic axis to compute magnetic anisotropy.

    Parameters:
    mpid (str): Material ID.
    comp (str): Compound name.
    dft  (str): DFT method. Available options: QE/VASP

    The function iterates over magnetic axis configurations specified in the saxis variable,
    updates the necessary input files for VASP calculations, and performs the VASP process.
    """
    input_data = config()
    if not os.path.isfile('mpid-magnetic.in'):
        entry = 0
    else:
        with open('mpid-magnetic.in', 'r') as mpid_read:
            lines = mpid_read.readlines()
        entry = len(lines)
    if dft in ('VASP', 'vasp'):
        # Define magnetic axis configurations
        saxis = input_data['magmom']['saxis']
        # config.json: magmom.force_theorem -- reuse the relaxed CHGCAR with
        # ICHARG = 11 instead of converging each axis from scratch
        force_theorem = bool(input_data['magmom'].get('force_theorem', False))
        # Process each magnetic axis configuration
        for i,axis in enumerate(saxis):
            sx = int(axis[0])
            sy = int(axis[1])
            sz = int(axis[2])
            # Create directory for the current magnetic axis configuration
            axis_dir = "R{}-saxis-{}{}{}-{}/relax".format(mpid,sx,sy,sz,comp)
            os.makedirs(axis_dir, exist_ok=True)
            source = "R{}-{}/relax".format(mpid,comp)
            # Copy necessary input files
            for name in ("INCAR", "POTCAR", "POSCAR", "KPOINTS"):
                shutil.copy(os.path.join(source, name), axis_dir)
            if os.path.isfile(os.path.join(source, "CHGCAR")):
                shutil.copy(os.path.join(source, "CHGCAR"), axis_dir)
            # FIX: this was a commented-out line, so the CHGCAR copied just
            # above was never actually used non-self-consistently and the
            # "force theorem" anisotropy run was in fact a fresh self-consistent
            # one.  It is now driven by config.json: magmom.force_theorem.
            if force_theorem:
                if os.path.isfile(os.path.join(axis_dir, "CHGCAR")):
                    rewrite_incar(os.path.join(axis_dir, "INCAR"),
                                  drop=("ICHARG",), add=(("ICHARG", 11),))
                else:
                    LOG.warning(
                        "magmom.force_theorem is set but %s has no CHGCAR; "
                        "running self-consistently instead", source)
            # Update INCAR with the new magnetic axis
            rewrite_incar(os.path.join(axis_dir, "INCAR"), drop=("SAXIS",),
                          add=(("SAXIS", "{} {} {}".format(sx,sy,sz)),))
            # Append to the mpid-magnetic.in file
            with open("mpid-magnetic.in", "a") as mpid_append:
                mpid_append.write("v{}".format(entry+1+i) + " " + mpid + "-saxis-{}{}{}".format(sx,sy,sz) + " " + comp + "\n")
            if os.path.isfile("config.json"):
                shutil.copy("config.json", axis_dir)
            if os.path.isfile("vasp.in"):
                shutil.copy("vasp.in", axis_dir)
            run_vasp_process("POSCAR", cwd=axis_dir)
    elif dft in ("QE", "qe"):
        # FIX: this branch printed "To be implemented" and returned, so
        # `mainprogram magenum` on a QE campaign exited 0 having written
        # nothing.  It now says so through the failure path.
        raise NotImplementedError(
            "magnetic-anisotropy enumeration (mainprogram magenum) is "
            "implemented for VASP only.  For Quantum ESPRESSO, set "
            "pwscf_in.magnetic = true in config.json and use the ordinary "
            "download/relax pipeline, which writes spin-polarised inputs.\n"
            "The outline of the missing QE implementation is kept in the "
            "comments below this message in htesp/magnetic.py.")
        # Perform scf calculation with magnetic configuration
        # extract (1)(2)... from starting_magnetization(1)(2)...
        # One idea is to grep "starting_magnetization" and
        # record its occurence. Now loop over and write following
        # in two files.
        # Define angle1(1)=0,angle2(1)=0 on one nscf.in
        # Define angle1(1)=90,angle2(1)=0 on second file.
        # Prepare 'nscf' calculation with lforcet = .true., nosym = .true., and startingpot='file'
        # Also need lspinorb = .true. and noncolin = .true. (can be included with obj.setting_qeinput())
        # Name this file nscf-{mpid}-{saxis_index}.in
        # Create similar folder as of vasp
        # copy scf.in and nscf.in files inside it
        # append other commands to run nscf.in file to run-scf.sh
        # Finally add another command to apply force theorem with projwfc.x
    else:
        raise ValueError(
            f"download.inp.calc is {dft!r}; expected 'QE' or 'VASP'")


def main():
    """
    Main function to control the workflow.

    Reads input data, determines the type of calculation and magnetic configuration,
    and performs the corresponding operations.

    Parameters:
    None

    Returns:
    None
    """
    input_data = config()
    # DFT calculation type
    dft = input_data['download']['inp']['calc']
    # Ordering or anisotropy ?
    mag_type = input_data['magmom']['type']
    # Read input.in to obtain list of materials
    with open("input.in","r") as read_in:
        lines = read_in.readlines()
    start = int(lines[0].split("\n")[0])
    end = int(lines[1].split("\n")[0])
    filename = lines[3].split("\n")[0]
    with open(filename,'r') as read_mpid:
        lines = read_mpid.readlines()
    lines = lines[start-1:end-1]
    # Initiate MpConnect object
    obj = MpConnect()
    config_mag = ['ferromagnetic','antiferromagnetic']
    # Loop over materials
    for line in lines:
        mpid = line.split(" ")[1]
        comp = line.split(" ")[2].split("\n")[0]
        # Creating different magnetic ordering
        if mag_type == "ordering":
            magnetic_structure(obj,mpid,comp,config_mag,dft)
        # Creating input files with different SAXIS
        elif mag_type == "anisotropy":
            changeaxis(mpid,comp,dft)
        else:
            print("Only ordering and anisotropy allowed\n")


if __name__ == "__main__":
    main()
