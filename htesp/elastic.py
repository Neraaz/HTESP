#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D. (tug11655@temple.edu)"""
"""Module to perform elastic calculations"""
import os
import shutil
import subprocess
import sys
import warnings
import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io import pwscf
from pymatgen.io.vasp.sets import Vasprun
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import structure
from pymatgen.analysis.elasticity import diff_fit,ElasticTensor,Stress
from pymatgen.analysis.elasticity import DeformedStructureSet,find_eq_stress
from htesp.cif_to_gsinput import pos_to_kpt
from htesp.write_potcar import poscar2potcar
from htesp.htepc import MpConnect
from htesp.check_json import config
from htesp.inputin import InputIn

# FIX(10): VASP reports stress in kBar while pymatgen's ElasticTensor works in
# GPa, so the conversion factor is 0.1, not 1.0; the sign is flipped because
# VASP's stress convention is the opposite of pymatgen's.
KBAR_TO_GPA = 0.1

# FIX(13): the CSV header listed 18 properties while every row skipped
# 'structure' and so had 17 values, shifting every column after
# debye_temperature.  'structure' is not a number and is dropped from both, and
# rows are written in this order rather than in dictionary order, so the header
# and the rows can no longer drift apart.
PROPNAME = ('trans_v', 'long_v', 'snyder_ac', 'snyder_opt', 'snyder_total',
            'clarke_thermalcond', 'cahill_thermalcond', 'debye_temperature',
            'k_voigt', 'k_reuss', 'k_vrh', 'g_voigt', 'g_reuss', 'g_vrh',
            'universal_anisotropy', 'homogeneous_poisson', 'y_mod')


def run_vasp_process(what, cwd=None):
    """Run ``vasp_process.py <what>`` without a shell, optionally in ``cwd``."""
    # FIX(36): ``vasp_process.py`` was an executable only while 1.x put src/ on
    # $PATH.  In 2.0 it is a module, so the bare name is a FileNotFoundError --
    # "[Errno 2] No such file or directory: 'vasp_process.py'" -- and every
    # command routed through here died.  workflow.py and oqmd_extract.py were
    # converted; these four call sites were missed.
    subprocess.run([sys.executable, "-m", "htesp.vasp_process", what],
                   check=True, cwd=cwd)


def force_isif(incar_path, value=2):
    """Rewrite ``incar_path`` so that it carries exactly one ``ISIF`` entry.

    Replaces ``sed -i '/ISIF/d' INCAR`` followed by ``echo >>``, which matched
    ISIF anywhere on a line instead of as the key.
    """
    try:
        with open(incar_path, "r") as handle:
            lines = handle.read().splitlines()
    except OSError:
        lines = []
    lines = [ln for ln in lines if ln.split("=")[0].strip().upper() != "ISIF"]
    lines.append("ISIF = {}".format(value))
    with open(incar_path, "w") as handle:
        handle.write("\n".join(lines) + "\n")


def qe_stress(scf_out):
    """Return the last stress tensor of a QE run, in GPa.

    pw.x prints ``total   stress`` followed by three rows of six numbers: the
    first three columns are the tensor in Ry/bohr^3 and the last three the same
    tensor in kbar.

    Parameters:
    -----------
    scf_out : str
        Path of the pw.x output file.

    Returns:
    --------
    numpy.ndarray
        3x3 stress tensor in GPa, in the sign convention pymatgen's
        ``ElasticTensor`` expects.
    """
    with open(scf_out, "r") as handle:
        lines = handle.read().splitlines()
    starts = [i for i, ln in enumerate(lines) if "total   stress" in ln]
    if not starts:
        raise ValueError("no 'total   stress' block in {}".format(scf_out))
    rows = []
    for line in lines[starts[-1] + 1:starts[-1] + 4]:
        fields = line.split()
        if len(fields) < 6:
            raise ValueError(
                "malformed stress block in {}: {!r}".format(scf_out, line))
        rows.append([float(value) for value in fields[3:6]])
    # FIX(11): the old code took columns [:, :3] -- the Ry/bohr^3 half -- and
    # multiplied by 21798.7, the Ry/angstrom^3 -> kbar factor, wrong by
    # (bohr/angstrom)^3 = 6.75, and it never applied the sign flip the VASP
    # branch applies.  Columns 3:6 already hold kbar; convert those to GPa with
    # the same sign convention as VASP.
    return -KBAR_TO_GPA * np.array(rows)


def deformation(mpid,obj,dft,orig_prefix,deformed_struc):
    """
    Function to create deformed structures utilizing pymatgen.analysis.elastic class.

    Parameters:
    -----------
    mpid : str
        Materials id.

    obj : object
        Object of MpConnect class.

    dft : str
        Density Functional Theory (DFT) method used, e.g., 'vasp', 'qe'.

    orig_prefix : str
        Prefix for the original undeformed structure.

    deformed_struc : object
        Object containing deformed structures.

    Returns:
    --------
    None

    Example:
    --------
    >>> from mpconnect import MpConnect
    >>> from pymatgen import Structure
    >>> from pymatgen.io.vasp import Poscar
    >>> obj = MpConnect()
    >>> mpid = "mp-1234"
    >>> orig_prefix = "Si2"
    >>> deformed_struc = ...  # Object containing deformed structures
    >>> # obtained with pymatgen.analysis.elasticity.DeformedStructureSet
    >>> deformation(mpid, obj, "vasp", orig_prefix, deformed_struc)
    """
    input_data = config()
    # Extracting deformed structures
    list_str = deformed_struc.deformed_structures
    os.makedirs("scf_dir", exist_ok=True)
    # Initialize index when mpid-deformed.in file not found
    if not os.path.isfile('mpid-deformed.in'):
        entry = 0
    else:
        with open('mpid-deformed.in', 'r') as mpid_read:
            lines = mpid_read.readlines()
        entry = len(lines)
    # Loop over deformed structures
    for i,struc in enumerate(list_str):
        if dft in ('vasp', 'VASP'):
            # Write VASP input files inside R{mpid}-{i+1}-{name}/relax
            obj.prefix = struc.composition.alphabetical_formula.replace(" ","")
            poscar = Poscar(structure=struc,comment=obj.prefix)
            relax_dir = "R{}-{}-{}/relax".format(mpid,i+1,obj.prefix)
            os.makedirs(relax_dir, exist_ok=True)
            shutil.copy("R{}-{}/relax/INCAR".format(mpid,orig_prefix), relax_dir)
            poscar.write_file(filename="{}/POSCAR".format(relax_dir))
            dwn = input_data['download']
            evenkpt = dwn['inp']['evenkpt']
            kptden = input_data['kptden']
            if evenkpt:
                print("Even kpoint mesh is utilized\n")
                pos_to_kpt("{}/POSCAR".format(relax_dir),kptden,True)
            else:
                pos_to_kpt("{}/POSCAR".format(relax_dir),kptden)
            shutil.move("KPOINTS", "{}/KPOINTS".format(relax_dir))
            structure_file = structure.Structure.from_file("{}/POSCAR".format(relax_dir))
            relax_set = MPRelaxSet(structure=structure_file)
            #relax_set.potcar.write_file("{}/POTCAR".format(relax_dir))
            if os.path.isfile("vasp.in"):
                shutil.copy("vasp.in", relax_dir)
            if os.path.isfile("config.json"):
                shutil.copy("config.json", relax_dir)
            # poscar2potcar() works on the current directory, so this one
            # chdir stays -- but inside try/finally, so an exception can no
            # longer leave the process inside the deformed directory.
            pwd = os.getcwd()
            try:
                os.chdir(relax_dir)
                # Write POTCAR from POSCAR
                poscar2potcar()
                # process INCAR file with vasp.in file
                run_vasp_process("POSCAR")
                force_isif("INCAR", 2)
            finally:
                os.chdir(pwd)
        else:
            # Read magnetic flag from pwscf_in dictionary
            magnetic = input_data['pwscf_in']['magnetic']
            # Create input files with FM ordering if magnetic flag is true
            if magnetic:
                default_magmoms = input_data['magmom']['magmom']
                struc.add_spin_by_element(default_magmoms)
                obj.structure = struc
            else:
                obj.structure = struc
            # Create scf.in file
            comp_list = []
            for composition in obj.structure.composition.elements:
                comp_list.append(composition.name)
            obj.comp_list = comp_list
            # FIX(12): getkpt() and setting_qeinput() both defaulted to
            # primitive=True, so the deformed cell was pushed back through
            # get_primitive_standard_structure().  A 0.5-1% strain sits well
            # inside the symmetry tolerance and was snapped away (and the frame
            # rotated), so every deformation wrote the same undeformed input.
            # Write the deformed structure as it is, like the VASP branch does.
            obj.getkpt(primitive=False)
            evenkpt = input_data['download']['inp']['evenkpt']
            if evenkpt:
                print("Utilizing even kpoint mesh\n")
                obj.getevenkpt()
            obj.maxecut_sssp_for_subs()
            obj.prefix = struc.composition.alphabetical_formula.replace(" ","")
            if magnetic:
                obj.setting_qeinput(calculation='relax',magnetic=True,pseudo_dir='../../pp',primitive=False)
            else:
                obj.setting_qeinput(calculation='relax',pseudo_dir='../../pp',primitive=False)
            shutil.move("scf-{}.in".format(obj.mpid),
                        "scf_dir/scf-{}-{}.in".format(mpid,i+1))
        print(f"Deformation: {i+1} ",mpid,obj.prefix)
        with open("mpid-deformed.in", "a") as mpid_append:
            mpid_append.write("v{}".format(entry+1+i) + " " + mpid + "-{}".format(i+1) + " " + obj.prefix + "\n")


def process_material(mpid,orig_prefix,mode):
    """
    Handle deformation and elastic computations for one material.

    Parameters:
    -----------
    mpid : str
        Materials id.

    orig_prefix : str
        Prefix for the undeformed structure.

    mode : str
        Either ``input`` (write the deformed inputs) or ``compute_elastic``
        (read the finished runs back and fit the elastic tensor).

    Returns:
    --------
    None
    """
    input_data = config()
    dft = input_data['download']['inp']['calc']
    obj = MpConnect()
    # Determine strain from the configuration
    strain = input_data.get('strain') or [-0.01,-0.005,0.005,0.01]
    try:
        # Try loading the structure from relaxed POSCAR
        obj.structure = structure.Structure.from_file("R{}-{}/relax/POSCAR".format(mpid,orig_prefix))
    except (OSError, ValueError):
        # Load structure from SCF input file if relaxed POSCAR doesn't exist
        obj.structure = pwscf.PWInput.from_file(f"scf_dir/scf-relax-{mpid}-{orig_prefix}.in").structure
    # Get conventional standardized structure
    structure_sym = SpacegroupAnalyzer(obj.structure, symprec=0.1).get_conventional_standard_structure()
    prefix_conv = structure_sym.composition.alphabetical_formula.replace(" ","")
    print(prefix_conv)
    # Generate deformed structure set
    deformed_struc = DeformedStructureSet(structure_sym,norm_strains=strain)
    # Handle different modes of operation
    if mode == "input":
        # Prepare inputs and submit calculations
        deformation(mpid,obj,dft,orig_prefix,deformed_struc)
    elif mode == "compute_elastic":
        # Compute elastic constants
        nstruc = 6*len(strain)
        deformation_mat = deformed_struc.deformations
        deform_list = []
        stress = []
        # Loop over inputs and extract stress from deformed structures after relaxation
        for istruc in range(nstruc):
            if dft in ('VASP','vasp'):
                data = Vasprun("R{}-{}-{}/relax/vasprun.xml".format(mpid,istruc+1,prefix_conv))
                # FIX(10): ionic_steps[0] is the *first*, unrelaxed step, and
                # the factor must be -0.1 (kBar -> GPa), not -1.0.
                ionic_steps = data.as_dict()['output']['ionic_steps']
                calc_stress = -KBAR_TO_GPA*np.array(ionic_steps[-1]['stress'])
            else:
                calc_stress = qe_stress(
                    "R{}-{}-{}/relax/scf.out".format(mpid,istruc+1,prefix_conv))
            # Creates stress object
            stress_obj = Stress(calc_stress)
            pk2stress = stress_obj.piola_kirchoff_2(deformation_mat[istruc])
            stress.append(pk2stress)
            # Extracts green lagrange strain
            deform_list.append(deformation_mat[istruc].green_lagrange_strain)
        stress = np.array(stress)
        strain = np.array(deform_list)
        # Fitting to obtain elastic tensor
        elastic_tens = diff_fit(strain,stress,order=2)[0]
        elastic_obj = ElasticTensor(elastic_tens)
        properties = elastic_obj.get_structure_property_dict(structure_sym)
        print(properties.keys())
        # Save elastic properties in elastic.csv file, in PROPNAME order so the
        # row always lines up with the header written by main()
        with open("elastic.csv","a") as write_elastic:
            write_elastic.write("{},{}".format(mpid,orig_prefix))
            for prop in PROPNAME:
                value = properties.get(prop)
                try:
                    write_elastic.write("," + str(round(float(value),2)))
                except (TypeError, ValueError):
                    write_elastic.write(",")
            write_elastic.write("\n")
    else:
        print("Only 2 mode is available, either input or compute_elastic\n")


def main(argv=None):
    """
    Entry point: ``elastic.py input`` or ``elastic.py compute_elastic``.

    Parameters:
    -----------
    argv : list, optional
        Command-line arguments without the program name.  Defaults to
        ``sys.argv[1:]``; the first token is the mode.

    Returns:
    --------
    None
    """
    # FIX(14): the driver used to live in ``if __name__ == "__main__":`` with
    # its own copy of the input.in parser, so htesp.workflow could not call it
    # (run_helper imports the module and calls main()).
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print("Only 2 mode is available, either input or compute_elastic\n")
        return
    mode = argv[0]
    settings = InputIn.load("input.in", config=config())
    with open(settings.track,'r') as read_mpid:
        lines = read_mpid.readlines()
    lines = lines[settings.start-1:settings.end-1]
    if not os.path.isfile("elastic.csv"):
        with open("elastic.csv","w") as write_elastic1:
            write_elastic1.write("materials_id,compound")
            for prop1 in PROPNAME:
                write_elastic1.write("," + prop1)
            write_elastic1.write("\n")
    for line in lines:
        fields = line.split()
        if len(fields) < 3:
            warnings.warn("ignoring malformed line {!r} in {}".format(
                line.strip(), settings.track), RuntimeWarning)
            continue
        process_material(fields[1],fields[2],mode)


if __name__ == "__main__":
    main()
