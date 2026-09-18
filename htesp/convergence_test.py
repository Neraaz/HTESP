#!/usr/bin/env python
#"""Written by Niraj K. Nepal, Ph.D.
"""
 Module to submit convergence tests.
 Runs with 'mainprogram convtest'"""

import contextlib
import os
import shutil
import subprocess
import sys
import time
import warnings
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.pwscf import PWOutput
from htesp.check_json import config
from htesp.write_potcar import stage_potcar
from htesp.inputin import InputIn

# FIX(20): conv_test.ecut is written verbatim as ``ecutwfc`` (Ry) for QE and as
# ``ENCUT`` (eV) for VASP, and both used to land in one ``param-energy.txt``
# with no record of which.  The numbers are still written exactly as given; the
# unit is recorded in the file header and in the collected file's name.
ECUT_UNIT = {'qe': 'Ry', 'vasp': 'eV'}


@contextlib.contextmanager
def pushd(path):
    """``cd path`` for the body, always coming back.

    FIX(22): every ``os.chdir(...) ... os.chdir("../../")`` chain in this module
    left the process in the wrong directory when anything in between raised,
    so the next material was set up inside the previous one's run directory.
    """
    previous = os.getcwd()
    os.chdir(path)
    try:
        yield path
    finally:
        os.chdir(previous)


def read_lines(path):
    """Return the lines of ``path`` without their trailing newline."""
    with open(path, "r") as handle:
        return handle.read().splitlines()


def write_lines(path, lines):
    """Write ``lines`` (no trailing newlines needed) to ``path``."""
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + ("\n" if lines else ""))


def run(command, cwd=None, check=True):
    """Run ``command`` (a list) without a shell, reporting a failure."""
    completed = subprocess.run(command, cwd=cwd, check=check)
    if completed.returncode:
        print("{} exited with {}".format(command[0], completed.returncode))
    return completed.returncode


def kpoint_vasp(kpoint,kconv=False):
    """
    Function to write kpoint.

    Parameters:
    - kpoint (list): List containing the k-point coordinates.
    - kconv (bool): True writes into R<kx>-<ky>-<kz>/KPOINTS, False rewrites
      KPOINTS in the current directory.

    Returns:
    - None

    Example:
    To write kpoints for a convergence calculation:
    >>> kpoint_vasp([2, 2, 2], kconv=True)

    To write kpoints for a regular calculation:
    >>> kpoint_vasp([4, 4, 4])
    """
    lines = read_lines("KPOINTS")
    # ``sed '4,5d'`` drops the mesh and the shift (1-based lines 4 and 5)
    head = lines[:3]
    body = ["{} {} {}".format(kpoint[0], kpoint[1], kpoint[2]), "0 0 0"]
    if kconv:
        # Generate KPOINTS with different k-point mesh in different folders
        target = "R{}-{}-{}".format(kpoint[0], kpoint[1], kpoint[2])
        os.makedirs(target, exist_ok=True)
        write_lines(os.path.join(target, "KPOINTS"), head + body)
    else:
        # Replace existing kpoints with new ones in KPOINTS file
        write_lines("KPOINTS", head + body)


def submission(param,jobscript,A,B,jj,cwd=None):
    """
    Submit a jobscript for computation based on parameters.

    Parameters:
    - param (str): Parameter type, 'ecut' or 'kpoint'.
    - jobscript (str): Name of the jobscript to be submitted.
    - A (str): Materials Id.
    - B (str): Name for the compound.
    - jj (str or tuple): The value of the parameter. If param is 'ecut', it's a string;
                          if param is 'kpoint', it's a list (e.g., (kx, ky, kz)).
    - cwd (str, optional): directory holding the jobscript. Defaults to the
      current one, which is what the callers used to arrange with os.chdir.

    Returns:
    None
    Example:

    To submit a jobscript named 'myjob.sh' for computation with parameter type 'ecut',
    Materials Id 'mp-123', compound name 'Si', and parameter value '500':
    >>> submission('ecut', 'myjob.sh', 'mp-123', 'Si', '500')
    """
    if param == 'kpoint':
        jj = f"{jj[0]}-{jj[1]}-{jj[2]}"
    elif param != 'ecut':
        print("Only kpoint and ecut are allowed\n")
    cwd = cwd or "."
    calc_visible_file = ""
    # Determine the name of file to be visible based on available options
    if os.path.isfile(os.path.join(cwd, "../../../CALC_VISIBLE_WITH_ID")):
        calc_visible_file = f"{A}-{jj}.sh"
    elif os.path.isfile(os.path.join(cwd, "../../../CALC_VISIBLE_WITH_NAME")):
        calc_visible_file = f"{B}-{jj}.sh"
    elif os.path.isfile(os.path.join(cwd, "../../../CALC_VISIBLE_WITH_ID-NAME")):
        calc_visible_file = f"{A}-{B}-{jj}.sh"
    # Move the jobscript to the visible file name and submit the job
    if calc_visible_file:
        shutil.move(os.path.join(cwd, jobscript),
                    os.path.join(cwd, calc_visible_file))
        jobscript = calc_visible_file
    # a queue that refuses the job must not abort the whole scan
    run(["sbatch", jobscript], cwd=cwd, check=False)
    time.sleep(1)


def incar_base(path, drop=("ISIF", "NSW", "ENCUT", "ENAUG"), add=(("NSW", 0),)):
    """Return the lines of ``path`` without ``drop``, plus the ``add`` pairs.

    The key is matched at the start of the line so that dropping ENCUT does not
    also drop ENCUTGW, which ``sed '/ENCUT/d'`` did.
    """
    drop = {key.upper() for key in drop}
    lines = []
    for line in read_lines(path):
        head = line.split("=", 1)[0].strip().upper() if "=" in line else ""
        if head in drop:
            continue
        lines.append(line)
    return lines + ["{} = {}".format(key, value) for key, value in add]


def main_vasp(file_name,parameter,start,end):
    """
    Execute convergence test calculations for kinetic energy cutoff.

    Parameters:
    - file_name (str): File holding materials id and name information.
    - parameter (str): Parameter to control the type of calculation.
    - start (int) : Starting index
    - end (int) : Ending index (not included)

    Returns:
    - None

    This function iterates through the specified range of indices and performs convergence test calculations.
    It modifies input files, creates directories, and submits jobs based on the provided parameter.

    Example:
    To perform a convergence test for kinetic energy cutoff ('ecut') from index 1 to 10:
    >>> main_vasp("mpid.in", "ecut", 1, 10)
    """
    input_data = config()

    print("---------------------------------------------------------------------------------------------------------------")
    print(f"Submitting convergence test calculations for {parameter}")
    print("---------------------------------------------------------------------------------------------------------------")

    for ii in range(start, end):
        with open(file_name, 'r') as file:
            for line in file:
                if f"v{ii} " not in line:
                    continue
                A = line.split()[1]
                B = line.split()[2]
                print(f"{A} {B}")
                # Create sub directories
                work = f"R{A}-{B}/{parameter}"
                os.makedirs(work, exist_ok=True)
                # Copy relevant input files.  FIX(38): the reference tree
                # cannot ship a licensed POTCAR, so build one from the POSCAR
                # when it is absent rather than dying on the copy.
                relax = f"R{A}-{B}/relax"
                # FIX(40): a POTCAR is built when it can be and reported when
                # it cannot; either way the scan carries on.  It is licensed,
                # so it is never in the reference tree and may be unavailable
                # altogether -- that is not a reason to abort input generation.
                wanted = ["INCAR", "KPOINTS", "POSCAR"]
                if stage_potcar(relax):
                    wanted.append("POTCAR")
                for name in wanted:
                    shutil.copy(f"{relax}/{name}", work)
                # Perform convergence test for kinetic energy cutoff ('ecut')
                if parameter == 'ecut':
                    ecut_list = input_data['conv_test']['ecut']
                    kpoint = input_data['conv_test']['kpoint'][0]
                    with pushd(work):
                        kpoint_vasp(kpoint,kconv=False)
                        base = incar_base("INCAR")
                        for jj in ecut_list:
                            run_dir = f"R{jj}"
                            os.makedirs(run_dir, exist_ok=True)
                            # FIX(19): the old loop moved INCAR1 into the first
                            # R<ecut> directory and then appended to a file
                            # that no longer existed, so every later cutoff got
                            # a two-line INCAR holding only its own ENCUT and
                            # ENAUG.  A complete INCAR is written every time.
                            # The value is the eV ENCUT, as given in conv_test.
                            write_lines(os.path.join(run_dir, "INCAR"),
                                        base + ["ENCUT = {}".format(jj),
                                                "ENAUG = {}".format(2*jj)])
                            for name in [n for n in ("KPOINTS", "POSCAR", "POTCAR")
                                         if os.path.isfile(n)]:
                                shutil.copy(name, run_dir)
                            shutil.copy("../../run-vasp.sh",
                                        os.path.join(run_dir, "run.sh"))
                            submission("ecut","run.sh",A,B,jj,cwd=run_dir)
                # Perform convergence test for k-points
                elif parameter == 'kpoint':
                    encut = input_data['conv_test']['ecut'][0]
                    enaug = encut * 2
                    kpoint_list = input_data['conv_test']['kpoint']
                    with pushd(work):
                        base = incar_base("INCAR") + [
                            "ENCUT = {}".format(encut),
                            "ENAUG = {}".format(enaug)]
                        for jj in kpoint_list:
                            run_dir = f"R{jj[0]}-{jj[1]}-{jj[2]}"
                            os.makedirs(run_dir, exist_ok=True)
                            # Creates KPOINTS file
                            kpoint_vasp(jj,True)
                            # Copy input files and submit job
                            shutil.copy("../../run-vasp.sh",
                                        os.path.join(run_dir, "run.sh"))
                            write_lines(os.path.join(run_dir, "INCAR"), base)
                            for name in [n for n in ("POSCAR", "POTCAR")
                                         if os.path.isfile(n)]:
                                shutil.copy(name, run_dir)
                            submission("kpoint","run.sh",A,B,jj,cwd=run_dir)
                else:
                    print("Only ecut and kpoint convergence tests can be performed\n")

    print("All done")


def kpoint_qe(kpoint,mpid):
    """
    Insert k-points into the SCF input file.

    Parameters:
    - kpoint (list): List containing the k-point coordinates.
    - mpid (str): The ID of the compound.

    Returns:
    - str: the path of the file written ('scf_dir/temp.in').

    Example:
    >>> kpoint_qe([0.5, 0.5, 0.5], "mp-12345")
    """
    source = "scf_dir/temp-{}.in".format(mpid)
    lines = read_lines(source)
    # Find the line index where K_POINTS is located (the last one, as before)
    idx = None
    for i,line in enumerate(lines):
        if 'K_POINTS' in line:
            idx = i
    if idx is None:
        raise ValueError("no K_POINTS card in {}".format(source))
    # Replace the mesh line following K_POINTS with the new one
    mesh = "{} {} {} 0 0 0".format(kpoint[0],kpoint[1],kpoint[2])
    if idx + 1 < len(lines):
        lines[idx + 1] = mesh
    else:
        lines.append(mesh)
    write_lines("scf_dir/temp.in", lines)
    return "scf_dir/temp.in"


def prepare_qe_template(source, dest):
    """Strip the cutoffs/pseudo_dir from an scf input and force ``calculation='scf'``.

    Replaces the five chained ``sed`` invocations that wrote and deleted four
    intermediate files.
    """
    drop = ('ecutrho', 'ecutwfc', 'calculation', 'pseudo_dir')
    kept = [ln for ln in read_lines(source)
            if not any(token in ln for token in drop)]
    out = []
    for line in kept:
        out.append(line)
        if '&CONTROL' in line:
            out.append("  calculation = 'scf',")
            out.append("  pseudo_dir = '../../../pp/',")
    write_lines(dest, out)
    return dest


def main_qe(file_name,parameter,start,end):
    """
    Execute convergence test calculations for kinetic energy cutoff.

    Parameters:
    - file_name (str): File holding materials id and name information.
    - parameter (str): Parameter to control the type of calculation.
    - start (int) : Starting index
    - end (int) : Ending index (not included)

    Returns:
    - list: the temporary files this run created, so the caller can remove
      exactly those.

    This function iterates through the specified range of indices and performs convergence test calculations.
    It modifies input files, creates directories, and submits jobs based on the provided parameter.

    Example:
    >>> main_qe("mpid.in", "ecut", 1, 5)
    """
    input_data = config()
    # FIX(21): ``rm scf_dir/temp*`` deleted the scratch inputs of every other
    # material (and of any other run sharing the directory).  The files this
    # run creates are collected here and only those are removed.
    created = []

    print("---------------------------------------------------------------------------------------------------------------")
    print(f"Submitting convergence test calculations for {parameter}")
    print("---------------------------------------------------------------------------------------------------------------")

    for ii in range(start, end):
        with open(file_name, 'r') as file:
            for line in file:
                if f"v{ii} " not in line:
                    continue
                A = line.split()[1]
                B = line.split()[2]
                print(f"{A} {B}")
                # Create directories if they don't exist
                os.makedirs(f"R{A}-{B}/{parameter}", exist_ok=True)
                # Updating ecutrho, ecutwfc, pseudo_dir, and calculation keywords of QE
                template = prepare_qe_template(f"scf_dir/scf-{A}.in",
                                               f"scf_dir/temp-{A}.in")
                created.append(template)
                # Perform convergence tests based on the parameter
                if parameter == 'ecut':
                    ecut_list = input_data['conv_test']['ecut']
                    kpoint = input_data['conv_test']['kpoint'][0]
                    created.append(kpoint_qe(kpoint,A))
                    # Iterate over ecut values
                    for jj in ecut_list:
                        run_dir = f"R{A}-{B}/{parameter}/R{jj}"
                        os.makedirs(run_dir, exist_ok=True)
                        kk = jj * 8
                        # FIX(20): conv_test.ecut is written straight through as
                        # ecutwfc, i.e. interpreted as Ry here and as eV by the
                        # VASP branch.  The value is deliberately not converted;
                        # the unit is recorded when the energies are collected.
                        lines = read_lines("scf_dir/temp.in")
                        out = []
                        for line in lines:
                            out.append(line)
                            if '&SYSTEM' in line:
                                out.append("  ecutwfc = {}, ecutrho = {},".format(jj, kk))
                        write_lines(os.path.join(run_dir, "scf.in"), out)
                        shutil.copy("run-scf.sh", run_dir)
                        submission('ecut',"run-scf.sh",A,B,jj,cwd=run_dir)
                elif parameter == 'kpoint':
                    ecut = input_data['conv_test']['ecut'][0]
                    ecutrho = ecut * 8
                    kpoint_list = input_data['conv_test']['kpoint']
                    lines = read_lines(template)
                    out = []
                    for line in lines:
                        out.append(line)
                        if '&SYSTEM' in line:
                            out.append("  ecutwfc = {}, ecutrho = {},".format(ecut, ecutrho))
                    write_lines(template, out)
                    # Iterate over kpoint values
                    for jj in kpoint_list:
                        run_dir = f"R{A}-{B}/{parameter}/R{jj[0]}-{jj[1]}-{jj[2]}"
                        os.makedirs(run_dir, exist_ok=True)
                        kpoint_qe(jj,A)
                        # Copy input files and submit job
                        shutil.copy("run-scf.sh", run_dir)
                        shutil.move("scf_dir/temp.in",
                                    os.path.join(run_dir, "scf.in"))
                        submission('kpoint',"run-scf.sh",A,B,jj,cwd=run_dir)
                else:
                    print("Only ecut and kpoint convergence tests can be performed\n")

    print("All done")
    return created


def extract(file_name,parameter,start,end,calc):
    """Collect the total energies of a finished convergence test.

    Parameters:
    - file_name (str): File holding materials id and name information.
    - parameter (str): 'ecut' or 'kpoint'.
    - start, end (int): range of tracking-file entries (end is exclusive).
    - calc (str): 'QE' or 'VASP'.

    Returns:
    - None

    Writes ``R<id>-<comp>/<parameter>/param-energy.txt`` for each material and
    copies it into ``convergence_result/``.
    """
    input_data = config()
    is_vasp = str(calc).lower() == 'vasp'
    dft = 'vasp' if is_vasp else 'qe'
    # FIX(20): the unit of conv_test.ecut depends on the code, and both used to
    # be written into an identically named file.
    unit = ECUT_UNIT[dft] if parameter == 'ecut' else 'grid'
    os.makedirs("convergence_result", exist_ok=True)
    with open(file_name, 'r') as read_mpid:
        lines = read_mpid.readlines()
    for line in lines[start-1:end-1]:
        fields = line.split()
        if len(fields) < 3:
            warnings.warn("ignoring malformed line {!r} in {}".format(
                line.strip(), file_name), RuntimeWarning)
            continue
        mpid, comp = fields[1], fields[2]
        folder = f"R{mpid}-{comp}/{parameter}"
        out_path = os.path.join(folder, "param-energy.txt")
        list_param = input_data['conv_test'][parameter]
        with open(out_path, "w") as write_energy:
            write_energy.write(f"# {parameter} unit: {unit}\n")
            for value in list_param:
                if parameter == 'ecut':
                    value1 = value
                else:
                    value1 = f"{value[0]}-{value[1]}-{value[2]}"
                run_dir = os.path.join(folder, f"R{value1}")
                if is_vasp:
                    data = Vasprun(os.path.join(run_dir, "vasprun.xml"))
                else:
                    data = PWOutput(os.path.join(run_dir, "scf.out"))
                energy = data.final_energy.real
                write_energy.write(f"{value1} {energy}" + "\n")
        shutil.copy(out_path,
                    f"convergence_result/{parameter}-{dft}-{mpid}-{comp}.txt")


def main(argv=None):
    """
    Entry point: ``convergence_test.py calculate`` or ``... extract``.

    Parameters:
    - argv (list, optional): command-line arguments without the program name.
      Defaults to ``sys.argv[1:]``; the first token is the mode.

    Returns:
    - None
    """
    # FIX(22): the driver used to live in ``if __name__ == "__main__":`` with
    # its own copy of the input.in parser, so htesp.workflow could not call it.
    argv = list(sys.argv[1:] if argv is None else argv)
    input_data = config()
    settings = InputIn.load("input.in", config=input_data)
    # Extract calculation type and parameter variable
    calc = input_data['download']['inp']['calc']
    parameter_variable = input_data['conv_test']['param']
    # Determine the mode of operation (calculate or extract results)
    mode = argv[0] if argv else "calculate"
    if mode == "calculate":
        if calc in ('QE','qe'):
            # Perform convergence test calculations for Quantum ESPRESSO
            created = main_qe(settings.track,parameter_variable,
                              settings.start,settings.end)
            for path in created:
                try:
                    os.remove(path)
                except FileNotFoundError:
                    continue
                except OSError as exc:
                    print("could not remove {}: {}".format(path, exc))
        else:
            # Perform convergence test calculations for VASP
            main_vasp(settings.track,parameter_variable,
                      settings.start,settings.end)
    else:
        print("Extracting results\n")
        extract(settings.track,parameter_variable,settings.start,
                settings.end,calc)


if __name__ == "__main__":
    main()
