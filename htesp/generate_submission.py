#!/usr/bin/env python
#"""Written by Niraj K. Nepal, Ph.D."""
"""Module to prepare job submission scripts"""
import os
import glob
import shutil
from htesp.ifermi_plot import ifermi
# FIX(22): this module used to `json.load(open("./config.json"))` directly,
# which bypassed the search-up-the-tree / packaged-default / $HTESP_CONFIG
# merge every other module gets from htesp.check_json.config().
from htesp.check_json import config

# FIX(21): `main` dispatched on the exact strings 'qe'/'QE', 'epw'/'EPW',
# 'wannier' (only), 'vasp'/'VASP' while the generator also accepted 'W' and
# 'WANNIER'.  A config saying which_calc: "WANNIER" therefore fell through
# every branch, `submission_files` stayed unbound and the next loop raised
# NameError.  which_calc is now normalised ONCE, here, and both the generator
# and main dispatch on the normalised value.
_CALC_ALIASES = {
    "qe": "qe",
    "epw": "epw",
    "wannier": "wannier", "w": "wannier", "wan": "wannier",
    "vasp": "vasp",
}


def normalise_calc(which_calc):
    """
    Map a `which_calc` value onto one of 'qe', 'epw', 'wannier', 'vasp'.

    Parameters
    ----------
    which_calc : str
        Value of the ``job_script.which_calc`` key, in any capitalisation.

    Returns
    -------
    str or None
        The canonical name, or None when the value is not recognised.
    """
    return _CALC_ALIASES.get(str(which_calc).strip().lower())


# FIX: every run line was built as `{parallel_command} -np {nproc} ...`.
# `-np` is mpirun syntax.  `srun` spells the same thing `-n` and rejects
# `-np`; `ibrun` takes no process count at all -- it runs the whole SLURM
# allocation -- and rejects both.  So on any site whose launcher is srun or
# ibrun (most TACC and many Cray machines) every generated run-*.sh failed at
# the first line, with a launcher usage message rather than anything from the
# DFT code.  The count flag now follows the launcher.
#: launcher basename -> the flag it spells "this many processes" with.
#: A launcher absent from this table keeps the historical ``-np``, which is
#: right for mpirun-family wrappers and is what sites have been using.
_PROCESS_FLAG = {
    "mpirun": "-np",
    "mpiexec": "-np",
    "mpiexec.hydra": "-np",
    "orterun": "-np",
    "srun": "-n",
    "aprun": "-n",
    "ibrun": None,          # no count: ibrun runs the whole allocation
}


def launch(parallel_command, nproc):
    """
    Build the launcher prefix for one run line, e.g. ``"mpirun -np 4"``.

    Parameters
    ----------
    parallel_command : str
        ``job_script.parallel_command``.  Either a bare launcher (``srun``)
        or one that already carries its own flags (``"srun --cpu-bind=cores
        -n 8"``), in which case it is passed through untouched -- a site that
        has spelled out its launcher means it.
    nproc : int or str
        ``job_script.nproc``.

    Returns
    -------
    str
        The prefix, with no trailing space, or ``""`` when no launcher is
        configured (the executable then runs serially rather than being
        handed to a command named "").
    """
    command = str(parallel_command or "").strip()
    if not command:
        return ""
    tokens = command.split()
    if any(token.startswith("-") for token in tokens):
        return command                      # already spelled out by the site
    flag = _PROCESS_FLAG.get(os.path.basename(tokens[0]).lower(), "-np")
    if flag is None:
        return command
    return "{} {} {}".format(command, flag, nproc)


def generate_submission_files(which_calc, parallel_command, nproc, command_list):
    """
    Generate submission files based on the calculation type, parallel command, number of processors,
    and a dictionary of commands.
    Parameters:
    - which_calc (str): The type of calculation.
      Valid options include 'qe', 'epw', 'wannier', and 'vasp'.
    - parallel_command (str): The parallel command to be used for execution, e.g., 'mpirun'.
    - nproc (str): The number of processors to be used.
    - command_list (dict): A dictionary containing command information for each calculation type.
      Keys are command names, and values are tuples containing the executable command,
      input file, and output file.

    Returns:
    - submission_files (dict): A dictionary containing generated submission files.
      Keys are command names, and values are formatted submission commands ready for execution.

    Raises:
    - ValueError: If an invalid calculation type is provided.

    Example:
    >>> command_list = {
    ...     'scf': ('pw.x', 'scf.in', 'scf.out'),
    ...     'nscf': ('pw.x', 'nscf.in', 'nscf.out')
    ... }
    >>> submission_files = generate_submission_files('qe', 'mpirun', '4', command_list)
    >>> print(submission_files)
    """
    npscf = int(nproc)
    run = launch(parallel_command, npscf)
    submission_files = {}

    # FIX(21): one normalisation point for every accepted spelling.
    which_calc = normalise_calc(which_calc)
    # Generate submission files based on the calculation type
    if which_calc == 'qe':
        for command_name, command_info in command_list.items():
            x_command, input_file, output_file = command_info
            submission_files[command_name] = f"{run} {x_command} -in {input_file} > {output_file}".lstrip()
            if x_command == 'bands.x':
                submission_files[command_name] += "\n"
                submission_files[command_name] += f"{run} projwfc.x -in projwfc.in > projwfc.out".lstrip()

    elif which_calc == 'epw':
        for command_name, command_info in command_list.items():
            x_command, input_file, output_file = command_info
            if command_name == 'epw':
                submission_files[command_name] = f"{run} {x_command} -npools {nproc} -i {input_file} > {output_file}".lstrip()
            else:
                submission_files[command_name] = f"{run} {x_command} -in {input_file} > {output_file}".lstrip()

    elif which_calc == 'wannier':
        for command_name, command_info in command_list.items():
            x_command, input_file, output_file = command_info
            if command_name in ('scf','nscf'):
                submission_files[command_name] = f"{run} {x_command} -in {input_file} > {output_file}".lstrip()
            elif command_name == 'pw2wannier90':
                submission_files[command_name] = f"{run} {x_command} -in {input_file} > {output_file}".lstrip()
            elif command_name == 'wannier_prepare':
                submission_files[command_name] = "wannier90.x -pp ex"
            elif command_name == 'wannier_band':
                submission_files[command_name] = "wannier90.x ex"
            else:
                print("command not found\n")

    elif which_calc == 'vasp':
        for command_name, command_info in command_list.items():
            x_command, input_file, output_file = command_info
            if command_name == 'ifermi' and os.path.isfile("ifermi.json"):
                    submission_files[command_name] = ifermi("plot","ifermi.json")
            elif command_name == 'wannier':
                submission_files[command_name] = "wannier90.x wannier90"
            else:
                submission_files[command_name] = f"{run} {x_command}".lstrip()
    else:
        print("Not valid options\n")

    return submission_files


def main():
    """
    Load parameters from a JSON file, generate submission files based on the specified calculation type,
    and create a script containing the submission commands.

    Reads parameters from 'config.json' and defines command dictionaries for different calculation types.
    Generates submission files based on the provided calculation type and prints the submission commands.
    Creates a script containing the submission commands named 'run-<calculation_type>.sh'.

    Parameters:
    None

    Returns:
    None
    """
    # FIX(22): was a direct json.load("./config.json") guarded by isfile, so
    # `parameters` stayed unbound whenever config.json was not in the cwd.
    # htesp.check_json.config() searches cwd + parents + $HTESP_CONFIG and
    # deep-merges the packaged default, so it always returns a full dict.
    parameters = config()['job_script']
    # Commands under different processes
    qe_elph_commands = {
        'scf': ('pw.x', 'scf.in', 'scf.out'),
        'band': ('pw.x', 'scf-band.in', 'scf-band.out'),
        'bandp': ('bands.x', 'band.in', 'band.out'),
        'dos': ('pw.x', 'scf-dos.in', 'scf-dos.out'),
        'dosp': ('dos.x', 'dos.in', 'dos.out'),
        'elph': ('ph.x', 'elph.in', 'elph.out'),
        'q2r': ('q2r.x', 'q2r.in', 'q2r.out'),
        'dynmat': ('dynmat.x', 'dynmat.in', 'dynmat.out'),
        'matdyn': ('matdyn.x', 'matdyn.in', 'matdyn.out'),
        'matdyn-dos': ('matdyn.x', 'matdyn-dos.in', 'matdyn-dos.out'),
        'lambda': ('lambda.x', 'lambda.in', 'lambda.out'),
        'pdos': ('projwfc.x', 'pdos.in', 'pdos.out')
    }

    epw_elph_commands = {
        'scf': ('pw.x', 'scf.in', 'scf.out'),
        'ph': ('ph.x', 'elph.in', 'elph.out'),
        'proj': ('pw.x', 'nscf-proj.in', 'nscf-proj.out'),
        'epw_nscf': ('pw.x', 'nscf_epw.in', 'nscf_epw.out'),
        'epw': ('epw.x', 'epw.in', 'epw.out')
    }

    wannier_band_commands = {
        'scf': ('pw.x', 'scf.in', 'scf.out'),
        'nscf': ('pw.x', 'nscf.in', 'nscf.out'),
        'wannier_prepare': ('wannier90.x','',''),
        'pw2wannier90': ('pw2wannier90.x', 'pw2wan.in', 'pw2wan.out'),
        'wannier_band': ('wannier90.x', '','')
    }

    vasp_commands = {
        'vasp': ('vasp_std', '', ''),
        'ifermi': ('ifermi', '', ''),
        'wannier': ('wannier','','')
    }


    # FIX(21): dispatch on the single normalised value so every accepted
    # spelling ('WANNIER', 'W', ' wannier ', ...) reaches a branch and
    # `submission_files` can never stay unbound.
    which_calc = normalise_calc(parameters['which_calc'])
    command_sets = {
        'qe': qe_elph_commands,
        'epw': epw_elph_commands,
        'wannier': wannier_band_commands,
        'vasp': vasp_commands,
    }
    if which_calc not in command_sets:
        raise ValueError(
            "job_script.which_calc = {!r} is not one of 'qe', 'epw', "
            "'wannier', 'vasp'".format(parameters['which_calc']))
    submission_files = generate_submission_files(
        which_calc, parameters['parallel_command'], parameters['nproc'],
        command_sets[which_calc])

    # Print the submission files dictionary directly
    # FIX: the cp/echo/mv shell calls interpolated the batch path and the
    # generated command line into a shell string; done with shutil/open now.
    commands = parameters['command_list']
    batch = parameters['batch']
    if parameters['command_combine'] is True:
        if not os.path.isfile("temp"):
            shutil.copy(batch, "temp")
        with open("temp", "a") as write_temp:
            for command in commands:
                write_temp.write(submission_files[command] + "\n")
        shutil.move("temp", "run-{}.sh".format(commands[-1]))
    else:
        for command in commands:
            shutil.copy(batch, "temp")
            with open("temp", "a") as write_temp:
                write_temp.write(submission_files[command] + "\n")
            shutil.move("temp", "run-{}.sh".format(command))

    calc_visible = parameters['calc_visible_with']
    # FIX: os.system("rm ...")/os.system("touch ...") replaced by os.remove and
    # a pure-Python touch.
    for stale in glob.glob("CALC_VISIBLE_WITH*"):
        os.remove(stale)
    marker = {"id": "CALC_VISIBLE_WITH_ID",
              "name": "CALC_VISIBLE_WITH_NAME",
              "id-name": "CALC_VISIBLE_WITH_ID-NAME"}.get(calc_visible)
    if marker:
        with open(marker, "a", encoding="utf-8"):
            pass


if __name__ == "__main__":
    main()
