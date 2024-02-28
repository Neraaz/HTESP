#!/usr/bin/env python
"""Written by Niraj K. Nepal, Ph.D."""
import os
import glob
import json


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
    """
    npscf = int(nproc)
    submission_files = {}

    if which_calc in ('qe','QE'):
        for command_name, command_info in command_list.items():
            x_command, input_file, output_file = command_info
            submission_files[command_name] = f"{parallel_command} -np {npscf} {x_command} < {input_file} > {output_file}"

    elif which_calc in ('epw','EPW'):
        for command_name, command_info in command_list.items():
            x_command, input_file, output_file = command_info
            if command_name == 'epw':
                submission_files[command_name] = f"{parallel_command} -np {nproc} {x_command} -npools {nproc} -i {input_file} > {output_file}"
            else:
                submission_files[command_name] = f"{parallel_command} -np {nproc} {x_command} < {input_file} > {output_file}"

    elif which_calc in ('wannier','W','WANNIER'):
        for command_name, command_info in command_list.items():
            x_command, input_file, output_file = command_info
            if command_name in ('scf','nscf'):
                submission_files[command_name] = f"{parallel_command} -np {npscf} {x_command} < {input_file} > {output_file}"
            elif command_name == 'pw2wannier90':
                submission_files[command_name] = f"{parallel_command} -np {npscf} {x_command} -in {input_file} > {output_file}"
            elif command_name == 'wannier_prepare':
                submission_files[command_name] = "wannier90.x -pp ex"
            elif command_name == 'wannier_band':
                submission_files[command_name] = "wannier90.x ex"
            else:
                print("command not found\n")

    elif which_calc in ('vasp','VASP'):
        for command_name, command_info in command_list.items():
            x_command, input_file, output_file = command_info
            if command_name == 'ifermi':
                submission_files[command_name] = "ifermi plot --property velocity --interpolation-factor 10 --property-colormap bwr"
            else:
                submission_files[command_name] = f"{parallel_command} -np {npscf} {x_command}"
    else:
        print("Not valid options\n")

    return submission_files


def main():
    """
    Load parameters from a JSON file, generate submission files based on the specified calculation type,
    and create a script containing the submission commands.

    Reads parameters from 'parameters.json' and defines command dictionaries for different calculation types.
    Generates submission files based on the provided calculation type and prints the submission commands.
    Creates a script containing the submission commands named 'run-<calculation_type>.sh'.

    Parameters:
    None

    Returns:
    None
    """
    # Load parameters from JSON file
    if os.path.isfile("htepc.json"):
        with open('htepc.json') as json_file:
            #parameters = json.load(json_file)
            parameters = json.load(json_file)['job_script']

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
        'ifermi': ('ifermi', '', '')
    }


    # Generate the submission files
    if parameters['which_calc'] in ('qe','QE'):
        submission_files = generate_submission_files(parameters['which_calc'], parameters['parallel_command'], parameters['nproc'], qe_elph_commands)
    elif parameters['which_calc'] in ('epw','EPW'):
        submission_files = generate_submission_files(parameters['which_calc'], parameters['parallel_command'], parameters['nproc'], epw_elph_commands)
    elif parameters['which_calc'] == 'wannier':
        submission_files = generate_submission_files(parameters['which_calc'], parameters['parallel_command'], parameters['nproc'], wannier_band_commands)
    elif parameters['which_calc'] in ('vasp','VASP'):
        submission_files = generate_submission_files(parameters['which_calc'], parameters['parallel_command'], parameters['nproc'], vasp_commands)
    else:
        print("Options not available\n")

    # Print the submission files dictionary directly
    commands = parameters['command_list']
    if parameters['command_combine'] is True:
        for command in commands:
            if not os.path.isfile("temp"):
                os.system("""cp {} temp""".format(parameters['batch']))
            os.system("""echo "{}" >> temp""".format(submission_files[command]))
        os.system("""mv temp run-{}.sh""".format(commands[-1]))
    else:
        for command in commands:
            os.system("""cp {} temp""".format(parameters['batch']))
            os.system("""echo "{}" >> temp""".format(submission_files[command]))
            os.system("""mv temp run-{}.sh""".format(command))

    calc_visible = parameters['calc_visible_with']
    check_calc_visible = glob.glob("CALC_VISIBLE_WITH*")
    if len(check_calc_visible) > 0:
        os.system("rm CALC_VISIBLE_WITH*")
    if calc_visible == "id":
        os.system("touch CALC_VISIBLE_WITH_ID")
    elif calc_visible == "name":
        os.system("touch CALC_VISIBLE_WITH_NAME")
    elif calc_visible == "id-name":
        os.system("touch CALC_VISIBLE_WITH_ID-NAME")
    else:
        pass


if __name__ == "__main__":
    main()
