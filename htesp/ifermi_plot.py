#!/usr/bin/env python
#"""Written by Niraj K. Nepal, Ph.D."""
"""Module to generate command for ifermi.

FIX(23): ``ifermi`` itself is an optional extra and is NOT imported here --
this module only assembles the ``ifermi ...`` command line from
``ifermi.json``.  :func:`check_ifermi` looks the executable up at the point of
use and says how to install it when it is missing, so importing this module
never requires the package.
"""
import os
import json
import shutil
import sys

def ifermi(command,input_file='ifermi.json'):
    """
    Construct a command string for the ifermi command based on the provided command
    and parameters from a JSON file.

    Parameters:
    - command (str): The command to execute (either 'info' or 'plot').

    Returns:
    - str: The constructed command string.

    If the 'ifermi.json' file exists, it reads the parameters from the file corresponding
    to the given command.
    If the parameter is a boolean value (True or False), it treats it as a boolean flag
    without a value in the command string.
    If the 'ifermi.json' file does not exist, it prints a message indicating the absence
    of the file.
    """
    # Read parameters from JSON file
    # ``command_params`` used to be left unbound when ifermi.json was absent,
    # so the loop below raised NameError instead of building a default command.
    command_params = {}
    if os.path.isfile(input_file):
        with open(input_file) as jsonfile:
            json_data = json.load(jsonfile)
            command_params = json_data.get(command, {})
    else:
        print("{} file not found; using bare 'ifermi {}'\n".format(input_file, command))

    # Construct command string
    command_str = f"ifermi {command}"

    for arg, value in command_params.items():
        # If the value is True or False, treat it as a boolean flag without a value
        if isinstance(value, bool):
            if value is True:
                command_str += f" {arg}"
        # Otherwise, treat it as a key-value pair
        else:
            command_str += f" {arg} {value}"

    return command_str
def check_ifermi():
    """Raise a clear message when the optional ``ifermi`` tool is absent.

    FIX(23).
    """
    if shutil.which("ifermi") is None:
        raise RuntimeError(
            "the 'ifermi' executable was not found on PATH.\n"
            "    pip install htesp[fermisurface]   (or: pip install ifermi)")
    return True
def main(command=None, input_file="ifermi.json", argv=None):
    """Print the assembled ifermi command line.

    FIX(21): a callable top-level entry point, so the module can be driven the
    same way as every other helper in this package.
    """
    if command is None:
        argv = list(sys.argv[1:] if argv is None else argv)
        command = argv[0] if argv else "plot"
        if len(argv) > 1:
            input_file = argv[1]
    print(ifermi(command, input_file))
    return 0
if __name__ == "__main__":
    sys.exit(main())
