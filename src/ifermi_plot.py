#!/usr/bin/env python
"""Written by Niraj K. Nepal, Ph.D."""
import os
import json

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
    if os.path.isfile(input_file):
        with open(input_file) as jsonfile:
            json_data = json.load(jsonfile)
            command_params = json_data.get(command, {})
    else:
        print("ifermi.json file not found\n")

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
