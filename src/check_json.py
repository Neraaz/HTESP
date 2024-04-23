#!/usr/bin/env python
#"""Written by Niraj K. Nepal, PhD."""
"""Loading config.json file"""
import os
import json

def config():
    """
    Read data from the config.json file.

    Returns:
    - input_data (dict): Dictionary containing data from the config.json file.
    """
    try:
        pwd = os.getcwd()
        if os.path.isfile(pwd + "/config.json"):
            json_file = pwd + "/config.json"
        else:
            json_file = "../../config.json"
        with open(json_file, "r") as readjson:
            input_data = json.load(readjson)
    except FileNotFoundError:
        print("config.json file not found\n")
        return None
    return input_data
