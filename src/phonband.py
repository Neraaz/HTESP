#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D."""
import sys

def phonband_in():
    """
    Prepare input file for processing phonon bandstructure.

    Usage:
    python script.py mpid compound prefix2

    Args:
    - mpid (str): Materials Project ID.
    - compound (str): Compound name.
    - prefix2 (str): Prefix string used for file naming.

    This function creates an input file named 'phonband-mpid-compound.in'
    within the 'scf_dir' directory.
    It writes the following contents into the file:
    1. The frequency file name (prefix2 without single quotes + '.freq').
    2. Frequency range: 0 to 5000.
    3. Output file for plotting the frequency: 'freq.plot'.
    4. Output postscript file for plotting: 'freq.ps'.
    5. Smearing parameter: 0.0.
    6. Broadening parameters: 100.0 and 0.0.

    This script is typically run within the 'create-inputs' bash script.
    """
    mpid = sys.argv[1]
    compound = sys.argv[2]
    prefix2=sys.argv[3]
    freq = prefix2.replace("'", "") + ".freq"
    with open("scf_dir/phonband-{}-{}.in".format(mpid,compound), 'w') as phon_process:
        phon_process.write(freq + "\n")
        phon_process.write("0 5000" + "\n")
        phon_process.write("freq.plot" + "\n")
        phon_process.write("freq.ps" + "\n")
        phon_process.write("0.0" + "\n")
        phon_process.write("100.0 0.0" + "\n")

if __name__ == "__main__":
    phonband_in()
