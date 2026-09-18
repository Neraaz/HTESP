#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to write input files to process phonon bandstructure"""
import sys

def phonband_in(mpid=None, compound=None, prefix2=None, argv=None):
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
    # FIX(21): ``phonband_in`` is the entry point htesp.workflow.run_helper()
    # calls; the three values may also be passed directly.
    if mpid is None or compound is None or prefix2 is None:
        argv = list(sys.argv[1:] if argv is None else argv)
        if len(argv) < 3:
            raise SystemExit("usage: phonband.py <mpid> <compound> <prefix>")
        mpid, compound, prefix2 = argv[0], argv[1], argv[2]
    freq = prefix2.replace("'", "") + ".freq"
    with open("scf_dir/phonband-{}-{}.in".format(mpid,compound), 'w') as phon_process:
        phon_process.write(freq + "\n")
        phon_process.write("0 5000" + "\n")
        phon_process.write("freq.plot" + "\n")
        phon_process.write("freq.ps" + "\n")
        phon_process.write("0.0" + "\n")
        phon_process.write("100.0 0.0" + "\n")

main = phonband_in
if __name__ == "__main__":
    phonband_in()
