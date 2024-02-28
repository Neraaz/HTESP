#!/usr/bin/env python
"""Writen by Niraj K. Nepal, Ph.D."""
import sys

def q2r_in():
    """
    Prepare input files for force constant in real space (q2r.in) and obtain atomic displacement
    for a particular mode (dynmat.in).

    This function is responsible for creating input files required for calculating force constants
    in real space (q2r.in) and obtaining atomic displacement for a specific mode (dynmat.in). The
    default names for the input files are q2r-mpid-compound.in and dynamt-mpid-compound.in. The
    function is typically invoked within the 'create-inputs' bash script.

    Parameters:
    None (relies on command-line arguments provided through sys.argv)

    Returns:
    None

    """
    mpid = sys.argv[1]
    comp = sys.argv[2]
    prefix2=sys.argv[3]
    dynmat = prefix2.replace("'", "") + ".dyn"
    dynmat1 = prefix2.replace("'", "") + ".dyn"
    frc = prefix2.replace("'", "") + ".fc"
    with open("q2r-{}-{}.in".format(mpid,comp), 'w') as q2r:
        q2r.write("&input" + "\n")
        q2r.write("zasr='simple'," + "\n")
        q2r.write("fildyn='{}',".format(dynmat) + "\n")
        q2r.write("flfrc='{}',".format(frc) + "\n")
        q2r.write("la2F=.true." + "\n")
        q2r.write("/" + "\n")
    with open("dynmat-{}-{}.in".format(mpid,comp), 'w') as dyn:
        dyn.write("&input" + "\n")
        dyn.write("asr='simple'," + "\n")
        dyn.write("fildyn='{}',".format(dynmat1) + "\n")
        dyn.write("/" + "\n")

if __name__ == "__main__":
    q2r_in()
