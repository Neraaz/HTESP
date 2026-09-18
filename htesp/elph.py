#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to prepare QE el-ph calculations"""
import sys
import os
import glob
import re
import numpy as np
from pymatgen.io.pwscf import PWInput
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from htesp.check_json import config

# FIX(24): the number of irreducible representations used to come out of
# `grep "irreducible representations" elph.out > irr.out` followed by
# `np.genfromtxt("irr.out")[:, 2]` - a fixed-name scratch file plus a fixed
# column index, which silently yields NaN if QE changes the wording or the
# line gains/loses a leading word.
_NIRR = re.compile(r"(\d+)\s+irreducible\s+representation")


def read_qpoint_mesh(candidates=("qpoint.dat", "qpoint.in")):
    """
    Read the q-mesh from the first of `candidates` that exists.

    FIX(24): both callers used to leave `qpt` unbound when no q-point file was
    present, so the next line raised NameError instead of naming what was
    missing.

    Parameters
    ----------
    candidates : sequence of str
        File names to try, in order.

    Returns
    -------
    numpy.ndarray
        The three q-mesh divisions.
    """
    for name in candidates:
        if os.path.isfile(name):
            return np.loadtxt(name)
    raise FileNotFoundError(
        "no q-point file found (looked for {})".format(", ".join(candidates)))


def count_irr(elph_out):
    """
    Count the irreducible representations QE reported for each q point.

    Parameters
    ----------
    elph_out : str
        Path to the ph.x output, e.g. ``R<mpid>-<comp>/calc/elph.out``.

    Returns
    -------
    list of int
        One entry per q point, in the order QE printed them.
    """
    with open(elph_out, "r") as read_out:
        return [int(match.group(1))
                for match in (_NIRR.search(line) for line in read_out)
                if match]


def elph_in(mpid,comp,prefix,elphmode="serial"):
    """

    Default filename: 'elph-mpid-compound.in'
    The file is saved inside the 'elph_dir' folder.
    These scripts are executed within the 'create-inputs' bash script.

    The function then generates the input file for electron-phonon coupling calculation.

    Parameters:
    -----------------
    mpid: (str) Materials IDs
    comp: (str) Compounds name
    prefix: (str) prefix to use in elph.in file

    Returns:
    -----------------
    None
    """
    # FIX(24): `mass` and `nion` used to stay unbound when mass.dat was absent,
    # so every write_file() call below raised NameError instead of simply
    # omitting the amass() lines.
    mass = None
    nion = 1
    mass_flag = False
    if os.path.isfile("mass.dat"):
        mass=np.loadtxt('mass.dat')
        if mass.ndim > 0:
            nion = mass.shape[0]
        else:
            nion = 1
            mass = [mass]
        mass_flag = True
    if elphmode == "serial":
        print("----Serial mode----\n")
        output=f"elph-{mpid}-{comp}.in"
        write_file(prefix,output,mass=mass,nion=nion,mass_flag=mass_flag)
    elif elphmode == "only_init":
        print("----only_init mode----\n")
        # FIX: os.system(f"rm {file_}") replaced by os.remove.
        for file_ in glob.glob("elph_dir/PARALLEL*"):
            os.remove(file_)
        output=f"elph-{mpid}-{comp}.in"
        write_file(prefix,output,mass=mass,nion=nion,mass_flag=mass_flag,only_init=True)
    elif elphmode == "parallel_q":
        print("q-point parallelization is on\n")
        len_q = irr_q(mpid,comp)
        os.makedirs("elph_dir", exist_ok=True)
        # FIX: os.system touch+echo replaced by a plain file write.
        if not os.path.isfile(f"elph_dir/PARALLEL_q-{mpid}-{comp}"):
            with open(f"elph_dir/PARALLEL_q-{mpid}-{comp}", "w") as write_par:
                write_par.write(f"{len_q}\n")
        for i in range(len_q):
            output = f"elph_dir/elph-{mpid}-{comp}-{i+1}.in"
            iq = i+1
            write_file(prefix,output,mass=mass,nion=nion,mass_flag=mass_flag,parallel_q_flag=True,i_q = iq)
    elif elphmode == "parallel_irr":
        len_q = irr_q(mpid,comp)
        print("Peform elph calculations with only_init keyword\n")
        print("parallel over q point and irreducible representations\n")
        os.makedirs("elph_dir", exist_ok=True)
        # FIX(24): grep into the fixed-name scratch file 'irr.out' and
        # np.genfromtxt(...)[:, 2] replaced by an explicit parse of elph.out;
        # the rm/touch/echo shell calls are plain file operations now.
        elph_out = f"R{mpid}-{comp}/calc/elph.out"
        if not os.path.isfile(elph_out):
            raise FileNotFoundError(
                f"{elph_out} not found; run the only_init el-ph step first")
        array_irr = count_irr(elph_out)
        if len(array_irr) < len_q:
            raise ValueError(
                f"{elph_out} reports {len(array_irr)} irreducible "
                f"representation counts but {len_q} q points are expected")
        with open(f"elph_dir/PARALLEL_irr-{mpid}-{comp}", "w") as write_par:
            for i in range(len_q):
                len_irr = int(array_irr[i])
                write_par.write(f"v{i+1} {len_irr}\n")
                for j in range(len_irr):
                    i_q = int(i+1)
                    i_irr = int(j+1)
                    output = f"elph_dir/elph-{mpid}-{comp}-{i_q}-{i_irr}.in"
                    write_file(prefix,output,mass=mass,nion=nion,mass_flag=mass_flag,i_q=i_q,parallel_ir_flag=True,i_irr=i_irr)
    else:
        print("elphmode available options are only_init, serial, parallel_q, and parallel_irr\n")

def write_file(prefix,output,mass=None,nion=1,mass_flag=False,parallel_q_flag=False,i_q = 1,parallel_ir_flag=False,i_irr=1,only_init=False):
    """
    Write input file for electron-phonon coupling calculation.

    Args:
        prefix (str): Prefix for output files.
        output (str): Output file path.
        mass (list): List of atomic masses. Default is None.
        nion (int): Number of ions. Default is 1.
        mass_flag (bool): Flag to include atomic masses. Default is False.
        parallel_q_flag (bool): Flag to enable parallel calculations
        over q points. Default is False.
        i_q (int): Starting q point index. Default is 1.
        parallel_ir_flag (bool): Flag to enable parallel calculations over q points
        and irreducible representations. Default is False.
        i_irr (int): Starting irreducible representation index. Default is 1.
        only_init (bool): Flag to include only initialization options. Default is False.
    """
    qpt = read_qpoint_mesh()
    print("q-mesh: {}".format(qpt))
    dynmat = prefix.replace("'", "") + ".dyn"
    with open(output, 'w') as elph_write:
        elph_write.write("electron phonon coupling \n")
        elph_write.write("&inputph" + "\n")
        elph_write.write("tr2_ph=1.0d-16," + "\n")
        elph_write.write("prefix={},".format(prefix) + "\n")
        elph_write.write("fildvscf='aldv'," + "\n")
        if only_init:
            elph_write.write("only_init = .true.\n")
            elph_write.write("lqdir = .true.\n")
            elph_write.write("outdir='./'," + "\n")
        elif parallel_q_flag:
            elph_write.write(f"start_q={i_q}," + "\n")
            elph_write.write(f"last_q={i_q}," + "\n")
            elph_write.write("outdir='./'," + "\n")
        elif parallel_ir_flag:
            elph_write.write("recover = .true.\n")
            elph_write.write(f"start_q={i_q}," + "\n")
            elph_write.write(f"last_q={i_q}," + "\n")
            elph_write.write(f"start_irr={i_irr}," + "\n")
            elph_write.write(f"last_irr={i_irr}," + "\n")
            elph_write.write(f"outdir='./{i_q}-{i_irr}'," + "\n")
            elph_write.write("low_directory_check=.true.,\n")
            elph_write.write("lqdir = .true.\n")
        else:
            elph_write.write("outdir='./'," + "\n")
        if mass_flag:
            for i in range(1,nion+1):
                elph_write.write("amass({})={},".format(i,mass[i-1]) + "\n")
        elph_write.write("fildyn='{}',".format(dynmat) + "\n")
        elph_write.write("electron_phonon='interpolated'," + "\n")
        elph_write.write("el_ph_sigma=0.005," + "\n")
        elph_write.write("el_ph_nsigma=10," + "\n")
        elph_write.write("trans=.true.," + "\n")
        elph_write.write("ldisp=.true." + "\n")
        elph_write.write("nq1={},nq2={},nq3={}".format(int(qpt[0]),int(qpt[1]),int(qpt[2])) + "\n")
        elph_write.write("/" + "\n")

def irr_q(mpid,comp):
    """
    Function to write irreducible qpoint in high_symm.in file,
    necessary for computing character table with phonopy
    for `mainprogram phono5` command.
    """
    qpt = read_qpoint_mesh()
    struc = PWInput.from_file(f"R{mpid}-{comp}/relax/scf.in").structure
    irrep = SpacegroupAnalyzer(struc,symprec=0.1)
    irrep = irrep.get_ir_reciprocal_mesh((int(qpt[0]),int(qpt[1]),int(qpt[2])))
    os.makedirs("elph_dir", exist_ok=True)
    with open(f"elph_dir/high_symm-{mpid}-{comp}.in", "w") as write_high:
        for i in range(len(irrep)):
            write_high.write(str(irrep[i][0][0]) + " ")
            write_high.write(str(irrep[i][0][1]) + " ")
            write_high.write(str(irrep[i][0][2]) + "\n")
    return len(irrep)

def parse_symmetry_analysis(file_name):
    """
    Parse symmetry analysis from symmetry_analysis.in file.
    from mainprogram phono5

    Args:
        file_name (str): The name of the input file containing symmetry analysis data.

    Returns:
        dict: A dictionary containing parsed symmetry analysis data.
    """
    # Read lines from the input file
    with open(file_name, "r") as f:
        lines = f.readlines()

    # Extract q-points, character indices, and summary indices
    #qpoints = [line for line in lines if "q-point" in line]
    qpoints = [line.strip().replace('[', '').replace(']', '').split() for line in lines if "q-point" in line]
    qpoints = [[float(coord) for coord in point[1:]] for point in qpoints]
    character_indices = [i for i, line in enumerate(lines) if "Character" in line]
    summary_indices = [i for i, line in enumerate(lines) if 'Summary' in line]

    # Parse symmetry analysis data
    nirr_dict = {}
    for i, (start_index, end_index) in enumerate(zip(character_indices, summary_indices)):
        char_lines = lines[start_index:end_index]
        nirr = sum(1 for line in char_lines if ":" in line)
        # FIX(24): `qpoints[i]` is a list of floats, so `.strip()` raised
        # AttributeError - this function could never return.
        qpt = qpoints[i]
        key = f"q{i+1}"
        value = [qpt, nirr - 1]
        nirr_dict[key] = value

    return nirr_dict

if __name__ == "__main__":
    #Extracting irreducible representations from phonopy
    #file_name = "symmetry_analysis.in"
    #result = parse_symmetry_analysis(file_name)
    mpid1 = sys.argv[1]
    comp1 = sys.argv[2]
    prefix1=sys.argv[3]
    input_data = config()
    mode1 = input_data.get("elph_mode", "serial")
    elph_in(mpid1,comp1,prefix1,elphmode=mode1)
    irr_q(mpid1,comp1)
