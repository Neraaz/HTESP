#!/usr/bin/env python
"""Written by Niraj K. Nepal, Ph.D.
 Submit convergence test for QE scf calculations.
 Runs with 'mainprogram convtest'"""

import os
import sys
import json
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.io.pwscf import PWOutput
try:
    PWD = os.getcwd()
    if os.path.isfile(PWD+"/htepc.json"):
        JSONFILE = PWD+"/htepc.json"
    else:
        JSONFILE = "../../htepc.json"
    with open(JSONFILE, "r") as readjson:
        input_data = json.load(readjson)
except FileNotFoundError:
    print("htepc.json file not found\n")

def kpoint_vasp(kpoint,kconv=False):
    """
    Function to write kpoint.

    Parameters:
    - kpoint (list): List containing the k-point coordinates.
    - mpid (str): The ID of the compound.
    - compound (str): The name of the compound.

    Returns:
    - None
    """
    if kconv:
        os.system(f"""sed '4,5d' KPOINTS > R{kpoint[0]}-{kpoint[1]}-{kpoint[2]}/KPOINTS""")
        os.system(f"""echo "{kpoint[0]} {kpoint[1]} {kpoint[2]}" >> R{kpoint[0]}-{kpoint[1]}-{kpoint[2]}/KPOINTS""")
        os.system(f"""echo "0 0 0" >> R{kpoint[0]}-{kpoint[1]}-{kpoint[2]}/KPOINTS""")
    else:
        os.system("""sed -i '4,5d' KPOINTS """)
        os.system(f"""echo "{kpoint[0]} {kpoint[1]} {kpoint[2]}" >> KPOINTS""")
        os.system("""echo "0 0 0" >> KPOINTS""")

def submission(param,jobscript,A,B,jj):
    """
    Submit a jobscript for computation based on parameters.

    Parameters:
    - param (str): Parameter type, 'ecut' or 'kpoint'.
    - jobscript (str): Name of the jobscript to be submitted.
    - A (str): Materials Id.
    - B (str): Name for the compound.
    - jj (str or tuple): The value of the parameter. If param is 'ecut', it's a string;
                          if param is 'kpoint', it's a list (e.g., (kx, ky, kz)).

    Returns:
    None
    """
    if param == 'ecut':
        jj = jj
    elif param == 'kpoint':
        jj = f"{jj[0]}-{jj[1]}-{jj[2]}"
    else:
        print("Only kpoint and ecut are allowed\n")
    calc_visible_file = ""
    if os.path.isfile("../../../CALC_VISIBLE_WITH_ID"):
        calc_visible_file = f"{A}-{jj}.sh"
    elif os.path.isfile("../../../CALC_VISIBLE_WITH_NAME"):
        calc_visible_file = f"{B}-{jj}.sh"
    elif os.path.isfile("../../../CALC_VISIBLE_WITH_ID-NAME"):
        calc_visible_file = f"{A}-{B}-{jj}.sh"

    if calc_visible_file:
        os.system(f"mv {jobscript} {calc_visible_file}")
        os.system(f"sbatch {calc_visible_file}")
        os.system("sleep 1s")
    else:
        os.system("sbatch {jobscript}")
        os.system("sleep 1s")
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

    """

    print("---------------------------------------------------------------------------------------------------------------")
    print(f"Submitting convergence test calculations for {parameter}")
    print("---------------------------------------------------------------------------------------------------------------")

    for ii in range(start, end):
        with open(file_name, 'r') as file:
            for line in file:
                if f"v{ii} " in line:
                    A = line.split()[1]
                    B = line.split()[2]
                    print(f"{A} {B}")

                    if not os.path.isdir(f"R{A}-{B}"):
                        os.mkdir(f"R{A}-{B}")

                    if not os.path.isdir(f"R{A}-{B}/{parameter}"):
                        os.mkdir(f"R{A}-{B}/{parameter}")

                    # Execute sed commands
                    os.system(f"""cp R{A}-{B}/relax/INCAR R{A}-{B}/{parameter}/""")
                    os.system(f"""cp R{A}-{B}/relax/KPOINTS R{A}-{B}/{parameter}/""")
                    os.system(f"""cp R{A}-{B}/relax/POSCAR R{A}-{B}/{parameter}/""")
                    os.system(f"""cp R{A}-{B}/relax/POTCAR R{A}-{B}/{parameter}/""")
                    if parameter == 'ecut':
                        ecut_list = input_data['conv_test']['ecut']
                        kpoint = input_data['conv_test']['kpoint'][0]
                        os.chdir(f"""R{A}-{B}/{parameter}""")
                        kpoint_vasp(kpoint,kconv=False)
                        os.system("sed '/ISIF/d' INCAR | sed '/NSW/d' | sed '/ENCUT/d' | sed '/ENAUG/d' > INCAR1")
                        os.system("""echo "NSW = 0" >> INCAR1""")
                        for jj in ecut_list:
                            if not os.path.isdir(f"R{jj}"):
                                os.mkdir(f"R{jj}")
                            # Write modified scf input file
                            os.system(f"""echo "ENCUT = {jj}" >> INCAR1""")
                            enaug = 2*jj
                            os.system(f"""echo "ENAUG = {enaug}" >> INCAR1""")
                            os.system(f"""mv INCAR1 R{jj}/INCAR""")
                            os.system(f"""cp KPOINTS R{jj}/""")
                            os.system(f"""cp POSCAR R{jj}/""")
                            os.system(f"""cp POTCAR R{jj}/""")
                            os.system(f"cp ../../run-vasp.sh R{jj}/run.sh")
                            os.chdir(f"R{jj}")
                            submission("ecut","run.sh",A,B,jj)
                            os.chdir("../")
                        os.chdir("../../")
                            # Copy shell script and submit job
                    elif parameter == 'kpoint':
                        encut = input_data['conv_test']['ecut'][0]
                        enaug = encut * 2
                        os.chdir(f"""R{A}-{B}/{parameter}""")
                        os.system("sed '/ISIF/d' INCAR | sed '/NSW/d' | sed '/ENCUT/d' | sed '/ENAUG/d' > INCAR1")
                        os.system("""echo "NSW = 0" >> INCAR1""")
                        os.system(f"""echo "ENCUT = {encut}" >> INCAR1""")
                        os.system(f"""echo "ENAUG = {enaug}" >> INCAR1""")
                        kpoint_list = input_data['conv_test']['kpoint']
                        for jj in kpoint_list:
                            if not os.path.isdir(f"R{jj[0]}-{jj[1]}-{jj[2]}"):
                                os.mkdir(f"R{jj[0]}-{jj[1]}-{jj[2]}")
                            kpoint_vasp(jj,True)
                            os.system(f"cp ../../run-vasp.sh R{jj[0]}-{jj[1]}-{jj[2]}/run.sh")
                            os.system(f"cp INCAR1 R{jj[0]}-{jj[1]}-{jj[2]}/INCAR")
                            os.system(f"cp POSCAR R{jj[0]}-{jj[1]}-{jj[2]}/")
                            os.system(f"cp POTCAR R{jj[0]}-{jj[1]}-{jj[2]}/")
                            os.chdir(f"R{jj[0]}-{jj[1]}-{jj[2]}")
                            submission("kpoint","run.sh",A,B,jj)
                            os.chdir("../")
                            # Copy shell script and submit job
                        os.chdir("../../")
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
    - None
    """
    with open("scf_dir/temp-{}.in".format(mpid),"r") as read_scf:
        lines = read_scf.readlines()
    for i,line in enumerate(lines):
        if 'K_POINTS' in line:
            idx = i+2
            idx_minus_1 = i+1
    os.system("""sed '{}d' scf_dir/temp-{}.in | sed '{} a\ {} {} {} 0 0 0' > scf_dir/temp.in""".format(idx,mpid,idx_minus_1,kpoint[0],kpoint[1],kpoint[2]))


def main_qe(file_name,parameter,start,end):
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

    """

    print("---------------------------------------------------------------------------------------------------------------")
    print(f"Submitting convergence test calculations for {parameter}")
    print("---------------------------------------------------------------------------------------------------------------")

    for ii in range(start, end):
        with open(file_name, 'r') as file:
            for line in file:
                if f"v{ii} " in line:
                    A = line.split()[1]
                    B = line.split()[2]
                    print(f"{A} {B}")

                    if not os.path.isdir(f"R{A}-{B}"):
                        os.mkdir(f"R{A}-{B}")

                    if not os.path.isdir(f"R{A}-{B}/{parameter}"):
                        os.mkdir(f"R{A}-{B}/{parameter}")

                    # Execute sed commands
                    os.system(f"sed '/ecutrho/d' scf_dir/scf-{A}.in | sed '/calculation/d' > scf_dir/scf-{A}-1.in")
                    os.system(f"sed '/ecutwfc/d' scf_dir/scf-{A}-1.in > scf_dir/scf-{A}-2.in")
                    os.system(f"rm scf_dir/scf-{A}-1.in")
                    os.system(f"sed '/pseudo_dir/d' scf_dir/scf-{A}-2.in > scf_dir/scf-{A}-1.in")
                    os.system(f"rm scf_dir/scf-{A}-2.in")
                    os.system(f"""sed \"/&CONTROL/a calculation = 'scf'," scf_dir/scf-{A}-1.in > scf_dir/scf-{A}-2.in""")
                    os.system(f"rm scf_dir/scf-{A}-1.in")
                    os.system(f"""sed \"/calculation = 'scf',/a pseudo_dir = '../../../pp/'," scf_dir/scf-{A}-2.in > scf_dir/temp-{A}.in""")
                    os.system(f"rm scf_dir/scf-{A}-2.in")
                    if parameter == 'ecut':
                        ecut_list = input_data['conv_test']['ecut']
                        kpoint = input_data['conv_test']['kpoint'][0]
                        kpoint_qe(kpoint,A)
                        os.system(f"rm scf_dir/temp-{A}.in")
                        for jj in ecut_list:
                            if not os.path.isdir(f"R{A}-{B}/{parameter}/R{jj}"):
                                os.mkdir(f"R{A}-{B}/{parameter}/R{jj}")
                            kk = jj * 8
                            # Write modified scf input file
                            os.system(f"sed \"/&SYSTEM/a   ecutwfc = {jj}, ecutrho = {kk},\" scf_dir/temp.in > R{A}-{B}/{parameter}/R{jj}/scf.in")
                            os.system(f"cp run-scf.sh R{A}-{B}/{parameter}/R{jj}/")
                            os.chdir(f"R{A}-{B}/{parameter}/R{jj}")
                            submission('ecut',"run-scf.sh",A,B,jj)
                            os.chdir("../../../")
                            # Copy shell script and submit job
                    elif parameter == 'kpoint':
                        ecut = input_data['conv_test']['ecut'][0]
                        ecutrho = ecut * 8
                        kpoint_list = input_data['conv_test']['kpoint']
                        os.system(f"sed \"/&SYSTEM/a   ecutwfc = {ecut}, ecutrho = {ecutrho},\" scf_dir/temp-{A}.in > scf_dir/temp-{A}-1.in")
                        os.system(f"mv scf_dir/temp-{A}-1.in scf_dir/temp-{A}.in")
                        for jj in kpoint_list:
                            if not os.path.isdir(f"R{A}-{B}/{parameter}/R{jj[0]}-{jj[1]}-{jj[2]}"):
                                os.mkdir(f"R{A}-{B}/{parameter}/R{jj[0]}-{jj[1]}-{jj[2]}")
                            kpoint_qe(jj,A)
                            os.system(f"cp run-scf.sh R{A}-{B}/{parameter}/R{jj[0]}-{jj[1]}-{jj[2]}/")
                            os.system(f"mv scf_dir/temp.in R{A}-{B}/{parameter}/R{jj[0]}-{jj[1]}-{jj[2]}/scf.in")
                            os.chdir(f"R{A}-{B}/{parameter}/R{jj[0]}-{jj[1]}-{jj[2]}")
                            submission('kpoint',"run-scf.sh",A,B,jj)
                            os.chdir("../../../")
                            # Copy shell script and submit job
                    else:
                        print("Only ecut and kpoint convergence tests can be performed\n")

    print("All done")

if __name__ == "__main__":
    #if len(sys.argv) != 4:
    #    print("Usage: script.py <start_index> <end_index> <filename>")
    #    sys.exit(1)
    #insertk([10,10,10],"mp-763")
    with open("input.in","r") as read_in:
        lines = read_in.readlines()
    START = int(lines[0].split("\n")[0])
    END = int(lines[1].split("\n")[0])
    FILENAME = lines[3].split("\n")[0]
    with open(FILENAME,'r') as read_mpid:
        lines = read_mpid.readlines()
    lines = lines[START-1:END-1]
    CALC = input_data['download']['inp']['calc']
    parameter_variable = input_data['conv_test']['param']  # Modify this according to your needs
    MODE = sys.argv[1]
    if MODE == "calculate":
        if CALC in ('QE','qe'):
            main_qe(FILENAME,parameter_variable,START,END)
            os.system("rm scf_dir/temp*")
        else:
            main_vasp(FILENAME,parameter_variable,START,END)
    else:
        print("Extracting results\n")
        if not os.path.isdir("convergence_result"):
            os.mkdir("convergence_result")
        for line in lines:
            MPID = line.split(" ")[1]
            COMP = line.split(" ")[2].split("\n")[0]
            FOLDER = f"R{MPID}-{COMP}/{parameter_variable}"
            os.chdir(FOLDER)
            if CALC in ('VASP','vasp'):
                with open("param-energy.txt", "w") as write_energy:
                    list_param = input_data['conv_test'][parameter_variable]
                    for value in list_param:
                        if parameter_variable == 'ecut':
                            value1 = value
                        else:
                            value1 = f"{value[0]}-{value[1]}-{value[2]}"
                        os.chdir(f"R{value1}")
                        data = Vasprun("vasprun.xml")
                        energy = data.final_energy.real
                        write_energy.write(f"{value1} {energy}" + "\n")
                        os.chdir("../")
                    os.chdir("../../")
                os.system(f"cp {FOLDER}/param-energy.txt convergence_result/{parameter_variable}-{MPID}-{COMP}.txt")
            else:
                with open("param-energy.txt", "w") as write_energy:
                    list_param = input_data['conv_test'][parameter_variable]
                    for value in list_param:
                        if parameter_variable == 'ecut':
                            value1 = value
                        else:
                            value1 = f"{value[0]}-{value[1]}-{value[2]}"
                        os.chdir(f"R{value1}")
                        data = PWOutput("scf.out")
                        energy = data.final_energy.real
                        write_energy.write(f"{value1} {energy}" + "\n")
                        os.chdir("../")
                    os.chdir("../../")
                os.system(f"cp {FOLDER}/param-energy.txt convergence_result/{parameter_variable}-{MPID}-{COMP}.txt")
