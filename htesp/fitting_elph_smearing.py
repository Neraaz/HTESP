#!/usr/bin/env python
#Writen by Niraj K. Nepal, Ph.D.
"""
This code extract data from lambda.out file and plot Tc vs sigma,
fit Tc = exp(-A*sigma**(1/3) + B) to obtain the value of decay
parameter "A".
"""
import sys
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def func(xdata,aval,bval):
    """
    Computes the value of a linear function at given points.

    Parameters:
    -----------------
    xdata : array-like
        The input data points where the function is evaluated.

    aval : float
        The coefficient 'a' in the linear function.

    bval : float
        The intercept 'b' in the linear function.

    Returns:
    -----------------
    ydata : array-like
        The computed values of the linear function at the input data points.
    """
    return aval*xdata + bval

def parse_lambda(filename):
    """
    Parses the lambda.out file obtained from lambda.x of Quantum ESPRESSO (QE).
    Creates a new file named "lambda.txt" containing columns for lambda, omega_log, and Tc.

    Parameters:
    -----------------
    filename : str
        The name of the file to parse, typically "lambda.out".

    Note:
    -----------------
    It extracts the rows that follow the 'omega_log' header line.

    Returns
    -------
    numpy.ndarray
        Array of shape (nsigma, 3) holding lambda, omega_log and Tc.
    """
    # FIX: `os.system("sed -n '/omega_log/,$p' ... | sed 1d > lambda.txt")`
    # replaced by a pure-Python parse - no shell, no scratch file to delete.
    with open(filename, "r") as read_lambda:
        lines = read_lambda.readlines()
    rows = []
    started = False
    for line in lines:
        if not started:
            if "omega_log" in line:
                started = True
            continue
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            rows.append([float(value) for value in fields[:3]])
        except ValueError:
            # a '*' or 'NaN' entry - caller's sanity check already covers it
            continue
    if not rows:
        raise ValueError("no lambda/omega_log/Tc rows found in {}".format(filename))
    return np.array(rows)

def fit_param(mpid,compound,ind,a_mgb2,begin_sigma,n_sigma):
    """
    Performs a fitting procedure to obtain exponential decay parameters and creates plots of the fitting.

    Parameters:
    -----------------
    mpid : str
        Materials ID.
    compound : str
        Compound name.
    ind : int
        Index of property. 1 for lambda and 2 for Tc.
    a_mgb2 : float
        Decay parameter obtained by fitting Tc vs sigma data of MgB2.
    begin_sigma : float
        Width of smearing used in double delta integration.
    n_sigma : int
        Number of sigma values.
    """
    filename="R{}-{}/calc/lambda.out".format(mpid,compound)
    if os.path.isfile(filename):
        # FIX(19): both greps wrote to the *same* file "file", so the second
        # overwrote the first and the '*' (overflow) check never ran.  Both
        # counts are now taken in Python and combined.
        with open(filename, "r") as read_lambda:
            lambda_text = read_lambda.read()
        n_star = sum(1 for line in lambda_text.splitlines() if "*" in line)
        n_nan = sum(1 for line in lambda_text.splitlines() if "NaN" in line)
        check = float(n_star + n_nan)
        if check < 1.0:
            try:
                prop = ["lambda","wlog", "Tc"]
                out = prop[ind-1]+"-{}-{}.png".format(mpid,compound)
                data1 = parse_lambda(filename)
                data = data1[:,ind-1]
                data[np.where(data == 0.0)] = 0.00001
                logy = np.log(data)
                factor = 1.0/3.0
                smlist = np.linspace(begin_sigma,begin_sigma*n_sigma,n_sigma)
                fit = np.polyfit(smlist**factor,logy,1)
                #print("{}: [A,B] = {}".format(prop[ind-1],fit))
                if data[3] < 0.1 and data[9] < 0.1:
                    with open("fitting_params_{}_zero.csv".format(prop[ind-1]), "a") as fwrite:
                        fwrite.write(mpid + "," + compound + ",")
                        fwrite.write(str(fit[0]) + "," + str(fit[1]) +",")
                        fwrite.write(str(data[0])+","+str(data[1])+",")
                        fwrite.write(str(data[2])+","+str(data[3])+",")
                        fwrite.write(str(data[4])+","+str(data[5])+",")
                        fwrite.write(str(data[6])+","+str(data[7])+",")
                        fwrite.write(str(data[8])+","+str(data[9]) + "\n")
                else:
                    if fit[0] < a_mgb2:
                        with open("fitting_params_{}_nonzero_noconv.csv".format(prop[ind-1]), "a") as fwrite:
                            fwrite.write(mpid + "," + compound + ",")
                            fwrite.write(str(fit[0]) + "," + str(fit[1]) +",")
                            fwrite.write(str(data[0])+","+str(data[1])+",")
                            fwrite.write(str(data[2])+","+str(data[3])+",")
                            fwrite.write(str(data[4])+","+str(data[5])+",")
                            fwrite.write(str(data[6])+","+str(data[7])+",")
                            fwrite.write(str(data[8])+","+str(data[9]) + "\n")
                    else:
                        with open("fitting_params_{}_nonzero_conv.csv".format(prop[ind-1]), "a") as fwrite:
                            fwrite.write(mpid + "," + compound + ",")
                            fwrite.write(str(fit[0]) + "," + str(fit[1]) +",")
                            fwrite.write(str(data[0])+","+str(data[1])+",")
                            fwrite.write(str(data[2])+","+str(data[3])+",")
                            fwrite.write(str(data[4])+","+str(data[5])+",")
                            fwrite.write(str(data[6])+","+str(data[7])+",")
                            fwrite.write(str(data[8])+","+str(data[9]) + "\n")
                y_fit = np.exp(fit[0]*smlist**factor + fit[1])
                plt.figure()
                plt.plot(smlist, data, 'r', lw=2.5)
                plt.plot(smlist, y_fit, 'k--o')
                plt.ylim(data.min()-0.1, data.max() + 0.1)
                os.makedirs("plots_fit", exist_ok=True)
                plt.savefig("plots_fit/" + out)
                plt.cla()
                plt.close()
            # FIX: a bare `except:` hid every failure, KeyboardInterrupt
            # included, and never said what went wrong.
            except (OSError, ValueError, IndexError, TypeError,
                    np.linalg.LinAlgError) as err:
                with open("problem_in_fit.in", "a") as fwrite:
                    fwrite.write("Problem in lambda or logomega, could be "
                                 "negative: {}-{} ({}: {})".format(
                                     mpid, compound, type(err).__name__, err) + "\n")
        else:
            with open("problem_in_fit.in", "a") as fwrite:
                fwrite.write("*** or NaN present on {}-{}".format(mpid,compound) + "\n")
        return
def main():
    """
    Main function to perform fitting procedure for electron-phonon coupling and handle input/output files.

    Usage: python fitting_elph_smearing.py result.csv AMgB2 smearing_begin number_of_smearing

    Parameters:
    -----------------
    result.csv : str
        CSV file containing results of electron-phonon coupling calculations.
    AMgB2 : float
        An exponential decay parameter for MgB2.
    smearing_begin : float
        Starting value of the smearing parameter.
    number_of_smearing : int
        Number of smearing parameter values.

    Returns:
    -----------------
    None
    """
    print("Usage: python fitting_elph_smearing.py result.csv")
    print(" AMgB2 smearing_begin number_of_smearing\n")
    id_file = pd.read_csv(sys.argv[1])
    a_mgb2 = -1.0*float(sys.argv[2])
    begin_sigma = float(sys.argv[3])
    n_sigma = int(sys.argv[4])
    id_file = id_file.dropna().reset_index(drop=True)
    id_file[id_file.Phonon_freq == "negative_freq"].to_csv("result_negative.csv",index=False)
    id_file = id_file[id_file.Phonon_freq == "nonegative_freq"].reset_index(drop=True)
    id_file.to_csv("result_non_negative.csv",index=False)
    idx = id_file["ID"]
    comp = id_file["compound"]
    n = id_file.shape[0]
    os.makedirs("Problematic", exist_ok=True)
    with open("fitting_params_Tc_nonzero_conv.csv", "w") as fwrite_file:
        fwrite_file.write("ID,compound,A,B,data1,data2,data3")
        fwrite_file.write(",data4,data5,data6,data7,data8,data9,data10\n")
    with open("fitting_params_Tc_nonzero_noconv.csv", "w") as fwrite_file:
        fwrite_file.write("ID,compound,A,B,data1,data2,data3")
        fwrite_file.write(",data4,data5,data6,data7,data8,data9,data10\n")
    with open("fitting_params_Tc_zero.csv", "w") as fwrite_file:
        fwrite_file.write("ID,compound,A,B,data1,data2,data3")
        fwrite_file.write(",data4,data5,data6,data7,data8,data9,data10\n")
    for i in range(n):
        #print(idx[i],comp[i])
        fit_param(idx[i],comp[i],3,a_mgb2,begin_sigma,n_sigma)
if __name__ == "__main__":
    main()
