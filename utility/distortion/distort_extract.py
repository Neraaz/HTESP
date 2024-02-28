"""Writen by Niraj K. Nepal, Ph.D. (tug11655@temple.edu)"""
import sys
import pandas as pd
import numpy as np
MPID = sys.argv[1]
COMP = sys.argv[2]
data = pd.read_csv('R{}-{}/Energy-mode.csv'.format(MPID,COMP))
data = data['Energy(eV)']
ENERGY = np.unique(data.round(decimals=4))
with open("distorted-energy.csv", "a") as distort:
    distort.write(MPID + "," + COMP + ",")
    for i in range(ENERGY.shape[0]):
        distort.write(str(ENERGY[i]) + ",")
    distort.write("\n")
