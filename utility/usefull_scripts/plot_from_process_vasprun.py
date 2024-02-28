# coding: utf-8
from matplotlib import pyplot as plt
import numpy as np
data = np.loadtxt("atom_dos_orb.txt")
data = data.reshape(5,2000,4)
dftot = np.zeros(2000)
for i in range(5):
    df = data[i]
    df1 = df[:,1]
    df2 = df[:,2]
    df3 = df[:,3]
    dftot += df1 + df2 + df3
    en = df[:,0]
    plt.plot(en,df1)
    plt.plot(en,df2)
    plt.plot(en,df3)
plt.plot(data[0][:,0],dftot)
    
plt.xlim(2,10)
plt.ylim(0,90)
plt.savefig("dos1.pdf")
