import sys
import pylab
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')
font = {'weight' : 'bold','size'   : 25}
matplotlib.rc('font', **font)
plt.rcParams["figure.autolayout"] = True



color = ['g', 'b', 'g']
marker = ['x', '+', 'o']
linestyle = ['--', 'dotted', '-.']
extreme_a2f = []
extreme_en = []
THZTOMEV = 4.135665538536
fig = plt.gcf()
fig.set_size_inches(8, 4)
def a2f_compare(filename,nfiles):
    """
    Function to plot isotropic Eliashberg spectral function
    parameters
    ------------
    filename : name of the output
    nfiles : Number of alpha2F files in alpha2F-i.dat format, i = 1, 2,...
    """
    for i in range(nfiles):
        data = np.loadtxt('alpha2F-{}.dat'.format(i+1))
        endata = data[:,0]
        a2f = data[:,2]
        delw = endata[1]
        endata[np.where(endata == 0.0)] = 0.0000000001
        a2fint = 2.0*delw*a2f/endata
        a2fint = a2fint.cumsum()
        b2f = a2f
        extreme_a2f.append(b2f.max())
        extreme_a2f.append(b2f.min())
        extreme_en.append(endata.min()*THZTOMEV)
        extreme_en.append(endata.max()*THZTOMEV)
        lmbd = round(a2fint[-1],4)
        plt.plot(endata,b2f,color[i],linestyle='-',lw=1.5)
        plt.plot(endata,a2fint,color[i],linestyle=linestyle[i],label=str(lmbd), lw=1.5)
    maxy = np.max(extreme_a2f) + 0.2
    maxx = 25.0
    miny = np.min(extreme_a2f) - 0.02
    plt.ylim((miny,maxy))
    plt.legend(fontsize=20)
    plt.xlim((0.0,maxx))
    plt.xlabel(r'$\omega$' + '(THz)', fontsize=20)
    plt.xticks([0,5,10,15,20,25],['0','5','10','15','20','25'],fontsize=20)
    plt.yticks([0,0.5,1.0,1.5,2.0],['0','0.5','1.0','1.5','2.0'],fontsize=20)
    plt.ylabel(r'$\alpha^2F$', fontsize=30)
    plt.savefig(filename+"-a2f.png")
    pylab.savefig(filename+ "-a2f.pdf", format='pdf')
if __name__ == "__main__":
    NFILES = int(sys.argv[1])
    a2f_compare('a2F-combine',NFILES)
