#!/usr/bin/env python
"""
USAGE: python bmfit.py energy-volume ncol-1

ncol is number of column in 'energy-volume' file

"""
'''Example of fitting the Birch-Murnaghan EOS to data'''

#from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
from pylab import * #this includes numpy as np!
from scipy.optimize import leastsq
import sys

# create some Equation of State (EOS) functions 
def Murnaghan(parameters,vol):
    '''
    given a vector of parameters and volumes, return a vector of energies.
    equation From PRB 28,5480 (1983)
    '''
    E0 = parameters[0]
    B0 = parameters[1]
    BP = parameters[2]
    V0 = parameters[3]
    
    E = E0 + B0*vol/BP*(((V0/vol)**BP)/(BP-1)+1) - V0*B0/(BP-1.)

    return E

def Birch_Murnaghan(parameters,vol):
    '''
    given a vector of parameters and volumes, return a vector of energies.
    '''
    E0  = parameters[0]
    B0  = parameters[1]
    BP  = parameters[2]
    V0  = parameters[3]
    eta = [V0 / x for x in vol]
    E   = []

    for ii in eta:
        E.append(E0 + B0*V0*(9./16.)*( BP*(ii**(2./3.) - 1.)**3. + (6.-4.*ii**(2./3.))*(ii**(2./3.)-1.)**2. ))

    return E
  
def SJEOSE(parameters, vol):
    E0  = parameters[0]
    B0  = parameters[1]
    BP  = parameters[2]
    V0  = parameters[3]
    eta = [V0 / x for x in vol]
    E   = []
    for ii in eta:
        E.append( -E0 + 9./2.*B0*V0*(4. - BP + (-3. + BP)*ii
         - (11. - 3.*BP)*(ii)**(1./3.) + (10. - 3.*BP)*(ii)**(2./3.)))
    return E

nmethod = int(sys.argv[2])
evec = []
at = np.loadtxt(str(sys.argv[1]))
v = np.array([x[0] for x in at])
for ii in range(1,nmethod+1):
  evec.append(np.array([x[ii] for x in at]))

resultsv = []
for e in evec:
  #make a vector of volumes to evaluate fits on with a lot of points so it looks smooth
  vfit = np.linspace(min(v),max(v),100)
  
  ### fit a parabola to the data
  # y = ax^2 + bx + c
  a,b,c = polyfit(v,e,2) #this is from pylab
  '''
  the parabola does not fit the data very well, but we can use it to get
  some analytical guesses for other parameters.
  
  V0 = minimum energy volume, or where dE/dV=0
  E = aV^2 + bV + c
  dE/dV = 2aV + b = 0
  V0 = -b/2a
  
  E0 is the minimum energy, which is:
  E0 = aV0^2 + bV0 + c
  
  B is equal to V0*d^2E/dV^2, which is just 2a*V0
  
  and from experience we know Bprime_0 is usually a small number like 4
  '''
  
  #now here are our initial guesses.
  v0 = -b/(2*a)
  e0 = a*v0**2 + b*v0 + c
  b0 = 2*a*v0
  bP = 4
  
  # and we define an objective function that will be minimized
  def objective(pars,y,x):
      #we will minimize this function
      err =  y - Birch_Murnaghan(pars,x)
      return err
  def objectivem(pars,y,x):
      #we will minimize this function
      err =  y - SJEOSE(pars,x)
      return err
  
  x0 = [e0, b0, bP, v0] #initial guesses in the same order used in the Murnaghan function
  
  # use leastsq fit to determine the fit parameters
  murnpars, ier1 = leastsq(objective, x0, args=(e,v)) #this is from scipy
  murnpars2, ier2 = leastsq(objectivem, x0, args=(e,v)) #this is from scipy
  
  #now we make a figure summarizing the results
  plot(v,e,'ro')
  plot(vfit, a*vfit**2 + b*vfit + c,'--',label='parabolic fit')
  plot(vfit, SJEOSE(murnpars2,vfit), '--', c='r', label='SJEOS fit')
  plot(vfit, Birch_Murnaghan(murnpars,vfit), c='k', label='Birch-Murnaghan fit')
  xlabel('Volume ($\AA^3$)')
  ylabel('Energy (eV)')
  legend(loc='best')
  
  #add some text to the figure in figure coordinates
  ax = gca()
  # BM params
  text(0.4,0.45,'BM Min volume = %1.8f $\AA^3$' % murnpars[3],
       transform = ax.transAxes)
  text(0.4,0.4,'BM Bulk modulus = %1.8f GPa' % (murnpars[1]*160.21773)
       , transform = ax.transAxes)
  # SJEOS params
  text(0.4,0.55,'M Min volume = %1.8f $\AA^3$' % murnpars2[3],
       transform = ax.transAxes)
  text(0.4,0.5,'M Bulk modulus = %1.8f GPa' % (murnpars2[1]*160.21773)
       , transform = ax.transAxes)
  savefig(sys.argv[1]+'.png')
  show()
  resultsv.append(murnpars)
  
  
  #print('BM lattice const : ', (murnpars[3]**(0.333333333333333)*2.0))
  #print('BM Bulk mod      : ', (murnpars[1]*160.21773))
resultsv = np.array(resultsv)
tags = ['E', 'B', 'dB/dP', 'V0']
for jj,row in enumerate(resultsv.T):
    print('  %6s  ' % tags[jj], end=" ")
    for col in row:
      if jj == 1:
        print('  %9.4f  ' % (col*160.21773), end=" ")
      else:
        print('  %9.4f  ' % col, end=" ")
    print('')
