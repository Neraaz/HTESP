"""
Equation of state, Niraj K. Nepal, Ph.D.
"""
import sys
import warnings
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
from pylab import polyfit
matplotlib.use('Agg')
warnings.filterwarnings("ignore")
def birch_murnaghan(parameters,vol):
    """
    Compute the energy for the Birch-Murnaghan equation of state.

    Parameters:
    -----------
    parameters : list
        A list of parameters [E0, B0, BP, V0] where:
        - E0: Energy at equilibrium volume
        - B0: Bulk modulus at equilibrium volume
        - BP: Pressure derivative of bulk modulus
        - V0: Equilibrium volume
    vol : list
        List of volumes at which to compute the energy.

    Returns:
    --------
    list
        A list of energies corresponding to the given volumes.
    """
    E0  = parameters[0]
    B0  = parameters[1]
    BP  = parameters[2]
    V0  = parameters[3]
    eta = [V0 / x for x in vol]
    E   = []
    for ii in eta:
        E.append(E0 + B0*V0*(9./16.)*( BP*(ii**(2./3.) - 1.)**3. + (6.-4.*ii**(2./3.))*(ii**(2./3.)-1.)**2. ))
    return E
def objective(pars,y,x):
    """
    Calculate the objective function for optimization.

    Parameters:
    -----------
    pars : list
        A list of parameters [E0, B0, BP, V0] used in the Birch-Murnaghan equation of state.
    y : list
        List of energy values.
    x : list
        List of volumes corresponding to the energy values.

    Returns:
    --------
    list
        A list of errors between the computed energies using the Birch-Murnaghan equation
        and the actual energy values.
    """
    err =  y - birch_murnaghan(pars,x)
    return err
def fit(v,e):
    """
    Perform fitting to determine the parameters of the Birch-Murnaghan equation of state.

    Parameters:
    -----------
    v : list
        List of volumes.
    e : list
        List of energy values corresponding to the volumes.

    Returns:
    --------
    tuple
        A tuple containing the parameters of the Birch-Murnaghan equation of state
        and an integer indicating the success of the optimization process.
    """
    a,b,c = polyfit(v,e,2)
    v0 = -b/(2*a)
    e0 = a*v0**2 + b*v0 + c
    b0 = 2*a*v0
    bP = 4
    x0 = [e0, b0, bP, v0] #initial guesses in the same order used in the Murnaghan function
  # use leastsq fit to determine the fit parameters
    murnpars, ier1 = leastsq(objective, x0, args=(e,v)) #this is from scipy
    return murnpars, ier1
def pressure(pars,volume):
    """
    Calculate the pressure according to the Birch-Murnaghan equation of state.

    Parameters:
    -----------
    pars : tuple
        A tuple containing the parameters of the Birch-Murnaghan equation of state.
        It includes E0, B0, BP, and V0.
    volume : float
        The volume at which pressure needs to be calculated.

    Returns:
    --------
    float
        The pressure calculated based on the Birch-Murnaghan equation of state.
    """
    _,b0,bp,v0 = pars
    p = (3.0*b0/2.0) * ((v0/volume)**(7.0/3.0) - (v0/volume)**(5.0/3.0)) * (1.0 + (3.0*(bp-4.0)/4.0)*((v0/volume)**(2.0/3.0)-1))
    return p
def solve(x,y,pars):
    """
    Solve the pressure-volume equation of state.

    This function calculates the difference between the observed pressure (y)
    and the pressure calculated using the Birch-Murnaghan equation of state
    with the given parameters (pars) and volumes (x).

    Parameters:
    -----------
    x : array_like
        An array of volumes.
    y : array_like
        An array of observed pressures.
    pars : tuple
        A tuple containing the parameters of the Birch-Murnaghan equation of state.
        It includes E0, B0, BP, and V0.

    Returns:
    --------
    array_like
        The difference between the observed pressure and the pressure calculated
        using the Birch-Murnaghan equation of state with the given parameters and volumes.
    """
    return y - pressure(pars,x)
def plot(ndata,outfile,type1,pstart,pend):
    """
    Plot function to visualize data.

    Parameters:
    -----------
    ndata : int
        Number of datasets to plot.
    outfile : str
        Output file name to save the plot.
    type1 : str
        Type of plot to generate. Should be one of ['ev', 'vp', 'hv', 'hp'].
    pstart : float
        Starting value of pressure.
    pend : float
        Ending value of pressure.

    Returns:
    --------
    None
    """
    _,ax_in = plt.subplots()
    color = ['r', 'g', 'b']
    label = ['phase1', 'phase2', 'phase3']
    #marker = ['s', '^', 'D']
    linestyle = ['solid', 'dotted', '--', '-.']
    xlabels_dict={"ev":r"V ($\AA^3$)", "vp":"P (GPa)", "hv":r"V ($\AA^3$)","hp":"P (GPa)"}
    ylabels_dict={"ev":"E (eV)", "vp":r"V ($\AA^3$)", "hv":"H (eV)","hp":"H (eV)"}
    ax_in.set_xlabel(xlabels_dict[type1], fontsize=20)
    ax_in.set_ylabel(ylabels_dict[type1], fontsize=20)
    ax_in.tick_params(axis='both', labelsize=15)
    for i in range(ndata):
        data = np.loadtxt("e-v-{}.dat".format(i+1))
        volume = data[:,0]
        energy = data[:,1]
        vfit = np.linspace(np.min(volume),np.max(volume),200)
        vol = vfit[::-1]
        press = np.linspace(pstart,pend,200)/160.217662
        pfit = np.linspace(pstart,pend,200)
        params,_ = fit(volume,energy)
        vol_final = []
        for j,press_x in enumerate(press):
            vol_sol = leastsq(solve, vol[j], args=(press_x,params))
            vol_final.append(vol_sol[0][0])
        vol_final = np.array(vol_final)
        #p = pressure(params,vfit)*160.217662
        efit = np.array(birch_murnaghan(params,vol_final))
        enthalpy = efit + pfit*vol_final*0.00624151
        #pv_term = p*vfit*0.00624151
        if i == 0:
            ebcc = enthalpy
        if type1 == "ev":
            efit_vfit = np.array(birch_murnaghan(params,vfit))
            ax_in.plot(vfit,efit_vfit, color[i], label=label[i],lw=1.5,linestyle=linestyle[i])
        elif type1 == "vp":
            ax_in.plot(pfit,vol_final, color[i], label=label[i],lw=1.5,linestyle=linestyle[i])
        elif type1 == "hv":
            vol_rev = vol_final[::-1]
            enthalpy_rev = enthalpy[::-1]
            ax_in.plot(vol_rev,enthalpy_rev, color[i], label=label[i],lw=1.5,linestyle=linestyle[i])
        elif type1 == "hp":
            ax_in.plot(pfit,enthalpy - ebcc, color[i], label=label[i],lw=1.5,linestyle=linestyle[i])
        else:
            print("Only ev, vp, hv, and hp type allowed\n")
        #ax_in.plot(vfit,efit, color[i], linestyle=linestyle[i], label=label[i], marker=marker[i],markerfacecolor='None',markersize=2,markeredgewidth=1)
    #if type1 == 'hp':
    #    ax_in.set_xlim(8,20)
    #    ax_in.set_ylim(-0.6,0.2)
    ax_in.legend(ncol=3,loc='best',fontsize='large',prop={'size': 15},frameon=False)
    plt.savefig(outfile, format='pdf', bbox_inches='tight')
if __name__ == "__main__":
    NDATA = sys.argv[1]
    if NDATA in ('help', 'H', 'h'):
        print("========================================================================================================\n")
        print("This script plots various equation of state plots from energy-volume data\n")
        print("Usage: python birch_murnaghan_enthalpy.py <number_of_energy_volume_data> <type_of_plot> <pressure_start> <pressure_end>\n")
        print("Options for <type_of_plot> are hp (enthalpy-pressure), vp (volume-pressure), ev (energy-volume), hv (enthalpy-volume)\n")
        print("Put your energy-volume data in e-v-{i}.dat format, i = 1,2,...\n")
        print("Make sure that pressure range you provide is compatible with volume range in your energy-volume data\n")
        print("Check minimum and maximum volumes, and associated external pressure in your (QE or VASP) output file\n")
        print("Adjust the limits of x-axis (xlim) and y-axis (ylim) according based on reliability of your data.")
        print("========================================================================================================\n")
    else:
        NDATA = int(NDATA)
        TYPE = sys.argv[2]
        PSTART = float(sys.argv[3])
        PEND = float(sys.argv[4])
        OUTFILE=TYPE + ".pdf"
        plot(NDATA,OUTFILE,TYPE,PSTART,PEND)
