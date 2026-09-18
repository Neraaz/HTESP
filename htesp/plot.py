#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module for plotting"""
import sys
import os
import glob
import re
import subprocess
import warnings
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pylab
import numpy as np
from pymatgen.io.vasp.inputs import Kpoints, Incar
from pymatgen.io.vasp.sets import Vasprun
from pymatgen.core import Composition
from pymatgen.electronic_structure.core import OrbitalType
from pymatgen.electronic_structure.plotter import DosPlotter
from pymatgen.electronic_structure.plotter import BSPlotter
from htesp.kpath import kpath
from htesp.check_json import config
warnings.filterwarnings("ignore")

font = {'weight' : 'bold','size'   : 20}
matplotlib.rc('font', **font)
plt.rcParams["figure.autolayout"] = True

_FLOAT = r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?"
# FIX(14): QE prints the Fermi level once per ionic step of a relaxation, and
# "Fermi-Dirac" also contains the word "Fermi".  Match the real statement only
# and take the LAST occurrence.
_QE_FERMI = re.compile(r"the Fermi energy is\s+(" + _FLOAT + r")")
_QE_FERMI_SPIN = re.compile(
    r"the spin up/dw Fermi energies are\s+(" + _FLOAT + r")\s+(" + _FLOAT + r")")
# FIX(14): VASP-5 writes "E-fermi :  x", VASP-6 writes "Fermi energy:  x".
# Both wordings are handled in ONE pass, and OUTCAR holds one such line per
# ionic step, so the LAST one is the converged value.
_VASP_FERMI = re.compile(
    r"(?:E-fermi\s*:|Fermi energy\s*:)\s*(" + _FLOAT + r")")
# FIX(18): the .dos header is "# E (eV)  dos(E) Int dos(E) EFermi = x eV" for
# nspin=1 but gains a column for nspin=2, so token [8] is the wrong field.
_DOS_EFERMI = re.compile(r"EFermi\s*=\s*(" + _FLOAT + r")")


def _to_float(text):
    """Convert a Fortran-style float ('1.0D-3') to a Python float."""
    return float(str(text).replace("D", "E").replace("d", "e"))


def fermi_from_qe_out(path):
    """
    Read the Fermi level from a QE ``scf.out``.

    FIX(14): the old code was ``grep Fermi scf.out | awk '{print $5}'`` and
    then took the FIRST line - the first ionic step of a relaxation - and it
    also matched "Fermi-Dirac".  The LAST genuine statement is used now, and a
    spin-polarised run (which prints two energies) is handled explicitly.

    Parameters
    ----------
    path : str
        Path to the QE output.

    Returns
    -------
    float
        The Fermi energy in eV.
    """
    with open(path, "r", errors="replace") as read_out:
        text = read_out.read()
    single = _QE_FERMI.findall(text)
    both = _QE_FERMI_SPIN.findall(text)
    if single and (not both or text.rfind("the Fermi energy is") >
                   text.rfind("the spin up/dw Fermi energies are")):
        return _to_float(single[-1])
    if both:
        up, down = (_to_float(value) for value in both[-1])
        print("spin-polarised run: Fermi energies up = {}, dw = {}; "
              "using the higher one\n".format(up, down))
        return max(up, down)
    raise ValueError("no Fermi energy found in {}".format(path))


def fermi_from_outcar(path):
    """
    Read the Fermi level from a VASP ``OUTCAR`` (VASP 5 or 6).

    FIX(14): the old code grepped only one of the two wordings, then called
    ``.item()`` on the loaded array - which raises for any relaxation with more
    than one ionic step, because OUTCAR then holds one such line per step.
    Both wordings are matched here and the LAST value is returned.

    Parameters
    ----------
    path : str
        Path to the OUTCAR.

    Returns
    -------
    float
        The Fermi energy in eV.
    """
    with open(path, "r", errors="replace") as read_out:
        text = read_out.read()
    matches = _VASP_FERMI.findall(text)
    if not matches:
        raise ValueError(
            "no 'E-fermi :' (VASP 5) or 'Fermi energy:' (VASP 6) line "
            "in {}".format(path))
    return _to_float(matches[-1])


def read_fermi(qe_out="scf.out", outcar="../relax/OUTCAR"):
    """
    Return the Fermi level from whichever of the two outputs is present.

    Parameters
    ----------
    qe_out : str
        QE output to try first. Default 'scf.out'.
    outcar : str
        VASP output to try second. Default '../relax/OUTCAR'.

    Returns
    -------
    float
        The Fermi energy in eV.
    """
    if os.path.isfile(qe_out):
        return fermi_from_qe_out(qe_out)
    if os.path.isfile(outcar):
        return fermi_from_outcar(outcar)
    raise FileNotFoundError(
        "neither {} nor {} is present".format(qe_out, outcar))


def nelect_from_outcar(path):
    """Read NELECT from a VASP OUTCAR (last occurrence)."""
    with open(path, "r", errors="replace") as read_out:
        values = re.findall(r"NELECT\s*=\s*(" + _FLOAT + r")", read_out.read())
    if not values:
        raise ValueError("no NELECT line in {}".format(path))
    return _to_float(values[-1])


def nelect_from_qe_out(path):
    """Read 'number of electrons' from a QE output (last occurrence)."""
    with open(path, "r", errors="replace") as read_out:
        values = re.findall(
            r"number of electrons\s*=\s*(" + _FLOAT + r")", read_out.read())
    if not values:
        raise ValueError("no 'number of electrons' line in {}".format(path))
    return _to_float(values[-1])


def read_lsorbit(incar="INCAR"):
    """
    Report whether spin-orbit coupling is switched on.

    FIX(15): the old test was ``grep LSORBIT INCAR | wc -l``, so an INCAR
    holding ``LSORBIT = .FALSE.`` counted as spin-orbit ON and the band index
    at the Fermi level was then left un-halved.  The value is parsed now.

    Parameters
    ----------
    incar : str
        Path to the INCAR. Default 'INCAR'.

    Returns
    -------
    bool
        True when LSORBIT is true.
    """
    if not os.path.isfile(incar):
        return False
    try:
        return bool(Incar.from_file(incar).get("LSORBIT", False))
    except (OSError, ValueError, KeyError) as err:
        print("could not parse {} ({}); assuming LSORBIT = .FALSE.\n".format(
            incar, err))
        return False


def parse_lambda_out(filename="lambda.out"):
    """
    Parse the ``lambda / omega_log / T_c`` table written by QE's lambda.x.

    FIX(17): lambda used to be read out of the hard-coded line 12 of
    lambda.out, i.e. one arbitrary smearing out of the ten el-ph smearings.
    The whole table is parsed here so the caller can choose a smearing.

    Parameters
    ----------
    filename : str
        Path to lambda.out. Default 'lambda.out'.

    Returns
    -------
    tuple
        ``(rows, header_index)`` - `rows` is an (nsigma, 3) array of
        lambda/omega_log/Tc and `header_index` is the 0-based line number of
        the 'omega_log' header, which is what the legacy default is measured
        against.
    """
    with open(filename, "r", errors="replace") as read_lambda:
        lines = read_lambda.readlines()
    header_index = None
    rows = []
    for i, line in enumerate(lines):
        if header_index is None:
            if "omega_log" in line:
                header_index = i
            continue
        fields = line.split()
        if len(fields) < 3:
            continue
        try:
            rows.append([float(value) for value in fields[:3]])
        except ValueError:
            continue
    if header_index is None or not rows:
        raise ValueError(
            "no lambda/omega_log/Tc table found in {}".format(filename))
    return np.array(rows), header_index


def select_smearing(rows, header_index, smearing_index=None,
                    legacy_line=12):
    """
    Choose which el-ph smearing to report, defaulting to today's behaviour.

    Parameters
    ----------
    rows : numpy.ndarray
        Table returned by :func:`parse_lambda_out`.
    header_index : int
        Header line index returned by :func:`parse_lambda_out`.
    smearing_index : int or None
        1-based smearing number.  None reproduces the historical choice
        (``lambda.out`` line `legacy_line`) so existing plots do not change.
    legacy_line : int
        The line the old code hard-coded. Default 12.

    Returns
    -------
    tuple
        ``(index, row)`` where `index` is 1-based.
    """
    if smearing_index is None:
        index = legacy_line - (header_index + 1)
        index = min(max(index, 0), len(rows) - 1)
    else:
        index = int(smearing_index) - 1
        if not 0 <= index < len(rows):
            raise IndexError(
                "smearing {} is out of range; lambda.out holds {} "
                "smearings".format(smearing_index, len(rows)))
    return index + 1, rows[index]


def kpt_labels(kpoints_files=("KPOINTS", "KPOINTS_OPT")):
    """
    Return the high-symmetry labels of a line-mode VASP k-point file.

    FIX(16): ``Kpoints.from_file("KPOINTS").labels`` is None whenever the run
    used ``kpt_opt: true`` (the shipped default), because the line-mode path
    then lives in KPOINTS_OPT and KPOINTS holds a plain mesh.  The old code
    fed that None straight into ``len()`` and every VASP band plot died with a
    TypeError.  KPOINTS_OPT, and then the ``label``/``lbpoint`` files the
    wannier band step writes, are tried in turn.

    Parameters
    ----------
    kpoints_files : sequence of str
        Candidate k-point files, in order of preference.

    Returns
    -------
    tuple
        ``(labels, num_kpts)``; `num_kpts` is None when the labels came from
        the `label` file rather than from a KPOINTS file.
    """
    for name in kpoints_files:
        if not os.path.isfile(name):
            continue
        try:
            kpt = Kpoints.from_file(name)
        except (OSError, ValueError, IndexError) as err:
            print("could not read {} ({})\n".format(name, err))
            continue
        if kpt.labels:
            return list(kpt.labels), int(kpt.num_kpts)
        print("{} carries no high-symmetry labels "
              "(not a line-mode file)\n".format(name))
    if os.path.isfile("label"):
        with open("label", "r") as read_label:
            labels = [line.split()[0] for line in read_label if line.split()]
        if labels:
            print("using the high-symmetry labels from the 'label' file\n")
            return labels, None
    raise RuntimeError(
        "no high-symmetry labels available: none of {} is a line-mode "
        "k-point file and no 'label' file was written.  Re-run the band step, "
        "or set kpt_opt correctly in config.json.".format(
            ", ".join(kpoints_files)))


def kptline():
    """
    Generate k-path for VASP-line format.

    Reads KPOINTS file and vasprun.xml file to generate the k-path for the VASP-line format.
    Returns:
    - nkpt (numpy.ndarray): Array containing distances between k-points.
    - sympoint (numpy.ndarray): Array containing the special k-points along the path.
    - symname (list): List of special k-point names formatted for LaTeX.

    This function computes the k-path for the VASP-line format based on the KPOINTS file and vasprun.xml file.
    It calculates the distances between k-points and identifies special k-points along the path.

    Example usage:
    nkpt, sympoint, symname = kptline()
    """
    # FIX(16): labels come from the first line-mode file that actually has
    # them (KPOINTS, then KPOINTS_OPT, then the 'label' file).
    symname, nline = kpt_labels()
    data = Vasprun("vasprun.xml")
    nkpt = np.array(data.get_band_structure("KPOINTS").distance)
    n = []
    npair = int(len(symname)/2)
    if nline is None:
        # labels came from the 'label' file - infer the segment length.
        nline = int(nkpt.shape[0] / npair) if npair else nkpt.shape[0]
    for i in range(npair):
        n.append(nline*i)
        n.append(nline*(i+1)-1)
    sympoint = nkpt[n]
    for i in range(sympoint.shape[0]):
        if sympoint[i-1] == sympoint[i]:
            if symname[i] != symname[i-1]:
                symadd = symname[i-1] + "|" + symname[i]
                symname[i-1] = ''
            else:
                symadd = symname[i]
            symname[i] = r'${}$'.format(symadd)
    symname[0] = r'${}$'.format("\\Gamma")
    return nkpt,sympoint,symname
def plot(plottype,file,comp,proj=None,read_kpoint=None,colormap=None):
    """
    Function to create quick plots.

    Parameters:
    -----------
    plottype : str
        Type of plot. 'band' for electronic bandstructure and DOS side by side,
        'phonband' for phonon bandstructure, 'gammaband' for phonon band with lambda projection,
        'a2f' for Eliashberg spectral function.
    file : str
        QE scf input file to generate k-point mesh.
    comp : str
        Compound name.

    Returns:
    ---------
    None

    This function creates various types of plots based on the input parameters.
    """
    input_data = config()
    if plottype == 'band':
        if proj:
            proj_data = np.loadtxt(proj)
        #nkpoint = int(sys.argv[4])
        #kcut = int(sys.argv[5])
        with open("../../input.in","r") as read_inputin:
            inputline = read_inputin.readlines()
        if read_kpoint:
            nkpoint = int(read_kpoint)
            kcut = 0
        else:
            nkpoint = int(inputline[2].split()[0])
            kcut = int(inputline[2].split()[1])
        vasp_line = False
        for line in inputline:
            if "vasp-line" in line:
                vasp_line = True
        if os.path.isfile("VASP_LINE"):
            vasp_line = True
        if os.path.isfile("KPT_OPT"):
            vasp_line = True
        if vasp_line:
            print("vasp-line plot found\n")
            pt_l,sympt,symlb = kptline()
            nkpoint = pt_l.shape[0]
        else:
            _,sympt,symlb,pt_l,_,_ = kpath(file,nkpoint,kcut)
        pt_l = np.array(pt_l)
        # FIX(14): the QE path used to grep "Fermi" (which also matches
        # "Fermi-Dirac") and take the FIRST match - the first ionic step of a
        # relaxation - while the VASP path grepped only "E-fermi" and then
        # called .item() on the result, raising for any relaxation with more
        # than one ionic step.  Both are now handled by one pair of helpers
        # that take the LAST match and accept VASP-5 and VASP-6 wordings.
        if os.path.isfile('scf.out'):
            fermi = fermi_from_qe_out('scf.out')
            band_fermi = int(nelect_from_qe_out('scf.out'))
        elif os.path.isfile('../relax/OUTCAR'):
            fermi = fermi_from_outcar('../relax/OUTCAR')
            band_fermi = int(nelect_from_outcar('../relax/OUTCAR'))
        else:
            raise FileNotFoundError("No scf.out and OUTCAR files present")
        print("Fermi level: {} eV, electrons: {}".format(fermi, band_fermi))
        # FIX(15): `grep LSORBIT INCAR | wc -l` counted "LSORBIT = .FALSE." as
        # spin-orbit ON; the value is parsed instead of merely detected.
        lsorbit = read_lsorbit("INCAR")
        if not lsorbit:
            if not band_fermi % 2 == 0:
                band_fermi = int((band_fermi + 1)/2)
            else:
                band_fermi = int(band_fermi/2)
        sympt = np.array(sympt)
        data=np.loadtxt('{}.dat.gnu'.format(comp))
        check_size = data.shape[1]
        nspin = 1
        #if check_size is 3 then nspin = 2.
        #print(check_size)
        if check_size == 3:
            nspin = 2
        fig,ax_p = plt.subplots()
        nband=int(data.shape[0]/nkpoint)
        print("nkpoints:{}-nband:{}".format(nkpoint,nband))
        with open("band_stat.csv", "w") as band_stat:
            band_stat.write("Band,Emin-Ef,Emax-Ef,Emin,Emax\n")
            for i in range(nband):
                if i == band_fermi - 1:
                    ax_p.plot(pt_l, data[0+nkpoint*i:nkpoint*(i+1)][:,1]-fermi,'r', lw = 0.75)
                    if nspin == 2:
                        ax_p.plot(pt_l, data[0+nkpoint*i:nkpoint*(i+1)][:,2]-fermi,'r--', lw = 0.75)
                    if proj:
                        projection = proj_data[0+nkpoint*i:nkpoint*(i+1)][:,1]
                        plt.scatter(pt_l,data[0+nkpoint*i:nkpoint*(i+1)][:,1]-fermi,c=projection,cmap=colormap,vmin=0,vmax=1)
                        #if nspin == 2:
                        #    projection = proj_data[0+nkpoint*i:nkpoint*(i+1)][:,2]
                        #    plt.scatter(pt_l,data[0+nkpoint*i:nkpoint*(i+1)][:,2]-fermi,c=projection,cmap='Blues',vmin=0,vmax=1)
                else:
                    # if nspin = 2, plot 3rd column with data[][:,2] with black dashed.
                    ax_p.plot(pt_l, data[0+nkpoint*i:nkpoint*(i+1)][:,1]-fermi,'k', lw = 0.75)
                    if nspin == 2:
                        ax_p.plot(pt_l, data[0+nkpoint*i:nkpoint*(i+1)][:,2]-fermi,'k--', lw = 0.75)
                    if proj:
                        projection = proj_data[0+nkpoint*i:nkpoint*(i+1)][:,1]
                        plt.scatter(pt_l,data[0+nkpoint*i:nkpoint*(i+1)][:,1]-fermi,c=projection,cmap=colormap,vmin=0,vmax=1)
                        #if nspin == 2:
                        #    projection = proj_data[0+nkpoint*i:nkpoint*(i+1)][:,2]
                        #    plt.scatter(pt_l,data[0+nkpoint*i:nkpoint*(i+1)][:,2]-fermi,c=projection,cmap='Blues',vmin=0,vmax=1)
                band_stat.write(str(i) + ",")
                band_stat.write(str(np.min(data[0+nkpoint*i:nkpoint*(i+1)][:,1])-fermi) + ",")
                band_stat.write(str(np.max(data[0+nkpoint*i:nkpoint*(i+1)][:,1])-fermi) + ",")
                band_stat.write(str(np.min(data[0+nkpoint*i:nkpoint*(i+1)][:,1])) + ",")
                band_stat.write(str(np.max(data[0+nkpoint*i:nkpoint*(i+1)][:,1])) + "\n")
        ax_p.set_xticks(sympt)
        ax_p.plot([sympt[0], sympt[-1]], [0, 0], 'k--', lw=0.1)
        ax_p.set_xticklabels(symlb)
        min_band = np.min(data[:,1])-fermi
        max_band = np.max(data[:,1])-fermi
        for i in range(sympt.shape[0]):
            ax_p.vlines(sympt[i],min_band-2,max_band+2,color='k',linestyles='dashed')
        ax_p.set_ylabel("E - Ef (eV)")
        #ax_p.set_ylim((min_band-2, max_band+2)) #change this depending on the energy levels.
        ylim = input_data['plot']['ylim']
        if ylim is not None:
            ax_p.set_ylim((ylim[0], ylim[1])) #change this depending on the energy levels.
        else:
            ax_p.set_ylim((-1, 1)) #change this depending on the energy levels.
        fig.set_figheight(8)
        fig.set_figwidth(12)
        #ax_p.set_ylim((-6, 4)) #change this depending on the energy levels.
        ax_p.set_title(comp)
        plt.savefig(comp + "-band.png")
        pylab.savefig(comp + '-band.pdf', format='pdf',bbox_inches='tight')

    if plottype == 'phonband':
        nkpoint = int(sys.argv[4])
        kcut = int(sys.argv[5])
        _,sympt,symlb,pt_l,_,_ = kpath(file,nkpoint,kcut)
        pt_l = np.array(pt_l)
        sympt = np.array(sympt)
        data=np.loadtxt('freq.plot')
        #dataph=np.loadtxt('phonon.dos')
        nband=int(data.shape[0]/nkpoint)
        print("nkpoints:{}-nband:{}".format(nkpoint,nband))
        percmtothz = 33.356
        #dataphx=dataph[:,0]/percmtothz
        fig,ax_p = plt.subplots()
        #color = ['r', 'b', 'g', 'k', 'm', 'y', 'c']
        #markers = ['o', 's', '+', '^', '*', 'D', 'p']
        for i in range(nband):
            ax_p.plot(pt_l, data[0+nkpoint*i:nkpoint*(i+1)][:,1]/percmtothz,'k', lw = 1)
        ax_p.set_ylabel(r'$\omega$'+'(THz)',fontsize=25)
        ax_p.set_xticks(sympt)
        ax_p.set_xticklabels(symlb)
        maxy = np.max(data[:,1])/percmtothz + 2.0
        miny = np.min(data[:,1])/percmtothz - 2.0
        for i in range(sympt.shape[0]):
            ax_p.vlines(sympt[i],0,maxy,color='k',linestyles='dashed')
        maxy = np.max(data[:,1])/percmtothz + 2.0
        miny = np.min(data[:,1])/percmtothz - 1.0
        ax_p.set_ylim((miny,maxy))
        ax_p.set_title(comp)
        #plt.savefig(comp+ "-phonon.png")
        pylab.savefig(comp + '-phonon.pdf', format='pdf',bbox_inches='tight')

    if plottype == 'gammaband':
        with open('lambda.dat') as lambd:
            lines = lambd.readlines()
        dosef=float(lines[2].split('\n')[0].split(' ')[-1])/3289.9146
        #dosef = float(lines[1].split()[11])/3289.9146
        nkpoint = int(sys.argv[4])
        kcut = int(sys.argv[5])
        _,sympt,symlb,pt_l,_,_ = kpath(file,nkpoint,kcut)
        pt_l = np.array(pt_l)
        sympt = np.array(sympt)
        data=np.loadtxt('gamma.plot')[:,1]/1000
        data2=np.loadtxt('freq.plot')[:,1]/33.356
        data3=data2**2.0
        #data3[np.where(data2 < -1)] *= -1
        np.seterr(invalid='ignore')
        lqu=data/(np.pi*dosef*data3)
        #lqu=np.nan_to_num(lqu)
        #lqu_norm = np.linalg.norm(lqu)
        lqu[np.where(data2 < 0.6)] = 0
        #lqu_min = lqu.max()*0.2
        #lqu = lqu/lqu.max()
        print("Setting minimum cutoff for plotting lambda as 10% of maximum value.\n")
        print("Change this if you need otherwise in src/plot.py 259 line\n")
        lqu_min = lqu.max()*0.1
        #print(dosef)
        #lqu = data2
        wph=data2
        nband=int(data.shape[0]/nkpoint)
        #if nband*n < data.shape[0]:
        #    n = 60
        #    nband=int(data.shape[0]/n)
        #    _,sympt,symlb,pt,_,_ = kpath(file,n,kcut)
        #    pt = np.array(pt)
        #    sympt = np.array(sympt)
        print("nkpoints:{}-nband:{}".format(nkpoint,nband))
        percmtothz = 33.356
        fig,ax_p = plt.subplots()
        for i in range(nband):
            x_data = pt_l
            y_data = wph[0+nkpoint*i:nkpoint*(i+1)]
            c_data = lqu[0+nkpoint*i:nkpoint*(i+1)]
            plot_gamma(x_data,y_data,c_data,lqu_min)
        #plt.xlabel("k-vector",fontsize=15)
        ax_p.set_ylabel(r'$\omega$'+'(THz)',fontsize=25)
        ax_p.set_xticks(sympt)
        ax_p.set_xticklabels(symlb)
        maxy = np.max(wph) + 1.0
        miny = np.min(wph) - 0.5
        for i in range(sympt.shape[0]):
            ax_p.vlines(sympt[i],0,maxy,color='k',linestyles='dashed')
        ax_p.set_ylim((miny,maxy))
        ax_p.set_title(comp)
        plt.savefig(comp+ "-gamma.png")
        pylab.savefig(comp + '-gamma.pdf', format='pdf',bbox_inches='tight')
    if plottype == 'a2f':
        fig,ax_p = plt.subplots()
        # FIX(17): lambda used to come from the hard-coded line 12 of
        # lambda.out and alpha2F from the hard-coded column 2 of alpha2F.dat -
        # one arbitrary smearing out of the ten el-ph smearings, and nothing
        # tied the two choices together.  The table is parsed properly now and
        # the smearing is selectable through plot.a2f_smearing (1-based).
        # Leaving that key unset reproduces exactly the old choice, so
        # existing plots do not change silently; the choice is always logged.
        rows, header_index = parse_lambda_out('lambda.out')
        smearing_key = input_data['plot'].get('a2f_smearing')
        ismear, row = select_smearing(rows, header_index, smearing_key)
        lam, omglog, t_c = row[0], row[1], row[2]
        data = np.loadtxt('alpha2F.dat')
        en_data = data[:,0]
        a2f_column = 2 if smearing_key is None else int(smearing_key)
        if a2f_column >= data.shape[1]:
            raise IndexError(
                "alpha2F.dat has {} columns; smearing {} would need column "
                "{}".format(data.shape[1], smearing_key, a2f_column))
        print("a2f: using smearing #{} of {} from lambda.out "
              "(lambda = {}, omega_log = {} K, Tc = {} K) and column {} of "
              "alpha2F.dat".format(ismear, len(rows), lam, omglog, t_c,
                                   a2f_column))
        a2f = data[:,a2f_column]
        #a2f = data[:,2]/en_data
        delw = en_data[1]
        en_data[np.where(en_data == 0.0)] = 0.0000000001
        plt.plot(en_data,a2f,'g',lw=2)
        a2fint = 2.0*delw*a2f/en_data
        a2fint = a2fint.cumsum()
        #lmbd = round(a2fint[-1],4)
        plt.plot(en_data,a2fint,'k--', lw=2)
        if np.max(a2f) > 2.0:
            maxy = np.max(a2f) + 1.0
        else:
            maxy = np.max(a2f) + 0.5
        maxx = np.max(en_data) + 2
        miny = np.min(a2f) - 0.1
        #for tick in ax_p.yaxis.get_major_ticks():
        #    tick.label.set_fontsize(20)
        #for tick in ax_p.yaxis.get_major_ticks():
        #    tick.label.set_fontsize(20)
        ax_p.tick_params(axis='both',labelsize=20)
        fig.set_figheight(6)
        fig.set_figwidth(8)
        ylim = input_data['plot']['ylim']
        xlim = input_data['plot']['xlim']
        if ylim is not None:
            plt.ylim((ylim[0], ylim[1])) #change this depending on the energy levels.
        else:
            plt.ylim((miny,maxy))
        if xlim is not None:
            plt.xlim((xlim[0], xlim[1])) #change this depending on the energy levels.
        else:
            plt.xlim((0.0,maxx))
        plt.xlabel(r'$\omega$' + '(THz)', fontsize=25)
        plt.ylabel(r'$\alpha^2F$'+r'($\omega$)', fontsize=30)
        plt.text(maxx-20,lam+0.5, r'$\lambda$' + "=" + str(round(lam,2)), fontsize=25,fontweight='bold')
        #plt.text(maxx-20,maxy-0.1,r'$\omega_{log}$' + "=" + str(int(omglog)) + " K", fontsize=20,fontweight='bold')
        #plt.text(maxx-20,maxy-1.5, r'$T_c$' + "=" + str(round(Tc,2))+" K", fontsize=20,fontweight='bold')
        plt.savefig(comp+"-a2f.png")
        pylab.savefig(comp+ "-a2f.pdf", format='pdf', bbox_inches='tight')
    if plottype == '':
        print("plot all")
def plot_gamma(xdata,ydata,color,min_):
    """
    Function to plot different sections with respect to different colors.

    Parameters:
    -----------
    xdata : array-like
        x-axis data.
    ydata : array-like
        y-axis data.
    color : array-like
        Projection data.
    min_ : float
        Minimum threshold for binary classification.

    Returns:
    --------
    None

    This function plots different sections of the data with respect to different colors.
    If the color value is less than or equal to the minimum threshold, the section is plotted in black.
    Otherwise, the section is plotted in green with marker sizes proportional to the color values.

    """
    ax_p = plt.gca()
    for i in np.arange(len(xdata) - 1):
        # Change marker size here
        marker_size = color[i]*5
        if color[i] <= min_:
            ax_p.plot([xdata[i],xdata[i+1]], [ydata[i], ydata[i+1]], lw=1.5, color='k')
        else:
            ax_p.plot([xdata[i],xdata[i+1]], [ydata[i], ydata[i+1]], lw=1.5, color='k')
            ax_p.plot([xdata[i],xdata[i+1]], [ydata[i], ydata[i+1]],linestyle='none',marker='o',color='green',markersize=marker_size,fillstyle='none')
def dos_plot(filedos,out='pdos.pdf'):
    """
    Function to plot density of states and partial density of states in different rows.

    Parameters:
    -----------
    filedos : str
        DOS file in .dos format obtained from QE calculations.
    out : str, optional
        Output plot in PDF format. Default is 'pdos.pdf'.

    Returns:
    --------
    None

    This function reads the DOS file and performs DOS and PDOS calculations.
    It requires a 'filedos.in' with ion and orbital contribution from different lines.
    It also uses 'sumpdos.sh' scripts as "sumpdos.sh element orbital" to create 'element-orbital.dat' files.
    For example, for boron and s orbital, element = B, orbital = s, B-s.dat file is created.
    The function then plots the density of states and partial density of states in different rows.

    """
    input_data = config()
    color = ['k', 'r', 'b', 'g','cyan','lightgreen','orange','yellow','lightblue']
    # FIX(18): the header token index [8] is only right for nspin = 1; the
    # nspin = 2 header carries an extra column and the Fermi energy moves.
    # Parse it by name.
    with open(filedos, "r") as dos:
        header = dos.readline()
    match = _DOS_EFERMI.search(header)
    if not match:
        raise ValueError(
            "no 'EFermi = <value>' in the header of {}: {!r}".format(
                filedos, header.strip()))
    fermi = _to_float(match.group(1))
    try:
        with open("filedos.in", "r") as p_dos:
            lines = p_dos.readlines()
    except FileNotFoundError:
        print("filedos.in not found\n")
        sys.exit()
    dos = np.loadtxt(filedos)
    fig,ax_p = plt.subplots(2,1)
    dict_ = {}
    for line in lines:
        dict_[line.split('\n')[0].split(' ')[0]] = line.split('\n')[0].split(' ')[1:]
    datalist = []
    datalist_name = []
    for key in dict_.keys():
        len_dict = len(dict_[key])
        for i in range(len_dict):
            # FIX: os.system("sumpdos.sh <elm> <orb>") interpolated an element
            # name into a shell string.  The helper is a shell script shipped
            # with the package, so it still runs, but through subprocess with
            # an argument list and a logged return code.
            sumpdos = subprocess.run(
                ["sumpdos.sh", str(key), str(dict_[key][i])], check=False)
            if sumpdos.returncode != 0:
                print("sumpdos.sh {} {} returned {}\n".format(
                    key, dict_[key][i], sumpdos.returncode))
            data = np.loadtxt("{}-{}.dat".format(key,dict_[key][i]))
            datalist.append(data)
            datalist_name.append('{} {}'.format(key,dict_[key][i]))
    en_l = datalist[0][:,0]
    for i,_ in enumerate(datalist):
        ax_p[0].plot(en_l-fermi,datalist[i][:,1],color=color[i],label=datalist_name[i],lw=3.0)
    ax_p[1].plot(dos[:,0]-fermi,dos[:,1],'k-',lw=3.0)
    ind = np.logical_and(dos[:,0]-fermi > -8,dos[:,0]-fermi < 4.1)
    dos_range = dos[ind][:,1]
    maxdos = dos_range.max()
    ax_p[0].plot([0,0], [0,maxdos], 'k-.', lw=0.75)
    ax_p[1].plot([0,0], [0,maxdos], 'k-.', lw=0.75)
    if len(datalist) < 5:
        ax_p[0].legend(loc="best",frameon=False,fontsize=20)
    else:
        ax_p[0].legend(ncol=2,loc="best",frameon=False,fontsize=20)
    ax_p[1].set_xlabel(r"E - E$_F$ (eV)", fontsize=30)
    ax_p[0].set_ylabel("PDOS (states/eV/cell)",fontsize=20)
    ax_p[1].set_ylabel("DOS (states/eV/cell)",fontsize=20)
    ax_p[0].set_xticklabels([])
    ax_p[0].tick_params(axis='x', labelsize=30)
    ax_p[0].tick_params(axis='y', labelsize=30)
    ax_p[1].tick_params(axis='x', labelsize=30)
    ax_p[1].tick_params(axis='y', labelsize=30)
    #for tick in ax_p[0].yaxis.get_major_ticks():
    #    tick.label.set_fontsize(30)
    #for tick in ax_p[1].yaxis.get_major_ticks():
    #    tick.label.set_fontsize(30)
    plt.xticks([-8,-6,-4,-2,0,2,4], ["-8", "-6", "-4", "-2", "0", "2", "4"], fontsize = 30)
    #plt.yticks([0, 10, 20, 30, 40], ["0", "10", "20", "30", "40"], fontsize=30)
    fig.set_figheight(8)
    fig.set_figwidth(12)
    plt.subplots_adjust(bottom=0.15)
    # FIX(18): 'ylim' is the ENERGY window; using it as a DOS-height limit is
    # meaningless.  A separate 'dos_ylim' key sets the DOS axis; when it is
    # absent the data-driven default is used, as before.
    dos_ylim = input_data['plot'].get('dos_ylim')
    xlim = input_data['plot']['xlim']
    if dos_ylim is not None:
        ax_p[0].set_ylim((dos_ylim[0], dos_ylim[1]))
        ax_p[1].set_ylim((dos_ylim[0], dos_ylim[1]))
    else:
        ax_p[0].set_ylim(0,maxdos)
        ax_p[1].set_ylim(0,maxdos)
    if xlim is not None:
        ax_p[0].set_xlim(xlim[0],xlim[1])
        ax_p[1].set_xlim(xlim[0],xlim[1])
    else:
        ax_p[0].set_xlim(-8,4)
        ax_p[1].set_xlim(-8,4)
    plt.savefig(out)
def band_wann_plot(fileout='plot.pdf'):
    """
    Function to plot Wannier interpolated bandstructure.

    Parameters:
    -----------
    fileband : str, optional
        Bandstructure file. Default is 'ex_band.dat'.
    fileout : str, optional
        Output file. Default is 'plot.pdf'.

    Returns:
    --------
    None

    This function reads the bandstructure file and plots the Wannier interpolated bandstructure.
    It extracts data from the provided files and plots the bandstructure accordingly.
    The function also checks for bands below the Fermi level and prints a message if found.

    """
    print("plotting wannier band\n")
    band_files = glob.glob("*_band.dat")
    if not band_files:
        raise FileNotFoundError("no *_band.dat file in the current directory")
    fileband = band_files[0]
    print(fileband)
    # FIX: the three `cat ... | awk ... > file` calls are a plain read of the
    # labelinfo file; the 'label' and 'lbpoint' files are still written,
    # because kpt_labels() falls back to 'label'.
    label_rows = []
    for info in sorted(glob.glob("*_band.labelinfo.dat")):
        with open(info, "r") as read_info:
            label_rows += [line.split() for line in read_info if line.split()]
    if not label_rows:
        raise FileNotFoundError("no *_band.labelinfo.dat file with content")
    symlb = [row[0] for row in label_rows]
    sympt = np.array([float(row[2]) for row in label_rows])
    nkpt = int(float(label_rows[-1][1]))
    with open("label", "w") as write_label:
        write_label.write("".join(name + "\n" for name in symlb))
    with open("lbpoint", "w") as write_point:
        write_point.write("".join(str(point) + "\n" for point in sympt))
    with open("n.dat", "w") as write_n:
        write_n.write(str(nkpt) + "\n")
    # FIX(14): same Fermi-level defect as in plot(): the QE branch took the
    # FIRST "Fermi" match (and matched "Fermi-Dirac"), and the VASP branch
    # here handled only the VASP-6 wording while plot() handled only VASP-5.
    fermi = read_fermi(qe_out="scf.out", outcar="OUTCAR")
    print("Fermi level: {} eV".format(fermi))
    color = ['r', 'b', 'r','g']
    linestyle = ['solid', 'dashed', 'dotted','-.']
    _,ax_p = plt.subplots()
    data=np.loadtxt(fileband)
    pt_k = data[:nkpt,0]
    nband=int(data.shape[0]/nkpt)
    print("nkpoints:{}-nband:{}".format(nkpt,nband))
    for i in range(nband):
        ax_p.plot(pt_k, data[0+nkpt*i:nkpt*(i+1)][:,1], lw = 0.5, linestyle=linestyle[0],color=color[0])
        if np.any(data[0+nkpt*i:nkpt*(i+1)][:,1] < fermi):
            print("Band below Fermi level: {} \n".format(i+1))
    ax_p.set_xticks(sympt)
    ax_p.set_xticklabels(symlb)
    #ax.set_ylim(fermi-1,fermi+1)
    #for i in range(sympt.shape[0]):
    #    ax.vlines(sympt[i],11.0,18.0,color='k',linestyles='dashed')
    ax_p.set_ylabel("E (eV)")
    ax_p.plot([sympt[0],sympt[-1]],[fermi,fermi],'k--', lw=0.5)
    pylab.savefig(fileout, format='pdf',bbox_inches='tight')

def plot_projection(scf_file,projection_file,phonon_freq,outfile,nkpt):
    """
    Function to plot atomic projection on phonon bandstructure.

    Parameters:
    -----------
    scf_file : str
        QE scf.in file to extract structure.
    projection_file : str
        File containing atomic projection data.
    phonon_freq : str
        File containing phonon dispersion data.
    outfile : str
        Output plot file.
    nkpt : int
        Number of q points.
    proj_cutoff : float, optional
        Cutoff to apply filter. Plotting only colors that have projection larger than the cutoff.
        If not provided, it will be calculated based on the number of atoms.

    Returns:
    --------
    None

    This function reads the necessary files and plots the atomic projection on the phonon bandstructure.
    It calculates the projection cutoff based on the number of atoms if not provided explicitly.
    The function generates a plot showing the atomic projection on the phonon dispersion.
    """
    input_data = config()
    _,sympt,symlb,pt_l,_,_ = kpath(scf_file,nkpt,0)
    pt_l = np.array(pt_l)
    sympt = np.array(sympt)
    #print(pt_l,sympt,symlb)
    data = np.loadtxt(phonon_freq)
    nband = data.shape[1]
    proj = np.loadtxt(projection_file)
    proj_list = []
    nat = int(proj.shape[0]/nkpt)
    # Define projection cutoff for systems with
    # different number of ions.
    atomproj = input_data['plot']['atomproj']
    if nat < 4:
        if 'atomproj' in input_data['plot'].keys() and atomproj is not None:
            proj_cutoff = atomproj
        else:
            print("atomproj key not found in config.json. Using default value\n")
            proj_cutoff = 0.6
    elif nat == 4:
        if 'atomproj' in input_data['plot'].keys() and atomproj is not None:
            proj_cutoff = atomproj
        else:
            print("atomproj key not found in config.json. Using default value\n")
            proj_cutoff = 0.5
    else:
        print("Plotting is available for systems with maximum 4 ions. Exiting\n")
        sys.exit()
    print("Cutoff filter for projection: {}\n".format(proj_cutoff))
    print("Change atomproj value in config.json for different cutoff\n")
    for i in range(nat):
        proj_list.append(proj[i*nkpt:(i+1)*nkpt,:])
    percmto_thz = 33.356
    proj_array = np.zeros_like(proj_list[0])
    color_list = ['red','blue','green','cyan','k']
    for i,proj_i in enumerate(proj_list):
        proj_array[np.where(proj_i > proj_cutoff)] = i
    _,ax_p = plt.subplots()
    for i in range(1,nband):
        xdata = pt_l
        ydata = data[:,i]/percmto_thz
        cdata = proj_array[:,i-1]
        for j in np.arange(len(xdata) - 1):
            ax_p.plot([xdata[j], xdata[j+1]], [ydata[j], ydata[j+1]], lw=1.5, color=color_list[int(cdata[j])])
    maxy = np.max(data)/percmto_thz + 2.0
    miny = np.min(data)/percmto_thz - 1.0
    plt.ylabel(r'$\omega$ (THz)',fontsize=20)
    for i in range(sympt.shape[0]):
        ax_p.vlines(sympt[i],0,maxy,color='k',linestyles='dashed')
    plt.xticks(sympt,symlb)
    plt.ylim(miny,maxy)
    plt.xlim(sympt[0],sympt[-1])
    plt.savefig(outfile)
def write_filedos(comp):
    """
    Function to write a 'filedos.in' file required for PDOS plot.

    Parameters:
    -----------
    comp : str
        Composition.

    Returns:
    --------
    None

    This function writes the 'filedos.in' file based on the composition provided.
    It extracts the elements and their electronic structure to determine the orbital contributions for PDOS.
    """
    if "-" in comp:
        comp = comp.split("-")[0]
    comp = Composition(comp)
    nelm = comp.elements
    with open('filedos.in', 'w') as write_pdos:
        for elm in nelm:
            els = elm.full_electronic_structure
            orb_list = []
            write_pdos.write(elm.symbol + " ")
            for orb in els:
                if 's' in orb and 's' not in orb_list:
                    orb_list.append('s')
                elif 'p' in orb and 'p' not in orb_list:
                    orb_list.append('p')
                elif 'd' in orb and 'd' not in orb_list:
                    orb_list.append('d')
                elif 'f' in orb and 'f' not in orb_list:
                    orb_list.append('f')
                else:
                    continue
            for i,orb in enumerate(orb_list):
                if i < len(orb_list) - 1:
                    write_pdos.write(orb + " ")
                else:
                    write_pdos.write(orb)
            write_pdos.write("\n")
def dos_plot_vasp(outfile="pdos.pdf"):
    """
    Function to plot DOS and partial DOS (pDOS) using the vasprun.xml file from VASP using the pymatgen package.

    Parameters:
    -----------
    outfile : str, optional
        Output plot file name (default is "pdos.pdf").

    Returns:
    --------
    None

    This function reads the 'filedos.in' file to determine which elements and orbitals to include in the pDOS plot.
    It uses the Vasprun object from the pymatgen package to parse the vasprun.xml file.
    The function plots the Total DOS and pDOS for the specified elements and orbitals.

    Note:
    -----
    'filedos.in' should be present in the current directory to specify the elements and orbitals for pDOS plotting.
    """
    input_data = config()
    try:
        with open("filedos.in", "r") as p_dos:
            lines = p_dos.readlines()
    except FileNotFoundError:
        print("filedos.in not found\n")
        sys.exit()
    try:
        result = Vasprun('vasprun.xml',parse_dos=True)
    # FIX: a bare `except:` here swallowed KeyboardInterrupt and every
    # unrelated failure; the retry is only meant to cover a POTCAR problem.
    except (OSError, ValueError, KeyError, IndexError) as err:
        print("Vasprun('vasprun.xml', parse_dos=True) failed ({}: {}); "
              "retrying without the POTCAR\n".format(type(err).__name__, err))
        result = Vasprun('vasprun.xml', parse_potcar_file=False)
    complete_dos = result.complete_dos
    nspin = len(complete_dos.densities.keys())
    plotter = DosPlotter()
    xlim1 = input_data['plot']['xlim']
    # FIX(18): 'ylim' is the ENERGY window - reusing it as a DOS height is
    # meaningless.  'dos_ylim' is the DOS-axis key; 'ylim' is no longer used
    # for the DOS height at all.
    dos_ylim = input_data['plot'].get('dos_ylim')
    # FIX(18): xlim1[0] / ylim1[0] used to be indexed BEFORE the
    # `is not None` tests below, so a null 'xlim' raised TypeError here.
    xmin = xmax = None
    if xlim1 is not None:
        xmin, xmax = xlim1[0], xlim1[1]
    ymin = ymax = None
    if dos_ylim is not None:
        ymin, ymax = dos_ylim[0], dos_ylim[1]
    if nspin == 1:
        plotter.add_dos('Total DOS', result.tdos)
        dict_ = {}
        for line in lines:
            dict_[line.split('\n')[0].split(' ')[0]] = line.split('\n')[0].split(' ')[1:]
        for key in dict_.keys():
            value = dict_[key]
            len_value = len(value)
            pdos_ion = complete_dos.get_element_spd_dos(key)
            for i in range(len_value):
                plotter.add_dos("{}({})".format(key,value[i]),pdos_ion[OrbitalType[value[i]]])
        if xlim1 is not None:
            plot_axis = plotter.get_plot(xlim=(xmin,xmax))
        else:
            plot_axis = plotter.get_plot(xlim=(-4,4))
        plot_axis.legend(loc="upper right")
        plot_axis.set_ylabel("Density of States (states/eV)",fontweight='bold',fontsize=25)
        plot_axis.set_xlabel(r"E - E$_F$ (eV)",fontweight='bold',fontsize=25)
        plt.savefig(outfile, dpi=500)
    elif nspin == 2:
        dos = complete_dos
        energies = dos.energies - dos.efermi
        spin = list(dos.densities.keys())
        dos_spin_up = dos.densities[spin[0]]
        dos_spin_down = -1*dos.densities[spin[1]]
        plt.figure(1)
        plt.plot(energies, dos_spin_up, label="Up", color='b')
        plt.plot(energies, dos_spin_down, label="Down", color='r')
        plt.axhline(0, color='black', linestyle='--', linewidth=0.5)
        plt.axvline(0, color='black', linestyle='--', linewidth=0.5)
        plt.legend(bbox_to_anchor=(1.05, 1),loc="upper right",frameon=False)
        plt.ylabel("DOS (states/eV)",fontweight='bold',fontsize=25)
        plt.xlabel(r"E - E$_F$ (eV)",fontweight='bold',fontsize=25)
        if xlim1 is not None:
            plt.xlim(xmin,xmax)
        else:
            plt.xlim(-8,4)
        # FIX(18): was driven by the energy 'ylim'; now by 'dos_ylim'.
        if dos_ylim is not None:
            ymax1 = ymax/2.0
            plt.ylim(-1*ymax1 - 1,ymax1 + 1)
        else:
            plt.ylim(-40,40)
        plt.savefig("pdos-spin-resolved.pdf", dpi=500)
        plt.figure(2)
        basecolor = list(mcolors.BASE_COLORS.keys())[:-1]
        color=basecolor+list(mcolors.TABLEAU_COLORS.keys())
        total_dos = result.complete_dos.get_densities()
        energies = result.complete_dos.energies
        efermi = result.efermi
        plt.plot(energies-efermi, total_dos, 'k-', lw=1.75, label="Total DOS")
        elm_list = result.complete_dos.structure.composition.elements
        ndos_data = 0
        for elm in elm_list:
            elm = str(elm)
            orb_list = list(result.complete_dos.get_element_spd_dos(elm).keys())
            for i, orb in enumerate(orb_list):
                x  = result.complete_dos.get_element_spd_dos(elm)[orb].get_densities()
                k = ndos_data
                plt.plot(energies-efermi, x, color=color[k], lw=1.75, label=f"{elm}-{orb}")
                ndos_data += 1
        # FIX(18): was driven by the energy 'ylim'; now by 'dos_ylim'.
        if dos_ylim is not None:
            plt.ylim(ymin,ymax)
        else:
            plt.ylim(-5,40)
        if xlim1 is not None:
            plt.xlim(xmin,xmax)
        else:
            plt.xlim(-8,4)
        # FIX(18): int(n/4) gives legend(ncol=0) - a ValueError - whenever
        # fewer than four curves were drawn.
        ncol = max(1, int(ndos_data/4))
        plt.legend(loc="upper left",ncol=ncol,fontsize='x-small',handletextpad=0.5,labelspacing=0.2)
        plt.ylabel("DOS (states/eV)",fontweight='bold',fontsize=25)
        plt.xlabel(r"E - E$_F$ (eV)",fontweight='bold',fontsize=25)
        plt.savefig(outfile, dpi=500)
    else:
        print("nspin should be 1 or 2\n")
def band_plot_vasp_line(mpid,compound):
    """
    Function to plot the bandstructure in line mode.

    Parameters:
    -----------
    mpid : str
        Materials ID.
    compound : str
        Name of the compound.

    Returns:
    --------
    None

    This function generates a PDF plot of the bandstructure in line mode using data from the vasprun.xml file.
    The plot is saved with the filename format: "<mpid>-<compound>-band.pdf".
    The band structure is obtained using the Vasprun object from the pymatgen package.
    The BSPlotter object is used to generate the band structure plot with specified settings.

    Note:
    -----
    Ensure that the vasprun.xml file containing band structure data is present in the current directory.
    """
    input_data = config()
    outfile = mpid + "-" + compound + "-" + "band.pdf"
    vaspout = Vasprun("vasprun.xml")
    bandstr = vaspout.get_band_structure(line_mode=True)
    plt1 = BSPlotter(bandstr)
    ylim1 = input_data['plot']['ylim']
    if ylim1 is not None:
        plt2 = plt1.get_plot(ylim=[ylim1[0],ylim1[1]],vbm_cbm_marker=True,zero_to_efermi=True)
    else:
        plt2 = plt1.get_plot(ylim=[-1.8,1.8],vbm_cbm_marker=True,zero_to_efermi=True)
    plt2.legend('',frameon=False)
    plt2.figure.savefig(outfile)
def main():
    """
    Main function to execute plotting tasks based on command-line arguments.

    This function reads command-line arguments to determine the type of plot and other required parameters.
    It performs different plotting tasks based on the specified plot type.

    Command-line Arguments:
    -----------------------
    plottype : str
        Type of plot to generate.
    mpid : str
        Materials ID.
    comp : str
        Name of the compound.

    Returns:
    --------
    None

    Plotting Tasks:
    ---------------
    - For 'pdos' plot type:
        - Checks if 'filedos.in' exists, creates one if not found.
        - Determines the type of calculation (QE or VASP).
        - Calls 'dos_plot_vasp' function if VASP calculation is detected, otherwise 'dos_plot' function.

    - For 'wann_band' plot type:
        - Calls 'band_wann_plot' function.

    - For 'phonproj' plot type:
        - Reads the number of k-points from command-line argument.
        - Calls 'plot_projection' function with appropriate parameters.

    - For other plot types:
        - Calls 'plot' function with specified plot type, file name, and compound name.

    Note:
    -----
    Ensure proper command-line arguments are provided.
    Only QE and VASP outputs are allowed for plotting tasks.
    """
    plottype = sys.argv[1]
    mpid = sys.argv[2]
    comp = sys.argv[3]
    if os.path.isfile('scf.in'):
        filename = 'scf.in'
    elif os.path.isfile('POSCAR'):
        filename = 'POSCAR'
    else:
        print("No scf.in and POSCAR exists\n")
    #input_file = espresso.read_espresso_in('scf.in')
    if plottype == 'pdos':
        if not os.path.isfile("filedos.in"):
            print("filedos.in file not found, creating one\n")
            write_filedos(comp)
        if os.path.isfile('POSCAR'):
            print("POSCAR present. searching vasprun.xml..\n")
            dos_plot_vasp("pdos-{}.pdf".format(comp))
        elif os.path.isfile('scf.in'):
            print("scf.in present. searching {}.dos and *pdos.pdos* files..\n".format(comp))
            dos_plot("{}.dos".format(comp),"pdos-{}.pdf".format(comp))
        else:
            print("Only QE and VASP outputs are allowed\n")
    elif plottype == 'wann_band':
        band_wann_plot()
    elif plottype == 'phonproj':
        nkpt = int(sys.argv[4])
        plot_projection(filename,"phonon-{}.proj.gp".format(comp),"{}.freq.gp".format(comp),"plot-proj-{}-{}.pdf".format(mpid,comp),nkpt)
    #elif plottype == 'line':
    #    band_plot_vasp_line(mpid,comp)
    else:
        plot(plottype,filename,comp)
if __name__ == "__main__":
    main()
