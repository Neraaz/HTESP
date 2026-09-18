#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to prepare QE input files"""
import copy
import os
import re
import warnings
from collections import OrderedDict
from pathlib import Path
import numpy as np
import scipy.linalg as alg
from ase.io import espresso,cif
from ase.cell import Cell
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io import pwscf
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.periodic_table import Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from htesp.kpath import kpath
from htesp.check_json import config
from htesp.config import api_key, require_api_key


def mprester(key):
    """Return an :class:`MPRester` for ``key``, importing ``mp_api`` on first use.

    FIX(20): ``from mp_api.client import MPRester`` at module scope made every
    importer of this module -- ``cif_to_gsinput``, ``elastic``, ``site_subs``,
    ``poscar_to_vasp``, ``qe_input`` -- load ``mp_api``, which in turn imports
    ``deltalake`` and its compiled extension.  Where that extension cannot be
    loaded the interpreter *aborts* (``Fatal Python error: Aborted``) rather
    than raising, so the whole process dies during the import -- including the
    test runner, and including commands that only ever write input files.

    Nothing but an actual Materials Project query needs ``mp_api``, so it is
    imported here, at the one point where a query is about to happen.
    """
    try:
        from mp_api.client import MPRester  # noqa: PLC0415 -- deliberately lazy
    except ImportError as exc:              # pragma: no cover - env dependent
        raise ImportError(
            "Materials Project access needs the mp_api package: "
            "pip install mp_api"
        ) from exc
    return MPRester(key)


class _LazyConfig(dict):
    """The configuration, resolved on first use rather than at import time.

    ``input_data = config()`` at module scope froze the configuration of
    whatever directory happened to be current when the module was first
    imported.  Under the per-material process pool -- and under any
    ``os.chdir`` -- that is the wrong directory.  This proxy resolves (and
    re-resolves, when the working directory changes) on the first attribute
    access instead, while keeping the ``input_data[...]`` spelling every
    function in this module already uses.
    """

    _resolved_in = None

    def _refresh(self):
        here = os.getcwd()
        if self._resolved_in != here:
            dict.clear(self)
            dict.update(self, config())
            self._resolved_in = here
        return self

    def __getitem__(self, key):
        return dict.__getitem__(self._refresh(), key)

    def get(self, key, default=None):
        return dict.get(self._refresh(), key, default)

    def __contains__(self, key):
        return dict.__contains__(self._refresh(), key)

    def keys(self):
        return dict.keys(self._refresh())

    def items(self):
        return dict.items(self._refresh())


input_data = _LazyConfig()

# FIX(5): ``getkpt`` symmetrised with symprec=0.1 and pymatgen's default
# ``international_monoclinic`` while ``setting_qeinput`` used symprec=0.01 with
# ``international_monoclinic=False``, so the k-mesh was computed for a
# different cell than the one written to scf-<mpid>.in.  Both go through these
# two constants now.
#: Summary fields that are *nested documents* rather than scalars.
#: FIX(31): ``setting()`` asked for ``available_fields[:-29]``, which still
#: includes these two.  emmet-core validates every requested sub-document, and
#: the payload the API returns for them no longer carries the fields its model
#: declares -- ``dos.elemental.Ru.total.1.task_id ... Field required``,
#: ``bandstructure.setyawan_curtarolo.equivalent_labels ... Field required``.
#: One unusable sub-document rejects the *whole* SummaryDoc, so
#: ``mainprogram download`` raised ValidationError for every material that has
#: band-structure or DOS data -- which is nearly all of them -- while the few
#: without it succeeded.  Nothing in HTESP reads either field: band structures
#: and densities of states are computed here, not fetched.
UNREQUESTABLE_FIELDS = ("bandstructure", "dos")

#: symmetry tolerance used everywhere a standard cell is derived in this module
SYMPREC = 0.01
#: 2010_SC-style monoclinic setting (alpha < 90, beta = gamma = 90); the
#: international setting puts the unique angle elsewhere, giving a cell that
#: does not match the k-mesh.
INTERNATIONAL_MONOCLINIC = False

# FIX(1c): ``starting_magnetization`` used to be written as 0 for every species
# of an undecorated structure under ``nspin = 2`` -- a non-magnetic calculation
# at twice the cost.  This is the last-resort value used when neither the
# structure nor ``config['magmom']['magmom']`` says anything about an element.
DEFAULT_MAGMOM = 0.5

# FIX(2): kept at module level so ``getecut_sssp`` can fall back to it when an
# element is missing from ``config['pseudo']['PSEUDO']``.
#: Wavefunction cutoffs in Ry, mirroring ``pseudo.PSEUDO`` of the packaged
#: config.json so the fallback can never contradict the shipped table (a
#: test asserts they are identical).  New elements are taken from SSSP 1.3.0
#: PBE efficiency; the entries that predate it are left as they were.
SSSP_EFFICIENCY = {
    'H': 60, 'He': 50, 'Li': 40, 'Be': 40, 'B': 35, 'C': 45, 'N': 60, 'O': 60,
    'F': 45, 'Ne': 50, 'Na': 40, 'Mg': 30, 'Al': 30, 'Si': 30, 'P': 30, 'S': 35,
    'Cl': 40, 'Ar': 60, 'K': 60, 'Ca': 30, 'Sc': 40, 'Ti': 35, 'V': 35, 'Cr': 40,
    'Mn': 65, 'Fe': 90, 'Co': 45, 'Ni': 45, 'Cu': 55, 'Zn': 40, 'Ga': 70, 'Ge': 40,
    'As': 35, 'Se': 30, 'Br': 30, 'Kr': 45, 'Rb': 30, 'Sr': 30, 'Y': 35, 'Zr': 30,
    'Nb': 40, 'Mo': 35, 'Tc': 30, 'Ru': 35, 'Rh': 35, 'Pd': 45, 'Ag': 50, 'Cd': 60,
    'In': 50, 'Sn': 60, 'Sb': 40, 'Te': 30, 'I': 35, 'Xe': 60, 'Cs': 30, 'Ba': 30,
    'La': 40, 'Ce': 50, 'Pr': 40, 'Nd': 40, 'Pm': 40, 'Sm': 40, 'Eu': 40, 'Gd': 40,
    'Tb': 40, 'Dy': 40, 'Ho': 40, 'Er': 40, 'Tm': 40, 'Yb': 40, 'Lu': 45, 'Hf': 50,
    'Ta': 45, 'W': 30, 'Re': 30, 'Os': 40, 'Ir': 55, 'Pt': 35, 'Au': 45, 'Hg': 50,
    'Tl': 50, 'Pb': 40, 'Bi': 45, 'Po': 30, 'At': 30, 'Rn': 30, 'Fr': 30, 'Ra': 30,
    'Ac': 30, 'Th': 30, 'Pa': 30, 'U': 30, 'Np': 30, 'Pu': 30, 'Am': 30, 'Cm': 30,
    'Bk': 60, 'Cf': 60, 'Es': 60, 'Fm': 60, 'Md': 60, 'No': 60, 'Lr': 60,
}


#: start of a QE "card" (everything that is not part of a namelist)
CARD_RE = re.compile(
    r"^\s*(ATOMIC_SPECIES|ATOMIC_POSITIONS|ATOMIC_VELOCITIES|ATOMIC_FORCES|"
    r"K_POINTS|CELL_PARAMETERS|OCCUPATIONS|CONSTRAINTS|SOLVENTS|HUBBARD)\b")

#: &SYSTEM entries that are rebuilt when a magnetic input is written
SYSTEM_KEY_RE = re.compile(
    r"^\s*(nspin|ntyp|nat|starting_magnetization\s*\(\s*\d+\s*\))\s*=")


def read_lines(path):
    """Return the lines of ``path`` with their trailing newline stripped."""
    return Path(path).read_text().splitlines()


def write_lines(path, lines):
    """Write ``lines`` (given without trailing newlines) to ``path``."""
    text = "\n".join(lines)
    Path(path).write_text(text + "\n" if text else "")


def section(lines, start, end):
    """``sed -n '/start/,/end/p' | sed '$d'`` without a shell.

    Returns the lines from the first one containing ``start`` up to, but not
    including, the first following line containing ``end``.  An empty list is
    returned when ``start`` never appears.
    """
    try:
        first = next(i for i, ln in enumerate(lines) if start in ln)
    except StopIteration:
        return []
    for j in range(first + 1, len(lines)):
        if end in lines[j]:
            return lines[first:j]
    return lines[first:]


def bare_element(label):
    """``'Fe,spin=5'`` / ``'Fe1'`` -> ``'Fe'``; anything else is returned as is."""
    label = str(label).split(",")[0].strip()
    match = re.match(r"^([A-Z][a-z]?)", label)
    return match.group(1) if match else label


def remove_files(*paths):
    """Delete each path that exists, quietly skipping the ones that do not."""
    for path in paths:
        try:
            Path(path).unlink()
        except FileNotFoundError:
            continue
        except OSError as exc:
            print("could not remove {}: {}".format(path, exc))
def pos_to_kpt(structure_filename,kpoint_density):
    """
    Obtain k-point mesh from a structure file.

    Parameters:
    ----------------------
    structure_filename : str
        Structure file (e.g., QE scf.in or VASP POSCAR).
    kpoint_density : float
        K-point density.

    Returns:
    ----------------------
    kmesh : list
        K-point mesh according to the k-point density.
    """
    kptsp = kpoint_density
    with open(structure_filename,"r") as read_struc:
        lines = read_struc.readlines()
    # Get cell vectors
    line = lines[1].split()
    latscl = float(line[0])
    ain = np.zeros((3, 3))
    amat = np.zeros((3, 3))
    bmat = np.zeros((3, 3))
    anorm = np.zeros(3)
    for i in range(3):
        line = lines[2 + i].split()
        for j in range(3):
            ain[i][j] = float(line[j]) * latscl
            amat[j][i] = ain[i][j]
        anorm[i] = np.sqrt(ain[i][0] ** 2 + ain[i][1] ** 2 + ain[i][2] ** 2)
    bmat = alg.inv(amat)
    bnorm = np.zeros(3)
    bnorm = alg.norm(bmat,axis=1)
    kratio = [bnorm[i] / bnorm[0] for i in range(3)]
    klat = bnorm[0] / kptsp
    kmesh = [int(kratio[i] * klat + 0.5) if int(kratio[i] * klat + 0.5) != 0 else 1 for i in range(3)]
    kratio = [bnorm[i] / bnorm[0] for i in range(3)]
    klat = bnorm[0] / kptsp
    kmesh = [int(kratio[i] * klat + 0.5) if int(kratio[i] * klat + 0.5) != 0 else 1 for i in range(3)]
    return kmesh
class MpConnect:
    """
    Class for connecting to Materials Project (MP) via MP API, extracting properties,
    and preparing input files for Quantum Espresso (QE) calculations.

    Parameters:
    --------------
    key : str, optional
        MP API KEY. Use your own key. If not provided, it should be available in config.json
        or ../../config.json file.

    Attributes:
    --------------
    prop : list
        List of available properties.
    mpid : str
        Materials ID.
    comp : str
        Compound name.
    data : dict
        Dictionary containing materials data.
    prefix : str
        Prefix for the compound.
    ecutwfc : float
        Kinetic energy cutoff for wavefunction.
    ecutrho : float
        Kinetic energy cutoff for charge density.
    kpt : list
        K-point grid.
    structure : Structure
        Structure object.
    pseudo_dir : str
        Path to pseudopotential files.
    outdir : str
        Directory for output files.
    evenkpt : tuple
        K-grid with an even number of points in all directions.
    kpshift : list
        K-point grid shifts.
    kptype : str
        K-point grid type.
    calc : str
        Type of calculation.
    smear : float
        Degauss value for smearing.
    smear_type : str
        Type of smearing.
    etot_conv_thr : float
        Total energy convergence threshold.
    forc_conv_thr : float
        Force convergence threshold.
    conv_thr : float
        Convergence threshold.
    dict_element : dict
        Dictionary containing kinetic energy cutoffs for different elements.
    comp_list : list
        List of elements in the compound.

    Example of using MpConnect class:
    >>> obj = MpConnect()  # Initialize MpConnect object
    >>> obj.setting('mp-763')  # Set Materials ID to 'mp-763'
    >>> obj.getkpt()  # Get k-points for the structure
    >>> obj.maxecut_sssp()  # Calculate maximum recommended energy cutoff using SSSP
    >>> obj.setting_qeinput()  # Set up Quantum Espresso input files based on the retrieved data
    """
    # Intialize the class
    def __init__(self):
        # The API key comes from $MP_API_KEY, then ~/.config/htesp/credentials,
        # then config.json -- never from an unguarded dictionary lookup that
        # left ``key`` unbound when no config.json was next to the caller.
        # ``api_key`` (not ``require_api_key``) is used here because most users
        # of this class -- elastic.py, site_subs.py, poscar_to_vasp.py -- only
        # want the input-writing half and never contact Materials Project;
        # ``setting`` raises through ``require_api_key`` when a key is needed.
        self.key = api_key(input_data)
        self.mpr = mprester(self.key) if self.key else None
        self.prop = None
        self.mpid = None
        self.comp = None
        self.data = None
        self.prefix = None
        self.ecutwfc = 0
        self.ecutrho = 0
        self.kpt = None
        self.structure = None
        self.pseudo_dir = ""
        self.outdir = ""
        self.evenkpt = None
        self.kpshift = None
        self.kptype = ""
        self.calc = ""
        self.smear = None
        self.smear_type = ""
        self.etot_conv_thr = None
        self.forc_conv_thr = None
        self.conv_thr = None
        self.dict_element = {}
        self.comp_list = None
    def setting(self,comp):
        """
        Initialize the process with a Materials ID.

        Parameters:
        ---------------------
        comp : str
            Materials ID.

        Returns:
        ---------------------
        comp : str
            Materials ID.

        Notes:
        ---------------------
        This function initializes the process with a Materials ID. It retrieves the data related to the given ID from the Materials Project (MP) database, sets up necessary parameters, and prepares the structure for further calculations.
        """
        # FIX(7): kpt, ecutwfc/ecutrho, comp_list and prefix survived from the
        # previous material when one process looped over several ids, so a
        # material whose own lookup failed silently inherited the last one's
        # k-mesh and cutoffs.  Reset them before anything else.
        self.kpt = None
        self.evenkpt = None
        self.kptype = ""
        self.kpshift = None
        self.ecutwfc = 0
        self.ecutrho = 0
        self.comp_list = None
        self.prefix = None
        self.structure = None
        self.data = None
        self.mpid = None
        # Set the Material ID
        self.comp=comp
        if self.mpr is None:
            # raises with the instructions for setting $MP_API_KEY
            self.key = require_api_key(input_data)
            self.mpr = mprester(self.key)
        # Retrieve available properties from Materials Project database
        # We only store properties except for last 29 elements
        # Mostly related to elastic properties, absent for many systems
        # FIX(31): drop the nested documents emmet-core cannot validate;
        # see UNREQUESTABLE_FIELDS for what that failure looks like.
        self.prop = [field
                     for field in self.mpr.materials.summary.available_fields[:-29]
                     if field not in UNREQUESTABLE_FIELDS]
        # ``config()`` always returns the full schema now, so the settings are
        # read unconditionally rather than only when a config.json sits in the
        # working directory.
        d = input_data["download"]
        # Append additional properties specified in the settings to the properties list
        for prop in d['element']['prop']:
            if prop not in self.prop:
                self.prop.append(prop)
        #self.data = self.mpr.materials.summary.get_data_by_id(self.comp,fields=self.prop).dict()
        # Fetch data for the given Materials ID from the Materials Project database
        # FIX(34): ``[0]`` on an empty result is "IndexError: list index out of
        # range", which says nothing about which id failed or why.  An id that
        # Materials Project does not know -- a typo, a deprecated entry, or an
        # identifier from another database (examples/QE/tutorial6 ships an
        # mpid-list.in full of ``oqmd-*`` ids) -- now names itself.
        found = self.mpr.materials.summary.search(material_ids=[self.comp],
                                                  fields=self.prop)
        if not found:
            raise LookupError(
                "Materials Project returned no entry for {!r}.  Check the id in "
                "the tracking file: ids from OQMD ('oqmd-...') or AFLOW "
                "('aflow:...') cannot be looked up here, and a deprecated "
                "mp-id has to be replaced with its successor.".format(self.comp))
        self.data = found[0]
        self.structure = self.data.structure
        self.data = self.data.dict()
        # Extract relevant information from the fetched data
        self.mpid = self.data['material_id']
        if self.mpid != self.comp:
            self.mpid = self.comp
        symbol_comp = ""
        elm = list(self.data['composition'].keys())
        count = list(self.data['composition'].values())
        # Construct a prefix based on elemental composition
        for i,_ in enumerate(elm):
            symbol_comp += elm[i]+str(int(count[i]))
        self.prefix = symbol_comp
        #print("******************************************\n")
        #print("Use property method to extract the following info\n")
        #print("property('name') or get_properties(['name1', 'name2', ..])\n")
        #print("******************************************\n")
        return comp
    def get_prop_list(self):
        """
        Print the list of properties available when downloading the data.

        Returns:
        ---------------------
        prop_list : dict_keys
            A list of keys representing the available properties.

        Notes:
        ---------------------
        This function retrieves the list of properties available for download from the Materials Project (MP) database for the current Materials ID (mpid). It returns a list of keys that represent the available properties that can be downloaded and accessed for analysis or further processing.
        """
        #return self.mpr.summary.get_data_by_id(self.mpid).dict().keys()
        return self.prop
        #return self.mpr.materials.summary.search(material_ids=[self.mpid],fields=self.prop).dict().keys()
    def download(self,filetype='cif'):
        """
        Download the structure file.

        Parameters:
        ---------------
        filetype : str, optional
            File format for downloading. Default is 'cif'.

        Notes:
        ---------------
        This function downloads the structure file for the current Materials ID (mpid) from the Materials Project (MP) database. It saves the file in the specified format, with the default format being CIF ('.cif'). Other supported formats may be specified as needed.
        """
        if filetype == 'cif':
            unsym_struc = self.structure
            CifWriter(unsym_struc, symprec=0.1).write_file('{}.cif'.format(self.mpid))
    def property(self,name):
        """
         Extract a particular property.

         Parameters:
         ------------------
         name : str
             The name of the property available from the list obtained from the get_prop_list() function.

         Returns:
         ------------------
         value
             The value of the specified property.

         Notes:
         ------------------
         This function extracts the value of a specific property identified by its name from the data retrieved for the current Materials ID (mpid). The property name should be one of the properties listed in the output of the get_prop_list() function.
        """
        return self.data[name]
    def getkpt(self,primitive=True):
        """
        Compute k-points based on k-point density.
        Parameters:
        --------------------
        primitive: logical
                 Use primitive standard structure
        Returns:
        ---------------------
        kpt : list
            The k-point grid.
        kptype : str
            The type of k-point grid.
        kptshift : list
            The shifts in the k-point grid, returns [0,0,0].

        Notes:
        ---------------------
        This function computes the k-point grid based on a specified k-point density,
        defaulting to 0.05 if unspecified. It returns the k-point grid, type, and shifts.
        """
        if primitive:
            # FIX(5): was symprec=0.1 with the international monoclinic
            # setting, i.e. a different cell from the one setting_qeinput
            # writes; both use the module constants now.
            struc = SpacegroupAnalyzer(
                self.structure, symprec=SYMPREC
            ).get_primitive_standard_structure(
                international_monoclinic=INTERNATIONAL_MONOCLINIC)
        else:
            struc = self.structure
        relax_set = MPRelaxSet(structure=struc)
        relax_set.poscar.write_file('POSCAR')
        try:
            kptden = input_data['kptden']
        except KeyError:
            print("Default kptden of 0.05 being utilized\n")
            kptden = 0.05
        self.kpt = pos_to_kpt("POSCAR",kptden)
        self.kptype = "automatic"
        self.kpshift = [0, 0, 0]
        remove_files("POSCAR")
        print("*********************************\n")
        print("KPOINT with kpoint density of {} \n".format(kptden))
        print("*********************************\n")
        return self.kpt,self.kptype,self.kpshift

    def getevenkpt(self):
        """
        Make the k-point mesh even.

        Returns:
        ---------------
        evenkpt : tuple
            K-grid with an even number of points in all directions.

        Notes:
        ---------------
        Ensures the generated k-point grid has even points in each dimension
        by incrementing odd components to the next even number.
        Returns the resulting even k-point grid as a tuple.
        """
        kptsize = len(self.kpt)
        kpoint_list = self.kpt
        for i in range(kptsize):
            if kpoint_list[i]%2 == 0:
                kpoint_list[i] = kpoint_list[i]
            else:
                kpoint_list[i] = kpoint_list[i] + 1
        self.evenkpt = tuple(kpoint_list)
        return self.evenkpt

    def getecut_sssp(self,element):
        """
        Obtain the kinetic energy cutoff for a specific element.

        Parameters:
        ---------------------
        element : str
            The type of element for which the kinetic energy cutoff is required.

        Returns:
        ----------------------
        float
            Kinetic energy cutoff for the particular element.

        Raises:
        ----------------------
        KeyError
            When the element is in neither ``config['pseudo']['PSEUDO']`` nor
            the module's :data:`SSSP_EFFICIENCY` table.

        Notes:
        ----------------------
        Gets cutoff values for elements from config.json and falls back to the
        SSSP efficiency set for anything the configuration does not cover.
        """
        psd_data = input_data.get("pseudo", {}).get("PSEUDO", {})
        self.dict_element = dict(psd_data) if psd_data else dict(SSSP_EFFICIENCY)
        # spin-decorated species arrive as 'Fe,spin=5'
        symbol = bare_element(element)
        # FIX(2): the old ``except KeyError`` path re-raised, so one element
        # missing from pseudo.PSEUDO aborted the whole scan.  Fall back to the
        # packaged SSSP table first and only then raise, naming the element.
        try:
            return self.dict_element[symbol]
        except KeyError:
            pass
        if symbol in SSSP_EFFICIENCY:
            warnings.warn(
                "element {!r} is missing from pseudo.PSEUDO in config.json; "
                "using the SSSP efficiency value of {} Ry".format(
                    symbol, SSSP_EFFICIENCY[symbol]),
                RuntimeWarning, stacklevel=2)
            return SSSP_EFFICIENCY[symbol]
        raise KeyError(
            "no plane-wave cutoff known for element {!r} (from {!r}): it is in "
            "neither pseudo.PSEUDO of config.json nor the built-in SSSP "
            "efficiency table.  Add \"{}\": <ecutwfc in Ry> to the "
            "pseudo.PSEUDO dictionary of your config.json.".format(
                symbol, element, symbol))
    def maxecut_sssp(self):
        """
        Choose the maximum kinetic energy cutoff among elements in a compound.

        Returns:
        -----------------
        tuple
            A tuple containing the maximum kinetic energy cutoff for waveFunction (ecutwfc)
            and the maximum kinetic energy cutoff for charge density (ecutrho).

        Notes:
        -----------------
        This function finds the maximum kinetic energy cutoff for compound elements,
        then returs the cutoffs for waveFunction and density.
        """
        # Try to retrieve elements from the fetched data.  The bare ``except:``
        # that used to follow swallowed KeyboardInterrupt; ``getecut_sssp``
        # now strips any ',spin=' decoration itself, so the third branch is no
        # longer needed.
        try:
            self.comp_list = self.data['elements']
        except (KeyError, TypeError):
            # If 'elements' is not in the fetched data, take them from the structure
            self.comp_list = [str(el) for el in self.structure.elements]
        # Get the kinetic energy cutoff for each element in the compound
        cutoff_list = [self.getecut_sssp(el) for el in self.comp_list]
        # Determine the maximum kinetic energy cutoff for waveFunction
        self.ecutwfc= max(cutoff_list)
        # Determine the maximum kinetic energy cutoff for charge density
        self.ecutrho = 8 * self.ecutwfc
        print("******************************************\n")
        print("SSSP tested K.E. cutoffs for waveFunction and density\n")
        print("******************************************\n")
        return self.ecutwfc,self.ecutrho
    def maxecut_sssp_for_subs(self):
        """
        Determine the maximum kinetic energy cutoff among elements in a compound during the substitution process.

        Returns:
        -----------------
        tuple
            A tuple containing the maximum kinetic energy cutoff for waveFunction (ecutwfc)
            and the maximum kinetic energy cutoff for charge density (ecutrho).

        Notes:
        -----------------
        Similar to maxecut_sssp but for substitution process.
        """
        # The bare ``except:`` that used to guard this swallowed
        # KeyboardInterrupt and hid a genuinely unknown element;
        # ``getecut_sssp`` strips any ',spin=' decoration itself now.
        if not self.comp_list:
            self.comp_list = [str(el) for el in self.structure.elements]
        cutoff_list = [self.getecut_sssp(el) for el in self.comp_list]
        self.ecutwfc= max(cutoff_list)
        self.ecutrho = 8 * self.ecutwfc
        return self.ecutwfc,self.ecutrho

    def ecut_set(self,ecutwfc=50.0,ecutrho=400.0):
        """
        Function to set kinetic energy cutoffs
        """
        self.ecutrho=ecutrho
        self.ecutwfc=ecutwfc

    def get_properties(self,property_name):
        """
        Function to extract multiple properties from MP.
        parameters
        -------------------
        property_name : list of properties: Default: ['material_id']
        Returns
        -------------------
        property_list : list of properties extracted.
        """
        property_list = []
        for prop in property_name:
            property_list.append(self.data[prop])
        property_list.insert(0,self.data['material_id'])
        with open('mpid.csv', 'a') as data:
            for prop in property_list:
                data.write(str(prop) + ",")
            data.write("\n")
        return property_list
    @staticmethod
    def _resolve_magnetization(magdict, pseudo_species):
        """Give every *undecorated* species a usable starting magnetization.

        Parameters:
        - magdict (OrderedDict): label -> sign produced by
          :func:`convert_species_list` (0 for a species carrying no spin).
        - pseudo_species (list): the original species strings, so a species
          that really is spin decorated can be told from one that is not.

        Returns:
        OrderedDict mapping each label to the value written as
        ``starting_magnetization``.
        """
        # FIX(1c): an undecorated structure got starting_magnetization(i) = 0
        # for every species under nspin = 2 -- a non-magnetic calculation at
        # twice the cost.  Take the per-element value from
        # config['magmom']['magmom'] instead, and warn + use DEFAULT_MAGMOM
        # when the element is not listed there.
        undecorated = {item for item in pseudo_species if ',' not in item}
        magmom_cfg = input_data.get('magmom', {}).get('magmom', {}) or {}
        resolved = OrderedDict()
        for label, value in magdict.items():
            if label not in undecorated:
                resolved[label] = value
                continue
            element = bare_element(label)
            if element in magmom_cfg:
                value = float(magmom_cfg[element])
            else:
                warnings.warn(
                    "no spin decoration and no magmom.magmom entry for element "
                    "{!r}; using starting_magnetization = {} rather than 0, "
                    "which would make this nspin = 2 run non-magnetic".format(
                        element, DEFAULT_MAGMOM),
                    RuntimeWarning, stacklevel=3)
                value = DEFAULT_MAGMOM
            if abs(value) > 1.0:
                warnings.warn(
                    "magmom.magmom[{!r}] = {} is a moment in Bohr magnetons, but "
                    "QE's starting_magnetization is a fraction in [-1, 1]; "
                    "clamping to {}".format(element, value,
                                            1.0 if value > 0 else -1.0),
                    RuntimeWarning, stacklevel=3)
                value = 1.0 if value > 0 else -1.0
            resolved[label] = value
        return resolved

    def _write_magnetic_input(self, scf_name, site_labels, site_coords, magdict):
        """Turn ``scf_name`` into a collinear spin-polarised QE input.

        Parameters:
        - scf_name (str): the QE input file to rewrite, in place.
        - site_labels (list): per-site type label ('Fe1', 'Fe2', 'Pd', ...),
          in the site order of the structure that was written.
        - site_coords (list): the matching fractional coordinates.
        - magdict (dict): label -> starting_magnetization.

        The ATOMIC_SPECIES block, ``ntyp``, ``nat``, ``nspin`` and the
        ``starting_magnetization(i)`` entries are rebuilt from those lists.
        """
        lines = read_lines(scf_name)
        try:
            spc_start = next(i for i, ln in enumerate(lines) if 'ATOMIC_SPECIES' in ln)
            pos_start = next(i for i, ln in enumerate(lines) if 'ATOMIC_POSITIONS' in ln)
        except StopIteration as exc:
            raise RuntimeError(
                "{} has no ATOMIC_SPECIES/ATOMIC_POSITIONS block to "
                "decorate".format(scf_name)) from exc
        natoms = len(site_labels)
        if pos_start + natoms >= len(lines):
            raise RuntimeError(
                "{} has fewer than the {} atomic positions the structure "
                "carries".format(scf_name, natoms))
        # distinct labels, in order of first appearance -- this is the QE type
        # index that starting_magnetization(i) refers to
        type_labels = list(OrderedDict.fromkeys(site_labels))
        # FIX(1b): one ATOMIC_SPECIES line per *label*, with the pseudopotential
        # named after the bare element.  The old code renamed line k of the
        # block for every spin-decorated species, which only lines up when
        # pymatgen happens to emit one line per decorated species, and left
        # ntyp describing the undecorated cell.
        species_block = ["ATOMIC_SPECIES"]
        for label in type_labels:
            element = bare_element(label)
            species_block.append("  {}  {:.4f} {}.upf".format(
                label, Element(element).atomic_mass, element))
        # FIX(1a): the labels and coordinates come from the standardised cell
        position_block = [lines[pos_start]]
        for j, label in enumerate(site_labels):
            coords = site_coords[j]
            position_block.append("  {} {:.10f} {:.10f} {:.10f}".format(
                label, coords[0], coords[1], coords[2]))
        # end of the existing ATOMIC_SPECIES block
        spc_end = spc_start + 1
        while (spc_end < len(lines) and lines[spc_end].strip()
               and not CARD_RE.match(lines[spc_end])):
            spc_end += 1
        new_lines = (lines[:spc_start] + species_block
                     + lines[spc_end:pos_start] + position_block
                     + lines[pos_start + 1 + natoms:])
        # drop whatever &SYSTEM already said about spin and type counts
        new_lines = [ln for ln in new_lines
                     if not SYSTEM_KEY_RE.match(ln)]
        try:
            sys_at = next(i for i, ln in enumerate(new_lines)
                          if '&SYSTEM' in ln.upper())
        except StopIteration as exc:
            raise RuntimeError(
                "{} has no &SYSTEM namelist".format(scf_name)) from exc
        additions = ["  nspin = 2,",
                     "  ntyp = {},".format(len(type_labels)),
                     "  nat = {},".format(natoms)]
        for index, label in enumerate(type_labels):
            additions.append("  starting_magnetization({}) = {},".format(
                index + 1, magdict.get(label, DEFAULT_MAGMOM)))
        new_lines = new_lines[:sys_at + 1] + additions + new_lines[sys_at + 1:]
        write_lines(scf_name, new_lines)

    def setting_qeinput(self,calculation='vc-relax',occupations='smearing',restart_mode='from_scratch',pseudo_dir='./',smearing=0.02,smearing_type='gauss',etot_conv_thr=1e-05,forc_conv_thr=1e-04,conv_thr=1e-16,ion_dynamics='bfgs',cell_dynamics='bfgs',magnetic=False,primitive=True):
        """
        Function to create input file for QE ground-state calculations.

        Parameters:
        - calculation (str): Type of calculation. Default: 'vc-relax'. Other options are 'relax' (ionic only), 'bands' for band structure, 'scf' for SCF calculations.
        - occupations (str): Occupation. Default: 'smearing'. Other options could be 'tetrahedra' and so on.
        - restart_mode (str): How to start the calculations.
        - pseudo_dir (str): Path to pseudopotential files. Default: './' (current directory).
        - smearing (float): Degauss value. Default: 0.02.
        - smearing_type (str): Type of smearing. Default: 'gauss'.
        - etot_conv_thr, forc_conv_thr, conv_thr (float): Convergence parameters of QE calculations.
        - ion_dynamics, cell_dynamics (str): Algorithm to perform relaxation. Default: 'bfgs'.
        - magnetic (logical): Magnetic flag if magnetic input need to be created
        - primitive (logical): if True, then use primitive structure

        Returns:
        Creates input files in scf-mpid.in format inside scf_dir/.
        """
        # lW  2010_SC has monoclinic conventional cell with alpha<90, beta=gamma=90,
        # different from international beta=90
        # FIX(1a): this used to run *after* the magnetic branch had read site
        # labels and fractional coordinates off self.structure, so the labels
        # and coordinates of the original cell were pasted onto the sites of
        # the standardised one.  Standardise first; everything below then
        # describes the cell that is actually written.
        if primitive:
            tmpanalyzer = SpacegroupAnalyzer(self.structure, symprec=SYMPREC)
            self.structure = tmpanalyzer.get_primitive_standard_structure(
                international_monoclinic=INTERNATIONAL_MONOCLINIC)
        # If magnetic calculation is enabled, handle pseudopotentials and magnetizations
        magdict = OrderedDict()
        site_labels = []
        site_coords = []
        if magnetic:
            elements_spin = self.structure.elements
            pseudo1 = {}
            for site in self.structure:
                site_labels.append(site.label)
                site_coords.append(site.frac_coords)
            site_labels, _ = convert_species_list(site_labels)
            pseudo_species = [str(elm_spin) for elm_spin in elements_spin]
            _, magdict = convert_species_list(pseudo_species)
            magdict = self._resolve_magnetization(magdict, pseudo_species)
            # FIX(1b)/FIX(33): key the dict exactly the way pymatgen looks it
            # up.  PWInput.__init__ does
            #     for species in self.structure.composition:
            #         if str(species) not in pseudo: raise PWInputError(...)
            # so for a spin-decorated structure the key is the *decorated*
            # string, 'Mg,spin=1.0'.  FIX(1b) had changed this to the bare
            # symbol on the belief that PWInput used site.specie.symbol; it
            # does not, and every magnetic QE input died with
            # "Missing Mg,spin=1.0 in pseudo specification!".  The file name
            # still uses the bare element, which is what the .upf is called.
            for species in self.structure.composition:
                element = bare_element(getattr(species, "symbol", species))
                pseudo1[str(species)] = element + ".upf"
                # tolerate an undecorated lookup too
                pseudo1.setdefault(element, element + ".upf")
        else:
            pseudo1 = {el:el+'.upf' for el in self.comp_list}
        # Set prefix for output files
        prefix = self.prefix
        # Set up control, system, and electrons parameters
        # FIX(6): these used to be the module-global config sub-dictionaries,
        # mutated in place below, so ecutwfc/ecutrho/prefix of one material
        # leaked into every later material handled by the same process.
        pwscf_in = input_data.get("pwscf_in") or {}
        if pwscf_in.get('control') and pwscf_in.get('system') and pwscf_in.get('electrons'):
            control = copy.deepcopy(pwscf_in['control'])
            system = copy.deepcopy(pwscf_in['system'])
            electrons = copy.deepcopy(pwscf_in['electrons'])
        else:
            control = {'calculation':calculation, 'nstep':300, 'restart_mode':restart_mode, 'pseudo_dir':pseudo_dir, 'outdir':'./', 'tprnfor':'.true.','tstress':'.true.', 'etot_conv_thr':etot_conv_thr, 'forc_conv_thr':forc_conv_thr}
            system = {'smearing':smearing_type, 'occupations':occupations, 'degauss':smearing}
            electrons = {'diagonalization':'david', 'mixing_mode':'plain', 'mixing_beta':0.7, 'conv_thr': conv_thr, 'electron_maxstep':300}
        # Set kinetic energy cutoffs
        system['ecutwfc'] = self.ecutwfc
        system['ecutrho'] = self.ecutrho
        control['prefix'] = prefix
        # Generate the input file based on the chosen calculation type
        if calculation == 'vc-relax':
            filename = pwscf.PWInput(self.structure, pseudo=pseudo1, control=control, system=system,electrons=electrons, kpoints_grid=self.kpt, ions={'ion_dynamics':ion_dynamics}, cell={'cell_dynamics':cell_dynamics,'press_conv_thr':0.05}, format_options={"coord_decimals": 10})
        elif calculation == 'relax':
            filename = pwscf.PWInput(self.structure, pseudo=pseudo1, control=control, system=system,electrons=electrons, kpoints_grid=self.kpt, ions={'ion_dynamics':ion_dynamics}, format_options={"coord_decimals": 10})
        else:
            filename = pwscf.PWInput(self.structure, pseudo=pseudo1, control=control, system=system,electrons=electrons, kpoints_grid=self.kpt, format_options={"coord_decimals": 10})
        # Write input file
        filename.write_file("temp.in")
        scf_name = "scf-{}.in".format(self.mpid)
        kpoint_name = "kpoint-{}.dat".format(self.mpid)
        kpath_name = "kpath-{}.dat".format(self.mpid)
        # Adjust the input file format and save it with the appropriate filename.
        # (was: sed "s/'.true.'/.true./")
        lines = [ln.replace("'.true.'", ".true.").replace("'.false.'", ".false.")
                 for ln in read_lines("temp.in")]
        write_lines(scf_name, lines)
        write_lines(kpoint_name,
                    section(lines, "K_POINTS automatic", "CELL_PARAMETERS angstrom"))
        obj = INPUTscf(scf_name)
        # Use symmetric structure
        obj.standardize(self.mpid, output="temp.dat")
        header = section(lines, "&CONTROL", "ATOMIC_SPECIES")
        species = section(lines, "ATOMIC_SPECIES", "ATOMIC_POSITIONS crystal")
        remove_files("temp.in")
        body = read_lines("temp.dat")
        if calculation == 'bands':
            # drop the automatic K_POINTS mesh and its following line
            trimmed = []
            skip = 0
            for line in body:
                if skip:
                    skip -= 1
                    continue
                if 'K_POINTS' in line:
                    skip = 1
                    continue
                trimmed.append(line)
            obj.generate_kpath(nqpoint=200, kcut=0, out=kpath_name)
            body = read_lines(kpath_name) + trimmed
            print("Input for band calculation. Provide nband = <n> within SYSTEM section. Use sufficient <n>.")
        write_lines(scf_name, header + species + body)
        # If magnetic calculation is enabled, modify the input file accordingly
        # to incorporate magnetic keywords in the input files
        if magnetic:
            self._write_magnetic_input(scf_name, site_labels, site_coords, magdict)
        # Clean up temporary files
        # FIX(9): the glob 'kpoint*' also deleted the kpoint files of every
        # other material sitting in the same working directory; only the files
        # this call created are removed now.
        remove_files("temp.dat", kpoint_name, kpath_name)
def reorder_dictionary(original_dict, new_order):
    """
    Reorder the keys of a dictionary based on a new list of keywords and return the list of keys.

    Parameters:
    - original_dict (dict): The original dictionary whose keys are to be reordered.
    - new_order (list): The new list of keywords specifying the desired order of keys.

    Returns:
    - list: The list of keys of the reordered dictionary based on the new order.

    Example:
    >>> original_dict = {'a': 1, 'b': 2, 'c': 3, 'd': 4}
    >>> new_order = ['c', 'a', 'd', 'b']
    >>> reorder_dictionary(original_dict, new_order)
    """
    # Create a new OrderedDict with the keys ordered according to the new list
    ordered_dict = OrderedDict((key, original_dict[key]) for key in new_order if key in original_dict)

    # Return the list of keys in the new order
    return ordered_dict


def extract_elements_with_spin(original_list):
    """
    Extract elements with their spin values from the given list.

    Parameters:
    - original_list (list): List containing strings representing elements with spin values.

    Returns:
    - dict: A dictionary where keys are elements and values are lists of spin values.
    """
    elements_with_spin = {}
    for item in original_list:
        if ',' in item:
            element, spin = item.split(',')
            spin_value = spin.split('=')[1]
            if element not in elements_with_spin:
                elements_with_spin[element] = [spin_value]
            elif spin_value not in elements_with_spin[element]:
                elements_with_spin[element].append(spin_value)
    return elements_with_spin

def parse_spin(label):
    """The spin out of a species label such as ``'Fe,spin=5'``.

    FIX(37): this was ``float(item.split('=')[-1])``, which numpy 2 broke.
    ``Species.__str__`` renders the spin with ``repr``, and numpy 2 changed the
    repr of a scalar from ``5.0`` to ``np.float64(5.0)`` -- so a structure
    decorated from ``config['magmom']['magmom']`` (whose values arrive as numpy
    floats through pymatgen) produced ``'Fe,spin=np.float64(5.0)'`` and
    ``mainprogram magenum`` died with

        ValueError: could not convert string to float: 'np.float64(5.0)'

    Pull the first number out of whatever the label says instead of assuming
    the text after ``=`` is already a bare float.  Returns 0.0 when the label
    carries no number at all.
    """
    text = str(label).split("=", 1)[-1].strip()
    # 'np.float64(5.0)' -> '5.0'.  Taking the first number in the raw text
    # would take the 64 out of "float64", so unwrap the call first.
    wrapped = re.search(r"\(([^()]*)\)\s*$", text)
    if wrapped:
        text = wrapped.group(1)
    match = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)
    return float(match.group(0)) if match else 0.0


def convert_species_list(original_list):
    """
    Rename elements in the original list with elements with different spins.

    Parameters:
    - original_list (list): List containing strings representing elements.

    Returns:
    - list: Updated list where elements are replaced with generic names for different spins.

    Example:
    >>> original_list = ['Fe,spin=5', 'Fe,spin=-5', 'pd', 'pd', 'I,spin=1', 'I,spin=-1']
    >>> convert_species_list(original_list)
    ['Fe1', 'Fe2', 'pd', 'pd', 'I1', 'I2']
    """
    elements_spin = extract_elements_with_spin(original_list)
    # Updated list
    updated_list = []
    magdict = {}
    # Iterate through the original list
    for item in original_list:
        if ',' in item:
            element = item.split(',')[0]
            spin = parse_spin(item)
            if len(elements_spin[element]) > 1 and spin > 0:
                updated_list.append(element + str(1))
                magdict[element + str(1)] = 1
            elif len(elements_spin[element]) > 1 and spin < 0:
                updated_list.append(element + str(2))
                magdict[element + str(2)] = -1
            else:
                updated_list.append(element)
                spin_s = spin
                mag_s = 1 if spin_s > 0 else -1 if spin_s < 0 else 0
                magdict[element] = mag_s
        else:
            updated_list.append(item)
            magdict[item] = 0
    return updated_list,magdict


def ase_cell_to_structure(ase_cell):
    # Extract lattice vectors from ASE Atoms object
    lattice_vectors = ase_cell.cell.tolist()
    
    # Extract atomic positions and species from ASE Atoms object
    atomic_positions = ase_cell.get_scaled_positions()
    species_symbols = ase_cell.get_chemical_symbols()
    
    # Create pymatgen Lattice object
    lattice = Lattice(lattice_vectors)
    
    # Create list of pymatgen Specie objects
    species = [Element(symbol) for symbol in species_symbols]
    
    # Create pymatgen Structure object
    structure = Structure(lattice, species, atomic_positions)
    
    return structure

class INPUTscf:
    """
    class to process QE input files
    parameters
    ------------------
    filename: (str) QE input file
    """
    def __init__(self,filename='scf.in'):
        self.filename = filename
        self.ase_cell = espresso.read_espresso_in(self.filename)
        # Convert the ASE cell.io.espresso object to a pymatgen Structure object
        structure = ase_cell_to_structure(self.ase_cell)
        # Use SpacegroupAnalyzer to get the symmetric structure
        symmetry_analyzer = SpacegroupAnalyzer(structure)
        symmetric_structure = symmetry_analyzer.get_symmetrized_structure()
        # Now, symmetric_structure contains the symmetrized version of the structure
        self.struc = symmetric_structure
        self.cell = symmetric_structure.lattice.matrix
        self.volume = symmetric_structure.lattice.volume
        #self.braiv_latt = Cell.get_bravais_lattice(self.cell)
        self.mpid = None
        self.comp = None
        self.prefix = None
        self.qpoint = []
        self.mass = []

    def scftocif(self,output='file.cif'):
        """
        Function to convert QE input file to structure file in .cif format
        parameters
        ---------------
        output : (str) name of the output file
        Returns
        ---------------
        file2 : output file object
        """
        self.struc.to(output)
        return self.struc

    def cellpar(self):
        """
        Function to calculate cell lengths and angles
        Returns
        -----------------------
        list of lenghts and angles
        """
        length = self.struc.lattice.abc
        angles = self.struc.lattice.angles
        return list(length) + list(angles)

    def standardize(self,mpid,output="standard.in"):
        """
        Function to get symmetrized structure.
        parameters
        -----------------
        mpid : (str) materials project ID
        output : (str) output file. Default: 'standard.in'
        """
        finalpos = self.struc.frac_coords
        finalcell = self.cell
        specieslist = self.struc.species
        with open("kpoint-{}.dat".format(mpid), "r") as kpoint:
            kplines = kpoint.readlines()
        with open(output, "w") as struc:
            struc.write("ATOMIC_POSITIONS crystal\n")
            for i,_ in enumerate(specieslist):
                struc.write(str(specieslist[i]) + " " + str(finalpos[i][0])+ " ")
                struc.write(str(finalpos[i][1]) + " " + str(finalpos[i][2]) + "\n")
            for i in range(2):
                struc.write(kplines[i])
            struc.write("CELL_PARAMETERS angstrom\n")
            for i in range(3):
                struc.write(str(finalcell[i][0]) + " " + str(finalcell[i][1]) + " " + str(finalcell[i][2]) + "\n")

    def generate_kpath(self,nqpoint=200,kcut=0,out='kpath.in',qualifier='crystal'):
        """
        Function to generate and write k-point mesh for bandstructure calculation to a file
        parameters
        --------------------
        nqpoint : (int) size of the k-point mesh. Default: 200
        kcut : (int) cutoff to the high-symmetry path of the Brillouin zone. Default: 0 for full Brillouin zone
        out : (str) output file. Default: 'kpath.in'
        qualifier : (str) K_POINTS unit. Default: 'crystal', matching
            :func:`htesp.kpath.printk`.  Use 'crystal_b' for the segment form.
        Returns
        --------------------
        kpts : numpy array of kpoints in linear axis after processing
        n : (int) size of kpts
        """
        kpts,_,_,_,_,_ = kpath(self.filename,nqpoint,kcut)
        nkpt = kpts.shape[0]
        with open(out, 'w') as kmesh:
            # FIX(3): a bare 'K_POINTS' makes QE read the list as tpiba (2pi/a
            # cartesian) while ASE's bandpath().kpts are fractional, so the
            # band structure was computed along the wrong path.  htesp.kpath
            # .printk writes 'K_POINTS crystal'; the two agree now.
            kmesh.write('K_POINTS {}\n'.format(qualifier))
            kmesh.write(str(nkpt) + '\n')
            for i in range(nkpt):
                kmesh.write(str(round(kpts[i][0],8)) + " " + str(round(kpts[i][1],8)) + " " + str(round(kpts[i][2],8)) + " " + str(0.0) + "\n")
        return nkpt,kpts

    def setting_input(self,comp,mass,qpoint):
        """
        Function to setup input
        parameters
        -------------
        comp : compound name
        mass: List of masses of elements
        qpoint : List of qpoint
        """
        self.comp=comp
        self.mass=mass
        self.qpoint=qpoint
        obj=MpConnect()
        obj.setting(self.comp)
        self.mpid = obj.mpid
        self.prefix = obj.prefix

    def create_elph(self,output='elph.in',tol=0.00000000000001,sigma=0.005):
        """
         Function similar to elph.py file inside src/
        """
        dynmat = self.prefix.replace("'", "") + ".dyn"
        nat = len(self.mass)
        with open(output, 'w') as elph:
            elph.write("electron phonon coupling \n")
            elph.write("&inputph" + "\n")
            elph.write("tr2_ph={},".format(tol) + "\n")
            elph.write("prefix={},".format(self.prefix) + "\n")
            elph.write("fildvscf='aldv'," + "\n")
            for i in range(1,nat+1):
                elph.write("amass({})={},".format(i,self.mass[i-1]) + "\n")
            elph.write("outdir='./'," + "\n")
            elph.write("fildyn='{}',".format(dynmat) + "\n")
            elph.write("electron_phonon='interpolated'," + "\n")
            elph.write("el_ph_sigma={},".format(sigma) + "\n")
            elph.write("el_ph_nsigma=10," + "\n")
            elph.write("trans=.true.," + "\n")
            elph.write("ldisp=.true." + "\n")
            elph.write("nq1={},nq2={},nq3={}".format(int(self.qpoint[0]),int(self.qpoint[1]),int(self.qpoint[2])) + "\n")
            elph.write("/" + "\n")

    def create_q2r(self,output='q2r.in'):
        """
         Function similar to q2r.py file inside src/
        """
        dynmat = self.prefix.replace("'", "") + ".dyn"
        frc = self.prefix.replace("'", "") + ".fc"
        with open(output, 'w') as q2r_write:
            q2r_write.write("&input" + "\n")
            q2r_write.write("zasr='simple'," + "\n")
            q2r_write.write("fildyn='{}',".format(dynmat) + "\n")
            q2r_write.write("flfrc='{}',".format(frc) + "\n")
            q2r_write.write("la2F=.true." + "\n")
            q2r_write.write("/" + "\n")

    def create_matdyn(self,out='matdyn.in',nqpt=60,kcut=0):
        """
         Function similar to matdyn.py file inside src/
        """
        freq = self.prefix.replace("'", "") + ".freq"
        frc = self.prefix.replace("'", "") + ".fc"
        eig = self.prefix.replace("'", "") + ".eig"
        nat = len(self.mass)
        nkpt,kpts = self.generate_kpath(nqpoint=nqpt,kcut=kcut)
        with open(out, 'w') as matdyn_write:
            matdyn_write.write("&input" + "\n")
            matdyn_write.write("asr='simple'," + "\n")
            for i in range(1,nat+1):
                matdyn_write.write("amass({})={},".format(i,self.mass[i-1]) + "\n")
            matdyn_write.write("flfrc='{}',".format(frc) + "\n")
            matdyn_write.write("flfrq='{}',".format(freq) + "\n")
            matdyn_write.write("fleig='{}',".format(eig) + "\n")
            matdyn_write.write("la2F=.true.," + "\n")
            matdyn_write.write("dos=.false.," + "\n")
            # FIX(4): the q-points written below come from
            # generate_kpath/ASE and are fractional, so matdyn.x must be told
            # to read them as crystal coordinates -- htesp/matdyn.py:41 already
            # does this and the two files must agree.
            matdyn_write.write("q_in_cryst_coord=.true." + "\n")
            matdyn_write.write("/" + "\n")
            matdyn_write.write(str(nkpt) + "\n")
            for i in range(nkpt):
                matdyn_write.write(str(round(kpts[i][0],8)) + " " + str(round(kpts[i][1],8)) + " ")
                matdyn_write.write(str(round(kpts[i][2],8)) + " " + str(0.0) + "\n")
    def create_phdos(self,kpts,out='phdos.in',ndos=200):
        """
         Function similar to matdyn_dos.py file inside src/
        """
        freq = self.prefix.replace("'", "") + "-dos.freq"
        frc = self.prefix.replace("'", "") + ".fc"
        nat = len(self.mass)
        with open(out, 'w') as phdos:
            phdos.write("&input" + "\n")
            phdos.write("asr='simple'," + "\n")
            for i in range(1,nat+1):
                phdos.write("amass({})={},".format(i,self.mass[i-1]) + "\n")
            phdos.write("flfrc='{}',".format(frc) + "\n")
            phdos.write("flfrq='{}',".format(freq) + "\n")
            phdos.write("la2F=.true.," + "\n")
            phdos.write("dos=.true.," + "\n")
            phdos.write("fldos='phonon.dos'," + "\n")
            phdos.write("nk1={},nk2={},nk3={},ndos={},".format(kpts[0],kpts[1],kpts[2],ndos) + "\n")
            phdos.write("/" + "\n")

    def create_dos(self,out1='dos.in',out2='pdos.in'):
        """
         Function similar to dos.py file inside src/
        """
        dynmat = self.prefix.replace("'", "") + ".dos"
        dynmat1 = self.prefix.replace("'", "") + ".pdos"
        with open(out1, 'w') as dos:
            dos.write("&dos" + "\n")
            dos.write("prefix={},".format(self.prefix) + "\n")
            dos.write("outdir='./'," + "\n")
            dos.write("fildos='{}',".format(dynmat) + "\n")
            dos.write("DeltaE=0.01" + "\n")
            dos.write("/" + "\n")
        with open(out2, 'w') as pdos:
            pdos.write("&projwfc" + "\n")
            pdos.write("prefix={},".format(self.prefix) + "\n")
            pdos.write("outdir='./'," + "\n")
            pdos.write("pfildos='{}',".format(dynmat1) + "\n")
            pdos.write("DeltaE=0.01" + "\n")
            pdos.write("/" + "\n")

    def post_band(self,out='band.in'):
        """
         Function similar to band.py file inside src/
        """
        dynmat = self.prefix.replace("'", "") + ".dat"
        with open(out, 'w') as band:
            band.write("&BANDS" + "\n")
            band.write("prefix={},".format(self.prefix) + "\n")
            band.write("outdir='./'," + "\n")
            band.write("filband='{}',".format(dynmat) + "\n")
            band.write("lsym=.true." + "\n")
            band.write("/" + "\n")

    def post_phband(self,out='phonband.in'):
        """
         Function similar to phonband.py file inside src/
        """
        freq = self.prefix.replace("'", "") + ".freq"
        with open(out, 'w') as phband:
            phband.write(freq + "\n")
            phband.write("0 5000" + "\n")
            phband.write("freq.plot" + "\n")
            phband.write("freq.ps" + "\n")
            phband.write("0.0" + "\n")
            phband.write("100.0 0.0" + "\n")

def scf_to_dos_scf(input_file='scf.in',out='scf-dos.in'):
    """
    Function to change QE input file from using smearing to 'tetrahedra' method, mainly for DOS calculation

    parameters
    ----------------------
    input_file : (str) input file for QE scf calculation
    out : (str) output file after modification

    """
    lines = [ln for ln in read_lines(input_file) if 'calculation' not in ln]
    write_lines(input_file, lines)
    insert(input_file=input_file,output='temp',keyword='&CONTROL',what="calculation = 'nscf',")
    lines = [ln.replace("'smearing'", "'tetrahedra'") for ln in read_lines('temp')]
    lines = [ln for ln in lines if 'degaus' not in ln and 'smearing' not in ln]
    write_lines(out, lines)
    remove_files('temp')
    print("Use denser k-mesh for dos calculations")

def insert(input_file='scf.in',output='scf-new.in',keyword=None,where="after",what=""):
    """
    Function to insert keyword in input file
    parameters
    -----------------------
    input_file : (str) input file. Default: 'scf.in'
    output : (str) output file. Default: 'scf-new.in'
    keyword : (str) keyword to look after
    where : (str) where to insert. Default: 'after', otherwise 'before'
    what : (str) what to insert. Default: ''
    """
    lines = read_lines(input_file)
    out_lines = []
    for line in lines:
        if keyword is not None and keyword in line:
            if where == "after":
                out_lines.append(line)
                out_lines.append("  " + what)
            else:
                out_lines.append("  " + what)
                out_lines.append(line)
        else:
            out_lines.append(line)
    write_lines(output, out_lines)
# Extract relaxed structure and update QE input.
class OUTPUTscf:
    """
    class to extract relax structure from QE output file and update input file
    parameters
    ---------------------
    filename : (str) QE output file. Default: 'scf.out'
    """
    def __init__(self,filename='scf.out'):
        self.filename = filename

    def extract_relax(self,output='relax.dat'):
        """
        Function to extract cell and positions of crystal structures from scf output file and write to a file
        parameters
        ----------------
        output : (str) output file. Default: 'relax.dat'
        Returns
        ----------------
        list of the lines written (empty when the run produced no coordinates)
        """
        # FIX(8): the old sed pipeline ('$d' then '1,4d' then '5d') hard-coded
        # the vc-relax layout -- the volume/density lines followed by
        # CELL_PARAMETERS.  For calculation='relax' there is no cell block, so
        # it deleted the header *and the first two atoms*, and for an
        # unconverged run (no final-coordinates block at all) it silently wrote
        # an empty file.  Both are parsed by keyword now; see
        # htesp/workflow.py QEText.final_coordinates for the same semantics.
        lines = read_lines(self.filename)
        block = []
        begin = next((i for i, ln in enumerate(lines)
                      if 'Begin final coordinates' in ln), None)
        if begin is not None:
            end = next((i for i in range(begin + 1, len(lines))
                        if 'End final coordinates' in lines[i]), len(lines))
            keep = False
            for line in lines[begin + 1:end]:
                head = line.strip().split()[0] if line.strip() else ""
                if head in ("CELL_PARAMETERS", "ATOMIC_POSITIONS"):
                    keep = True
                    line = (line.replace("(angstrom)", "angstrom")
                                .replace("(crystal)", "crystal")
                                .replace("(alat=", "alat="))
                if keep:
                    block.append(line)
        else:
            # unconverged or still running: fall back to the last geometry
            starts = [i for i, ln in enumerate(lines) if 'ATOMIC_POSITIONS' in ln]
            if starts:
                last = starts[-1]
                block.append(lines[last].replace("(crystal)", "crystal")
                                        .replace("(angstrom)", "angstrom")
                                        .replace("(alat=", "alat="))
                for line in lines[last + 1:]:
                    if not line.strip() or CARD_RE.match(line):
                        break
                    if len(line.split()) < 4:
                        break
                    block.append(line)
        while block and not block[-1].strip():
            block.pop()
        if not block:
            warnings.warn(
                "{} contains neither a 'Begin final coordinates' block nor an "
                "ATOMIC_POSITIONS block: the run did not converge (or did not "
                "start).  {} is left empty.".format(self.filename, output),
                RuntimeWarning, stacklevel=2)
        write_lines(output, block)
        return block
    def update_scf(self,input_file='scf.in',output='scf-new.in',what='structure'):
        """
        Function to update QE input file and update with new structure
        parameters
        ------------
        input_file : (str) QE scf input file to update. Default: 'scf.in'
        output : (str) QE scf output file, updated with new structure. Default: 'scf-new.in'
        what : (str) what to update. Default: 'structure'. Other, not implemented yet !
        """
        lines = read_lines(input_file)
        kpoint_block = section(lines, "K_POINTS automatic", "CELL_PARAMETERS angstrom")
        header = section(lines, "&CONTROL", "ATOMIC_POSITIONS crystal")
        if what == 'structure':
            relax = self.extract_relax('relax.in')
            write_lines(output, header + kpoint_block + relax)
        else:
            print('do nothing\n')
        remove_files("relax.in")

#def extract_energy(self,input='scf.out'):
#    os.system("""Rytoev=`echo "scale=6;13.605698" | bc`""")
#    os.system("""en=`grep "!    total energy              =   " {} | tail -n1 | awk '{print $5}'` """.format(input))
#    os.system("""e=`echo "scale=6; $en * $Rytoev  " | bc`""")
#    os.system("""echo "the total energy: $e" """)
#    return
#class for calculations. In progress........
#class calculation:
#    def __init__(self,batch_header='batch.header'):
#        self.batch_header = batch_header
#def create_submission_scripts(self):
#def relax_structure(self,input='scf.in',out='scf.out'):
#def calculate_Tc(self,mu=0.16,degaussq=0.12,ngaussq=0,phfile='elph.out',phsigma=0.005):
