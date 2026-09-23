#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to extract data from OQMD database"""
import sys
import os
import re
import shutil
import subprocess
import time
from ase.io.vasp import read_vasp
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.core import structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import contextlib
import socket

import qmpy_rester as qr
from htesp.cif_to_gsinput import pos_to_kpt
from htesp.write_potcar import poscar2potcar
from htesp.htepc import MpConnect
from htesp.check_json import config
from htesp.inputin import InputIn

#: seconds a single OQMD socket operation may block before it is abandoned.
#:
#: FIX: ``qmpy_rester`` builds a bare ``requests.Session()`` and passes no
#: timeout to it, so a connection that stalls blocks for ever.  It did:
#: ``mainprogram oqmd-download`` was found nineteen minutes in, alive, with
#: two open sockets and not a line of output, and nothing short of a kill
#: would have ended it.  A socket timeout is the right instrument here
#: because it distinguishes the two cases a wall-clock budget cannot -- a
#: *stalled* connection is abandoned, while a query that is merely slow (OQMD
#: searches have taken anywhere from 35 to 100 seconds) keeps going as long
#: as data is still arriving.
OQMD_SOCKET_TIMEOUT = 30.0


@contextlib.contextmanager
def socket_timeout(seconds: float = OQMD_SOCKET_TIMEOUT):
    """Apply a default socket timeout for the duration of the block.

    Scoped rather than set at import, deliberately.
    ``socket.setdefaulttimeout`` is process-global and
    ``htesp/aflow_extract.py`` imports this module, so setting it on import
    would quietly put a timeout on AFLOW's queries too -- and on everything
    ``mainprogram data-combine`` does, which runs all three front ends in one
    process.  Only sockets opened inside the block are affected, and the
    previous value is restored even when the query raises.
    """
    previous = socket.getdefaulttimeout()
    socket.setdefaulttimeout(seconds)
    try:
        yield
    finally:
        socket.setdefaulttimeout(previous)


#: how many times search() retries an OQMD query with a smaller limit
MAX_SEARCH_RETRIES = 8
#: base delay, in seconds, between those retries (multiplied by the attempt)
RETRY_BACKOFF_SECONDS = 2.0

def poscar_to_input(calc_type,mpid,compound,keven):
    """
    Function to convert POSCAR into ground-state input files.

    Parameters
    ----------
    calc_type : str
        Calculation type. "QE" or "VASP".
    mpid : str
        Material id.
    compound : str
        Compound name.
    keven : bool
        If even k-mesh to use.

    Returns
    -------
    str
        Compound name.

    Notes
    -----
    This function reads a POSCAR file and generates input files
    required for ground-state calculations
    using either Quantum Espresso (QE) or VASP software.

    It checks if 'config.json' exists, and if so, it retrieves settings such as k-point density and
    download configurations.

    Based on the calculation type, it creates input files and directories accordingly. For VASP,
    it prepares INCAR, POSCAR, and POTCAR files and
    potentially runs VASP processing scripts if available.
    For QE, it creates input files with appropriate settings.

    Finally, it moves input files to the designated directories.

    """
    input_data = config()
    # check_json.config() always returns a fully populated dict now, so the
    # "config.json exists?" guard that used to leave ``d`` unbound is gone.
    d = input_data['download']
    evenkpt = d['inp']['evenkpt']
    kptden = input_data['kptden']
    # Read POSCAR with ASE
    data = read_vasp("POSCAR")
    symbol = list(data.symbols)
    relax_dir = os.path.join("R{}-{}".format(mpid, compound), "relax")
    # Creates input for VASP
    if calc_type in ('VASP','vasp'):
        # FIX(all): os.makedirs(..., exist_ok=True) instead of isdir-then-mkdir,
        # and no shell interpolation of mpid/compound.
        os.makedirs(relax_dir, exist_ok=True)
        # Write VASP input files
        structure_file = structure.Structure.from_file("POSCAR")
        structure_file = SpacegroupAnalyzer(structure_file, symprec=0.1).get_primitive_standard_structure()
        relax_set = MPRelaxSet(structure=structure_file)
        relax_set.poscar.write_file("POSCAR")
        # FIX(8): the k-mesh used to be computed from the raw cell, three
        # statements before MPRelaxSet replaced POSCAR with the standardised
        # primitive cell -- so KPOINTS described a different cell from the
        # POSCAR it shipped with.  It is computed here, from the POSCAR that
        # is actually written.
        if evenkpt:
            print("Even kpoint mesh is utilized\n")
            k_mesh = pos_to_kpt("POSCAR",kptden,True)
        else:
            k_mesh = pos_to_kpt("POSCAR",kptden)
        poscar2potcar()
        relax_set.incar.write_file("INCAR")
        # Copy files to R{mpid}-{compound}/relax/ folder
        # FIX(all): os.system("mv ...") -> shutil.move
        for name in ("KPOINTS", "POSCAR", "INCAR", "POTCAR"):
            if os.path.isfile(name):
                shutil.move(name, os.path.join(relax_dir, name))
        # Update INCAR with vasp.in
        if os.path.isfile("vasp.in"):
            shutil.copy("vasp.in", relax_dir)
            # FIX(all): os.chdir()/os.system() -> subprocess.run(cwd=...)
            proc = subprocess.run([sys.executable, "-m", "htesp.vasp_process",
                                   "POSCAR"], cwd=relax_dir, check=False)
            if proc.returncode != 0:
                print("vasp_process returned {} in {}\n".format(proc.returncode,
                                                               relax_dir))
    else:
        # FIX(8): QE keeps the raw cell, so the mesh is computed from it here.
        if evenkpt:
            print("Even kpoint mesh is utilized\n")
            k_mesh = pos_to_kpt("POSCAR",kptden,True)
        else:
            k_mesh = pos_to_kpt("POSCAR",kptden)
        if os.path.isfile("KPOINTS"):
            os.remove("KPOINTS")
    if keven:
        for i in range(3):
            if k_mesh[i]%2 == 0:
                k_mesh[i] = k_mesh[i]
            else:
                k_mesh[i] = k_mesh[i] + 1
    if calc_type in ('QE','qe'):
        obj = MpConnect()
        structure_file = structure.Structure.from_file("POSCAR")
        magnetic = input_data['pwscf_in']['magnetic']
        # Magnetic (FM) structure generation if magnetic flag is true
        if magnetic:
            # Reading magnetic moments from input file
            default_magmoms = input_data['magmom']['magmom']
            # Obtaining Ferromagnetic structure
            structure_file.add_spin_by_element(default_magmoms)
            obj.structure = structure_file
        else:
            obj.structure = structure_file
        # Get composition list
        comp_list = []
        for composition in obj.structure.composition.elements:
            comp_list.append(composition.name)
        obj.comp_list = comp_list
        # Getting k-point mesh
        obj.getkpt()
        evenkpt = input_data['download']['inp']['evenkpt']
        if evenkpt:
            print("Utilizing even kpoint mesh\n")
            obj.getevenkpt()
        # Getting maximum kinetic energy cutoff among elements
        obj.maxecut_sssp()
        obj.prefix = compound
        obj.mpid = mpid
        # Writing QE input scf.in file
        if magnetic:
            obj.setting_qeinput(magnetic=True,pseudo_dir='../../pp')
        else:
            obj.setting_qeinput(pseudo_dir='../../pp')
        os.makedirs("scf_dir", exist_ok=True)
        scf_name = "scf-{}.in".format(mpid)
        if os.path.isfile(scf_name):
            shutil.move(scf_name, os.path.join("scf_dir", scf_name))
def search(kwargs,properties):
    """
    Function to search compounds.

    Parameters
    ----------
    kwargs : dict
        Dictionary object for query to the OQMD website.
    properties : list
        Properties to extract.

    Notes
    -----
    This function uses QMPYRester to retrieve data from the OQMD (Open Quantum Materials Database) website.
    It attempts to get data with the provided query parameters. If the limit is reached, it decreases the limit
    and retries until successful.

    It then writes the retrieved data to 'download.csv' file,
    including IDs and properties specified in the list
    'properties'. It also creates 'mpid-list.in' file containing
    IDs, names, and compounds for each entry.
    Example
    -------
    >>> kwargs = {'element_set': '(Fe-Mn),B', 'stability': '<-0.1', 'natom': '<10', 'limit': 100}
    >>> properties = ['entry_id', 'name', 'enthalpy', 'composition']
    >>> search(kwargs, properties)
    """
    # Initialize QMPYRester object
    obj = qr.QMPYRester()
    limit = int(kwargs['limit'])
    orig_limit = limit
    # FIX(5): the old loop decremented ``limit`` by ``int(0.1*orig_limit)``,
    # which is 0 for orig_limit < 10 -- an infinite loop -- and, when it did
    # terminate, left ``data`` unbound so the next statement raised NameError.
    # The step is now at least 1, the number of attempts is bounded, and
    # exhaustion raises instead of falling through.
    step = max(1, int(orig_limit * 0.1))
    data = None
    last_error = None
    for attempt in range(MAX_SEARCH_RETRIES):
        if limit <= 0:
            break
        kwargs['limit'] = limit
        try:
            # Attempt to retrieve data from OQMD with provided query parameters.
            # FIX(32): verbose=True is qmpy_rester's default, and it does not
            # merely print -- it calls input('Proceed? [Y/n]:').  HTESP is
            # always run non-interactively (a batch script, the tutorial
            # runner, a process pool), so that prompt either blocks forever on
            # a terminal nobody is watching -- `htesp-tutorials` sat on QE/4
            # for 33 minutes until Ctrl-C -- or raises "EOF when reading a
            # line" when stdin is closed, which the retry loop below then
            # burned all its attempts on.  It also returns None for any answer
            # that is not Y/y/Yes/yes, so the next line would raise TypeError.
            with socket_timeout():
                response = obj.get_oqmd_phases(verbose=False, **kwargs)
            if response is None:
                raise RuntimeError(
                    "OQMD returned no response for limit={}".format(limit))
            # Extract data from response
            data = response['data']
            # Exit loop if data retrieval is successful
            break
        except Exception as exc:          # qmpy_rester raises bare urllib errors
            last_error = exc
            print("OQMD query failed with limit={} ({}); retrying with a "
                  "smaller limit\n".format(limit, exc))
            limit -= step
            if attempt + 1 < MAX_SEARCH_RETRIES and limit > 0:
                time.sleep(RETRY_BACKOFF_SECONDS * (attempt + 1))
    if data is None:
        raise RuntimeError(
            "OQMD search failed after {} attempts (limit went from {} down to "
            "{}).  Last error: {}".format(MAX_SEARCH_RETRIES, orig_limit,
                                          max(limit, 0), last_error))
    # Counter for entries
    ind = 1
    # FIX(all): os.makedirs(..., exist_ok=True) instead of isdir-then-mkdir
    os.makedirs("download", exist_ok=True)
    # Write data to CSV file and create 'mpid-list.in' file
    with open("download/download-oqmd.csv", "w") as write_download:
        write_download.write("ID,")
        for i,prop in enumerate(properties):
            if prop == 'composition':
                write_download.write("compound,")
            else:
                if i < len(properties) - 1:
                    write_download.write(prop + ",")
                else:
                    write_download.write(prop)
        write_download.write("\n")
    with open("mpid-list.in", "w") as write_mpid:
        for subdata in data:
            strings = subdata['name']
            strings = re.findall(r'[A-Z][a-z]*', strings)
            elements = kwargs['element_set']
            elements = re.findall(r'[A-Za-z]+', elements)
            exists_in_string = any(string not in elements for string in strings)
            # Write entry to 'mpid-list.in' and CSV file if it meets the criteria
            if not exists_in_string:
                write_mpid.write("v{}".format(ind) + " " + "oqmd-"+str(subdata['entry_id']) + " " + subdata['name'] + "\n")
                ind += 1
                with open("download/download-oqmd.csv", "a") as write_download:
                    write_download.write("oqmd-" + str(subdata['entry_id']) + ",")
                    for i,prop in enumerate(properties):
                        if prop == 'composition':
                            write_download.write(subdata[prop].replace(" ", "") + ",")
                        else:
                            if i < len(properties) - 1:
                                write_download.write(str(subdata[prop]) + ",")
                            else:
                                write_download.write(str(subdata[prop]))
                    write_download.write("\n")
def download(calc_type,start,end):
    """
    Function to create input files.

    Parameters
    ----------
    calc_type : str
        Type of calculations, QE or VASP.
    start : int
        Start index of the compounds to download.
    end : int
        End index of the compounds to download.

    Notes
    -----
    This function utilizes QMPYRester to retrieve data from the Open Quantum Materials Database (OQMD).
    It reads material IDs and compounds from 'mpid-list.in' file and downloads the corresponding data.
    For each compound, it creates a POSCAR file containing the compound structure.
    Then, it invokes the 'poscar_to_input' function to generate input files for QE or VASP calculations.
    Finally, it updates 'mpid.in' file with downloaded compounds.

    """
    #kwargs = {
    #    'element_set': '(Fe-Mn),B',      # composition include (Fe OR Mn) AND O
    #    'stability': '<-0.1',            # hull distance smaller than -0.1 eV
    #    'natom': '<10',                  # number of atoms less than 10
    #    }
    obj = qr.QMPYRester()
    # Read the material IDs and compounds from 'mpid-list.in' file
    with open("mpid-list.in","r") as read_mpid:
        mpid_data = read_mpid.readlines()
    mpid_data = mpid_data[start-1:end-1]
    # Loop through the specified range of compounds
    for mpid in mpid_data:
        oqmd_id = int(mpid.split(" ")[1].split("-")[1])
        with socket_timeout():
            subd = obj.get_entry_by_id(oqmd_id)
        try:
            # Extract necessary data from the OQMD entry
            oqmd_id = subd['id']
            oqmd_id = "oqmd-" + str(oqmd_id)
            compound = subd['name']
            # Write POSCAR file with the compound structure
            with open('POSCAR', 'w') as write_poscar:
                write_poscar.write(compound + '\n')
                write_poscar.write("1.0 \n")
                unit_cell = subd['unit_cell']
                for lattice in unit_cell:
                    write_poscar.write(str(lattice[0]) + " " + str(lattice[1]) + " " + str(lattice[2]) + "\n")
                basis = []
                elements = []
                for i,elm in enumerate(subd['sites']):
                    elements.append(subd['sites'][i].split("@")[0].split(" ")[0])
                    basis.append(subd['sites'][i].split("@")[1].split(" ")[1:])
                # FIX(6): REGRESSION GUARD -- the counts line was written from
                # Counter(elements), which groups by species, while the
                # positions were written in the original OQMD site order.
                # OQMD does not guarantee that sites arrive grouped by
                # species, so for e.g. sites [Mg, B, Mg] VASP read the header
                # "Mg B / 2 1" against positions (Mg, B, Mg) and assigned the
                # second Mg position to boron.  The sites are now emitted
                # grouped by species, in first-appearance order, so the counts
                # and the positions describe the same structure.
                species_order = []
                grouped = {}
                for i, elm in enumerate(elements):
                    if elm not in grouped:
                        species_order.append(elm)
                        grouped[elm] = []
                    grouped[elm].append(basis[i])
                write_poscar.write(" ".join(species_order) + " \n")
                write_poscar.write(
                    " ".join(str(len(grouped[elm])) for elm in species_order) + " \n")
                write_poscar.write("Direct\n")
                for elm in species_order:
                    for bas in grouped[elm]:
                        write_poscar.write(bas[0] + " " + bas[1] + " " + bas[2] + " " + elm + "\n")
            # Generate input files for QE or VASP calculations
            poscar_to_input(calc_type,oqmd_id,compound,False)
            # FIX(all): os.system("mv ...") -> os.replace, no shell
            if os.path.isfile("POSCAR"):
                os.replace("POSCAR", "POSCAR-{}".format(oqmd_id))
            # Update 'mpid.in' file with downloaded compounds
            if not os.path.isfile('mpid.in'):
                k_ind = 0
                with open("mpid.in", "w") as write_mpid:
                    write_mpid.write("v{}".format(k_ind+1) + " " + oqmd_id + " " + compound + "\n")
            else:
                with open('mpid.in', 'r') as read_mpid:
                    lines = read_mpid.readlines()
                k_ind = len(lines)
                if not any(oqmd_id in line for line in lines):
                    with open("mpid.in", "a") as write_mpid:
                        write_mpid.write("v{}".format(k_ind+1) + " " + oqmd_id + " " + compound + "\n")
        # FIX(all): bare except -> a named exception with the reason logged
        except (KeyError, IndexError, TypeError, ValueError, OSError) as exc:
            print("Structure data not found for {} ({}: {})\n".format(
                oqmd_id, type(exc).__name__, exc))
            continue
def build_kwargs(input_data=None):
    """Assemble the OQMD query dictionary and the property list.

    Returns ``(kwargs, properties, calc_type)``.  This used to live inside
    ``if __name__ == "__main__":``, which is why none of it could be tested or
    reused -- see FIX(9).
    """
    if input_data is None:
        input_data = config()
    d = input_data['download']
    doqmd = d['oqmd']
    limit = doqmd['limit']
    ntype = doqmd['ntype_constraint']
    elm_list = doqmd['entries']
    # element_set syntax: '(A-B)' means A OR B, ',' joins AND-ed terms.
    elms = '({})'.format("-".join(str(e) for e in elm_list))
    must_include = doqmd['must_include']
    # FIX(7): the old loop wrote '(A-B),MgB,C' for must_include ['Mg','B','C'] --
    # its 'j < len-1 and j > 0' test never fired for the first element, so no
    # comma was emitted after it and two species were glued into one token.
    if must_include:
        elms = elms + "," + ",".join(str(e) for e in must_include)
    properties = doqmd['prop']
    metal = doqmd['metal']
    band_gap = 0.001 if metal else None
    neg_fe = doqmd['FE']
    form_e = 0.0001 if neg_fe else None
    thermo_stable = doqmd['thermo_stable']
    ebh = 0.001 if thermo_stable else None
    nsites = doqmd['size_constraint']
    kwargs = {
               'element_set': '{}'.format(elms),
               'band_gap': '<{}'.format(band_gap),
               'delta_e': '<{}'.format(form_e),
               'stability': '<{}'.format(ebh),
               'natoms': '<{}'.format(nsites),
               'ntypes': '<{}'.format(ntype),
               'limit' : limit,
              }
    if not metal:
        kwargs.pop('band_gap')
    if not neg_fe:
        kwargs.pop('delta_e')
    if not thermo_stable:
        kwargs.pop('stability')
    calc_type = d['inp']['calc']
    return kwargs, properties, calc_type


def main(argv=None):
    """Entry point for ``oqmd_extract search`` and ``oqmd_extract download``.

    FIX(9): the whole driver -- reading ``input.in``, reading ``config.json``
    and building KWARGS -- used to sit in the ``if __name__ == "__main__":``
    block, so importing this module gave no way to run it and the names it
    defined (LIMIT, CALC_TYPE, ...) were module globals.  ``__main__`` now
    only calls this function; the ``search`` / ``download`` argument contract
    is unchanged.
    """
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print("wrong options. Allowed options are either search or download\n")
        return 1
    condition = argv[0]
    input_data = config()
    # FIX(9): input.in is parsed by htesp.inputin.InputIn rather than by hand.
    inp = InputIn.load("input.in", input_data)
    start = inp.start
    end = inp.end
    kwargs, properties, calc_type = build_kwargs(input_data)
    print(calc_type)
    if condition == 'search':
        search(kwargs, properties)
    elif condition == 'download':
        download(calc_type, start, end)
    else:
        print("wrong options. Allowed options are either search or download\n")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
