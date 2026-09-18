#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to write QE/VASP input files from .cif files.

FIX(19) -- DUPLICATED FUNCTION.  ``pos_to_kpt`` below is a diverged copy of
``htesp.htepc.pos_to_kpt``.  They are NOT interchangeable:

* this copy takes a third argument ``evenkpt``; when true it rounds every odd
  mesh division up to the next even number,
* this copy writes a ``KPOINTS`` file in the current directory as a side
  effect and returns the mesh, whereas ``htepc.pos_to_kpt`` only returns it,
* ``htepc.pos_to_kpt`` computes ``kratio``/``klat``/``kmesh`` twice (a
  copy-paste leftover; the second computation is identical to the first).

The numerical core is otherwise the same.  They are deliberately left
unmerged in this pass; ``htesp/htepc.py`` is maintained separately.  Callers
that need ``KPOINTS`` on disk or an even mesh must use this one.
"""
import sys
import os
import shutil
import subprocess
import tempfile
import warnings
import glob
import scipy.linalg as alg
from ase.io import vasp
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import structure
import numpy as np
from htesp.write_potcar import poscar2potcar
from htesp.htepc import MpConnect
from htesp.check_json import config
def pos_to_kpt(structure_filename,kpoint_density,evenkpt=False):
    """
    Function to obtain a k-point mesh from a structure file.

    .. note:: FIX(19) -- diverged duplicate of :func:`htesp.htepc.pos_to_kpt`;
       that one has no ``evenkpt`` argument and does not write ``KPOINTS``.

    Parameters:
    - structure_filename (str): Path to the structure file (QE scf.in or VASP POSCAR).
    - kpoint_density (float): Desired k-point density.
    - evenkpt (bool): Flag indicating whether to enforce an even number of k-points along each axis.
      Default is False.

    Returns:
    - kmesh (list): K-point mesh according to the k-point density.
    Example:
        >>> # Generating a k-point mesh for a VASP POSCAR with a k-point density of 0.05
        >>> pos_to_kpt("POSCAR", 0.05, evenkpt=True)
    """
    kptsp = kpoint_density
    # Read the structure file
    with open(structure_filename,"r") as read_struc:
        lines = read_struc.readlines()
    # Get cell vectors and compute reciprocal lattice vectors
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
    #bnorm = np.zeros(3)
    # Compute the length of reciprocal lattice vectors
    bnorm = alg.norm(bmat,axis=1)
    kratio = [bnorm[i] / bnorm[0] for i in range(3)]
    # Compute the k-point mesh
    klat = bnorm[0] / kptsp
    kmesh = [int(kratio[i] * klat + 0.5) if int(kratio[i] * klat + 0.5) != 0 else 1 for i in range(3)]
    # Enforce even number of kpoints if evenkpt flag is true
    if evenkpt:
        for i in range(3):
            if kmesh[i]%2 == 0:
                kmesh[i] = kmesh[i]
            else:
                kmesh[i] = kmesh[i] + 1
    # Write KPOINTS file
    with open("KPOINTS", "w") as write_kpt:
        write_kpt.write("KPOINTS automatic" + "\n")
        write_kpt.write("0"+"\n")
        write_kpt.write("Gamma\n")
        write_kpt.write(f"{kmesh[0]} {kmesh[1]} {kmesh[2]}"+"\n")
        write_kpt.write("0 0 0\n")
    return kmesh
def read_mpid_entries(path="mpid.in"):
    """Return ``[(mpid, compound), ...]`` from a ``v<N> <mpid> <compound>`` file.

    A missing or unreadable file is an empty list (FIX(18)).
    """
    entries = []
    try:
        with open(path, "r") as read_track:
            lines = read_track.readlines()
    except OSError:
        return entries
    for line in lines:
        parts = line.split()
        if len(parts) >= 3:
            entries.append((parts[1], parts[2]))
        elif len(parts) == 2:
            entries.append((parts[1], ""))
    return entries
def find_mpid(mpid, path="mpid.in"):
    """1-based index of ``mpid`` in the tracking file, or ``None`` (FIX(18))."""
    for index, (existing, _) in enumerate(read_mpid_entries(path)):
        if existing == str(mpid):
            return index + 1
    return None
def register_mpid(mpid, compound, path="mpid.in"):
    """Record ``mpid`` in the ``v<N> <mpid> <compound>`` tracking file.

    FIX(18): the original read-then-append pattern (repeated verbatim in
    ``qe_input.py``, ``vasp_input.py``, ``oqmd_extract.py`` and
    ``aflow_extract.py``) numbered the new entry ``v<len(lines)+1>``, so a file
    that had ever been hand-edited, or two processes appending at once,
    produced duplicate or skipped ``v<N>`` numbers -- and every consumer looks
    an entry up with ``grep "v$ii "``.  This helper

    * is idempotent -- an mpid already present is left alone and its existing
      index returned,
    * renumbers ``v<N>`` densely from 1 on every write, and
    * is atomic -- the whole file is written to a temporary file in the same
      directory and ``os.replace``d into place.

    Returns the 1-based index of ``mpid`` in the file.
    """
    entries = read_mpid_entries(path)
    for index, (existing, _) in enumerate(entries):
        if existing == str(mpid):
            return index + 1
    entries.append((str(mpid), str(compound)))
    directory = os.path.dirname(os.path.abspath(path)) or "."
    handle, tmp_path = tempfile.mkstemp(dir=directory, prefix=".mpid-", suffix=".tmp")
    try:
        with os.fdopen(handle, "w") as write_track:
            for index, (entry_mpid, entry_comp) in enumerate(entries):
                write_track.write("v{} {} {}\n".format(index + 1, entry_mpid, entry_comp))
        os.replace(tmp_path, path)
    except BaseException:
        if os.path.isfile(tmp_path):
            os.remove(tmp_path)
        raise
    return len(entries)
def pymatgen_cif(infile):
    """
    Function to convert (ICSD) CIF files to pymatgen format.

    Parameters:
    - infile (str): Path to the input CIF file.

    Returns:
    - structure (Structure): Pymatgen Structure object representing the CIF structure.
    Example:
        >>> # Convert a CIF file named 'example.cif' to pymatgen format
        >>> pymatgen_cif("example.cif")
    """
    # Parse the CIF file using CifParser
    cif_parser = CifParser(infile)
    structure = cif_parser.parse_structures()[0]  # Assuming there's only one structure in the CIF
    oxidation = False
    # Check the oxidation states and merge sites if present
    for elem in structure.elements:
        if '+' in str(elem):
            oxidation = True
    if oxidation:
        structure.merge_sites()
        structure = structure.remove_oxidation_states()
    # Write the CIF file using CifWriter
    cif_writer = CifWriter(structure, symprec=0.1)
    cif_writer.write_file(infile)
    # Get the primitive standard structure using SpacegroupAnalyzer
    structure = SpacegroupAnalyzer(structure=structure,symprec=0.1).get_primitive_standard_structure()
    return structure
def ciftoscf(calc_type,mpid,cif2cell=True,keven=False):
    """
    Function to convert a CIF file to a Quantum ESPRESSO (QE) input file 'scf.in'.

    Parameters:
    - calc_type (str): The type of calculation, either 'VASP' or 'QE'.
    - mpid (str): The Material ID.
    - cif2cell (bool): If True, uses cif2cell to process the CIF file. Default is True.
    - keven (bool): If True, even k-mesh is used. Default is False.

    Returns:
    - compound (str): The compound name.

    Writes the 'scf-mpid.in' file.
    Example:
        >>> # Convert a CIF file to a Quantum ESPRESSO input file
        >>> ciftoscf("QE", "mp-1234", cif2cell=False, keven=False)
    """
    input_data = config()
    # Get current directory
    pwd = os.getcwd()
    # Check if CIF file exists in the current directory
    # glob() returns an empty list rather than raising, so the old
    # try/except FileNotFoundError left ``file_path`` unbound.
    matches = glob.glob(os.path.join(pwd, "{}.cif".format(mpid)))
    if not matches:
        raise FileNotFoundError("File {}.cif not found in {}".format(mpid, pwd))
    file_path = matches[0]
    # Process CIF file using cif2cell if use_cif2cell flag is true
    # Otherwise, it uses pymatgen
    if not cif2cell:
        # Use pymatgen to process CIF file and obtain structure
        struc_poscar = pymatgen_cif(file_path)
        cell = struc_poscar.lattice.matrix
        compound = str(struc_poscar.composition).replace(" ", "")
        pos = struc_poscar.frac_coords
        symbol = []
        for site in struc_poscar.sites:
            symbol.append(str(site.specie))
        # Write POSCAR file
        with open("POSCAR", "w") as write_poscar:
            write_poscar.write("Structure for {}".format(mpid) + "\n")
            write_poscar.write("1.0\n")
            for i in range(3):
                write_poscar.write("{} {} {}".format(cell[i][0], cell[i][1], cell[i][2]) + "\n")
            dict_symbol = {}
            for i,sym in enumerate(symbol):
                if sym not in dict_symbol:
                    dict_symbol[sym] = 1
                else:
                    dict_symbol[sym] += 1
            for symb in dict_symbol:
                write_poscar.write(symb + " ")
            write_poscar.write("\n")
            for symb in dict_symbol:
                write_poscar.write(str(dict_symbol[symb]) + " ")
            write_poscar.write("\n")
            write_poscar.write("Direct\n")
            for i,_ in enumerate(symbol):
                write_poscar.write(str(pos[i][0]) + " ")
                write_poscar.write(str(pos[i][1]) + " " + str(pos[i][2]) + " " + symbol[i] +"\n")
    else:
        # FIX(all): os.system -> subprocess.run with an argument list, so the
        # CIF path is never interpreted by a shell.
        subprocess.run(["cif2cell", file_path, "-p", "vasp",
                        "--vasp-cartesian-lattice-vectors"], check=True)
        with open("POSCAR", "r") as read_poscar:
            lines = read_poscar.readlines()
        list_first = lines[0].split("\n")[0].split(" ")
        index = None
        for i,element in enumerate(list_first):
            if "order:" in element:
                index = i
        if index is None:
            raise ValueError("cif2cell POSCAR header has no 'order:' field: "
                             + lines[0].rstrip())
        ion_name = list_first[index+1:-1]
        insert_text = " ".join(ion_name)
        # FIX(all): sed/mv via the shell replaced by an in-memory line insert
        # (sed '5 a <text>' appends after line 5, i.e. before old line 6).
        lines.insert(5, insert_text + "\n")
        with open("POSCAR", "w") as write_poscar:
            write_poscar.writelines(lines)
        symbol = ion_name
        data = vasp.read_vasp('POSCAR')
        compound = str(data.symbols)
    # Obtain k-point mesh
    evenkpt = input_data['download']['inp']['evenkpt']
    kptden = input_data['kptden']
    if evenkpt:
        k_mesh = pos_to_kpt("POSCAR",kptden,True)
    else:
        k_mesh = pos_to_kpt("POSCAR",kptden)
    # Generate VASP input set
    relax_dir = os.path.join("R{}-{}".format(mpid, compound), "relax")
    if calc_type in ('VASP','vasp'):
        # FIX(all): os.makedirs(..., exist_ok=True), no shell
        os.makedirs(relax_dir, exist_ok=True)
        structure_file = structure.Structure.from_file("POSCAR")
        # Generate MPRelaxSet object
        relax_set = MPRelaxSet(structure=structure_file)
        relax_set.poscar.write_file("POSCAR")
        # Generate POTCAR from POSCAR file
        poscar2potcar()
        relax_set.incar.write_file("INCAR")
        # FIX(all): os.system("mv ...") -> shutil.move
        for name in ("KPOINTS", "POSCAR", "INCAR", "POTCAR"):
            if os.path.isfile(name):
                shutil.move(name, os.path.join(relax_dir, name))
        print(compound)
        # Copy vasp.in file if it exists
        if os.path.isfile("vasp.in"):
            shutil.copy("vasp.in", relax_dir)
            # Process VASP input files reading vasp.in file if it exists
            # Otherwise, it will simply print INCAR from MPRelaxSet
            # FIX(all): os.chdir()/os.system() -> subprocess.run(cwd=...)
            proc = subprocess.run([sys.executable, "-m", "htesp.vasp_process",
                                   "POSCAR"], cwd=relax_dir, check=False)
            if proc.returncode != 0:
                print("vasp_process returned {} in {}\n".format(proc.returncode,
                                                               relax_dir))
    else:
        if os.path.isfile("KPOINTS"):
            os.remove("KPOINTS")
    if keven:
        for i in range(3):
            if k_mesh[i]%2 == 0:
                k_mesh[i] = k_mesh[i]
            else:
                k_mesh[i] = k_mesh[i] + 1
    # QE calculation setup
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
    return compound
def main(calc_type=None, argv=None):
    """
    The main function of the script.

    This function performs the following operations:
    - Retrieves a list of CIF files in the current directory.
    - Sets the calculation type.
    - Sets the flag for using cif2cell.
    - Iterates over CIF files, converts them to QE input files, and writes data to 'mpid.in'.

    Parameters:
    - calc_type (str): Calculation type ('VASP' or 'QE').  FIX(20): this is now
      a real argument; it used to be read from ``sys.argv[1]`` in the middle of
      the function, so the module could not be driven from Python.  When it is
      omitted the value still comes from the command line, keeping the
      ``cif_to_gsinput.py <calc>`` contract.
    - argv (list): command line used when ``calc_type`` is omitted.

    Returns: None
    """
    input_data = config()
    warnings.filterwarnings('ignore')
    if calc_type is None:
        argv = list(sys.argv[1:] if argv is None else argv)
        if not argv:
            raise SystemExit("usage: cif_to_gsinput.py <QE|VASP>")
        calc_type = argv[0]
    # Get a list of CIF files in the current directory
    list_cif = glob.glob("*.cif",recursive=True)
    # CIF2CELL if True, uses cif2cell package to create POSCAR from given .cif files.
    # if False, It uses pymatgen cifparser to read and produce cif output, which then explicitely
    # read and uses to write POSCAR.
    cif2cell = input_data['download']['inp']['use_cif2cell']
    if cif2cell:
        print("CIF2CELL is True. Using cif2cell package....\n")
    # Iterate over CIF files
    for cif in list_cif:
        mpid = cif.split('.')[0] # Extract Material ID from file name
        print(mpid)
        # Convert CIF to QE input file
        compound = ciftoscf(calc_type,mpid,cif2cell,False)
        # FIX(18): read-then-append replaced by the atomic, idempotent,
        # densely renumbering helper above.
        register_mpid(mpid, compound)
if __name__ == '__main__':
    # FIX(20): the __main__ guard is present and calls main() with no
    # hand-rolled argument parsing of its own.
    main()
