#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to process vasp input files"""
import os
import shutil
import sys
import warnings
import numpy as np
from ase.io import vasp,espresso
from ase.cell import Cell
from pymatgen.core import structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp.sets import Incar
from htesp.cif_to_gsinput import pos_to_kpt
from htesp.check_json import config

# FIX(25): the MAGMOM branch fires for every ISPIN = 2 run and used to index
# config['magmom']['magmom'] directly, so any element the user had not listed
# raised KeyError and aborted the run.  Elements that are absent get this
# documented non-magnetic guess instead (VASP's own default is 1.0 per atom;
# 0.6 is small enough not to force a moment but large enough to break the
# symmetry of a spin-polarised start).
NONMAGNETIC_MAGMOM = 0.6


def incar_key(line):
    """Return the INCAR key a line sets, upper-cased, or ``None``.

    A line sets a key when the text before the first ``=`` is a single word.
    This is what ``sed -i '/KEY/d'`` was meant to match.
    """
    if "=" not in line:
        return None
    head = line.split("=", 1)[0].strip()
    if not head or len(head.split()) != 1:
        return None
    return head.upper()


def incar_value(line):
    """Return the text a line assigns to its key, stripped."""
    return line.split("=", 1)[1].strip() if "=" in line else ""


def read_incar_lines(path="INCAR"):
    """Return the lines of an INCAR, or an empty list when it is missing."""
    try:
        with open(path, "r") as handle:
            return handle.read().splitlines()
    except OSError:
        return []


def write_incar_lines(lines, path="INCAR"):
    """Write ``lines`` back to an INCAR."""
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + ("\n" if lines else ""))


def drop_incar_keys(keys, path="INCAR"):
    """Remove every line of ``path`` that sets one of ``keys``.

    FIX(23): ``sed -i '/ENCUT/d' INCAR`` is a substring match, so deleting
    ENCUT also deleted ENCUTGW and deleting NELM also deleted NELMIN and
    NELMDL.  The key is matched here at the start of the line, followed by
    optional whitespace and ``=``.
    """
    wanted = {key.strip().upper() for key in keys}
    lines = read_incar_lines(path)
    kept = [ln for ln in lines if incar_key(ln) not in wanted]
    write_incar_lines(kept, path)
    return len(lines) - len(kept)


def parse_vasp_in(path="vasp.in"):
    """Read ``vasp.in`` into an ordered list of ``(key, value)`` pairs.

    ``value`` is ``None`` for a key-only line, which means "delete this key
    from INCAR".

    FIX(24): the old parser appended keys and values to two independent lists
    and then paired them by index, so a key-only (delete) line anywhere but at
    the end shifted every later value onto the wrong key.
    """
    # FIX(24): keys and values went into two independent lists and were paired
    # by index, so a key-only (delete) line shifted every later value onto the
    # wrong key; one pass builds the pairs directly.
    pairs = []
    try:
        with open(path, 'r') as read_vasp:
            lines = read_vasp.readlines()
    except OSError:
        return pairs
    for line in lines:
        fields = line.split()
        if not fields:
            continue
        key = fields[0]
        value = " ".join(fields[1:]) if len(fields) > 1 else None
        pairs.append((key, value))
    return pairs


def encut_check():
    """
    Function to check ENCUT in INCAR file and replace
    with 1.3xENCUTMAX of POTCAR, when ENCUT is less than that.
    """
    # Obtain ENMAX from POTCAR (was: grep ENMAX | awk '{print $3}' | sed 's/;//')
    encut = []
    try:
        with open("POTCAR", "r") as read_potcar:
            for line in read_potcar:
                if "ENMAX" not in line:
                    continue
                fields = line.split()
                if len(fields) >= 3:
                    try:
                        encut.append(float(fields[2].rstrip(";")))
                    except ValueError:
                        continue
    except OSError as exc:
        print("cannot read POTCAR ({}); ENCUT left untouched\n".format(exc))
        return
    if not encut:
        print("No ENMAX found in POTCAR; ENCUT left untouched\n")
        return
    # Calculate 1.3 times of maximum ENCUT found in POTCAR
    encutmax = max(encut)*1.3
    # Get ENCUT from current INCAR, set 0 if not found
    current = None
    for line in read_incar_lines():
        # FIX(23): ``grep ENCUT INCAR`` also matched ENCUTGW
        if incar_key(line) == "ENCUT":
            try:
                current = float(incar_value(line).split()[0])
            except (IndexError, ValueError):
                current = None
            break
    if current is None:
        print("No ENCUT parameters found in INCAR or vasp.in\n")
        print("1.3 times ENMAX is used from POTCAR\n")
        current = 0.0
    # If ENCUT in INCAR is less than encutmax, update it with ecutmax
    if current < encutmax:
        print("------------------------------------------------------\n")
        print("ENCUT less than 1.3 times ENCUTMAX in POTCAR found\n")
        print("Adjusting ENCUT in INCAR .........................\n")
        print("------------------------------------------------------\n")
        drop_incar_keys(["ENCUT"])
        write_incar_lines(read_incar_lines() + ["ENCUT = {}".format(encutmax)])


def vasp_process():
    """
    Processes the INCAR file using the 'vasp.in' file.
    Look for 'vasp.in' file in the utility folder as a demo.

    If 'vasp.in' is found, it reads the file and replaces corresponding keys in the 'INCAR' file with the provided values.
    Existing keys such as ISPIN, MAGMOM, and LORBIT in 'INCAR' are removed to avoid conflicts.
    If magnetic moment values are provided in 'config.json', it updates the 'MAGMOM' keyword accordingly.
    If 'METAGGA' is present, 'LMIXTAU = .TRUE.' is added to 'INCAR'.
    To remove unwanted keywords from 'INCAR', put those keywords at the bottom in 'vasp.in' file after the keyword that has values.

    Returns:
        None
    """
    input_data = config()
    # Get the type of magnetic enumeration
    magenum = input_data['magmom']['type']
    # Check if 'vasp.in' file exists
    if os.path.isfile("vasp.in"):
        pairs = parse_vasp_in("vasp.in")
        if magenum == 'anisotropy':
            print("anisotropy type found in magmom dictionary\n")
            # If magnetic enumeration is anisotropy, ignore NSW from 'vasp.in'.
            # This used to rewrite the user's vasp.in in place with
            # ``sed -i '/NSW/d'``; dropping it from the parsed pairs has the
            # same effect on INCAR without destroying the input file.
            pairs = [(key, value) for key, value in pairs
                     if key.strip().upper() != "NSW"]
        keys = [key for key, _ in pairs]
        key_set = {key.strip().upper() for key in keys}
        # Remove all the keys in INCAR which are found in vasp.in
        drop_incar_keys(keys)
        newincar = read_incar_lines()
        backupkey = {'ISPIN':1,'MAGMOM':False,'LORBIT':False}
        for oldkey in newincar:
            keyx = incar_key(oldkey)
            if keyx == 'ISPIN':
                try:
                    backupkey[keyx] = int(incar_value(oldkey).split()[0])
                except (IndexError, ValueError):
                    backupkey[keyx] = 1
            elif keyx == 'LORBIT':
                backupkey['LORBIT'] = True
            elif keyx == 'MAGMOM':
                backupkey['MAGMOM'] = True
        # Remove MAGMOM keyword if exists and magenum is not 'anisotropy'
        if backupkey['MAGMOM'] and magenum != 'anisotropy':
            drop_incar_keys(["MAGMOM"])
        # Remove NSW keyword if magenum is 'anisotropy'
        if magenum == 'anisotropy':
            drop_incar_keys(["NSW"])
        if backupkey['LORBIT']:
            drop_incar_keys(["LORBIT"])
        # Write new keys and values to 'INCAR'
        metagga = False
        spinval = 1
        added = []
        for key, value in pairs:
            if value is None:
                # key-only line: the key was deleted from INCAR above
                continue
            if key == "ISPIN":
                spinval = value
            if key == "METAGGA":
                metagga = True
            added.append(key + " = " + str(value))
        try:
            spin_is_two = int(str(spinval).split()[0]) == 2
        except (IndexError, ValueError):
            spin_is_two = False
        # If ISPIN=2, update MAGMOM keyword from config.json
        if ("ISPIN" in key_set and spin_is_two) or backupkey['ISPIN'] == 2:
            magmoms = input_data.get('magmom', {}).get('magmom', {}) or {}
            if magmoms:
                struc = structure.Structure.from_file("POSCAR")
                sites = struc.sites
                pieces = []
                for site in sites:
                    element = str(site.specie)
                    # FIX(25): ``m[element]`` raised KeyError for any element
                    # the user had not listed under magmom.magmom.
                    if element in magmoms:
                        moment = magmoms[element]
                    else:
                        warnings.warn(
                            "no magmom.magmom entry for element {!r}; using the "
                            "non-magnetic default of {}".format(
                                element, NONMAGNETIC_MAGMOM), RuntimeWarning)
                        moment = NONMAGNETIC_MAGMOM
                    # if LSORBIT key found, rewrite MAGMOM in mx my mz format
                    if 'LSORBIT' not in key_set:
                        pieces.append(str(moment))
                    else:
                        pieces.append("0 0 " + str(moment))
                magmom_string = " ".join(pieces)
                if magenum != 'anisotropy':
                    added.append("MAGMOM = {}".format(magmom_string))
                    added.append("LORBIT = 11")
                else:
                    print("type is anisotropy, therefore doesn't update MAGMOM keyword\n")
            else:
                print("magmom values not provided in config.json\n")
                print("Provide magnetic moment values as dictionary magmom={'A':2, 'B':3}\n")
        # Add LMIXTAU = .TRUE. if METAGGA is present in INCAR
        if metagga:
            added.append("LMIXTAU = .TRUE.")
            added.append("LASPH = .TRUE.")
        write_incar_lines(read_incar_lines() + added)
    # Check and adjust ENCUT
    encut_check()


def eigen_process():
    """
    Process the EIGENVAL file generated by VASP.

    Reads the 'input.in' file to determine if the 'vasp-line' keyword is present, indicating
    a line-mode calculation.

    Reads the 'EIGENVAL' file to extract band structure data and writes it to 'band.dat'.

    If 'vasp-line' is present or the weight of a k-point is less than 0.000001, it writes the k-point index
    and corresponding band energies to 'band.dat'. Otherwise, it skips the k-point and its associated data.

    Additionally, it creates 'band.dat.gnu', a file suitable for plotting band structures in GNUPlot.

    Returns:
        None
    """
    with open("../../input.in","r") as read_inputin:
        inputline = read_inputin.readlines()
    vasp_line = False
    for line in inputline:
        if "vasp-line" in line:
            vasp_line = True
    if os.path.isfile("KPT_OPT"):
        vasp_line = True
    # TO DO: parse EIGENVAL for spin-polarized calculations
    incar = Incar.from_file("INCAR")
    if 'LSORBIT' not in incar.keys():
        incar['LSORBIT'] = False
    if 'ISPIN' not in incar.keys():
        incar['ISPIN'] = 1
    if 'LNONCOLLINEAR' not in incar.keys():
        incar['LNONCOLLINEAR'] = False
    if incar['LSORBIT'] or incar['LNONCOLLINEAR']:
        nspin = 4
    elif incar['ISPIN'] == 2 and not incar['LNONCOLLINEAR']:
        nspin = 2
    else:
        nspin = 1
    with open("EIGENVAL", "r") as read_eig:
        for i in range(5):
            read_eig.readline()
        _, nkpoints, nbands = [int(eig) for eig in read_eig.readline().split()]
        read_eig.readline()
        with open('band.dat', 'w') as write_eig:
            k_ind = 1
            for i in range(nkpoints):
                _, _, _, weight = [float(eig) for eig in read_eig.readline().split()]
                if weight < 0.000001 or vasp_line:
                    write_eig.write(str(k_ind) + " ")
                if weight < 0.000001 or vasp_line:
                    for j in range(nbands):
                        fields = read_eig.readline().split()
                        if j < nbands - 1:
                            if nspin == 2:
                                write_eig.write(str(fields[1]) + " " + str(fields[2]) + " ")
                            else:
                                write_eig.write(str(fields[1]) + " ")
                        else:
                            if nspin == 2:
                                write_eig.write(str(fields[1]) + " " + str(fields[2]) + "\n")
                            else:
                                write_eig.write(str(fields[1]) + "\n")
                    read_eig.readline()
                    k_ind += 1
                else:
                    for j in range(nbands):
                        fields = read_eig.readline().split()
                    read_eig.readline()
    with open("band.dat.gnu", "w") as write_dat_gnu:
        data = np.loadtxt('band.dat')
        nrow,_ = data.shape
        if nspin == 2:
            band_value = data[:,1:]
            band_value_1 = band_value[:, ::2]
            band_value_2 = band_value[:, 1::2]
            band_value = band_value_1
        else:
            band_value = data[:,1:]
        for i in range(band_value.shape[1]):
            for j in range(band_value.shape[0]):
                if j < nrow - 1:
                    if nspin == 2:
                        write_dat_gnu.write(str(j) + " " + str(band_value_1[j,i]) + " " + str(band_value_2[j,i]) + "\n")
                    else:
                        write_dat_gnu.write(str(j) + " " + str(band_value[j,i]) + "\n")
                else:
                    if nspin == 2:
                        write_dat_gnu.write(str(j) + " " + str(band_value_1[j,i]) + " " + str(band_value_2[j,i]) + "\n")
                        write_dat_gnu.write("\n")
                    else:
                        write_dat_gnu.write(str(j) + " " + str(band_value[j,i]) + "\n")
                        write_dat_gnu.write("\n")


def split_path_labels(path_string):
    """Split an ASE band-path string into its labels.

    ``'GXWK1G'`` -> ``['G', 'X', 'W', 'K1', 'G']``: a digit belongs to the
    label in front of it.
    """
    labels = []
    for char in path_string:
        if char.isdigit() and labels:
            labels[-1] += char
        else:
            labels.append(char)
    return labels


def band_phonopy(filename):
    """
    Write k-points of high-symmetry points for Phonopy band structure plot using 'band.conf' file.

    Parameters:
    -----------
    filename : str
        The name of the file containing the crystal structure, either 'POSCAR' for VASP or 'scf.in' for Espresso.

    Returns:
    --------
    None
    """
    if filename == 'POSCAR':
        data = vasp.read_vasp('POSCAR')
    else:
        data = espresso.read_espresso_in('scf.in')
    cell_bp = Cell.bandpath(data.cell)
    # FIX(26): the ',' in an ASE path marks a discontinuity (e.g. 'GXWK,UX').
    # Stripping it joined the two branches into one, so phonopy interpolated a
    # spurious segment between the end of one branch and the start of the next.
    # Each comma-separated branch becomes its own phonopy segment instead.
    segments = [split_path_labels(chunk)
                for chunk in cell_bp.path.split(',') if chunk]
    npt = sum(len(segment) for segment in segments)
    band = data.cell.bandpath(path=cell_bp.path,npoints=npt)
    specialpt = band.special_points
    with open("band_phonopy.in", "w") as band_phon:
        band_phon.write("BAND = ")
        rendered = []
        for segment in segments:
            points = []
            for special in segment:
                coords = specialpt[special]
                points.append("{} {} {}".format(coords[0], coords[1], coords[2]))
            rendered.append("\t".join(points))
        # phonopy separates disconnected branches of BAND/BAND_LABELS with ','
        band_phon.write(", ".join(rendered))
        band_phon.write("\n")
        band_phon.write("BAND_LABELS = ")
        band_phon.write(", ".join("\t".join(segment) for segment in segments))
    with open("high_symm.in", "w") as high_sym:
        special = list(Cell.bandpath(data.cell).special_points.values())
        for j,_ in enumerate(special):
            high_sym.write(str(special[j][0]) + " ")
            high_sym.write(str(special[j][1]) + " " + str(special[j][2]) + "\n")


def main():
    """
    main function
    """
    input_data = config()
    filename = sys.argv[1]
    kptden = input_data["kptden"]
    if filename == 'POSCAR':
        band_phonopy(filename)
        vasp_process()
        # remove EIGENVAL if you only update INCAR and KPOINTS
        if os.path.isfile("EIGENVAL"):
            shutil.copy("KPOINTS", "KPOINTS_band")
            eigen_process()
        pos_to_kpt(filename,kptden)
    # Change POSCAR to conventional unit cell
    # and update KPOINTS
    elif filename == 'conventional':
        data = structure.Structure.from_file("POSCAR")
        data = SpacegroupAnalyzer(data,symprec=0.1).get_conventional_standard_structure()
        data.to("POSCAR")
        pos_to_kpt("POSCAR",kptden)
    elif filename == 'eigen':
        # FIX(39): KPOINTS_band is written by the 'POSCAR' branch above, but
        # only when EIGENVAL is present -- that is, only after a real VASP run.
        # Moving it unconditionally turned "the band run has not happened yet"
        # into a bare "[Errno 2] No such file or directory: 'KPOINTS_band'".
        if not os.path.isfile("KPOINTS_band"):
            raise FileNotFoundError(
                "KPOINTS_band is missing in {}.  It is created from KPOINTS by "
                "'vasp_process.py POSCAR' once EIGENVAL exists, so run the "
                "band-structure calculation first (mainprogram 13) and only "
                "then post-process it (mainprogram 15).".format(os.getcwd()))
        shutil.move("KPOINTS_band", "KPOINTS")
    # Symmetrize the POSCAR
    elif filename == 'symmetrize':
        print("Symmetrizing the primitive structure\n")
        data = structure.Structure.from_file("POSCAR")
        data = SpacegroupAnalyzer(data,symprec=0.1).get_primitive_standard_structure()
        data.to("POSCAR")
        pos_to_kpt("POSCAR",kptden)
    else:
        band_phonopy(filename)


if __name__ == "__main__":
    main()
