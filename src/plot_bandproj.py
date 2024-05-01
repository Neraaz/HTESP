#!/usr/bin/env python
"""
Script to process projection data
"""
import os
import sys
import glob
from collections import Counter
import numpy as np
from pymatgen.io.pwscf import PWInput
from pymatgen.io.vasp.sets import Incar
from pymatgen.core import structure
from check_json import config
from plot import plot

class PROCARProcessor:
    """
    A class to process PROCAR file and extract band information.

    Parameters:
    - procar_file (str): Path to the PROCAR file.
    - band_file (str): Path to the band file.

    Attributes:
    - nkpoint (int): Number of k-points.
    - nband (int): Number of bands.
    - nion (int): Number of ions.
    - ncol (int): Number of columns in the PROCAR data.
    - col (list): List of column labels.
    - data (numpy.ndarray): Processed PROCAR data.
    - nspin (int): Number of spin channels.
    - elements (list): List of elements parsed from structure file.
    - nkpoint_extra (int): Additional k-points in PROCAR file not in band file.

    Methods:
    - write_band_orbital(orb): Writes band orbital data to the specified file.
    - write_band_element(elm): Writes band element data to the specified file.
    - extract_data(elm, orb): Extracts data for a specific element and orbital.
    - extract_data_2(elm, orb): Extracts data for a specific element and orbital, summing over suborbitals.
    - print_info(): Prints available options for elements and orbitals.
    - clean(): Cleans up temporary files created during processing.

    Example:
    ```python
    # Instantiate PROCARProcessor
    >>> procar_processor = PROCARProcessor("PROCAR", "BANDS")

    # Print available options
    >>> procar_processor.print_info()

    # Extract data for a specific element and orbital
    >>> procar_processor.extract_data('Fe', 's')

    # Write band orbital data to a file
    >>> procar_processor.write_band_orbital('s')

    # Write band element data to a file
    >>> procar_processor.write_band_element('Fe')

    # Clean up temporary files
    >>> procar_processor.clean()
    ```
    """
    def __init__(self, procar_file, band_file):
        """
        Initialize PROCARProcessor with the provided PROCAR file and band file.

        Parameters:
        - procar_file (str): Path to the PROCAR file.
        - band_file (str): Path to the band file.
        """
        # Assign the provided file paths
        self.procar_file = procar_file
        self.band_file = band_file
        # Parse header information from the band file and assign to attributes
        self.nkpoint, self.nband, self.nion, self.ncol, self.col = self._parse_procar_header(self.band_file)
        # Parse data from the PROCAR file and assign to attributes
        self.data, self.nspin = self._parse_procar_data()
        # Parse structure information and assign to elements attribute
        self.elements = self._parse_structure()
        # Initialize additional k-points attribute
        # Useful for METAGGA/Hybrid calculations with 2 sets of k-mesh in KPOINTS file
        # uniform k-mesh and high-symmetry k-mesh
        self.nkpoint_extra = None

    def _parse_procar_header(self,band_file):
        """
        Parse header information from the PROCAR file.

        Parameters:
        - band_file (str): Path to the band file.

        Returns:
        - nkpoint (int): Number of k-points.
        - nband (int): Number of bands.
        - nion (int): Number of ions.
        - ncol (int): Number of columns in the PROCAR data.
        - col (list): List of column labels.
        """
        # Read PROCAR file
        with open(self.procar_file, "r") as f:
            procar = f.readlines()
        # Determine number of columns
        ncol = len(procar[7].split())
        col = procar[7].split()[1:-1]
        # Load band file to get number of k-points
        orig_band = np.loadtxt(self.band_file)
        nkpoint = int(orig_band.shape[0])
        # Extract total number of k-points
        nkpoint_total = int(procar[1].split(":")[1].split()[0])
        self.nkpoint_extra = nkpoint_total - nkpoint
        # Extract number of bands and ions
        nband = int(procar[1].split(":")[2].split()[0])
        nion = int(procar[1].split(":")[3].split()[0])
        # Calculate number of skipped lines for MGGA/Hybrid calculations
        # with 2 sets of k-points in KPOINTS
        self.skip = int(nband * self.nkpoint_extra)
        return nkpoint, nband, nion, ncol, col
    def _parse_structure(self):
        """
        Parse structure information from the POSCAR file.

        Returns:
        - elements (list): List of element symbols.
        """
        # Load structure from POSCAR file
        struc = structure.Structure.from_file("POSCAR")
        # Extract element symbols
        elements = []
        for elm in struc:
            elements.append(elm.species_string)
        return elements

    def print_info(self):
        """
        Print available options for elements and orbitals.
        """
        # Find unique elements
        unique_elements = list(set(self.elements))
        orbital_names = self.col + ["p", "d"]
        # Print information
        print(f"Available options: elements {unique_elements}, orbitals {orbital_names}")

    def clean(self):
        """
        Clean up temporary files created during processing.
        """
        for elm in self.elements:
            clean_files = glob.glob(f"{elm}-*.dat.gnu")
            if len(clean_files) > 0:
                for file_ in clean_files:
                    os.system(f"rm {file_}")

    def _parse_procar_data(self):
        """
        Parse data from the PROCAR file.

        Returns:
        - data (numpy.ndarray): Processed PROCAR data.
        - nspin (int): Number of spin channels.
        """
        # Load INCAR file to determine spin configuration
        incar = Incar.from_file("INCAR")
        if 'LSORBIT' not in incar.keys():
            incar['LSORBIT'] = False
        if 'ISPIN' not in incar.keys():
            incar['ISPIN'] = 1
        if 'LNONCOLLINEAR' not in incar.keys():
            incar['LNONCOLLINEAR'] = False
        if incar['LSORBIT'] or incar['LNONCOLLINEAR']:
            self.nspin = 4
        elif incar['ISPIN'] == 2 and not incar['LNONCOLLINEAR']:
            self.nspin = 2
        else:
            self.nspin = 1
        # Read PROCAR file
        with open(self.procar_file, "r") as f:
            procar = f.readlines()
        # Initialize ncount variable
        ncount = 0
        # Write processed PROCAR data to a temporary file
        with open("procar_process.txt", "w") as procar_write:
            if self.nspin in (1, 4):
                # Loop over lines in PROCAR file
                for i, line in enumerate(procar):
                    if 'ion      s' in line:
                        # Check if the line is within the range to skip
                        if ncount >= self.skip:
                            for j in range((self.nion+1)*self.nspin):
                                procar_write.write(procar[i+j+1])
                        ncount += 1
            else:
                len_s = int(len(procar[1:])/2)
                procar_up = procar[1:len_s+1]
                procar_dn = procar[len_s+1:]
                # Loop over lines in PROCAR file for spin-polarized calculation
                for i, line in enumerate(procar_up):
                    if 'ion      s' in line:
                        if ncount >= self.skip:
                            for j in range(self.nion+1):
                                procar_write.write(procar_up[i+j+1])
                        ncount += 1
                ncount = 0
                for i, line in enumerate(procar_dn):
                    if 'ion      s' in line:
                        if ncount >= self.skip:
                            for j in range(self.nion+1):
                                procar_write.write(procar_dn[i+j+1])
                        ncount += 1
        # Read processed PROCAR data from temporary file
        data = np.genfromtxt("procar_process.txt")
        # Reshape data according to spin configuration
        if self.nspin == 2:
            data = data.reshape(self.nspin, self.nkpoint, self.nband, self.nion+1, self.ncol)
        else:
            data = data.reshape(self.nkpoint, self.nband, self.nspin, self.nion+1, self.ncol)
        # For noncollinear calculation, select only the first spin channel
        if self.nspin == 4:
            data = data[:,:,0,:,:]
        return data, self.nspin
        #return data.reshape(self.nkpoint, self.nband, self.nspin, self.nion+1, self.ncol)

    def extract_data(self,elm,orb):
        """
        Extract data for a specific element and orbital.

        Parameters:
        - elm (str): Element symbol.
        - orb (str): Orbital type.
        """
        # To Do: Handle for ISPIN = 2 only
        # Find indices of the orbital
        indices = find_all_indices(self.col,orb)
        # Extract subdata based on spin configuration
        if self.nspin == 2:
            subdata = self.data[:,:,:,find_all_indices(self.elements,elm),indices]
        else:
            subdata = self.data[:,:,find_all_indices(self.elements,elm),indices]
        # Write extracted data to a file
        with open(f"{elm}-{orb}.dat.gnu", "w") as write_procar:
            for i in range(self.nband):
                if self.nspin == 2:
                    subdata_iband = subdata[:,:,i,:]
                    subdata_iband = np.sum(subdata_iband,axis=2)
                else:
                    subdata_iband = subdata[:,i,:]
                    subdata_iband = np.sum(subdata_iband,axis=1)
                for j in range(self.nkpoint):
                    if self.nspin == 2:
                        procar_jband = subdata_iband[:,j]
                        write_procar.write(str(j+1) + " " + str(procar_jband[0]) + " " + str(procar_jband[1]) +  "\n")
                    else:
                        procar_jband = subdata_iband[j]
                        write_procar.write(str(j+1) + " " + str(procar_jband) + "\n")
                write_procar.write("\n")
    def extract_data_2(self,elm,orb):
        """
        Extract data for a specific element and orbital, summing over suborbitals.

        Parameters:
        - elm (str): Element symbol.
        - orb (str): Orbital type.
        """
        # Determine suborbitals for each orbital
        if orb == 'p':
            orb_list = ['px', 'py', 'pz']
        elif orb == 'd':
            orb_list = ['dxy', 'dyz', 'dxz', 'dz2', 'x2-y2']
        else:
            print("Supports p or d\n")
        # Extract data for each suborbital
        for orbi in orb_list:
            self.extract_data(elm,orbi)
        data0 = np.loadtxt(f"{elm}-{orb_list[0]}.dat.gnu")
        data = np.zeros_like(data0)
        data[:,0] = data0[:,0]
        # Sum data over suborbitals
        for orbi in orb_list:
            data_i = np.loadtxt(f"{elm}-{orbi}.dat.gnu")
            data[:,1] += data_i[:,1]
            if self.nspin == 2:
                data[:,2] += data_i[:,2]
        # Save summed data to a file
        np.savetxt(f"{elm}-{orb}.dat.gnu",data)

    def write_band_orbital(self, orb):
        """
        Write orbital data to the specified file.

        Parameters:
        - output_file (str): Path to the output file.
        """
        # Open file for writing orbital data
        with open(f"{orb}.dat.gnu", "w") as write_procar:
            # Loop over bands
            for i in range(self.nband):
                # Extract data for the current band and orbital
                if self.nspin == 2:
                    procar_iband = self.data[:, :, i, self.nion, 1:self.ncol]
                    procar_iband = procar_iband[:, :,find_all_indices(self.col,orb)]
                else:
                    procar_iband = self.data[:, i, self.nion, 1:self.ncol]
                    procar_iband = procar_iband[:, find_all_indices(self.col,orb)]
                # Loop over k-points
                for j in range(self.nkpoint):
                    if self.nspin == 2:
                        #procar_jband = procar_iband[j]
                        write_procar.write(str(j+1) + " " + str(procar_jband[0,j]) + " " + str(procar_jband[1,j]) + "\n")
                    else:
                        procar_jband = procar_iband[j][0]
                        write_procar.write(str(j+1) + " " + str(procar_jband) + "\n")
                write_procar.write("\n")

    def write_band_element(self, elm):
        """
        Write element data to the specified file.

        Parameters:
        - output_file (str): Path to the output file.
        """
        with open(f"{elm}.dat.gnu", "w") as write_procar:
            for i in range(self.nband):
                # Extract data for the current band and element
                if self.nspin == 2:
                    procar_iband = self.data[:, :, i, find_all_indices(self.elements,elm), -1]
                    procar_iband = np.sum(procar_iband,axis=2)
                else:
                    procar_iband = self.data[:, i, find_all_indices(self.elements,elm), -1]
                    procar_iband = np.sum(procar_iband,axis=1)
                for j in range(self.nkpoint):
                    if self.nspin == 2:
                        #procar_jband = procar_iband[j]
                        write_procar.write(str(j+1) + " " + str(procar_jband[0,j]) + " " +str(procar_jband[1,j]) + "\n")
                    else:
                        procar_jband = procar_iband[j]
                        write_procar.write(str(j+1) + " " + str(procar_jband) + "\n")
                write_procar.write("\n")

def find_all_indices(lst, element):
    """
    Find all indices of a specific element in a list.

    Parameters:
    - lst (list): List to search.
    - element (str): Element to find.

    Returns:
    - indices (list): List of indices where the element is found.
    """
    indices = []
    index = -1
    while True:
        try:
            index = lst.index(element, index + 1)
            indices.append(index)
        except ValueError:
            break
    return indices


class DataProcessor:
    """
    A class to process projection data in "proj.out*" format from QE calculations.

    Attributes:
        projdata (str): The projection data.
        nkpoint (int): The number of k-points.
        nband (int): The number of bands.
        nrow (int): The number of rows in the data.
        dict_ (dict): A dictionary for storing projection data.
        dict_2 (dict): Another dictionary for storing data with keys
                       as elements and values as list storing orbital info.
        dict_1 (dict): Similar to dict_1, but only with useful orbital info.
        ref_data (str): Reference data to initialize numpy array.
        proj_dict (dict): Dictionary for projection data.
        final_dict (dict): Dictionary for final data for unique element-orbital combination.
        In a compound like A2B3, when combining orbital contributions for elements A and B,
        we need to consider the stoichiometry.
        For example, if we're summing contributions from the pz orbital:

        - For element A, the contribution from the pz orbital should be summed 2 times
          (since there are 2 A atoms).
        - For element B, the contribution from the pz orbital should be summed 3 times
          (since there are 3 B atoms).
        This ensures that the overall contribution is correctly weighted by the number
        of atoms of each element in the compound.
    Example:
    >>> processor = DataProcessor(structure_file="scf.in", proj_file="proj.out*")
    """
    def __init__(self, structure_file="scf.in", proj_file="proj.out*"):
        """
        Initialize DataProcessor with structure and projection files.

        Args:
            structure_file (str): The filename of the structure input file.
            proj_file (str): The filename pattern for the projection output files.
        """
        self.projdata = None
        self.nkpoint = None
        self.nband = None
        self.nrow = None
        self.dict_ = None
        self.dict_1 = None
        self.dict_2 = None
        self.ref_data = None
        self.proj_dict = {}
        self.final_dict = {}
        self._parse_proj_file(structure_file, proj_file)

    def convert_gnu(self, filename):
        """
        Convert data file to GNU plot format.

        Args:
            filename (str): The filename of the data file to be converted.
        """
        # Load data from the file as output of projwfc.x format
        data = np.loadtxt(filename)
        # Extract the third column of the data
        data = data[:,2]
        # Reshape the data into a 2D array based on nkpoint and nband
        data = data.reshape(self.nkpoint, self.nband)
        # Write the data to a GNU plot format file
        with open(f"{filename}.gnu", "w") as write_gnu:
            # Loop over each band
            for i in range(data.shape[1]):
                datak = data[:,i]
                # Write data for each k-point in the band
                for j in range(datak.shape[0]):
                    write_gnu.write(str(j) + " " + str(datak[j]) + "\n")
                # Add a newline after writing data for each band
                write_gnu.write("\n")

    def combine_orbital(self, key, orb, if_gnu=False):
        """
        Combine projection data for a specific orbital and write to files.

        Args:
            key (str): Key for accessing data in dictionaries.
            orb (str): Orbital type to combine.
            if_gnu (bool): Whether to convert the combined data to GNU plot format.
        """
        # Write projection data for each element/orbital combination
        self.write_proj_elm_orbi()
        # Iterate over each orbital name in the final_dict for the given key
        for name in self.final_dict[key]:
            # Retrieve projection array for the current orbital name
            proj_array = self.final_dict[key][name]
            # Save projection array to a file with appropriate formatting
            np.savetxt(f"{key}-{name}.dat",proj_array,fmt=['%d', '%d', '%.4f'])
            # Convert to GNU plot format if specified
            if if_gnu:
                self.convert_gnu(f"{key}-{name}.dat")
        # Combine data for the specified orbital type
        sdata = np.zeros((self.nrow, 3))
        sfile = glob.glob(f"{key}-*{orb}*.dat")
        s0 = np.loadtxt(sfile[0])
        sdata[:,0] = s0[:,0]
        sdata[:,1] = s0[:,1]
        len_s = len(sfile)
        for i in range(len_s):
            idata = np.loadtxt(sfile[i])
            #sdata[:,2] += idata[:,2] / len_s
            sdata[:,2] += idata[:,2]
        np.savetxt(f"{key}-{orb}.dat", sdata, fmt=['%d', '%d', '%.4f'])
        self.convert_gnu(f"{key}-{orb}.dat")

    def clean(self,delete_file):
        """
        Deleting file
        """
        os.system(f"rm {delete_file}")

    def _parse_proj_file(self, structure_file, proj_file):
        """
        Parse projection file to extract relevant data.

        Args:
            structure_file (str): The filename of the structure input file.
            proj_file (str): The filename pattern for the projection output files.
        """
        # Get the elements from the structure file
        struc = PWInput.from_file(structure_file).structure
        element = []
        for elm in struc.elements:
            element.append(str(elm))
        # To Do: For spin pol
        #if spin == 2:
        #projfile = glob.glob(proj_file)
        # Find the projection file matching the pattern
        projfile = glob.glob(proj_file)[0]
        #if spin == 2:
        # Add outer loop
        #for proj in projfile:
        # Extract orbital informations for each element from the projection file
        for elm in element:
            os.system(f"grep {elm} {projfile} | sed '1d' > proj-{elm}.dat")
        # Read the proj-{elm}.dat file and extract necessary information
        with open(f"{projfile}", "r") as read_projfile:
            self.projdata = read_projfile.readlines()
        _, self.nkpoint, self.nband = self.projdata[12].split()
        self.nkpoint = int(self.nkpoint)
        self.nband = int(self.nband)
        self.nrow = int(self.nkpoint * self.nband)
        self.dict_ = {}
        self.dict_1 = {}
        self.dict_2 = {}
        # For each element, store data, line with orbital names in short form, and full form
        # as lists
        for elm in element:
            sublist = []
            typelist = []
            namelist = []
            with open(f"proj-{elm}.dat", "r") as read_elm:
                lines = read_elm.readlines()
                for line in lines:
                    line = line.rstrip()
                    namelist.append(line)
                    linex = " ".join([line.split()[k] for k in [2, 3, 6]])
                    typelist.append(linex)
                    for i in range(len(self.projdata)):
                        if line in self.projdata[i]:
                            sublist.append(i)
            self.dict_[elm] = sublist
            self.dict_1[elm] = typelist
            self.dict_2[elm] = namelist

    def write_proj_elm_orbi(self):
        """
        Write projection data for each element/orbital combination.
        """
        # Iterate over each element key
        for key in self.dict_.keys():
            value = self.dict_[key]
            nameval = self.dict_1[key]
            namecount = Counter(nameval)
            lenv = len(value)
            dict_value = {}
            # Iterate over each index and value pair
            for i, idx in enumerate(value):
                ncount = namecount[nameval[i]]
                subname = nameval[i].split()
                # Define the orbital name
                if "S" in nameval[i]:
                    name = subname[1]
                elif "P 1" in nameval[i]:
                    name = subname[1] + "z"
                elif "P 2" in nameval[i]:
                    name = subname[1] + "x"
                elif "P 3" in nameval[i]:
                    name = subname[1] + "y"
                elif "D 1" in nameval[i]:
                    name = subname[1] + "z2"
                elif "D 2" in nameval[i]:
                    name = subname[1] + "xz"
                elif "D 3" in nameval[i]:
                    name = subname[1] + "yz"
                elif "D 4" in nameval[i]:
                    name = subname[1] + "x2-y2"
                elif "D 5" in nameval[i]:
                    name = subname[1] + "xy"
                else:
                    continue
                name = name.lower()
                # extract subdata from projection data
                subdata = self.projdata[idx + 1:idx + self.nrow + 1]
                subdata = self.convert_to_numpy(subdata)
                # Update the subdata if it has more than one occurence for a given element type.
                if name in dict_value:
                    dict_value[name] += subdata
                    #dict_value[name] += subdata / ncount
                else:
                    dict_value[name] = subdata
                    #dict_value[name] = subdata / ncount
                # Similar to dict_1, but with new defined orbital names.
                if key in self.proj_dict:
                    if name not in self.proj_dict[key]:
                        self.proj_dict[key].append(name)
                else:
                    self.proj_dict[key] = [name]
            # Similar to dict_, but with summed over similar atom type.
            self.final_dict[key] = dict_value

    def print_info(self):
        """
        Printing available options
        """
        self.write_proj_elm_orbi()
        print(f"Other Available options: {self.proj_dict}\n")

    def extract_data(self,elm,name):
        """
        Extract projection data for a specific element and orbital and save it to a file.

        Args:
            elm (str): Element for which projection data is extracted.
            name (str): Name of the orbital for which projection data is extracted.
        """
        # Write projection data for each element/orbital combination
        self.write_proj_elm_orbi()
        # Retrieve projection array for the specified element and orbital
        proj_array = self.final_dict[elm][name]
        # Save projection array to a file with appropriate formatting
        np.savetxt(f"{elm}-{name}.dat", proj_array, fmt=['%d', '%d', '%.4f'])
        # Convert to GNU plot format
        self.convert_gnu(f"{elm}-{name}.dat")

    def extract_data_2(self,elm,name):
        """
        Extract projection data for a specific element and orbital and save it to a file,
        while summing contributions from multiple suborbitals if applicable.
        Here, orbitals are "p", "d" instead of suborbitals "px", "py", "pz", etc.

        Args:
            elm (str): Element for which projection data is extracted.
            name (str): Name of the orbital or orbital type for which projection data is extracted.
        """
        # Determine the orbital type and list of orbitals to extract
        norb,lorb = name[:-1],name[-1]
        if lorb == 'p':
            orb_list = ['px', 'py', 'pz']
        elif lorb == 'd':
            orb_list = ['dxy', 'dxz', 'dyz', 'dz2', 'x2-y2']
        else:
            orb_list = [name]
        # Extract projection data for each orbital in the list
        for orbi in orb_list:
            self.extract_data(elm, norb+orbi)
        # Combine data from multiple suborbitals if applicable
        files = glob.glob(f"{elm}-{name}*.dat")
        for i, file_ in enumerate(files):
            self.ref_data = file_
            if i == 0:
                data = np.loadtxt(file_)
            else:
                newdata = np.loadtxt(file_)
                data[:, 2] += newdata[:, 2]
        np.savetxt(f"{elm}-{name}.dat", data, fmt=['%d', '%d', '%.4f'])
        self.convert_gnu(f"{elm}-{name}.dat")


    def combine_element(self,elm):
        """
        Combine projection data for a specific element and save it to a file.

        Args:
            elm (str): Element for which projection data is combined.
        """
        # Find files matching the element pattern
        files = glob.glob(f"{elm}-*.dat")
        # Iterate over files and combine data
        for i, file_ in enumerate(files):
            self.ref_data = file_
            if i == 0:
                data = np.loadtxt(file_)
            else:
                newdata = np.loadtxt(file_)
                data[:, 2] += newdata[:, 2]
        np.savetxt(f"{elm}.dat", data, fmt=['%d', '%d', '%.4f'])
        self.convert_gnu(f"{elm}.dat")

    def combine_orbital(self,orb):
        """
        Combine projection data for a specific orbital and save it to a file.

        Args:
            orb (str): Orbital type for which projection data is combined.
        """
        # Find files matching the orbital pattern
        files = glob.glob(f"*-{orb}.dat")
        # Iterate over files and combine data
        for i, file_ in enumerate(files):
            self.ref_data = file_
            if i == 0:
                data = np.loadtxt(file_)
            else:
                newdata = np.loadtxt(file_)
                data[:, 2] += newdata[:, 2]
        np.savetxt(f"{orb}.dat", data, fmt=['%d', '%d', '%.4f'])
        self.convert_gnu(f"{orb}.dat")

    def convert_to_numpy(self, input_data):
        """
        Convert input data to a NumPy array.

        Args:
            input_data (list): List of strings representing data.

        Returns:
            numpy.ndarray: NumPy array containing the converted data.
        """
        # Remove newline characters and convert each string to a list of floats
        data_list = [list(map(float, line.strip().split())) for line in input_data]
        # Convert the list of lists to a NumPy array
        np_array = np.array(data_list)
        return np_array

def main():
    """
    Main function to execute the program
    """
    # Loading configuration file and compound name
    input_data = config()
    comp = sys.argv[1]
    # Looking for structure file
    if os.path.isfile('scf.in'):
        filename = '../relax/scf.in'
    elif os.path.isfile('POSCAR'):
        filename = 'POSCAR'
    else:
        print("No scf.in and POSCAR exists\n")
    # Reading parameters from config.json file
    plot_info = input_data["plot"]["bandproj"]
    proj_type = plot_info["proj_type"]
    proj = plot_info["proj"]
    color = plot_info["colormap"]
    dft = input_data["download"]["inp"]["calc"]
    # If DFT is "VASP"
    if dft in ("VASP", "vasp"):
        # Initialize PROCARProcessor instance
        procar_processor = PROCARProcessor("PROCAR","band.dat")
        # Clean if there exists projection files
        procar_processor.clean()
        # Print available options
        print(f"projection type: {proj_type}, projections: {proj}")
        procar_processor.print_info()
        if proj_type == "element":
            # For element mode, loop over list of elements
            for elm in proj:
                # Write projection files for each elements
                procar_processor.write_band_element(elm)
        elif proj_type == "orbital":
            # Similar for orbital mode
            for orb in proj:
                procar_processor.write_band_orbital(orb)
        elif proj_type == "element-orbital":
            # For element-orbital mode, dictionary of key, value pair is used.
            # Keys are elements, and values are either single orbital ("p" or "px")
            # or list of orbitals (could be ["p", "d"] or ["px", "py", ...])
            # Doesn't work for f-orbitals for now.
            for elm_orb in proj:
                elm = elm_orb
                orb = proj[elm]
                if isinstance(orb,list):
                    for orbi in orb:
                        if orbi in ['p', 'd']:
                            procar_processor.extract_data_2(elm,orbi)
                        else:
                            procar_processor.extract_data(elm,orbi)
                else:
                    if orb in ['p', 'd']:
                        procar_processor.extract_data_2(elm,orb)
                    else:
                        procar_processor.extract_data(elm,orb)
    elif dft in ("QE", "qe"):
        # For spin polarized we may need to loop over up and down spin channel.
        # Initialize DataProcessor instance for QE proj.out* file.
        dp = DataProcessor(structure_file=filename)
        # Print available options
        print(f"projection type: {proj_type}, projections: {proj}")
        dp.print_info()
        # Clean files if already present
        clean_files = glob.glob("*.dat")
        for file_ in clean_files:
            if file_ not in f"{comp}.dat":
                dp.clean(file_)
        clean_files = glob.glob("*.dat.gnu")
        for file_ in clean_files:
            if file_ not in f"{comp}.dat.gnu":
                dp.clean(file_)
        # First creates all keys/values pair data
        # and combine according to plot_type.
        if proj_type == "element":
            # Creates projection files based on elements
            for elm in proj:
                proj_list = dp.proj_dict[elm]
                for projection in proj_list:
                    dp.extract_data(elm,projection)
                dp.combine_element(elm)
                for projection in proj_list:
                    os.system(f"rm {elm}-{projection}.dat.gnu")
        elif proj_type == "orbital":
            # Creates projection files based on orbitals
            for elm in dp.proj_dict:
                for projection in proj:
                    if projection in dp.proj_dict[elm]:
                        dp.extract_data(elm,projection)
            for orb in proj:
                dp.combine_orbital(orb)
            for elm in dp.proj_dict:
                for projection in proj:
                    if projection in dp.proj_dict[elm]:
                        os.system(f"rm {elm}-{projection}.dat.gnu")
        elif proj_type == "element-orbital":
            # Creates projection files based on element-orbital combination
            # For element-orbital mode, dictionary of key, value pair is used.
            # Keys are elements, and values are either single orbital ("2p" or "2px")
            # or list of orbitals (could be ["2p", "3d"] or ["2px", "2py", ...])
            # Doesn't work for f-orbitals for now.
            # Unlike VASP, here one has to define 2p, 3d instead of p, d.
            for elm_orb in proj:
                elm = elm_orb
                orb = proj[elm]
                if isinstance(orb,list):
                    for orbi in orb:
                        try:
                            dp.extract_data(elm,orbi)
                        except:
                            dp.extract_data_2(elm,orbi)
                else:
                    try:
                        dp.extract_data(elm,orb)
                    except:
                        dp.extract_data_2(elm,orb)
    else:
        print("Either QE/qe or VASP/vasp allowed\n")
    # Plot band structure with projections overlay.
    #color_maps = ['Reds', 'Blues', 'Greens', 'Purples', 'Greys', 'Oranges',
    #                  'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
    #                  'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
    if dft in ("VASP", "vasp"):
        proj_files = glob.glob("*.dat.gnu")
        print("-----------------------------------\n")
        print(f"Band file: {comp}.dat.gnu")
        print("\n")
        print("Processing PROCAR file\n")
        print("-----------------------------------\n")
        for _,proj_file in enumerate(proj_files):
            if proj_file not in f"{comp}.dat.gnu":
                print(proj_file)
                plot("band",filename,comp,proj_file,procar_processor.nkpoint,colormap=color)
                proj_name = proj_file.split(".dat.gnu")[0]
                os.system(f"mv {comp}-band.pdf {comp}-{proj_name}-band.pdf")
    elif dft in ("QE", "qe"):
        proj_files = glob.glob("*.dat.gnu")
        print("-----------------------------------\n")
        print(f"Band file: {comp}.dat.gnu")
        print("\n")
        print("Processing proj.out* file\n")
        print("-----------------------------------\n")
        for _,proj_file in enumerate(proj_files):
            if proj_file not in f"{comp}.dat.gnu":
                print(proj_file)
                plot("band",filename,comp,proj_file,colormap=color)
                proj_name = proj_file.split(".dat.gnu")[0]
                os.system(f"mv {comp}-band.pdf {comp}-{proj_name}-band.pdf")
    else:
        print("Only VASP/vasp or QE/qe options available\n")

if __name__ == "__main__":
    main()
