"""Written by Niraj K. Nepal,Ph.D."""
# coding: utf-8
import yaml
import numpy as np
from pymatgen.core import structure

class DisplacePhonopy:
    """
    A class to handle displacement based on VASP+PHONOPY band.yaml file.
    To obtain eigenvectors on band.yaml, switch on eigenvector tag to true in band.conf file.

    Parameters
    ----------
    filename : str
        The name of the YAML file containing the displacement data, usually band.conf.

    Attributes
    ----------
    filename : str
        The name of the YAML file containing the displacement data.
    data : dict
        The parsed data from the YAML file.
    eigen : list
        List containing eigenvalues.
    """
    def __init__(self,filename):
        self.filename = filename
        self.data = None
        self.eigen = None

    def parse(self):
        """
        Parses the YAML file specified by 'filename'.

        Returns
        -------
        None
        """
        with open(self.filename, "r") as stream:
            try:
                self.data = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)

    def print_qpt(self):
        """
        Prints the q-point positions from the parsed data.

        Returns
        -------
        None
        """
        for i,l in enumerate(self.data['phonon']):
            print(i,l['q-position'])

    def print_mode(self,iqpt):
        """
        Prints the energy modes for a given q-point index.

        Parameters
        ----------
        iqpt : int
            Index of the q-point for which energy modes are to be printed.

        Returns
        -------
        None
        """
        energy = self.data['phonon'][iqpt]['band']
        nmode = len(energy)
        for i,en in enumerate(energy):
            print(i+1,en['frequency'])

    def negative_freq(self):
        """
        Identifies q-points with negative frequencies and returns them as a dictionary.

        Constructs a dictionary where each key represents a q-point and its associated
        value is a list containing the q-point index, band index, and frequency with
        negative values. The q-point index and band index are combined in the key.

        Returns
        -------
        dict
            Dictionary with q-point as key and a list containing q-point index, band index,
            and negative frequency as values.
        """
        dict_imag = {}
        print("Returns dictionary with qpoint as key, and qpoint index, bandindex and frequency as values combined in list\n")
        for i,en in enumerate(self.data['phonon']):
            data_q = self.data['phonon'][i]['band']
            for j, freq in enumerate(data_q):
                if freq['frequency'] < 0:
                    dict_imag['{}-{}'.format(en['q-position'],i)] = [j,freq['frequency']]
        return dict_imag

    def geteigenvec(self,iqpt,imode):
        """
        Retrieves the eigenvectors for a specified q-point index and mode index.

        Parameters
        ----------
        iqpt : int
            Index of the q-point.
        imode : int
            Index of the mode.

        Returns
        -------
        None
        """
        self.eigen = self.data['phonon'][iqpt]['band'][imode]['eigenvector']


    def applydisplace(self):
        """
        Applies displacement to the structure based on the eigenvectors.

        Retrieves the structure from the "POSCAR" file using pymatgen and applies
        displacement to each site based on the corresponding eigenvector values.
        The modified structure is then written to a new VASP POSCAR file named
        "displace_poscar.vasp".

        Returns
        -------
        None
        """
        struc = structure.Structure.from_file("POSCAR")
        for i,site in enumerate(struc.sites):
            displace_eigen = self.eigen[i]
            site.coords += np.array([displace_eigen[0][0],displace_eigen[1][0],displace_eigen[2][0]])
        struc.to("displace_poscar.vasp",fmt="poscar")
