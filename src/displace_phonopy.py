#!/usr/bin/env python
#"""Written by Niraj K. Nepal,Ph.D."""
"""Module to process band.yaml file from phonopy calculations to extract eigenvectors"""
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

    Example
    -------
    >>> displace = DisplacePhonopy("band.conf") #create an object
    >>> displace.parse() #parsing band.conf
    >>> displace.extract_mass() #Extracting mass and store in self.mass
    >>> displace.print_qpt() #print qpoints and its index
    0 [0.0, 0.0, 0.0]
    1 [0.1, 0.1, 0.1]
    ...
    >>> displace.print_mode(0) #print band energies for a qpoint
    1 0.0
    2 0.1
    ...
    >>> displace.negative_freq() #returns dictionary of qpoints, band index, and phonon frequency
    >>> #Easy to locate qpoint (iq) and band index (inu) for negative frequencies
    >>> displace.geteigenvec(0, 0) #get eigen vector corresponding to particular iq and inu.
    >>> displace.applydisplace() #Apply distortion to the structure, provided an undistorted
    >>> #structure exists as POSCAR file
    """
    def __init__(self,filename):
        self.filename = filename
        self.data = None
        self.eigen = None
        self.mass = None

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

    def mass_extract(self):
        """
        Returns list of masses of the ions
        """
        mass = []
        point_list = self.data['points']
        for point in point_list:
            mass.append(point['mass'])
        self.mass = mass
            
        

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
            displace_eigen = self.eigen[i]/np.sqrt(self.mass[i])
            site.coords += np.array([displace_eigen[0][0],displace_eigen[1][0],displace_eigen[2][0]])
        struc.to("displace_poscar.vasp",fmt="poscar")
