#!/usr/bin/env python
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

    Attributes:
    - nkpoint (int): Number of k-points.
    - nband (int): Number of bands.
    - nion (int): Number of ions.
    - data (numpy.ndarray): Processed PROCAR data.

    Methods:
    - write_band_orbital(output_file): Writes band orbital data to the specified file.
    - write_band_element(output_file): Writes band element data to the specified file.
    """

    def __init__(self, procar_file, band_file):
        """
        Initialize PROCARProcessor with the provided PROCAR file.
        """
        self.procar_file = procar_file
        self.band_file = band_file
        self.nkpoint, self.nband, self.nion, self.ncol, self.col = self._parse_procar_header(self.band_file)
        self.data, self.nspin = self._parse_procar_data()
        self.elements = self._parse_structure()
        self.nkpoint_extra = None

    def _parse_procar_header(self,band_file):
        """
        Parse header information from the PROCAR file.
        """
        with open(self.procar_file, "r") as f:
            procar = f.readlines()
        ncol = len(procar[7].split())
        col = procar[7].split()[1:-1]
        orig_band = np.loadtxt(self.band_file)
        nkpoint = int(orig_band.shape[0])
        nkpoint_total = int(procar[1].split(":")[1].split()[0])
        self.nkpoint_extra = nkpoint_total - nkpoint
        nband = int(procar[1].split(":")[2].split()[0])
        nion = int(procar[1].split(":")[3].split()[0])
        self.skip = int(nband * self.nkpoint_extra)
        return nkpoint, nband, nion, ncol, col
    def _parse_structure(self):
        struc = structure.Structure.from_file("POSCAR")
        elements = []
        for elm in struc:
            elements.append(elm.species_string)
        return elements
    def print_info(self):
        unique_elements = list(set(self.elements))
        orbital_names = self.col + ["p", "d"]
        print(f"Available options: elements {unique_elements}, orbitals {orbital_names}")
    def clean(self):
        for elm in self.elements:
            clean_files = glob.glob(f"{elm}-*.dat.gnu")
            if len(clean_files) > 0:
                for file_ in clean_files:
                    os.system(f"rm {file_}")

    def _parse_procar_data(self):
        """
        Parse data from the PROCAR file.
        """
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
        with open(self.procar_file, "r") as f:
            procar = f.readlines()
        ncount = 0
        with open("procar_process.txt", "w") as procar_write:
            if self.nspin in (1, 4):
                for i, line in enumerate(procar):
                    if 'ion      s' in line:
                        if ncount >= self.skip:
                            for j in range((self.nion+1)*self.nspin):
                                procar_write.write(procar[i+j+1])
                        ncount += 1
            else:
                len_s = int(len(procar[1:])/2)
                procar_up = procar[1:len_s+1]
                procar_dn = procar[len_s+1:]
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
                        
        data = np.genfromtxt("procar_process.txt")
        if self.nspin == 2:
            data = data.reshape(self.nspin, self.nkpoint, self.nband, self.nion+1, self.ncol)
        else:
            data = data.reshape(self.nkpoint, self.nband, self.nspin, self.nion+1, self.ncol)
        # For noncollinear
        if self.nspin == 4:
            data = data[:,:,0,:,:] 
        return data, self.nspin
        #return data.reshape(self.nkpoint, self.nband, self.nspin, self.nion+1, self.ncol)

    def extract_data(self,elm,orb):
        # To Do: Handle for ISPIN = 2 only
        indices = find_all_indices(self.col,orb)
        if self.nspin == 2: 
            subdata = self.data[:,:,:,find_all_indices(self.elements,elm),indices]
        else:
            subdata = self.data[:,:,find_all_indices(self.elements,elm),indices]
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
        # To Do: Handle for ISPIN = 2 only
        if orb == 'p':
            orb_list = ['px', 'py', 'pz']
        elif orb == 'd':
            orb_list = ['dxy', 'dyz', 'dxz', 'dz2', 'x2-y2']
        else:
            print("Supports p or d\n")
        for orbi in orb_list:
            self.extract_data(elm,orbi)
        data0 = np.loadtxt(f"{elm}-{orb_list[0]}.dat.gnu")
        data = np.zeros_like(data0)
        data[:,0] = data0[:,0]
        for orbi in orb_list:
            data_i = np.loadtxt(f"{elm}-{orbi}.dat.gnu")
            data[:,1] += data_i[:,1]
            if self.nspin == 2:
                data[:,2] += data_i[:,2]
        np.savetxt(f"{elm}-{orb}.dat.gnu",data) 

    def write_band_orbital(self, orb):
        """
        Write band orbital data to the specified file.

        Parameters:
        - output_file (str): Path to the output file.
        """
        with open(f"{orb}.dat.gnu", "w") as write_procar:
            for i in range(self.nband):
                if self.nspin == 2:
                    procar_iband = self.data[:, :, i, self.nion, 1:self.ncol]
                    procar_iband = procar_iband[:, :,find_all_indices(self.col,orb)]
                else:
                    procar_iband = self.data[:, i, self.nion, 1:self.ncol]
                    procar_iband = procar_iband[:, find_all_indices(self.col,orb)]
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
        Write band element data to the specified file.

        Parameters:
        - output_file (str): Path to the output file.
        """
        with open(f"{elm}.dat.gnu", "w") as write_procar:
            for i in range(self.nband):
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
    def __init__(self, structure_file="scf.in", proj_file="proj.dat*"):
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
        data = np.loadtxt(filename)
        data = data[:,2]
        data = data.reshape(self.nkpoint, self.nband)
        with open(f"{filename}.gnu", "w") as write_gnu:
            for i in range(data.shape[1]):
                datak = data[:,i]
                for j in range(datak.shape[0]):
                    write_gnu.write(str(j) + " " + str(datak[j]) + "\n")
                write_gnu.write("\n")

    def combine_orbital(self, key, orb, if_gnu=False):
        self.write_proj_elm_orbi()
        for name in self.final_dict[key]:
            proj_array = self.final_dict[key][name]
            np.savetxt(f"{key}-{name}.dat",proj_array,fmt=['%d', '%d', '%.4f'])
            if if_gnu:
                self.convert_gnu(f"{key}-{name}.dat")
        sdata = np.zeros((self.nrow, 3))
        sfile = glob.glob(f"{key}-*{orb}*.dat")
        s0 = np.loadtxt(sfile[0])
        sdata[:,0] = s0[:,0]
        sdata[:,1] = s0[:,1]
        len_s = len(sfile)
        for i in range(len_s):
            idata = np.loadtxt(sfile[i])
            sdata[:,2] += idata[:,2] / len_s
        np.savetxt(f"{key}-{orb}.dat", sdata, fmt=['%d', '%d', '%.4f'])
        self.convert_gnu(f"{key}-{orb}.dat")

    def clean(self,onlydat = False):
        if onlydat:
            os.system("rm *.dat")
        else:
            os.system("rm *.dat")
            os.system("rm *.dat.gnu")

    def _parse_proj_file(self, structure_file, proj_file):
        struc = PWInput.from_file(structure_file).structure
        element = []
        for elm in struc.elements:
            element.append(str(elm))
        # To Do: For spin pol
        #if spin == 2:
        #projfile = glob.glob(proj_file)
        projfile = glob.glob(proj_file)[0]
        #if spin == 2:
        # Add outer loop 
        #for proj in projfile:
        for elm in element:
            os.system(f"grep {elm} {projfile} | sed '1d' > proj-{elm}.dat")
        with open(f"{projfile}", "r") as read_projfile:
            self.projdata = read_projfile.readlines()
        _, self.nkpoint, self.nband = self.projdata[12].split()
        self.nkpoint = int(self.nkpoint)
        self.nband = int(self.nband)
        self.nrow = int(self.nkpoint * self.nband)
        self.dict_ = {}
        self.dict_1 = {}
        self.dict_2 = {}
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
        for key in self.dict_.keys():
            value = self.dict_[key]
            nameval = self.dict_1[key]
            namecount = Counter(nameval)
            lenv = len(value)
            dict_value = {}
            for i, idx in enumerate(value):
                ncount = namecount[nameval[i]]
                subname = nameval[i].split()
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
                    name = subname[1] + "zx"
                elif "D 3" in nameval[i]:
                    name = subname[1] + "zy"
                elif "D 4" in nameval[i]:
                    name = subname[1] + "x2-y2"
                elif "D 5" in nameval[i]:
                    name = subname[1] + "xy"
                else:
                    continue
                name = name.lower()
                subdata = self.projdata[idx + 1:idx + self.nrow + 1]
                subdata = self.convert_to_numpy(subdata)
                if name in dict_value:
                    dict_value[name] += subdata
                    #dict_value[name] += subdata / ncount
                else:
                    dict_value[name] = subdata
                    #dict_value[name] = subdata / ncount
                if key in self.proj_dict:
                    if name not in self.proj_dict[key]:
                        self.proj_dict[key].append(name)
                else:
                    self.proj_dict[key] = [name]
            self.final_dict[key] = dict_value

    def print_info(self):
        self.write_proj_elm_orbi()
        print(f"Other Available options: {self.proj_dict}\n")
        
    def extract_data(self,elm,name):
        self.write_proj_elm_orbi()
        proj_array = self.final_dict[elm][name]
        np.savetxt(f"{elm}-{name}.dat", proj_array, fmt=['%d', '%d', '%.4f'])
        self.convert_gnu(f"{elm}-{name}.dat")

    def combine_element(self,elm):
        files = glob.glob(f"{elm}-*.dat")
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
        files = glob.glob(f"*-{orb}.dat")
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
        # Remove newline characters and convert each string to a list of floats
        data_list = [list(map(float, line.strip().split())) for line in input_data]
        # Convert the list of lists to a NumPy array
        np_array = np.array(data_list)
        return np_array

def main():
    input_data = config()
    comp = sys.argv[1]
    if os.path.isfile('scf.in'):
        filename = 'scf.in'
    elif os.path.isfile('POSCAR'):
        filename = 'POSCAR'
    else:
        print("No scf.in and POSCAR exists\n")
    plot_info = input_data["plot"]["bandproj"]
    proj_type = plot_info["proj_type"]
    proj = plot_info["proj"]
    dft = input_data["download"]["inp"]["calc"]
    # Command for plotting
    # plot("band",filename,comp,proj)
    # Copy adding orbital_projection
    # cp comp + '-band.pdf' to comp+proj_name+'-band.pdf'
    # proj_name.dat.gnu.split(".dat.gnu")[0]
    if dft in ("VASP", "vasp"):
        procar_processor = PROCARProcessor("PROCAR","band.dat")
        procar_processor.clean()
        print(f"projection type: {proj_type}, projections: {proj}")
        procar_processor.print_info()
        if proj_type == "element":
            for elm in proj:
                procar_processor.write_band_element(elm)
        elif proj_type == "orbital":
            for orb in proj:
                procar_processor.write_band_orbital(orb)
        elif proj_type == "element-orbital":
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
        dp = DataProcessor(structure_file=filename)
        print(f"projection type: {proj_type}, projections: {proj}")
        dp.print_info()
        dp.clean()
        if proj_type == "element":
            for elm in proj:
                proj_list = dp.proj_dict[elm]
                for projection in proj_list:
                    dp.extract_data(elm,projection)
                dp.combine_element(elm)
                for projection in proj_list:
                    os.system(f"rm {elm}-{projection}.dat.gnu")
        elif proj_type == "orbital":
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
            for elm_orb in proj:
                elm = elm_orb
                orb = proj[elm]
                dp.extract_data(elm,orb)
        dp.clean(onlydat=True)
    else:
        print("Either QE/qe or VASP/vasp allowed\n")
    if dft in ("VASP", "vasp"):
        proj_files = glob.glob("*.dat.gnu")
        print("-----------------------------------\n")
        print(f"Band file: {comp}.dat.gnu")
        print("\n")
        print("Processing PROCAR file\n")
        print("-----------------------------------\n")
        for proj_file in proj_files:
            if proj_file not in f"{comp}.dat.gnu":
                print(proj_file)
                plot("band",filename,comp,proj_file,procar_processor.nkpoint)
                proj_name = proj_file.split(".dat.gnu")[0]
                os.system(f"mv {comp}-band.pdf {comp}-{proj_name}-band.pdf")
    elif dft in ("QE", "qe"):
        proj_files = glob.glob("*.dat.gnu")
        print(f"Band file: {comp}.dat.gnu")
        for proj_file in proj_files:
            if proj_file not in f"{comp}.dat.gnu":
                print(proj_file)
                plot("band",filename,comp,proj_file)
                proj_name = proj_file.split(".dat.gnu")[0]
                os.system(f"mv {comp}-band.pdf {comp}-{proj_name}-band.pdf")
    else:
        print("Only VASP/vasp or QE/qe options available\n")

if __name__ == "__main__":
    main()
