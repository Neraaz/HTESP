# coding: utf-8
import os
import sys
import numpy as np
import itertools
from pymatgen.core import structure
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.vasp.sets import MPRelaxSet
from pymatgen.io import pwscf
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
from cif_to_gsinput import pos_to_kpt
from write_potcar import poscar2potcar
from htepc import MpConnect
class layer_rotation:
    def __init__(self):
        self.filename = None
        self.data = None
        self.cm_frame_xy = None
        self.zcoord_to_rotate = None

    def read_structure(self,filename):
        self.filename = filename
        self.data = structure.Structure.from_file(self.filename)

    def find_zcord_to_rotate(self,layer_to_rotate):
        zcoord = self.data.cart_coords[:,2]
        zcoord=np.round(zcoord,4)
        layer_coord = np.unique(zcoord)
        #layer_to_rotate = np.array([2,4,6,8])
        self.zcoord_to_rotate = layer_coord[layer_to_rotate - 1]

    def find_cm_frame(self):
        self.cm_frame_xy = np.matmul(self.data.lattice.matrix.T,np.array([0.5,0.5,0]))

    def layer_rotate(self,angle,direction_to_rotate):
        for i,site in enumerate(self.data.sites):
            idx = np.where(self.zcoord_to_rotate == round(site.coords[2],4))
            if idx[0].shape[0] > 0:
                #cm_frame = np.array([0,0,zcoord_to_rotate[idx[0][0]]])
                cm_frame = self.cm_frame_xy + [0,0,self.zcoord_to_rotate[idx[0][0]]] 
                direction_index = idx[0][0]
                coords = np.array(site.coords)
                site.coords = rotation(coords,angle,direction_to_rotate[direction_index],cm_frame) + cm_frame
        self.data.to("rotate_poscar.vasp",fmt='poscar')
def rotation(position,angle,direction,cm_frame):
    radian = angle * np.pi/180.0
    position = position - cm_frame
    if direction == "-":
        rotation_matrix = np.array([[np.cos(radian),-1.0*np.sin(radian),0],[np.sin(radian),np.cos(radian),0],[0,0,1]])
    elif direction == "+":
        rotation_matrix = np.array([[np.cos(radian),np.sin(radian),0],[-1.0*np.sin(radian),np.cos(radian),0],[0,0,1]])
    else:
        rotation_matrix = np.array([[1,0,0],[0,1,0],[0,0,1]])
    #print("rotation matrix: ",rotation_matrix)
    new_pos = np.matmul(rotation_matrix,position)
    return new_pos
#center = np.zeros(4,3)
if __name__ == "__main__":
    mpid = sys.argv[1]
    compound = sys.argv[2]
    orig_prefix = compound
    layer = np.array([2,4,6,8])
    #direction_to_rotate = ['+','-','-','+']
    angle = 10
    dft = "VASP"
    filename = "R{}-{}/relax/POSCAR".format(mpid,compound)
    obj = layer_rotation()
    obj.read_structure(filename)
    obj.find_zcord_to_rotate(layer)
    obj.find_cm_frame()
    characters = ['+', '-', '0']
    combinations = itertools.product(characters, repeat=4)
    list_struc = []
    for i,combination in enumerate(combinations):
        combination = list(combination)
        obj.layer_rotate(angle,combination)
        struc = structure.Structure.from_file("rotate_poscar.vasp")
        list_struc.append(struc)
        os.system("mv rotate_poscar.vasp rotate_poscar_{}.vasp".format(i+1))
    print("Initial structure combinations: {}".format(len(list_struc)))
    matcher = StructureMatcher()
    unique_structures = matcher.group_structures(list_struc)
    print("Final unique structures after applying StructureMatcher: {}".format(len(unique_structures)))
    if not os.path.isdir('scf_dir'):
        os.mkdir("scf_dir")
    if not os.path.isfile('mpid-rotate.in'):
        entry = 0
    else:
        with open('mpid-rotate.in', 'r') as mpid_read:
            lines = mpid_read.readlines()
        entry = len(lines)
    for i,struclist in enumerate(unique_structures):
        struc = struclist[0]
        struc.to("poscar{}.vasp".format(i+1),fmt='poscar')
        if dft in ('vasp', 'VASP'):
            obj.prefix = struc.composition.alphabetical_formula.replace(" ","")
            poscar = Poscar(structure=struc,comment=obj.prefix)
            if not os.path.isdir("R{}-{}-{}".format(mpid,i+1,obj.prefix)):
                os.mkdir("R{}-{}-{}".format(mpid,i+1,obj.prefix))
            if not os.path.isdir("R{}-{}-{}/relax".format(mpid,i+1,obj.prefix)):
                os.mkdir("R{}-{}-{}/relax".format(mpid,i+1,obj.prefix))
            pwd=os.getcwd()
            os.system("cp R{}-{}/relax/INCAR {}/R{}-{}-{}/relax/".format(mpid,orig_prefix,pwd,mpid,i+1,obj.prefix))
            poscar.write_file(filename="R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix))
            #os.system("cp R{}-{}/relax/KPOINTS {}/R{}-{}-{}/relax/".format(mpid,orig_prefix,pwd,mpid,i+1,obj.prefix))
            if os.path.isfile("download.py"):
                import download as d
            evenkpt = d.inp['evenkpt']
            if evenkpt:
                print("Even kpoint mesh is utilized\n")
                pos_to_kpt("R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix),0.025,True)
            else:
                pos_to_kpt("R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix),0.025)
            os.system("mv KPOINTS R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
            structure_file = structure.Structure.from_file("R{}-{}-{}/relax/POSCAR".format(mpid,i+1,obj.prefix))
            relax_set = MPRelaxSet(structure=structure_file)
            #relax_set.potcar.write_file("R{}-{}-{}/relax/POTCAR".format(mpid,i+1,obj.prefix))
            if os.path.isfile("vasp.in"):
                os.system("cp vasp.in R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
            if os.path.isfile("magmom.py"):
                os.system("cp magmom.py R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
            if os.path.isfile("pseudo.py"):
                os.system("cp pseudo.py R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
            os.chdir("R{}-{}-{}/relax/".format(mpid,i+1,obj.prefix))
            poscar2potcar()
            os.system("vasp_process.py POSCAR")
            os.chdir("../../")
        else:
            obj.structure = struc
            comp_list = []
            for composition in obj.structure.composition.elements:
                comp_list.append(composition.name)
            obj.comp_list = comp_list
            obj.getkpt()
            obj.getevenkpt()
            obj.maxecut_sssp_for_subs()
            obj.prefix = struc.composition.alphabetical_formula.replace(" ","")
            obj.setting_qeinput(pseudo_dir='../../pp')
            os.system("""mv scf-None.in""" + """ scf_dir/scf-{}-{}.in""".format(mpid,i+1))
        print(mpid,obj.prefix)
        with open("mpid-rotate.in", "a") as mpid_append:
            mpid_append.write("v{}".format(entry+1+i) + " " + mpid + "-{}".format(i+1) + " " + obj.prefix + "\n")
    if not os.path.isdir("structures"):
        os.mkdir("structures")
        os.system("mv *.vasp structures")
    else:
        os.system("mv *.vasp structures")
