"""
Input_output script for berkeleyGW code, primariliy with Quantum Espresso (QE)
written by Niraj K. Nepal, Ph.D.
Email: tug11655@temple.edu
"""
import sys
import os
import numpy as np
from pymatgen.io.pwscf import PWInput

class bgw:
    def __init__(self,mpid,compound):
        self.mpid = mpid
        self.comp = compound

    def pw2bgw_in(self,kpoint,vxc,complx=True,wfng_flag=True,wfng_file="WFN",wfng_kgrid=True,out="pp_in"):
        nkpts = kpoint[0]
        kshift = kpoint[1]
        nk1_list = ["wfng_nk1","wfng_nk2","wfng_nk3"] 
        dk1_list = ["wfng_dk1","wfng_dk2","wfng_dk3"]
        os.system("""grep prefix R{}-{}/relax/scf.in > prefix.dat""".format(self.mpid,self.comp))
        with open("prefix.dat", "r") as read_prefix:
            prefix = read_prefix.readlines()[0].rsplit(",")[0].rsplit()[-1]
        with open(out,"w") as write_pp:
            write_pp.write("&input_pw2bgw\n")
            write_pp.write("prefix = {}".format(prefix) + "\n")
            write_pp.write("outdir = './'" + "\n")
            if complx:
                write_pp.write("real_or_complex = 2" + "\n")
            else:
                write_pp.write("real_or_complex = 1" + "\n")
            if wfng_flag:
                write_pp.write("wfng_flag = .true." + "\n")
            else:
                write_pp.write("wfng_flag = .false." + "\n")
            if wfng_kgrid:
                write_pp.write("wfng_kgrid = .true." + "\n")
            else:
                write_pp.write("wfng_kgrid = .false." + "\n")
            if complx:
                write_pp.write("wfng_file = 'wfn.cplx'" + "\n")
            else:
                write_pp.write("wfng_file = '{}.real'" + "\n")
            for i,_ in enumerate(nk1_list):
                write_pp.write("{} = {}".format(nk1_list[i],int(nkpts[i])) + "\n")
            for i,_ in enumerate(dk1_list):
                if wfng_file == "WFN_co":
                    write_pp.write("{} = {}".format(dk1_list[i],0) + "\n")
                else:
                    write_pp.write("{} = {}".format(dk1_list[i],kshift[i]/2.0) + "\n")
            if wfng_file == "WFN_co":
                write_pp.write("rhog_flag = .true.\n")
                write_pp.write("vxcg_flag = .false.\n")
                if complx:
                    write_pp.write("rhog_file = 'rho.cplx'\n")
                    write_pp.write("vxcg_file = 'vxc.cplx'\n")
                else:
                    write_pp.write("rhog_file = 'rho.real'\n")
                    write_pp.write("vxcg_file = 'vxc.real'\n")
                write_pp.write("vxc_flag = .true.\n")
                write_pp.write("vxc_file = 'vxc.dat'\n")
                write_pp.write("vxc_diag_nmin = {}\n".format(vxc[0][0]))
                write_pp.write("vxc_diag_nmax = {}\n".format(vxc[0][1]))
                write_pp.write("vxc_offdiag_nmin = {}\n".format(vxc[1][0]))
                write_pp.write("vxc_offdiag_nmax = {}\n".format(vxc[1][1]))
            write_pp.write("/" + "\n")

    def write_wfn_grid(self,kpoint,q_shift,fft_size,time_rev,out="WFN.in"):
        with open(out,"w") as wfn_write:
            nkpts = kpoint[0]
            kshift = kpoint[1]
            wfn_write.write(str(int(nkpts[0])) + " " + str(int(nkpts[1])))
            wfn_write.write(" "+ str(int(nkpts[2])) + "\n")
            wfn_write.write(str(kshift[0]/2.0) + " " + str(kshift[1]/2.0))
            wfn_write.write(" "+ str(kshift[2]/2.0) + "\n")
            if len(q_shift) == 3:
                wfn_write.write(str(q_shift[0]) + " " + str(q_shift[1]) + " " + str(q_shift[2]) + "\n")
            wfn_write.write("\n")
            #os.system("""sed -n '/crystal axes:/,/reciprocal axes:/p' R{}-{}/relax/scf.out | sed '$d' | sed '$d' | sed '1d' > lattice.txt""".format(self.mpid,self.comp))
            os.system("""grep -A 3 'crystal axes:' R{}-{}/relax/scf.out | tail -n 3 > lattice.txt""".format(self.mpid,self.comp))
            
            #os.system("""sed -n '/crystal axes:/,/reciprocal axes:/p' R{}-{}/relax/scf.out | sed '$d' | sed '$d' | sed '1d' | awk '{print $4 " " $5 " " $6}' > lattice.txt""".format(self.mpid,self.comp))
            lattice = np.genfromtxt("lattice.txt")[:,[3,4,5]]
            for lat in lattice:
                wfn_write.write(str(lat[0]) + " " + str(lat[1]) + " " + str(lat[2]) + "\n")
            try:
                os.system("""grep nat R{}-{}/relax/scf.in | awk '{print $3}' > nat.txt""".format(self.mpid,self.comp))
            except:
                os.system("""grep 'number of atoms/cell' R{}-{}/relax/scf.out | tail -n 1 > nat.txt""".format(self.mpid,self.comp))
            nat = int(np.genfromtxt("nat.txt")[-1])
            wfn_write.write(str(int(nat)) + "\n")
            os.system("""grep -A {} 'positions (alat units)' R{}-{}/relax/scf.out | tail -n {} > position.txt""".format(nat,self.mpid,self.comp,nat))
            data = PWInput.from_file("R{}-{}/relax/scf.in".format(self.mpid,self.comp)).structure
            new_dict = {}
            comp_list = data.composition.elements
            for i,elm in enumerate(comp_list):
                new_dict[elm.name] = i+1
            with open("position.txt", "r") as read_pos:
                positions = read_pos.readlines()
            for line1 in positions:
                line1 = line1.rsplit()
                line = [line1[i] for i in [1,6,7,8]]
                #line = line.split("\n")[0].split(" ")
                elm_type = new_dict[line[0]]
                wfn_write.write(str(elm_type) + " " + str(line[1]) + " " + str(line[2]) + " " + str(line[3]) + "\n")
            wfn_write.write(str(fft_size[0]) + " " + str(fft_size[1]) + " " + str(fft_size[2]) + "\n")
            wfn_write.write(str(time_rev) + "\n")
            #.false.                ! use time-reversal symmetry. Set to false for BerkeleyGW

if __name__ == "__main__":
    MPID = sys.argv[1]
    COMP = sys.argv[2]
    CALC_TYPE = sys.argv[3]
    obj = bgw(MPID,COMP)
    if os.path.isfile("kpoint.in"):
        kpoint = np.loadtxt("kpoint.in")
        if kpoint.ndim == 1:
            kpoint = kpoint.reshape(2,3)
    else:
        os.system("""sed -n '/K_POINTS automatic/,/CELL_PARAMETERS angstrom/p' R{}-{}/relax/scf.in | sed '$d' | sed '1d' > kpoint.in""".format(MPID,COMP))
        kpoint = np.loadtxt("kpoint.in")
        kpoint = kpoint.reshape(2,3)
    if os.path.isfile("bgw_input.py"):
        import bgw_input as bgw_in
        fft_size = bgw_in.fft
        dqshift = bgw_in.qshift['qshift']
        qshift_i = bgw_in.qshift['index']
        random_kshift = bgw_in.random_kshift
        time_rev = bgw_in.time_rev
        nband = bgw_in.nband
        vxc_diag = bgw_in.vxc_diag
        vxc_offdiag = bgw_in.vxc_offdiag
    else:
        print("bgw_input.py not found, one is created with default values\n")
        print("Please edit the file according to your need\n")
        fft_size = [24,24,24]
        dqshift = 0.001
        qshift_i = 2
        time_rev = ".false."
        random_kshift = []
        nband = {'nbnd':2,'cbnd':1,'vbnd':1}
        vxc_diag = [1,2]
        vxc_offdiag=[0,0]
        with open("bgw_input.py", "w") as write_bgw:
            write_bgw.write("fft={}".format(fft_size) + "\n")
            write_bgw.write("random_kshift={}".format(random_kshift) + "\n")
            write_bgw.write("qshift={}".format({'qshift':qshift,'index':qshift_i}) + "\n")
            write_bgw.write("time_rev={}".format(time_rev) + "\n")
            write_bgw.write("nband={}".format(nband) + "\n")
            write_bgw.write("vxc_diag={}".format(vxc_diag) + "\n")
            write_bgw.write("vxc_offdiag={}".format(vxc_offdiag) + "\n")
    if CALC_TYPE == "kgrid":
        with open("kpt.in", "w") as write_kpt:
            write_kpt.write("K_POINTS automatic\n")
            offset = np.ceil(kpoint[1])
            write_kpt.write(str(int(kpoint[0][0])) + " " + str(int(kpoint[0][1])))
            write_kpt.write(" " + str(int(kpoint[0][2])) + " " + str(int(offset[0])))
            write_kpt.write(" " + str(int(offset[1])) + " " + str(int(offset[2])) + "\n")
        qshift_1 = [0.,0.,0.]
        qshift_2 = [0.,0.,0.]
        qshift_2[qshift_i] = qshift_1[qshift_i] + dqshift
        obj.write_wfn_grid(kpoint,q_shift=qshift_2,fft_size=fft_size,time_rev=time_rev,out="WFNq.in")
        obj.write_wfn_grid(kpoint,q_shift=qshift_1,fft_size=fft_size,time_rev=time_rev,out="WFN.in")
        os.system("""sed '2d' WFN.in | sed '2 i 0 0 0' > WFN_co.in""")
        kpoint[0] = kpoint[0]*2
        if len(random_kshift) > 0:
            kpoint[1] = random_kshift
        obj.write_wfn_grid(kpoint,q_shift=qshift_1,fft_size=fft_size,time_rev=time_rev,out="WFN_fi.in")
        kpoint[1][2] += dqshift
        obj.write_wfn_grid(kpoint,q_shift=qshift_1,fft_size=fft_size,time_rev=time_rev,out="WFNq_fi.in")
        if not os.path.isdir("R{}-{}/bgw".format(MPID,COMP)):
            os.mkdir("R{}-{}/bgw".format(MPID,COMP))
        if not os.path.isdir("R{}-{}/bgw/00-kgrid".format(MPID,COMP)):
            os.mkdir("R{}-{}/bgw/00-kgrid".format(MPID,COMP))
        os.system("mv kpt.in WFN.in WFNq.in WFN_co.in WFN_fi.in WFNq_fi.in R{}-{}/bgw/00-kgrid/".format(MPID,COMP))
        os.chdir("R{}-{}/bgw/00-kgrid".format(MPID,COMP))
        WFN_LIST = ["WFN.in","WFNq.in","WFN_co.in","WFN_fi.in","WFNq_fi.in"]
        for WFN in WFN_LIST:
            os.system("kgrid.x ./{} ./{}.out ./{}.log".format(WFN,WFN,WFN))
        os.chdir("../../../")
    elif CALC_TYPE == "input":
        obj.pw2bgw_in(kpoint,vxc=[vxc_diag,vxc_offdiag])
        #obj.pw2bgw_in(kpoint,vxc=[vxc_diag,vxc_offdiag],wfng_file="WFN_co")
        bandval = nband.values()
        with open("bgw_band.in", "w") as write_gwband:
            for bands in bandval:
                write_gwband.write("{} ".format(int(bands)))
            write_gwband.write("\n")
        
    elif CALC_TYPE == "scf":
        pass
    elif CALC_TYPE == "scfall":
        pass
    else:
        print("Please use provided options\n")
