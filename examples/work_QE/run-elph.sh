#!/bin/bash

#SBATCH --partition=dense --ntasks=2 --cpus-per-task=12 --time=2-0

# record start date and host
date
hostname

source /shared/intel/bin/compilervars.sh intel64
source /shared/intel/impi/2019.0.117/intel64/bin/mpivars.sh

#mpirun -np 24 plotband.x < plotband.in
#mpirun -np 12 matdyn.x < matdyn.in.dos > matdyn.out.dos
mpirun -np 24 ph.x < elph.in > elph.out
# record end date
date
