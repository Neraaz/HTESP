#!/bin/bash

#SBATCH --partition=dense --ntasks=1 --cpus-per-task=12 --time=0-1

# record start date and host
date
hostname


source /shared/intel/bin/compilervars.sh intel64
source /shared/intel/impi/2019.0.117/intel64/bin/mpivars.sh
#mpirun -np 24 plotband.x < plotband.in
#mpirun -np 12 matdyn.x < matdyn.in.dos > matdyn.out.dos
mpirun -np 12 matdyn.x < matdyn-dos.in > matdyn-dos.out
# record end date
date
