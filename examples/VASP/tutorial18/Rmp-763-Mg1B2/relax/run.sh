#!/bin/bash

#SBATCH --partition=dense --ntasks=1 --cpus-per-task=48 --time=1-0

# record start date and host


mpirun -np 24 /home/nnepal/bin/vasp.6.3.0-wannier/bin/vasp_std
