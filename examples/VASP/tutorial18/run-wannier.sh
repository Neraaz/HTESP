#!/bin/bash

#SBATCH --partition=dense --ntasks=1 --cpus-per-task=48 --time=1-0

# record start date and host


mpirun -np 48 vasp_std
wannier90.x wannier90
