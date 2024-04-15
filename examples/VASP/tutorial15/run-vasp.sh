#!/bin/bash

#SBATCH -x dense001 --partition=dense --ntasks=1 --cpus-per-task=48 --time=1-0

# record start date and host


mpirun -np 24 vasp_std
