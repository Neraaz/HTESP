#!/bin/bash

#SBATCH --partition=dense --ntasks=1 --cpus-per-task=48 --time=1-0

# record start date and host


mpirun -np 24 pw.x -in scf.in > scf.out
mpirun -np 24 pw.x -in nscf.in > nscf.out
wannier90.x -pp ex
mpirun -np 24 pw2wannier90.x -in pw2wan.in > pw2wan.out
wannier90.x ex
