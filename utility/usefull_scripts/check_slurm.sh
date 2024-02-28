#!/bin/bash


ncalc=$(squeue | grep "nnepal" | awk '{ print $1 }')

#ncalc=$(echo $ncalc + 1 | bc)

for id in $ncalc; do
 find * -name slurm* | grep $id
done
