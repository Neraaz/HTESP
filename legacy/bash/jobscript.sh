#!/bin/bash
SLEEP=$(echo "scale=2; 0.01" | bc)
vasprun() {
    A="$1"
    B="$2"
    if [[ -f ../../CALC_VISIBLE_WITH_ID ]]; then
     mv run.sh $A.sh
     sbatch $A.sh
    elif [[ -f ../../CALC_VISIBLE_WITH_NAME ]]; then
     mv run.sh $B.sh
     sbatch $B.sh
    elif [[ -f ../../CALC_VISIBLE_WITH_ID-NAME ]]; then
     mv run.sh $A-$B.sh
     sbatch $A-$B.sh
    else
     sbatch run.sh
    fi
}


qerun() {
    A="$1"
    B="$2"
    qerun="$3"
     if [[ -f ../../CALC_VISIBLE_WITH_ID ]]; then
       mv run-$qerun.sh $A-$qerun.sh
       sbatch $A-$qerun.sh
     elif [[ -f ../../CALC_VISIBLE_WITH_NAME ]]; then
       mv run-$qerun.sh $B-$qerun.sh
       sbatch $B-$qerun.sh
     elif [[ -f ../../CALC_VISIBLE_WITH_ID-NAME ]]; then
       mv run-$qerun.sh $A-$B-$qerun.sh
       sbatch $A-$B-$qerun.sh
     else
       sbatch run-$qerun.sh
     fi
}
