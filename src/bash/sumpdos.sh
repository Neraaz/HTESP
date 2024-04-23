#!/bin/bash
# Writen by Niraj K. Nepal, Ph.D.
# script to extract partial DOS using sumpdos.x script
# Usage: sumpdos.sh atom_name orbital_name

atom=$1
orb=$2
x=$atom
y=$orb
sumpdos.x *\($atom\)* > $atom-tot.dat
sumpdos.x *\($orb\)* > $orb-tot.dat
sumpdos.x *\($atom\)*\($orb\) > $atom-$orb.dat
sumpdos.x *\($atom\)*\($orb\)* > $atom-$orb.dat
