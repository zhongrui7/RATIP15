#!/bin/bash


# Example 1: Double Ionization of Neon Atoms
# for the 2017 Joint ICTP-IAEA School on Atomic Processes in Plasmas
# 
# Authors: Stephan Fritzsche, Randolf Beerwerth
#          S.Fritzsche@gsi.de, Randolf.Beerwerth@uni-jena.de


# Define a neon nucleus with Z = 10 and M = 20
# Default estimate for the nuclear radius is used, 
# and nuclear moments remain undefined
$g2k/rnucleus <<EOF
10
20
n
20
0
0
0

EOF
