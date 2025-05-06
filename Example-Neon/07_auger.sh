#!/bin/bash


# Example 1: Double Ionization of Neon Atoms
# for the 2017 Joint ICTP-IAEA School on Atomic Processes in Plasmas
# 
# Authors: Stephan Fritzsche, Randolf Beerwerth
#          S.Fritzsche@gsi.de, Randolf.Beerwerth@uni-jena.de


# remove the previous result, if already existing
rm ratip_auger.sum


# Utilized RATIP's Auger component to compute the Auger transition rates
$ratip/xauger <<EOL
ratip_auger.sum
isodata
eV
1000
y
n
y
n
n
n
n
n
n
y
n
y
n
0
0

n
n
ne-1+-hole.c
ne-2+-A.noblock.c
ne-1+-hole.cm
ne-2+-A.cm
ne-1+-hole.w
ne-2+-A.w
EOL

