#!/bin/bash


# Example 1: Double Ionization of Neon Atoms
# for the 2017 Joint ICTP-IAEA School on Atomic Processes in Plasmas
# 
# Authors: Stephan Fritzsche, Randolf Beerwerth
#          S.Fritzsche@gsi.de, Randolf.Beerwerth@uni-jena.de


# remove the previous result, if existing
rm ratip_photo.sum

# Utilized RATIP's photo component to compute the photo ionization cross 
# sections in an energy range from 870 to 1500 eV in 10eV steps
$ratip/xphoto <<EOL
ratip_photo.sum
isodata
E1
eV
1000
y
870 1500 10
y
n
y
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

n
n
ne-0+.c
ne-1+-hole.c
ne-0+.cm
ne-1+-hole.cm
ne-0+.w
ne-1+-hole.w
EOL


