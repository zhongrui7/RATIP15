#!/bin/bash


# Example 1: Double Ionization of Neon Atoms
# for the 2017 Joint ICTP-IAEA School on Atomic Processes in Plasmas
# 
# Authors: Stephan Fritzsche, Randolf Beerwerth
#          S.Fritzsche@gsi.de, Randolf.Beerwerth@uni-jena.de



# remove the previously existing results
rm ne-1+-hole-cesd.sum
rm ne-1+-hole.xpn

# expand the initial wave function in Slater determinants
$ratip/xcesd <<EOL
ne-1+-hole-cesd.sum
ne-1+-hole.c
y
n
n
ne-1+-hole.cm
ne-1+-hole.xpn
EOL

# remove the previously existing results
rm ne-1+-gs-cesd.sum
rm ne-1+-gs.xpn

# expand the final wave function in Slater determinants
$ratip/xcesd <<EOL
ne-1+-gs-cesd.sum
ne-1+-gs.noblock.c
y
n
n
ne-1+-gs.cm
ne-1+-gs.xpn
EOL


# remove the existing result for the transition rates
rm ratip_reos.sum

# Compute the transition rates utilizing the REOS component of RATIP
# and the previously generated expansions in Slater determinants
$ratip/xreos <<EOL
ratip_reos.sum
n
isodata
E1 M1
eV
y
n
n
y
n
n
n
n
y
n
n
n
n
ne-1+-hole.xpn
ne-1+-gs.xpn
ne-1+-hole.w
ne-1+-gs.w
EOL

