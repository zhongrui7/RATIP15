#!/bin/bash


# Example 1: Double Ionization of Neon Atoms
# for the 2017 Joint ICTP-IAEA School on Atomic Processes in Plasmas
# 
# Authors: Stephan Fritzsche, Randolf Beerwerth
#          S.Fritzsche@gsi.de, Randolf.Beerwerth@uni-jena.de


# define the name for the wave function files to be generated
name="ne-2+-A"

# generate the CSL for the three configurations in doubly-ionized neon that
# can be reached by Auger processes from K-shell ionized neon ions
$g2k/csl <<EOF
y

1s 2s 2p
n
n
1s(2) 2s(0) 2p(6)

1s(2) 2s(1) 2p(5)

1s(2) 2s(2) 2p(4)


0
EOF

# As before, save the original CSL file for RATIP
cp rcsl.out ${name}.noblock.c
mv rcsl.out rcsf.inp

# Then perform the separation into blocks of equal angular momentum
$g2k/rcsfblock <<EOF
y
EOF

mv rcsf.out rcsf.inp


# Perform the angular integration, same as before
$g2k/rangular <<EOF
y
EOF

# Generate an initial estimate for the radial orbitals
$g2k/rwfnestimate <<EOF
y
2
*
EOF


# Iterative solution of the radial wave functions with the MCDF method
# This time a total of 10 CSFs/fine-structure levels is given that need
# to be optimized
$g2k/rmcdhf <<EOF
y
3
1
1
2
2
1
5
*
*
100
EOF

# Store the resulting wave functions under the given name
$g2k/rsave ${name}


# Remove the temporary files resulting from the angular integration
rm mcp.*

# Remove, if existing, the previous results
rm ${name}.cm
rm ${name}-ci.sum

# Solve the configuration interaction problem using RATIP
# now computation of the 10 wave functions with serial numbers
# 1-10 (last line)
$ratip/xrelci <<EOF
${name}-ci.sum
${name}.noblock.c
isodata
y
y
n
n
n
y
y
n
eV
n
n
${name}.w
${name}.cm
1-10
EOF

