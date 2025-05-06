#!/bin/bash


# Example 1: Double Ionization of Neon Atoms
# for the 2017 Joint ICTP-IAEA School on Atomic Processes in Plasmas
# 
# Authors: Stephan Fritzsche, Randolf Beerwerth
#          S.Fritzsche@gsi.de, Randolf.Beerwerth@uni-jena.de


# define the name for the wave function files to be generated
name="ne-1+-gs"


# generate the configuration state list for the ground (and first excited)
# configuration of singly-ionized neon
$g2k/csl <<EOF
y

1s 2s 2p
n
n
1s(2) 2s(2) 2p(5)

1s(2) 2s(1) 2p(6)


0
EOF

# Here we need a copy of the resulting CSL that is not seperated into blocks of equal 
# total angular momentum for the RATIP components
cp rcsl.out ${name}.noblock.c

# This copies the CSL to the input of the next GRASP program
mv rcsl.out rcsf.inp

# This program separates the CSL into blocks of equal total angular momentum
$g2k/rcsfblock <<EOF
y
EOF

# again move the resulting CSL to the input for the angular integration
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
# Now we have three CSFs and hence three fine-structure levels with J = 1/2+, 1/2-, 3/2-
# These need to be stated seperately for the optimization
$g2k/rmcdhf <<EOF
y
1
1
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
# RATIP does not use the same separation in blocks of equal J, hence the 
# levels that are to be computed are given in the form 1-3, meaning, that all 
# three wave functions should be computed.
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
1-3
EOF

