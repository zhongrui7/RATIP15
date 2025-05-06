#!/bin/bash


# Example 1: Double Ionization of Neon Atoms
# for the 2017 Joint ICTP-IAEA School on Atomic Processes in Plasmas
# 
# Authors: Stephan Fritzsche, Randolf Beerwerth
#          S.Fritzsche@gsi.de, Randolf.Beerwerth@uni-jena.de


# define the name for the wave function files to be generated
name="ne-0+"


# generate the configuration state list (CSL) for the ground configuration
# of neutral neon
$g2k/csl <<EOF
y

1s 2s 2p
n
n
1s(2) 2s(2) 2p(6)


0
EOF

# The resulting output file needs to be copied to the input of the next program
mv rcsl.out rcsf.inp


# Perform the angular integration
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
$g2k/rmcdhf <<EOF
y
1
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
# remark: The mcdf program (rmcdhf) does in principle solve the
# eigenvalue problem. However, the resulting mixing coefficients file (*.m, *.cm)
# is not compatible with RATIP's components, therefore we use
# RATP to generate a compatible output file

$ratip/xrelci <<EOF
${name}-ci.sum
${name}.c
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
1
EOF

