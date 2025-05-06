#!/bin/bash

# Example 1: Double Ionization of Neon Atoms
# for the 2017 Joint ICTP-IAEA School on Atomic Processes in Plasmas
# 
# Authors: Stephan Fritzsche, Randolf Beerwerth
#          S.Fritzsche@gsi.de, Randolf.Beerwerth@uni-jena.de


# This scripts executes all scripts that are necessary to perform the computations
# for example 1


source 00_path.sh

./01_iso.sh

./02_wave_0+.sh

./03_wave_1+-hole.sh

./04_wave_1+-gs.sh

./05_wave_2+.sh

./06_photo.sh

./07_auger.sh

./08_fluorescence.sh




