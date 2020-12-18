#!/bin/bash

INPUT1="/home/sci/ssane/data/RedSeaEddy/velocityMagnitude/mean/meanVelocityMagnitude.vtk"
INPUT2="/home/sci/ssane/data/RedSeaEddy/velocityMagnitude/dev/devVelocityMagnitude.vtk"
INPUT3="/home/sci/ssane/data/RedSeaEddy/temperature/mean/meanTemperature.vtk"
INPUT4="/home/sci/ssane/data/RedSeaEddy/temperature/dev/devTemperature.vtk"
OUTPUT="/home/sci/ssane/data/RedSeaEddy/UFLS_RedSeaEddy.vtk"
ATT1_MIN="0.15"
ATT1_MAX="0.25"
ATT2_MIN="20"
ATT2_MAX="30"

ATTR1="meanVelocityMagnitude"
ATTR2="devVelocityMagnitude"
ATTR3="meanTemperature"
ATTR4="devTemperature"

./FLS $INPUT1 $INPUT2 $INPUT3 $INPUT4 $OUTPUT $ATT1_MIN $ATT1_MAX $ATT2_MIN $ATT2_MAX $ATTR1 $ATTR2 $ATTR3 $ATTR4 


