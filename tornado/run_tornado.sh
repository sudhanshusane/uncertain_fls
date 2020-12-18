#!/bin/bash

INPUT="/home/sci/ssane/data/Tornado/Tornado_07100.vtk"
OUTPUT="OutputTornado_CI.vtk"
ATT1_MIN="0.4"
ATT1_MAX="2.2"

ATT2_MIN="-73"
ATT2_MAX="-20"

ATTR1="vortmag"
ATTR2="prespert"

./FLS $INPUT $OUTPUT $ATT1_MIN $ATT1_MAX $ATT2_MIN $ATT2_MAX $ATTR1 $ATTR2


