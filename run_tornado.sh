#!/bin/bash

INPUT="/home/sci/ssane/data/Tornado/Tornado_07100.vtk"
OUTPUT="OutputTornado_CI.vtk"
#ATT2_MIN="-20.0"
#ATT2_MAX="-19.0"
#ATT2_MIN="-0.5"
#ATT2_MAX="0.5"
ATT1_MIN="0.4"
ATT1_MAX="2.2"

ATT2_MIN="-73"
ATT2_MAX="-20"

ATTR1="vortmag"
ATTR2="prespert"

STDDEV1="0.2"
STDDEV2="5"

./FLS $INPUT $OUTPUT $ATT1_MIN $ATT1_MAX $ATT2_MIN $ATT2_MAX $ATTR1 $ATTR2 $STDDEV1 $STDDEV2


