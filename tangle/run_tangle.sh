#!/bin/bash

INPUT1="/home/sci/ssane/data/Tangle/vtk/meanTangle.vtk"
INPUT2="/home/sci/ssane/data/Tangle/vtk/devTangle.vtk"
OUTPUT="/home/sci/ssane/data/Tangle/UFLS_Tangle2.vtk"
ATT1_MIN="0"
ATT1_MAX="62"

ATTR1="meanTangle"
ATTR2="devTangle"

./FLS $INPUT1 $INPUT2 $OUTPUT $ATT1_MIN $ATT1_MAX $ATTR1 $ATTR2


