#!/bin/bash

INPUT="/home/sci/ssane/data/SyntheticCylinderFlow_Unsteady/synthetic_cylinderflow_unsteady.vtk"
FIELD=1
OUTPUT="OutputSCF_CI_FLS"$FIELD".vtk"
ATT1_MIN="-0.01"
ATT1_MAX="0.01"

ATT2_MIN="-0.01"
ATT2_MAX="0.01"

ATTR1="u"
ATTR2="v"

STDDEV1="0.01"
STDDEV2="0.01"

./FLS $INPUT $OUTPUT $ATT1_MIN $ATT1_MAX $ATT2_MIN $ATT2_MAX $ATTR1 $ATTR2 $STDDEV1 $STDDEV2


