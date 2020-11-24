#!/bin/bash

INPUT="/home/sci/ssane/data/KarmanVortexStreet_Unsteady/cylinder2d_unsteady.vtk"
FIELD=1
OUTPUT="OutputKVS_CI_FLS"$FIELD".vtk"
ATT1_MIN="1.5"
ATT1_MAX="1.834"

ATT2_MIN="-0.01"
ATT2_MAX="0.01"

ATTR1="u"
ATTR2="v"

STDDEV1="0.05"
STDDEV2="0.05"

./FLS $INPUT $OUTPUT $ATT1_MIN $ATT1_MAX $ATT2_MIN $ATT2_MAX $ATTR1 $ATTR2 $STDDEV1 $STDDEV2


