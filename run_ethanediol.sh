#!/bin/bash

INPUT="/home/sci/ssane/data/EthaneDiol/ethane_diol.vtk"
FIELD=4
OUTPUT="OutputEthane_CI_FLS"$FIELD".vtk"
ATT1_MIN="2"
ATT1_MAX="2.2"

ATT2_MIN="-0.5"
ATT2_MAX="0.0"

ATTR1="log_lp_Rho_rp_"
ATTR2="log_lp_s_rp_"

STDDEV1="0.2"
STDDEV2="0.2"

./FLS $INPUT $OUTPUT $ATT1_MIN $ATT1_MAX $ATT2_MIN $ATT2_MAX $ATTR1 $ATTR2 $STDDEV1 $STDDEV2


