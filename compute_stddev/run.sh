#!/bin/bash

FILE="/home/sci/ssane/data/EthaneDiol/ethane_diol.vtk"
FIELDS=2
NEIGHBORHOOD=3
OUTPUT="/home/sci/ssane/data/EthaneDiol/uncertain_ethane_diol.vtk"

./StdDevGen $FILE $FIELDS $NEIGHBORHOOD $OUTPUT
