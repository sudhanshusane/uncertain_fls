#!/bin/bash

FILE="/home/sci/ssane/data/Tornado/Tornado_07100.vtk"
#FILE="/home/sci/ssane/data/EthaneDiol/ethane_diol.vtk"
FIELDS=4
NEIGHBORHOOD=5
OUTPUT="/home/sci/ssane/data/Tornado/uncertain_Tornado.vtk"
#OUTPUT="/home/sci/ssane/data/EthaneDiol/uncertain_ethane_diol.vtk"

./StdDevGen $FILE $FIELDS $NEIGHBORHOOD $OUTPUT
