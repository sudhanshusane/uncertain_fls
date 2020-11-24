#!/bin/bash

FILE="/home/sci/ssane/data/EthaneDiol/ethane_diol.vtk"
FIELDS=2
NEIGHBORHOOD=3

./StdDevGen $FILE $FIELDS $NEIGHBORHOOD
