#!/bin/bash
#
# JBR : 1/18/21
#
# This script cycles through the directories and compiles all required FORTRAN 
# libraries *.a
#
# CURRENTLY ASSUMES GFORTRAN!
#

echo "Building FORTRAN libraries"

###### libcip ########################
echo "libcip"
cd ./libcip
make -f makefile_gfortran
cd ..

###### libtau ########################
echo "libtau"
cd ./libtau
make -f makefile_gfortran
cd ..

###### libharm ########################
echo "libharm"
cd ./libharm
make -f makefile_gfortran
cd ..

###### libutil ########################
echo "libutil"
cd ./libutil
make -f makefile_gfortran
cd ..
