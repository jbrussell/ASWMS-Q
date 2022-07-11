#!/bin/bash
#
# JBR : 8/19/18
# Script to quickly compile MINEOS and idagrn6 fortran codes
#
# COMPILER : specify the compiler... I've only tested gfortran
# LIBPATH : path to libraries used
# BINPATH : path to output binaries
#

COMPILER="gfortran"
LIBPATH="$PWD/libgfortran"
BINPATH="$PWD/bin"

echo "Compiler: $COMPILER"

##########################################################
#                         MINEOS                         #
##########################################################
echo "Building MINEOS"

###### plot_wk ########################
echo "plot_wk"
MAKEPATH="./MINEOS/prog_plotwk/"

MAKEFILE="1_plot_wk"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o

###### frechet ########################
echo "frechet"
MAKEPATH="./MINEOS/prog_frechet/"

MAKEFILE="1_frechet_all"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o

MAKEFILE="2_frechet_ACFLN_love"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o

MAKEFILE="3_frechet_Q"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o

###### mineos ########################
echo "mineos"
MAKEPATH="./MINEOS/prog_mineos/"

MAKEFILE="1_eig_recover"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o

MAKEFILE="2_mineos_qcorrectphv"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o

MAKEFILE="3_mineos_strip"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o

MAKEFILE="4_mineos_table"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o

MAKEFILE="5_get_eigfxn_grvelo"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o

MAKEFILE="5_get_eigfxn_grvelo_int"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o

MAKEPATH="./MINEOS/prog_mineos/Mineos/"
MAKEFILE="1_mineos_nohang"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o

###### eigen ########################
echo "eigen"
MAKEPATH="./MINEOS/prog_eigen/"

MAKEFILE="1_eigenST_asc"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o

##########################################################
#                        IDAGRN6                         #
##########################################################
echo "idagrn6"
MAKEPATH="./idagrn6/prog/"

MAKEFILE="1_idagrn6_sac"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o

MAKEFILE="2_idagrn6_sac_excite"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o

MAKEFILE="3_idagrn6_excite"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o

###### mask ########################
echo "idagrn6_mask"
MAKEPATH="./idagrn6/prog_mask/"
MAKEFILE="1_idagrn6_mask"
sed -e "s/FCxreplace/${COMPILER}/ ; s+BINxreplace+${BINPATH}+ ; s+LIBxreplace+${LIBPATH}+" ${MAKEPATH}${MAKEFILE}.mkin > ${MAKEPATH}${MAKEFILE}.mk
(cd ${MAKEPATH} && make clean -f ${MAKEFILE}.mk)
(cd ${MAKEPATH} && make -f ${MAKEFILE}.mk)
rm -r ${MAKEPATH}/*.o