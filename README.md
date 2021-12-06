# GLASCOW
#======================================================================

#-Archive-description--------------------------------------------------

This archive is divided in 5 main directories.

BIN contains a static binary of Glascow

MAN contains GLASCOW's pdf manual

SRC contains source code and a Makefile

EX contains 4 directories each with examples datasets.

.CHK is for developpement purpose and run tests on sample data.

#-Installation--------------------------------------------------------

By default, GLASCOW will be compiled with ifort and LAPACK (mkl).  Use
of IFORT and LAPACK are strongly recommended because they greatly
improve speed of computations. Nevertheless, compilation flag for
gfortran use are also provided. You'll have to change compiler
definition in the Makefile, some hints are given in it.

The binary provided is statically compiled and should therefore work
as is on x86_64 linux platform. 
