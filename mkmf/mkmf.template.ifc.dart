# template for Intel 7.1 Fortran Compiler on a RedHat 7.3 architecture 
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next two lines are maintained by CVS, please do not edit>
# $Id$
# $Name$
#
# typical use with mkmf
# mkmf -t mkmf.template.xxxx -c"-Duse_netCDF" ...
#
# NETCDF and LIBS needs to be customized for your site
#       We had the netcdf libraries built right in: /opt/Intel/[lib,include]
#
# MPICH libraries are only needed for Bgrid model and only when you turn on
#       the -Duse_libMPI preprocessor flag.
#
# FFLAGS   used all the time
#    -fpp
#
# FFLAGS   for benchmarking
#    -O0
#    -pc64
#
# FFLAGS   for debugging
#    -g
#
# FFLAGS   for production
#    -O2
#
# FFLAGS   no longer needed as of Hawaii release
#    -i4
#    -r8
#
# LDFLAGS
#     -static -KPIC    needed for relocatable code
#
FC = ifc
LD = ifc
NETCDF = /opt/Intel
INCS = -I$(NETCDF)/include
FFLAGS = -i4 -r8 -fpp -O2 $(INCS)
LIBS = -L/usr/local/lib -lPEPCF90 -L$(NETCDF)/lib -ludunits -lnetcdf -L/usr/mpi/lib 
LDFLAGS = $(LIBS)
