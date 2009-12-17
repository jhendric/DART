#!/bin/sh 
#
# DART software - Copyright � 2004 - 2010 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# This script will make executables which extract data
# from ADP BUFR input files, and place the data into a basic text file.
# prepbufr.x:  used to extract data from prepbufr files
# ** Make sure the "ar" command location has been set in your path
# environment variable.  Type "which ar" to check if this is done. **
#
# if you get a link or runtime error about bort being undefined, go into
# the lib directory and read the README_BUFRLIB file about how to fix it.
# you can do something like:  cc='cc -DUNDERSCORE' to set the flag for all
# compiles at once.

set -eua
 
#  ------------------------------------------------------------------------
#  CPLAT - platform type (linux, sgi, aix, sun, macosx)
#  ------------------------------------------------------------------------
 
#CPLAT=macosx
#CPLAT=linux
CPLAT=aix

#  set up the compilers to use
#  -----------------------------------------------------
 
if [ $CPLAT = sgi ]
then
   cc=cc; ff=f77
elif [ $CPLAT = linux ]
then
#   cc=cc; ff=pgf90
   cc=cc; ff=ifort
elif [ $CPLAT = aix ]
then
   cc='cc -O'; ff='f77 -O'
elif [ $CPLAT = sun ]
then
   cc=cc; ff=f77
elif [ $CPLAT = macosx ]
then
   cc=gcc; ff=gfortran
fi

#  Compile and archive the Bufr Library
#  ------------------------------------

echo 'Compiling the Bufr library'

cd lib
$ff -c *.f
$cc -c *.c
ar crv bufrlib.a *.o
rm *.o
 
#  Compile and link the decode programs
#  ---------------------------------------

echo 'Compiling the prepbufr programs'

cd ../src
$ff prepbufr.f     ../lib/bufrlib.a -o ../exe/prepbufr.x 
$ff prepbufr_03Z.f ../lib/bufrlib.a -o ../exe/prepbufr_03Z.x 

 
# ONLY IF NEEDED, compile the 2 auxiliary conversion programs
#  (comment one or both of these sections in)

#  Compile the binary format converter program
#  ---------------------------------------
 
echo 'Compiling the grabbufr converter program'
 
cd ../convert_bufr
$ff grabbufr.f ../lib/bufrlib.a -o ../exe/grabbufr.x
 
 
#  Compile the block/unblock converter program
#  ---------------------------------------
 
#echo 'Compiling the blk/ublk converter program'
 
#cd ../blk_ublk
#./cwordsh.make
 

#  clean up
#  --------

echo 'Finished making executables'
cd ..

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

