#!/bin/csh
#
# DART software - Copyright � 2004 - 2010 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id: quickbuild.csh 4504 2010-09-24 21:05:09Z nancy $
#
# This script compiles all executables in this directory.

set MODEL = "system_simulation"

@ n = 0

#----------------------------------------------------------------------
# Build all the single-threaded targets
#----------------------------------------------------------------------

foreach TARGET ( mkmf_* )

   set PROG = `echo $TARGET | sed -e 's#mkmf_##'`

   @ n = $n + 1
   echo
   echo "---------------------------------------------------"
   echo "${MODEL} build number ${n} is ${PROG}" 
   \rm -f ${PROG}
   csh $TARGET || exit $n
   make        || exit $n

end

# clean up.  comment this out if you want to keep the .o and .mod files around
\rm -f *.o *.mod input.nml.*_default

echo "Success: All DART programs compiled."  

exit 0

# <next few lines under version control, do not edit>
# $URL: https://proxy.subversion.ucar.edu/DAReS/DART/trunk/observations/gps/work/quickbuild.csh $
# $Revision: 4504 $
# $Date: 2010-09-24 15:05:09 -0600 (Fri, 24 Sep 2010) $

