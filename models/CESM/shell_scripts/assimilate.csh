#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$

# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.

echo "`date` -- BEGIN CESM ASSIMILATE"

./cam_assimilate.csh

# for now, just call the cam script
#./pop_assimilate.csh
#./clm_assimilate.csh

echo "`date` -- END CESM ASSIMILATE"

exit 0

