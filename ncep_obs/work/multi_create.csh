#!/bin/csh  -v
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2007, Data Assimilation Research Section
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# test version - convert multiple days in parallel
#
#--------------------------------------------------------------
# DESCRIPTION:
#   Run create_obs_seq for the timespan described in input.nml
#   Up to a month at a time, in parallel.
#
#--------------------------------------------------------------

# you might need other options set for LSF on your system; 
#  e.g. -P project_charge_code, a different queue name, etc.

#BSUB -J createobs
#BSUB -o createobs.out
#BSUB -e createobs.err
#BSUB -q economy
#BSUB -P 00000000
#BSUB -W 3:30
#BSUB -n 1

# Set the year, month, and day range here.
# These 4 vars are the only user inputs to this script.
set     year = 2008
set    month =    0
set startday =    1
set   endday =   31


# Loop over days, first to convert, then to clean up
set day = $startday
set lastday = $endday

while ( $day <= $lastday )

   mkdir tempdir_${day}
   cd tempdir_${day}

   set mo = `printf %02d $month`
   cp -f ../create_real_obs .
   sed -e "s/YYYY/${year}/" \
       -e "s/MM/${mo}/"     \
       -e "s/NN/${day}/"    < ../input.nml.template >! input.nml

   echo starting day $day
   ./create_real_obs >! obs.log &

   cd ..

   @ day++
end

wait

echo all obs_seq files created at `date`

set day = $startday
set lastday = $endday

# remove the various bits but leave the newly created
# obs_seq file there.  the script could move it automatically
# but for now i am going to look at them by hand.  if they run
# ok for a couple trials, then i will move them to where they
# belong.  must then set a base dir to move them into.
# could then remove the dir, although might miss errors that way.
while ( $day <= $lastday )
 rm tempdir_${day}/input.nml tempdir_${day}/create_real_obs tempdir_${day}/dart_log.nml
 #rmdir tempdir_${day} 
 @ day++
end

exit 0

# <next few lines under version control, do not edit>
# $URL: http://subversion.ucar.edu/DAReS/DART/trunk/ncep_obs/work/multi_create.csh $
# $Id: prepbufr.csh 3042 2007-08-02 22:03:02Z thoar $
# $Revision: 3042 $
# $Date: 2007-08-02 16:03:02 -0600 (Thu, 02 Aug 2007) $
