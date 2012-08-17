#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# This is an example script for how to stage the files in CENTRALDIR
# in preparation for an assimilation.
#
#==============================================================================
# Set the commands so we can avoid problems with aliases, etc.
#==============================================================================

set   MOVE = '/usr/local/bin/mv -fv'
set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
set   LINK = '/usr/local/bin/ln -fvs'
set REMOVE = '/usr/local/bin/rm -fr'

set   MOVE = '/bin/mv -fv'
set   COPY = '/bin/cp -fvp'
set   LINK = '/bin/ln -fvs'
set REMOVE = '/bin/rm -fr'

#==============================================================================
# Stage all the required files in CENTRALDIR
#
# CENTRALDIR is where 'filter' will run, each model advance takes place
# in a subdirectory created and populated by 'advance_model.csh'
#
# The files may exist from a perfect model setup, use them. If not, get them.
#==============================================================================

set CENTRALDIR = `pwd`

set NOAHDIR = /Users/thoar/svn/DART/devel/models/noah/src/hrldas-v3.3
set DARTDIR = /Users/thoar/svn/DART/devel/models/noah

foreach FILE ( Noah_hrldas_beta SOILPARM.TBL VEGPARM.TBL GENPARM.TBL URBPARM.TBL )
   if ( -e ${FILE} )  then
      echo "Using existing $FILE"
   else
      ${COPY} ${NOAHDIR}/Run/${FILE} . || exit 1
   endif
end

if ( -e wrfinput ) then
   echo "Using existing wrfinput"
else
   ${COPY} ${NOAHDIR}/Run/wrfinput.template wrfinput  || exit 1
endif

foreach FILE ( namelist.hrldas obs_seq.out input.nml filter dart_to_noah noah_to_dart restart_file_tool )
   if ( -e ${FILE} )  then
      echo "Using existing $FILE"
   else
      echo "$FILE needs to be copied."
      ${COPY} ${DARTDIR}/work/${FILE} . || exit 1
   endif
end

foreach FILE ( run_filter.csh advance_model.csh )
   if ( -e ${FILE} )  then
      echo "Using existing $FILE"
   else
      ${COPY} ${DARTDIR}/shell_scripts/${FILE} . || exit 1
   endif
end

#==============================================================================
# need a set of noah restart files to define the initial ensemble
#
# 1) point to a directory full of noah restart files and pick N of them;
# 2) convert them to DART format;
# 3) make sure the time in each DART file is 'identical'
# 4) the original noah restart files are also needed to start up
#    advance_model.csh for the first model advance. The tricky part is
#    that the Time variable in those files is all wrong.
#
# NOTE : This is for testing the machinery ONLY! If you try to publish a paper
# with this ensemble I will REJECT IT AND EVERY OTHER PAPER IN YOUR CAREER.
#
# COME UP WITH YOUR OWN ENSEMBLE.
#==============================================================================

set ENSEMBLESOURCE = /Users/thoar/svn/DART/devel/models/noah/src/hrldas-v3.3/Run/hourly_output

set   nfiles = `ls -1 ${ENSEMBLESOURCE}/RESTART*DOMAIN* | wc -l`
set filelist = `ls -1 ${ENSEMBLESOURCE}/RESTART*DOMAIN*`

@ ifile = 1
@ ensemble_member = 0
while ($ifile <= $nfiles)

   @ ensemble_member = $ensemble_member + 1
   set fext = `printf %04d $ensemble_member`

   ${COPY} $filelist[$ifile] restart.nc

   # change the time here ... and simplify life

   ./noah_to_dart                     || exit 3
   ${MOVE} dart_ics filter_ics.$fext  || exit 4

   @ ifile = $ifile + 10
end

#==============================================================================
# Now we make sure all the initial ensemble files have the same time.
#==============================================================================

echo "Make sure the earliest time in the obs_seq.out file is at or after"
echo "the time we are inserting in the initial ensemble."
echo "DART can advance the model states to the observation time."
echo "DART cannot move the model state back in time."
 
set MODEL_DAY    = 147192
set MODEL_SECOND = 0

ex input.nml <<ex_end
/restart_file_tool_nml
g;ens_size ;s;= .*;= 1,;
g;single_restart_file_in ;s;= .*;= .true.,;
g;single_restart_file_out ;s;= .*;= .true.,;
g;write_binary_restart_files ;s;= .*;= .false.,;
g;overwrite_data_time ;s;= .*;= .true.,;
g;new_data_days ;s;= .*;= $MODEL_DAY,;
g;new_data_secs ;s;= .*;= $MODEL_SECOND,;
g;input_is_model_advance_file ;s;= .*;= .false.,;
g;output_is_model_advance_file ;s;= .*;= .false.,;
g;gregorian_cal ;s;= .*;= .true.;
wq
ex_end

foreach FILE ( filter_ics.* )

   ln -svf ${FILE} filter_restart
   ./restart_file_tool 
   ${MOVE} filter_updated_restart ${FILE}

   ${REMOVE} filter_restart
end

# Since we have some knowledge of the ensemble size, 
# provide reasonable default values.

ex input.nml <<ex_end
/filter
g;ens_size ;s;= .*;= $ensemble_member,;
g;num_output_state_members ;s;= .*;= $ensemble_member,;
g;num_output_obs_members ;s;= .*;= $ensemble_member,;
g;single_restart_file_in ;s;= .*;= .false.,;
g;single_restart_file_out ;s;= .*;= .false.,;
wq
ex_end

#==============================================================================
# Finish up.
#==============================================================================

echo
echo "CENTRALDIR is ${CENTRALDIR}"
echo "Configure the ${CENTRALDIR}/input.nml"
echo "Configure the ${CENTRALDIR}/namelist.hrldas"
echo "execute       ${CENTRALDIR}/run_filter.csh"
echo "Configure     ${CENTRALDIR}/wrfinput"
echo

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

