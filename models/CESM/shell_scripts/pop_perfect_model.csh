#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# This block is an attempt to localize all the machine-specific
# changes to this script such that the same script can be used
# on multiple platforms. This will help us maintain the script.
# Search below for TIMECHECK to see what times this script will
# run.

echo "`date` -- BEGIN GENERATE POP TRUE STATE"

set nonomatch       # suppress "rm" warnings if wildcard does not match anything

# The FORCE options are not optional.
# The VERBOSE options are useful for debugging though
# some systems don't like the -v option to any of the following 
switch ("`hostname`")
   case be*:
      # NCAR "bluefire"
      set   MOVE = '/usr/local/bin/mv -fv'
      set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
      set   LINK = '/usr/local/bin/ln -fvs'
      set REMOVE = '/usr/local/bin/rm -fr'

      set BASEOBSDIR = /glade/proj3/image/Observations/WOD09
   breaksw

   case ys*:
      # NCAR "yellowstone"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set BASEOBSDIR = /glade/p/image/Observations/WOD09
   breaksw

   default:
      # NERSC "hopper"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set BASEOBSDIR = /scratch/scratchdirs/nscollin/ACARS
   breaksw
endsw

#-------------------------------------------------------------------------
# Determine time of model state ... from file name of first member
# of the form "./${CASE}.pop.r.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#-------------------------------------------------------------------------

set FILE = `head -n 1 rpointer.ocn.restart`
set FILE = $FILE:t
set FILE = $FILE:r
set OCN_DATE_EXT = `echo $FILE:e`
set OCN_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
set OCN_YEAR     = `echo $OCN_DATE[1] | bc`
set OCN_MONTH    = `echo $OCN_DATE[2] | bc`
set OCN_DAY      = `echo $OCN_DATE[3] | bc`
set OCN_SECONDS  = `echo $OCN_DATE[4] | bc`
set OCN_HOUR     = `echo $OCN_DATE[4] / 3600 | bc`

echo "valid time of model is $OCN_YEAR $OCN_MONTH $OCN_DAY $OCN_SECONDS (seconds)"
echo "valid time of model is $OCN_YEAR $OCN_MONTH $OCN_DAY $OCN_HOUR (hours)"

#-------------------------------------------------------------------------
# Determine if current time is 0Z - if so, generate perfect obs.
# If not, return before doing anything.
#-------------------------------------------------------------------------

## TIMECHECK:
if ( $OCN_HOUR == 0 ) then
   echo "Hour is $OCN_HOUR so we are generating perfect obs for the ocean"
else
   echo "Hour is not 0Z so we are skipping generating perfect obs for the ocean"
   echo "`date` -- END   GENERATE POP TRUE STATE"
   exit 0
endif

#-------------------------------------------------------------------------
# Create temporary working directory for the perfect model and go there
#-------------------------------------------------------------------------

set temp_dir = pmo_pop
echo "temp_dir is $temp_dir"

if ( -d $temp_dir ) then
   ${REMOVE} $temp_dir/*
else
   mkdir -p $temp_dir
endif
cd $temp_dir

#-----------------------------------------------------------------------------
# Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of POP.
#-----------------------------------------------------------------------------

set YYYYMM   = `printf %04d%02d                ${OCN_YEAR} ${OCN_MONTH}`
set OBSFNAME = `printf obs_seq.0Z.%04d%02d%02d ${OCN_YEAR} ${OCN_MONTH} ${OCN_DAY}`
set OBS_FILE = ${BASEOBSDIR}/${YYYYMM}/${OBSFNAME}

if (  -e   ${OBS_FILE} ) then
   ${LINK} ${OBS_FILE} obs_seq.in
else
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... no observation file $OBS_FILE"
   exit -1
endif

#=========================================================================
# Block 1: Populate a run-time directory with the input needed to run DART,
# and Block 2: Convert 1 POP restart file to DART initial conditions file.
# At the end of the block, we have DART initial condition file  perfect_ics
# that came from pointer file ../rpointer.ocn.restart
#
# REQUIRED DART namelist settings:
# &perfect_model_obs_nml:  restart_in_file_name    = 'perfect_ics'
# &pop_to_dart_nml:        pop_to_dart_output_file = 'dart_ics'
#=========================================================================

echo "`date` -- BEGIN POP-TO-DART"

   if ( ! -e   ${CASEROOT}/pop_input.nml ) then
      echo "ERROR ... DART required file ${CASEROOT}/pop_input.nml not found ... ERROR"
      echo "ERROR ... DART required file ${CASEROOT}/pop_input.nml not found ... ERROR"
      exit -2
   endif

   # make sure there are no old output logs hanging around
   $REMOVE output.pop_to_dart

   set OCN_RESTART_FILENAME = ../${CASE}.pop.r.${OCN_DATE_EXT}.nc
   set     OCN_NML_FILENAME = ../pop2_in
   set     DART_IC_FILENAME = perfect_ics

   sed -e "s#dart_ics#${DART_IC_FILENAME}#" < ${CASEROOT}/pop_input.nml >! input.nml

   ${LINK} $OCN_RESTART_FILENAME pop.r.nc
   ${LINK} $OCN_NML_FILENAME     pop_in

   ${EXEROOT}/pop_to_dart >! output.pop_to_dart

   if ($status != 0) then
      echo "ERROR ... DART died in 'pop_to_dart' ... ERROR"
      echo "ERROR ... DART died in 'pop_to_dart' ... ERROR"
      exit -3
   endif

   cd ..

echo "`date` -- END POP-TO-DART"

#=========================================================================
# Block 3: Run perfect_model_obs and harvest the synthetic observations
# and diagnostic files.
#
# DART namelist settings required:
# &perfect_model_obs_nml:           async                  = 0,
# &perfect_model_obs_nml:           adv_ens_command        = "no_advance_script",
# &perfect_model_obs_nml:           output_restart         = .false.,
# &perfect_model_obs_nml:           restart_in_file_name   = 'perfect_ics'
# &perfect_model_obs_nml:           restart_out_file_name  = 'not_created'
# &perfect_model_obs_nml:           obs_sequence_in_name   = 'obs_seq.in'
# &perfect_model_obs_nml:           obs_sequence_out_name  = 'obs_seq.perfect'
# &perfect_model_obs_nml:           init_time_days         = -1,
# &perfect_model_obs_nml:           init_time_seconds      = -1,
# &perfect_model_obs_nml:           first_obs_days         = -1,
# &perfect_model_obs_nml:           first_obs_seconds      = -1,
# &perfect_model_obs_nml:           last_obs_days          = -1,
# &perfect_model_obs_nml:           last_obs_seconds       = -1,
#
#=========================================================================


echo "`date` -- BEGIN POP PERFECT_MODEL_OBS"
${EXEROOT}/perfect_model_obs_pop || exit -7
echo "`date` -- END   POP PERFECT_MODEL_OBS"


${MOVE} True_State.nc      ../pop_True_State.${OCN_DATE_EXT}.nc
${MOVE} Prior_Diag.nc      ../pop_Prior_Diag.${OCN_DATE_EXT}.nc
${MOVE} Posterior_Diag.nc  ../pop_Posterior_Diag.${OCN_DATE_EXT}.nc
${MOVE} obs_seq.perfect    ../pop_obs_seq.${OCN_DATE_EXT}.perfect
${MOVE} dart_log.out       ../pop_dart_log.${OCN_DATE_EXT}.out

#=========================================================================
# Block 4: Update the pop restart files.
#=========================================================================

# not needed ... perfect_model_obs does not update the model state.


#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

echo "`date` -- END   GENERATE POP TRUE STATE"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

