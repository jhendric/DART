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

echo "`date` -- BEGIN POP_ASSIMILATE"

switch ("`hostname`")
   case be*:
      # NCAR "bluefire"
      # The FORCE options are not optional.
      # the VERBOSE options are useful for debugging.
      set   MOVE = '/usr/local/bin/mv -fv'
      set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
      set   LINK = '/usr/local/bin/ln -fvs'
      set REMOVE = '/usr/local/bin/rm -fr'

      set BASEOBSDIR = /glade/proj3/image/Observations/WOD09
      set DARTDIR    = ${HOME}/svn/DART/dev
      set LAUNCHCMD  = mpirun.lsf
   breaksw

   case ys*:
      # NCAR "yellowstone"
      # The FORCE options are not optional.
      # the VERBOSE options are useful for debugging.
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set BASEOBSDIR = /glade/p/image/Observations/WOD09
      set DARTDIR    = ${HOME}/svn/DART/dev
      set LAUNCHCMD  = mpirun.lsf
   breaksw

   default:
      # NERSC "hopper"
      set   MOVE = 'mv -fv'
      set   COPY = 'cp -fv --preserve=timestamps'
      set   LINK = 'ln -fvs'
      set REMOVE = 'rm -fr'

      set BASEOBSDIR = /scratch/scratchdirs/nscollin/ACARS
      set DARTDIR    = ${HOME}/devel
      set LAUNCHCMD  = "aprun -n $NTASKS"
   breaksw
endsw

set ensemble_size = ${NINST_OCN}

# Create temporary working directory for the assimilation
set temp_dir = assimilate_pop
echo "temp_dir is $temp_dir"

# Create a clean temporary directory and go there
if ( -d $temp_dir ) then
   ${REMOVE} $temp_dir/*
else
   mkdir -p $temp_dir
endif
cd $temp_dir

#-------------------------------------------------------------------------
# Determine time of model state ... from file name of first member
# of the form "./${CASE}.pop_${ensemble_member}.r.2000-01-06-00000.nc"
#
# Piping stuff through 'bc' strips off any preceeding zeros.
#-------------------------------------------------------------------------

set FILE = `head -1 ../rpointer.ocn_0001.restart`
set FILE = $FILE:t
set FILE = $FILE:r
set MYCASE = `echo $FILE | sed -e "s#\..*##"`
set OCN_DATE_EXT = `echo $FILE:e`
set OCN_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
set OCN_YEAR     = `echo $OCN_DATE[1] | bc`
set OCN_MONTH    = `echo $OCN_DATE[2] | bc`
set OCN_DAY      = `echo $OCN_DATE[3] | bc`
set OCN_SECONDS  = `echo $OCN_DATE[4] | bc`
set OCN_HOUR     = `echo $OCN_DATE[4] / 3600 | bc`

echo "valid time of model is $OCN_YEAR $OCN_MONTH $OCN_DAY $OCN_SECONDS (seconds)"
echo "valid time of model is $OCN_YEAR $OCN_MONTH $OCN_DAY $OCN_HOUR (hours)"

#-----------------------------------------------------------------------------
# Get observation sequence file ... or die right away.
# The observation file names have a time that matches the stopping time of POP.
#-----------------------------------------------------------------------------

set YYYYMMDD = `printf %04d%02d%02d ${OCN_YEAR} ${OCN_MONTH} ${OCN_DAY}`
set YYYYMM   = `printf %04d%02d     ${OCN_YEAR} ${OCN_MONTH}`
set OBSFNAME = obs_seq.0Z.${YYYYMMDD}
set OBS_FILE = ${BASEOBSDIR}/${YYYYMM}/${OBSFNAME}

if (  -e   ${OBS_FILE} ) then
   ${LINK} ${OBS_FILE} obs_seq.out
else
   echo "ERROR ... no observation file $OBS_FILE"
   echo "ERROR ... no observation file $OBS_FILE"
   exit -1
endif

#=========================================================================
# Block 1: Populate a run-time directory with the input needed to run DART.
#=========================================================================

echo "`date` -- BEGIN COPY BLOCK"

if (  -e   ${CASEROOT}/input.nml ) then
   ${COPY} ${CASEROOT}/input.nml .
else
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   echo "ERROR ... DART required file ${CASEROOT}/input.nml not found ... ERROR"
   exit -2
endif

# Modify the DART input.nml such that
# the DART ensemble size matches the CESM number of instances
# WARNING: the output files contain ALL enemble members ==> BIG

ex input.nml <<ex_end
g;ens_size ;s;= .*;= $ensemble_size;
g;num_output_obs_members ;s;= .*;= $ensemble_size;
wq
ex_end

# COULD also  ERROR OUT if the important setting is not the same.
# problem is ... there are multiple instances of 'ens_size' in input.nml
#
echo "`date` -- END COPY BLOCK"

#=========================================================================
# Block 2: Stage the files needed for SAMPLING ERROR CORRECTION
#
# The sampling error correction is a lookup table.
# The tables are stored in the DART distribution.
# Each ensemble size has its own (static) file.
# It is only needed if
# input.nml:&assim_tools_nml:sampling_error_correction = .true.,
#=========================================================================

set  MYSTRING = `grep sampling_error_correction input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set SECSTRING = `echo $MYSTRING[2] | tr '[:upper:]' '[:lower:]'`

if ( $SECSTRING == true ) then
   set SAMP_ERR_FILE = ${DARTDIR}/system_simulation/final_full_precomputed_tables/final_full.${ensemble_size}
   if (  -e   ${SAMP_ERR_FILE} ) then
      ${COPY} ${SAMP_ERR_FILE} .
   else
      echo "ERROR: no sampling error correction file for this ensemble size."
      echo "ERROR: looking for ${SAMP_ERR_FILE}"
      exit -3
   endif
else
   echo "Sampling Error Correction not requested for this assimilation."
endif

#=========================================================================
# Block 3: DART INFLATION
# This stages the files that contain the inflation values.
# The inflation values change through time and should be archived.
#
# This file is only relevant if 'inflation' is turned on -
# i.e. if inf_flavor(1) /= 0 AND inf_initial_from_restart = .TRUE.
#
# filter_nml
# inf_flavor                  = 2,                       0,
# inf_initial_from_restart    = .true.,                  .false.,
# inf_in_file_name            = 'prior_inflate_ics',     'post_inflate_ics',
# inf_out_file_name           = 'prior_inflate_restart', 'post_inflate_restart',
# inf_diag_file_name          = 'prior_inflate_diag',    'post_inflate_diag',
#
# NOTICE: the archiving scripts more or less require the names of these
# files to be as listed above. When being archived, the filenames get a
# unique extension (describing the assimilation time) appended to them.
#
# The inflation file is essentially a duplicate of the model state ...
# it is slaved to a specific geometry. The initial files are created
# offline with values of unity. For the purpose of this script, they are
# thought to be the output of a previous assimilation, so they should be
# named something like prior_inflate_restart.YYYY-MM-DD-SSSSS
#
# The first inflation file can be created with 'fill_inflation_restart'
# which can be built in the usual DART manner.
#
# The strategy is to use the LATEST inflation file from the CESM 'rundir'.
# After an assimilation, the new inflation values/files will be moved to
# the CESM rundir to be used for subsequent assimilations. If the short-term
# archiver has worked correctly, only the LATEST files will available. Of
# course, it is not required to have short-term archiving turned on, so ...
#=========================================================================

set  MYSTRING = `grep inf_flavor input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  PRIOR_INF = $MYSTRING[2]
set  POSTE_INF = $MYSTRING[3]

set  MYSTRING = `grep inf_initial_from_restart input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  PRIOR_TF = `echo $MYSTRING[2] | tr '[:upper:]' '[:lower:]'`
set  POSTE_TF = `echo $MYSTRING[3] | tr '[:upper:]' '[:lower:]'`

# its a little tricky to remove both styles of quotes from the string.

set  MYSTRING = `grep inf_in_file_name input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set  PRIOR_INF_IFNAME = $MYSTRING[2]
set  POSTE_INF_IFNAME = $MYSTRING[3]

set  MYSTRING = `grep inf_out_file_name input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set  PRIOR_INF_OFNAME = $MYSTRING[2]
set  POSTE_INF_OFNAME = $MYSTRING[3]

set  MYSTRING = `grep inf_diag_file_name input.nml`
set  MYSTRING = `echo $MYSTRING | sed -e "s#[=,'\.]# #g"`
set  MYSTRING = `echo $MYSTRING | sed -e 's#"# #g'`
set  PRIOR_INF_DIAG = $MYSTRING[2]
set  POSTE_INF_DIAG = $MYSTRING[3]

# IFF we want PRIOR inflation:

if ( $PRIOR_INF > 0 ) then

   if ($PRIOR_TF == false) then
      echo "inf_flavor(1) = $PRIOR_INF, using namelist values."
   else
      # Look for the output from the previous assimilation
      (ls -rt1 ../pop_${PRIOR_INF_OFNAME}.* | tail -1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      # If one exists, use it as input for this assimilation
      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest ${PRIOR_INF_IFNAME}
      else
         echo "ERROR: Requested PRIOR inflation but specified no incoming inflation file."
         echo "ERROR: expected something like ../pop_${PRIOR_INF_OFNAME}.YYYY-MM-DD-SSSSS"
         exit -4
      endif

   endif
else
   echo "Prior Inflation           not requested for this assimilation."
endif

# POSTERIOR: We look for the 'newest' and use it - IFF we need it.

if ( $POSTE_INF > 0 ) then

   if ($POSTE_TF == false) then
      echo "inf_flavor(2) = $POSTE_INF, using namelist values."
   else

      # Look for the output from the previous assimilation
      (ls -rt1 ../pop_${POSTE_INF_OFNAME}.* | tail -1 >! latestfile) > & /dev/null
      set nfiles = `cat latestfile | wc -l`

      # If one exists, use it as input for this assimilation
      if ( $nfiles > 0 ) then
         set latest = `cat latestfile`
         ${LINK} $latest ${POSTE_INF_IFNAME}
      else
         echo "ERROR: Requested POSTERIOR inflation but specified no incoming inflation file."
         echo "ERROR: expected something like ../pop_${POSTE_INF_OFNAME}.YYYY-MM-DD-SSSSS"
         exit -5
      endif
   endif
else
   echo "Posterior Inflation       not requested for this assimilation."
endif

#=========================================================================
# Block 4: Convert N POP restart files to DART initial condition files.
# pop_to_dart is serial code, we can do all of these at the same time
# as long as we can have unique namelists for each of them.
#
# At the end of the block, we have DART initial condition files  filter_ics.[1-N]
# that came from pointer files ../rpointer.ocn_[1-N].restart
#
# REQUIRED DART namelist settings:
# &filter_nml:           restart_in_file_name    = 'filter_ics'
#                        restart_out_file_name   = 'filter_restart'
# &ensemble_manager_nml: single_restart_file_in  = '.false.'
# &pop_to_dart_nml:      pop_to_dart_output_file = 'dart_ics',
# &dart_to_pop_nml:      dart_to_pop_input_file  = 'dart_restart',
#                        advance_time_present    = .false.
#=========================================================================

echo "`date` -- BEGIN POP-TO-DART"

set member = 1
while ( ${member} <= ${ensemble_size} )

   # Each member will do its job in its own directory.
   # That way, we can do N of them simultaneously -

   set MYTEMPDIR = member_${member}
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   # make sure there are no old output logs hanging around
   $REMOVE output.${member}.pop_to_dart

   set OCN_RESTART_FILENAME = `printf ../../${MYCASE}.pop_%04d.r.${OCN_DATE_EXT}.nc  ${member}`
   set     OCN_NML_FILENAME = `printf ../../pop2_in_%04d  ${member}`
   set     DART_IC_FILENAME = `printf filter_ics.%04d     ${member}`
   set    DART_RESTART_FILE = `printf filter_restart.%04d ${member}`

   sed -e "s/dart_ics/..\/${DART_IC_FILENAME}/" \
       -e "s/dart_restart/..\/${DART_RESTART_FILE}/" < ../input.nml >! input.nml

   ${LINK} $OCN_RESTART_FILENAME pop.r.nc
   ${LINK}     $OCN_NML_FILENAME pop_in

   echo "starting pop_to_dart for member ${member} at "`date`
   ${EXEROOT}/pop_to_dart >! output.${member}.pop_to_dart &

   cd ..

   @ member++
end

wait

set nsuccess = `fgrep 'Finished ... at YYYY' member*/output.[0-9]*.pop_to_dart | wc -l`
if (${nsuccess} != ${ensemble_size}) then
   echo "ERROR ... DART died in 'pop_to_dart' ... ERROR"
   echo "ERROR ... DART died in 'pop_to_dart' ... ERROR"
   exit -6
endif

echo "`date` -- END POP-TO-DART for all ${ensemble_size} members."

#=========================================================================
# Block 5: Actually run the assimilation.
# Will result in a set of files : 'filter_restart.xxxx'
#
# DART namelist settings required:
# &filter_nml:           async                  = 0,
# &filter_nml:           adv_ens_command        = "./no_model_advance.csh",
# &filter_nml:           restart_in_file_name   = 'filter_ics'
# &filter_nml:           restart_out_file_name  = 'filter_restart'
# &filter_nml:           obs_sequence_in_name   = 'obs_seq.out'
# &filter_nml:           obs_sequence_out_name  = 'obs_seq.final'
# &filter_nml:           init_time_days         = -1,
# &filter_nml:           init_time_seconds      = -1,
# &filter_nml:           first_obs_days         = -1,
# &filter_nml:           first_obs_seconds      = -1,
# &filter_nml:           last_obs_days          = -1,
# &filter_nml:           last_obs_seconds       = -1,
# &ensemble_manager_nml: single_restart_file_in = .false.
#
#=========================================================================

# POP always needs a pop_in and a pop.r.nc to start.
# Lots of ways to get the filename

set OCN_RESTART_FILENAME = `head -1 ../rpointer.ocn_0001.restart`

${LINK} ../$OCN_RESTART_FILENAME pop.r.nc
${LINK} ../pop2_in_0001          pop_in

# On yellowstone, you can explore task layouts with the following:
if ( $?LSB_PJL_TASK_GEOMETRY ) then
   setenv ORIGINAL_LAYOUT "${LSB_PJL_TASK_GEOMETRY}"

   # setenv GEOMETRY_32_1NODE \
   #    "{(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31)}";
   # setenv LSB_PJL_TASK_GEOMETRY "${GEOMETRY_32_1NODE}"
endif

echo "`date` -- BEGIN FILTER"
${LAUNCHCMD} ${EXEROOT}/filter || exit -7
echo "`date` -- END FILTER"

if ( $?LSB_PJL_TASK_GEOMETRY ) then
   setenv LSB_PJL_TASK_GEOMETRY "${ORIGINAL_LAYOUT}"
endif

${MOVE} Prior_Diag.nc      ../pop_Prior_Diag.${OCN_DATE_EXT}.nc
${MOVE} Posterior_Diag.nc  ../pop_Posterior_Diag.${OCN_DATE_EXT}.nc
${MOVE} obs_seq.final      ../pop_obs_seq.${OCN_DATE_EXT}.final
${MOVE} dart_log.out       ../pop_dart_log.${OCN_DATE_EXT}.out

# Accomodate any possible inflation files
# 1) rename file to reflect current date
# 2) move to CENTRALDIR so the DART INFLATION BLOCK works next time and
#    that they can get archived.

foreach FILE ( ${PRIOR_INF_OFNAME} ${POSTE_INF_OFNAME} ${PRIOR_INF_DIAG} ${POSTE_INF_DIAG} )
   if ( -e ${FILE} ) then
      ${MOVE} ${FILE} ../pop_${FILE}.${OCN_DATE_EXT}
   else
      echo "No ${FILE} for ${OCN_DATE_EXT}"
   endif
end

#=========================================================================
# Block 6: Update the POP restart files ... simultaneously ...
#
# Each member will do its job in its own directory, which already exists
# and has the required input files remaining from 'Block 4'
#=========================================================================

echo "`date` -- BEGIN DART-TO-POP"
set member = 1
while ( $member <= $ensemble_size )

   cd member_${member}

   ${REMOVE} output.${member}.dart_to_pop

   echo "starting dart_to_pop for member ${member} at "`date`
   ${EXEROOT}/dart_to_pop >! output.${member}.dart_to_pop &

   cd ..

   @ member++
end

wait

set nsuccess = `fgrep 'Finished ... at YYYY' member*/output.[0-9]*.dart_to_pop | wc -l`
if (${nsuccess} != ${ensemble_size}) then
   echo "ERROR ... DART died in 'dart_to_pop' ... ERROR"
   echo "ERROR ... DART died in 'dart_to_pop' ... ERROR"
   exit -8
endif

echo "`date` -- END DART-TO-POP for all ${ensemble_size} members."

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

echo "`date` -- END POP_ASSIMILATE"

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

