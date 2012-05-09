#!/usr/local/bin/tcsh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$

# The FORCE options are not optional.
# the VERBOSE options are useful for debugging.
set   MOVE = '/usr/local/bin/mv -fv'
set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
set   LINK = '/usr/local/bin/ln -fvs'
set REMOVE = '/usr/local/bin/rm -fr'

set ensemble_size = ${NINST_LND}

# Create temporary working directory for the assimilation
set temp_dir = assimilate_dir
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
# of the form "./${CASE}.clm2_${ensemble_member}.r.2000-01-06-00000.nc"
#-------------------------------------------------------------------------

set FILE = `head -1 ../rpointer.lnd_0001`
set FILE = $FILE:t
set FILE = $FILE:r
set MODEL_DATE_EXT = `echo $FILE:e`
set MODEL_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
set MODEL_YEAR     = $MODEL_DATE[1]
set MODEL_MONTH    = $MODEL_DATE[2]
set MODEL_DAY      = $MODEL_DATE[3]
set MODEL_SECONDS  = $MODEL_DATE[4]

echo "valid time of model is $MODEL_YEAR $MODEL_MONTH $MODEL_DAY $MODEL_SECONDS"

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------

set DARTDIR = ${HOME}/DART/models/clm/work

set DART_OBS_DIR = ${MODEL_YEAR}${MODEL_MONTH}
set  OBSDIR = /ptmp/dart/Obs_sets/clm/${DART_OBS_DIR}

#-------------------------------------------------------------------------
# DART COPY BLOCK
# Populate a run-time directory with the bits needed to run DART
# Either get them from the CCSM 'run' directory or some stock repository
# The grid files are absolute paths ... so they need not move.
#-------------------------------------------------------------------------

foreach FILE ( input.nml filter clm_to_dart dart_to_clm )

   if ( -e ${CASEROOT}/${FILE} ) then
      ${COPY}   ${CASEROOT}/${FILE} .
   else if ( -e    ../${FILE} ) then
      ${COPY} ../${FILE} .
   else if ( -e ${DARTDIR}/${FILE} ) then
      ${COPY}   ${DARTDIR}/${FILE} .
   else
      echo "DART required file $FILE not found ... ERROR"
      exit 1
   endif

end

#-------------------------------------------------------------------------
# This is the file for the sampling error correction.
# Each ensemble size has its own file.
# It is static - it does not need to be archived, etc.
# It is only needed if
# input.nml:&assim_tools_nml:sampling_error_correction = .true.,
#-------------------------------------------------------------------------

set SAMP_ERR_FILE = ${DARTDIR}/system_simulation/final_full.${ensemble_size}

if ( -e ${SAMP_ERR_FILE}/ ) then
   ${COPY} ${SAMP_ERR_FILE} .
else
   echo "WARNING: no sampling error correction file for this ensemble size."
   echo "warning: looking for system_simulation/final_full.${ensemble_size}"
endif

#-------------------------------------------------------------------------
# DART INFLATION BLOCK
# This file is only relevant if 'inflation' is turned on -
# i.e. if inf_flavor(1) /= 0 - AND we are in a 'restart' mode.
#
# filter_nml
# inf_flavor                  = 2,                       0,
# inf_initial_from_restart    = .true.,                  .false.,
# inf_in_file_name            = 'prior_inflate_ics',     'post_inflate_ics',
#
# This is a 'test' configuration for this script. We are simply
# assuming that the namelist values are set such that we need this file,
# and that it is called 'prior_inflate_ics'. Since the inflation file is
# essentially a duplicate of the model state ... it is slaved to a specific
# geometry. I created the file offline.  The inflation values are all unity.
#
# The strategy is to use the LATEST inflation file from CENTRALDIR if one exists -
#
# After an assimilation, the output file will be copied back to CENTRALDIR
# to be used for subsequent assimilations.
#-------------------------------------------------------------------------

foreach FILE ( prior post )

   # These files may or may not exist. This causes some complexity.
   # So - we look for the 'newest' and use it. And Pray.

   (ls -rt1 ../${FILE}_inflate.*.restart.* | tail -1 >! latestfile) > & /dev/null
   set nfiles = `cat latestfile | wc -l`

   if ( $nfiles > 0 ) then
      set latest = `cat latestfile`
      ${LINK} $latest ${FILE}_inflate_ics
   else
      # MUST HAVE inf_initial_from_restart = .false.
      echo "WARNING: no incoming ${FILE}_inflate.YYYY-MM-DD-00000.restart.endiansuffix"
   endif

end

#-------------------------------------------------------------------------
# Block 1: convert N clm restart files to DART initial conditions file(s)
# clm_to_dart is serial code, we can do all of these at the same time
# and just wait for them to finish IFF it were not for the fact we'd have
# to have unique namelists for all of them.
#
# At the end of the block, we have DART restart files  filter_ics.[1-N]
# that came from pointer files ../rpointer.lnd_[1-N]
#
# DART namelist settings appropriate/required:
# &filter_nml:           restart_in_file_name    = 'filter_ics'
# &ensemble_manager_nml: single_restart_file_in  = '.false.'
# &clm_to_dart_nml:      clm_to_dart_output_file = 'dart.ud',
#-------------------------------------------------------------------------

set member = 1
while ( ${member} <= ${ensemble_size} )

   # Each member will do its job in its own directory.
   # That way, we can do N of them simultaneously -
   # they all read their OWN 'input.nml' ... the output
   # filenames must inserted into the appropriate input.nml

   set MYTEMPDIR = member_${member}
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   set POINTER_FILENAME = `printf rpointer.lnd_%04d ${member}`
   set MODEL_RESTART_FILENAME = `head -1 ../../${POINTER_FILENAME}`
   set MODEL_HISTORY_FILENAME = `echo ${MODEL_RESTART_FILENAME} | sed "s/\.r\./\.h0\./"`
   ${LINK} ../../$MODEL_RESTART_FILENAME clm_restart.nc
   ${LINK} ../../$MODEL_HISTORY_FILENAME clm_history.nc

   # patch the CLM restart files to ensure they have the proper
   # _FillValue and missing_value attributes.
   ncatted -O -a    _FillValue,frac_sno,o,d,1.0e+36   clm_restart.nc
   ncatted -O -a missing_value,frac_sno,o,d,1.0e+36   clm_restart.nc
   ncatted -O -a    _FillValue,DZSNO,o,d,1.0e+36      clm_restart.nc
   ncatted -O -a missing_value,DZSNO,o,d,1.0e+36      clm_restart.nc
   ncatted -O -a    _FillValue,H2OSOI_LIQ,o,d,1.0e+36 clm_restart.nc
   ncatted -O -a missing_value,H2OSOI_LIQ,o,d,1.0e+36 clm_restart.nc
   ncatted -O -a    _FillValue,H2OSOI_ICE,o,d,1.0e+36 clm_restart.nc
   ncatted -O -a missing_value,H2OSOI_ICE,o,d,1.0e+36 clm_restart.nc
   ncatted -O -a    _FillValue,T_SOISNO,o,d,1.0e+36   clm_restart.nc
   ncatted -O -a missing_value,T_SOISNO,o,d,1.0e+36   clm_restart.nc

   # the slash in the filename screws up 'sed' ... unless
   set DART_IC_FILE = `printf ..\\/filter_ics.%04d ${member}`

   sed -e "s/dart.ud/${DART_IC_FILE}/" < ../input.nml >! input.nml

   echo "starting clm_to_dart for member ${member} at "`date`
   ../clm_to_dart >! output.${member}.clm_to_dart &
   echo "finished clm_to_dart for member ${member} at "`date`

   cd ..

   @ member++
end

wait

#-------------------------------------------------------------------------
# Block 2: Actually run the assimilation.
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
# &ensemble_manager_nml: single_restart_file_in = '.false.'
#
#-------------------------------------------------------------------------

# clm always needs a clm_restart.nc, and a clm_history.nc to start.

set MODEL_RESTART_FILENAME = `head -1 ../rpointer.lnd_0001`
set MODEL_HISTORY_FILENAME = `echo ${MODEL_RESTART_FILENAME} | sed "s/\.r\./\.h0\./"`

${LINK} ../$MODEL_RESTART_FILENAME clm_restart.nc
${LINK} ../$MODEL_HISTORY_FILENAME clm_history.nc

# Determine proper observation sequence file.

set OBSFNAME = obs_seq.0Z.${MODEL_YEAR}${MODEL_MONTH}${MODEL_DAY}
set OBS_FILE = ${OBSDIR}/${OBSFNAME}

${LINK} ${OBS_FILE}   obs_seq.out

# FIXME: special for trying out non-monotonic task layouts.
# FIXME setenv ORG_PATH "${PATH}"
# FIXME setenv LSF_BINDIR /contrib/lsf/tgmpatch
# FIXME setenv PATH ${LSF_BINDIR}:${PATH}
# FIXME setenv ORG_TASK_GEOMETRY "${LSB_PJL_TASK_GEOMETRY}"

# layout: flat
# FIXME setenv NANCY_GEOMETRY_54_1NODE \
# FIXME 	"{(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53)}";

# FIXME setenv LSB_PJL_TASK_GEOMETRY "${NANCY_GEOMETRY_54_1NODE}"

# FIXME which mpirun.lsf

mpirun ./filter || exit 2

${MOVE} Prior_Diag.nc      ../Prior_Diag.${MODEL_DATE_EXT}.nc
${MOVE} Posterior_Diag.nc  ../Posterior_Diag.${MODEL_DATE_EXT}.nc
${MOVE} obs_seq.final      ../obs_seq.${MODEL_DATE_EXT}.final
${MOVE} dart_log.out       ../dart_log.${MODEL_DATE_EXT}.out

# Accomodate any possible inflation files

foreach INFLATION ( prior post )

   if ( -e ${INFLATION}_inflate_restart ) then
      # 1) rename file to reflect current date
      # 2) move to CENTRALDIR so the DART INFLATION BLOCK works next time

      ${MOVE} ${INFLATION}_inflate_restart ../${INFLATION}_inflate.${MODEL_DATE_EXT}.restart.be
   else
      echo "No ${INFLATION}_inflate_restart for ${MODEL_DATE_EXT}"
   endif

   if ( -e ${INFLATION}_inflate_diag ) then
      ${MOVE} ${INFLATION}_inflate_diag ../${INFLATION}_inflate.${MODEL_DATE_EXT}.diag
   else
      echo "No ${INFLATION}_inflate_diag for ${MODEL_DATE_EXT}"
   endif

end

# FIXME: special for trying out non-monotonic task layouts.
# FIXME setenv PATH "${ORG_PATH}"
# FIXME setenv LSB_PJL_TASK_GEOMETRY "${ORG_TASK_GEOMETRY}"

#-------------------------------------------------------------------------
# Block 3: Update the clm restart files ... simultaneously ...
#
# DART namelist settings required:
# &filter_nml:           restart_out_file_name  = 'filter_restart'
# &ensemble_manager_nml: single_restart_file_in = '.false.'
# &dart_to_clm_nml:      dart_to_clm_input_file = 'dart.ic',
# &dart_to_clm_nml:      advance_time_present   = .false.
#-------------------------------------------------------------------------

set member = 1
while ( ${member} <= ${ensemble_size} )

   # Each member will do its job in its own directory.
   # That way, we can do N of them simultaneously -
   # they all read their OWN 'input.nml' ... the output
   # filenames must inserted into the appropriate input.nml

   set MYTEMPDIR = member_${member}
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   set DART_RESTART_FILE = `printf filter_restart.%04d ${member}`
   ${LINK} ../$DART_RESTART_FILE dart.ic

   set POINTER_FILENAME = `printf rpointer.lnd_%04d ${member}`
   set MODEL_RESTART_FILENAME = `head -1 ../../${POINTER_FILENAME}`
   set MODEL_HISTORY_FILENAME = `echo ${MODEL_RESTART_FILENAME} | sed "s/\.r\./\.h0\./"`
   ${LINK} ../../$MODEL_RESTART_FILENAME clm_restart.nc
   ${LINK} ../../$MODEL_HISTORY_FILENAME clm_history.nc

   echo "starting dart_to_clm for member ${member} at "`date`
   ../dart_to_clm >! output.${member}.dart_to_clm &
   echo "finished dart_to_clm for member ${member} at "`date`

   cd ..

   @ member++
end

wait

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

ls -lrt

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

