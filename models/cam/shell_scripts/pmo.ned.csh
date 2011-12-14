#!/usr/bin/tcsh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id: assimilate.hopper.csh 5403 2011-11-16 15:20:05Z nancy $

# The FORCE options are not optional. 
# the VERBOSE options are useful for debugging.
set   MOVE = 'mv -fv'
set   COPY = 'cp -fv --preserve=timestamps'
set   LINK = 'ln -fvs'
set REMOVE = 'rm -fr'

# for perfect model we only have a single state
set ensemble_size = 1

# Create temporary working directory for the assimilation
set temp_dir = assimilate_dir
echo "temp_dir is $temp_dir"
mkdir -p $temp_dir
cd $temp_dir

#-------------------------------------------------------------------------
# Determine time of model state ... from file name of first member
# of the form "./${CASE}.cam_${ensemble_member}.r.2000-01-06-00000.nc"
#-------------------------------------------------------------------------

set FILE = `head -1 ../rpointer.atm_0001`
set FILE = $FILE:t
set FILE = $FILE:r
set MYCASE = $FILE:ar
set MODEL_DATE_EXT = `echo $FILE:e`
set MODEL_DATE     = `echo $FILE:e | sed -e "s#-# #g"`
set MODEL_YEAR     = $MODEL_DATE[1]
set MODEL_MONTH    = $MODEL_DATE[2]
set MODEL_DAY      = $MODEL_DATE[3]
set MODEL_SECONDS  = $MODEL_DATE[4]
set MODEL_HOUR     = `echo $MODEL_DATE[4] / 3600 | bc`

echo "valid time of model is $MODEL_YEAR $MODEL_MONTH $MODEL_DAY $MODEL_SECONDS (seconds)"
echo "valid time of model is $MODEL_YEAR $MODEL_MONTH $MODEL_DAY $MODEL_HOUR (hours)"

#-----------------------------------------------------------------------------
# Set variables containing various directory names where we will GET things
#-----------------------------------------------------------------------------

# FIXME: different for everyone
set DARTROOT = ${HOME}/devel
set DARTDIR = ${DARTROOT}/models/cam/work

# FIXME: different on hopper
set DART_OBS_DIR = UVT_skeleton_12H
set  OBSDIR = /scratch/scratchdirs/nscollin/Synthetic/${DART_OBS_DIR}

#-------------------------------------------------------------------------
# DART COPY BLOCK
# Populate a run-time directory with the bits needed to run DART 
#-------------------------------------------------------------------------

foreach FILE ( input_pmo.nml perfect_model_obs cam_to_dart dart_to_cam )
   if (  -e   ${DARTDIR}/${FILE} ) then
      ${COPY} ${DARTDIR}/${FILE} .
   else
      echo "DART required file ${DARTDIR}/${FILE} not found ... ERROR"
      exit 1
   endif
end

${COPY} input_pmo.nml input.nml

# surface height data
${COPY} $HOME/cam_phis.nc .

#-------------------------------------------------------------------------
# Block 1: convert 1 cam restart files to DART initial conditions file
#
# At the end of the block, we have DART restart files  perfect_ics
# that came from pointer files ../rpointer.atm_[1]
#
# DART namelist settings appropriate/required:
# &perfect_model_obs_nml: restart_in_file_name   = 'perfect_ics'
# &perfect_model_obs_nml: restart_out_file_name  = 'perfect_restart'
# &ensemble_manager_nml: single_restart_file_in  = '.true.'
# &cam_to_dart_nml:      cam_to_dart_output_file = 'perfect_ics',
#-------------------------------------------------------------------------

# only a single member, no need for subdirs

set member = 1

set POINTER_FILENAME = `printf rpointer.atm_%04d ${member}`
set MODEL_RESTART_FILENAME = `head -1 ../${POINTER_FILENAME}`
set MODEL_INITIAL_FILENAME = `echo ${MODEL_RESTART_FILENAME} | sed "s#\.r\.#\.i\.#"`
${LINK} ../$MODEL_INITIAL_FILENAME caminput.nc

# TJH can we use a .h0. file instead of some arbitrary cam_phis.nc

set DART_IC_FILE = "filter_ics"

echo "starting cam_to_dart for member ${member} at "`date`
./cam_to_dart >! output.${member}.cam_to_dart 
echo "finished cam_to_dart for member ${member} at "`date`


#-------------------------------------------------------------------------
# Block 2: Actually run the assimilation.
# Will result in a set of files : 'perfect_restart'
#
# DART namelist settings required:
# &perfect_model_obs_nml:  async                  = 0,
# &perfect_model_obs_nml:  adv_ens_command        = "no_model_advance",
# &perfect_model_obs_nml:  restart_in_file_name   = 'perfect_ics'
# &perfect_model_obs_nml:  restart_out_file_name  = 'perfect_restart'
# &perfect_model_obs_nml:  obs_sequence_in_name   = 'obs_seq.in'
# &perfect_model_obs_nml:  obs_sequence_out_name  = 'obs_seq.out'
# &perfect_model_obs_nml:  init_time_days         = -1,
# &perfect_model_obs_nml:  init_time_seconds      = -1,
# &perfect_model_obs_nml:  first_obs_days         = -1,
# &perfect_model_obs_nml:  first_obs_seconds      = -1,
# &perfect_model_obs_nml:  last_obs_days          = -1,
# &perfect_model_obs_nml:  last_obs_seconds       = -1,
# &ensemble_manager_nml:   single_restart_file_in = .true.
#
#-------------------------------------------------------------------------

# cam always needs a cam_initial.nc and a cam_history.nc to start.

set MODEL_RESTART_FILENAME = `head -1 ../rpointer.atm_0001`
set MODEL_INITIAL_FILENAME = `echo ${MODEL_RESTART_FILENAME} | sed "s#\.r\.#\.i\.#"`
set MODEL_HISTORY_FILENAME = `echo ${MODEL_RESTART_FILENAME} | sed "s#\.r\.#\.h0\.#"`

${LINK} ../$MODEL_INITIAL_FILENAME caminput.nc
#${LINK} ../$MODEL_RESTART_FILENAME cam_restart.nc
#${LINK} ../$MODEL_HISTORY_FILENAME cam_history.nc

# Determine proper observation sequence file.

set OBSFNAME = `printf obs_seq${MODEL_YEAR}${MODEL_MONTH}${MODEL_DAY}%02d ${MODEL_HOUR}`
set OBS_FILE = ${OBSDIR}/${OBSFNAME} 

${LINK} ${OBS_FILE} obs_seq.in


# the work happens here

echo starting perfect_model_obs executable now with 1 task
./perfect_model_obs || exit 7


${MOVE} obs_seq.out        ../obs_seq.${MODEL_DATE_EXT}.out
${MOVE} dart_log.out       ../dart_log.${MODEL_DATE_EXT}.out
${MOVE} True_State.nc      ../True_State.${MODEL_DATE_EXT}.nc
${MOVE} Prior_Diag.nc      ../Prior_Diag.${MODEL_DATE_EXT}.nc


#-------------------------------------------------------------------------
# Block 3: Update the cam restart file
#
# DART namelist settings required:
# &perfect_model_obs_nml: restart_out_file_name  = 'perfect_restart'
# &ensemble_manager_nml:  single_restart_file_in = '.true.'
# &dart_to_cam_nml:       dart_to_cam_input_file = 'perfect_restart',
# &dart_to_cam_nml:       advance_time_present   = .false.
# &atm_in_xxxx:ncdata = 'cam_initial_x.nc'
#-------------------------------------------------------------------------

set member = 1

set DART_RESTART_FILE = "perfect_restart"

set ATM_POINTER_FILENAME = `printf rpointer.atm_%04d ${member}`
set LND_POINTER_FILENAME = `printf rpointer.lnd_%04d ${member}`
set ICE_POINTER_FILENAME = `printf rpointer.ice_%04d ${member}`

set ATM_RESTART_FILENAME = `head -1 ../${ATM_POINTER_FILENAME}`
set LND_RESTART_FILENAME = `echo ${ATM_RESTART_FILENAME} | sed "s#\.cam_#\.clm2_#"`
set ICE_RESTART_FILENAME = `echo ${ATM_RESTART_FILENAME} | sed "s#\.cam_#\.cice_#"`

set ATM_INITIAL_FILENAME = `echo ${ATM_RESTART_FILENAME} | sed "s#\.r\.#\.i\.#"`

#  set ATM_HISTORY_FILENAME = `echo ${ATM_RESTART_FILENAME} | sed "s#\.r\.#\.h0\.#"`
#  ${LINK} ../$ATM_RESTART_FILENAME cam_restart.nc
#  ${LINK} ../$ATM_HISTORY_FILENAME cam_history.nc
   ${LINK} ../$ATM_INITIAL_FILENAME caminput.nc

echo "starting dart_to_cam for member ${member} at "`date`
./dart_to_cam >! output.${member}.dart_to_cam
echo "finished dart_to_cam for member ${member} at "`date`

# The initial filenames are static and come from the atm_in_xxxx namelist.
# We must copy the updated initial files to the static names.

${COPY} ../$ATM_INITIAL_FILENAME ../cam_initial_${member}.nc
${COPY} ../$LND_RESTART_FILENAME ../clm_restart_${member}.nc
${COPY} ../$ICE_RESTART_FILENAME ../ice_restart_${member}.nc


#-------------------------------------------------------------------------
# Now that everything is staged, we have to communicate the current 
# model time to the drv_in&seq_timemgr_inparm namelist 
# which is built from CASEROOT/user_nl_drv by the *.run script
#-------------------------------------------------------------------------

ex ${CASEROOT}/Buildconf/cpl.buildnml.csh << ex_end
g; start_ymd;s;=[ ]*.*;= ${MODEL_YEAR}${MODEL_MONTH}${MODEL_DAY};
g; start_tod;s;=[ ]*.*;= $MODEL_SECONDS;
wq
ex_end

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

exit 0

# <next few lines under version control, do not edit>
# $URL: https://proxy.subversion.ucar.edu/DAReS/DART/branches/development/models/cam/shell_scripts/assimilate.hopper.csh $
# $Revision: 5403 $
# $Date: 2011-11-16 07:20:05 -0800 (Wed, 16 Nov 2011) $

