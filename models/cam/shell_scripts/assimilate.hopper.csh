#!/usr/bin/tcsh
# FIXME: where it might be on bluefire but not hopper:
#!/usr/local/bin/tcsh
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$

# The FORCE options are not optional. 
# the VERBOSE options are useful for debugging.
set   MOVE = 'mv -fv'
set   COPY = 'cp -fv --preserve=timestamps'
set   LINK = 'ln -fvs'
set REMOVE = 'rm -fr'
# FIXME: not on hopper:
#set   MOVE = '/usr/local/bin/mv -fv'
#set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
#set   LINK = '/usr/local/bin/ln -fvs'
#set REMOVE = '/usr/local/bin/rm -fr'

set ensemble_size = ${NINST_ATM}

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
set DARTDIR = ${HOME}/devel/models/cam/work

# FIXME: different on hopper
set DART_OBS_DIR = ${MODEL_YEAR}${MODEL_MONTH}_6H
set  OBSDIR = /scratch/scratchdirs/nscollin/ACARS/${DART_OBS_DIR}

#-------------------------------------------------------------------------
# DART COPY BLOCK
# Populate a run-time directory with the bits needed to run DART 
#-------------------------------------------------------------------------

foreach FILE ( input.nml filter cam_to_dart dart_to_cam )
   if (  -e   ${DARTDIR}/${FILE} ) then
      ${COPY} ${DARTDIR}/${FILE} .
   else
      echo "DART required file ${DARTDIR}/${FILE} not found ... ERROR"
      exit 1
   endif
end

# FIXME: doesn't exist on hopper
#${COPY} /glade/proj3/DART/raeder/FV1deg_4.0/cam_phis.nc .
${COPY} $HOME/cam_phis.nc .

# Make sure the DART ensemble size matches CESM instances

ex input.nml <<ex_end
g;ens_size ;s;= .*;= $ensemble_size;
g;num_output_state_members ;s;= .*;= $ensemble_size;
g;num_output_obs_members ;s;= .*;= $ensemble_size;
wq
ex_end


#-------------------------------------------------------------------------
# This is the file for the sampling error correction.
# Each ensemble size has its own file.
# It is static - it does not need to be archived, etc.
# It is only needed if 
# input.nml:&assim_tools_nml:sampling_error_correction = .true.,
#-------------------------------------------------------------------------

# FIXME: only do this is SEC is true
set SAMP_ERR_FILE = ${DARTDIR}/system_simulation/final_full.${ensemble_size}

if ( -e ${SAMP_ERR_FILE}/ ) then
   ${COPY} ${SAMP_ERR_FILE} .
else
   echo "WARNING: no sampling error correction file for this ensemble size."
   echo "warning: looking for ${SAMP_ERR_FILE}"
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
# Block 1: convert N cam restart files to DART initial conditions file(s)
# cam_to_dart is serial code, we can do all of these at the same time
# and just wait for them to finish IFF it were not for the fact we'd have
# to have unique namelists for all of them.
#
# At the end of the block, we have DART restart files  filter_ic_old.[1-N]
# that came from pointer files ../rpointer.atm_[1-N]
#
# DART namelist settings appropriate/required:
# &filter_nml:           restart_in_file_name    = 'filter_ic_old'
# &ensemble_manager_nml: single_restart_file_in  = '.false.'
# &cam_to_dart_nml:      cam_to_dart_output_file = 'dart_ics',
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

   set POINTER_FILENAME = `printf rpointer.atm_%04d ${member}`
   set MODEL_RESTART_FILENAME = `head -1 ../../${POINTER_FILENAME}`
   set MODEL_INITIAL_FILENAME = `echo ${MODEL_RESTART_FILENAME} | sed "s#\.r\.#\.i\.#"`
   ${LINK} ../../$MODEL_INITIAL_FILENAME caminput.nc
   ${LINK} ../cam_phis.nc .

   # TJH can we use a .h0. file instead of some arbitrary cam_phis.nc

   set DART_IC_FILE = `printf ../filter_ic_old.%04d ${member}`

   sed -e "s#dart_ics#${DART_IC_FILE}#" < ../input.nml >! input.nml

   echo "starting cam_to_dart for member ${member} at "`date`
   ../cam_to_dart >! output.${member}.cam_to_dart &
   echo "finished cam_to_dart for member ${member} at "`date`

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
# &filter_nml:           restart_in_file_name   = 'filter_ic_old'
# &filter_nml:           restart_out_file_name  = 'filter_ic_new'
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

${LINK} ${OBS_FILE} obs_seq.out

# special for trying out non-monotonic task layouts.
# FIXME: HOPPER does not use LSF, so none of these are defined
#setenv ORG_PATH "${PATH}"
#setenv LSF_BINDIR /contrib/lsf/tgmpatch
#setenv PATH ${LSF_BINDIR}:${PATH}
#setenv ORG_TASK_GEOMETRY "${LSB_PJL_TASK_GEOMETRY}"

# FIXME: hopper uses aprun not mpirun or mpirun.lsf
#if (bluefire) then
#which mpirun.lsf
#
#mpirun.lsf ./filter || exit 2
# else
# FIXME: this needs to come from the calling script
setenv NTASKS 24
aprun -n $NTASKS ./filter || exit 2
#endif


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

# special for trying out non-monotonic task layouts.
# FIXME:  not for hopper
#setenv PATH "${ORG_PATH}"
#setenv LSB_PJL_TASK_GEOMETRY "${ORG_TASK_GEOMETRY}"

#-------------------------------------------------------------------------
# Block 3: Update the cam restart files ... simultaneously ...
#
# DART namelist settings required:
# &filter_nml:           restart_out_file_name  = 'filter_ic_new'
# &ensemble_manager_nml: single_restart_file_in = '.false.'
# &dart_to_cam_nml:      dart_to_cam_input_file = 'temp_ic',
# &dart_to_cam_nml:      advance_time_present   = .false.
# &atm_in_xxxx:ncdata = 'cam_initial_x.nc'
#-------------------------------------------------------------------------

set member = 1
while ( ${member} <= ${ensemble_size} )

   # Each member will do its job in its own directory.
   # Cannot do these simultaneously -

   set MYTEMPDIR = member_${member}
   mkdir -p $MYTEMPDIR
   cd $MYTEMPDIR

   set DART_RESTART_FILE = `printf filter_ic_new.%04d ${member}`
   ${LINK} ../$DART_RESTART_FILE temp_ic

   set ATM_POINTER_FILENAME = `printf rpointer.atm_%04d ${member}`
   set LND_POINTER_FILENAME = `printf rpointer.lnd_%04d ${member}`
   set ICE_POINTER_FILENAME = `printf rpointer.ice_%04d ${member}`

   set ATM_RESTART_FILENAME = `head -1 ../../${ATM_POINTER_FILENAME}`
   set LND_RESTART_FILENAME = `echo ${ATM_RESTART_FILENAME} | sed "s#\.cam_#\.clm2_#"`
   set ICE_RESTART_FILENAME = `echo ${ATM_RESTART_FILENAME} | sed "s#\.cam_#\.cice_#"`

   set ATM_INITIAL_FILENAME = `echo ${ATM_RESTART_FILENAME} | sed "s#\.r\.#\.i\.#"`

#  set ATM_HISTORY_FILENAME = `echo ${ATM_RESTART_FILENAME} | sed "s#\.r\.#\.h0\.#"`
#  ${LINK} ../../$ATM_RESTART_FILENAME cam_restart.nc
#  ${LINK} ../../$ATM_HISTORY_FILENAME cam_history.nc
   ${LINK} ../../$ATM_INITIAL_FILENAME caminput.nc

   echo "starting dart_to_cam for member ${member} at "`date`
   ../dart_to_cam >! output.${member}.dart_to_cam
   echo "finished dart_to_cam for member ${member} at "`date`

   # The initial filenames are static and come from the atm_in_xxxx namelist.
   # We must rename the updated initial files to the static names.

   ${MOVE} ../../$ATM_INITIAL_FILENAME ../../cam_initial_${member}.nc
   ${MOVE} ../../$LND_RESTART_FILENAME ../../clm_restart_${member}.nc
   ${MOVE} ../../$ICE_RESTART_FILENAME ../../ice_restart_${member}.nc

   # This is a safety precaution. We want to make sure CESM uses the static
   # files for initialization as opposed to the pointer files and their
   # contents. Basically, we are preventing CESM from using restarts instead
   # of the initial files. Move all of the files for this instance to
   # this directory

   set INSTANCE = `printf %04d ${member}`
   ${MOVE} ../../$ATM_POINTER_FILENAME         .
   ${MOVE} ../../$LND_POINTER_FILENAME         .
   ${MOVE} ../../$ICE_POINTER_FILENAME         .
   ${MOVE} ../../${MYCASE}.*_${INSTANCE}.*.nc  .

   cd ..

   @ member++
end

wait

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

mkdir safehouse
${MOVE} *pointer* FD* safehouse

#-------------------------------------------------------------------------
# Cleanup
#-------------------------------------------------------------------------

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

