#!/bin/csh
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$

# ---------------------
# Purpose
# ---------------------
#
# This script is designed to configure and build a multi-instance CESM model
# that has XXX,YYY,ZZZZ as active components
# and will use DART to assimilate observations at regular intervals.
# This script does not build DART. It works best if the appropriate DART
# executables have been built, however.
#
# This script relies heavily on the information in:
# http://www.cesm.ucar.edu/models/cesm1.1/cesm/doc/usersguide/book1.html
#
# ---------------------
# How to set up the script
# ---------------------
# -- Either edit and run this script in the $DART/models/CESM/shell_scripts
#    directory where it first appears,
#    or copy it to somewhere that it will be preserved and run it there.
#    It will create a 'case' directory, where the model will be built,
#    and an execution directory, where each forecast and assimilation will
#    take place.  The short term archiver will use a third directory for
#    storage of model output until it can be moved to long term storage (HPSS)
# -- Examine the whole script to identify things to change for your experiments.
# -- Provide any initial files needed by your run:
#       inflation
#       sampling error correction
#       CESM initial ensemble
#       ...
# -- Run this script.
# -- Edit the DART input.nml that appears in the $CASEROOT directory.
# -- Submit the job using $CASEROOT/${case}.submit
#
# ---------------------
# Important features
# ---------------------
#
# If you want to change something in your case other than the runtime
# settings, it is safest to delete everything and start the run from scratch.
# For the brave, read
#
# http://www.cesm.ucar.edu/models/cesm1.1/cesm/doc/usersguide/x1142.html
#
# and you may be able to salvage something with
# ./cesm_setup -clean
# ./cesm_setup
# ./${case}.clean_build
# ./${case}.build
#
# ==============================================================================
# ====  Set case options
# ==============================================================================

# case will be used many ways;
#    directory and file names, both locally and on HPSS, and
#    script names; so consider it's length and information content.
# num_instances:  Number of ensemble members

setenv case                 cesm_startup
setenv compset              B_2000_CAM5
setenv resolution           0.9x1.25_gx1v6
setenv cesmtag              cesm1_1_1
setenv num_instances        2

# ==============================================================================
# define machines and directories
#
# mach            Computer name
# cesm_datadir    Root path of the public CESM data files
# cesmroot        Location of the cesm code base
#                 For cesm1_1 on yellowstone
# DARTroot        Location of DART code tree.
#                    Executables, scripts and input in $DARTroot/models/dev/...
# caseroot        Your (future) cesm case directory, where this CESM+DART will be built.
#                    Preferably not a frequently scrubbed location.
#                    This script will delete any existing caseroot, so this script,
#                    and other useful things should be kept elsewhere.
# rundir          (Future) Run-time directory; scrubbable, large amount of space needed.
# exeroot         (Future) directory for executables - scrubbable, large amount of space needed.
# archdir         (Future) Short-term archive directory
#                    until the long-term archiver moves it to permanent storage.
# ==============================================================================

setenv mach         yellowstone
setenv cesm_datadir /glade/p/cesm/cseg/inputdata
setenv cesmroot     /glade/p/cesm/cseg/collections/$cesmtag
setenv caseroot     /glade/p/work/${USER}/cases/${case}
setenv exeroot      /glade/scratch/${USER}/${case}/bld
setenv rundir       /glade/scratch/${USER}/${case}/run
setenv archdir      /glade/scratch/${USER}/archive/${case}

setenv DARTroot     /glade/u/home/${USER}/svn/DART/dev

set RTM_stagedir = /glade/scratch/thoar/DART_POP_RESTARTS/2004-01-01-00000
set CLM_stagedir = /glade/scratch/thoar/DART_POP_RESTARTS/CLM_2004-01-01-00000/cesm_test
set CAM_stagedir = /glade/p/cesm/cseg/inputdata/atm/cam/inic/fv
set POP_stagedir = /glade/p/work/aliciak/DART_IC/CCSM4_ensembles/rest/2004-01-01-00000

# ==============================================================================
# configure settings ... run_startdate format is yyyy-mm-dd
# ==============================================================================

setenv refyear     2004
setenv refmon      01
setenv refday      01
setenv run_reftod  00000
setenv run_refdate $refyear-$refmon-$refday

# ==============================================================================
# runtime settings --  How many assimilation steps will be done after this one
#
# stop_option   Units for determining the forecast length between assimilations
# stop_n        Number of time units in the forecast
# assim_n       Number of time units between assimilations
# ==============================================================================

setenv resubmit      0
setenv stop_option   nhours
setenv stop_n        72
setenv assim_n       24

# ==============================================================================
# job settings
#
# timewall   can be changed during a series by changing the ${case}.run
# queue      can be changed during a series by changing the ${case}.run
#
# TJH: How many T62_gx1v6 CESM instances can fit on 1 node?
# ==============================================================================

setenv ACCOUNT      P8685nnnn
setenv timewall     0:30
setenv queue        regular
setenv ptile        15

# ==============================================================================
# set these standard commands based on the machine you are running on.
# ==============================================================================

switch ("`hostname`")
   case be*:
      # NCAR "bluefire"
      # The FORCE options are not optional.
      # the VERBOSE options are useful for debugging.
      set   MOVE = '/usr/local/bin/mv -fv'
      set   COPY = '/usr/local/bin/cp -fv --preserve=timestamps'
      set   LINK = '/usr/local/bin/ln -fvs'
      set REMOVE = '/usr/local/bin/rm -fr'

   breaksw
   default:
      # NERSC "hopper", NWSC "yellowstone"
      set   MOVE = '/bin/mv -fv'
      set   COPY = '/bin/cp -fv --preserve=timestamps'
      set   LINK = '/bin/ln -fvs'
      set REMOVE = '/bin/rm -fr'
      set nonomatch

   breaksw
endsw

# ==============================================================================
# some simple error checking before diving into the work
# ==============================================================================

# make sure these directories exist
set musthavedirs = "cesm_datadir cesmroot DARTroot"
foreach VAR ( $musthavedirs )
   # VAR is the shell variable name, DIR is the value
   set DIR = `eval echo \${$VAR}`
   if ( ! -d $DIR ) then
      echo "ERROR: directory '$DIR' not found"
      echo " In the setup script check the setting of: $VAR "
      exit -1
   endif
end

# make sure there is a filter in these dirs
set musthavefiles = "cam POP clm"
foreach MODEL ( $musthavefiles )
   set target = $DARTroot/models/$MODEL/work/filter
   if ( ! -x $target ) then
      echo "ERROR: executable file 'filter' not found"
      echo " Looking for: $target "
      echo " Make sure all DART assimilation executables have "
      echo " been compiled before running this setup script."
      exit -1
   endif
end

# ==============================================================================
# Create the case.
#
# For list of the pre-defined cases: ./create_newcase -list
# To create a variant case, see the CESM documentation and carefully
# incorporate any needed changes into this script.
# ==============================================================================

   echo "removing old files from ${caseroot}"
   echo "removing old files from ${exeroot}"
   echo "removing old files from ${rundir}"
   ${REMOVE} ${caseroot}
   ${REMOVE} ${exeroot}
   ${REMOVE} ${rundir}
   ${cesmroot}/scripts/create_newcase -case ${caseroot} -mach ${mach} \
                   -res ${resolution} -compset ${compset}

   if ( $status != 0 ) then
      echo "ERROR: Case could not be created."
      exit -1
   endif

# ==============================================================================
# Configure the case - this creates the CASEROOT directory.
# ==============================================================================

cd ${caseroot}

# Save a copy for debug purposes
foreach FILE ( *xml )
   if ( ~ -e        ${FILE}.original ) then
      ${COPY} $FILE ${FILE}.original
   endif
end

   # This is only for the purpose of debugging the code.
   # A more efficient layout must be found
   @ atm_pes = $ptile * $num_instances * 4
   @ ocn_pes = $ptile * $num_instances * 4
   @ lnd_pes = $ptile * $num_instances * 4
   @ ice_pes = $ptile * $num_instances * 1
   @ glc_pes = $ptile * $num_instances
   @ rof_pes = $ptile * $num_instances
   @ cpl_pes = $ptile * 4

   @ glc_root = $lnd_pes + $ice_pes
   @ rof_root = $lnd_pes + $ice_pes + $glc_pes

#echo "task partitioning ... atm+ocn // lnd+ice+glc+rof"
echo ""
echo "ATM  gets $atm_pes"
echo "CPL  gets $cpl_pes"
echo "ICE  gets $ice_pes"
echo "LND  gets $lnd_pes"
echo "GLC  gets $glc_pes"
echo "DROF gets $rof_pes"
echo "OCN  gets $ocn_pes"
echo ""

#   total number of hw pes = 240
#   cpl hw pe range ~ from 0 to 59
#   cam hw pe range ~ from 0 to 119
#   pop2 hw pe range ~ from 120 to 239
#   clm hw pe range ~ from 0 to 59
#   cice hw pe range ~ from 60 to 89
#   sglc hw pe range ~ from 90 to 119
#   rtm hw pe range ~ from 120 to 149
#   TJH FIXME ... CLM could use a lot more processors.

./xmlchange NTHRDS_CPL=1,NTASKS_CPL=$cpl_pes
./xmlchange NTHRDS_GLC=1,NTASKS_GLC=$glc_pes,NINST_GLC=1
./xmlchange NTHRDS_ATM=1,NTASKS_ATM=$atm_pes,NINST_ATM=$num_instances
./xmlchange NTHRDS_LND=1,NTASKS_LND=$lnd_pes,NINST_LND=$num_instances
./xmlchange NTHRDS_ICE=1,NTASKS_ICE=$ice_pes,NINST_ICE=$num_instances
./xmlchange NTHRDS_ROF=1,NTASKS_ROF=$rof_pes,NINST_ROF=$num_instances
./xmlchange NTHRDS_OCN=1,NTASKS_OCN=$ocn_pes,NINST_OCN=$num_instances
./xmlchange ROOTPE_ATM=0
./xmlchange ROOTPE_OCN=0
./xmlchange ROOTPE_CPL=0
./xmlchange ROOTPE_LND=0
./xmlchange ROOTPE_ICE=0
./xmlchange ROOTPE_GLC=0
./xmlchange ROOTPE_ROF=0

# http://www.cesm.ucar.edu/models/cesm1.1/cesm/doc/usersguide/c1158.html#run_start_stop
# "A hybrid run indicates that CESM is initialized more like a startup, but uses
# initialization datasets from a previous case. This is somewhat analogous to a
# branch run with relaxed restart constraints. A hybrid run allows users to bring
# together combinations of initial/restart files from a previous case (specified
# by $RUN_REFCASE) at a given model output date (specified by $RUN_REFDATE).
# Unlike a branch run, the starting date of a hybrid run (specified by $RUN_STARTDATE)
# can be modified relative to the reference case. In a hybrid run, the model does not
# continue in a bit-for-bit fashion with respect to the reference case. The resulting
# climate, however, should be continuous provided that no model source code or
# namelists are changed in the hybrid run. In a hybrid initialization, the ocean
# model does not start until the second ocean coupling (normally the second day),
# and the coupler does a "cold start" without a restart file.

./xmlchange RUN_TYPE=startup
./xmlchange RUN_STARTDATE=$run_refdate
./xmlchange START_TOD=$run_reftod
./xmlchange RUN_REFDATE=$run_refdate
./xmlchange RUN_REFTOD=$run_reftod
./xmlchange GET_REFCASE=FALSE
./xmlchange EXEROOT=${exeroot}

./xmlchange CALENDAR=GREGORIAN

./xmlchange STOP_OPTION=$stop_option
./xmlchange STOP_N=$stop_n
./xmlchange CONTINUE_RUN=FALSE
./xmlchange RESUBMIT=$resubmit
./xmlchange PIO_TYPENAME=pnetcdf

./xmlchange CLM_CONFIG_OPTS='-bgc cn'

./xmlchange DOUT_S=FALSE
./xmlchange DOUT_S_ROOT=${archdir}
./xmlchange DOUT_S_SAVE_INT_REST_FILES=FALSE
./xmlchange DOUT_L_MS=FALSE
./xmlchange DOUT_L_MSROOT="csm/${case}"
./xmlchange DOUT_L_HTAR=FALSE

# level of debug output, 0=minimum, 1=normal, 2=more, 3=too much, valid values: 0,1,2,3 (integer)

./xmlchange DEBUG=FALSE
./xmlchange INFO_DBUG=0

# ==============================================================================
# Set up the case.
# This creates the EXEROOT and RUNDIR directories.
# ==============================================================================

./cesm_setup

if ( $status != 0 ) then
   echo "ERROR: Case could not be set up."
   exit -2
endif

# ==============================================================================
# Modify namelist templates for each instance.
# ==============================================================================

@ inst = 1
while ($inst <= $num_instances)

   set instance  = `printf %04d $inst`
   set instance2 = `printf %02d $inst`

   # ===========================================================================
   set fname = "user_nl_cam_${instance}"
   # ===========================================================================
   # For a HOP TEST ... empty_htapes = .false.
   # For a HOP TEST ... use a default fincl1 

   echo " inithist      = 'DAILY'"                      >> ${fname}
   echo " ncdata        = 'cam_initial_${instance}.nc'" >> ${fname}
   echo " empty_htapes  = .true. "                      >> ${fname}
   echo " fincl1        = 'PHIS:I' "                    >> ${fname}
   echo " nhtfrq        = -$assim_n "                   >> ${fname}
   echo " mfilt         = 1 "                           >> ${fname}

   # ===========================================================================
   set fname = "user_nl_pop2_${instance}"
   # ===========================================================================

   # POP Namelists
   # init_ts_suboption = 'data_assim'   for non bit-for-bit restarting (assimilation mode)
   # init_ts_suboption = 'null'         for 'perfect' restarting/forecasting
   # For a HOP TEST (untested)... tavg_file_freq_opt = 'nmonth' 'nday' 'once'"

   echo "init_ts_file = 'b40.20th.005_ens${instance2}.pop.r.2004-01-01-00000'" >> $fname
   echo "init_ts_suboption = 'spunup'" >> $fname

   # ===========================================================================
   set fname = "user_nl_cice_${instance}"
   # ===========================================================================
   # CICE Namelists
   
   echo "ice_ic = 'b40.20th.005_ens${instance2}.cice.r.2004-01-01-00000.nc'" >> $fname

   # ===========================================================================
   set fname = "user_nl_clm_${instance}"
   # ===========================================================================
   
   # Customize the land namelists
   # The history tapes are a work in progress. If you write out the instantaneous
   # flux variables every 30 minutes to the .h1. file, the forward observation
   # operators for these fluxes should just read them from the .h1. file rather
   # than trying to create them from the (incomplete DART) CLM state.
   # For a HOP TEST ... hist_empty_htapes = .false.
   # For a HOP TEST ... use a default hist_fincl1 
   #
   # old ... stagedir = /glade/scratch/afox/bptmp/MD_40_PME/run/MD_40_PME

   echo "finidat = '${CLM_stagedir}.clm2_${instance}.r.${run_refdate}-${run_reftod}.nc'" >> $fname
   echo "hist_empty_htapes = .true."                 >> $fname
   echo "hist_fincl1 = 'NEP'"                        >> $fname
   echo "hist_fincl2 = 'NEP','FSH','EFLX_LH_TOT_R'"  >> $fname
   echo "hist_nhtfrq = -$assim_n,1,"                 >> $fname
   echo "hist_mfilt  = 1,48"                         >> $fname
   echo "hist_avgflag_pertape = 'A','A'"             >> $fname

   # ===========================================================================
   set fname = "user_nl_rtm_${instance}"
   # ===========================================================================
   # RIVER RUNOFF CAN START FROM AN OLD CLM RESTART FILE

   echo "finidat_rtm = 'b40.20th.005_ens${instance2}.clm2.r.${run_refdate}-${run_reftod}.nc'" >> $fname

   @ inst ++
end

./preview_namelists

# ==============================================================================
# Update source files if need be
#    Ideally, using DART will not require any source mods.
#    Until then, this script accesses source mods from a hard-wired location below.
#    Those may eventually be packaged into the DART repository.
#    If you have additional source mods, they will need to be merged into any DART
#    mods and put in the SourceMods subdirectory found in the 'case' directory.
# ==============================================================================

if (    -d     ~/${cesmtag}/SourceMods ) then
   ${COPY} -r  ~/${cesmtag}/SourceMods/* ${caseroot}/SourceMods/
else
   echo "ERROR - No SourceMods for this case."
   echo "ERROR - No SourceMods for this case."
   echo "DART requires modifications to several src.pop2/ files."
   echo "These files can be downloaded from:"
   echo "http://www.image.ucar.edu/pub/DART/CESM/DART_SourceMods_cesm1_1_1.tar"
   echo "untar these into your HOME directory - they will create a"
   echo "~/cesm_1_1_1  directory with the appropriate SourceMods structure."
   exit -4
endif

# ==============================================================================
# build
# ==============================================================================

echo ''
echo 'Building the case'
echo ''

./${case}.build

if ( $status != 0 ) then
   echo "ERROR: Case could not be built."
   exit -5
endif

# ==============================================================================
# Stage the restarts now that the run directory exists
# I ran this compset once without setting user_nl_atm_00nn to see where the
# initial files come from.
# ==============================================================================

cd ${rundir}

echo ''
echo "Linking the restart files from the staging directories"
echo 'into the CESM run directory.'
echo ''

@ inst = 1
while ($inst <= $num_instances)
   set n4 = `printf %04d $inst`
   set n2 = `printf %02d $inst`

   echo ''
   echo "Staging restarts for instance $inst of $num_instances"

   ${LINK} ${CAM_stagedir}/cami-mam3_0000-01-01_0.9x1.25_L30_c100618.nc      cam_initial_${n4}.nc

   ${LINK} ${POP_stagedir}/b40.20th.005_ens${n2}.pop.r.2004-01-01-00000      .
   ${LINK} ${POP_stagedir}/b40.20th.005_ens${n2}.pop.r.2004-01-01-00000.hdr  .
   ${LINK} ${POP_stagedir}/b40.20th.005_ens${n2}.cice.r.2004-01-01-00000.nc  .
   ${LINK} ${RTM_stagedir}/b40.20th.005_ens${n2}.clm2.r.2004-01-01-00000.nc  .

#  ${LINK} ${POP_stagedir}/rpointer.ocn_${n4}.restart                       ${rundir}
#  ${LINK} ${POP_stagedir}/rpointer.ocn_${n4}.ovf                           ${rundir}
#  ${LINK} ${POP_stagedir}/rpointer.ice_${n4}                               ${rundir}
#  ${LINK} ${POP_stagedir}/b40.20th.005_ens${n2}.pop.ro.2004-01-01-00000    ${rundir}
#
   @ inst ++
end

# ==============================================================================
# Edit the run script to reflect project, queue, and wallclock
# ==============================================================================

cd ${caseroot}

echo ''
echo 'Updating the run script to set wallclock and queue.'
echo ''

if ( ~ -e  ${case}.run.original ) then
   ${COPY} ${case}.run ${case}.run.original
endif

source Tools/ccsm_getenv
set BATCH = `echo $BATCHSUBMIT | sed 's/ .*$//'`
switch ( $BATCH )
   case bsub*:
      # NCAR "bluefire", "yellowstone"
      set TIMEWALL=`grep BSUB ${case}.run | grep -e '-W' `
      set    QUEUE=`grep BSUB ${case}.run | grep -e '-q' `
      sed -e "s/ptile=[0-9][0-9]*/ptile=$ptile/" \
          -e "s/$TIMEWALL[3]/$timewall/" \
          -e "s/$QUEUE[3]/$queue/" < ${case}.run >! temp.$$
      ${MOVE} temp.$$  ${case}.run
   breaksw

   default:

   breaksw
endsw

# ==============================================================================
# The *.run script must be modified to call the DART assimilate script.
# The modifications are contained in a "here" document that MUST NOT
# expand the wildcards etc., before it is run. This is achieved by
# double-quoting the characters used to delineate the start/stop of
# the "here" document. No kidding. It has to be "EndOfText",
# not 'EndOfText' or EndOfText.
# ==============================================================================

echo ''
echo 'Adding the call to assimilate.csh to the *.run script.'
echo ''

cat << "EndOfText" >! add_to_run.txt

# -------------------------------------------------------------------------
# START OF DART: if CESM finishes correctly (pirated from ccsm_postrun.csh);
# perform an assimilation with DART.

set CplLogFile = `ls -1t cpl.log* | head -n 1`
if ($CplLogFile == "") then
   echo 'ERROR: Model did not complete - no cpl.log file present - exiting.'
   echo 'ERROR: Assimilation will not be attempted.'
   exit -1
endif

grep 'SUCCESSFUL TERMINATION' $CplLogFile
if ( $status == 0 ) then
#  ${CASEROOT}/assimilate.csh

   if ( $status == 0 ) then
      echo "`date` -- DART HAS FINISHED"
   else
      echo "`date` -- DART FILTER ERROR - ABANDON HOPE"
      setenv LSB_PJL_TASK_GEOMETRY "{(0)}"
      mpirun.lsf "exit -3"
   endif
else
   echo 'ERROR: Model did not complete successfully - exiting.'
   echo 'ERROR: Assimilation will not be attempted.'
   setenv LSB_PJL_TASK_GEOMETRY "{(0)}"
   mpirun.lsf "exit -2"
endif

# END OF DART BLOCK
# -------------------------------------------------------------------------
"EndOfText"

# Now that the "here" document is created,
# determine WHERE to insert it -- ONLY IF it is not already there.

grep "ABANDON HOPE" ${case}.run
set STATUSCHECK = $status

if ( ${STATUSCHECK} == 0 ) then
   echo "DART block already present in ${case}.run"
else if ( ${STATUSCHECK} == 1 ) then

   set MYSTRING = `grep --line-number "CSM EXECUTION HAS FINISHED" ${case}.run`
   set MYSTRING = `echo $MYSTRING | sed -e "s#:# #g"`

   @ origlen = `cat ${case}.run | wc -l`
   @ keep = $MYSTRING[1]
   @ lastlines = $origlen - $keep

   head -n $keep      ${case}.run    >! temp.$$
   cat                add_to_run.txt >> temp.$$
   tail -n $lastlines ${case}.run    >> temp.$$

   ${MOVE} temp.$$ ${case}.run
   ${REMOVE} add_to_run.txt

else
   echo "ERROR in grep of ${case}.run: aborting"
   echo "status was ${STATUSCHECK}"
   exit -6
endif

chmod 0744 ${case}.run

# ==============================================================================
# Stage the required parts of DART in the CASEROOT directory.
# ==============================================================================

# The standard CESM short-term archiving script may need to be altered
# to archive addtional or subsets of things, or to reduce the amount of
# data that is sent to the long-term archive.  Put a version of st_archive.sh
# in  ${DARTroot}/models/CESM/shell_scripts when/if necessary

if ( ~ -e  Tools/st_archive.sh.original ) then
   ${COPY} Tools/st_archive.sh Tools/st_archive.sh.original
else
   echo "a Tools/st_archive.sh backup copy already exists"
endif

# ${COPY} ${DARTroot}/models/CESM/shell_scripts/st_archive.sh   Tools/ TJH DEBUG
${COPY} ${DARTroot}/models/CESM/shell_scripts/assimilate.csh          assimilate.csh
${COPY} ${DARTroot}/models/CESM/shell_scripts/cam_assimilate.csh  cam_assimilate.csh
${COPY} ${DARTroot}/models/CESM/shell_scripts/pop_assimilate.csh  pop_assimilate.csh
${COPY} ${DARTroot}/models/CESM/shell_scripts/clm_assimilate.csh  clm_assimilate.csh

# ==============================================================================
# Stage the DART executables in the CESM execution root directory: EXEROOT
# ==============================================================================

${COPY} ${DARTroot}/models/cam/work/cam_to_dart   ${exeroot}/.
${COPY} ${DARTroot}/models/cam/work/dart_to_cam   ${exeroot}/.
${COPY} ${DARTroot}/models/cam/work/filter        ${exeroot}/filter_cam
${COPY} ${DARTroot}/models/cam/work/filter        ${exeroot}/filter
${COPY} ${DARTroot}/models/cam/work/input.nml                cam_input.nml

${COPY} ${DARTroot}/models/clm/work/clm_to_dart   ${exeroot}/.
${COPY} ${DARTroot}/models/clm/work/dart_to_clm   ${exeroot}/.
${COPY} ${DARTroot}/models/clm/work/filter        ${exeroot}/filter_clm
${COPY} ${DARTroot}/models/clm/work/input.nml                clm_input.nml

${COPY} ${DARTroot}/models/POP/work/pop_to_dart   ${exeroot}/.
${COPY} ${DARTroot}/models/POP/work/dart_to_pop   ${exeroot}/.
${COPY} ${DARTroot}/models/POP/work/filter        ${exeroot}/filter_pop
${COPY} ${DARTroot}/models/POP/work/input.nml                pop_input.nml

${COPY} ${DARTroot}/models/CESM/work/cesm_to_dart ${exeroot}/.
${COPY} ${DARTroot}/models/CESM/work/dart_to_cesm ${exeroot}/.
${COPY} ${DARTroot}/models/CESM/work/filter       ${exeroot}/filter_cesm
${COPY} ${DARTroot}/models/CESM/work/input.nml               input.nml

# ==============================================================================
# What to do next
# ==============================================================================

echo ''
echo "Time to check the case."
echo ''
echo "cd into ${caseroot}"
echo "Modify what you like in input.nml, make sure the observation directory"
echo "names set in assimilate.csh match those on your system, and submit"
echo "the CESM job by running:"
echo "./${case}.submit"
echo ''
echo "For continued submissions after the initial (hybrid) startup,"
echo "make the following changes to the env_run variables:"
echo ''
echo "  ./xmlchange -file env_run.xml -id STOP_N        -val 24"
echo "  ./xmlchange -file env_run.xml -id CONTINUE_RUN  -val TRUE"
echo "  ./xmlchange -file env_run.xml -id RESUBMIT      -val <your_favorite_number>"
echo ''
echo "Once you get to 2004-01-04, there are more things to do ... "
echo "in the CASEROOT directory, uncomment the ncdata in the user_nl_cam* files ..."
echo "ncdata       = 'cam_initial_${instance}.nc'" >> ${fname}
echo "in the run directory, link the current cam initial files ..."
echo "make sure the history tapes for cam,clm are being created at the right frequency"
echo ''
echo "Check the streams listed in the streams text files.  If more or different"
echo 'dates need to be added, then do this in the $CASEROOT/user_*files*'
echo "then invoke 'preview_namelists' so you can check the information in the"
echo "CaseDocs or ${rundir} directories."
echo ''

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

