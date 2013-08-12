#!/bin/csh
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$

# ---------------------
# Purpose
# ---------------------
#
# This script is designed to set up, stage, and build a multi-instance run of CESM
# using a B compset where CAM, POP, and CLM are all active. The initial states for
# the models come from as single reference case so a CESM hybrid setup is used.
#
# POP: uses the result of the 'b40.20th.005_ens${instance}' experiments. The POP
#      restart files were saved as binary files, which is somewhat problematic for
#      data assimilation. Consequently, the entire model state must be advanced for
#      several days before a viable netCDF restart file can be produces. We advance
#      for 72 hours.
#
# RTM: uses the result of the 'b40.20th.005_ens${instance}' experiments. These
#      restart files are actually the CLM restart files from this experiment because
#      the RTM was part of CLM when the experiment was run. The standalone version of
#      RTM can read old CLM restart files.
#
# CICE: uses the result of the 'b40.20th.005_ens${instance}' experiments.
#
# CLM: uses the result of the 'b40.20th.005_ens${instance}' experiments - sort of.
#      CLM has changed since then and the CLM restart files needed to be converted
#      to the new format. I ran 'interpinic' to (essentially) reformat the files
#      and changed the CASENAME to 'cesm_test' and used the multi-instance naming
#      convention for the new files.
#
# CAM: We want to use a newer version of CAM than was used for the b40.20th.005
#      experiments. Consequently, we took a SINGLE instance of CAM and replicated
#      it for all the desired instances. After just a few days, the differences in
#      the ocean and land states will induce enough variability in the CAM states.
#
# Much of the complexity comes from ensuring compatibility between the namelists
# for each instance and staging of the files. The original experiments were run
# before the multi-instance capacity was developed and the naming convention decided.
# Consequently, there is a lot of manipulation of the 'instance' portion of the
# filenames.
#
# This script results in a viable setup for a CESM multi-instance experiment. You
# are STRONGLY encouraged to run the multi-instance CESM a few times and experiment
# with different settings BEFORE you try to assimilate observations. The amount of
# data volume is quite large and you should become comfortable using CESM's restart
# capability to re-stage files in your RUN directory
#
# To perform assimilation, it is required only to insert a few dozen lines into the
# CASEROOT/*.run script. This, and the required setup, is performed in
# CESM_DART_config - which can be run at a later date. e.g. you can use this
# script to advance this ensemble from 2004-01-01 to 2004-02-01 and then run
# CESM_DART_config to augment the existing run script, modify STOP_N to 24 hours,
# and start assimilating observationg at midnight from 2004-02-01 on ...
#
# This script relies heavily on the information in:
# http://www.cesm.ucar.edu/models/cesm1.1/cesm/doc/usersguide/book1.html
#
# ---------------------
# How to set up the script
# ---------------------
#
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

setenv case                 cesm_continue
setenv compset              B_2000_CAM5
setenv resolution           0.9x1.25_gx1v6
setenv cesmtag              cesm1_1_1
setenv num_instances        30

# ==============================================================================
# define machines and directories
#
# mach            Computer name
# cesmroot        Location of the cesm code base
#                 For cesm1_1 on yellowstone
# caseroot        Your (future) cesm case directory, where this CESM+DART will be built.
#                    Preferably not a frequently scrubbed location.
#                    This script will delete any existing caseroot, so this script,
#                    and other useful things should be kept elsewhere.
# rundir          (Future) Run-time directory; scrubbable, large amount of space needed.
# exeroot         (Future) directory for executables - scrubbable, large amount of space needed.
# archdir         (Future) Short-term archive directory
#                    until the long-term archiver moves it to permanent storage.
# dartroot        Location of _your_ DART installation
#                    This is passed on to the CESM_DART_config script.
# ==============================================================================

setenv mach         yellowstone
setenv cesmroot     /glade/p/cesm/cseg/collections/$cesmtag
setenv caseroot     /glade/p/work/${USER}/cases/${case}
setenv exeroot      /glade/scratch/${USER}/${case}/bld
setenv rundir       /glade/scratch/${USER}/${case}/run
setenv archdir      /glade/scratch/${USER}/archive/${case}
setenv dartroot     /glade/u/home/${USER}/svn/DART/trunk

# ==============================================================================
# configure settings ... run_startdate format is yyyy-mm-dd
# ==============================================================================

setenv run_refcase cesm_hybrid
setenv refyear     2004
setenv refmon      01
setenv refday      10
setenv run_reftod  00000
setenv run_refdate $refyear-$refmon-$refday

# THIS IS THE LOCATION of the 'reference case'.

set stagedir = /glade/scratch/thoar/archive/${run_refcase}/rest/${run_refdate}-${run_reftod}

# ==============================================================================
# runtime settings --  How many assimilation steps will be done after this one
#                      plus archiving options
#
# resubmit      How many job steps to run on continue runs (will be 0 initially)
# stop_option   Units for determining the forecast length between assimilations
# stop_n        Number of time units in the first forecast
# assim_n       Number of time units between assimilations
#
# If the long-term archiver is off, you get a chance to examine the files before
# you
#
# ==============================================================================

setenv resubmit            12
setenv stop_option         nhours
setenv stop_n              6
setenv assim_n             6
setenv short_term_archiver on
setenv long_term_archiver  off

# ==============================================================================
# job settings
#
# queue      can be changed during a series by changing the ${case}.run
# timewall   can be changed during a series by changing the ${case}.run
#
# TJH: 30 instances with 900 pes took about 30 minutes on yellowstone.
# ==============================================================================

setenv ACCOUNT      P86850054
setenv queue        economy
setenv timewall     0:20

# ==============================================================================
# set these standard commands based on the machine you are running on.
# ==============================================================================

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

   breaksw
   default:
      # NERSC "hopper", NWSC "yellowstone"
      set   MOVE = '/bin/mv -fv'
      set   COPY = '/bin/cp -fv --preserve=timestamps'
      set   LINK = '/bin/ln -fvs'
      set REMOVE = '/bin/rm -fr'

   breaksw
endsw

# ==============================================================================
# Make sure the CESM directories exist.
# VAR is the shell variable name, DIR is the value
# ==============================================================================

foreach VAR ( stagedir cesmroot dartroot )
   set DIR = `eval echo \${$VAR}`
   if ( ! -d $DIR ) then
      echo "ERROR: directory '$DIR' not found"
      echo " In the setup script check the setting of: $VAR "
      exit -1
   endif
end

# ==============================================================================
# Create the case - this creates the CASEROOT directory.
#
# For list of the pre-defined cases: ./create_newcase -list
# To create a variant case, see the CESM documentation and carefully
# incorporate any needed changes into this script.
# ==============================================================================

# fatal idea to make caseroot the same dir as where this setup script is
# since the build process removes all files in the caseroot dir before
# populating it.  try to prevent shooting yourself in the foot.

if ( $caseroot == `dirname $0` ) then
   echo "ERROR: the setup script should not be located in the caseroot"
   echo "directory, because all files in the caseroot dir will be removed"
   echo "before creating the new case.  move the script to a safer place."
   exit -1
endif

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
# Record the DARTROOT directory and copy the DART setup script to CASEROOT.
# CESM_DART_config can be run at some later date if desired, but it presumes
# to be run from a CASEROOT directory. If CESM_DART_config does not exist locally,
# then it better exist in the expected part of the DARTROOT tree.
# ==============================================================================

if ( ! -e CESM_DART_config ) then
   ${COPY} ${dartroot}/models/CESM/shell_scripts/CESM_DART_config .
endif

if (   -e CESM_DART_config ) then
   sed -e "s#BOGUS_DART_ROOT_STRING#$dartroot#" < CESM_DART_config >! temp.$$
   ${MOVE} temp.$$ ${caseroot}/CESM_DART_config
   chmod 755       ${caseroot}/CESM_DART_config
else
   echo "WARNING: the script to configure for data assimilation is not available."
   echo "         CESM_DART_config should be present locally or in"
   echo "         ${dartroot}/models/CESM/shell_scripts/"
   echo "         You can stage this script later, but you must manually edit it"
   echo "         to reflect the location of the DART code tree."
endif

# ==============================================================================
# Configure the case.
# ==============================================================================

cd ${caseroot}

source ./Tools/ccsm_getenv || exit -2

@ ptile = $MAX_TASKS_PER_NODE / 2
@ nthreads = 1

# Save a copy for debug purposes
foreach FILE ( *xml )
   if ( ! -e        ${FILE}.original ) then
      ${COPY} $FILE ${FILE}.original
   endif
end

if ( $num_instances < 10) then

   # This is only for the purpose of debugging the code.
   @ atm_tasks = $ptile * $num_instances * 4
   @ ocn_tasks = $ptile * $num_instances * 4
   @ lnd_tasks = $ptile * $num_instances * 4
   @ ice_tasks = $ptile * $num_instances * 1
   @ glc_tasks = $ptile * $num_instances
   @ rof_tasks = $ptile * $num_instances
   @ cpl_tasks = $ptile * 4

else

   # This works, but a more efficient layout should be used
   @ atm_tasks = $ptile * $num_instances * 2
   @ ocn_tasks = $ptile * $num_instances * 2
   @ lnd_tasks = $ptile * $num_instances * 2
   @ ice_tasks = $ptile * $num_instances
   @ glc_tasks = $ptile * $num_instances
   @ rof_tasks = $ptile * $num_instances
   @ cpl_tasks = $ptile * $num_instances

endif

# echo "task partitioning ... perhaps ... atm // ocn // lnd+ice+glc+rof"
# presently, all components run 'serially' - one after another.
echo ""
echo "ATM  gets $atm_tasks"
echo "CPL  gets $cpl_tasks"
echo "ICE  gets $ice_tasks"
echo "LND  gets $lnd_tasks"
echo "GLC  gets $glc_tasks"
echo "DROF gets $rof_tasks"
echo "OCN  gets $ocn_tasks"
echo ""

./xmlchange NTHRDS_CPL=$nthreads,NTASKS_CPL=$cpl_tasks
./xmlchange NTHRDS_GLC=$nthreads,NTASKS_GLC=$glc_tasks,NINST_GLC=1
./xmlchange NTHRDS_ATM=$nthreads,NTASKS_ATM=$atm_tasks,NINST_ATM=$num_instances
./xmlchange NTHRDS_LND=$nthreads,NTASKS_LND=$lnd_tasks,NINST_LND=$num_instances
./xmlchange NTHRDS_ICE=$nthreads,NTASKS_ICE=$ice_tasks,NINST_ICE=$num_instances
./xmlchange NTHRDS_ROF=$nthreads,NTASKS_ROF=$rof_tasks,NINST_ROF=$num_instances
./xmlchange NTHRDS_OCN=$nthreads,NTASKS_OCN=$ocn_tasks,NINST_OCN=$num_instances
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
# and the coupler does a "cold start" without a restart file."

# TJH:
# DART's CAM implementation causes a bit more complexity. DART only uses CAM _initial_
# files, not RESTART files, so there are sourcemods to force a hybrid start for CAM to
# read initial files - even when CONTINUE_RUN = TRUE.

./xmlchange RUN_TYPE=hybrid
./xmlchange RUN_STARTDATE=$run_refdate
./xmlchange START_TOD=$run_reftod
./xmlchange RUN_REFCASE=$run_refcase
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

# This sets the ocean coupling time to 6 hours.
# All related namelist settings are based on this value.
./xmlchange OCN_NCPL=4

./xmlchange CLM_CONFIG_OPTS='-bgc cn'

if ($short_term_archiver == 'off') then
   ./xmlchange DOUT_S=FALSE
else
   ./xmlchange DOUT_S=TRUE
   ./xmlchange DOUT_S_ROOT=${archdir}
   ./xmlchange DOUT_S_SAVE_INT_REST_FILES=FALSE
endif
if ($long_term_archiver == 'off') then
   ./xmlchange DOUT_L_MS=FALSE
else
   ./xmlchange DOUT_L_MS=TRUE
   ./xmlchange DOUT_L_MSROOT="csm/${case}"
   ./xmlchange DOUT_L_HTAR=FALSE
endif

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
# Edit the run script to reflect queue and wallclock
# ==============================================================================

echo ''
echo 'Updating the run script to set wallclock and queue.'
echo ''

if ( ! -e  ${case}.run.original ) then
   ${COPY} ${case}.run ${case}.run.original
endif

source Tools/ccsm_getenv
set BATCH = `echo $BATCHSUBMIT | sed 's/ .*$//'`
switch ( $BATCH )
   case bsub*:
      # NCAR "bluefire", "yellowstone"
      set TIMEWALL=`grep BSUB ${case}.run | grep -e '-W' `
      set    QUEUE=`grep BSUB ${case}.run | grep -e '-q' `
      sed -e "s/$TIMEWALL[3]/$timewall/" \
          -e "s/ptile=[0-9][0-9]*/ptile=$ptile/" \
          -e "s/$QUEUE[3]/$queue/" < ${case}.run >> temp.$$
          ${MOVE} temp.$$ ${case}.run
          chmod 755       ${case}.run
   breaksw

   default:

   breaksw
endsw

# ==============================================================================
# Modify namelist templates for each instance. This is a bit of a nuisance in
# that we are pulling in restart and initial files from 'all over the place'
# and each model component has a different strategy.
#
# In a hybrid run with CONTINUE_RUN = FALSE (i.e. just starting up):
#
# CAM has been forced to read initial files - specified by namelist var:ncdata
# POP reads from pointer files
# CICE reads from namelist variable 'ice_ic'
# CLM builds its own 'finidat' value from the REFCASE variables but in CESM1_1_1
#     it does not use the instance string. There is a patch for clm.buildnml.csh
# RTM reads from namelist variable 'finidat_rtm', but rtm.buildnml.csh also is buggy.
#
# All of these must later on be staged with these same filenames.
# OR - all these namelists can be changed to match whatever has been staged.
# MAKE SURE THE STAGING SECTION OF THIS SCRIPT MATCHES THESE VALUES.
# ==============================================================================

@ inst = 1
while ($inst <= $num_instances)

   # instance now includes the leading underscore
   set instance  = `printf _%04d $inst`

   # ===========================================================================
   set fname = "user_nl_cam${instance}"
   # ===========================================================================
   # For a HOP TEST ... empty_htapes = .false.
   # For a HOP TEST ... use a default fincl1
   # FIXME ... add documentation for configuring CAM history files

   echo " inithist      = '6-HOURLY'"                   >> ${fname}
   echo " ncdata        = 'cam_initial${instance}.nc'"  >> ${fname}
   echo " empty_htapes  = .true. "                      >> ${fname}
   echo " fincl1        = 'PHIS:I' "                    >> ${fname}
   echo " nhtfrq        = -$assim_n "                   >> ${fname}
   echo " mfilt         = 1 "                           >> ${fname}

   # ===========================================================================
   set fname = "user_nl_clm${instance}"
   # ===========================================================================

   # Customize the land namelists
   # The filename is built using the REFCASE/REFDATE/REFTOD information.
   #
   # This is the time to consider how DART and CESM will interact.  If you intend
   # on assimilating flux tower observations (nominally at 30min intervals),
   # then it is required to create a .h1. file with the instantaneous flux
   # variables every 30 minutes. Despite being in a namelist, these values
   # HAVE NO EFFECT once CONTINUE_RUN = TRUE so now is the time to set these.
   #
   # DART's forward observation operators for these fluxes just reads them
   # from the .h1. file rather than trying to create them from the subset of
   # CLM variables that are available in the DART state vector.
   #
   # For a HOP TEST ... hist_empty_htapes = .false.
   # For a HOP TEST ... use a default hist_fincl1
   #
   # FIXME ... add documentation for configuring CLM history files

   @ thirtymin = $assim_n * 2

   echo "hist_empty_htapes = .true."                 >> $fname
   echo "hist_fincl1 = 'NEP'"                        >> $fname
   echo "hist_fincl2 = 'NEP','FSH','EFLX_LH_TOT_R'"  >> $fname
   echo "hist_nhtfrq = -$assim_n,1,"                 >> $fname
   echo "hist_mfilt  = 1,$thirtymin"                 >> $fname
   echo "hist_avgflag_pertape = 'A','A'"             >> $fname

   # ===========================================================================
   set fname = "user_nl_pop2${instance}"
   # ===========================================================================

   # POP Namelists
   # GIVEN: init_ts_option = 'ccsm_hybrid'  when  RUN_TYPE=hybrid, then
   #
   # init_ts_suboption = 'data_assim'   --> non bit-for-bit restarting (assimilation mode)
   # init_ts_suboption = 'rest'         --> default behavior
   #
   # No matter the setting of init_ts_suboption, POP always uses the
   # 'data_assim' behavior when RUN_TYPE=hybrid and CONTINUE_RUN=FALSE.
   # As soon as CONTINUE_RUN=TRUE, the value becomes important.
   #
   # For a HOP TEST ... Would like to have restart files every day, not just for end.
   # For a HOP TEST (untested)... tavg_file_freq_opt = 'nmonth' 'nday' 'once'"

   echo "init_ts_suboption = 'rest'" >> $fname

   # ===========================================================================
   set fname = "user_nl_cice${instance}"
   # ===========================================================================
   # CICE Namelists

   echo "ice_ic = '${run_refcase}.cice${instance}.r.${run_refdate}-${run_reftod}.nc'" >> $fname

   # ===========================================================================
   set fname = "user_nl_rtm${instance}"
   # ===========================================================================
   # RIVER RUNOFF CAN START FROM AN OLD CLM RESTART FILE
   # you can specify the RTM filename here and override the settings from
   # RUN_REFCASE/RUN_REFDATE/RUN_REFTOD (something you cannot do with CLM).

   echo "finidat_rtm = '${run_refcase}.rtm${instance}.r.${run_refdate}-${run_reftod}.nc'" >> $fname

   @ inst ++
end

# ==============================================================================
# Update source files.
#    Ideally, using DART would not require any modifications to the model source.
#    Until then, this script accesses sourcemods from a hardwired location.
#    If you have additional sourcemods, they will need to be merged into any DART
#    mods and put in the SourceMods subdirectory found in the 'case' directory.
# ==============================================================================

if (    -d     ~/${cesmtag}/SourceMods ) then
   ${COPY} -r  ~/${cesmtag}/SourceMods/* ${caseroot}/SourceMods/
else
   echo "ERROR - No SourceMods for this case."
   echo "ERROR - No SourceMods for this case."
   echo "DART requires modifications to several src files."
   echo "These files can be downloaded from:"
   echo "http://www.image.ucar.edu/pub/DART/CESM/DART_SourceMods_cesm1_1_1.tar"
   echo "untar these into your HOME directory - they will create a"
   echo "~/cesm_1_1_1  directory with the appropriate SourceMods structure."
   exit -4
endif

# The CESM multi-instance capability is relatively new and still has a few
# implementation bugs. These are known problems and will be fixed soon.
# this should be removed when the files are fixed:

echo "REPLACING BROKEN CESM FILES HERE - SHOULD BE REMOVED WHEN FIXED"
echo caseroot is ${caseroot}
if ( -d ~/${cesmtag} ) then

   # preserve the original version of the files
   if ( ! -e  ${caseroot}/Buildconf/clm.buildnml.csh.original ) then
      ${MOVE} ${caseroot}/Buildconf/clm.buildnml.csh \
              ${caseroot}/Buildconf/clm.buildnml.csh.original
   endif
   if ( ! -e  ${caseroot}/Buildconf/rtm.buildnml.csh.original ) then
      ${MOVE} ${caseroot}/Buildconf/rtm.buildnml.csh \
              ${caseroot}/Buildconf/rtm.buildnml.csh.original
   endif
   if ( ! -e  ${caseroot}/preview_namelists.original ) then
      ${MOVE} ${caseroot}/preview_namelists \
              ${caseroot}/preview_namelists.original
   endif

   # patch/replace the broken files
   ${COPY} ~/${cesmtag}/clm.buildnml.csh  ${caseroot}/Buildconf/.
   ${COPY} ~/${cesmtag}/rtm.buildnml.csh  ${caseroot}/Buildconf/.
   ${COPY} ~/${cesmtag}/preview_namelists ${caseroot}/.

endif

./preview_namelists

# ==============================================================================
# Stage the restarts now that the run directory exists
# THIS IS THE STAGING SECTION - MAKE SURE THIS MATCHES THE NAMELISTS.
# POP/CAM/CICE read from pointer files. The others use namelist values initially.
# ==============================================================================

cat << EndOfText >! stage_initial_cesm_files
#!/bin/sh

cd ${rundir}

echo ''
echo 'Copying the restart files from the staging directory.'
echo ''

${COPY} ${stagedir}/* .

let inst=1
while ((\$inst <= $num_instances)); do
   inst_string=\`printf _%04d \$inst\`

   ${LINK} ${stagedir}/${run_refcase}.cam\${inst_string}.i.${run_refdate}-${run_reftod}.nc cam_initial\${inst_string}.nc

   let inst+=1
done

exit 0

EndOfText
chmod 0755 stage_initial_cesm_files

./stage_initial_cesm_files

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
# What to do next
# ==============================================================================

echo ""
echo "Time to check the case."
echo ""
echo "1) cd ${rundir}"
echo "   and check the compatibility between the namelists/pointer"
echo "   files and the files that were staged."
echo ""
echo "2) cd ${caseroot}"
echo "   (on yellowstone) If the ${case}.run script still contains:"
echo '   #BSUB -R "select[scratch_ok > 0]"'
echo "   around line 9, delete it."
echo ""
echo "3) If you want to assimilate 'right away', configure and execute"
echo "   the ${caseroot}/CESM_DART_config script."
echo ""
echo "3) Verify the contents of env_run.xml and submit the CESM job:"
echo "   ./${case}.submit"
echo ""
echo "4) After the job has run, make sure it worked and that "
echo "   POP is creating netCDF restart files."
echo ""
echo "5) To extend the run in $assim_n '"$stop_option"' steps,"
echo "   change the env_run.xml variables:"
echo ""
echo "  ./xmlchange CONTINUE_RUN=TRUE"
echo "  ./xmlchange RESUBMIT=<number_of_cycles_to_run>"
echo "  ./xmlchange STOP_N=$assim_n"
echo ""

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

