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
# that has CAM, CLM, and CICE as active components over a single data ocean,
# and will use DART to assimilate observations at regular intervals.
# This script does not build DART.
#
# This script relies heavily on the information in:
# http://www.cesm.ucar.edu/models/cesm1.1/cesm/doc/usersguide/book1.html
#
# ---------------------
# How to set up the script
# ---------------------
# -- Either edit and run this script in the $DART/models/cam/shell_scripts
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
#       CAM/CLM/CICE initial ensemble
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
# ==============================================================================
# === IMPORTANT modifications to the distribution code before ANYTHING
# ==============================================================================

# Had to edit the following to remove the LSB_PJL_... word too long error
# cesm1_1_beta04/scripts/ccsm_utils/Machines/mkbatch.bluefire .
# On bluefire for cesm1_1_beta04 it's not enough to use the hard-wired
# SourceMods below; mkbatch.bluefire cannot be modified via SourceMods.
# So the cesm1_1_beta04 version in Tim's directory must be used.
# As long as OMP_NUM_THREADS == 1 ... the default is fine.
# This may also be a problem if the number of nodes is increased, e.g.
# if the regular memory nodes are used, or if the resolution is increased,
# because the number of entries in BIND_THRD_GEOMETRY will increase.
#
# The cesm1_1_beta04 lt_archive script did not create parent dirs
# if they did not already exist.  To fix the script, edit:
#   ${cesmroot}/scripts/ccsm_utils/Tools/lt_archive.csh
# and change 'mkdir' to 'mkdir -p'.  This is fixed in more recent
# versions of the code.
#
# For these reasons, it is required that you use ~thoar/cesm1_1_beta04
# as the source of the CESM code.

# ==============================================================================
# ====  Set case options
# ==============================================================================

# case will be used many ways;
#    directory and file names, both locally and on HPSS, and
#    script names; so consider it's length and information content.
# num_instances:  Number of ensemble members

setenv case                 cam_test_1
setenv compset              F_AMIP_CAM5
setenv cesmtag              cesm1_1_1
setenv resolution           f09_f09
setenv num_instances        4

# ==============================================================================
# define machines and directories
#
# mach            Computer name
# cesm_datadir    Root path of the public CESM data files
# cesmroot        Location of the cesm code base
#                 For cesm1_1 on yellowstone
# DARTroot        Location of DART code tree.
#                    Executables, scripts and input in $DARTroot/models/cam/...
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

# ==============================================================================
# configure settings
# ==============================================================================
# yyyy-mm-dd
setenv run_startdate 2008-10-31

setenv sst_dataset ${cesm_datadir}/atm/cam/sst/sst_HadOIBl_bc_0.9x1.25_1850_2008_c100127.nc
setenv year_start  1850
setenv year_end    2010

# TJH FIXME ... the sst_dataset above does not appear in that directory._CAM5
# Try this:  ${cesm_datadir}/atm/cam/sst/sst_HadOIBl_bc_0.9x1.25_1850_2011_c110307.nc

# ==============================================================================
# runtime settings --  How many assimilation steps will be done after this one
#
# stop_option   Units for determining the forecast length between assimilations
#               Changing stop_option requires changes to user_nl_cam below.
# stop_n        Number of time units in the forecast
# ==============================================================================

setenv resubmit      0
setenv stop_option   nhours
setenv stop_n        6

# ==============================================================================
# job settings
#
# timewall   can be changed during a series by changing the ${case}.run
# queue      can be changed during a series by changing the ${case}.run
#         lrg_ queues are used in order to fit more instances on each node.
#         FV 1-degree can comfortably fit 4 instances on 1 lrg_ node (~60 gbyte)
#         On bluefire the regular queue (or higher) is probably necessary,
#         because it appears that there are not separate queues for the lrg memory
#         and the regular memory nodes.  So economy jobs requesting smaller numbers
#         of processors seem to prevent this lrg_economy job (20 nodes for 1-degree)
#         from running for long periods.
# ==============================================================================

setenv ACCOUNT      P86850054
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

   breaksw
endsw

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
   ${COPY} $FILE ${FILE}.original
end

# num_tasks_per_instance = #tasks_node / #instances_node
# Bluefire: #tasks_node = 64 using SMT
#           #instances_node = 1-degree: 4 on lrg_ nodes, 2 on standard.
#                             2-degree: 16 on lrg_ nodes, 8 on standard
#           CAM5; to be determined, but no more than listed for CAM4
set num_tasks_per_node = 15
set num_tasks_per_instance = 8
set num_threads = 1

# This is hard-wiring for the current (1/17/2011) multi-instance CESM restriction
# that all instances must be advanced simultaneously.  Work is underway to relax that.
# @ total_nt = $num_instances * $num_tasks_per_instance
# echo "total MPI tasks requested = $total_nt"

# Atm gets all the nodes and runs.
# The other components divide them up.
# ? ? ?  What about sglc?  It's a stub, and doesn't matter what pes are assigned to it.
# This algorithm figures out whether there are enough processors requested
# to run each component on whole nodes, or the components need to share some nodes.
# The hard-coded numbers (and ratios between them) are estimates; change them if you
# know better.
#@ atm_pes  = $total_nt
#@ large_small = $total_nt / (8 * $num_tasks_per_node)
#if ($large_small > 0) then
#   # Large_small > 0 means there are at least 8 nodes requested.
#   # Allot whole nodes to the major components.
#   @ cpl_pes = ( $total_nt / (3 * $num_tasks_per_node) ) * $num_tasks_per_node
#   @ ice_pes = ( $total_nt / (3 * $num_tasks_per_node) ) * $num_tasks_per_node
#   @ glc_pes = ( $total_nt / (3 * $num_tasks_per_node) ) * $num_tasks_per_node
#   @ ocn_pes = $num_tasks_per_node
#   @ rof_pes = $num_tasks_per_node
#   @ lnd_pes = $total_nt - ($cpl_pes + $ice_pes + $ocn_pes + $glc_pes + $rof_pes)
#else
   # dummy layout
   @ cpl_pes = 8
   @ glc_pes = 8
   @ ocn_pes = 8
   @ ice_pes = 32
   @ rof_pes = 32
   @ lnd_pes = 32
   @ atm_pes = 32
#endif

echo "task layout"
echo ""
echo "CPL gets $cpl_pes"
echo "ATM gets $atm_pes"
echo "ICE gets $ice_pes"
echo "LND gets $lnd_pes"
echo "ROF gets $rof_pes"
echo "GLC gets $glc_pes"
echo "OCN gets $ocn_pes"
echo ""

./xmlchange NTHRDS_CPL=$num_threads,NTASKS_CPL=$cpl_pes
./xmlchange NTHRDS_ATM=$num_threads,NTASKS_ATM=$atm_pes,NINST_ATM=$num_instances
./xmlchange NTHRDS_ICE=$num_threads,NTASKS_ICE=$ice_pes,NINST_ICE=$num_instances
./xmlchange NTHRDS_LND=$num_threads,NTASKS_LND=$lnd_pes,NINST_LND=$num_instances
./xmlchange NTHRDS_ROF=$num_threads,NTASKS_ROF=$rof_pes,NINST_ROF=$num_instances
./xmlchange NTHRDS_GLC=$num_threads,NTASKS_GLC=$glc_pes,NINST_GLC=1
./xmlchange NTHRDS_OCN=$num_threads,NTASKS_OCN=$ocn_pes,NINST_OCN=1
./xmlchange ROOTPE_OCN=$ice_pes

./xmlchange RUN_TYPE=startup
./xmlchange RUN_STARTDATE=$run_startdate
./xmlchange EXEROOT=${exeroot}

./xmlchange DOUT_S_ROOT=${archdir}
./xmlchange DOUT_S=FALSE
./xmlchange DOUT_S_SAVE_INT_REST_FILES=TRUE
./xmlchange DOUT_L_MS=FALSE
./xmlchange DOUT_L_MSROOT="csm/${case}"
./xmlchange DOUT_L_HTAR=FALSE

./xmlchange SSTICE_DATA_FILENAME=$sst_dataset
./xmlchange SSTICE_YEAR_ALIGN=$year_start
./xmlchange SSTICE_YEAR_START=$year_start
./xmlchange SSTICE_YEAR_END=$year_end

# Do not change the CALENDAR or the CONTINUE_RUN
# DOUT_S     is to turn on/off the short-term archiving
# DOUT_L_MS  is to store to the HPSS (formerly "MSS")

./xmlchange CALENDAR=GREGORIAN
./xmlchange STOP_OPTION=$stop_option
./xmlchange STOP_N=$stop_n
./xmlchange CONTINUE_RUN=FALSE
./xmlchange RESUBMIT=$resubmit
./xmlchange PIO_TYPENAME=pnetcdf

# The river transport model ON is useful only when using an active ocean or
# land surface diagnostics. Setting ROF_GRID to 'null' turns off the RTM.

./xmlchange ROF_GRID='r05'

# level of debug output, 0=minimum, 1=normal, 2=more, 3=too much, valid values: 0,1,2,3 (integer)

./xmlchange DEBUG=TRUE
./xmlchange INFO_DBUG=2

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
#
# hist_empty_htapes = .true.     suppresses the creation of all history files
# hist_fincl1 = 'TG',            except the first one, which will have one variable
# hist_nhtfrq = -$stop_n,        create one every $stop_n HOURS
# hist_mfilt  =  1,              with precisely one day in it
# hist_avgflag_pertape = 'I'     use instantaneous values - no average
#
# ==============================================================================

set chem_datapath = "${cesm_datadir}/atm/cam/chem/trop_mozart_aero"

@ inst = 1
while ($inst <= $num_instances)

   set instance  = `printf %04d $inst`

   # ===========================================================================
   set fname = "user_nl_cam_$instance"
   # ===========================================================================
   # A lot of the files specified here are because the 'default' files only
   # contain data through 2005 and we are interested in timeframes after that.

#  echo " iradae               = -$stop_n "                                                      >> ${fname}

   echo " inithist = 'ENDOFRUN' "                                                                >> ${fname}
   echo " ncdata = 'cam_initial_${instance}.nc' "                                                >> ${fname}
   echo " avgflag_pertape      = 'A','A' "                                                       >> ${fname}
   echo " div24del2flag        = 4 "                                                             >> ${fname}
   echo " empty_htapes         = .true. "                                                        >> ${fname}
   echo " fincl1               = 'PHIS:I' "                                                      >> ${fname}
   echo " fincl2               = 'PS', 'Z3', 'T', 'U', 'V', 'FLNT', 'PSL', 'OMEGA' "             >> ${fname}
   echo " nhtfrq               = -$stop_n,-1 "                                                   >> ${fname}
   echo " mfilt                = 1,1 "                                                           >> ${fname}

   echo " aerodep_flx_datapath  = '${chem_datapath}/aero' "                                      >> ${fname}
   echo " aerodep_flx_file      = 'aerosoldep_monthly_1849-2006_1.9x2.5_c090803.nc' "            >> ${fname}
   echo " aerodep_flx_cycle_yr  = 2000 "                                                         >> ${fname}
   echo " aerodep_flx_type      = 'CYCLICAL' "                                                   >> ${fname}
   echo " bndtvghg              = '${cesm_datadir}/atm/cam/ggas/ghg_hist_1765-2009_c100902.nc' " >> ${fname}
   echo " prescribed_ozone_file = 'ozone_1.9x2.5_L26_1850-2015_rcp45_c101108.nc' "               >> ${fname}
   echo " tracer_cnst_file      = 'oxid_1.9x2.5_L26_1850-2015_rcp45_c101108.nc' "                >> ${fname}

   echo " ext_frc_specifier     = "                                                              >> ${fname}
   echo " 'SO2    -> ${chem_datapath}/emis/ar5_mam3_so2_elev_1850-2010_c20100902_v12.nc',"       >> ${fname}
   echo " 'bc_a1  -> ${chem_datapath}/emis/ar5_mam3_bc_elev_1850-2010_c20100902_v12.nc',"        >> ${fname}
   echo " 'num_a1 -> ${chem_datapath}/emis/ar5_mam3_num_a1_elev_1850-2010_c20100902_v12.nc',"    >> ${fname}
   echo " 'num_a2 -> ${chem_datapath}/emis/ar5_mam3_num_a2_elev_1850-2010_c20100902_v12.nc',"    >> ${fname}
   echo " 'pom_a1 -> ${chem_datapath}/emis/ar5_mam3_oc_elev_1850-2010_c20100902_v12.nc',"        >> ${fname}
   echo " 'so4_a1 -> ${chem_datapath}/emis/ar5_mam3_so4_a1_elev_1850-2010_c20100902_v12.nc',"    >> ${fname}
   echo " 'so4_a2 -> ${chem_datapath}/emis/ar5_mam3_so4_a2_elev_1850-2010_c20100902_v12.nc' "    >> ${fname}

   echo " srf_emis_specifier    = "                                                              >> ${fname}
   echo " 'DMS    -> ${chem_datapath}/emis/aerocom_mam3_dms_surf_1849-2010_c20100902.nc',"       >> ${fname}
   echo " 'SO2    -> ${chem_datapath}/emis/ar5_mam3_so2_surf_1850-2010_c20100902_v12.nc',"       >> ${fname}
   echo " 'SOAG   -> ${chem_datapath}/emis/ar5_mam3_soag_1.5_surf_1850-2010_c20100902_v12.nc',"  >> ${fname}
   echo " 'bc_a1  -> ${chem_datapath}/emis/ar5_mam3_bc_surf_1850-2010_c20100902_v12.nc',"        >> ${fname}
   echo " 'num_a1 -> ${chem_datapath}/emis/ar5_mam3_num_a1_surf_1850-2010_c20100902_v12.nc',"    >> ${fname}
   echo " 'num_a2 -> ${chem_datapath}/emis/ar5_mam3_num_a2_surf_1850-2010_c20100902_v12.nc',"    >> ${fname}
   echo " 'pom_a1 -> ${chem_datapath}/emis/ar5_mam3_oc_surf_1850-2010_c20100902_v12.nc',"        >> ${fname}
   echo " 'so4_a1 -> ${chem_datapath}/emis/ar5_mam3_so4_a1_surf_1850-2010_c20100902_v12.nc',"    >> ${fname}
   echo " 'so4_a2 -> ${chem_datapath}/emis/ar5_mam3_so4_a2_surf_1850-2010_c20100902_v12.nc' "    >> ${fname}

   echo " solar_data_file = '${cesm_datadir}/atm/cam/solar/spectral_irradiance_Lean_1610-2009_ann_c100405.nc' " >> ${fname}


   # ===========================================================================
   set fname = "user_nl_clm_$instance"
   # ===========================================================================
   # faerdep does not make its way to the final namelist in the F_AMIP_CAM5 case:

#  echo "faerdep  = '${cesm_datadir}/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_rcp4.5_monthly_1849-2104_0.9x1.25_c100407.nc' " >> ${fname}
   echo "hist_empty_htapes = .true. "                       >> ${fname}
   echo "finidat           = 'clm_restart_${instance}.nc' " >> ${fname}


   # ===========================================================================
   set fname = "user_nl_cice_$instance"
   # ===========================================================================

   echo  " ice_ic = 'ice_restart_${instance}.nc' " >> ${fname}
  
   @ inst ++
end

# ==============================================================================
# to create custom streamfiles ...
# "To modify the contents of a stream txt file, first use preview_namelists to
#  obtain the contents of the stream txt files in CaseDocs, and then place a copy
#  of the modified stream txt file in $CASEROOT with the string user_ prepended."
#
# -or-
#
# we copy a template stream txt file from the
# $DARTroot/models/POP/shell_scripts directory and modify one for each instance.
#
# ==============================================================================

# none needed for CAM

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
   exit -3
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
   exit -4
endif

# ==============================================================================
# Stage the restarts now that the run directory exists
# ==============================================================================

# CAM5 ... /glade/p/work/raeder/Exp/Pincus_Cam5/archive/Exp3/rest/2008-11-01-00000

set stagedir = /glade/p/work/raeder/Exp/Pincus_Cam5/archive/Exp3/rest/2008-11-01-00000

echo ''
echo "Copying the restart files from ${stagedir}"
echo 'into the CESM run directory.'
echo ''

@ n = 1
while ($n <= $num_instances)
   set instance  = `printf %04d $n`
   echo "Staging restarts for instance $n of $num_instances"
   ${COPY} ${stagedir}/cam_initial_${n}.nc ${rundir}/cam_initial_${instance}.nc
   ${COPY} ${stagedir}/clm_restart_${n}.nc ${rundir}/clm_restart_${instance}.nc
   ${COPY} ${stagedir}/ice_restart_${n}.nc ${rundir}/ice_restart_${instance}.nc
   @ n++
end

if (  -e   ${stagedir}/prior_inflate_restart* ) then
   ${COPY} ${stagedir}/prior_inflate_restart* ${rundir}/.
endif

# ==============================================================================
# IMPORTANT: All resubmits must be type 'startup'.
# Change Tools/ccsm_postrun.csh line 83 to CONTINUE_RUN -val FALSE'
# ==============================================================================

echo ''
echo 'Changing Tools/ccsm_postrun.csh such that all the resubmits are "startup",'
echo 'which means CONTINUE_RUN should be FALSE in ccsm_postrun.csh'
echo ''

ex Tools/ccsm_postrun.csh <<ex_end
/use COMP_RUN_BARRIERS as surrogate for timing run logical/
/CONTINUE_RUN/
s;TRUE;FALSE;
wq
ex_end

# ==============================================================================
# Edit the run script to reflect project, queue, and wallclock
# ==============================================================================

echo ''
echo 'Updating the run script to set wallclock and queue.'
echo ''

${COPY} ${case}.run ${case}.run.orig

source Tools/ccsm_getenv
set BATCH = `echo $BATCHSUBMIT | sed 's/ .*$//'`
switch ( $BATCH )
   case bsub*:
      # NCAR "bluefire", "yellowstone"
      set TIMEWALL=`grep BSUB ${case}.run | grep -e '-W' `
      set    QUEUE=`grep BSUB ${case}.run | grep -e '-q' `
      sed -e "s/ptile=32/ptile=$ptile/" \
          -e "s/$TIMEWALL[3]/$timewall/" \
          -e "s/$QUEUE[3]/$queue/" < ${case}.run >! temp.$$
      ${MOVE} temp.$$  ${case}.run
   breaksw

   default:

   breaksw
endsw

chmod 0755 ${case}.run

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
   ${CASEROOT}/assimilate.csh

   if ( $status == 0 ) then
      echo "`date` -- DART HAS FINISHED"
   else
      echo "`date` -- DART FILTER ERROR - ABANDON HOPE"
      exit -3
   endif
else
   echo 'ERROR: Model did not complete successfully - exiting.'
   echo 'ERROR: Assimilation will not be attempted.'
   exit -2
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

   ${MOVE}   temp.$$ ${case}.run
   ${REMOVE} add_to_run.txt

else
   echo "ERROR in grep of ${case}.run: aborting"
   echo "status was ${STATUSCHECK}"
   exit -5
endif

# ==============================================================================
# Stage the required parts of DART in the CASEROOT directory.
# ==============================================================================

# The cesm1_1_beta04 release had an error in that it did not
# provide the lt_archive.sh script, and the one in the repos
# did not have the -p flag, which is a good idea. So, for now ...

if ( ~ -e  Tools/st_archive.sh.orig ) then
   ${COPY} Tools/st_archive.sh      Tools/st_archive.sh.orig
else
   echo "a Tools/st_archive.sh backup copy already exists"
endif

# ${COPY} ${DARTroot}/models/cam/shell_scripts/st_archive.sh  Tools/

${COPY} ${DARTroot}/models/cam/shell_scripts/assimilate.csh  .
${COPY} ${DARTroot}/models/cam/work/input.nml                .

# ==============================================================================
# Stage the DART executables in the CESM execution root directory: EXEROOT
# ==============================================================================

foreach FILE ( filter cam_to_dart dart_to_cam )
   ${COPY} ${DARTroot}/models/cam/work/${FILE} ${exeroot}/
   if ( $status != 0 ) then
      echo "ERROR: ${DARTroot}/models/cam/work/${FILE} not copied to ${exeroot}"
      echo "ERROR: ${DARTroot}/models/cam/work/${FILE} not copied to ${exeroot}"
      exit -6
   endif
end

# This script will copy an existing prior_inflate_restart to the run dir if found.
if (-f ${stagedir}/DART/prior_inflate_restart) then
   ${COPY} ${stagedir}/DART/prior_inflate_restart ${exeroot}/run/prior_inflate_restart.${run_startdate}-00000
   echo ''
   echo "${stagedir}/DART/prior_inflate_restart has been copied to "
   echo "${exeroot}/run/prior_inflate_restart.${run_startdate}-00000"
   echo 'If that has the wrong state vector, you will need to replace it before running.'
   echo ''
else
   echo ''
   echo 'If using inflation in DART you may need to copy an inflation restart file'
   echo "to ${exeroot}/run/prior_inflate_restart.${run_startdate}-00000"
   echo 'before running.  It must include the exact fields as your DART state vector.'
   echo "You can make one with ${DARTroot}/models/cam/work/fill_inflation_restart"
   echo ''
endif


# only warn people if a precomputed final_full for this number of instances
# does not already exist.
if (! -f ${DARTroot}/system_simulation/final_full_precomputed_tables/final_full.${num_instances}) then
   echo ''
   echo 'If you are using the DART sampling error correction feature'
   echo 'the assimilate.csh script will expect to copy this file:'
   echo "${DARTroot}/system_simulation/final_full_precomputed_tables/final_full.${num_instances}"
   echo 'and it does not exist for your number of ensemble members.'
   echo "Generate one by building and running ${DARTroot}/system_simulation/work/full_error"
   echo 'with the namelist set to your ensemble size.'
   echo ''
endif

# ==============================================================================
# What to do next
# ==============================================================================

echo ''
echo "Time to check the case."
echo ''
echo "cd into ${caseroot}"
echo 'Modify what you like in input.nml, make sure the observation directory'
echo 'names set in assimilate.csh match those on your system, and submit'
echo 'the CESM job by running:'
echo "./$case.submit"
echo ''

