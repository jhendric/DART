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
# -- Examin the whole script to identify things to change for your experiments.
# -- Provide any initial files needed by your run:
#       inflation
#       sampling error correction
#       CAM/CLM/CICE initial ensemble
#       ...
# -- Run this script.
# -- Edit the DART input.nml that appears in the $caseroot directory.
# -- Submit the job using $caseroot/${case}.${mach}.submit
#      ($mach may not be needed for cesm releases after cesm1_1_beta04)
#
# ---------------------
# Important features
# ---------------------
#
# If you want to change something in your case other than the runtime
# settings, you need to delete everything and start the run from scratch.
#
# ./${CASENAME}.*.clean_build
# ./configure -cleanall
#
# ====================================================================
# === IMPORTANT modifications to the distribution code before ANYTHING
# ====================================================================

# Had to edit the following to remove the LSB_PJL_... word too long error
# cesm1_1_beta04/scripts/ccsm_utils/Machines/mkbatch.bluefire .
# On bluefire for cesm1_1_bet04 it's not enough to use the hard-wired
# SourceMods below; mkbatch.bluefire cannot be modified via SourceMods.
# So the cesm1_1_beta04 version in Tim's directory must be used.
# As long as OMP_NUM_THREADS == 1 ... the default is fine.
# This may also be a problem if the number of nodes is increased, e.g.
# if the regular memory nodes are used, or if the resolution is increased,
# because the number of entries in BIND_THRD_GEOMETRY will increase.
#
# The cesm1_1_beta04 lt_archive script also did not create parent dirs
# sed -e "s#mkdir#mkdir -p#" < ${cesmroot}/scripts/ccsm_utils/Tools/lt_archive.csh \
#        >! Tools/lt_archive.csh

# For these reasons, it is required that you use ~thoar/cesm1_1_beta04
# as the source of the CESM code.

# ====================================================================
# ====  Set case options
# ====================================================================

# case will be used many ways;
#    directory and file names, both locally and on HPSS, and
#    script names; so consider it's length and information content.
# num_instances:  Number of ensemble members
# reuse_existing_case:
#    false; Remove $caseroot and $exeroot and rebuild
#    true;  configure -cleannamelist

setenv case                 Exp1
setenv compset              F_2000
setenv cesmtag              cesm1_1_beta04
setenv resolution           f09_f09
setenv num_instances        4
setenv reuse_existing_case  false

# ====================================================================
# define machines and directories
#
# mach            Computer name
# cesm_datadir    Root path of the public CESM data files
# cesmroot        Location of the cesm code base
#                 For cesm1_1_beta04 on bluefire, MUST use 'thoar' value provided.
# DARTroot        Location of DART code tree.
#                    Executables, scripts and input in $DARTroot/models/cam/...
# caseroot        Your (future) cesm case directory, where this CESM+DART will be built.
#                    Preferably not a frequently scrubbed location.
#                 caseroot will be deleted if reuse_existing_case is false (below)
#                 So this script, and other useful things should be kept elsewhere.
# exeroot         (Future) Run-time directory; scrubbable, large amount of space needed.
# archdir         (Future) Short-term archive directory
#                    until the long-term archiver moves it to permanent storage.
# ====================================================================

setenv mach         bluefire
setenv cesm_datadir /glade/proj3/cseg/inputdata
setenv cesmroot     /glade/home/thoar/${cesmtag}

setenv DARTroot     /glade/home/${USER}/svn/DART/dev

setenv caseroot     /glade/user/${USER}/cases/${case}
setenv exeroot      /ptmp/${USER}/${case}
setenv archdir      /ptmp/${USER}/archive/${case}
set    CESM_setup_dir = `pwd`

# ====================================================================
# configure settings
# ====================================================================

setenv run_startdate 2008-10-31

setenv sst_dataset ${cesm_datadir}/atm/cam/sst/sst_HadOIBl_bc_0.9x1.25_1850_2011_c110307.nc
setenv year_start  1850
setenv year_end    2010

# ====================================================================
# runtime settings --  How many assimilation steps will be done after this one
#
# stop_option   Units for determining the forecast length between assimilations
#               Changing stop_option requires changes to user_nl_cam below.
# stop_n        Number of time units in the forecast
# ====================================================================

setenv resubmit      0
setenv stop_option   nhours
setenv stop_n        6

# ====================================================================
# job settings
#
# timewall    can be changed during a series by changing the ${case}.${mach}.run
# queue   can be changed during a series by changing the ${case}.${mach}.run
#         lrg_ queues are used in order to fit more instances on each node.
#         FV 1-degree can comfortably fit 4 instances on 1 lrg_ node (~60 gbyte)
#         On bluefire the regular queue (or higher) is probably necessary,
#         because it appears that there are not separate queues for the lrg memory
#         and the regular memory nodes.  So economy jobs requesting smaller numbers
#         of processors seem to prevent this lrg_economy job (20 nodes for 1-degree)
#         from running for long periods.
# ====================================================================

setenv proj         12345678
setenv timewall     2:00
setenv queue        lrg_regular

# ====================================================================
# Create the case.
#
# For list of the pre-defined cases: ./create_newcase -list
# To create a variant case, see the CESM documentation and carefully
# incorporate any needed changes into this script.
# ====================================================================

if ("${reuse_existing_case}" == "false") then
   echo "removing old files from ${caseroot} and ${exeroot}"
   \rm -fr ${caseroot}
   \rm -fr ${exeroot}
   ${cesmroot}/scripts/create_newcase -case ${caseroot} -mach ${mach} \
                   -res ${resolution} -compset ${compset} -skip_rundb

   if ( $status != 0 ) then
      echo "ERROR: Case could not be created."
      exit 1
   endif
else
   cd ${caseroot}
   ./configure  -cleannamelist
endif

# ====================================================================
# Configure the case.
# ====================================================================

cd ${caseroot}

./xmlchange -file env_build.xml    -id EXEROOT        -val ${exeroot}
./xmlchange -file env_build.xml    -id USE_ESMF_LIB   -val TRUE
#./xmlchange -file env_build.xml    -id ESMF_LIBDIR    -val ${nancy_scratch}/esmf-mpi

# num_tasks_per_instance = #tasks_node / #instances_node
# Bluefire: #tasks_node = 64 using SMT
#           #instances_node = 1-degree: 4 on lrg_ nodes, 2 on standard.
#                             2-degree: 16 on lrg_ nodes, 8 on standard
#           CAM5; to be determined, but no more than listed for CAM4
set num_tasks_per_node = 64
set num_tasks_per_instance = 16

# This is hard-wiring for the current (1/17/2011) multi-instance CESM restriction
# that all instances must be advanced simultaneously.  Work is underway to relax that.
@ total_nt = $num_instances * $num_tasks_per_instance
echo "total MPI tasks requested = $total_nt"

# Atm gets all the nodes and runs.
# The other components divide them up.
# ? ? ?  What about sglc?  It's a stub, and doesn't matter what pes are assigned to it.
# This algorithm figures out whether there are enough processors requested
# to run each component on whole nodes, or the components need to share some nodes.
# The hard-coded numbers (and ratios between them) are estimates; change them if you
# know better.
@ atm_pes  = $total_nt
@ large_small = $total_nt / (8 * $num_tasks_per_node)
if ($large_small > 0) then
   # Large_small > 0 means there are at least 8 nodes requested.
   # Allot whole nodes to the major components.
   @ cpl_pes  = ( $total_nt / (3 * $num_tasks_per_node) ) * $num_tasks_per_node
   @ cice_pes = ( $total_nt / (3 * $num_tasks_per_node) ) * $num_tasks_per_node
   @ docn_pes = $num_tasks_per_node
   @ lnd_pes  = $total_nt - ($cpl_pes + $cice_pes + $docn_pes)
else
   # 1/4 cpl,  1/4 cice, 1/8 docn, 3/8 lnd.  These will occupy fractions of nodes.
   @ cpl_pes  = $total_nt / 4
   @ cice_pes = $total_nt / 4
   @ docn_pes = $total_nt / 8
   @ lnd_pes  = $total_nt - ($cpl_pes + $cice_pes + $docn_pes)
endif

# first pe of each component (counted from 0)
@ atm_rootpe  = 0
@ cpl_rootpe  = 0
@ cice_rootpe = $cpl_rootpe  + $cpl_pes
@ docn_rootpe = $cice_rootpe + $cice_pes
@ lnd_rootpe  = $docn_rootpe + $docn_pes

echo "check pe counting; last_pe should = total_nt - 1"
@ last_pe = ( $lnd_rootpe + $lnd_pes ) - 1
echo $last_pe

echo "task layout"
echo "[$atm_rootpe ......................... ATM ............................. $atm_pes]"
echo "[$cpl_rootpe ... CPL ... $cice_rootpe ... ICE ... $docn_rootpe ... OCN ... $lnd_rootpe ... LND ... $total_nt]"
echo ""
echo "ATM gets $atm_pes"
echo "ICE gets $cice_pes"
echo "LND gets $lnd_pes"
echo "CPL gets $cpl_pes"
echo "OCN gets $docn_pes"
echo ""

./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $atm_pes
./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val $atm_rootpe
./xmlchange -file env_mach_pes.xml -id  NINST_ATM -val $num_instances

./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $cpl_pes
./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val $cpl_rootpe

./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $cice_pes
./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val $cice_rootpe
./xmlchange -file env_mach_pes.xml -id  NINST_ICE -val $num_instances

./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $docn_pes
./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val $docn_rootpe

./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $lnd_pes
./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val 1
./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val $lnd_rootpe
./xmlchange -file env_mach_pes.xml -id  NINST_LND -val $num_instances

./xmlchange -file env_conf.xml -id RUN_TYPE                -val startup
./xmlchange -file env_conf.xml -id RUN_STARTDATE           -val $run_startdate
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_FILENAME   -val $sst_dataset
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_YEAR_START -val $year_start
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_YEAR_END   -val $year_end
./xmlchange -file env_conf.xml -id CLM_CONFIG_OPTS         -val '-rtm off'
# The river transport model ON is useful only when using an active ocean or
# land surface diagnostics.

# Substantial archiving changes exist in the Tools/st_archive.sh script.
# Do not change the CALENDAR or the CONTINUE_RUN
# DOUT_S     is to turn on/off the short-term archiving
# DOUT_L_MS  is to store to the HPSS (formerly "MSS")

./xmlchange -file env_run.xml -id CONTINUE_RUN               -val FALSE
./xmlchange -file env_run.xml -id RESUBMIT                   -val $resubmit
./xmlchange -file env_run.xml -id STOP_OPTION                -val $stop_option
./xmlchange -file env_run.xml -id STOP_N                     -val $stop_n
./xmlchange -file env_run.xml -id CALENDAR                   -val GREGORIAN
./xmlchange -file env_run.xml -id DOUT_S_ROOT                -val ${archdir}
./xmlchange -file env_run.xml -id DOUT_S                     -val TRUE
./xmlchange -file env_run.xml -id DOUT_S_SAVE_INT_REST_FILES -val TRUE
./xmlchange -file env_run.xml -id DOUT_L_MS                  -val TRUE
./xmlchange -file env_run.xml -id DOUT_L_MSROOT              -val "csm/${case}"
./xmlchange -file env_run.xml -id DOUT_L_HTAR                -val FALSE

# ====================================================================
# Create namelist template: user_nl_cam, user_nl_clm
# ====================================================================

cd ${caseroot}

cat <<EOF >! user_nl_cam
&camexp
 inithist             = 'ENDOFRUN'
 div24del2flag        = 4
 empty_htapes         = .true.
 fincl1               = 'PHIS:I'
 nhtfrq               = -$stop_n
 iradae               = -$stop_n
 aerodep_flx_datapath = '${cesm_datadir}/atm/cam/chem/trop_mozart_aero/aero'
 aerodep_flx_file     = 'aerosoldep_monthly_1849-2006_1.9x2.5_c090803.nc'
 aerodep_flx_cycle_yr = 2000
 aerodep_flx_type     = 'CYCLICAL'
/
EOF

cat <<EOF >! user_nl_clm
&clmexp
  fatmgrid = '${cesm_datadir}/lnd/clm2/griddata/griddata_0.9x1.25_070212.nc'
  faerdep  = '${cesm_datadir}/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_rcp4.5_monthly_1849-2104_0.9x1.25_c100407.nc'
  outnc_large_files = .true.
  hist_empty_htapes = .true.
/
EOF

# ====================================================================
# Update source files if need be
#    Ideally, using DART will not require any source mods.
#    Until then, this script accesses source mods from a hard-wired location below.
#    Those may eventually be packaged into the DART repository.
#    If you have additional source mods, they will need to be merged into any DART
#    mods and put in the SourceMods subdirectory found in the 'case' directory.
# ====================================================================

\cp -rf ~thoar/${cesmtag}/SourceMods/* ${caseroot}/SourceMods/
if ( $status == 0) then
   echo "FYI - Local Source Modifications used for this case:"
   ls -lr ${caseroot}/SourceMods/*
else
   echo "FYI - No SourceMods for this case"
endif

# ====================================================================
# Configure
# ====================================================================

cd ${caseroot}

./configure -case

if ( $status != 0 ) then
   echo "ERROR: Case could not be configured."
   exit 2
endif

# ====================================================================
# Stage the required parts of DART in the caseroot directory.
# ====================================================================

cd ${caseroot}

if ("${reuse_existing_case}" == "false") then
   \mv Tools/st_archive.sh Tools/st_archive.sh.orig
endif
\cp -f ${DARTroot}/models/cam/shell_scripts/st_archive.sh Tools/.

# The cesm1_1_beta04 release had an error in that it did not
# provide the lt_archive.csh script, and the one in the repos
# did not have the -p flag, which is a good idea. So, for now ...
\cp -f ${cesmroot}/scripts/ccsm_utils/Tools/lt_archive.csh Tools/.

\cp -f ${DARTroot}/models/cam/shell_scripts/assimilate.csh .
\cp -f ${DARTroot}/models/cam/work/input.nml                .

# ====================================================================
# Update the scripts that build the namelists.
# The active components scripts need to support the multi-instance naming.
# ====================================================================

echo ''
echo 'Editing the Buildconf/{cam,cice,clm}.buildnml.csh files'
echo ''

cd ${caseroot}/Buildconf

cp -f  cam.buildnml.csh  cam.buildnml.csh.orig
cp -f cice.buildnml.csh cice.buildnml.csh.orig
cp -f  clm.buildnml.csh  clm.buildnml.csh.orig

# The CAM buildnml script only needs changing in one place.

ex cam.buildnml.csh <<ex_end
/cam_inparm/
/ncdata/
s;= '.*';= "cam_initial_\${atm_inst_counter}.nc";
wq
ex_end

# The CICE buildnml script only needs changing in one place.

ex cice.buildnml.csh <<ex_end
/setup_nml/
/ice_ic/
s;= '.*';= "ice_restart_\${ice_inst_counter}.nc";
wq
ex_end

# The CLM buildnml script needs changing in MULTIPLE places.

@ n = 1
while ($n <= $num_instances)
   set inst = `printf "%04d" $n`
   ex clm.buildnml.csh <<ex_end
/lnd_in_$inst/
/finidat/
s;= '.*';= "clm_restart_${n}.nc";
wq
ex_end
   @ n++
end

# ====================================================================
# The *.run script must be modified to call the DART assimilate script.
# The modifications are contained in a "here" document that MUST NOT
# expand the wildcards etc., before it is run. This is achieved by
# double-quoting the characters used to delineate the start/stop of
# the "here" document. No kidding. It has to be "EndOfText",
# not 'EndOfText' or EndOfText.
# ====================================================================

cd ${caseroot}

echo ''
echo 'Adding the call to assimilate.csh to the *.run script.'
echo ''

cat << "EndOfText" >! add_to_run.txt

# -------------------------------------------------------------------------
# START OF DART: if CESM finishes correctly (pirated from ccsm_postrun.csh);
# perform an assimilation with DART.
# -------------------------------------------------------------------------

set CplLogFile = `ls -1t cpl.log* | head -n 1`
if ($CplLogFile == "") then
 echo 'ERROR: Model did not complete - no cpl.log file present - exiting'
 echo 'ERROR: Assimilation will not be attempted.'
 exit -4
endif

grep 'SUCCESSFUL TERMINATION' $CplLogFile
if ( $status == 0 ) then
  ${CASEROOT}/assimilate.csh

  if ( $status == 0 ) then
     echo "`date` -- DART HAS FINISHED"
  else
     echo "`date` -- DART FILTER ERROR - ABANDON HOPE"
     exit -5
  endif
endif

# END OF DART BLOCK
# -------------------------------------------------------------------------

"EndOfText"

# Now that the "here" document is created,
# determine WHERE to insert it -- ONLY IF it is not already there.

grep "ABANDON HOPE" ${case}.${mach}.run
set STATUSCHECK = $status

if ( ${STATUSCHECK} == 0 ) then
   echo "DART block already present in ${case}.${mach}.run"
else if ( ${STATUSCHECK} == 1 ) then

   set MYSTRING = `grep --line-number "CSM EXECUTION HAS FINISHED" ${case}.${mach}.run`
   set MYSTRING = `echo $MYSTRING | sed -e "s#:# #g"`

   @ origlen = `cat ${case}.${mach}.run | wc -l`
   @ keep = $MYSTRING[1]
   @ lastlines = $origlen - $keep

   mv ${case}.${mach}.run ${case}.${mach}.run.orig

   head -n $keep      ${case}.${mach}.run.orig >! ${case}.${mach}.run
   cat                add_to_run.txt           >> ${case}.${mach}.run
   tail -n $lastlines ${case}.${mach}.run.orig >> ${case}.${mach}.run

else
   echo "ERROR in grep of ${case}.${mach}.run: aborting"
   echo "status was ${STATUSCHECK}"
   exit 10
endif

# ====================================================================
# Edit the run script to reflect project, queue, and wallclock
# ====================================================================

set PROJ=`grep BSUB $case.$mach.run | grep -e '-P' `
sed s/$PROJ[3]/$proj/ < $case.$mach.run >! temp
/bin/mv temp  $case.$mach.run

set TIMEWALL=`grep BSUB $case.$mach.run | grep -e '-W' `
sed s/$TIMEWALL[3]/$timewall/ < $case.$mach.run >! temp
/bin/mv temp  $case.$mach.run

set QUEUE=`grep BSUB $case.$mach.run | grep -e '-q' `
sed s/$QUEUE[3]/$queue/ < $case.$mach.run >! temp
/bin/mv temp  $case.$mach.run

chmod 0744 $case.$mach.run

# ====================================================================
# IMPORTANT: All resubmits must be type 'startup'.
# Change Tools/ccsm_postrun.csh line 83 to CONTINUE_RUN -val FALSE'
# ====================================================================

cd ${caseroot}/Tools

echo ''
echo 'Changing Tools/ccsm_postrun.csh such that all the resubmits are "startup",'
echo 'which means CONTINUE_RUN should be FALSE in ccsm_postrun.csh'
echo ''

ex ccsm_postrun.csh <<ex_end
/use COMP_RUN_BARRIERS as surrogate for timing run logical/
/CONTINUE_RUN/
s;TRUE;FALSE;
wq
ex_end

# ====================================================================
# build
# ====================================================================

cd ${caseroot}

echo ''
echo 'Building the case'
echo ''

./$case.$mach.build

if ( $status != 0 ) then
   echo "ERROR: Case could not be built."
   exit 3
endif

# ====================================================================
# Stage the required parts of DART in the execution root directory,
# now that EXEROOT exists.
# ====================================================================

foreach FILE ( filter cam_to_dart dart_to_cam )
   \cp -f ${DARTroot}/models/cam/work/${FILE} ${exeroot}/.
   if ( $status != 0 ) then
      echo "ERROR: ${DARTroot}/models/cam/work/${FILE} not copied to ${exeroot}"
      echo "ERROR: ${DARTroot}/models/cam/work/${FILE} not copied to ${exeroot}"
      exit 3
   endif
end

# ====================================================================
# Stage the restarts now that the run directory exists
# ====================================================================

# 20080801 ... /glade/proj3/DART/raeder/FV1deg_4.0/Exp1/obs_0000
# 20081031 ... /ptmp/thoar/restarts

set stagedir = /ptmp/thoar/restarts

echo ''
echo "Staging the restarts from {$stagedir}"
echo ''

@ n = 1
while ($n <= $num_instances)
  echo "Staging restarts for instance $n of $num_instances"
  cp --preserve=timestamps ${stagedir}/CAM/caminput_${n}.nc ${exeroot}/run/cam_initial_${n}.nc
  cp --preserve=timestamps ${stagedir}/CLM/clminput_${n}.nc ${exeroot}/run/clm_restart_${n}.nc
  cp --preserve=timestamps ${stagedir}/ICE/iceinput_${n}.nc ${exeroot}/run/ice_restart_${n}.nc
  @ n++
end

# This script may stage a prior_inflate_restart
# from $CESM_setup_dir  to  ${exeroot}/[prior,pos]_inflate_restart.YYYY-MM-DD-SSSSS
if (-f ${CESM_setup_dir}/prior_inflate_restart) then
   cp  ${CESM_setup_dir}/prior_inflate_restart \
     ${exeroot}/run/prior_inflate_restart.${run_startdate}-00000
   echo "${CESM_setup_dir}/prior_inflate_restart has been copied to "
   echo "${exeroot}/run/prior_inflate_restart.${run_startdate}-00000"
   echo "If that has the wrong state vector, you will need to replace it before running"
else
   echo "May need to stage a ${exeroot}/run/prior_inflate_restart.${run_startdate}-00000"
   echo "appropriate for this state vector in ${exeroot}/run."
   echo "You may make one with ${DARTroot}/models/cam/work/fill_inflation_restart"
endif

echo '================================================================================'
echo " If you're using DART's sampling error correction,"
echo " the assimilate.csh script will try to copy one from"
echo " ${DARTroot}/system_simulation/final_full_precomputed_tables/final_full.${num_instances}"
echo '================================================================================'

# ====================================================================
# What to do next
# ====================================================================

echo ''
echo 'Time to check the case.'
echo "cd ${caseroot}"
echo "Modify what you like in input.nml, make sure the observation directory"
echo 'names built in assimilate.csh match those on your system, submit the job:'
echo "./$case.$mach.submit"
echo ''

