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
# The real purpose of this set of notes is to record what is needed to configure
# and build a CESM instance that has CAM, CLM, and CICE as active components
# in a multi-instance configuration over a single data ocean ... etc.
# Despite looking like a script, it might best be used as a set of notes.
# ---------------------
# How to set the script
# ---------------------
# -- Copy this script into your directory
# -- Choose a case name (by changing "setenv case" ) and save the script as $case.csh
# -- Set the case options at the top of the script
# -- If you have source mods, the script assumes that your mods are in: mods_$case.
#    So, if you have source mods, create a subdirectory mods_$case that contains your mods.
#    If you don t have any source mods, the script creates an empty subdirectory mods_$case.
# -- If you have namelist mods, you need to add them to the namelist template: user_nl_cam
#    Set your namelist variables over there (without modifying the syntax to create user_nl_cam
# -- Now, you are ready to go. Save your script and submit your run with the command: ./$case.csh
#    The script creates your case, configure, compile and submit your job.
# -- The script also creates a subdirectory (nml_$case) that contains your namelists.
#
# ---------------------
# Important features
# ---------------------
#
# If you want to change something in your case other than the runtime settings,
# you need to delete everything and start the run from scratch.
#
# ./${CASENAME}.*.clean_build
# ./configure -cleanall
#
# ====================================================================
# === IMPORTANT modifications to the distribution code before ANYTHING 
# ====================================================================

# had to edit the following to remove the LSB_PJL_... word too long error
# cesm1_1_beta04/scripts/ccsm_utils/Machines/mkbatch.bluefire
# as long as OMP_NUM_THREADS == 1 ... the default is fine.

# ====================================================================
# ====  Set case options
# ====================================================================

setenv case                 startup0
setenv compset              F_2000
setenv ccsmtag              cesm1_1_beta04
setenv resolution           f09_f09
setenv num_instances        4
setenv reuse_existing_case  false

# ================================
# define machines and directories
# ================================
#
# mach            computer name
# cesm_datadir    location of public CESM data files
# cesm_public     location of public CESM code distributions
# caseroot        your (future) cesm case directory
# rundir          (future) run-time directory
# archdir         (future) short-term archive directory
# ccsmroot        location of the cesm code base
# DARTdir         location of DART executables, scripts and input

setenv mach         bluefire
setenv cesm_datadir /glade/proj3/cseg/inputdata
setenv cesm_public  /glade/proj3/cseg
setenv caseroot     /glade/user/${USER}/cases/${case}
setenv rundir       /glade/scratch/${USER}/${case}
setenv archdir      /glade/scratch/${USER}/archive/${case}
setenv DARTdir      /glade/home/thoar/svn/DART/dev

setenv ccsmroot     /glade/home/thoar/${ccsmtag}
#setenv ccsmroot     ${cesm_public}/collections/${ccsmtag}

# ======================
# configure settings
# ======================

setenv run_startdate 2008-10-31

setenv sst_dataset ${cesm_datadir}/atm/cam/sst/sst_HadOIBl_bc_0.9x1.25_1850_2011_c110307.nc
setenv year_start  1850
setenv year_end    2010

# ======================
# runtime settings
# ======================

setenv resubmit      0
setenv stop_n        6
setenv stop_option   nhours

# ======================
# job settings
# ======================

setenv proj         93300315
setenv timewall     1:00
setenv queue        lrg_regular

# ====================================================================
# Create the case.
# For list of the cases: ./create_newcase -list
# ====================================================================

# if reuse_existing_case is false and the directory does not exist, ...

if ("${reuse_existing_case}" == "false") then
   echo "removing old files from ${caseroot} and ${rundir}"
   \rm -fr ${caseroot}
   \rm -fr ${rundir}
   ${ccsmroot}/scripts/create_newcase -case ${caseroot} -mach ${mach} \
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

./xmlchange -file env_build.xml    -id EXEROOT        -val ${rundir}
./xmlchange -file env_build.xml    -id USE_ESMF_LIB   -val TRUE
#./xmlchange -file env_build.xml    -id ESMF_LIBDIR    -val ${nancy_scratch}/esmf-mpi

# The game here is that - for memory reasons - we want 4 ATMs running
# on a bluefire node ... capable of 64 tasks.
set num_tasks_per_atm_instance = 16
set nthreads = 1

# PE LAYOUT: 
#   total number of tasks  = 1280 
#   maximum threads per task = 1 
#   cpl ntasks = 320  nthreads=1 rootpe=  0 ninst=1 
#   cam ntasks =1280  nthreads=1 rootpe=  0 ninst=80 
#   clm ntasks = 576  nthreads=1 rootpe=704 ninst=80 
#   cice ntasks= 320  nthreads=1 rootpe=320 ninst=80 
#   docn ntasks=  64  nthreads=1 rootpe=640 ninst=1 
#   sglc ntasks=  64  nthreads=1 rootpe=  0 ninst=1 
#   
#   total number of hw pes = 1280 
#     cpl hw pe range ~ from 0 to 319 
#     cam hw pe range ~ from 0 to 1279 
#     clm hw pe range ~ from 704 to 1279 
#     cice hw pe range ~ from 320 to 639 
#     docn hw pe range ~ from 640 to 703 
#     sglc hw pe range ~ from 0 to 63 

if ( $num_instances == 4 ) then
   @ total_atm = 128
   @ total_ice = 32
   @ total_cpl = 32
   @ total_ocn = 32
   @ ice_start = $total_cpl
   @ ocn_start = $total_cpl + $total_ice
   @ lnd_start = $ocn_start + $total_ocn
   @ total_lnd = $total_atm - $lnd_start
else
   @ total_atm = $num_tasks_per_atm_instance * $num_instances
   @ total_ice = $num_tasks_per_atm_instance * $num_instances / 4
   @ total_cpl = $num_tasks_per_atm_instance * $num_instances / 4
   @ total_ocn = $num_tasks_per_atm_instance * 4
   @ ice_start = $total_cpl
   @ ocn_start = $total_cpl + $total_ice
   @ lnd_start = $ocn_start + $total_ocn
   @ total_lnd = $total_atm - $lnd_start
endif

echo "Node layout"
echo "[0 ......................... ATM ..................... $total_atm]"
echo "[0 ... CPL ... $total_cpl ... ICE ... $ocn_start ... OCN ... $lnd_start ... LND ... $total_atm]"

# Tony: "if f19_g16, i'd recommend using 320 for cice, 32 for
# docn, 320 for cpl, and then the rest for clm as a starting point."
# I'm changing the ratios ... not sure about using half a node
# [CPL][ICE][DOCN][LND] .... assuming we're using 20 nodes ...
# 1 node DOCN
# 5 nodes CPL
# 5 nodes ICE
# 9 nodes LND

./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $total_atm
./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $nthreads
./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val 0
./xmlchange -file env_mach_pes.xml -id  NINST_ATM -val $num_instances

./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $total_lnd
./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $nthreads
./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val $lnd_start
./xmlchange -file env_mach_pes.xml -id  NINST_LND -val $num_instances

./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $total_ice
./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $nthreads
./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val $ice_start
./xmlchange -file env_mach_pes.xml -id  NINST_ICE -val $num_instances

./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $total_cpl
./xmlchange -file env_mach_pes.xml -id NTHRDS_CPL -val $nthreads
./xmlchange -file env_mach_pes.xml -id ROOTPE_CPL -val 0

./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $total_ocn
./xmlchange -file env_mach_pes.xml -id NTHRDS_OCN -val $nthreads
./xmlchange -file env_mach_pes.xml -id ROOTPE_OCN -val $ocn_start
./xmlchange -file env_mach_pes.xml -id  NINST_OCN -val 1

./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $total_ocn
./xmlchange -file env_mach_pes.xml -id NTHRDS_GLC -val $nthreads
./xmlchange -file env_mach_pes.xml -id ROOTPE_GLC -val 0
./xmlchange -file env_mach_pes.xml -id  NINST_GLC -val 1

./xmlchange -file env_mach_pes.xml -id TOTALPES   -val $total_atm

./xmlchange -file env_conf.xml -id RUN_TYPE                -val startup
./xmlchange -file env_conf.xml -id RUN_STARTDATE           -val $run_startdate
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_FILENAME   -val $sst_dataset
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_YEAR_START -val $year_start
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_YEAR_END   -val $year_end
./xmlchange -file env_conf.xml -id CLM_CONFIG_OPTS         -val '-rtm off'

# The river transport model ON is useful only when using an active ocean or
# land surface diagnostics.

./xmlchange -file env_run.xml      -id CONTINUE_RUN     -val FALSE
./xmlchange -file env_run.xml      -id RESUBMIT         -val $resubmit
./xmlchange -file env_run.xml      -id STOP_OPTION      -val $stop_option
./xmlchange -file env_run.xml      -id STOP_N           -val $stop_n
./xmlchange -file env_run.xml      -id CALENDAR         -val GREGORIAN

# Substantial archiving changes exist in the Tools/st_archive.sh script.
# DOUT_S     is to turn on/off the short-term archiving
# DOUT_L_MS  is to store to the HPSS (formerly "MSS")
./xmlchange -file env_run.xml -id DOUT_S_ROOT                -val ${archdir}
./xmlchange -file env_run.xml -id DOUT_S                     -val TRUE
./xmlchange -file env_run.xml -id DOUT_S_SAVE_INT_REST_FILES -val TRUE
./xmlchange -file env_run.xml -id DOUT_L_MS                  -val FALSE 
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

# these cause problems when running with the full cesm:
#  nhtfrq                       = -12
# faerdep ... perhaps only used if CN modeling on ... 

cat <<EOF >! user_nl_clm
&clmexp
  fatmgrid = '${cesm_datadir}/lnd/clm2/griddata/griddata_0.9x1.25_070212.nc'
  faerdep  = '${cesm_datadir}/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_rcp4.5_monthly_1849-2104_0.9x1.25_c100407.nc'
  outnc_large_files = .true.
  hist_empty_htapes = .true.
/
EOF

# these cause problems when running with the full cesm:
#  hist_nhtfrq = -12

# ====================================================================
# Update source files if need be
# ====================================================================

\cp -rf ~thoar/${ccsmtag}/SourceMods/* ${caseroot}/SourceMods/
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

\mv Tools/st_archive.sh Tools/st_archive.sh.org
\cp -f ${DARTdir}/models/cam/shell_scripts/st_archive.sh Tools/st_archive.sh
# only needed for beta04 - fixed in more recent versions
\cp -f ${ccsmroot}/scripts/ccsm_utils/Tools/lt_archive.csh Tools/.

\cp -f ${DARTdir}/models/cam/shell_scripts/assimilate.startup.csh assimilate.csh
\cp -f ${DARTdir}/models/cam/work/input.nml    .
\cp -f ${DARTdir}/models/cam/work/filter       .
\cp -f ${DARTdir}/models/cam/work/cam_to_dart  .
\cp -f ${DARTdir}/models/cam/work/dart_to_cam  .

# ====================================================================
# Update the scripts that build the namelists.
# The active components scripts need to support the multi-instance naming.
# ====================================================================

echo ''
echo 'Editing the Buildconf/{cam,cice,clm}.buildnml.csh files'
echo ''

cd ${caseroot}/Buildconf

cp -f  cam.buildnml.csh  cam.buildnml.csh.org
cp -f cice.buildnml.csh cice.buildnml.csh.org
cp -f  clm.buildnml.csh  clm.buildnml.csh.org

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
# the "here" document. No kidding. It has to be "EndOfText"
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

set CplLogFile = `ls -1t cpl.log* | head -1`
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

grep --line-number "ABANDON HOPE" ${case}.${mach}.run
if ( $status == 0 ) then
   echo "DART block already present in ${case}.${mach}.run"
else

   set MYSTRING = `grep --line-number "CSM EXECUTION HAS FINISHED" ${case}.${mach}.run`
   set MYSTRING = `echo $MYSTRING | sed -e "s#:# #g"`

   @ orglen = `cat ${case}.${mach}.run | wc -l`
   @ keep = $MYSTRING[1]
   @ lastlines = $orglen - $keep 

   mv ${case}.${mach}.run ${case}.${mach}.run.orig

   head -$keep      ${case}.${mach}.run.orig >! ${case}.${mach}.run
   cat              add_to_run.txt           >> ${case}.${mach}.run
   tail -$lastlines ${case}.${mach}.run.orig >> ${case}.${mach}.run

endif

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
# Stage the restarts now that the run directory exists
# ====================================================================

# 20081031 ... /ptmp/thoar/restarts
# 20080801 ... /glade/proj3/DART/raeder/FV1deg_4.0/Exp1/obs_0000
set stagedir = /glade/proj3/DART/raeder/FV1deg_4.0/Exp1/obs_0000
set stagedir = /ptmp/thoar/restarts

echo ''
echo "Staging the restarts from {$stagedir}"
echo ''

@ n = 1
while ($n <= $num_instances)

   echo "Staging restarts for instance $n of $num_instances"
   cp --preserve=timestamps ${stagedir}/CAM/caminput_${n}.nc ${rundir}/run/cam_initial_${n}.nc
   cp --preserve=timestamps ${stagedir}/CLM/clminput_${n}.nc ${rundir}/run/clm_restart_${n}.nc
   cp --preserve=timestamps ${stagedir}/ICE/iceinput_${n}.nc ${rundir}/run/ice_restart_${n}.nc

 @ n++
end

echo 'If inflation is being used ... '
echo "must stage a ${rundir}/[prior,pos]_inflate_restart.YYYY-MM-DD-SSSSS"

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
# Submit job
# ====================================================================

echo ''
echo 'Time to check the case.'
echo "cd ${caseroot}"
echo "Modify what you like in input.nml, make sure the observation directory"
echo 'names built in assimilate.csh match those on your system, submit the job:'
echo "./$case.$mach.submit"
echo ''

