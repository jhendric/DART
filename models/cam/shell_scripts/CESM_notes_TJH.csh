#!/bin/csh
# ====================================================================
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# ====================================================================
#
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
# -- If you have namelist mods, you need to add them to the namelist template: user_nl_cam_$case
#    Set your namelist variables over there (without modifying the syntax to create user_nl_cam_$case
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

# ====================================================================
# ====  Set case options
# ====================================================================

# ======================
# case settings
# ======================

#setenv case          F2000
#setenv compset       F_2000
setenv case          F_AMIP_CN
setenv compset       F_AMIP_CN
setenv ccsmtag       cesm1_1_beta04
setenv resolution    f09_f09
setenv num_instances 20

# ================================
# define machines and directories
# ================================

setenv mach bluefire                                ;# machine
setenv DARTdir /glade/home/thoar/svn/DART/dev       ;# DART executables, scripts and input
setenv cesm_public /glade/proj3/cseg/inputdata      ;# aka DIN_LOC_ROOT
setenv ccsmroot /glade/home/thoar/${ccsmtag}        ;# this is where the cesm code resides
setenv caseroot /glade/user/thoar/cases/${case}     ;# your (future) cesm case directory
setenv rundir /ptmp/thoar/${case}                   ;# (future) run-time directory
setenv archdir /glade/scratch/thoar/archive/${case} ;# (future) short-term archive directory

# ======================
# clear out previous builds
# ======================

echo "removing old files from ${caseroot} and ${rundir}"
\rm -fr ${caseroot}
\rm -fr ${rundir}

# ======================
# configure settings
# ======================

setenv run_startdate 2008-10-31

setenv sst_dataset ${cesm_public}/atm/cam/sst/sst_HadOIBl_bc_0.9x1.25_1850_2011_c110307.nc
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
setenv timewall     2:00
setenv queue        regular

# ======================
# namelist variables
# ======================
# Create namelist templates that get copied once the case has been configured.

setenv this_dir `pwd`

cat <<EOF >! user_nl_cam_${case}
&camexp
 inithist                     = 'ENDOFRUN'
 div24del2flag                = 4
 bndtvghg                     = '${cesm_public}/atm/cam/ggas/ghg_rcp45_1765-2500_c100405.nc'
 prescribed_ozone_datapath    = '${cesm_public}/atm/cam/ozone'
 prescribed_ozone_file        = 'ozone_1.9x2.5_L26_1850-2015_rcp45_c101108.nc'
 prescribed_ozone_name        = 'O3'
 prescribed_ozone_type        = 'INTERP_MISSING_MONTHS'
 prescribed_volcaero_datapath = '${cesm_public}/atm/cam/volc'
 prescribed_volcaero_file     = 'CCSM4_volcanic_1850-2011_prototype1.nc'
 solar_data_file              = '${cesm_public}/atm/cam/solar/SOLAR_TSI_Lean_1610-2140_annual_c100301.nc'
 prescribed_aero_datapath     = '${cesm_public}/atm/cam/chem/trop_mozart_aero/aero'
 prescribed_aero_file         = 'aero_rcp45_v1_1.9x2.5_L26_1995-2105_c100316.nc'
 aerodep_flx_datapath         = '${cesm_public}/atm/cam/chem/trop_mozart_aero/aero'
 aerodep_flx_file             = 'aerosoldep_rcp4.5_monthly_1849-2104_0.9x1.25_c100407.nc'
/
EOF

cat <<EOF >! user_nl_clm_${case}
&clmexp
  fpftdyn  = '${cesm_public}/lnd/clm2/surfdata/surfdata.pftdyn_0.9x1.25_rcp4.5_simyr1850-2100_c100406.nc'
  faerdep  = '${cesm_public}/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_rcp4.5_monthly_1849-2104_0.9x1.25_c100407.nc'
/
EOF

# ====================================================================
# ====  End of case options
# ====================================================================

# ====================================================================
# Create a new case, configure, and build.
# For list of the cases: ./create_newcase -list
# ====================================================================

cd  ${ccsmroot}/scripts
./create_newcase -case ${caseroot} -mach $mach -res $resolution -compset $compset -skip_rundb

if ( $status != 0 ) then
   echo "ERROR: Case could not be created."
   exit 1
endif

cd ${caseroot}

./xmlchange -file env_build.xml    -id USE_ESMF_LIB   -val TRUE
#./xmlchange -file env_build.xml    -id ESMF_LIBDIR    -val ${nancy_scratch}/esmf-mpi

set num_tasks_per_instance = 32
set nthreads = 1
@ total_nt = $num_instances * $num_tasks_per_instance

./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $total_nt
./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val $nthreads
./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val 0
./xmlchange -file env_mach_pes.xml -id  NINST_ATM -val $num_instances

./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $total_nt
./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val $nthreads
./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val 0
./xmlchange -file env_mach_pes.xml -id  NINST_LND -val $num_instances

./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $total_nt
./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val $nthreads
./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val 0
./xmlchange -file env_mach_pes.xml -id  NINST_ICE -val $num_instances

./xmlchange -file env_conf.xml -id RUN_TYPE                -val startup
./xmlchange -file env_conf.xml -id RUN_STARTDATE           -val $run_startdate
./xmlchange -file env_conf.xml -id RUN_REFDATE             -val $run_startdate
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_FILENAME   -val $sst_dataset
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_YEAR_START -val $year_start
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_YEAR_END   -val $year_end
./xmlchange -file env_conf.xml -id CLM_CONFIG_OPTS         -val '-rtm off'
#./xmlchange -file env_conf.xml -id CLM_CONFIG_OPTS         -val '-bgc cn -rtm off'

./xmlchange -file env_run.xml -id RESUBMIT    -val $resubmit
./xmlchange -file env_run.xml -id STOP_OPTION -val $stop_option
./xmlchange -file env_run.xml -id STOP_N      -val $stop_n
./xmlchange -file env_run.xml -id CALENDAR    -val GREGORIAN

# Substantial archiving changes exist in the Tools/st_archive.sh script.
# DOUT_S     is to turn on/off the short-term archiving
# DOUT_L_MS  is to store to the HPSS (formerly "MSS")
./xmlchange -file env_run.xml -id DOUT_S_ROOT                -val ${archdir}
./xmlchange -file env_run.xml -id DOUT_S                     -val TRUE
./xmlchange -file env_run.xml -id DOUT_S_SAVE_INT_REST_FILES -val TRUE
./xmlchange -file env_run.xml -id DOUT_L_MS                  -val TRUE
./xmlchange -file env_run.xml -id DOUT_L_HTAR                -val TRUE

# ====================================================================
# Create namelist template: user_nl_cam
# ====================================================================

\mv -f ${this_dir}/user_nl_cam_{$case} ${caseroot}/user_nl_cam
\mv -f ${this_dir}/user_nl_clm_{$case} ${caseroot}/user_nl_clm

# ====================================================================
# Update source files if need be
# ====================================================================

\cp ~/seq*F90 ${caseroot}/SourceMods/src.drv 
if ( $status == 0) then
   echo "FYI - Local Source Modifications used for this case:"
   ls -l ${caseroot}/SourceMods/src.drv/*
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
# Stage a copy of the DART assimilate.csh script HERE
# ====================================================================

cd ${caseroot}

\mv Tools/st_archive.sh Tools/st_archive.sh.org

\cp -f ${DARTdir}/models/cam/shell_scripts/assimilate.csh .
\cp -f ${DARTdir}/models/cam/shell_scripts/st_archive.sh Tools/st_archive.sh

# ====================================================================
# update the namelists and scripts as needed
# ====================================================================

echo ''
echo 'Editing the Buildconf/{cam,cice,clm}.buildnml.csh files'
echo ''

cd ${caseroot}/Buildconf

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

# The CLM buildnml script needs changing in MANY places.

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
# determine WHERE to insert it.

set MYSTRING = `grep --line-number "CSM EXECUTION HAS FINISHED" ${case}.${mach}.run`
set MYSTRING = `echo $MYSTRING | sed -e "s#:# #g"`

@ orglen = `cat ${case}.${mach}.run | wc -l`
@ keep = $MYSTRING[1]
@ lastlines = $orglen - $keep 

mv ${case}.${mach}.run ${case}.${mach}.run.orig

head -$keep      ${case}.${mach}.run.orig >! ${case}.${mach}.run
cat              add_to_run.txt           >> ${case}.${mach}.run
tail -$lastlines ${case}.${mach}.run.orig >> ${case}.${mach}.run

# ====================================================================
# IMPORTANT: All resubmits must be coldstarts.
# Change Tools/ccsm_postrun.csh line 83 to CONTINUE_RUN -val FALSE'
# ====================================================================

cd ${caseroot}/Tools

echo ''
echo 'Changing Tools/ccsm_postrun.csh such that all the resubmits are coldstarts.'
echo ''

ex Tools/ccsm_postrun.csh <<ex_end
/use COMP_RUN_BARRIERS as surrogate for timing run logical/
/CONTINUE_RUN/
s;TRUE;FALSE;
wq
ex_end

grep CONTINUE_RUN Tools/ccsm_postrun.csh
echo 'CONTINUE_RUN should be FALSE'

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

echo ''
echo 'Staging the restarts'
echo ''

@ n = 1
while ($n <= $num_instances)

   echo "Staging restarts for instance $n of $num_instances"
   cp /ptmp/thoar/restarts/CAM/caminput_${n}.nc ${rundir}/run/cam_initial_${n}.nc
   cp /ptmp/thoar/restarts/CLM/clminput_${n}.nc ${rundir}/run/clm_restart_${n}.nc
   cp /ptmp/thoar/restarts/ICE/iceinput_${n}.nc ${rundir}/run/ice_restart_${n}.nc

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
echo 'case is now ready to submit'
echo "cd into ${caseroot} and run: ./$case.$mach.submit"
echo ''

# exit

cd ${caseroot}
./$case.$mach.submit

