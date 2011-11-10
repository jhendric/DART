#!/bin/csh
# ====================================================================
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# ---------------------
# Purpose
# ---------------------
#
# The real purpose of this set of notes is to record what is needed to configure
# and build a CESM instance that has CAM, CLM, and CICE as active components
# in a multi-instance configuration over a single data ocean ... etc.
# Despite looking like a script, it might best be used as a set of notes.
# 
# ---------------------
# Important features
# ---------------------
#
# If you want to change something in your case other than the runtime settings,
# you need to delete everything and start the run from scratch. 
#
#
# ./${CASENAME}.*.clean_build
# ./configure -cleanall

# ====================================================================
# ====  Set case options
# ==================================================================== 

# ======================
# case settings 
# ======================

setenv ccsmtag       cesm1_1_beta04
setenv case          Farchive
setenv compset       F
setenv compset       F_AMIP_CN
setenv resolution    f09_f09   

# ======================
# configure settings
# ======================

setenv run_startdate 2008-10-31

setenv sst_dataset /glade/proj3/cseg/inputdata/atm/cam/sst/sst_HadOIBl_bc_0.9x1.25_1850_2011_c110307.nc
setenv year_start  1850
setenv year_end    2010

# ======================
# runtime settings
# ======================

setenv resubmit      0
setenv stop_n        6
setenv stop_option   nhours

# ======================
# namelist variables
# ======================
# We create namelist templates that get copied once the case has been configured.

setenv this_dir `pwd`

cat <<EOF >! user_nl_cam_${case}
&camexp
 inithist                   = 'ENDOFRUN'
 div24del2flag              = 4
 bndtvghg                   = '/fis/cgd/cseg/csm/inputdata/atm/cam/ggas/ghg_rcp45_1765-2500_c100405.nc'
 prescribed_ozone_datapath  = '/fis/cgd/cseg/csm/inputdata/atm/cam/ozone'
 prescribed_ozone_file      = 'ozone_1.9x2.5_L26_1850-2015_rcp45_c101108.nc'
 prescribed_ozone_name      = 'O3'
 prescribed_ozone_type      = 'INTERP_MISSING_MONTHS'
 prescribed_volcaero_file   = 'CCSM4_volcanic_1850-2011_prototype1.nc'
 solar_data_file            = '/fis/cgd/cseg/csm/inputdata/atm/cam/solar/SOLAR_TSI_Lean_1610-2140_annual_c100301.nc'
 prescribed_aero_datapath   = '/fis/cgd/cseg/csm/inputdata/atm/cam/chem/trop_mozart_aero/aero'
 prescribed_aero_file       = 'aero_rcp45_v1_1.9x2.5_L26_1995-2105_c100316.nc'
 aerodep_flx_datapath       = '/fis/cgd/cseg/csm/inputdata/atm/cam/chem/trop_mozart_aero/aero'
 aerodep_flx_file           = 'aerosoldep_rcp4.5_monthly_1849-2104_0.9x1.25_c100407.nc'
/
EOF

 
cat <<EOF >! user_nl_clm_${case}
&clmexp
  fpftdyn  = '/fis/cgd/cseg/csm/inputdata/lnd/clm2/surfdata/surfdata.pftdyn_0.9x1.25_rcp4.5_simyr1850-2100_c100406.nc'
  faerdep  = '/fis/cgd/cseg/csm/inputdata/atm/cam/chem/trop_mozart_aero/aero/aerosoldep_rcp4.5_monthly_1849-2104_0.9x1.25_c100407.nc'
/
EOF

# ======================
# job settings
# ======================

# project numbers available 
setenv proj         93300315
setenv timewall     6:00
setenv queue        premium

# ================================
#  define machines and directories 
# ================================

setenv mach bluefire                           ;# machine
setenv ptmp /ptmp/thoar                        ;# your ptmp on this machine
setenv ccsmroot /glade/home/thoar/${ccsmtag}   ;# this is where the cesm code lives
setenv caseroot /glade/user/thoar/cases/$case  ;# your cesm case directory
setenv lockroot /glade/user/thoar/locked_cases ;# locked cases
setenv rundir ${ptmp}/${case}

# ====================================================================
# ====  End of case options
# ==================================================================== 

# ====================================================================
# Create a new case, configure, and build 
# for list of the cases: ./create_newcase -list
# ====================================================================

cd  $ccsmroot/scripts
./create_newcase -case $caseroot -mach $mach -res $resolution -compset $compset -skip_rundb  

cd $caseroot

./xmlchange -file env_build.xml    -id USE_ESMF_LIB   -val TRUE

./xmlchange -file env_mach_pes.xml -id NINST_ATM -val 4
./xmlchange -file env_mach_pes.xml -id NINST_LND -val 4
./xmlchange -file env_mach_pes.xml -id NINST_ICE -val 4

#  ./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val '256'
#  ./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val '1'
#  ./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val '0'

#  ./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val '256'
#  ./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val '1'
#  ./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val '0'

#  ./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val '256'
#  ./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val '1'
#  ./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val '0'
   

./xmlchange -file env_conf.xml -id RUN_TYPE                -val startup
./xmlchange -file env_conf.xml -id RUN_STARTDATE           -val $run_startdate
./xmlchange -file env_conf.xml -id RUN_REFDATE             -val $run_startdate
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_FILENAME   -val $sst_dataset
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_YEAR_START -val $year_start    
./xmlchange -file env_conf.xml -id DOCN_SSTDATA_YEAR_END   -val $year_end  
./xmlchange -file env_conf.xml -id CLM_CONFIG_OPTS         -val '-bgc cn -rtm off'


./xmlchange -file env_run.xml -id RESUBMIT    -val $resubmit
./xmlchange -file env_run.xml -id STOP_OPTION -val $stop_option 
./xmlchange -file env_run.xml -id STOP_N      -val $stop_n 
./xmlchange -file env_run.xml -id CALENDAR    -val GREGORIAN

# Substantial archiving changes ...
# DOUT_S     is to turn on/off the short-term archiving
# DOUT_L_MS  is to store to the HPSS (formerly "MSS")
./xmlchange -file env_run.xml -id DOUT_S      -val 'TRUE' 
./xmlchange -file env_run.xml -id DOUT_S_SAVE_INT_REST_FILES      -val 'TRUE' 

#------------------
#  Create namelist template: user_nl_cam 
#------------------

\mv -f ${this_dir}/user_nl_cam_{$case} ${caseroot}/user_nl_cam
\mv -f ${this_dir}/user_nl_clm_{$case} ${caseroot}/user_nl_clm

#------------------
# configure
#------------------

  cd $caseroot
  ./configure -case

#------------------
#  Stage a copy of the DART assimilate.csh script HERE
#------------------

  cd $caseroot

  ln -s /glade/home/thoar/svn/DART/dev/models/cam/shell_scripts/assimilate.csh .
  ln -s /glade/home/thoar/svn/DART/dev/models/cam/shell_scripts/st_archive.sh Tools/st_archive.sh

#------------------
#  TJH FIXME EARLY EXIT
#------------------

echo ''
echo 'Edit the Buildconf/{cam,cice,clm}.buildnml.csh files'
echo '  cam.buildnml.csh      ncdata = "cam_initial_${atm_inst_counter}.nc"'
echo ' cice.buildnml.csh      ice_ic = "ice_restart_${ice_inst_counter}.nc"'
echo '  clm.buildnml.csh  hand edit all the clm_restart_x.nc instances in all HERE docs'
echo ''
echo 'Add the call to assimilate.csh to the *.run script.'
echo ''
echo 'IMPORTANT: Change Tools/ccsm_postrun.csh:83 to CONTINUE_RUN -val FALSE'
echo 'for all the resubmits to be coldstarts.'
echo ''

exit

#------------------
#  build 
#------------------
   
  cd $caseroot
  ./$case.$mach.build 

#------------------
# Stage the restarts now that the run directory exists
#------------------

@ n = 1
set Nens = 4
while ($n <= $Nens)

   echo "Staging restarts for instance $n of $Nens"
   cp /ptmp/thoar/restarts/CAM/caminput_${n}.nc ${rundir}/run/cam_initial_${n}.nc
   cp /ptmp/thoar/restarts/CLM/clminput_${n}.nc ${rundir}/run/clm_restart_${n}.nc
   cp /ptmp/thoar/restarts/ICE/iceinput_${n}.nc ${rundir}/run/ice_restart_${n}.nc

 @ n++
end

# ====================================================================
#  Continue run
# ====================================================================

# if ($create_build_case == false) then
#   echo "Found:" {$lockroot}/{$case}.locked
#   echo "set parameters for continue run"
#   cd $caseroot 
#   ./xmlchange  -file env_run.xml -id RESUBMIT     -val "$resubmit"
#   ./xmlchange  -file env_run.xml -id STOP_OPTION  -val $stop_option 
#   ./xmlchange  -file env_run.xml -id STOP_N       -val $stop_n 
#   ./xmlchange  -file env_run.xml -id CONTINUE_RUN -val TRUE  
# endif 

# ====================================================================
# Edit the run script
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

# ====================================================================
# Submit job  
# ====================================================================

  cd $caseroot
  ./$case.$mach.submit

