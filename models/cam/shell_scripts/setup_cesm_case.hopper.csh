#!/bin/csh
## ====================================================================
##
## ---------------------
## Purpose
## ---------------------
## 
## ---------------------
## How to set the script
## ---------------------
## -- Copy this script in your directory
## -- Choice a case name (by changing "setenv case" ) and save the script as $case.csh
## -- Set the case options at the top of the script 
## -- If you have source mods, the script assumes that your mods are in: mods_$case. 
##    So, if you have source mods, create a subdirectory mods_$case that contains your mods.
##    If you don t have any source mods, the script creates an empty subdirectory mods_$case. 
## -- If you have namelist mods, you need to add them to the namelist template: user_nl_cam_$case
##    Set your namelist variables over there (without modifying the syntax to create user_nl_cam_$case
## -- Now, you are ready to go. Save your script and submit your run with the command: ./$case.csh  
##    The script creates your case, configure, compile and submit your job. 
## -- The script also creates a subdirectory (nml_$case) that contains your namelists.
## 
## ---------------------
## Important features
## ---------------------
## -- If you want to change something in your case other than the runtime settings above,
##    you need to delete everything and start the run from scratch. 

## ====================================================================
## ====  Set case options
## ==================================================================== 

## ======================
## case settings 
## ======================

setenv ccsmtag       cesm1_1_beta04
setenv case          F_AMIP8
#setenv compset       F
#setenv compset       F2000
setenv compset       F_AMIP_CN
setenv resolution    f09_f09   

setenv num_instances    80

## ======================
## configure settings
## ======================


# my home, and temp/scratch space
setenv my_home $HOME
setenv my_scratch $SCRATCH
setenv cesm_public /project/projectdirs/ccsm1

# for data files not in the public cems area
setenv my_datadir ${my_scratch}/cesm_datafiles

setenv run_startdate 2008-10-31
setenv sst_dataset $my_datadir/sst_HadOIBl_bc_0.9x1.25_1850_2011_c110307.nc
setenv year_start  1850
setenv year_end    2010


## ================================
##  define machines and directories 
## ================================

setenv mach hopp2                                       ;# machine - must match cesm names
setenv mytmp ${SCRATCH}                                 ;# your scratch/tmp space on this machine
setenv ccsmroot ${my_home}/${ccsmtag}                   ;# this is where the cesm code lives
#setenv ccsmroot ${cesm_public}/collections/${ccsmtag}  ;# this is where the cesm code lives
setenv caseroot ${my_home}/cases/$case                  ;## your cesm case directory
setenv lockroot ${my_home}/locked_cases                 ;## locked cases
setenv rundir ${mytmp}/${case}

## ======================
## runtime settings
## ======================

# each 'job' is 6 hours, so 8 resubmits should be 2 full days
setenv resubmit      8  
setenv stop_n        6
setenv stop_option   nhours

## ======================
## namelist variables
## ======================
## We create namelist templates that get copied once the case has been configured.

setenv this_dir `pwd`

cat <<EOF >! user_nl_cam_${case}
&camexp
 inithist                   = 'ENDOFRUN'
 div24del2flag              = 4
 bndtvghg                   = '${cesm_public}/inputdata/atm/cam/ggas/ghg_rcp45_1765-2500_c100405.nc'
 prescribed_ozone_datapath  = '${cesm_public}/inputdata/atm/cam/ozone'
 prescribed_ozone_file      = 'ozone_1.9x2.5_L26_1850-2015_rcp45_c101108.nc'
 prescribed_ozone_name      = 'O3'
 prescribed_ozone_type      = 'INTERP_MISSING_MONTHS'
 prescribed_volcaero_datapath = "$my_datadir"
 prescribed_volcaero_file   = 'CCSM4_volcanic_1850-2011_prototype1.nc'
 solar_data_file            = '${cesm_public}/inputdata/atm/cam/solar/SOLAR_TSI_Lean_1610-2140_annual_c100301.nc'
 prescribed_aero_datapath   = '${cesm_public}/inputdata/atm/cam/chem/trop_mozart_aero/aero'
 prescribed_aero_file       = 'aero_rcp45_v1_1.9x2.5_L26_1995-2105_c100316.nc'
 aerodep_flx_datapath       = "$my_datadir"
 aerodep_flx_file           = 'aerosoldep_rcp4.5_monthly_1849-2104_0.9x1.25_c100407.nc'
/
EOF

 #aerodep_flx_datapath       = '${cesm_public}/inputdata/atm/cam/chem/trop_mozart_aero/aero'
 #prescribed_volcaero_file   = 'CCSM4_volcanic_1850-2011_prototype1.nc'
 
cat <<EOF >! user_nl_clm_${case}
&clmexp
  fsurdat  = '${cesm_public}/inputdata/lnd/clm2/surfdata/surfdata_0.9x1.25_simyr1850_c091006.nc'
  fpftdyn  = '${cesm_public}/inputdata/lnd/clm2/surfdata/surfdata.pftdyn_0.9x1.25_rcp2.6_simyr1850-2100_c100323.nc'
  faerdep  = "$my_datadir/aerosoldep_rcp4.5_monthly_1849-2104_0.9x1.25_c100407.nc"
/
EOF

  #fsurdat  = '${cesm_public}/inputdata/lnd/clm2/surfdata/surfdata_0.9x1.25_simyr2000_c100505.nc'


## ======================
## job settings
## ======================

## project numbers available 
#setenv proj         93300315
#setenv timewall     6:00
#setenv queue        premium

## ====================================================================
## ====  End of case options
## ==================================================================== 

## ====================================================================
## Create a new case, configure, and build 
## for list of the cases: ./create_newcase -list
## ====================================================================

cd  $ccsmroot/scripts
./create_newcase -case $caseroot -mach $mach -res $resolution -compset $compset -skip_rundb  

cd $caseroot

./xmlchange -file env_build.xml    -id USE_ESMF_LIB   -val TRUE
./xmlchange -file env_build.xml    -id ESMF_LIBDIR    -val ${my_scratch}/esmf-mpi

# number of instances == ensemble size
./xmlchange -file env_mach_pes.xml -id NINST_ATM -val $num_instances
./xmlchange -file env_mach_pes.xml -id NINST_LND -val $num_instances
./xmlchange -file env_mach_pes.xml -id NINST_ICE -val $num_instances

# was 12 when nthreads was 1
set num_tasks_per_instance = 2
set nthreads = 6
@ total_nt = $num_instances * $num_tasks_per_instance

./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val "$total_nt"
./xmlchange -file env_mach_pes.xml -id NTHRDS_ATM -val "$nthreads"
./xmlchange -file env_mach_pes.xml -id ROOTPE_ATM -val '0'

./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val "$total_nt"
./xmlchange -file env_mach_pes.xml -id NTHRDS_LND -val "$nthreads"
./xmlchange -file env_mach_pes.xml -id ROOTPE_LND -val '0'

./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val "$total_nt"
./xmlchange -file env_mach_pes.xml -id NTHRDS_ICE -val "$nthreads"
./xmlchange -file env_mach_pes.xml -id ROOTPE_ICE -val '0'
   

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

# Have not sorted out the archiving just yet ... turning it off for now.
# DOUT_S     is to turn on/off the short-term archiving
# DOUT_L_MS  is to store to the HPSS (formerly "MSS")
./xmlchange -file env_run.xml -id DOUT_S      -val 'TRUE' 
./xmlchange -file env_run.xml -id DOUT_L_MS   -val 'FALSE' 
./xmlchange -file env_run.xml -id DOUT_S_SAVE_INT_REST_FILES      -val 'TRUE'

#------------------
##  Create namelist template: user_nl_cam 
#------------------

\mv -f ${this_dir}/user_nl_cam_{$case} ${caseroot}/user_nl_cam
\mv -f ${this_dir}/user_nl_clm_{$case} ${caseroot}/user_nl_clm

# FIXME: updated source files
cp ~/seq*F90 $caseroot/SourceMods/src.drv
echo updated source files:
ls -l $caseroot/SourceMods/src.drv/*

#------------------
## configure
#------------------

  cd $caseroot
  ./configure -case

#------------------
##  Stage a copy of the DART assimilate.csh script HERE
#------------------

  cd $caseroot

  setenv SCRIPTHOME $HOME/devel/models/cam/shell_scripts/
  ln -sf $SCRIPTHOME/assimilate.csh .
  ln -sf $SCRIPTHOME/st_archive.sh Tools/


#------------------
## update the namelists and scripts as needed
#------------------


echo
echo 'Editing the Buildconf/{cam,cice,clm}.buildnml.csh files'
echo

cd Buildconf
mv cam.buildnml.csh cam.buildnml.csh.orig
sed -e '/ ncdata/c  ncdata = "cam_initial_${atm_inst_counter}.nc" ' cam.buildnml.csh.orig >! cam.buildnml.csh
chmod 0755 cam.buildnml.csh

mv cice.buildnml.csh cice.buildnml.csh.orig
sed -e '/ ice_ic/c  ice_ic = "ice_restart_${ice_inst_counter}.nc" ' cice.buildnml.csh.orig >! cice.buildnml.csh
chmod 0755 cice.buildnml.csh

#mv clm.buildnml.csh clm.buildnml.csh.orig
#sed -e '/ finidat/c  finidat = "clm_restart_X.nc" ' clm.buildnml.csh.orig >! clm.buildnml.csh
cp -f clm.buildnml.csh clm.buildnml.csh.orig
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

chmod 0755 clm.buildnml.csh

cd ..

echo
echo 'Adding the call to assimilate.csh to the *.run script.'
echo 

# save original file minus the last 4 lines; insert the
# call to the assimilate script in front of those 4 lines.
# could also do this with an ex script

cat <<"EOF" >! add_to_run.csh

# -------------------------------------------------------------------------
# START OF DART STUFF:
# See if CSM finishes correctly (pirated from ccsm_postrun.csh) and if so,
# perform an assimilation with DART.
# -------------------------------------------------------------------------

set CplLogFile = `ls -1t cpl.log* | head -1`
if ($CplLogFile == "") then
 echo "ERROR: Model did not complete - no cpl.log file present - exiting"
 echo "ERROR: Assimilation will not be attempted."
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

# END OF DART STUFF
# -------------------------------------------------------------------------


"EOF"

mv ${case}.${mach}.run ${case}.${mach}.run.orig
set len = `cat ${case}.${mach}.run.orig | wc -l`
set lastlines = `grep -A 500 "CSM EXECUTION HAS FINISHED" ${case}.${mach}.run.orig | wc -l`
@ lastlines -= 2
@ keep = $len - $lastlines

head -$keep ${case}.${mach}.run.orig > ${case}.${mach}.run
cat ./add_to_run.csh >> ${case}.${mach}.run
tail -$lastlines ${case}.${mach}.run.orig >> ${case}.${mach}.run


cd Tools
cp ccsm_postrun.csh ccsm_postrun.orig
sed -e '/CONTINUE_RUN/s/TRUE/FALSE/' ccsm_postrun.orig >! ccsm_postrun.csh
chmod 0755 ccsm_postrun.csh 
cd ..

echo 'Changed Change Tools/ccsm_postrun.csh:83 to CONTINUE_RUN -val FALSE'
echo 'so all the resubmits are coldstarts.'


# if start date was not set above, set:
#   Buildconf/cpl.buildnml.csh: start_ymd, start_tod

#------------------
##  build 
#------------------
   
echo
echo 'Building the case'
echo

  cd $caseroot
  ./$case.$mach.build 

#------------------
##  restarts  - rundir does not exist until build runs
#------------------

echo
echo 'Staging the restarts'
echo

@ n = 1
while ($n <= $num_instances)

   echo "Staging restarts for instance $n of $num_instances"
   cp ${my_scratch}/tim_datafiles/CAM/caminput_${n}.nc ${rundir}/run/cam_initial_${n}.nc
   cp ${my_scratch}/tim_datafiles/CLM/clminput_${n}.nc ${rundir}/run/clm_restart_${n}.nc
   cp ${my_scratch}/tim_datafiles/ICE/iceinput_${n}.nc ${rundir}/run/ice_restart_${n}.nc

 @ n++
end

cp ${my_scratch}/tim_datafiles/DART/p*inflate_restart.* ${rundir}/run/

#------------------
##  Save another copy of the original 'initial' namelist so we 
## can tell if build changes anything?
#------------------
   

## ====================================================================
##  Continue run
## ====================================================================

# if ($create_build_case == false) then
#   echo "Found:" {$lockroot}/{$case}.locked
#   echo "set parameters for continue run"
#   cd $caseroot 
#   ./xmlchange  -file env_run.xml -id RESUBMIT     -val "$resubmit"
#   ./xmlchange  -file env_run.xml -id STOP_OPTION  -val $stop_option 
#   ./xmlchange  -file env_run.xml -id STOP_N       -val $stop_n 
#   ./xmlchange  -file env_run.xml -id CONTINUE_RUN -val TRUE  
# endif 

## ====================================================================
## Edit the run script - this is only for bluefire.
## ====================================================================

#set PROJ=`grep BSUB $case.$mach.run | grep -e '-P' `
#sed s/$PROJ[3]/$proj/ < $case.$mach.run >! temp
#/bin/mv temp  $case.$mach.run
#
#set TIMEWALL=`grep BSUB $case.$mach.run | grep -e '-W' `
#sed s/$TIMEWALL[3]/$timewall/ < $case.$mach.run >! temp
#/bin/mv temp  $case.$mach.run
#
#set QUEUE=`grep BSUB $case.$mach.run | grep -e '-q' `
#sed s/$QUEUE[3]/$queue/ < $case.$mach.run >! temp
#/bin/mv temp  $case.$mach.run

## ====================================================================
## Submit job  
## ====================================================================

echo 
echo 'Job ready to be submitted here'
echo "cd into $caseroot and run: ./$case.$mach.submit"
echo 

#  cd $caseroot
#  ./$case.$mach.submit

