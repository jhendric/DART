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
# -- Either edit and run this script in the $DART/models/clm/shell_scripts
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
# ====  Set case options
# ====================================================================

# case will be used many ways;
#    directory and file names, both locally and on HPSS, and
#    script names; so consider it's length and information content.
# num_instances:  Number of ensemble members
# reuse_existing_case:
#    false; Remove $caseroot and $exeroot and rebuild
#    true;  configure -cleannamelist

setenv case                 clm_cesm1_1
setenv compset              I_2000_CN
setenv cesmtag              cesm1_1
setenv resolution           f19_f19
setenv num_instances        4
setenv reuse_existing_case  false

# ====================================================================
# define machines and directories
#
# mach            Computer name
# cesm_datadir    Root path of the public CESM data files
# cesmroot        Location of the cesm code base
#                 For cesm1_1_beta04 on yellowstone, MUST use 'thoar' value provided.
# DARTroot        Location of DART code tree.
#                    Executables, scripts and input in $DARTroot/models/dev/...
# caseroot        Your (future) cesm case directory, where this CESM+DART will be built.
#                    Preferably not a frequently scrubbed location.
#                 caseroot will be deleted if reuse_existing_case is false (below)
#                 So this script, and other useful things should be kept elsewhere.
# exeroot         (Future) Run-time directory; scrubbable, large amount of space needed.
# archdir         (Future) Short-term archive directory
#                    until the long-term archiver moves it to permanent storage.
# ====================================================================

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
# ====================================================================

setenv stream_year_first 2000
setenv stream_year_last  2000
setenv stream_year_align 2000
setenv refyear     2000
setenv refmon      01
setenv refday      31
setenv run_reftod  00000
setenv run_refdate $refyear-$refmon-$refday

# ====================================================================
# runtime settings --  How many assimilation steps will be done after this one
#
# stop_option   Units for determining the forecast length between assimilations
#               Changing stop_option requires changes to user_nl_clm below.
# stop_n        Number of time units in the forecast
# ====================================================================

setenv resubmit      4
setenv stop_option   nhours
setenv stop_n        24

# ====================================================================
# job settings
#
# timewall   can be changed during a series by changing the ${case}.run
# queue      can be changed during a series by changing the ${case}.run
#
# TJH How many f19_f19 CLM instances can fit on 1 'regular' node?
# ====================================================================

setenv proj         P86850054
setenv timewall     0:29
setenv queue        premium

# ====================================================================
# set these standard commands based on the machine you are running on.
# ====================================================================

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


# ====================================================================
# Create the case.
#
# For list of the pre-defined cases: ./create_newcase -list
# To create a variant case, see the CESM documentation and carefully
# incorporate any needed changes into this script.
# ====================================================================

if ("${reuse_existing_case}" == "false") then
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
      exit 1
   endif
else
   cd ${caseroot}
   ./configure -cleannamelist
   ./configure -cleanmach

endif

# ====================================================================
# Configure the case.
# ====================================================================

cd ${caseroot}

# This is just for debug purposes
foreach FILE ( *xml )
   cp $FILE ${FILE}.original
end


# The river transport model ON is useful only when using an active ocean or
# land surface diagnostics. The biogeochemistry should also be turned (back) on.
# If you have 'CN' in the compset name, CLM_CONFIG_OPTS defaults properly.
# since we are turning off the RTM, we need to turn back on "the right thing".
# ./xmlchange -file env_build.xml -id CLM_CONFIG_OPTS -val '-rtm off -bgc cn'
#
./xmlchange CLM_CONFIG_OPTS='-bgc cn'
./xmlchange ROF_GRID='null'

@ total_nt = 128
@ atm_pes  = $total_nt
@ cpl_pes  = $total_nt / 8
@ ice_pes  = $total_nt / 8
@ ocn_pes  = $total_nt / 8
@ rof_pes  = $total_nt / 8
@ lnd_pes  = $total_nt - ($cpl_pes + $ice_pes + $ocn_pes + $rof_pes)

echo "task partitioning ..."
echo ""
echo "ATM gets $atm_pes"
echo "CPL gets $cpl_pes"
echo "ICE gets $ice_pes"
echo "OCN gets $ocn_pes"
echo "LND gets $lnd_pes"
echo ""

./xmlchange NTHRDS_ATM=1,NTASKS_ATM=$atm_pes,NINST_ATM=$num_instances
./xmlchange NTHRDS_ICE=1,NTASKS_ICE=$ice_pes,NINST_ICE=1
./xmlchange NTHRDS_OCN=1,NTASKS_OCN=$ocn_pes,NINST_OCN=1
./xmlchange NTHRDS_GLC=1,NTASKS_GLC=$ocn_pes,NINST_GLC=1
./xmlchange NTHRDS_LND=1,NTASKS_LND=$lnd_pes,NINST_LND=$num_instances
./xmlchange NTHRDS_ROF=1,NTASKS_ROF=$rof_pes,NINST_ROF=$num_instances
./xmlchange NTHRDS_CPL=1,NTASKS_CPL=$cpl_pes

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
./xmlchange RUN_REFCASE=${case}
./xmlchange RUN_REFDATE=$run_refdate
./xmlchange RUN_REFTOD=$run_reftod
./xmlchange BRNCH_RETAIN_CASENAME=TRUE
./xmlchange GET_REFCASE=FALSE

./xmlchange STOP_OPTION=$stop_option
./xmlchange STOP_N=$stop_n
./xmlchange CONTINUE_RUN=FALSE
./xmlchange RESUBMIT=$resubmit

./xmlchange DOUT_S=FALSE
./xmlchange DOUT_S_ROOT=${archdir}
./xmlchange DOUT_S_SAVE_INT_REST_FILES=FALSE
./xmlchange DOUT_L_MS=FALSE
./xmlchange DOUT_L_MSROOT="csm/${case}"
./xmlchange DOUT_L_HTAR=FALSE

./xmlchange DATM_MODE=CPLHIST3HrWx
./xmlchange DATM_CPLHIST_CASE=$case
./xmlchange DATM_CPLHIST_YR_ALIGN=$refyear
./xmlchange DATM_CPLHIST_YR_START=$refyear
./xmlchange DATM_CPLHIST_YR_END=$refyear

#./xmlchange USE_ESMF_LIB=TRUE
#./xmlchange ESMF_LIBDIR=${nancy_scratch}/esmf-mpi

./xmlchange CALENDAR=GREGORIAN
./xmlchange EXEROOT=${exeroot}

# ====================================================================
# Configure
# ====================================================================

./cesm_setup

# ====================================================================
# Create namelist template: user_nl_clm
# Example user_nl_clm namelist adding and removing fields on primary history file
# hist_fincl1 = 'COSZEN', 'DECL'
# hist_fexcl1 = 'TG', 'TV', 'TSOI', 'H2OSOI'
# DART needs the lon,lat,levgrnd,lonatm,latatm,lonrof,latrof DIMENSION
# information from the .h0. history file - nothing else.
#
# hist_empty_htapes = .true.     suppresses the creation of all history files
# hist_fincl1 = 'TG',            except the first one, which will have one variable
# hist_nhtfrq = -$stop_n,        create one every $stop_n HOURS
# hist_mfilt  =  1,              with precisely one day in it
# hist_avgflag_pertape = 'I'     use instantaneous values - no average
#
# The fincl2 history tape has the half-hourly flux tower observations.
# The observation operators in obs_def_tower_mod.f90
# are going to read from the .h1. history file for these values.
# ====================================================================

@ inst = 1
while ($inst <= $num_instances)

   set instance  = `printf %04d $inst`

   set fname = "user_nl_datm_$instance"

   echo "dtlimit = 1.5, 1.5, 1.5"                    >> $fname
   echo "fillalgo = 'nn', 'nn', 'nn'"                >> $fname
   echo "fillmask = 'nomask','nomask','nomask'"      >> $fname
   echo "mapalgo = 'bilinear','bilinear','bilinear'" >> $fname
   echo "mapmask = 'nomask','nomask','nomask'"       >> $fname
   echo "streams = 'datm.streams.txt.CPLHIST3HrWx.Solar_$instance             $stream_year_align $stream_year_first $stream_year_last'," >> $fname
   echo "          'datm.streams.txt.CPLHIST3HrWx.Precip_$instance            $stream_year_align $stream_year_first $stream_year_last'," >> $fname
   echo "          'datm.streams.txt.CPLHIST3HrWx.nonSolarNonPrecip_$instance $stream_year_align $stream_year_first $stream_year_last'"  >> $fname
   echo "taxmode = 'cycle','cycle','cycle'"          >> $fname
   echo "tintalgo = 'coszen','coszen','coszen'"      >> $fname
   echo "restfils = 'unset'"                         >> $fname
   echo "restfilm = 'unset'"                         >> $fname

   # Customize the land namelists
   # Initially, just use the default restart files ... then move to
   # something ... specify in finidat?

   set fname = "user_nl_clm_$instance"

   set stagedir = /glade/scratch/afox/bptmp/MD_40_PME/run/MD_40_PME


#  echo "finidat = '${stagedir}.clm2_$instance.r.${run_refdate}-${run_reftod}.nc'" >> $fname
   echo "hist_empty_htapes = .false."                >> $fname
   echo "hist_fincl1 = 'NEP'"                        >> $fname
   echo "hist_fincl2 = 'NEP','FSH','EFLX_LH_TOT_R'"  >> $fname
   echo "hist_nhtfrq = -$stop_n,1,"                  >> $fname
   echo "hist_mfilt  = 1,48"                         >> $fname
   echo "hist_avgflag_pertape = 'A','A'"             >> $fname

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

./preview_namelists

# This gives us the right number of stream files

foreach FILE (CaseDocs/*streams*)
   set FNAME = $FILE:t

   switch ( ${FNAME} )
      case *presaero*:
         echo "Skipping prescribed aerosol stream txt file." 
         breaksw
      default:  
         cp -v $FILE user_$FNAME
         chmod 644   user_$FNAME
         breaksw
   endsw

end

# This gives us the right number of stream files

foreach FILE (user*streams*)
   set FNAME = $FILE
   set name_parse = `echo $FNAME:q | sed 's/\_/ /g'`

   if ($#name_parse == 3) then

      set instance = $name_parse[3]
      if (-e $DARTroot/models/clm/shell_scripts/*$name_parse[2]*template) then
         echo "Copying over DART template for $FNAME and changing instances"
         cp $DARTroot/models/clm/shell_scripts/*$name_parse[2]*template $FNAME
         sed s/ninst/$instance/g $FNAME >! out
         mv out $FNAME
      else
         echo "Looking for multi-instance template for $FNAME"
      endif

   else 

      echo "Looking for template for $FNAME"
      if (-e $DARTroot/models/clm/shell_scripts/${FNAME}_template) then
         echo "Copying over the DART template for $FNAME"
         cp  $DARTroot/models/clm/shell_scripts/${FNAME}_template $FNAME
      else
         echo "WARNING: cannot find DART template for $FNAME"
         echo "WARNING: cannot find DART template for $FNAME"
         exit 7
      endif

   endif
end

./preview_namelists

# ====================================================================
# Update source files if need be
#    Ideally, using DART will not require any source mods.
#    Until then, this script accesses source mods from a hard-wired location below.
#    Those may eventually be packaged into the DART repository.
#    If you have additional source mods, they will need to be merged into any DART
#    mods and put in the SourceMods subdirectory found in the 'case' directory.
# ====================================================================

# this one needs a recursive copy to get all the files in the subdirs
#  ${COPY} -r  ~thoar/${cesmtag}/SourceMods/* ${caseroot}/SourceMods/
#if ( $status == 0) then
#   echo "FYI - Local Source Modifications used for this case:"
#   ls -lr ${caseroot}/SourceMods/*
#else
#   echo "FYI - No SourceMods for this case"
#endif

# ====================================================================
# build
# ====================================================================

echo ''
echo 'Building the case'
echo ''

./$case.build

if ( $status != 0 ) then
   echo "ERROR: Case could not be built."
   exit 3
endif

# ====================================================================
# Stage the required parts of DART in the caseroot directory.
# ====================================================================

if ("${reuse_existing_case}" == "false") then
   ${MOVE} Tools/st_archive.sh Tools/st_archive.sh.orig
endif
${COPY} ${DARTroot}/models/clm/shell_scripts/st_archive.sh Tools/
# ${COPY} ${DARTroot}/models/clm/shell_scripts/datm.buildnml.csh Buildconf/

${COPY} ${DARTroot}/models/clm/work/input.nml                .

# ====================================================================
# The *.run script must be modified to call the DART assimilate script.
# The modifications are contained in a "here" document that MUST NOT
# expand the wildcards etc., before it is run. This is achieved by
# double-quoting the characters used to delineate the start/stop of
# the "here" document. No kidding. It has to be "EndOfText",
# not 'EndOfText' or EndOfText.
# ====================================================================

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

   mv ${case}.run ${case}.run.orig

   head -n $keep      ${case}.run.orig >! ${case}.run
   cat                add_to_run.txt   >> ${case}.run
   tail -n $lastlines ${case}.run.orig >> ${case}.run

endif

# ====================================================================
# Edit the run script to reflect project, queue, and wallclock
# ====================================================================

echo ''
echo 'Updating the run script to set project number, wallclock, and queue.'
echo ''

set PROJ=`grep BSUB $case.run | grep -e '-P' `
sed s/$PROJ[3]/$proj/ < $case.run >! temp
${MOVE} temp  $case.run

set TIMEWALL=`grep BSUB $case.run | grep -e '-W' `
sed s/$TIMEWALL[3]/$timewall/ < $case.run >! temp
${MOVE} temp  $case.run

set QUEUE=`grep BSUB $case.run | grep -e '-q' `
sed s/$QUEUE[3]/$queue/ < $case.run >! temp
${MOVE} temp  $case.run

chmod 0744 $case.run

# ====================================================================
# Stage the required parts of DART in the execution root directory,
# now that EXEROOT exists.
# ====================================================================

foreach FILE ( filter clm_to_dart dart_to_clm input.nml )
   ${COPY} ${DARTroot}/models/clm/work/${FILE} ${exeroot}/
   if ( $status != 0 ) then
      echo "ERROR: ${DARTroot}/models/clm/work/${FILE} not copied to ${exeroot}"
      echo "ERROR: ${DARTroot}/models/clm/work/${FILE} not copied to ${exeroot}"
      exit 3
   endif
end

${COPY} ${DARTroot}/models/clm/shell_scripts/assimilate.csh  assimilate.csh

# ====================================================================
# Stage the restarts now that the run directory exists
# ====================================================================
#
# obs sequences files: /ptmp/yfzhang/Obs_seqs

# 20000501 ... /ptmp/afox/MD_40_PME/run
set stagedir = /glade/scratch/afox/bptmp/MD_40_PME/run

# 20021101 ... /ptmp/yfzhang/inputdata_cam/lnd/clm2/initdata
# set stagedir = /ptmp/yfzhang/inputdata_cam/lnd/clm2/initdata

echo ''
echo "Copying the restart files from ${stagedir}"
echo ''

@ n = 1
while ($n <= $num_instances)

#   echo "Staging restarts for instance $n of $num_instances"

#  set LANDFILE = `printf ${stagedir}/init1998.clm2_%04d.r.2002-11-01-00000.nc $n`
   set LANDFILE = `printf ${stagedir}/MD_40_PME.clm2_%04d.r.2000-01-31-00000.nc $n`
   set LND_RESTART_FILENAME = `printf "${case}.clm2_%04d.r.%04d-%02d-%02d-%05d.nc" $n $refyear $refmon $refday $run_reftod`

#   ${COPY} ${LANDFILE} ${exeroot}/run/${LND_RESTART_FILENAME}

 @ n++
end

# ====================================================================
# What to do next
# ====================================================================

echo ''
echo 'Time to check the case.'
echo "cd into ${caseroot}"
echo 'Modify what you like in input.nml, make sure the observation directory'
echo 'names set in assimilate.csh match those on your system, and submit'
echo 'the CESM job by running:'
echo "./$case.submit"
echo ''

