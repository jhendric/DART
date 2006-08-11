#!/bin/csh
#
# Data Assimilation Research Testbed -- DART
# Copyright 2004-2006, Data Assimilation Research Section 
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#
# <next three lines automatically updated by CVS, do not edit>
# $Id$
# $Source$
# $Name$
#
# Standard script for use in assimilation applications
# where the model advance is executed as a separate process.

# This script copies the necessary files into the temporary directory
# for a model run. It assumes that there is ${CENTRALDIR}/WRF directory
# where boundary conditions files reside.

# If the ensemble mean assim_model_state_ic_mean is present in the CENTRALDIR,
# it is converted to a WRF netCDF format.
# It is then used in update_wrf_bc the calculate the deviation from the mean.
# This deviation from the mean is then added at the end of the interval to
# calculate new boundary tendencies. The magnitude of the perturbation added
# at the end of the interval is controled by infl. The purpose is to increase
# time correlation at the lateral boundaries.

set infl = 0.0

set days_in_month = ( 31 28 31 30 31 30 31 31 30 31 30 31 )

set      myname = $0
set  CENTRALDIR = $1
set     element = $2
set    temp_dir = $3

set REMOVE = 'rm -rf'
set   COPY = 'cp -p'
set   MOVE = 'mv -f'

# Shell script to run the WRF model from DART input.
# where the model advance is executed as a separate process.

echo "starting ${myname} for ens member $element at "`date`
echo "CENTRALDIR is ${CENTRALDIR}"
echo "temp_dir is ${temp_dir}"

# Create a clean temporary directory and go there
${REMOVE} ${temp_dir}
mkdir -p  ${temp_dir}
cd        ${temp_dir}


# Copy or link the required files to the temp directory

ln -s ${CENTRALDIR}/input.nml .

ln -s ${CENTRALDIR}/RRTM_DATA .
ln -s ${CENTRALDIR}/LANDUSE.TBL .
ln -s ${CENTRALDIR}/VEGPARM.TBL .
ln -s ${CENTRALDIR}/SOILPARM.TBL .
ln -s ${CENTRALDIR}/GENPARM.TBL .
ln -s ${CENTRALDIR}/wrf.exe .
ln -s ${CENTRALDIR}/gribmap.txt .

# nfile is required when using mpi to run wrf.exe
# nfile is machine specific; ideally, it should be
# constructed by the script advance_ens.csh

hostname >! nfile
hostname >> nfile
###ln -s  ${CENTRALDIR}/nfile$element nfile
${COPY} ${CENTRALDIR}/wrfinput_d0? .
                   # Provides auxilliary info not avail. from DART state vector

if (  -e ${CENTRALDIR}/assim_model_state_ic_mean ) then
   ln -s ${CENTRALDIR}/assim_model_state_ic_mean dart_wrf_vector
   echo ".true." | ${CENTRALDIR}/dart_tf_wrf >& out.dart_to_wrf_mean
   ${COPY} wrfinput_d01 wrfinput_mean
endif

${MOVE} ${CENTRALDIR}/assim_model_state_ic$element dart_wrf_vector # ICs for run

# Convert DART to wrfinput

echo ".true." | ${CENTRALDIR}/dart_tf_wrf >& out.dart_to_wrf

${REMOVE} dart_wrf_vector

# The program dart_tf_wrf has created the file wrf.info.
# Time information is extracted from wrf.info.

set secday = `head -1 wrf.info`
set targsecs = $secday[1]
set targdays = $secday[2]
set targkey = `echo "$targdays * 86400 + $targsecs" | bc`

set secday = `head -2 wrf.info | tail -1`
set wrfsecs = $secday[1]
set wrfdays = $secday[2]
set wrfkey = `echo "$wrfdays * 86400 + $wrfsecs" | bc`

# If model blew up in the previous cycle, the member is now likely an outlier.
# Set infl = 0. to avoid further deterioration of the ensemble member.

if ( -e ${CENTRALDIR}/blown_${wrfdays}_${wrfsecs}.out ) then
   set MBLOWN = `cat ${CENTRALDIR}/blown_${wrfdays}_${wrfsecs}.out`
   set NBLOWN = `cat ${CENTRALDIR}/blown_${wrfdays}_${wrfsecs}.out | wc -l`
   set BLOWN = 0
   set imem = 1
   while ( $imem <= $NBLOWN )
      if ( $MBLOWN[$imem] == $element ) then
         @ BLOWN ++
      endif
      @ imem ++
   end
   if ( $BLOWN > 0 ) then
      set infl = 0.0
   endif
endif

# Find all BC's file available and sort them with "keys".

#--1st, check if LBCs are "specified" (in which case wrfbdy files are req'd)
set SPEC_BC = `grep specified ${CENTRALDIR}/namelist.input | grep true | cat | wc -l`

if ($SPEC_BC > 0) then
   ls ${CENTRALDIR}/WRF/wrfbdy_*_$element > bdy.list
else
   echo ${CENTRALDIR}/WRF/wrfbdy_${targdays}_${targsecs}_$element > bdy.list
endif

echo ${CENTRALDIR}/WRF/wrfbdy_ > str.name
sed 's/\//\\\//g' < str.name > str.name2
set STRNAME = `cat str.name2`
set COMMAND = s/`echo ${STRNAME}`//

sed $COMMAND < bdy.list > bdy.list2
sed 's/_/ /g' < bdy.list2 > bdy.list
set num_files = `cat bdy.list | wc -l`
set items = `cat bdy.list`
set ifile = 1
set iday = 1
set isec = 2
while ( $ifile <= $num_files )
   set key = `echo "$items[$iday] * 86400 + $items[$isec]" | bc`
   echo $key >> keys
   @ ifile ++
   set iday = `expr $iday \+ 3`
   set isec = `expr $isec \+ 3`
end
set keys = `sort keys`

set cal_date    = `head -3 wrf.info | tail -1`
set START_YEAR  = $cal_date[1]
set START_MONTH = $cal_date[2]
set START_DAY   = $cal_date[3]
set START_HOUR  = $cal_date[4]
set START_MIN   = $cal_date[5]
set START_SEC   = $cal_date[6]

set END_YEAR    = $cal_date[1]
set END_MONTH   = $cal_date[2]
set END_DAY     = $cal_date[3]
set END_HOUR    = $cal_date[4]
set END_MIN     = $cal_date[5]
set END_SEC     = $cal_date[6]

set MY_NUM_DOMAINS    = `head -4 wrf.info | tail -1`
set ADV_MOD_COMMAND   = `head -5 wrf.info | tail -1`

if ( `expr $END_YEAR \% 4` == 0 ) then
   set days_in_month[2] = 29
endif
if ( `expr $END_YEAR \% 100` == 0 ) then
   if ( `expr $END_YEAR \% 400` == 0 ) then
      set days_in_month[2] = 29
   else
      set days_in_month[2] = 28
   endif
endif

set ifile = 1
# Find the next BC's file available.

while ( $keys[${ifile}] <= $wrfkey )
   if ($ifile < $num_files ) then
      @ ifile ++
   else
      echo No boundary file available to move beyond
      echo ${START_YEAR}-${START_MONTH}-${START_DAY}_${START_HOUR}:${START_MIN}:${START_SEC}
      exit
   endif
end

###############################################################
# Advance the model with new BC until target time is reached. #
###############################################################

while ( $wrfkey < $targkey )

   set iday = `echo "$keys[$ifile] / 86400" | bc`
   set isec = `echo "$keys[$ifile] % 86400" | bc`

   # Copy the boundary condition file to the temp directory.

   ${COPY} ${CENTRALDIR}/WRF/wrfbdy_${iday}_${isec}_$element wrfbdy_d01

   ${COPY} ${CENTRALDIR}/WRF/wrflowinp_d01_${iday}_${isec} wrflowinp_d01

   if ( $targkey > $keys[$ifile] ) then
      set INTERVAL_SS = `echo "$keys[$ifile] - $wrfkey" | bc`
   else
      set INTERVAL_SS = `echo "$targkey - $wrfkey" | bc`
   endif
   set RUN_HOURS   = `expr $INTERVAL_SS \/ 3600`
   set REMAIN      = `expr $INTERVAL_SS \% 3600`
   set RUN_MINUTES = `expr $REMAIN \/ 60`
   set RUN_SECONDS = `expr $REMAIN \% 60`
   set INTERVAL_MIN = `expr $INTERVAL_SS \/ 60`

   @ END_SEC = $END_SEC + $INTERVAL_SS
   while ( $END_SEC >= 60 )
      @ END_SEC = $END_SEC - 60
      @ END_MIN ++
      if ($END_MIN >= 60 ) then
         @ END_MIN = $END_MIN - 60
         @ END_HOUR ++
      endif
      if ($END_HOUR >= 24 ) then
         @ END_HOUR = $END_HOUR - 24
         @ END_DAY ++
      endif
      if ($END_DAY > $days_in_month[$END_MONTH] ) then
         set END_DAY = 1
         @ END_MONTH ++
      endif
      if ($END_MONTH > 12 ) then
         set END_MONTH = 1
         @ END_YEAR ++

         if ( `expr $END_YEAR \% 4` == 0 ) then
            set days_in_month[2] = 29
         endif
         if ( `expr $END_YEAR \% 100` == 0 ) then
            if ( `expr $END_YEAR \% 400` == 0 ) then
               set days_in_month[2] = 29
            else
               set days_in_month[2] = 28
            endif
         endif

      endif
   end

   set END_SEC = `expr $END_SEC \+ 100`
   set END_SEC = `echo $END_SEC | cut -c2-3`
   set END_MIN = `expr $END_MIN \+ 100`
   set END_MIN = `echo $END_MIN | cut -c2-3`
   set END_HOUR = `expr $END_HOUR \+ 100`
   set END_HOUR = `echo $END_HOUR | cut -c2-3`
   set END_DAY = `expr $END_DAY \+ 100`
   set END_DAY = `echo $END_DAY | cut -c2-3`
   set END_MONTH = `expr $END_MONTH \+ 100`
   set END_MONTH = `echo $END_MONTH | cut -c2-3`

#-----------------------------------------------------------------------
# Update time control entries in the WRF namelist.input:
#-----------------------------------------------------------------------

${REMOVE} script.sed
cat > script.sed << EOF
 /run_hours/c\
 run_hours                  = ${RUN_HOURS}
 /run_minutes/c\
 run_minutes                = ${RUN_MINUTES}
 /run_seconds/c\
 run_seconds                = ${RUN_SECONDS}
 /start_year/c\
 start_year                 = ${START_YEAR}, ${START_YEAR}, ${START_YEAR}
 /start_month/c\
 start_month                = ${START_MONTH}, ${START_MONTH}, ${START_MONTH}
 /start_day/c\
 start_day                  = ${START_DAY}, ${START_DAY}, ${START_DAY}
 /start_hour/c\
 start_hour                 = ${START_HOUR}, ${START_HOUR}, ${START_HOUR}
 /start_minute/c\
 start_minute               = ${START_MIN}, ${START_MIN}, ${START_MIN}
 /start_second/c\
 start_second               = ${START_SEC}, ${START_SEC}, ${START_SEC}
 /end_year/c\
 end_year                   = ${END_YEAR}, ${END_YEAR}, ${END_YEAR}
 /end_month/c\
 end_month                  = ${END_MONTH}, ${END_MONTH}, ${END_MONTH}
 /end_day/c\
 end_day                    = ${END_DAY}, ${END_DAY}, ${END_DAY}
 /end_hour/c\
 end_hour                   = ${END_HOUR}, ${END_HOUR}, ${END_HOUR}
 /end_minute/c\
 end_minute                 = ${END_MIN}, ${END_MIN}, ${END_MIN}
 /end_second/c\
 end_second                 = ${END_SEC}, ${END_SEC}, ${END_SEC}
#  dart_tf_wrf is expecting only a single time per file
 /frames_per_outfile/c\
 frames_per_outfile         = 1, 1, 1,
EOF

 sed -f script.sed \
    ${CENTRALDIR}/namelist.input > namelist.input

# Update boundary conditions

   echo $infl | ${CENTRALDIR}/update_wrf_bc >& out.update_wrf_bc

   if ( -e rsl.out.integration ) then
      ${REMOVE} rsl.*
   endif

   ${ADV_MOD_COMMAND} >>& rsl.out.integration

   sleep 1

   set SUCCESS = `grep "wrf: SUCCESS COMPLETE WRF" rsl.* | cat | wc -l`
   if ($SUCCESS == 0) then
      echo $element >> ${CENTRALDIR}/blown_${targdays}_${targsecs}.out
   endif

if ( -e ${CENTRALDIR}/extract ) then
   if ( $element == 1 ) then
      ls wrfout_d0${MY_NUM_DOMAINS}_* > wrfout.list
      if ( -e ${CENTRALDIR}/psfc.nc ) then
         ${COPY} ${CENTRALDIR}/psfc.nc .
      endif
      echo `cat wrfout.list | wc -l` | ${CENTRALDIR}/extract
      ${MOVE} psfc.nc ${CENTRALDIR}/.
   endif
endif

   set dn = 1
   while ( $dn <= $MY_NUM_DOMAINS )
      ${MOVE} wrfout_d0${dn}_${END_YEAR}-${END_MONTH}-${END_DAY}_${END_HOUR}:${END_MIN}:${END_SEC} wrfinput_d0${dn}
      @ dn ++
   end

   ${REMOVE} wrfout*

   set START_YEAR  = $END_YEAR
   set START_MONTH = $END_MONTH
   set START_DAY   = $END_DAY
   set START_HOUR  = $END_HOUR
   set START_MIN   = $END_MIN
   set START_SEC   = $END_SEC
   set wrfkey = $keys[$ifile]
   @ ifile ++

end

##############################################
# At this point, the target time is reached. #
##############################################

# create new input to DART (taken from "wrfinput")
echo ".false." | ${CENTRALDIR}/dart_tf_wrf >& out.wrf_to_dart

${MOVE} dart_wrf_vector ${CENTRALDIR}/assim_model_state_ud$element

# Change back to working directory and get rid of temporary directory
cd ${CENTRALDIR}
#${REMOVE} ${temp_dir}
