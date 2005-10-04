! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
PROGRAM wrf_3dvar_tf_dart

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

use         types_mod, only : r8, missing_r8, missing_data, DEG2RAD, earth_radius
use     utilities_mod, only : open_file, close_file, file_exist, initialize_utilities, &
                              finalize_utilities, register_module, logfileunit, E_MSG, E_ERR, &
                              error_handler
use  obs_sequence_mod, only : obs_type, obs_sequence_type, init_obs_sequence, insert_obs_in_seq, &
                              set_copy_meta_data, set_qc_meta_data, write_obs_seq, assignment(=), &
                              init_obs, static_init_obs_sequence, set_obs_def, set_obs_values, set_qc
use       obs_def_mod, only : set_obs_def_location, set_obs_def_error_variance, &
                              set_obs_def_kind, set_obs_def_time, set_obs_def_key, &
                              obs_def_type, DOPPLER_RADIAL_VELOCITY, RADAR_REFLECTIVITY
use      location_mod, only : location_type, set_location
use  time_manager_mod, only : time_type, set_date, set_calendar_type, GREGORIAN
use obs_def_radar_mod, only : set_rad_vel

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: seq
type(obs_type)          :: obs, prev_obs
type(obs_def_type)      :: obs_def
type(location_type)     :: location
type(time_type)         :: time

INTEGER           :: iunit, ierr, iost, io

character(len=19) :: radar_name
real(r8)          :: radar_long, radar_lat, radar_elev
real(r8)          :: orientation(3)
character(len=19) :: start_scan_date
integer           :: num_prof, max_levels

character(len=80) :: dummy

character(len=12) :: platform_name, err_string, nml_string
integer           :: year, month, day, hours, minutes, seconds
real(r8)          :: lat,lon,elv
integer           :: levels

integer           :: ii, num_Radar, rv_qc, rf_qc, key, it
integer           :: num_obs, num_copies, num_qc, max_num_obs

real(r8)          :: height, rv_inv, rf_inv, rv_error, rf_error
real(r8)          :: obs_value(1), rstatus(1,1)
real(r8)          :: h, spath, x, y, rad_lon, rad_lat, obs_lon, obs_lat
real(r8)          :: rgate, raz, elev_rad, elev_obs, ae

!-----------------------------------------------------------------------------
! Namelist with default values
!
character(len = 129) :: wrf_3dvar_file        = 'qc_radr_3dvar_2002083100.dat', &
                        obs_seq_out_file_name = 'obs_seq.out'
integer              :: calendar_type         = GREGORIAN

namelist /wrf_3dvar_tf_dart_nml/ wrf_3dvar_file, obs_seq_out_file_name, calendar_type

!------------------------------------------------------------------------------

call initialize_utilities('wrf_3dvar_tf_dart')
call register_module(source, revision, revdate)
!write(logfileunit,*)'STARTING wrf_3dvar_tf_dart ...'
!call error_handler(E_MSG,'wrf_3dvar_tf_dart','STARTING ...',source,revision,revdate)

! Begin by reading the namelist input
if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   read(iunit, nml = wrf_3dvar_tf_dart_nml, iostat = io)
   if(io /= 0) then
      ! A non-zero return means a bad entry was found for this namelist
      ! Reread the line into a string and print out a fatal error message.
      BACKSPACE iunit
      read(iunit, '(A)') nml_string
      write(err_string, *) 'INVALID NAMELIST ENTRY: ', trim(adjustl(nml_string))
      call error_handler(E_ERR, 'wrf_3dvar_tf_dart:&wrf_3dvar_tf_dart_nml problem', &
                         err_string, source, revision, revdate)
   endif
   call close_file(iunit)
endif

! Record the namelist values used for the run ...
call error_handler(E_MSG,'wrf_3dvar_tf_dart','wrf_3dvar_tf_dart_nml values are',' ',' ',' ')
write(logfileunit, nml=wrf_3dvar_tf_dart_nml)
write(     *     , nml=wrf_3dvar_tf_dart_nml)

call set_calendar_type(calendar_type)

iunit = open_file(wrf_3dvar_file, action = 'read')

! -------------------------------------------------------------------
! Initialize the counters:

num_Radar = 0
num_obs = 0
key = 0

!-----------------------------------------------------------------------------!
! Read the header of a MM5 3D-VAR 2.0 Radar observation file
!-----------------------------------------------------------------------------!

READ (UNIT = iunit, IOSTAT = iost, &
     FMT = '(A19,F8.3,2X,F8.3,F10.1,2X,A19,2I6)' ) &
     radar_name, radar_long, radar_lat, radar_elev, &
     start_scan_date, &
     num_prof, max_levels

READ (UNIT = iunit, IOSTAT = iost, FMT = '(A80)' ) dummy
READ (UNIT = iunit, IOSTAT = iost, FMT = '(A80)' ) dummy

call static_init_obs_sequence()

max_num_obs = num_prof*max_levels*2.0_r8
num_copies = 1
num_qc = 1

! Initialize an obs_sequence structure
call init_obs_sequence(seq, num_copies, num_qc, max_num_obs)

call set_copy_meta_data(seq, 1, 'observations')
call set_qc_meta_data(seq, 1, 'missing')

call init_obs(obs, num_copies, num_qc)

!  READ FORMATS
!  ------------

!  LOOP OVER RECORDS
!  -----------------

reports: &
     DO
!     READ STATION GENERAL INFO
!     =============================

READ (UNIT = iunit, IOSTAT = iost, &
     FMT = '(A12,3X,I4,5(A1,I2),2X,2(F12.3,2X),F8.1,2X,I6)' ) &
     platform_name,  &
     year, dummy, month, dummy, day, dummy, hours, dummy, minutes, dummy, seconds, &
     lat,       &
     lon,       &
     elv,       &
     levels

time = set_date(year, month, day, hours, minutes, seconds)

IF (iost /= 0) THEN
   WRITE (0,'(/,A,I3,/)') ' END OF UNIT: ',iunit
   WRITE (0,'(A,I3)')     ' IOSTAT == ',iost
   EXIT reports
ENDIF

!     READ EACH LEVELS
!     ----------------

loop_level: DO ii = 1, levels

   READ (UNIT = iunit, FMT = '( 3X, F12.1, 2(F12.3,I4,F12.3,2X) )' ) &
        height,         &
        rv_inv,         &
        rv_qc,          &
        rv_error,       &
        rf_inv,         &
        rf_qc,          &
        rf_error

   if (rv_inv /= missing_r8) then

      num_obs = num_obs + 1
      key = key + 1

      location = set_location(lon, lat, height, 3)
      call set_obs_def_location(obs_def, location)

      location = set_location(radar_long, radar_lat, radar_elev, 3)
      call set_obs_def_key(obs_def, key)

      obs_lat = lat*DEG2RAD
      rad_lat = radar_lat*DEG2RAD
      obs_lon = lon*DEG2RAD
      rad_lon = radar_long*DEG2RAD

      ae = 1000.0_r8 * earth_radius
      x = ae * cos((obs_lat + rad_lat)/2.0_r8) * (obs_lon - rad_lon)
      y = ae * (obs_lat - rad_lat)
      raz = atan(x/y)
      spath = sqrt(x*x + y*y)
      h = height - radar_elev

      ae = 4000.0_r8 * earth_radius / 3.0_r8
      rgate = spath

      do it=1,10
         elev_rad = asin((h*h + 2.0_r8*ae*h - rgate*rgate)/(2.0_r8*ae*rgate))
         rgate = (ae+h)*sin(spath/ae)/cos(elev_rad)
      enddo

      elev_obs = sqrt(1.5_r8 * h / (1000.0_r8 * earth_radius) + elev_rad*elev_rad)

      orientation(1) = sin(raz)*cos(elev_obs)
      orientation(2) = cos(raz)*cos(elev_obs)
      orientation(3) = sin(elev_obs)

      call set_rad_vel(key, location, orientation)

      call set_obs_def_time(obs_def, time)

      call set_obs_def_kind(obs_def, DOPPLER_RADIAL_VELOCITY)

      if (rv_error /= missing_r8) then
         call set_obs_def_error_variance(obs_def, rv_error*rv_error)
      else
         call set_obs_def_error_variance(obs_def, 2.0_r8)
      endif

      call set_obs_def(obs, obs_def)

      obs_value(1) = rv_inv
      call set_obs_values(obs, obs_value, 1)

      if (rv_inv == missing_r8 .or. &
           rv_error == missing_r8 ) then

         rv_qc = missing_data

      end if

      rstatus(1,1) = rv_qc
      call set_qc(obs, rstatus(1,:), 1)

      if(num_obs == 1) then
         call insert_obs_in_seq(seq, obs)
      else
         call insert_obs_in_seq(seq, obs, prev_obs)
      endif

      prev_obs = obs

   endif

   if (rf_inv /= missing_r8) then

      num_obs = num_obs + 1

      call set_obs_def_kind(obs_def, RADAR_REFLECTIVITY)

      if (rv_error /= missing_r8) then
         call set_obs_def_error_variance(obs_def, rf_error*rf_error)
      else
         call set_obs_def_error_variance(obs_def, 2.0_r8)
      endif

      call set_obs_def(obs, obs_def)

      obs_value(1) = rf_inv
      call set_obs_values(obs, obs_value, 1)

      if (rf_inv == missing_r8 .or. &
           rf_error == missing_r8 ) then

         rf_qc = missing_data

      end if
 
      rstatus(1,1) = rf_qc
      call set_qc(obs, rstatus(1,:), 1)

      if(num_obs == 1) then
         call insert_obs_in_seq(seq, obs)
      else
         call insert_obs_in_seq(seq, obs, prev_obs)
      endif

      prev_obs = obs

   endif

ENDDO loop_level

num_Radar = num_Radar + 1
   
ENDDO reports

call close_file(iunit)                                                        

!  PRINT OUT
!  =============
 
write(unit=*, fmt='(5x,a,i6,a)') &
     'Read:  ', num_Radar, ' Radar reports,'

! Write out the sequence
call write_obs_seq(seq, obs_seq_out_file_name)

write(logfileunit,*)'FINISHED wrf_3dvar_tf_dart.'
write(logfileunit,*)

call finalize_utilities ! closes the log file.
 
END PROGRAM wrf_3dvar_tf_dart
