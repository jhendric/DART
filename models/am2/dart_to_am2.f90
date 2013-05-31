! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
 
program dart_to_am2

!----------------------------------------------------------------------
! purpose: interface between AM2 and DART
!
! method: Read DART state vector ("proprietary" format), but not time(s).
!         Reform state vector back into AM2 fields.
!         Replace those fields on the AM2 initial file with the new values,
!         preserving all other information on the file.
!
! author: Patrick Hofmann 3/7/08
!         based on prog_var_to_vector and vector_to_prog_var by Robert Pincus
! mod:    to read temp_ic (assim_model_state_ic; 2 times) or temp_ud (1 time) and put
!         the fields into the AM2 initial files
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : get_unit, file_exist, open_file, logfileunit, &
                             initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read, close_file
use        model_mod, only : model_type, static_init_model, init_model_instance, &
                             write_model_init, vector_to_prog_var 
use  assim_model_mod, only : get_model_size, aread_state_restart, &
                             open_restart_read, close_restart
use time_manager_mod, only : time_type, read_time, get_time, &
                             print_time, print_date, operator(-)

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type(model_type)       :: var
type(time_type)        :: adv_to_time, model_time, deltat
real(r8), allocatable  :: statevector(:)
integer                :: x_size, sec, day, iunit, io

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 128) :: dart_to_am2_input_file = 'dart_restart'
character (len = 128) :: restart_file           = 'fv_rst.res.nc'
character (len = 128) :: tracer_file            = 'atmos_tracers.res.nc'
logical               :: advance_time_present   = .false.

namelist /dart_to_am2_nml/ dart_to_am2_input_file, &
                           restart_file, tracer_file, &
                           advance_time_present

!------------------------------------------------------------------
! Namelist coupler_nml default values
!------------------------------------------------------------------

integer               :: months = 0, days = 0, hours = 0
integer               :: dt_atmos = 1800, dt_ocean = 21600, dt_cpld = 21600
integer, dimension(6) :: current_date = (/2007,1,1,0,0,0/) 
character (len = 8)   :: calendar = "'julian'"
!character (len = 15)  :: current_date = '2006,1,1,0,0,0,'

namelist /coupler_nml/ &
   months, days, hours, dt_atmos, dt_ocean, dt_cpld, current_date, calendar

!------------------------------------------------------------------

call initialize_utilities('dart_to_am2')

call find_namelist_in_file("input.nml", "dart_to_am2_nml", iunit)
read(iunit, nml = dart_to_am2_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_am2_nml")

call static_init_model()
x_size = get_model_size()
allocate(statevector(x_size))

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_to_am2_input_file)
if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! Update the current AM2 restart file
!----------------------------------------------------------------------

! Allocate the instance of the am2 model type
call init_model_instance(var)

! decompose vector back into AM2 fields
call vector_to_prog_var (statevector, var)
deallocate (statevector)

! if any of the tracer fields are negative, change them to zero
! also, restrain CF to be between 0 and 1
where(var%tracers < 0) var%tracers = 0
where(var%tracers(:,:,:,3) > 1) var%tracers(:,:,:,3) = 1 

! write fields to the netCDF initial file
call write_model_init(restart_file, tracer_file, var)

!----------------------------------------------------------------------
! Write a new coupler namelist with advance-to-time if need be.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_am2: AM2  model date')
call print_time( model_time,'dart_to_am2: DART model time')
call print_date( model_time,'dart_to_am2: AM2  model date',logfileunit)
call print_time( model_time,'dart_to_am2: DART model time',logfileunit)

if (advance_time_present) then

   deltat = adv_to_time - model_time
   call get_time(deltat,sec,day)

   ! Read append.nml's coupler_nml values
   call find_namelist_in_file("append.nml","coupler_nml",iunit)
   read(iunit, nml = coupler_nml, iostat = io)
   call check_namelist_read(iunit, io, "coupler_nml")

   ! Change days and hours to advance
   days = day
   hours = sec/3600

   !Write newappend.nml with new days and hours variables
   iunit = get_unit()
   open(unit=iunit,file="newappend.nml",action='write')
   write(iunit,nml = coupler_nml, iostat = io)
   call close_file(iunit)

   call print_time(adv_to_time,'dart_to_am2:advance_to time')
   call print_date(adv_to_time,'dart_to_am2:advance_to date')
   call print_time(adv_to_time,'dart_to_am2:advance_to time',logfileunit)
   call print_date(adv_to_time,'dart_to_am2:advance_to date',logfileunit)

endif

call finalize_utilities('dart_to_am2')

end program dart_to_am2

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

