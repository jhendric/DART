! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program dart_to_noah1D

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between DART and the NOAH model
!
! method: Read DART state vector and overwrite values in a noah restart file.
!         If the DART state vector has an 'advance_to_time' present, 
!         it is read ... but nothing happens with it at this time.
!         DART is NEVER expected to advance noah.
!
!         The dart_to_noah_nml namelist setting for advance_time_present 
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!
! author: Tim Hoar 12 July 2011
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, timestamp, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file, &
                             error_handler, E_MSG
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, print_time, print_date, operator(-), get_time
use        model_mod, only : static_init_model, dart_vector_to_model_file, &
                             get_model_size

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 128) :: dart_to_noah_input_file = 'dart_restart'
logical               :: advance_time_present   = .false.

namelist /dart_to_noah_nml/ dart_to_noah_input_file, &
                           advance_time_present

!----------------------------------------------------------------------

character(len=20)     :: noah_restart_filename = 'noah_input.nml'
integer               :: iunit, io, x_size
type(time_type)       :: model_time, adv_to_time
real(r8), allocatable :: statevector(:)
logical               :: verbose              = .FALSE.

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_noah', output_flag=verbose)

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the NOAH namelist
! to set location and state vector
!----------------------------------------------------------------------

call static_init_model()

x_size = get_model_size()
allocate(statevector(x_size))

! Read the namelist to get the input filename. 

call find_namelist_in_file("input.nml", "dart_to_noah_nml", iunit)
read(iunit, nml = dart_to_noah_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_noah_nml")

write(*,*)
write(*,'(''dart_to_noah:converting DART file '',A, &
      &'' to NOAH input namelist '',A)') &
     trim(dart_to_noah_input_file), trim(noah_restart_filename)

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_to_noah_input_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! write out the new namelist ...
!----------------------------------------------------------------------

if ( advance_time_present ) then
   call dart_vector_to_model_file(statevector, noah_restart_filename, model_time, adv_to_time)
else
   call dart_vector_to_model_file(statevector, noah_restart_filename, model_time)
endif

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_noah:noah  model date')
call print_time( model_time,'dart_to_noah:DART    model time')
call print_date( model_time,'dart_to_noah:noah  model date',logfileunit)
call print_time( model_time,'dart_to_noah:DART    model time',logfileunit)

if ( advance_time_present ) then
   call print_time(adv_to_time,'dart_to_noah:advance_to time')
   call print_date(adv_to_time,'dart_to_noah:advance_to date')
   call print_time(adv_to_time,'dart_to_noah:advance_to time',logfileunit)
   call print_date(adv_to_time,'dart_to_noah:advance_to date',logfileunit)
endif

! When called with 'end', timestamp will call finalize_utilities()
call timestamp(string1=source, pos='end')

end program dart_to_noah1D
