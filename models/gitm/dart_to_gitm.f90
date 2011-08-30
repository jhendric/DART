! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program dart_to_gitm

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between DART and the gitm model
!
! method: Read DART state vector and overwrite values in a gitm restart file.
!         If the DART state vector has an 'advance_to_time' present, a
!         file called gitm_in.DART is created with a time_manager_nml namelist 
!         appropriate to advance gitm to the requested time.
!
!         The dart_to_gitm_nml namelist setting for advance_time_present 
!         determines whether or not the input file has an 'advance_to_time'.
!         Typically, only temporary files like 'assim_model_state_ic' have
!         an 'advance_to_time'.
!
! author: Tim Hoar 25 Jun 09, revised 12 July 2010
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, timestamp, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, open_file, close_file
use  assim_model_mod, only : open_restart_read, aread_state_restart, close_restart
use time_manager_mod, only : time_type, print_time, print_date, operator(-), get_time
use        model_mod, only : static_init_model, sv_to_restart_file, &
                             get_model_size, get_base_time, get_gitm_restart_dirname

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!------------------------------------------------------------------
! The namelist variables
!------------------------------------------------------------------

character (len = 128) :: dart_to_gitm_input_file = 'dart.ic'
logical               :: advance_time_present       = .false.

namelist /dart_to_gitm_nml/ dart_to_gitm_input_file, &
                           advance_time_present

!----------------------------------------------------------------------

character(len=256)    :: gitm_restart_dirname
integer               :: iunit, io, x_size, diff1, diff2
type(time_type)       :: model_time, adv_to_time, base_time
real(r8), allocatable :: statevector(:)
logical               :: verbose              = .FALSE.

!----------------------------------------------------------------------

call initialize_utilities(progname='dart_to_gitm', output_flag=verbose)

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the gitm namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------

call static_init_model()

x_size = get_model_size()
allocate(statevector(x_size))

! Read the namelist to get the input dirname. 

call find_namelist_in_file("input.nml", "dart_to_gitm_nml", iunit)
read(iunit, nml = dart_to_gitm_nml, iostat = io)
call check_namelist_read(iunit, io, "dart_to_gitm_nml")

call get_gitm_restart_dirname( gitm_restart_dirname )

write(*,*)
write(*,'(''dart_to_gitm:converting DART file '',A, &
      &'' to gitm restart file '',A)') &
     trim(dart_to_gitm_input_file), trim(gitm_restart_dirname)

!----------------------------------------------------------------------
! Reads the valid time, the state, and the target time.
!----------------------------------------------------------------------

iunit = open_restart_read(dart_to_gitm_input_file)

if ( advance_time_present ) then
   call aread_state_restart(model_time, statevector, iunit, adv_to_time)
else
   call aread_state_restart(model_time, statevector, iunit)
endif
call close_restart(iunit)

!----------------------------------------------------------------------
! update the current gitm state vector
! Convey the amount of time to integrate the model ...
! time_manager_nml: stop_option, stop_count increments
!----------------------------------------------------------------------

call sv_to_restart_file(statevector, gitm_restart_dirname, model_time)

if ( advance_time_present ) then
   base_time = get_base_time(gitm_restart_dirname)
   call get_time((model_time  - base_time), diff1)
   call get_time((adv_to_time - base_time), diff2)
   iunit = open_file('times', action='write')
   write(iunit, '(I8, I8)') diff1, diff2
   call close_file(iunit)
endif

!----------------------------------------------------------------------
! Log what we think we're doing, and exit.
!----------------------------------------------------------------------

call print_date( model_time,'dart_to_gitm:gitm  model date')
call print_time( model_time,'dart_to_gitm:DART model time')
call print_date( model_time,'dart_to_gitm:gitm  model date',logfileunit)
call print_time( model_time,'dart_to_gitm:DART model time',logfileunit)

if ( advance_time_present ) then
call print_time(adv_to_time,'dart_to_gitm:advance_to time')
call print_date(adv_to_time,'dart_to_gitm:advance_to date')
call print_time(adv_to_time,'dart_to_gitm:advance_to time',logfileunit)
call print_date(adv_to_time,'dart_to_gitm:advance_to date',logfileunit)
endif

! When called with 'end', timestamp will call finalize_utilities()
call timestamp(string1=source, pos='end')

end program dart_to_gitm
