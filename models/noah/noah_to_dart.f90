! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program noah_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between NOAH and DART
!
! method: Read noah "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The noah filename is read from the noah_in namelist
!         <edit noah_to_dart_output_file in input.nml:noah_to_dart>
!         noah_to_dart
!
! author: Tim Hoar 11 July 2012
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, timestamp, &
                             find_namelist_in_file, check_namelist_read, &
                             logfileunit, do_output
use        model_mod, only : get_model_size, noah_to_dart_vector, &
                             get_noah_restart_filename
use  assim_model_mod, only : awrite_state_restart, open_restart_write, close_restart
use time_manager_mod, only : time_type, print_time, print_date

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!-----------------------------------------------------------------------
! namelist parameters with default values.
!-----------------------------------------------------------------------

character(len=128) :: noah_to_dart_output_file  = 'dart_ics'

namelist /noah_to_dart_nml/ noah_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

logical               :: verbose = .TRUE.
integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)
character(len=256)    :: noah_restart_filename

!======================================================================

call initialize_utilities(progname='noah_to_dart', output_flag=verbose)

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "noah_to_dart_nml", iunit)
read(iunit, nml = noah_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "noah_to_dart_nml") ! closes, too.

call get_noah_restart_filename( noah_restart_filename )

if (do_output()) then
   write(*,*)
   write(*,'(''noah_to_dart:converting noah restart file '',A, &
      &'' to DART file '',A)') &
       trim(noah_restart_filename), trim(noah_to_dart_output_file)
   write(logfileunit,*)
   write(logfileunit,'(''noah_to_dart:converting noah restart file '',A, &
      &'' to DART file '',A)') &
       trim(noah_restart_filename), trim(noah_to_dart_output_file)
endif

!----------------------------------------------------------------------
! get to work
!----------------------------------------------------------------------

x_size = get_model_size()
allocate(statevector(x_size))

call noah_to_dart_vector(noah_restart_filename, statevector, model_time) 

iunit = open_restart_write(noah_to_dart_output_file)

call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

!----------------------------------------------------------------------
! When called with 'end', timestamp will call finalize_utilities()
!----------------------------------------------------------------------

if (do_output()) then
   call print_date(model_time, str='noah_to_dart:DART model date',iunit=logfileunit)
   call print_date(model_time, str='noah_to_dart:DART model date')
   call print_time(model_time, str='noah_to_dart:DART model time')
   call print_time(model_time, str='noah_to_dart:DART model time',iunit=logfileunit)
endif

call timestamp(string1=source, pos='end')

end program noah_to_dart

