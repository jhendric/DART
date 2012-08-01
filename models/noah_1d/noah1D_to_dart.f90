! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program noah1D_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between noah1D and DART
!
! method: Read noah1D "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The noah1D filename is read from the noah1D_in namelist
!         <edit noah1D_to_dart_output_file in input.nml:noah1D_to_dart>
!         noah1D_to_dart
!
! author: Tim Hoar 11 July 2012
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, timestamp, &
                             find_namelist_in_file, check_namelist_read
use        model_mod, only : get_model_size, noah1d_to_dart_vector, &
                             get_noah1D_restart_filename
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

character(len=128) :: noah1D_to_dart_output_file  = 'dart_ics'

namelist /noah1D_to_dart_nml/ noah1D_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

logical               :: verbose = .TRUE.
integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)
character(len=256)    :: noah1D_restart_filename

!======================================================================

call initialize_utilities(progname='noah1D_to_dart', output_flag=verbose)

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "noah1D_to_dart_nml", iunit)
read(iunit, nml = noah1D_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "noah1D_to_dart_nml") ! closes, too.

call get_noah1D_restart_filename( noah1D_restart_filename )

write(*,*)
write(*,'(''noah1D_to_dart:converting noah1D restart file '',A, &
      &'' to DART file '',A)') &
       trim(noah1D_restart_filename), trim(noah1D_to_dart_output_file)

!----------------------------------------------------------------------
! get to work
!----------------------------------------------------------------------

x_size = get_model_size()
allocate(statevector(x_size))

call noah1d_to_dart_vector(noah1D_restart_filename, statevector, model_time) 

iunit = open_restart_write(noah1D_to_dart_output_file)

call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

!----------------------------------------------------------------------
! When called with 'end', timestamp will call finalize_utilities()
!----------------------------------------------------------------------

call print_date(model_time, str='noah1D_to_dart:noah1D  model date')
call print_time(model_time, str='noah1D_to_dart:DART    model time')
call timestamp(string1=source, pos='end')

end program noah1D_to_dart

