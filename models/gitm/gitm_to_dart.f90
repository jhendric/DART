! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program gitm_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between gitm and DART
!
! method: Read gitm "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The gitm dirname is read from the gitm_in namelist
!         <edit gitm_to_dart_output_file in input.nml:gitm_to_dart_nml>
!         gitm_to_dart
!
! author: Tim Hoar 6/24/09
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, timestamp, &
                             find_namelist_in_file, check_namelist_read
use        model_mod, only : get_model_size, restart_file_to_sv, &
                             get_gitm_restart_dirname
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

character(len=128) :: gitm_to_dart_output_file  = 'dart.ud'

namelist /gitm_to_dart_nml/ gitm_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

logical               :: verbose = .TRUE.
integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)
character(len=256)    :: gitm_restart_dirname

!======================================================================

call initialize_utilities(progname='gitm_to_dart', output_flag=verbose)

!----------------------------------------------------------------------
! Read the namelist to get the output dirname.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "gitm_to_dart_nml", iunit)
read(iunit, nml = gitm_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "gitm_to_dart_nml") ! closes, too.

write(*,*)
write(*,'(''gitm_to_dart:converting gitm restart file '',A, &
      &'' to DART file '',A)') &
       trim(gitm_restart_dirname), trim(gitm_to_dart_output_file)

!----------------------------------------------------------------------
! get to work
!----------------------------------------------------------------------

x_size = get_model_size()
allocate(statevector(x_size))

call get_gitm_restart_dirname( gitm_restart_dirname )

call restart_file_to_sv(gitm_restart_dirname, statevector, model_time) 

iunit = open_restart_write(gitm_to_dart_output_file)

call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

!----------------------------------------------------------------------
! When called with 'end', timestamp will call finalize_utilities()
!----------------------------------------------------------------------

call print_date(model_time, str='gitm_to_dart:gitm  model date')
call print_time(model_time, str='gitm_to_dart:DART model time')
call timestamp(string1=source, pos='end')

end program gitm_to_dart

