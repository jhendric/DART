! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program naaps_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
 
!----------------------------------------------------------------------
! purpose: interface between naaps and DART
!
! method: Read naaps "restart" files of model state
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The naaps filename is read from the naaps_in namelist
!         <edit naaps_to_dart_output_file in input.nml:naaps_to_dart_nml>
!         naaps_to_dart
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read
use        model_mod, only : get_model_size, analysis_file_to_statevector, &
                             get_naaps_restart_path, get_naaps_metadata,   &
                             get_naaps_dtg, get_naaps_ensemble_member,     &
                             static_init_model
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

character(len=128) :: naaps_to_dart_output_file = 'dart.ud' 

namelist /naaps_to_dart_nml/ naaps_to_dart_output_file 

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

LOGICAL               :: verbose = .TRUE.
INTEGER               :: io, iunit, x_size
INTEGER               :: nx, ny, nz, ns, member
TYPE(time_type)       :: model_time
REAL(r8), allocatable :: statevector(:)
CHARACTER(len=256)    :: naaps_restart_path 
CHARACTER(len=10)     :: dtg

!======================================================================

call initialize_utilities(progname='naaps_to_dart', output_flag=verbose)

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "naaps_to_dart_nml", iunit)
read(iunit, nml = naaps_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "naaps_to_dart_nml")

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the model namelist
! to set grid size, paths, etc.
!----------------------------------------------------------------------

call static_init_model()
call get_naaps_restart_path( naaps_restart_path )
call get_naaps_dtg( dtg )
call get_naaps_metadata( naaps_restart_path, dtg, model_time )
call get_naaps_ensemble_member( member )

write(*,*)
write(*,'(''naaps_to_dart:converting naaps restart dir '',A, &
      &'' member '',i4,'' to DART file '',A)') &
       trim(naaps_restart_path), member, trim(naaps_to_dart_output_file)

!----------------------------------------------------------------------
! get to work
!----------------------------------------------------------------------

x_size = get_model_size()
allocate(statevector(x_size))

call analysis_file_to_statevector( naaps_restart_path, statevector, member, model_time ) 
iunit = open_restart_write(naaps_to_dart_output_file)

call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

!----------------------------------------------------------------------
! When called with 'end', timestamp will call finalize_utilities()
!----------------------------------------------------------------------

call print_date(model_time, str='naaps_to_dart:naaps  model date')
call print_time(model_time, str='naaps_to_dart:DART model time')
call finalize_utilities()

end program naaps_to_dart
