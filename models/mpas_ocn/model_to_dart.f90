! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program model_to_dart

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between model and DART
!
! method: Read MPAS "history" files of model state.
!         Reform fields into a DART state vector (control vector).
!         Write out state vector in "proprietary" format for DART.
!         The output is a "DART restart file" format.
! 
! USAGE:  The model filename is read from the model_in namelist
!         <edit model_to_dart_output_file in input.nml:model_to_dart_nml>
!         model_to_dart
!
! author: Tim Hoar 12 Sep 2011
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : initialize_utilities, finalize_utilities, &
                             find_namelist_in_file, check_namelist_read
use        model_mod, only : get_model_size, analysis_file_to_statevector, &
                             get_model_analysis_filename, static_init_model
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

character(len=128) :: model_to_dart_output_file  = 'dart.ud'

namelist /model_to_dart_nml/    &
     model_to_dart_output_file

!----------------------------------------------------------------------
! global storage
!----------------------------------------------------------------------

logical               :: verbose = .TRUE.
integer               :: io, iunit, x_size
type(time_type)       :: model_time
real(r8), allocatable :: statevector(:)
character(len=256)    :: model_analysis_filename

!======================================================================

call initialize_utilities(progname='model_to_dart')
! if verbose is false, E_MSG won't get printed.  i'm not sure we want that.
!call initialize_utilities(progname='model_to_dart', output_flag=verbose)

!----------------------------------------------------------------------
! Read the namelist to get the output filename.
!----------------------------------------------------------------------

call find_namelist_in_file("input.nml", "model_to_dart_nml", iunit)
read(iunit, nml = model_to_dart_nml, iostat = io)
call check_namelist_read(iunit, io, "model_to_dart_nml") ! closes, too.

!----------------------------------------------------------------------
! Call model_mod:static_init_model() which reads the model namelists
! to set grid sizes, etc.
!----------------------------------------------------------------------
call static_init_model()

call get_model_analysis_filename(model_analysis_filename)

x_size = get_model_size()
allocate(statevector(x_size))

write(*,*)
write(*,*) 'model_to_dart: converting model analysis file ', &
           "'"//trim(model_analysis_filename)//"'" 
write(*,*) ' to DART file ', "'"//trim(model_to_dart_output_file)//"'"

!----------------------------------------------------------------------
! Read the valid time and the state from the MPAS netcdf file
!----------------------------------------------------------------------
call analysis_file_to_statevector(model_analysis_filename, statevector, model_time) 

!----------------------------------------------------------------------
! Write the valid time and the state to the dart restart file
!----------------------------------------------------------------------
iunit = open_restart_write(model_to_dart_output_file)

call awrite_state_restart(model_time, statevector, iunit)
call close_restart(iunit)

!----------------------------------------------------------------------
! finish up
!----------------------------------------------------------------------

call print_date(model_time, str='model_to_dart:model model date')
call print_time(model_time, str='model_to_dart:DART model time')
call finalize_utilities()

end program model_to_dart

