! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
program trans_sv_pv

! <next three lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$

!----------------------------------------------------------------------
! purpose: interface between CAM and DART
!
! method: Read DART state vector ("proprietary" format)
!         Reform state vector back into CAM fields.
!         Replace those fields on the CAM initial file with the new values,
!         preserving all other information on the file.
!
! author: Kevin Raeder 2/21/03
!         based on prog_var_to_vector and vector_to_prog_var by Jeff Anderson
!
!----------------------------------------------------------------------

use       types_mod, only : r8
use   utilities_mod, only : get_unit
use       model_mod, only : model_type, init_model_instance, write_cam_init, &
   vector_to_prog_var
use assim_model_mod, only : assim_model_type, static_init_assim_model, &
   init_assim_model, get_model_size, get_model_state_vector, read_state_restart, &
   open_restart_read, close_restart
! Guam clean out advance_model
use time_manager_mod, only : time_type, read_time

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(assim_model_type) :: x
type(model_type)       :: var
! Guam clean out advance_model
type(time_type)        :: adv_to_time
real(r8), allocatable  :: x_state(:)
integer                :: file_unit, x_size
character (len = 128)  :: file_name = 'caminput.nc', file_in = 'temp_ic'
character (len = 16)   :: file_form

! Static init assim model calls static_init_model
PRINT*,'static_init_assim_model in trans_sv_pv'
call static_init_assim_model()
call init_assim_model(x)

! Allocate the instance of the cam model type for storage
call init_model_instance(var)

! Get file for DART vector input
! debug file_unit = 17
file_unit = open_restart_read(file_in)
PRINT*,'In trans_sv_pv file_in unit  = ',file_unit
PRINT*,' '

! Guam clean out advance_model
! Read in time to which CAM must advance.  
! Neither this, nor time in x is used in this program
! read in state vector from DART
call read_state_restart(x, file_unit, adv_to_time)
call close_restart(file_unit)

! Get the state part of the assim_model type x
x_size = get_model_size()
allocate(x_state(x_size))
x_state = get_model_state_vector(x)

! decompose vector back into CAM fields
call vector_to_prog_var (x_state, var)
deallocate (x_state)

! write fields to the netCDF initial file
call write_cam_init(file_name, var)

end program trans_sv_pv
