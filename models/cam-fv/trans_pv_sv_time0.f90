! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
program trans_pv_sv_time0

! <next five lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
! $Name$

!----------------------------------------------------------------------
! purpose: interface between CAM and DART
!
! method: Read CAM 'initial' file (netCDF format).
!         Reform fields into a state vector.
!         Write out state vector in "proprietary" format for DART
!
! author: Kevin Raeder 2/21/03
!         based on prog_var_to_vector and vector_to_prog_var by Jeff Anderson
!
!----------------------------------------------------------------------

use        types_mod, only : r8
use    utilities_mod, only : logfileunit, error_handler, E_ERR, E_MSG, find_namelist_in_file, &
                             check_namelist_read
use        model_mod, only : model_type, init_model_instance, read_cam_init, &
                             prog_var_to_vector

use  assim_model_mod, only : assim_model_type, static_init_assim_model, &
   init_assim_model, get_model_size , set_model_state_vector, write_state_restart, &
   set_model_time, open_restart_read, open_restart_write, close_restart, &
   aread_state_restart
! Guam; move time stripping from advance_model to here
use time_manager_mod, only : time_type, read_time, set_time

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

character (len = 128) :: file_name = 'caminput.nc', file_out = 'temp_ud'
! Hawaii;                          file_time = 'temp_ic'
! trans_pv_sv_time0 should get it's time from the namelist, 
! not from temp_ic, which came from filter_ics, which will not exist for a new
! set of fields comprising the state vector.

! Temporary allocatable storage to read in a native format for cam state
type(assim_model_type) :: x
type(model_type)       :: var
type(time_type)        :: model_time
real(r8), allocatable  :: x_state(:), x_temp(:)
integer                :: file_unit, x_size, iunit, io

!-------------------------
! define junk for reading in whole namelist.
logical :: start_from_restart = .false., output_restart = .false.
integer :: async = 0
integer :: init_time_days = 0, init_time_seconds = 0, output_interval = 1
character(len = 129) :: restart_in_file_name  = 'perfect_ics',     &
                        restart_out_file_name = 'perfect_restart', &
                        obs_seq_in_file_name  = 'obs_seq.in',      &
                        obs_seq_out_file_name = 'obs_seq.out',     &
                        adv_ens_command       = './advance_ens.csh'


! namelist /filter_nml/init_time_days, init_time_seconds
namelist /perfect_model_obs_nml/ async, adv_ens_command, obs_seq_in_file_name, &
   obs_seq_out_file_name, start_from_restart, output_restart, &
   restart_in_file_name, restart_out_file_name, init_time_days, init_time_seconds, &
   output_interval

!-------------------------

! Static init assim model calls static_init_model
PRINT*,'static_init_assim_model in trans_pv_sv'
call static_init_assim_model()

! Initialize the assim_model instance
call init_assim_model(x)

! Allocate the local state vector
x_size = get_model_size()
allocate(x_state(x_size), x_temp(x_size))

! Allocate the instance of the cam model type for storage
call init_model_instance(var)

! Read the file cam state fragments into var
call read_cam_init(file_name, var)

! transform fields into state vector for DART
call prog_var_to_vector(var, x_state)

! Put this in the structure
call set_model_state_vector(x, x_state)

! Read the namelist entry
call find_namelist_in_file("input.nml", "perfect_model_obs_nml", iunit)
read(iunit, nml = perfect_model_obs_nml, iostat = io)
call check_namelist_read(iunit, io, "perfect_model_obs_nml")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'trans_pv_sv_time0','perfect_model_obs_nml values are',' ',' ',' ')
write(logfileunit, nml=perfect_model_obs_nml)
write(     *     , nml=perfect_model_obs_nml)


call filter_set_initial_time

call set_model_time (x, model_time)
call close_restart(file_unit)

! Get channel for output 
! debug file_unit = 13
file_unit = open_restart_write(file_out)
PRINT*,'In trans_pv_sv file_out unit = ',file_unit
PRINT*,' '
! write out state vector in "proprietary" format
call write_state_restart(x, file_unit)
call close_restart(file_unit)

contains

! WARNING: THERE IS SOME DANGER IN USING THESE SCOPED SUBROUTINES
!==========================================================================

!-------------------------------------------------------------------------
subroutine filter_set_initial_time

if(init_time_days >= 0) then
   model_time = set_time(init_time_seconds, init_time_days)
else
   model_time = set_time(0, 0)
endif

end subroutine filter_set_initial_time

!-------------------------------------------------------------------------
end program trans_pv_sv_time0

