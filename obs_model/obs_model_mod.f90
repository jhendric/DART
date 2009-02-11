! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use types_mod,            only : r8
use utilities_mod,        only : register_module, error_handler,     &
                                 E_ERR, E_MSG, E_WARN,               &
                                 get_unit, file_exist
use assim_model_mod,      only : aget_closest_state_time_to, get_model_time_step, &
                                 open_restart_write, open_restart_read,           &
                                 awrite_state_restart, close_restart, adv_1step,  &
                                 aread_state_restart
use obs_sequence_mod,     only : obs_sequence_type, obs_type, get_obs_from_key,      &
                                 get_obs_def, init_obs, destroy_obs, get_num_copies, &
                                 get_num_qc, get_first_obs, get_next_obs_from_key,   &
                                 get_obs_time_range
use obs_def_mod,          only : obs_def_type, get_obs_def_time
use time_manager_mod,     only : time_type, set_time, get_time, print_time,   &
                                 operator(/=), operator(>), operator(-),      &
                                 operator(/), operator(+), operator(<), operator(==), &
                                 operator(<=), operator(>=)
use ensemble_manager_mod, only : get_ensemble_time, ensemble_type
use mpi_utilities_mod,    only : my_task_id, task_sync, task_count, block_task, &
                                 sum_across_tasks

implicit none
private

public :: move_ahead, advance_state

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical :: module_initialized = .false.
logical :: print_time_details = .true.

! Module storage for writing error messages
character(len = 129) :: errstring

contains

!======================================================================

subroutine initialize_module

call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module

!-------------------------------------------------------------------------

subroutine move_ahead(ens_handle, ens_size, seq, last_key_used, &
   key_bounds, num_obs_in_set, async, adv_ens_command)

! Advances all ensemble copies of the state to the time appropriate for
! the next unused observation in the observation sequence. Note that there
! may be more than ens_size copies in the ens_handle storage. Copies other
! than the ensemble copies need not be advanced in time.

implicit none

type(ensemble_type),     intent(inout) :: ens_handle
integer,                 intent(in)    :: ens_size
type(obs_sequence_type), intent(in)    :: seq
integer,                 intent(in)    :: last_key_used, async
integer,                 intent(out)   :: key_bounds(2), num_obs_in_set
character(len = *),      intent(in)    :: adv_ens_command

type(time_type)    :: next_time, time2, start_time, end_time, delta_time, ens_time
type(obs_type)     :: observation
type(obs_def_type) :: obs_def
logical            :: is_this_last, is_there_one, out_of_range
integer            :: my_first_copy, leaving_early

! Initialize if needed
if(.not. module_initialized) then
   call initialize_module
   module_initialized = .true.
endif

! Prepare values for 'early' returns
num_obs_in_set  =   0
key_bounds(1:2) = -99
leaving_early   =   0

! If none of my copies are regular ensemble members no need to advance. 
! This is true either if I have no copies or if the copy #s I have are
! all larger than the ensemble count.  This does assume that ensemble 
! copies are in the first ens_size global copies of ensemble and that 
! global indices are monotonically increasing in each pes local copies
! (note: Violating the private boundaries on ens_handle by direct access)
! This was the original code:
!   if((ens_handle%my_num_copies < 1) .or. (ens_handle%my_copies(1) > ens_size)) then

if(ens_handle%my_num_copies >= 1) then
   my_first_copy = ens_handle%my_copies(1)
else
   my_first_copy = 0
endif


! if this task has no ensembles to advance, return directly from here.
! this routine has code to sync up with the other ensembles which do
! have an advance to do - it will not return until the others are done.
if ((ens_handle%my_num_copies < 1) .or. (my_first_copy > ens_size)) then
   call wait_if_needed(async)
   return
endif

! if you get here, this task has at least one ensemble to try to advance;
! it is possible we are at the end of the observations and there in fact
! is no need to advance.  if so, we need to tell any sleepers to give up.

! Initialize a temporary observation type to use
call init_obs(observation, get_num_copies(seq), get_num_qc(seq))

! Get the next observation in the sequence that hasn't already been assimilated
if(last_key_used > 0) then
   call get_next_obs_from_key(seq, last_key_used, observation, is_this_last)
   ! If the end of the observation sequence has been passed, return
   if(is_this_last) leaving_early = 1
else
   is_there_one = get_first_obs(seq, observation)
   if(.not. is_there_one) leaving_early = 1
endif

if (leaving_early > 0) then
   call wait_if_needed(async)
   ! need to destroy obs here before returning
   call destroy_obs(observation)
   return
endif

! Get the time of this observation
call get_obs_def(observation, obs_def)
next_time = get_obs_def_time(obs_def)

! Get the time of the ensemble, assume consistent across all ensemble copies
call get_ensemble_time(ens_handle, 1, ens_time)

! Compute the model time step and center a window around the closest time
delta_time = get_model_time_step()

! print out current window, if requested
if (print_time_details) then
   if(delta_time / 2 > ens_time) then
      start_time = set_time(0, 0)
   else
      start_time = ens_time - delta_time / 2 + set_time(1, 0)
   endif
   end_time = ens_time + delta_time / 2
   call timechat(ens_time,    'move_ahead', .true.,   'Model data timestamp currently:        ')
   call timechat(start_time,  'move_ahead', .false.,  'Current assimilation window starts at: ')
   call timechat(end_time,    'move_ahead', .false.,  'Current assimilation window ends at:   ')
   !call timechat(delta_time,  'move_ahead', .false., 'Width of assimilation window:          ')
endif

! now recompute for the next window, so the code below can remain unchanged.
! Figure out time to which to advance model 

! More control over time window use of observations would come in here
time2 = aget_closest_state_time_to(ens_time, next_time)

! WATCH OUT FOR USING BOUNDARY TIME OBS TWICE; add one second to bottom time
! ALSO, avoid having a negative time for the start (for low-order models in general)
if(delta_time / 2 > time2) then
   start_time = set_time(0, 0)
else
   start_time = time2 - delta_time / 2 + set_time(1, 0)
endif
end_time = time2 + delta_time / 2

! Output very brief current start and end time at the message level
if (.not. print_time_details) then
   call timechat(start_time,  'move_ahead', .true.,   'Next assimilation window starts at: ')
   call timechat(end_time,    'move_ahead', .false.,  'Next assimilation window ends at:   ')
endif

! If the next observation is not in the window, then have an error
if(next_time < start_time .or. next_time > end_time .or. print_time_details) then
   ! Is this test still really needed?  If you comment out that test, you can
   ! simply skip early obs, but for now, leave the code alone and print out a
   ! better set of error messages.

   if (next_time < start_time .or. next_time > end_time) then
      call error_handler(E_MSG, ' ', ' ')
      call error_handler(E_MSG, 'move_ahead', 'Inconsistent model state/observation times: ')
   endif


   if (time2 /= ens_time) then
      call timechat(next_time,   'move_ahead', .true.,  'Next available observation time:       ')
      call timechat(time2,       'move_ahead', .false., 'New data time should be:               ', &
         'Not within current window, model will be called to advance state.')
   else 
      if (next_time >= start_time .and. next_time <= end_time) then
         call timechat(next_time,   'move_ahead', .true.,  'Next available observation time:       ', &
            'Within current assimilation window, model does not need advance.')
      else 
         call timechat(next_time,   'move_ahead', .true.,  'Next available observation time:       ', &
            'Next obs outside current assimilation window.')
      endif
   endif
   call error_handler(E_MSG, ' ', ' ')

   if (next_time < start_time .or. next_time > end_time) then
      call timechat(start_time,  'move_ahead', .true.,  'New assimilation window starts at:     ')
      call timechat(end_time,    'move_ahead', .false., 'New assimilation window ends at:       ')
      if (next_time < start_time) then
         call error_handler(E_MSG, 'move_ahead', &
            'Next observation cannot be earlier than start of new time window')
         call error_handler(E_MSG, 'move_ahead', &
            'If this is the start of the obs_seq file, ')
         call error_handler(E_MSG, 'move_ahead', &
            'use filter namelist to set first obs or data init time.')
      else
         call error_handler(E_MSG, 'move_ahead', &
            'Next observation is later than end of new time window')
         call error_handler(E_MSG, 'move_ahead', &  
            'should not happen; code has miscomputed how far to advance, somehow')
      endif
      call error_handler(E_MSG, ' ', ' ')
      call error_handler(E_ERR, 'move_ahead', &
        'Inconsistent model state/observation times, cannot continue', source, revision, revdate)
   endif
endif

! WANT BETTER WINDOW CONTROL, TIME INTERPOLATION, TOO.
! Get all the observations that are in the observation window
call get_obs_time_range(seq, start_time, end_time, key_bounds, num_obs_in_set, &
   out_of_range, observation)

! ok, not really a time detail, but if turned off, the output is pretty much
! the same as the original.
if (print_time_details) then
   write (errstring, '(A,I8,A)') 'Next window contains up to ', num_obs_in_set, ' observations'
   call error_handler(E_MSG, 'move_ahead', errstring)
endif

! Advance all ensembles to the time for the assimilation
if(time2 /= ens_time) then
   call advance_state(ens_handle, ens_size, time2, async, adv_ens_command)

   ! Get the time of the ensemble, assume consistent across all ensemble copies
   call get_ensemble_time(ens_handle, 1, ens_time)
   
   if (ens_time /= time2 .or. print_time_details) then
      call timechat(start_time, 'move_ahead', .false., 'New assimilation window starts at:     ')
      call timechat(end_time,   'move_ahead', .false., 'New assimilation window ends at:       ')

      ! error out if model state did not advance to when requested.
      if (ens_time /= time2) then
         call timechat(next_time,  'move_ahead', .true.,  'Model advance time NOT what requested: ')
         call error_handler(E_ERR, 'move_ahead', 'Model advance complete but model time not correct')
      else
         call error_handler(E_MSG, 'move_ahead', 'Model advance complete, model time updated to requested time')
         call error_handler(E_MSG, ' ', ' ')
      endif
   endif

else
   call wait_if_needed(async)
endif

! Release the storage associated with the observation temp variable
call destroy_obs(observation)

end subroutine move_ahead


!------------------------------------------------------------------------------

subroutine advance_state(ens_handle, ens_size, target_time, async, adv_ens_command) 

! Advances the model extended state until time is equal to the target_time.

implicit none

! ens_size is the number of ensembles that are model state and need to be advanced
type(ensemble_type), intent(inout) :: ens_handle
type(time_type),     intent(in)    :: target_time
integer,             intent(in)    :: ens_size, async
character(len=*),    intent(in)    :: adv_ens_command

character(len = 129), dimension(ens_handle%num_copies) :: ic_file_name, ud_file_name 
character(len = 129)                                   :: control_file_name
! WARNING: ARE THESE LENGTHS OKAY?
character(len = 129)                                   :: system_string

type(time_type) :: time_step
integer         :: is1, is2, id1, id2, my_num_state_copies, global_ens_index
integer         :: i, control_unit, ic_file_unit, ud_file_unit


! Initialize if needed
if(.not. module_initialized) then
   call initialize_module
   module_initialized = .true.
endif

! Determine model time_step
time_step = get_model_time_step()

! Loop through each model state and advance; Count number of copies that are state ensembles
my_num_state_copies = 0
ENSEMBLE_MEMBERS: do i = 1, ens_handle%my_num_copies
   global_ens_index = ens_handle%my_copies(i)

   ! If global_ens_index is not one of the model state ensembles, don't need to advance
   ! Assumes that state ensemble copies have smallest global index and are stored 
   ! contiguously on local processor
   if(global_ens_index > ens_size) exit ENSEMBLE_MEMBERS

   ! Increment number of ensemble member copies I have
   my_num_state_copies = my_num_state_copies + 1

   ! If no need to advance return
   if(ens_handle%time(i) == target_time) then
      call wait_if_needed(async)
      return
   endif

   ! Check for time error
   if(ens_handle%time(i) > target_time) then
      call get_time(ens_handle%time(i), is1, id1)
      call get_time(target_time, is2, id2)
      write(errstring,*)'target time ',is2,id2,' is before model_time ',is1,id1
      call error_handler(E_ERR,'advance_state', errstring, source, revision, revdate)
   endif

   !------------- Block for subroutine callable adv_1step interface -----------
   if(async == 0) then

      do while(ens_handle%time(i) < target_time)
         call adv_1step(ens_handle%vars(:, i), ens_handle%time(i))
         ens_handle%time(i) = ens_handle%time(i) + time_step
      end do

   !-------------- End single subroutine callable adv_1step interface ---------
   else
   !-------------- Block for calling shell to advance model -------------------

      ! Can only handle up to 10000 ensemble members
      if(global_ens_index >= 10000) then
         write(errstring,*)'Trying to use ', ens_size,' model states -- too many.'
         call error_handler(E_WARN,'advance_state',errstring,source,revision,revdate)
         call error_handler(E_ERR,'advance_state','Use less than 10000 member ensemble.', &
            source,revision,revdate)
      endif

      ! Create file names for input and output state files for this copy
      write(ic_file_name(i), '("assim_model_state_ic.", i4.4)') global_ens_index
      write(ud_file_name(i), '("assim_model_state_ud.", i4.4)') global_ens_index


      ! Open a restart file and write state following a target time
      ! Force writing to binary to have bit-wise reproducing advances
      ic_file_unit = open_restart_write(trim(ic_file_name(i)), "unformatted")
      call awrite_state_restart(ens_handle%time(i), ens_handle%vars(:, i), ic_file_unit, target_time)
      call close_restart(ic_file_unit)

   endif

   !-------------- End of block for calling shell to advance model -------------

end do ENSEMBLE_MEMBERS


! Following is for async options that use shell to advance model
SHELL_ADVANCE_METHODS: if(async /= 0) then
   ! Get a unique name for the control file; use process id
   if(my_task_id() >= 10000) call error_handler(E_ERR, 'advance_state', &
      'Can only have 10000 processes', source, revision, revdate)
   write(control_file_name, '("filter_control", i5.5)') my_task_id()

   if (async == 2) then

      ! Everyone writes the ic and ud file names to the control file
      control_unit = get_unit()
      open(unit = control_unit, file = trim(control_file_name))
      ! Write out the file names to a control file
      ! Also write out the global ensemble member number for the script 
      ! if it needs to create unique filenames for its own use.
      do i = 1, my_num_state_copies
         write(control_unit, '(i5)') ens_handle%my_copies(i)
         write(control_unit, '(a)' ) trim(ic_file_name(i))
         write(control_unit, '(a)' ) trim(ud_file_name(i))
      end do
      close(control_unit)

      ! Arguments to advance model script are unique id and number of copies
      write(system_string, '(i10, 1x, i10)') my_task_id(), my_num_state_copies
     
      ! Issue a system command with arguments my_task, my_num_copies, and control_file
      call system(trim(adv_ens_command)//' '//trim(system_string)//' '//trim(control_file_name))
   
      ! if control file is still here, the advance failed
      if(file_exist(control_file_name)) then
        write(errstring,*)'control file for task ',my_task_id(),' still exists; model advance failed'
        call error_handler(E_ERR,'advance_state', errstring, source, revision, revdate)
      endif

   else if (async == 4) then

      ! Only task 0 writes all ens names to file.
      if (my_task_id() == 0) then
         ! Write the ic and ud file names to the control file
         control_unit = get_unit()
         open(unit = control_unit, file = trim(control_file_name))
         ! Write out the file names to a control file
         ! Also write out the global ensemble member number for the script 
         ! if it needs to create unique filenames for its own use.
         do i = 1, ens_size
            write(control_unit, '(i5)') i
            write(control_unit, '("assim_model_state_ic.", i4.4)') i
            write(control_unit, '("assim_model_state_ud.", i4.4)') i
         end do
         close(control_unit)
      endif

      ! make sure all tasks have finished writing their initial condition
      ! files before letting process 0 start the model advances.
      call task_sync()

      ! PE0 tells the run script to all the model advance script here, then sleeps.  
      ! the other tasks just sleep.  TODO: this code is too intertwined with the
      ! run script - it assumes that the script has already created 2 fifo/pipe 
      ! files with specific names:  mkfifo filter_to_model.lock model_to_filter.lock
      ! eventually it should probably be yet another separate executable which
      ! keeps the code all together in the mpi module.

      ! tell any waiting tasks that we are indeed doing a model advance
      ! here and they need to block.
      call sum_across_tasks(1, i)

      ! block, not spin, until the model advance is done.
      ! the code does not resume running until the separate program 'wakeup_filter'
      ! calls the companion 'resume_task()'.
      call block_task()


      ! if control file is still here, the advance failed
      if(file_exist(control_file_name)) then
         write(errstring,*)'control file for task ',my_task_id(),' still exists; model advance failed'
         call error_handler(E_ERR,'advance_state', errstring, source, revision, revdate)
      endif

   else
      ! Unsupported option for async error
      write(errstring,*)'input.nml - async is ',async,' must be 0, 2, or 4'
      call error_handler(E_ERR,'advance_state', errstring, source, revision, revdate)
   endif

   ! All should be done, read in the states and proceed
   do i = 1, my_num_state_copies
      ud_file_unit = open_restart_read(trim(ud_file_name(i)))
      call aread_state_restart(ens_handle%time(i), ens_handle%vars(:, i), ud_file_unit)
      call close_restart(ud_file_unit)
   end do

end if SHELL_ADVANCE_METHODS

end subroutine advance_state

!--------------------------------------------------------------------

subroutine wait_if_needed(async) 

! returns true if this task has at least one ensemble to advance.
! if not, this routine blocks until it is safe to return in sync
! with all other tasks with work to do.

integer,             intent(in)    :: async

integer :: should_block


! if async = 0, 2, or 5, no synchronizing needed.  return now.
if ((async == 0) .or. (async == 2) .or. (async == 5)) return


! If only task 0 is handling the advances for all tasks, we have to 
! be sure they are all done writing out the ic files, and then afterwards
! that they do not try to read the ud files before they are done.  
! Also, this task needs to block, not spin, while the model advances.

if (async == 4) then

   ! make sure all tasks have finished writing their initial condition
   ! files before letting process 0 start the model advances.
   call task_sync()

   ! collective call to see if the code is planning to block or not.
   ! do not block only on this tasks behalf.  if no one needs to block, 
   ! return here.
   call sum_across_tasks(0, should_block)
   if (should_block == 0) return

   ! sleep (not spin) until the model advance has finished.
   call block_task()
  
   return

endif

! should not reach here with valid values of async.
write(errstring,*) 'input.nml - async is ',async,' must be 0, 2, 4 or 5'
call error_handler(E_ERR, 'work_to_do', errstring, '', '', '')

end subroutine wait_if_needed

!--------------------------------------------------------------------

subroutine timechat(a_time, label, blank, string1, string2, string3)
 type(time_type),  intent(in)           :: a_time
 character(len=*), intent(in)           :: label, string1
 logical,          intent(in)           :: blank
 character(len=*), intent(in), optional :: string2, string3

 
integer :: a_day, a_sec

call get_time(a_time, a_sec, a_day)
write (errstring, '(A,A,I8,A,I6)') string1, ' day=', a_day, ' sec=', a_sec

if (blank) call error_handler(E_MSG, ' ', ' ')
call error_handler(E_MSG, label, trim(errstring))

if (present(string2)) then
   call error_handler(E_MSG, label, string2)
endif

if (present(string3)) then
   call error_handler(E_MSG, label, string3)
endif

end subroutine timechat
!--------------------------------------------------------------------

end module obs_model_mod
