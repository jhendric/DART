! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_sequence_mod

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$

! WARNING OPERATOR OVERLOAD FOR EQUIVALENCE???
! FURTHER WARNING: Compiler problems exist with the use of assignment(=) in
! use only statement. First, can only use it at the level above if the internals
! of the type are not private. Second, if I inherit assignment(=) from obs_def
! and also define one in obs_sequence, I get an error if I try to make it public
! to a module that uses obs_sequence but not obs_def with the intel compiler. No
! obvious workaround exists. For now, make modules at higher levels use explicit
! copy subruoutines. USERS MUST BE VERY CAREFUL TO NOT DO DEFAULT ASSIGNMENT
! FOR THESE TYPES THAT HAVE COPY SUBROUTINES.

use types_mod, only : r8
use obs_def_mod, only : obs_def_type, get_obs_def_time, read_obs_def, &
   write_obs_def, destroy_obs_def, interactive_obs_def, copy_obs_def
use time_manager_mod, only : time_type, operator(>), operator(<), operator(>=), &
   operator(/=)
use utilities_mod, only : get_unit, open_file, close_file, file_exist, check_nml_error, &
   register_module, error_handler, E_ERR, E_MSG

implicit none

private

interface assignment(=)
   module procedure copy_obs
end interface

! Public interfaces for obs sequences
public obs_sequence_type, init_obs_sequence, interactive_obs_sequence, &
   get_num_copies, get_num_qc, get_num_obs, get_max_num_obs, get_copy_meta_data, &
   get_qc_meta_data, get_next_obs, get_prev_obs, insert_obs_in_seq, &
   delete_obs_from_seq, set_copy_meta_data, set_qc_meta_data, get_first_obs, &
   get_last_obs, add_copies, add_qc, write_obs_seq, read_obs_seq,  &
   append_obs_to_seq, get_obs_from_key, get_obs_time_range, set_obs, get_time_range_keys, &
   get_num_times, static_init_obs_sequence

! Public interfaces for obs
public obs_type, init_obs, destroy_obs, get_obs_def, set_obs_def, &
   get_obs_values, set_obs_values, get_qc, set_qc, write_obs, read_obs, &
   interactive_obs, copy_obs, assignment(=)

! Public interfaces for obs covariance modeling
public obs_cov_type

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type obs_sequence_type
   private
   integer :: num_copies, num_qc, num_obs, max_num_obs
   character(len = 129), pointer :: copy_meta_data(:)
   character(len = 129), pointer :: qc_meta_data(:)
   integer :: first_time, last_time
!   integer :: first_avail_time, last_avail_time
   type(obs_type), pointer :: obs(:)
! What to do about groups
end type obs_sequence_type

type obs_type
   private
! The key is needed to indicate the element number in the storage for the obs_sequence
! Do I want to enforce the identity of the particular obs_sequence?
   integer :: key
   type(obs_def_type) :: def
   real(r8), pointer :: values(:)
   real(r8), pointer :: qc(:)
! Put sort indices directly into the data structure
   integer :: prev_time, next_time
   integer :: cov_group
end type obs_type

type obs_cov_type
   private    
   integer :: num_cov_groups
! ??????
end type obs_cov_type

!-------------------------------------------------------------
! Namelist with default values
! binary_restart_files  == .true.  -> use unformatted file format.
!                                     Full precision, faster, smaller,
!                                     but not as portable.
! binary_restart_files  == .false.  -> use ascii file format.
!                                     Portable, but loses precision,
!                                     slower, and larger.

logical  :: read_binary_obs_sequence = .false.
logical  :: write_binary_obs_sequence = .false.

namelist /obs_sequence_nml/ read_binary_obs_sequence, write_binary_obs_sequence

!-------------------------------------------------------------      


contains

!--------------------------------------------------------------

subroutine static_init_obs_sequence

! reads namelist and registers module
! Read the namelist input                                                                                             

integer :: iunit, ierr, io

call register_module(source, revision, revdate)

if(file_exist('input.nml')) then
   iunit = open_file(fname = 'input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(iunit, nml = obs_sequence_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'obs_sequence_nml')
   enddo
 11 continue
   call close_file(iunit)
endif             

write(*, *) 'in static_init_obs_sequence'
write(*, *) 'read write', read_binary_obs_sequence, write_binary_obs_sequence

end subroutine static_init_obs_sequence

!--------------------------------------------------------------

!WHAT ABOUT PASS THROUGHS TO THE OBS_DEF???
! WhAT ABOUT copy_obs_sequence similar to read.
!-------------------------------------------------
subroutine init_obs_sequence(seq, num_copies, num_qc, &
   expected_max_num_obs)

! Constructor for an obs_sequence

type(obs_sequence_type), intent(out) :: seq
integer, intent(in) :: num_copies, num_qc, expected_max_num_obs

integer :: i

seq%num_copies = num_copies
seq%num_qc = num_qc
seq%num_obs = 0
seq%max_num_obs = expected_max_num_obs
allocate(seq%copy_meta_data(num_copies), &
   seq%qc_meta_data(num_qc), seq%obs(expected_max_num_obs))
do i = 1, num_copies
   seq%copy_meta_data(i) = 'Copy metadata not initialized'
end do
do i = 1, num_qc
   seq%qc_meta_data(i) = 'QC metadata not initialized'
end do
! Initialize the pointers to allocated but pointing to zero space???
do i = 1, expected_max_num_obs
   allocate(seq%obs(i)%values(num_copies), seq%obs(i)%qc(num_qc))
end do
seq%first_time = -1
seq%last_time = -1
!seq%first_avail_time = -1
!seq%last_avail_time = -1

end subroutine init_obs_sequence

!--------------------------------------------------------------

function interactive_obs_sequence()

! Interactive creation of an observation sequence
type(obs_sequence_type) :: interactive_obs_sequence

type(obs_type) :: obs
integer :: max_num_obs, num_copies, num_qc, i, end_it_all


write(*, *) 'Input upper bound on number of observations in sequence'
read(*, *) max_num_obs

write(*, *) 'Input number of copies of data (0 for just a definition)'
read(*, *) num_copies

write(*, *) 'Input number of quality control values per field (0 or greater)'
read(*, *) num_qc

! Initialize an obs_sequence structure
call init_obs_sequence(interactive_obs_sequence, num_copies, num_qc, max_num_obs)

do i = 1, num_copies
   write(*, *) 'input meta data for data copy ', i
   read(*, *) interactive_obs_sequence%copy_meta_data(i)
end do

do i = 1, num_qc
   write(*, *) 'input meta data for qc field ', i
   read(*, *) interactive_obs_sequence%qc_meta_data(i)
end do

! Initialize the obs variable
call init_obs(obs, num_copies, num_qc)

! Loop to initialize each observation in turn; terminate by ???
do i = 1, max_num_obs
   write(*, *) 'input a -1 if there are no more obs'
   read(*, *) end_it_all
   if(end_it_all == -1) exit
   call interactive_obs(num_copies, num_qc, obs)
   if(i == 1) then
      write(*, *) 'calling insert obs in sequence'
      call insert_obs_in_seq(interactive_obs_sequence, obs)
      write(*, *) 'back from insert obs in sequence'
   else
      call insert_obs_in_seq(interactive_obs_sequence, obs, interactive_obs_sequence%obs(i - 1))
   endif
end do

end function interactive_obs_sequence

!---------------------------------------------------------

function get_num_copies(seq)


type(obs_sequence_type), intent(in) :: seq
integer :: get_num_copies

get_num_copies = seq%num_copies

end function get_num_copies

!-------------------------------------------------

function get_num_qc(seq)


type(obs_sequence_type), intent(in) :: seq
integer :: get_num_qc

get_num_qc= seq%num_qc

end function get_num_qc

!-------------------------------------------------

function get_num_obs(seq)


type(obs_sequence_type), intent(in) :: seq
integer :: get_num_obs

get_num_obs = seq%num_obs

end function get_num_obs

!-------------------------------------------------

function get_max_num_obs(seq)


type(obs_sequence_type), intent(in) :: seq
integer :: get_max_num_obs

get_max_num_obs = seq%max_num_obs

end function get_max_num_obs
!-------------------------------------------------

function get_copy_meta_data(seq, copy_num)


type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: copy_num
character(len=129) :: get_copy_meta_data

! Should have an error check for copy_num range
get_copy_meta_data = seq%copy_meta_data(copy_num)

end function get_copy_meta_data

!-------------------------------------------------
function get_qc_meta_data(seq, qc_num)


type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: qc_num
character(len=129) :: get_qc_meta_data

! Should have an error check for qc_num range
get_qc_meta_data = seq%qc_meta_data(qc_num)

end function get_qc_meta_data

!-------------------------------------------------

subroutine get_next_obs(seq, obs, next_obs, is_this_last)


type(obs_sequence_type), intent(in) :: seq
type(obs_type), intent(in) :: obs
type(obs_type), intent(out) :: next_obs
logical, intent(out) :: is_this_last

integer :: next_index

! Get index of the next observation
next_index = obs%next_time
if(next_index == -1) then
   is_this_last = .true.
   return
else
   is_this_last = .false.
   next_obs = seq%obs(next_index)
endif

end subroutine get_next_obs

!-------------------------------------------------

subroutine get_prev_obs(seq, obs, prev_obs, is_this_first)


type(obs_sequence_type), intent(in) :: seq
type(obs_type), intent(in) :: obs
type(obs_type), intent(out) :: prev_obs
logical, intent(out) :: is_this_first

integer :: prev_index

! Get index of the next observation
prev_index = obs%prev_time
if(prev_index == -1) then
   is_this_first= .true.
   return
else
   is_this_first= .false.
   prev_obs = seq%obs(prev_index)
endif

end subroutine get_prev_obs

!-------------------------------------------------------------

subroutine get_obs_from_key(seq, key, obs)

type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: key
type(obs_type) :: obs

obs = seq%obs(key)

end subroutine get_obs_from_key

!-----------------------------------------------------------------

subroutine set_obs(seq, obs, key_in)

! Copies the obs into the key element of sequence where key is the key field
! in obs. If the integer argument key is present, the obs is copied into
! the key-th element of the sequence.

type(obs_sequence_type), intent(inout) :: seq
type(obs_type), intent(in) :: obs
integer, intent(in), optional :: key_in

integer :: key

! Get the key to copy into
if(present(key_in)) then 
   key = key_in
else
   key = obs%key
endif

seq%obs(key) = obs

! Make sure the key in sequence is set propoerly
seq%obs(key)%key = key

end subroutine set_obs

!-------------------------------------------------------------------

subroutine get_obs_time_range(seq, time1, time2, key_bounds, num_keys, out_of_range, obs)

! Add other options for getting the first time to minimize searh
type(obs_sequence_type),  intent(in) :: seq
type(time_type),          intent(in) :: time1, time2
integer,                 intent(out) :: key_bounds(2)
integer,                 intent(out) :: num_keys
logical,                 intent(out) :: out_of_range
type(obs_type), intent(in), optional :: obs

type(time_type)    :: cur_time
type(obs_def_type) :: obs_def
integer            :: current, last_key

! Returns the first key and last key of sequence of obs between time1 and
! time2 along with the total number.
! A complete list of the keys can be obtained by call to get_time_range_keys
! Logical out_of_range is true if the time range is all past the end of sequence times

num_keys = 0
out_of_range = .false.

! The optional argument obs says the search can be started at this observation

! Figure out where to begin search
if(present(obs)) then
   current = obs%key
else
   current = seq%first_time
endif

! Find the first element in the time window
do while(current /= -1)
   call get_obs_def(seq%obs(current), obs_def)
   cur_time = get_obs_def_time(obs_def)
   if(cur_time >= time1) goto 10
   current = seq%obs(current)%next_time
end do
! Falling off the end means there are no times greater than time1
out_of_range = .true.
return

10 continue
! current is pointer to first

! First pass, count the keys for storage requirements
key_bounds(1) = current
last_key = current
do while(current /= -1)
   call get_obs_def(seq%obs(current), obs_def)
   cur_time = get_obs_def_time(obs_def)
   if(cur_time > time2) goto 20
! Found a time in the range
   num_keys = num_keys + 1
   last_key = current
   current = seq%obs(current)%next_time
end do

20 continue
key_bounds(2) = last_key

end subroutine get_obs_time_range

!---------------------------------------------------------------

subroutine get_time_range_keys(seq, key_bounds, num_keys, keys)

! Given bounds from get_obs_time_range and an array keys big enough to hold
! all the keys in the range, returns the keys in the range

type(obs_sequence_type), intent(in) :: seq
integer, intent(in) :: key_bounds(2), num_keys
integer, intent(out) :: keys(num_keys)

integer :: current, i

! Now loop through again to get these keys
current = key_bounds(1)
do i = 1, num_keys
   keys(i) = seq%obs(current)%key
   current = seq%obs(current)%next_time
end do

end subroutine get_time_range_keys


!-------------------------------------------------

subroutine insert_obs_in_seq(seq, obs, prev_obs)

type(obs_sequence_type), intent(inout) :: seq
type(obs_type), intent(inout) :: obs
type(obs_type), intent(in), optional :: prev_obs

type(time_type) :: obs_time, current_time
integer :: prev, next, current

! Inserts an observation into a sequence, optional argument
! prev_obs says that this was the predecessor in time.
! This avoids time search in cases where one is building
! a sequence from scratch.


! Make sure there is room, fail for now if not
if(seq%num_obs == seq%max_num_obs) then
! Later do an increase of space and copy
   write(*, *) 'ran out of room in obs_sequence'
   stop
endif

! Set the key for the observation
obs%key = seq%num_obs + 1
seq%num_obs = seq%num_obs + 1

! Get the time for the observation
obs_time = get_obs_def_time(obs%def)

if(present(prev_obs)) then
    prev = prev_obs%key
    next = prev_obs%next_time
else

! Have to search through the linked list to find last member
! already in with a time less than or equal to obs time
   prev = -1
   current = -1
   next = seq%first_time
   do while(next /= -1)
      prev = current
      current = next
      next = seq%obs(current)%next_time
      current_time = get_obs_def_time(seq%obs(current)%def)
! If the time of the observation in the sequence is >, stop
      if(current_time > obs_time) then 
! The observation that will follow the one being inserted is current
         next = current
         goto 10 
      endif
   end do

! Falling off the end means that next is -1, so current should be previous for insertion
   prev = current
endif

! If the time check occured, previous is already pointing to previous
10 continue

! prev now holds the key of the previous observation, next holds the one after

! Link into the foward moving pointer chain
! If prev is -1, new observation goes at the start
if(prev == -1) then
   obs%next_time = seq%first_time
   obs%prev_time = -1
   seq%first_time = obs%key
else
   obs%prev_time = prev
   obs%next_time = next
   seq%obs(prev)%next_time = obs%key
endif   

! Link into the backward moving pointer chain
if(next == -1) then
   obs%prev_time = seq%last_time
   obs%next_time = -1
   seq%last_time = obs%key
else
   seq%obs(next)%prev_time = obs%key
endif

! Finally, copy this obs structure into the sequence
seq%obs(obs%key) = obs

end subroutine insert_obs_in_seq

!----------------------------------------------------------------------

subroutine append_obs_to_seq(seq, obs)

! Appends an observation to an existing sequence; Error if new obs is 
! not later than time of last obs already in seq

type(obs_sequence_type), intent(inout) :: seq
type(obs_type), intent(inout) :: obs

type(obs_type) :: last_obs
type(time_type) :: obs_time, last_time

! Initialize obs_type before using
call init_obs(last_obs, 0, 0)

! If this is first, just put it in
if(.not. get_last_obs(seq, last_obs)) then
   call insert_obs_in_seq(seq, obs)
else

! Otherwise, get last obs from sequence and do insert with it as
! the previous after checking times

! Get the time for the observation
   obs_time = get_obs_def_time(obs%def)
   last_time = get_obs_def_time(last_obs%def)
   if(obs_time < last_time) then
      write(*, *) 'Error: tried to append an obs to sequence with bad time'
      stop
   endif

!!!   call insert_obs_in_seq(seq, obs)
!!!   if(1 == 1) return

! Make sure there is room, fail for now if not
   if(seq%num_obs == seq%max_num_obs) then
! Later do an increase of space and copy
      write(*, *) 'ran out of room in obs_sequence'
      stop
   endif

! Set the key for the observation
   obs%key = seq%num_obs + 1
   seq%num_obs = seq%num_obs + 1
! Link into the pointer chains
! Previous last points to this one, this one points back to previous last
   obs%prev_time = seq%last_time
   seq%obs(seq%last_time)%next_time = obs%key
   seq%last_time = obs%key

! Put this obs into the sequence's last slot
   seq%obs(seq%num_obs) = obs

endif

end subroutine append_obs_to_seq

!---------------------------------------------------------------

!subroutine insert_obs_group_in_seq(seq, obs_grp, prev_obs)

! Insert a group of observations from the same time into a sequence
!type(obs_sequence_type), intent(inout) :: seq
!type(obs_type), intent(inout) :: obs
!type(obs_type), intent(in), optional :: prev_obs
!
!end subroutine insert_obs_group_in_seq

!-------------------------------------------------

subroutine delete_obs_from_seq(seq, obs)

! Removes this observation from the sequence, does not free storage in this implementation
type(obs_sequence_type), intent(inout) :: seq
type(obs_type), intent(inout) :: obs

integer :: prev, next

prev = obs%prev_time
next = obs%next_time

! Previous should now point to next; if it deleted was first update sequence first_time
if(prev /= -1) then
   seq%obs(prev)%next_time = next
else
   seq%first_time = next
endif

! Next should point to previous
seq%obs(next)%prev_time = prev

! If deleted obs was first, update first pointer in sequence
if(prev == -1) seq%first_time = next

! If deleted was last, update last pointer in sequence
if(next == -1) seq%last_time = prev

end subroutine delete_obs_from_seq

!-------------------------------------------------
subroutine set_copy_meta_data(seq, copy_num, meta_data)

! Need all sorts of error checking to avoid silly stuff eventually

type(obs_sequence_type), intent(inout) :: seq
integer, intent(in) :: copy_num
character(len = 129), intent(in) :: meta_data

seq%copy_meta_data(copy_num) = meta_data

end subroutine set_copy_meta_data

!-------------------------------------------------

subroutine set_qc_meta_data(seq, qc_num, meta_data)

! Need error checks
type(obs_sequence_type), intent(inout) :: seq
integer, intent(in) :: qc_num
character(len = 129), intent(in) :: meta_data

seq%qc_meta_data(qc_num) = meta_data

end subroutine set_qc_meta_data

!-------------------------------------------------

function get_first_obs(seq, obs)

type(obs_sequence_type), intent(in) :: seq
type(obs_type),         intent(out) :: obs
logical                             :: get_first_obs

if(seq%num_obs == 0) then
   get_first_obs = .false.
else
   get_first_obs = .true.
   obs = seq%obs(seq%first_time)
endif

end function get_first_obs

!-------------------------------------------------

function get_last_obs(seq, obs)

type(obs_sequence_type), intent(in) :: seq
type(obs_type), intent(out) :: obs
logical :: get_last_obs

if(seq%num_obs == 0) then
   get_last_obs = .false.
   return
else
   get_last_obs = .true.
   obs = seq%obs(seq%last_time)
endif

end function get_last_obs

!-------------------------------------------------

subroutine add_copies(seq, num_to_add)

! This requires a complete recreation of the entire obs sequence
! Add additional copyies to an observation sequence. This increases
! the space for copy meta_data and goes through the whole string of
! observations deallocating and allocating (yuck), to add space.
! In the long run, may want a smoother way to do this globally.

type(obs_sequence_type), intent(inout) :: seq
integer, intent(in) :: num_to_add

character(len = 129) :: meta_temp(seq%num_copies)
real(r8) :: values_temp(seq%num_copies)
integer :: i, old_num

old_num = seq%num_copies
seq%num_copies = old_num + num_to_add

! Copy the old copy metadata to temp storage, reallocate and copy
if(old_num > 0) then
   meta_temp = seq%copy_meta_data
endif

! Deallocate and reallocate with enhanced length
deallocate(seq%copy_meta_data)
allocate(seq%copy_meta_data(old_num + num_to_add))
seq%copy_meta_data(1:old_num) = meta_temp
seq%copy_meta_data(old_num+1 : old_num + num_to_add) = 'Copy metadata not initialized'

! Loop through all the observations, copy and increase size
!??? WHAT IS THE STORY WITH NUM_OBS WHEN A DELETION IS DONE???
do i = 1, seq%max_num_obs

! Copy the existing values
   if(old_num > 0 .and. i < seq%num_obs) then
      values_temp = seq%obs(i)%values
   end if

! Deallocate, reallocate and copy
   deallocate(seq%obs(i)%values)
   allocate(seq%obs(i)%values(old_num + num_to_add))
   seq%obs(i)%values(1:old_num) = values_temp
   
end do

end subroutine add_copies

!-------------------------------------------------

subroutine add_qc(seq, num_to_add)

! This requires a complete recreation of the entire obs sequence
! Add additional copies to an observation sequence. This increases
! the space for copy meta_data and goes through the whole string of
! observations deallocating and allocating (yuck), to add space.
! In the long run, may want a smoother way to do this globally.

type(obs_sequence_type), intent(inout) :: seq
integer,                    intent(in) :: num_to_add

character(len = 129) :: qc_temp(seq%num_copies)
real(r8)             :: values_temp(seq%num_copies)
integer              :: i, old_num

old_num = seq%num_qc
seq%num_qc = old_num + num_to_add

! Copy the old copy metadata to temp storage, reallocate and copy
if(old_num > 0) then
   qc_temp = seq%qc_meta_data
endif

! Deallocate and reallocate with enhanced length
deallocate(seq%qc_meta_data)
allocate(seq%qc_meta_data(old_num + num_to_add))
seq%qc_meta_data(1:old_num) = qc_temp
seq%qc_meta_data(old_num+1 : old_num + num_to_add) = 'QC metadata not initialized'

! Loop through all the observations, copy and increase size
!??? WHAT IS THE STORY WITH NUM_OBS WHEN A DELETION IS DONE???
do i = 1, seq%max_num_obs

! Copy the existing values
   if(old_num > 0 .and. i < seq%num_obs) then
      values_temp = seq%obs(i)%qc
   end if

! Deallocate, reallocate and copy
   deallocate(seq%obs(i)%qc)
   allocate(seq%obs(i)%qc(old_num + num_to_add))
   seq%obs(i)%qc(1:old_num) = values_temp
   
end do

end subroutine add_qc

!------------------------------------------------------------------

subroutine write_obs_seq(seq, file_name)

type(obs_sequence_type), intent(in) :: seq
character(len = 129), intent(in) :: file_name

integer :: i, file_id

! Open the file
file_id = get_unit()
if(write_binary_obs_sequence) then
   write(*, *) 'opening file in write_obs-seq unformatted'
   open(unit = file_id, file = file_name, form = "unformatted")
else
   write(*, *) 'opening file in write_obs-seq regular'
   open(unit = file_id, file = file_name)
endif

! First inefficient ugly pass at writing an obs sequence, need to update for storage size
if(write_binary_obs_sequence) then
   write(file_id) seq%num_copies, seq%num_qc, seq%num_obs, seq%max_num_obs
else
   write(file_id, *) seq%num_copies, seq%num_qc, seq%num_obs, seq%max_num_obs
endif 

do i = 1, seq%num_copies
   if(write_binary_obs_sequence) then
      write(file_id) seq%copy_meta_data(i)
   else
      write(file_id, '(a129)') seq%copy_meta_data(i)
   endif
end do

do i = 1, seq%num_qc
   if(write_binary_obs_sequence) then
      write(file_id) seq%qc_meta_data(i)
   else
      write(file_id, '(a129)') seq%qc_meta_data(i)
   endif
end do

if(write_binary_obs_sequence) then
   write(file_id) seq%first_time, seq%last_time
else
   write(file_id, *) seq%first_time, seq%last_time
endif

! write(*, *) first_avail_time, last_avail_time
do i = 1, seq%num_obs
   call write_obs(seq%obs(i), file_id, seq%num_copies, seq%num_qc)
end do

! Close up the file
call close_file(file_id)

end subroutine write_obs_seq

!------------------------------------------------------------------

subroutine read_obs_seq(file_name, add_copies, add_qc, add_obs, seq)

! Be able to increase size at read in time for efficiency

character(len = 129), intent(in) :: file_name
integer, intent(in) :: add_copies, add_qc, add_obs
type(obs_sequence_type), intent(out) :: seq

integer :: i, num_copies, num_qc, num_obs, max_num_obs, file_id

! Open the file
file_id = get_unit()
if(read_binary_obs_sequence) then
   open(unit = file_id, file = file_name, form = "unformatted")
else
   open(unit = file_id, file = file_name)
endif

! First, determine the size that was written out
if(read_binary_obs_sequence) then
   read(file_id) num_copies, num_qc, num_obs, max_num_obs
else
   read(file_id, *) num_copies, num_qc, num_obs, max_num_obs
endif

call init_obs_sequence(seq, num_copies + add_copies, &
   num_qc + add_qc, num_obs + add_obs)

! Set the number of obs available at present
seq%num_obs = num_obs

! Get the available copy_meta_data
do i = 1, num_copies
   if(read_binary_obs_sequence) then
      read(file_id) seq%copy_meta_data(i)
   else
      read(file_id, '(a129)') seq%copy_meta_data(i)
   endif
end do

! Get the available qc_meta_data
do i = 1, num_qc
   if(read_binary_obs_sequence) then
      read(file_id) seq%qc_meta_data(i)
   else
      read(file_id, '(a129)') seq%qc_meta_data(i)
   endif
end do

! Read the first and last avail_time pointers
if(read_binary_obs_sequence) then
   read(file_id) seq%first_time, seq%last_time
   !read(file_id) seq%first_avail_time, seq%last_avail_time
else
   read(file_id, *) seq%first_time, seq%last_time
   !read(file_id, *) seq%first_avail_time, seq%last_avail_time
endif

! Now read in all the previously defined observations
do i = 1, num_obs
   call read_obs(file_id, num_copies, add_copies, num_qc, add_qc, i, seq%obs(i))
! Also set the key in the obs
   seq%obs(i)%key = i
end do

! Close up the file
call close_file(file_id)

end subroutine read_obs_seq

!-------------------------------------------------

!=================================================

! Functions for the obs_type
!-------------------------------------------------
subroutine init_obs(obs, num_copies, num_qc)

! Sort of a constructor for obs_type
! Should this be public or private just for sequence?

integer, intent(in) :: num_copies, num_qc
type(obs_type), intent(out) :: obs

obs%key = -1
allocate(obs%values(num_copies), &
   obs%qc(num_qc))
obs%prev_time = -1
obs%next_time = -1
obs%cov_group = -1

end subroutine init_obs

!-----------------------------------------------------

subroutine copy_obs(obs1, obs2)

! To be overloaded with =

type(obs_type), intent(out) :: obs1
type(obs_type), intent(in) :: obs2

obs1%key = obs2%key
call copy_obs_def(obs1%def, obs2%def)

!write(*, *) 'in copy obs'
!write(*, *) 'size of obs1, obs2 ', size(obs1%values), size(obs2%values)
if(size(obs1%values) /= size(obs2%values) .or. &
      size(obs1%qc) /= size(obs2%qc)) then
   deallocate(obs1%values)
   deallocate(obs1%qc)
!   write(*, *) 'allocating in copy_obs'
   allocate(obs1%values(size(obs2%values)), obs1%qc(size(obs2%qc)))
endif
obs1%values = obs2%values
obs1%qc = obs2%qc
obs1%prev_time = obs2%prev_time
obs1%next_time = obs2%next_time
obs1%cov_group = obs2%cov_group

!write(*, *) 'done with copy_obs'

end subroutine copy_obs

!-------------------------------------------------

subroutine destroy_obs(obs)

! Free up allocated storage in an observation type
type(obs_type), intent(inout) :: obs

deallocate(obs%values, obs%qc)
call destroy_obs_def(obs%def)

end subroutine destroy_obs

!-------------------------------------------------
subroutine get_obs_def(obs, obs_def)

type(obs_type), intent(in) :: obs
type(obs_def_type), intent(out) :: obs_def

! WARNING: NEED TO DEFINE A COPY ROUTINE FOR OBS_DEF !!!
call copy_obs_def(obs_def, obs%def)

end subroutine get_obs_def

!-------------------------------------------------
subroutine set_obs_def(obs, obs_def)

type(obs_type), intent(out) :: obs
type(obs_def_type), intent(in) :: obs_def

call copy_obs_def(obs%def, obs_def)

end subroutine set_obs_def
!-------------------------------------------------
subroutine get_obs_values(obs, values, copy_indx)


type(obs_type), intent(in) :: obs
real(r8), intent(out) :: values(:)
integer, optional, intent(in) :: copy_indx

if(present(copy_indx)) then
   values(1) = obs%values(copy_indx)
else
   values = obs%values
endif

end subroutine get_obs_values

!-------------------------------------------------

subroutine set_obs_values(obs, values, copy_indx)


type(obs_type), intent(out) :: obs
real(r8), intent(in) :: values(:)
integer, optional, intent(in) :: copy_indx

if(present(copy_indx)) then
   obs%values(copy_indx) = values(1)
else
   obs%values = values
endif

end subroutine set_obs_values

!-------------------------------------------------
subroutine get_qc(obs, qc, qc_indx)


type(obs_type),    intent(in) :: obs
real(r8),         intent(out) :: qc(:)
integer, optional, intent(in) :: qc_indx

if(present(qc_indx)) then
   qc(1) = obs%qc(qc_indx)
else
   qc = obs%qc
endif

end subroutine get_qc

!-------------------------------------------------
subroutine set_qc(obs, qc, qc_indx)

type(obs_type),   intent(out) :: obs
real(r8),          intent(in) :: qc(:)
integer, optional, intent(in) :: qc_indx

if(present(qc_indx)) then
   obs%qc(qc_indx) = qc(1)
else
   obs%qc = qc
endif

end subroutine set_qc

!-------------------------------------------------

subroutine write_obs(obs, file_id, num_copies, num_qc)

! Write out an observation to file, inefficient

type(obs_type), intent(in) :: obs
integer, intent(in) :: file_id, num_copies, num_qc

integer :: i

do i = 1, num_copies
   if(write_binary_obs_sequence) then
      write(file_id) obs%values(i)
   else
      write(file_id, *) obs%values(i)
   endif
end do

do i = 1, num_qc
   if(write_binary_obs_sequence) then
      write(file_id) obs%qc(i)
   else
      write(file_id, *) obs%qc(i)
   endif
end do

if(write_binary_obs_sequence) then
   write(file_id) obs%prev_time, obs%next_time, obs%cov_group
   call write_obs_def(file_id, obs%def, 'unformatted')
else
   write(file_id, *) obs%prev_time, obs%next_time, obs%cov_group
   call write_obs_def(file_id, obs%def)
endif

end subroutine write_obs

!-------------------------------------------------

subroutine read_obs(file_id, num_copies, add_copies, num_qc, add_qc, key, obs)

! Read in observation from file, watch for allocation of storage

integer, intent(in) :: file_id, num_copies, add_copies, num_qc, add_qc, key
type(obs_type), intent(inout) :: obs

integer :: i

! Read in values and qc
if(num_copies > 0) then
   if(read_binary_obs_sequence) then
      do i = 1, num_copies
         read(file_id) obs%values(i)
      end do
   else
      read(file_id, *) obs%values(1:num_copies)
   endif
endif

if(num_qc > 0) then
   if(read_binary_obs_sequence) then
      do i = 1, num_qc
         read(file_id) obs%qc(i)
      end do
   else
      read(file_id, *) obs%qc(1:num_qc)
   endif
endif

if(read_binary_obs_sequence) then
   read(file_id) obs%prev_time, obs%next_time, obs%cov_group
   call read_obs_def(file_id, obs%def, 'unformatted')
else
   read(file_id, *) obs%prev_time, obs%next_time, obs%cov_group
   call read_obs_def(file_id, obs%def)
endif


end subroutine read_obs

!------------------------------------------------------------------------------

subroutine interactive_obs(num_copies, num_qc, obs)

integer,           intent(in) :: num_copies, num_qc
type(obs_type), intent(inout) :: obs

integer :: i

! Does interactive initialization of an observation type

call interactive_obs_def(obs%def)
do i = 1, num_copies
   write(*, *) 'Enter value ', i, 'for this observation'
   read(*, *) obs%values(i)
end do

do i = 1, num_qc
   write(*, *) 'Enter quality control value ', i, 'for this observation'
   read(*, *) obs%qc(i)
end do

! WHAT ABOUT THE COVARIANCE GROUPING???

end subroutine interactive_obs


!---------------------------------------------------------

function get_num_times(seq)

! Returns number of different times for observations in sequence
! Could also be computed as sequence is built?

type(obs_sequence_type), intent(in) :: seq
integer :: get_num_times

integer :: next
type(obs_def_type) :: obs_def
type(time_type) :: time, prev_time

! Just loop through the time sorted sequence and look for different times
get_num_times = 0
next = seq%first_time

do while (next /= -1)
   call get_obs_def(seq%obs(next), obs_def)
   time = get_obs_def_time(obs_def)
   if(get_num_times == 0) then
      get_num_times = 1
   else if(time /= prev_time) then
      get_num_times = get_num_times + 1
   endif
   prev_time = time
   next = seq%obs(next)%next_time
end do
   
end function get_num_times


!-------------------------------------------------
!subroutine get_cov_group
!-------------------------------------------------
!subroutine set_cov_group ???

!=================================================

end module obs_sequence_mod
