module obs_set_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!

use types_mod
use set_def_list_mod, only : set_def_list_type, get_total_num_obs
use time_manager_mod, only : time_type, read_time, write_time

private

public obs_set_type, init_obs_set, get_obs_set_time, get_obs_values,&
   set_obs_values, set_obs_missing, set_obs_set_time, &
   contains_data, obs_value_missing, &
   read_obs_set, write_obs_set

type obs_set_type
   private
   real(r8), pointer :: obs(:, :)
   logical, pointer :: missing(:, :)
   integer :: num_copies
   integer :: num_obs
   type(time_type) :: time
   integer :: def_index
!   type(obs_set_def_type) :: def
end type obs_set_type


contains

!=======================================================


function init_obs_set(set_def_list, index, num_copies_in)
!--------------------------------------------------------------------
!
! Creates and initializes storage for an obs_set associated
! with the obs_set_def index in obs_set_def_list and with num_copies of the
! observations associated (default is 1). 

type(obs_set_type) :: init_obs_set
type(set_def_list_type), intent(in) :: set_def_list
integer, intent(in) :: index
integer, optional, intent(in) :: num_copies_in

integer :: num_copies, num_obs

! Begin by getting the number of copies
num_copies = 1
if(present(num_copies_in)) num_copies = num_copies_in
if(num_copies < 0) then
   write(*, *) 'Error: Negative num_copies in init_obs_set'
   stop
endif

! Set the number of copies
init_obs_set%num_copies = num_copies

! Assign the set definition
init_obs_set%def_index = index

! Get the total number of obs in the set
num_obs = get_total_num_obs(set_def_list, index)
init_obs_set%num_obs = num_obs

! Allocate storage for obs; need to verify that the 0 size allocation
! is legal F90
allocate(init_obs_set%obs(num_obs, num_copies), init_obs_set%missing(num_obs, num_copies))

end function init_obs_set



function get_obs_set_time(set)
!-------------------------------------------------------
!
! Returns the time associated with this observation set

implicit none

type(time_type) :: get_obs_set_time
type(obs_set_type), intent(in) :: set

get_obs_set_time = set%time

end function get_obs_set_time




subroutine get_obs_values(set, obs, index_in)
!-----------------------------------------------------------------
!
! Returns the values of observations from the set. If the set has
! multiple observations per definition, the optional argument index
! selects which set to return.

implicit none

type(obs_set_type), intent(in) :: set
real(r8), intent(out) :: obs(:)
integer, optional, intent(in) :: index_in

integer :: index

! Get the appropriate index
index = 1
if(present(index_in)) index = index_in
if(index < 1 .or. index > set%num_copies) then
   write(*, *) 'Error: Out of range index in init_obs_set'
   stop
endif

! Next make sure there's enough room in obs
if(size(obs) < set%num_obs) then
   write(*, *) 'Error: obs array too small in get_obs_values'
   stop
endif

!Copy the obs
obs = set%obs(:, index)

end subroutine get_obs_values



subroutine set_obs_values(set, obs, index_in)
!--------------------------------------------------------------------
!
! Sets the obs values; index is optional with default 1.

implicit none

type(obs_set_type), intent(inout) :: set
real(r8), intent(in) :: obs(:)
integer, optional, intent(in) :: index_in

integer :: index

! Get the appropriate index
index = 1
if(present(index_in)) index = index_in
if(index < 1 .or. index > set%num_copies) then
   write(*, *) 'Error: Out of range index in init_obs_set'
   stop
endif

! Make sure the obs array is the right size
if(size(obs) /= set%num_obs) then
   write(*, *) 'Error: obs array wrong size in set_obs_values'
   stop
endif

!Copy the obs
set%obs(:, index) = obs

end subroutine set_obs_values




subroutine set_obs_set_missing(set, missing, index_in)
!--------------------------------------------------------------------
!
! Sets the obs values; index is optional with default 1.

implicit none

type(obs_set_type), intent(inout) :: set
logical, intent(in) :: missing(:)
integer, optional, intent(in) :: index_in

integer :: index

! Get the appropriate index
index = 1
if(present(index_in)) index = index_in
if(index < 1 .or. index > set%num_copies) then
   write(*, *) 'Error: Out of range index in init_obs_set'
   stop
endif

! Make sure the missing array is the right size
if(size(missing) /= set%num_obs) then
   write(*, *) 'Error: missing array wrong size in set_obs_set_missing'
   stop
endif

! Set the data missing
set%missing(:, index) = missing

end subroutine set_obs_set_missing



subroutine set_obs_set_time(set, time)
!----------------------------------------------------
!
! Set the time for the obs_set

implicit none

type(obs_set_type), intent(inout) :: set
type(time_type), intent(in) :: time

set%time = time

end subroutine set_obs_set_time




function contains_data(set)
!--------------------------------------------------------------------
!
! Returns true if the number of associated data copies with the set is
! not 0.

implicit none

logical :: contains_data
type(obs_set_type), intent(in) :: set

contains_data = set%num_copies > 0

end function contains_data




subroutine obs_value_missing(set, missing, index_in)
!--------------------------------------------------------------------
!
! Returns true if the data associated with the copy is missing. Default
! for index_in is 1.

implicit none

type(obs_set_type), intent(in) :: set
logical, intent(out) :: missing(:)
integer, optional, intent(in) :: index_in

integer :: index

! Get the appropriate index
index = 1
if(present(index_in)) index = index_in
if(index < 1 .or. index > set%num_copies) then
   write(*, *) 'Error: Out of range index in init_obs_set'
   stop
endif

! Next make sure there's enough room in obs
if(size(missing) < set%num_obs) then
   write(*, *) 'Error: missing array too small in obs_value_missing'
   stop
endif

! Copy the missing data
missing = set%missing(:, index)

end subroutine obs_value_missing


function read_obs_set(file_id)
!------------------------------------------------------------------------
!
! Reads and obs_set from a file

implicit none

type(obs_set_type) :: read_obs_set
integer, intent(in) :: file_id

character*5 :: header
integer :: num_obs, num_copies, i

! Read the header and verify 
read(file_id, *) header
if(header /= 'obset') then
   write(*, *) 'Error: Expected "obset" in header in read_obs_set' 
   stop
end if

! Read the obs_set def index
read(file_id, *) read_obs_set%def_index

! Read the number of obs and the number of copies
read(file_id, *) num_obs, num_copies
read_obs_set%num_obs = num_obs
read_obs_set%num_copies = num_copies

! Allocate space
allocate(read_obs_set%obs(num_obs, num_copies), read_obs_set%missing(num_obs, num_copies))

! Read the data for each copy in turn
do i = 1, num_copies
   read(file_id, *) read_obs_set%obs(:, i)
end do

! Read the missing fields for each copy in turn
do i = 1, num_copies
   read(file_id, *) read_obs_set%missing(:, i)
end do

! Read the time 
read_obs_set%time = read_time(file_id)

end function read_obs_set



subroutine write_obs_set(file_id, set)
!-------------------------------------------------------------------------
!
! Writes an obs_set with all copies of its observations to file, currently
! represented by integer unit number.

implicit none

type(obs_set_type), intent(in) :: set
integer, intent(in) :: file_id

integer :: i

! First write ascii header saying set is coming
write(file_id, *) 'obset'

! Write the obs_set_def_index
write(file_id, *) set%def_indeX

! The number of obs followed by the number of copies
write(file_id, *) set%num_obs, set%num_copies

! Write out the data for each copy in turn
do i = 1, set%num_copies
   write(file_id, *) set%obs(:, i)
end do

! Write out the missing fields for each copy in turn
do i = 1, set%num_copies
   write(file_id, *) set%missing(:, i)
end do

! Write out the time associated with this observation
call write_time(file_id, set%time)

end subroutine write_obs_set




end module obs_set_mod
