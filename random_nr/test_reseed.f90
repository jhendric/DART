! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program test_reseed

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use     types_mod, only : r4, r8, digits12, i8
use utilities_mod, only : register_module, error_handler, E_ERR, &
                          initialize_utilities, finalize_utilities
use time_manager_mod, only : time_type, operator(+), set_time, get_time, &
                             set_calendar_type, print_time, print_date
use random_seq_mod, only : random_seq_type, init_random_seq, &
                           random_uniform, random_gaussian

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type (random_seq_type) :: seq
type (time_type) :: state_time, delta_time
integer :: i, l, rep
integer :: maxreps = 1000
integer, parameter :: ntests = 5
integer :: nloops(ntests) = (/ 1, 10, 100, 1000, 1000000 /)

real(r8) :: next_val, mean
integer :: nextseed


call initialize_utilities('test_reseed')
call register_module(source,revision,revdate)


call set_calendar_type('GREGORIAN')

state_time = set_time(0, 148000)
delta_time = set_time(43200)

print *, ' '
call print_date(state_time, 'setting start time')
call print_time(delta_time, 'delta time')

print *, ' '
print *, 'doing ', maxreps, ' repeats of ', ntests, ' counts:'
do i=1, ntests
   print *, 'mean for ', nloops(i), ' uniform distribution'
enddo
do i=1, ntests
   print *, 'mean for ', nloops(i), ' gaussian distribution'
enddo
print *, ' '

do rep=1, maxreps

   nextseed = generate_seed(state_time)

   print *, ' '
   print *, 'new seed: ', nextseed

! ---------

   do l=1, ntests

      call init_random_seq(seq, nextseed)

      mean = 0.0_r8
      do i=1, nloops(l)
    
         next_val = random_uniform(seq)
         mean = mean + next_val
      
      enddo

      print *, mean/nloops(l), nloops(l), ' uniform'

   enddo

! ---------

   print *, ' '

! ---------

   do l=1, ntests

      call init_random_seq(seq, nextseed)

      mean = 0.0_r8
      do i=1, nloops(l)
 
         next_val = random_gaussian(seq, 0.0_r8, 1.0_r8)
         mean = mean + next_val
   
      enddo

      print *, mean/nloops(l), nloops(l), ' gaussian'

   enddo

! ---------

   print *, ' '

! ---------

   state_time = state_time + delta_time

enddo

call finalize_utilities()

contains

!-------------------------------------------------------------------------

function generate_seed(state_time)
! use the state time to set the seed for the (repeatable) random sequence 

type(time_type), intent(in) :: state_time
integer                     :: generate_seed

integer :: days,seconds,bigint
integer(kind=i8) :: bigtime,bigone,bigtwo

call get_time(state_time, seconds, days)

bigtime       = int(days,i8)*100000_i8 + int(seconds,i8)
bigint        = huge(seconds)        ! biggest 32bit integer
bigone        = int(bigint,i8)       ! coerce to 64bit integer
bigtwo        = mod(bigtime,bigone)  ! modulo arith on 64bit integers
generate_seed = -1 * int(bigtwo)     ! coerce back to 32bit integer

!write(*,*)'TJH DEBUG generate_seed ',bigtime,bigint,generate_seed

end function generate_seed

end program test_reseed
