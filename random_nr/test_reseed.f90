! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program test_reseed

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r4, r8, i8
use    utilities_mod, only : register_module, error_handler, E_ERR, &
                             initialize_utilities, finalize_utilities, &
                             open_file, close_file, find_namelist_in_file, &
                             check_namelist_read
use time_manager_mod, only : time_type, operator(+), set_time, get_time, &
                             set_calendar_type, print_time, print_date
use   random_seq_mod, only : random_seq_type, init_random_seq, &
                             random_uniform, random_gaussian

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type (random_seq_type) :: seq
type (time_type) :: state_time, delta_time

integer :: i, l, rep, nextseed
integer :: maxreps = 1000
integer, parameter :: ntests = 5
integer :: nloops(ntests) = (/ 1, 10, 100, 1000, 1000000 /)

real(r8), allocatable, dimension(:) :: ranarray
real(r8) :: mean, var

integer :: iunit1,iunit2,iunit3, io
character(len=32) :: filename1,filename2,filename3

! Namelist with default values

integer  :: time_days = 148000
integer  :: time_seconds = 0
integer  :: time_step_seconds = 3600

namelist /test_reseed_nml/ time_days, time_seconds, time_step_seconds

!---------------------------------------------------------------

call initialize_utilities('test_reseed')
call register_module(source,revision,revdate)
call set_calendar_type('GREGORIAN')

! Read the namelist entry
call find_namelist_in_file("input.nml", "test_reseed_nml", iunit1)
read(iunit1, nml = test_reseed_nml, iostat = io)
call check_namelist_read(iunit1, io, "test_reseed_nml")

state_time = set_time(time_seconds, time_days)
delta_time = set_time(time_step_seconds)

print *, ' '
call print_date(state_time, 'setting start time')
call print_time(delta_time, 'delta time')

print *, ' '
print *, 'defining r4 to be ',r4
print *, 'defining r8 to be ',r8
print *, 'doing ', maxreps, ' repeats of ', ntests, ' counts:'
do i=1, ntests
   print *, 'mean[,std] for ', nloops(i), ' uniform  distribution'
enddo
do i=1, ntests
   print *, 'mean[,std] for ', nloops(i), ' gaussian distribution'
enddo
print *, ' '

write(filename1,'(''ran_test_r'',I1,''_'',I5.5,''_compact.out'')') r8, time_step_seconds
write(filename2,'(''ran_test_r'',I1,''_'',I5.5,''_full.out''   )') r8, time_step_seconds
write(filename3,'(''ran_test_r'',I1,''_'',I5.5,''_numbers.out'')') r8, time_step_seconds

iunit1 = open_file(filename1, form='formatted', action='write')
iunit2 = open_file(filename2, form='formatted', action='write')
iunit3 = open_file(filename3, form='formatted', action='write')

write(*,*)'creating compact output file ',trim(filename1)
write(*,*)'creating         output file ',trim(filename2)
write(*,*)'creating         output file ',trim(filename3)

do rep=1, maxreps

   nextseed = generate_seed(state_time)

   write(   *  ,'('' trial '',i5,'' new seed: '', i14)')rep, nextseed
   write(iunit1,'('' trial '',i5,'' new seed: '', i14)')rep, nextseed
   write(iunit2,'('' trial '',i5,'' new seed: '', i14)')rep, nextseed

! ---------

   do l=1, ntests

      call init_random_seq(seq, nextseed)

      allocate(ranarray(nloops(l)))

      mean = 0.0_r8
      var  = 0.0_r8
      do i=1, nloops(l)
         ranarray(i) = random_uniform(seq)
         mean = mean + ranarray(i)
      enddo

      mean = mean/real(nloops(l),r8)

      do i=1, nloops(l)
         var = var + (mean - ranarray(i))**2
      enddo

      if (l /= 1) then
         write(iunit1,'(''uniform  '',i7,1x,f9.6,  1x,f9.6  )') nloops(l), mean, sqrt(var/(nloops(l)-1))
         write(iunit2,'(''uniform  '',i7,1x,f19.16,1x,f19.16)') nloops(l), mean, sqrt(var/(nloops(l)-1))
      else
         write(iunit1,'(''uniform  '',i7,1x,f9.6            )') nloops(l), mean
         write(iunit2,'(''uniform  '',i7,1x,f19.16          )') nloops(l), mean
      endif

      deallocate(ranarray)

   enddo

! ---------

   write(  *   ,*)
   write(iunit1,*)
   write(iunit2,*)

! ---------

   do l=1, ntests

      call init_random_seq(seq, nextseed)

      allocate(ranarray(nloops(l)))

      mean = 0.0_r8
      var  = 0.0_r8
      do i=1, nloops(l)
         ranarray(i) = random_gaussian(seq, 0.0_r8, 1.0_r8)
         mean = mean + ranarray(i)
      enddo

      mean = mean/real(nloops(l),r8)

      do i=1, nloops(l)
         var = var + (mean - ranarray(i))**2
      enddo

      if (l /= 1) then
         write(iunit1,'(''gaussian '',i7,1x,f9.6,  1x,f9.6  )') nloops(l), mean, sqrt(var/(nloops(l)-1))
         write(iunit2,'(''gaussian '',i7,1x,f19.16,1x,f19.16)') nloops(l), mean, sqrt(var/(nloops(l)-1))
      else
         write(iunit1,'(''gaussian '',i7,1x,f9.6            )') nloops(l), mean
         write(iunit2,'(''gaussian '',i7,1x,f19.16          )') nloops(l), mean
      endif

      if ( nloops(l) == 100 ) then
         do i=1, nloops(l)
            write(iunit3,*) ranarray(i)
         enddo
      endif

      deallocate(ranarray)

   enddo

! ---------

   write(  *   ,*)
   write(iunit1,*)
   write(iunit2,*)

! ---------

   state_time = state_time + delta_time

enddo

call close_file(iunit1)
call close_file(iunit2)
call close_file(iunit3)

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
