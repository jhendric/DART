program location_test

! <next three lines automatically updated by CVS, do not edit>
!  $Source$
!  $Revision$
!  $Date$

! Simple test program to exercise twod_sphere location module.

use location_mod
use types_mod,     only : r8, PI
use utilities_mod, only : get_unit

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(location_type) :: loc1, loc2
integer  :: iunit, i
real(r8) :: loc2_val, lon, lat

! Test distribution of random locations
call interactive_location(loc2)
do i = 1, 10
   loc1 = loc2
   call interactive_location(loc2)
   write(*, *) 'location is ', get_location(loc2)
   write(*, *) 'distance to previous is ', get_dist(loc1, loc2)
   write(*, *) 'distance to previous is ', get_dist(loc2, loc1)
end do

! Open an output file
iunit = get_unit()
open(iunit, file = 'location_test_file')

! Set the first location
call interactive_location(loc1)
call interactive_location(loc2)

! Write this location to the file
call write_location(iunit, loc1)
call write_location(iunit, loc2)

close(iunit)

! Now read them back in and compute the distances from loc1
open(iunit, file = 'location_test_file')
loc1 = read_location(iunit)
loc2 = read_location(iunit)

write(*, *) 'location 1 is ', get_location(loc1)
write(*, *) 'location 2 is ', get_location(loc2)
write(*, *) 'distance is   ', get_dist(loc1, loc2)

close(iunit)

end program location_test

