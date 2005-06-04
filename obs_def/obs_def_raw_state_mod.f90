! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module obs_def_raw_state_mod

use        types_mod, only : r8, missing_i, missing_r8, RAD2DEG
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, file_exist, &
                             open_file, check_nml_error, logfileunit, close_file
use     location_mod, only : location_type, set_location, get_location 
!WRF use     location_mod, only : query_location
use time_manager_mod, only : time_type, read_time, write_time, set_time, set_time_missing, &
                             interactive_time
use  assim_model_mod, only : get_state_meta_data, interpolate
use cov_cutoff_mod,   only : comp_cov_factor

implicit none

public write_1d_integral, read_1d_integral, interactive_1d_integral, get_expected_1d_integral

! Storage for the special information required for observations of this type
integer, parameter                      :: max_1d_integral_obs = 100
integer                                 :: num_1d_integral_obs = 0
real(r8)                                :: half_width(max_1d_integral_obs)
integer, dimension(max_1d_integral_obs) :: num_points, localization_type

! For now, read in all info on first read call, write all info on first write call
logical :: already_read = .false., already_written = .false.

contains

!----------------------------------------------------------------------

subroutine write_1d_integral(key, ifile, fileformat)

integer, intent(in)             :: key, ifile
character(len=32), intent(in)   :: fileformat

integer :: i

! Philosophy, dump ALL information about this special obs_type at once???
! For now, this means you can only write ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this writing
if(.not. already_written) then
   already_written = .true.   
   ! Write out the number of 1d_integral obs descriptions
   SELECT CASE (fileformat)
      CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
         write(ifile) num_1d_integral_obs
      CASE DEFAULT
         write(ifile, *) num_1d_integral_obs
   END SELECT

   ! Write out the half_width, num_points, and localization_type for each  
   do i = 1, num_1d_integral_obs
      SELECT CASE (fileformat)
         CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
            write(ifile) half_width(i), num_points(i), localization_type(i) 
         CASE DEFAULT
            write(ifile, *) half_width(i), num_points(i), localization_type(i) 
      END SELECT
   end do
endif

! Write out the obs_def key for this observation
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      write(ifile) key
   CASE DEFAULT
      write(ifile, *) key
END SELECT

end subroutine write_1d_integral

!----------------------------------------------------------------------

subroutine read_1d_integral(key, ifile, fileformat)

integer, intent(out)            :: key
integer, intent(in)             :: ifile
character(len=32), intent(in)   :: fileformat

integer :: i

! Philosophy, read ALL information about this special obs_type at once???
! For now, this means you can only read ONCE (that's all we're doing 3 June 05)
! Toggle the flag to control this reading
if(.not. already_read) then
   already_read = .true.   
   ! Read the number of 1d_integral obs descriptions
   SELECT CASE (fileformat)
      CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
         read(ifile) num_1d_integral_obs
      CASE DEFAULT
         read(ifile, *) num_1d_integral_obs
   END SELECT
   
   ! Read the half_width, num_points, and localization_type for each  
   do i = 1, num_1d_integral_obs
      SELECT CASE (fileformat)
         CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
            read(ifile) half_width(i), num_points(i), localization_type(i) 
         CASE DEFAULT
            read(ifile, *) half_width(i), num_points(i), localization_type(i) 
      END SELECT
   end do
endif

! Read in the key for this particular observation
SELECT CASE (fileformat)
   CASE ("unf", "UNF", "unformatted", "UNFORMATTED")
      read(ifile) key
   CASE DEFAULT
      read(ifile, *) key
END SELECT

end subroutine read_1d_integral

!----------------------------------------------------------------------

subroutine interactive_1d_integral(key)

integer, intent(out) :: key

! Initializes the specialized part of a 1d_integral observation
! Passes back up the key for this one

! Make sure there's enough space, if not die for now (clean later)
if(num_1d_integral_obs >= max_1d_integral_obs) then
   ! PUT IN ERROR HANDLER CALL
   stop
endif

! Increment the index
num_1d_integral_obs = num_1d_integral_obs + 1
key = num_1d_integral_obs

! Otherwise, prompt for input for the three required beasts
write(*, *) 'Creating an interactive_1d_integral observation'
write(*, *) 'Input half width of integral '
read(*, *) half_width(num_1d_integral_obs)
write(*, *) 'Input the number of evaluation points (??? recommended) '
read(*, *) num_points(num_1d_integral_obs)
write(*, *) 'Input localization type: 1=Gaspari-Cohn; 2=Boxcar; 3=Ramped Boxcar'
read(*, *) localization_type(num_1d_integral_obs)

end subroutine interactive_1d_integral

!----------------------------------------------------------------------

subroutine get_expected_1d_integral(state, location, key, val, istatus)

real(r8), intent(in)            :: state(:)
type(location_type), intent(in) :: location
integer, intent(in)             :: key
real(r8), intent(out)           :: val
integer, intent(out)            :: istatus

integer :: i
real(r8) :: range, loc, bottom, dx, x, sum, dist, weight, weight_sum
type(location_type) :: location2

! Figure out the total range of the integrated funtion (1 is max)
range = 4.0_r8 * half_width(key)
if(range > 1.0_r8) range = 1.0_r8
!write(*, *) 'range is ', range

! Get the location value
loc = get_location(location)
!write(*, *) 'loc base is ', loc

! Compute the bottom and top of the range
bottom = loc - range / 2.0_r8
if(bottom < 0.0_r8) bottom = bottom + 1.0_r8
!write(*, *) 'bottom is ', bottom

! Next figure out where to put all the points
dx = range / (num_points(key) - 1)
!write(*, *) 'dx is ', dx

! Loop to compute the value at each point, then multiply by localization
! to get weighted integral
sum = 0.0
weight_sum = 0.0
do i = 1, num_points(key)
   x = bottom + (i - 1) * dx
   if(x > 1.0_r8) x = x - 1.0_r8
!write(*, *) 'location for int ', i, 'is ', x
   location2 = set_location(x)
   call interpolate(state, location2, 1, val, istatus)
   dist = abs(loc - x)
   if(dist > 0.5_r8) dist = 1.0_r8 - dist
!write(*, *) 'dist ', i, dist
   weight = comp_cov_factor(dist, half_width(key), localization_type(key))
!write(*, *) 'weight ', i, weight
   sum = sum + weight * val
   weight_sum = weight_sum + weight
end do

val = sum / weight_sum

!write(*, *) 'get_expected_1d_integral key is ', key
!write(*, *) half_width(key), num_points(key), localization_type(key)

end subroutine get_expected_1d_integral

!----------------------------------------------------------------------

end module obs_def_raw_state_mod
