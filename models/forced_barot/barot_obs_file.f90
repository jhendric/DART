! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
program barot_obs_file

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 
!
! Generates a regular set, currently hardcoded to 40 by 30,
! of observation file locations for the forced barot model.

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

write(*, *) 1200
do i = 1, 40
   do j = 1, 30
      lon = ((i - 0.5) / 40.0) * 360.0 
! Don't want to have to interpolate lats
      lat = -86.0 + (j / 30.0) * 170.0
      write(*, *) lon, lat, (1e6)**2
   end do
end do
end program barot_obs_file
