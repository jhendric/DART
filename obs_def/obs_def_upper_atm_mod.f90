! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! This module supports the observation types from the AIRS instruments.
! http://winds.jpl.nasa.gov/missions/quikscat/index.cfm

! BEGIN DART PREPROCESS KIND LIST
! SAT_TEMPERATURE,           KIND_TEMPERATURE,           COMMON_CODE
! SAT_TEMPERATURE_ELECTRON,  KIND_TEMPERATURE_ELECTRON,  COMMON_CODE
! SAT_TEMPERATURE_ION,       KIND_TEMPERATURE_ION,       COMMON_CODE
! SAT_DENSITY_NEUTRAL_O3P,   KIND_DENSITY_NEUTRAL_O3P,   COMMON_CODE
! SAT_DENSITY_NEUTRAL_O2,    KIND_DENSITY_NEUTRAL_O2,    COMMON_CODE
! SAT_DENSITY_NEUTRAL_N2,    KIND_DENSITY_NEUTRAL_N2,    COMMON_CODE
! SAT_DENSITY_NEUTRAL_N4S,   KIND_DENSITY_NEUTRAL_N4S,   COMMON_CODE
! SAT_DENSITY_NEUTRAL_NO,    KIND_DENSITY_NEUTRAL_NO,    COMMON_CODE
! SAT_DENSITY_NEUTRAL_N2D,   KIND_DENSITY_NEUTRAL_N2D,   COMMON_CODE
! SAT_DENSITY_NEUTRAL_N2P,   KIND_DENSITY_NEUTRAL_N2P,   COMMON_CODE
! SAT_DENSITY_NEUTRAL_H,     KIND_DENSITY_NEUTRAL_H,     COMMON_CODE
! SAT_DENSITY_NEUTRAL_HE,    KIND_DENSITY_NEUTRAL_HE,    COMMON_CODE
! SAT_DENSITY_NEUTRAL_CO2,   KIND_DENSITY_NEUTRAL_CO2,   COMMON_CODE
! SAT_DENSITY_NEUTRAL_O1D,   KIND_DENSITY_NEUTRAL_O1D,   COMMON_CODE
! SAT_DENSITY_ION_O4SP,      KIND_DENSITY_ION_O4SP,      COMMON_CODE
! SAT_DENSITY_ION_O2P,       KIND_DENSITY_ION_O2P,       COMMON_CODE
! SAT_DENSITY_ION_N2P,       KIND_DENSITY_ION_N2P,       COMMON_CODE
! SAT_DENSITY_ION_NP,        KIND_DENSITY_ION_NP,        COMMON_CODE
! SAT_DENSITY_ION_NOP,       KIND_DENSITY_ION_NOP,       COMMON_CODE
! SAT_DENSITY_ION_O2DP,      KIND_DENSITY_ION_O2DP,      COMMON_CODE
! SAT_DENSITY_ION_O2PP,      KIND_DENSITY_ION_O2PP,      COMMON_CODE
! SAT_DENSITY_ION_HP,        KIND_DENSITY_ION_HP,        COMMON_CODE
! SAT_DENSITY_ION_HEP,       KIND_DENSITY_ION_HEP,       COMMON_CODE
! SAT_DENSITY_ION_E,         KIND_DENSITY_ION_E,         COMMON_CODE
! SAT_VELOCITY_U,            KIND_VELOCITY_U,            COMMON_CODE
! SAT_VELOCITY_V,            KIND_VELOCITY_V,            COMMON_CODE
! SAT_VELOCITY_W,            KIND_VELOCITY_W,            COMMON_CODE
! SAT_VELOCITY_U_ION,        KIND_VELOCITY_U_ION,        COMMON_CODE
! SAT_VELOCITY_V_ION,        KIND_VELOCITY_V_ION,        COMMON_CODE
! SAT_VELOCITY_W_ION,        KIND_VELOCITY_W_ION,        COMMON_CODE
! SAT_VELOCITY_VERTICAL_O3P, KIND_VELOCITY_VERTICAL_O3P, COMMON_CODE
! SAT_VELOCITY_VERTICAL_O2,  KIND_VELOCITY_VERTICAL_O2,  COMMON_CODE
! SAT_VELOCITY_VERTICAL_N2,  KIND_VELOCITY_VERTICAL_N2,  COMMON_CODE
! SAT_VELOCITY_VERTICAL_N4S, KIND_VELOCITY_VERTICAL_N4S, COMMON_CODE
! SAT_VELOCITY_VERTICAL_NO,  KIND_VELOCITY_VERTICAL_NO,  COMMON_CODE
! SAT_F107,                  KIND_1D_PARAMETER,          COMMON_CODE
! SAT_RHO,                   KIND_DENSITY,               COMMON_CODE
! GPS_PROFILE,               KIND_ELECTRON_DENSITY,      COMMON_CODE
! GND_GPS_VTEC,		     KIND_GND_GPS_VTEC
! CHAMP_DENSITY,             KIND_DENSITY
! MIDAS_TEC,                 KIND_VERTICAL_TEC
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_upper_atm_mod, only : get_expected_upper_atm_density
!  use obs_def_upper_atm_mod, only : get_expected_gnd_gps_vtec
!  use obs_def_upper_atm_mod, only : get_expected_vtec
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
! case(CHAMP_DENSITY) 
!      call get_expected_upper_atm_density(state, location, obs_val, istatus)
! case(MIDAS_TEC) 
!      call get_expected_vtec(state, location, obs_val, istatus)
! case(GND_GPS_VTEC)
!      call get_expected_gnd_gps_vtec(state, location, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
! case(CHAMP_DENSITY) 
!      continue
! case(MIDAS_TEC) 
!        continue
! case(GND_GPS_VTEC)
!      continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
! case(CHAMP_DENSITY) 
!      continue
! case(MIDAS_TEC) 
!        continue
! case(GND_GPS_VTEC)
!      continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
! case(CHAMP_DENSITY) 
!      continue
! case(MIDAS_TEC) 
!        continue
! case(GND_GPS_VTEC)
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_upper_atm_mod

use        types_mod, only : r8, MISSING_R8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type, get_location, set_location, &
                             VERTISHEIGHT, VERTISLEVEL
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_ATOMIC_OXYGEN_MIXING_RATIO, &
                             KIND_MOLEC_OXYGEN_MIXING_RATIO, &
                             KIND_TEMPERATURE, &
                             KIND_PRESSURE, &
                             KIND_DENSITY_ION_E, KIND_GND_GPS_VTEC, &
                             KIND_GEOPOTENTIAL_HEIGHT

implicit none
private
public :: get_expected_upper_atm_density, &
          get_expected_gnd_gps_vtec, &
          get_expected_vtec

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"


real(r8), PARAMETER       :: universal_gas_constant = 8314.0_r8 ! [J/K/kmol]
logical, save             :: module_initialized = .false.

contains

subroutine initialize_module
!-----------------------------------------------------------------------------
call register_module(source, revision, revdate)
module_initialized = .true.

end subroutine initialize_module


subroutine get_expected_upper_atm_density(x, location, obs_val, istatus)
!-----------------------------------------------------------------------------
!Given DART state vector and a location, 
!it computes thermospheric neutral density [Kg/m3] 
!The istatus variable should be returned as 0 unless there is a problem
!
real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus
real(r8)                        :: mmro1, mmro2 ! mass mixing ratio 
real(r8)                        :: pressure, temperature 

if ( .not. module_initialized ) call initialize_module

call interpolate(x, location, KIND_ATOMIC_OXYGEN_MIXING_RATIO, mmro1, istatus)
if (istatus /= 0) then
   obs_val = MISSING_R8
   return
endif
call interpolate(x, location, KIND_MOLEC_OXYGEN_MIXING_RATIO, mmro2, istatus)
if (istatus /= 0) then
   obs_val = MISSING_R8
   return
endif
call interpolate(x, location, KIND_PRESSURE, pressure, istatus)
if (istatus /= 0) then
   obs_val = MISSING_R8
   return
endif
call interpolate(x, location, KIND_TEMPERATURE, temperature, istatus)
if (istatus /= 0) then
   obs_val = MISSING_R8
   return
endif

!density [Kg/m3] =  pressure [N/m2] * M [g/mol] / temperature [K] / R [N*m/K/kmol] 
obs_val          =  pressure &
                 /(mmro1/16.0_r8+mmro2/32.0_r8+(1-mmro1-mmro2)/28.0_r8) &
                 /temperature/universal_gas_constant 

end subroutine get_expected_upper_atm_density


subroutine get_expected_gnd_gps_vtec(state_vector, location, obs_val, istatus)
!-----------------------------------------------------------------------------
!Given DART state vector and a location, 
!it computes ground GPS vertical total electron content
!The istatus variable should be returned as 0 unless there is a problem
!
real(r8),            intent(in) :: state_vector(:)
type(location_type), intent(in) :: location
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

! Given a location and the state vector from one of the ensemble members, 
! compute the model-predicted total electron content that would be in the
! integrated column from an instrument looking straight down at the tangent point.
! 'istatus' is the return code.  0 is success; any positive value signals an
! error (different values can be used to indicate different error types).
! Negative istatus values are reserved for internal use only by DART.

integer  :: nLons, nLats, nAlts, iAlt
real(r8), allocatable :: LON(:), LAT(:), ALT(:), IDensityS_ie(:) 
real(r8) :: loc_vals(3)
real(r8) :: tec
type(location_type) :: probe

if ( .not. module_initialized ) call initialize_module

istatus = 36 !initially bad return code
obs_val = MISSING_R8

! something larger than the expected number of vert levels in the model
allocate(ALT(500), IDensityS_ie(500))

loc_vals = get_location(location)

nAlts = 0
LEVELS: do iAlt=1, size(ALT)+1
   ! loop over levels.  if we get to one more than the allocated array size,
   ! this model must have more levels than we expected.  increase array sizes,
   ! recompile, and try again.
   if (iAlt > size(ALT)) then
      call error_handler(E_ERR, 'get_expected_gnd_gps_vtec', 'more than 500 levels in model', &
           source, revision, revdate, &
           text2='increase ALT, IDensityS_ie array sizes in code and recompile')
   endif

   ! At each altitude interpolate the 2D IDensityS_ie to the lon-lat where data 
   ! point is located. After this loop we will have a column centered at the data 
   ! point's lon-lat and at all model altitudes.
   probe = set_location(loc_vals(1), loc_vals(2), real(iAlt, r8), VERTISLEVEL) !probe is where we have data 
   call interpolate(state_vector, probe, KIND_DENSITY_ION_E, IDensityS_ie(iAlt), istatus) 
   if (istatus /= 0) exit LEVELS
   call interpolate(state_vector, probe, KIND_GEOPOTENTIAL_HEIGHT, ALT(iAlt), istatus) 
   if (istatus /= 0) exit LEVELS
   nAlts = nAlts+1
enddo LEVELS

if (nAlts == 0) return

tec=0.0_r8 !start with zero for the summation

do iAlt = 1, nAlts-1 !approximate the integral over the altitude as a sum of trapezoids
   !area of a trapezoid: A = (h2-h1) * (f2+f1)/2
   tec = tec + ( ALT(iAlt+1)-ALT(iAlt) )  * ( IDensityS_ie(iAlt+1)+IDensityS_ie(iAlt) ) /2.0_r8
enddo
obs_val = tec * 10.0**(-16) !units of TEC are "10^16" #electron/m^2 instead of just "1" #electron/m^2

deallocate(ALT, IDensityS_ie)

! Good return code. 
istatus = 0

end subroutine get_expected_gnd_gps_vtec


subroutine get_expected_vtec(x, location, obs_val, istatus)
!-----------------------------------------------------------------------------
!Given DART state vector and a location, 
!it computes thermospheric neutral density [Kg/m3] 
!The istatus variable should be returned as 0 unless there is a problem
!
real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus
real(r8)                        :: mmro1, mmro2 ! mass mixing ratio 
real(r8)                        :: pressure, temperature 

if ( .not. module_initialized ) call initialize_module

call error_handler(E_ERR, 'get_expected_vtec', 'routine needs to be written', &
           source, revision, revdate)

end subroutine get_expected_vtec

end module obs_def_upper_atm_mod
! END DART PREPROCESS MODULE CODE      

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
