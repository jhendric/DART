! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! This module supports the observation types from the AIRS instruments.
! http://winds.jpl.nasa.gov/missions/quikscat/index.cfm

! "BOB" is simply a placeholder for any kind of observation platform ...
! "SAT" would be more professional ...

! BEGIN DART PREPROCESS KIND LIST
! SAT_TEMPERATURE,           KIND_TEMPERATURE, COMMON_CODE
! SAT_TEMPERATURE_ELECTRON,  KIND_TEMPERATURE_ELECTRON, COMMON_CODE
! SAT_TEMPERATURE_ION,       KIND_TEMPERATURE_ION, COMMON_CODE
! SAT_DENSITY_NEUTRAL_O3P,   KIND_DENSITY_NEUTRAL_O3P, COMMON_CODE
! SAT_DENSITY_NEUTRAL_O2,    KIND_DENSITY_NEUTRAL_O2, COMMON_CODE
! SAT_DENSITY_NEUTRAL_N2,    KIND_DENSITY_NEUTRAL_N2, COMMON_CODE
! SAT_DENSITY_NEUTRAL_N4S,   KIND_DENSITY_NEUTRAL_N4S, COMMON_CODE
! SAT_DENSITY_NEUTRAL_NO,    KIND_DENSITY_NEUTRAL_NO, COMMON_CODE
! SAT_DENSITY_NEUTRAL_N2D,   KIND_DENSITY_NEUTRAL_N2D, COMMON_CODE
! SAT_DENSITY_NEUTRAL_N2P,   KIND_DENSITY_NEUTRAL_N2P, COMMON_CODE
! SAT_DENSITY_NEUTRAL_H,     KIND_DENSITY_NEUTRAL_H, COMMON_CODE
! SAT_DENSITY_NEUTRAL_HE,    KIND_DENSITY_NEUTRAL_HE, COMMON_CODE
! SAT_DENSITY_NEUTRAL_AR,    KIND_DENSITY_NEUTRAL_AR, COMMON_CODE
! SAT_DENSITY_NEUTRAL_O1D,   KIND_DENSITY_NEUTRAL_O1D, COMMON_CODE
! SAT_DENSITY_ION_O4SP,      KIND_DENSITY_ION_O4SP, COMMON_CODE
! SAT_DENSITY_ION_O2P,       KIND_DENSITY_ION_O2P, COMMON_CODE
! SAT_DENSITY_ION_N2P,       KIND_DENSITY_ION_N2P, COMMON_CODE
! SAT_DENSITY_ION_NP,        KIND_DENSITY_ION_NP, COMMON_CODE
! SAT_DENSITY_ION_NOP,       KIND_DENSITY_ION_NOP, COMMON_CODE
! SAT_DENSITY_ION_O2DP,      KIND_DENSITY_ION_O2DP, COMMON_CODE
! SAT_DENSITY_ION_O2PP,      KIND_DENSITY_ION_O2PP, COMMON_CODE
! SAT_DENSITY_ION_HP,        KIND_DENSITY_ION_HP, COMMON_CODE
! SAT_DENSITY_ION_HEP,       KIND_DENSITY_ION_HEP, COMMON_CODE
! SAT_DENSITY_ION_E,         KIND_DENSITY_ION_E, COMMON_CODE
! SAT_VELOCITY_U,            KIND_VELOCITY_U, COMMON_CODE
! SAT_VELOCITY_V,            KIND_VELOCITY_V, COMMON_CODE
! SAT_VELOCITY_W,            KIND_VELOCITY_W, COMMON_CODE
! SAT_VELOCITY_U_ION,        KIND_VELOCITY_U_ION, COMMON_CODE
! SAT_VELOCITY_V_ION,        KIND_VELOCITY_V_ION, COMMON_CODE
! SAT_VELOCITY_W_ION,        KIND_VELOCITY_W_ION, COMMON_CODE
! SAT_VELOCITY_VERTICAL_O3P, KIND_VELOCITY_VERTICAL_O3P, COMMON_CODE
! SAT_VELOCITY_VERTICAL_O2,  KIND_VELOCITY_VERTICAL_O2, COMMON_CODE
! SAT_VELOCITY_VERTICAL_N2,  KIND_VELOCITY_VERTICAL_N2, COMMON_CODE
! SAT_VELOCITY_VERTICAL_N4S, KIND_VELOCITY_VERTICAL_N4S, COMMON_CODE
! SAT_VELOCITY_VERTICAL_NO,  KIND_VELOCITY_VERTICAL_NO, COMMON_CODE
! CHAMP_DENSITY,		     KIND_DENSITY
! GPS_PROFILE,               KIND_ELECTRON_DENSITY, COMMON_CODE
! END DART PREPROCESS KIND LIST

! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_upper_atm_mod, only : get_expected_upper_atm_density
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE

! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
! case(CHAMP_DENSITY) 
!      call get_expected_upper_atm_density(state, location, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF

! BEGIN DART PREPROCESS READ_OBS_DEF
! case(CHAMP_DENSITY) 
!      continue
! END DART PREPROCESS READ_OBS_DEF

! BEGIN DART PREPROCESS WRITE_OBS_DEF
! case(CHAMP_DENSITY) 
!      continue
! END DART PREPROCESS WRITE_OBS_DEF

! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
! case(CHAMP_DENSITY) 
!      continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF

! BEGIN DART PREPROCESS MODULE CODE
module obs_def_upper_atm_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r8, missing_r8
use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG
use     location_mod, only : location_type
use  assim_model_mod, only : interpolate
use     obs_kind_mod, only : KIND_ATOMIC_OXYGEN_MIXING_RATIO, &
                             KIND_MOLEC_OXYGEN_MIXING_RATIO, &
                             KIND_TEMPERATURE, &
                             KIND_PRESSURE
implicit none
private
public                    :: get_expected_upper_atm_density

! version controlled file description for error handling, do not edit
character(len=128) :: &
source   = "$URL$", &
revision = "$Revision$", &
revdate  = "$Date$"


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
   obs_val = missing_r8
   return
endif
call interpolate(x, location, KIND_MOLEC_OXYGEN_MIXING_RATIO, mmro2, istatus)
if (istatus /= 0) then
   obs_val = missing_r8
   return
endif
call interpolate(x, location, KIND_PRESSURE, pressure, istatus)
if (istatus /= 0) then
   obs_val = missing_r8
   return
endif
call interpolate(x, location, KIND_TEMPERATURE, temperature, istatus)
if (istatus /= 0) then
   obs_val = missing_r8
   return
endif

!density [Kg/m3] =  pressure [N/m2] * M [g/mol] / temperature [K] / R [N*m/K/kmol] 
obs_val          =  pressure &
                 /(mmro1/16.0_r8+mmro2/32.0_r8+(1-mmro1-mmro2)/28.0_r8) &
                 /temperature/universal_gas_constant 

end subroutine get_expected_upper_atm_density


end module obs_def_upper_atm_mod
! END DART PREPROCESS MODULE CODE      
