! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

! BEGIN DART PREPROCESS KIND LIST
!SOIL_TEMPERATURE,               KIND_SOIL_TEMPERATURE,         COMMON_CODE
!LAYER_LIQUID_WATER,             KIND_LIQUID_WATER,             COMMON_CODE
!LAYER_ICE,                      KIND_ICE,                      COMMON_CODE
!SNOW_THICKNESS,                 KIND_SNOW_THICKNESS,           COMMON_CODE
!SNOW_WATER,                     KIND_SNOW_WATER,               COMMON_CODE
!MODIS_SNOWCOVER_FRAC,           KIND_SNOWCOVER_FRAC,           COMMON_CODE
!LEAF_CARBON,                    KIND_LEAF_CARBON,              COMMON_CODE
!LEAF_AREA_INDEX,                KIND_LEAF_AREA_INDEX,          COMMON_CODE
!TOWER_AIR_TEMPERATURE,          KIND_TEMPERATURE,              COMMON_CODE
!TOWER_SOIL_TEMPERATURE,         KIND_TEMPERATURE,              COMMON_CODE
!TOWER_U_WIND_COMPONENT,         KIND_U_WIND_COMPONENT,         COMMON_CODE
!TOWER_V_WIND_COMPONENT,         KIND_V_WIND_COMPONENT,         COMMON_CODE
!TOWER_LATENT_HEAT_FLUX,         KIND_LATENT_HEAT_FLUX
!TOWER_SENSIBLE_HEAT_FLUX,       KIND_SENSIBLE_HEAT_FLUX
!TOWER_GLOBAL_RADIATION,         KIND_RADIATION,                COMMON_CODE
!TOWER_NETC_ECO_EXCHANGE,        KIND_NET_CARBON_PRODUCTION
!TOWER_NET_CARBON_FLUX,          KIND_NET_CARBON_FLUX,          COMMON_CODE
! END DART PREPROCESS KIND LIST

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!  use obs_def_tower_mod, only : get_expected_latent_heat_flux
!  use obs_def_tower_mod, only : get_expected_sensible_heat_flux
!  use obs_def_tower_mod, only : get_expected_net_C_production
! END DART PREPROCESS USE OF SPECIAL OBS_DEF MODULE
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!  case(TOWER_LATENT_HEAT_FLUX)
!     call get_expected_latent_heat_flux(state, state_time, ens_index, location, obs_def%key, obs_val, istatus)
!  case(TOWER_SENSIBLE_HEAT_FLUX)
!     call get_expected_sensible_heat_flux(state, state_time, ens_index, location, obs_def%key, obs_val, istatus)
!  case(TOWER_NETC_ECO_EXCHANGE)
!     call get_expected_net_C_production(state, state_time, ens_index, location, obs_def%key, obs_val, istatus)
! END DART PREPROCESS GET_EXPECTED_OBS_FROM_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS READ_OBS_DEF
!    case(TOWER_LATENT_HEAT_FLUX,TOWER_SENSIBLE_HEAT_FLUX,TOWER_NETC_ECO_EXCHANGE)
!       continue
! END DART PREPROCESS READ_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS WRITE_OBS_DEF
!    case(TOWER_LATENT_HEAT_FLUX,TOWER_SENSIBLE_HEAT_FLUX,TOWER_NETC_ECO_EXCHANGE)
!       continue
! END DART PREPROCESS WRITE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS INTERACTIVE_OBS_DEF
!    case(TOWER_LATENT_HEAT_FLUX,TOWER_SENSIBLE_HEAT_FLUX,TOWER_NETC_ECO_EXCHANGE)
!       continue
! END DART PREPROCESS INTERACTIVE_OBS_DEF
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! BEGIN DART PREPROCESS MODULE CODE
module obs_def_tower_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r8, missing_r8, PI, deg2rad
use     location_mod, only : location_type
use time_manager_mod, only : time_type
use    utilities_mod, only : register_module, E_ERR, error_handler

implicit none
private

public :: get_expected_latent_heat_flux,   &
          get_expected_sensible_heat_flux, &
          get_expected_net_C_production

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

logical :: module_initialized = .false.

character(len=129) :: string1,string2,string3

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Start of executable routines
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine initialize_module

! Called once to set values and allocate space

! integer :: iunit, io, rc

! Prevent multiple calls from executing this code more than once.
if (module_initialized) return

module_initialized = .true.

! Log the version of this source file.
call register_module(source, revision, revdate)

! ! Read the namelist entry.
! call find_namelist_in_file("input.nml", "obs_def_tower_nml", iunit)
! read(iunit, nml = obs_def_ocean_nml, iostat = io)
! call check_namelist_read(iunit, io, "obs_def_ocean_nml")

! ! Record the namelist values used for the run ... 
! if (do_nml_file()) write(nmlfileunit, nml=obs_def_ocean_nml)
! if (do_nml_term()) write(     *     , nml=obs_def_ocean_nml)

end subroutine initialize_module


subroutine get_expected_latent_heat_flux(state, state_time, ens_index, location, obs_key, obs_val, istatus)
! the routine must return values for:
! obs_val -- the computed forward operator value
! istatus -- return code: 0=ok, > 0 is error, < 0 reserved for system use

real(r8),            intent(in)  :: state(:)
type(time_type),     intent(in)  :: state_time
integer,             intent(in)  :: ens_index
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_key
real(r8),            intent(out) :: obs_val
integer,             intent(out) :: istatus

obs_val = MISSING_R8
istatus = 1

call error_handler(E_ERR, 'get_expected_latent_heat_flux', &
            'not implemented yet.', &
             source, revision, revdate)

end subroutine get_expected_latent_heat_flux



subroutine get_expected_sensible_heat_flux(state, state_time, ens_index, location, obs_key, obs_val, istatus)
! the routine must return values for:
! obs_val -- the computed forward operator value
! istatus -- return code: 0=ok, > 0 is error, < 0 reserved for system use
real(r8),            intent(in)  :: state(:)
type(time_type),     intent(in)  :: state_time
integer,             intent(in)  :: ens_index
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_key
real(r8),            intent(out) :: obs_val
integer,             intent(out) :: istatus

obs_val = MISSING_R8
istatus = 1

call error_handler(E_ERR, 'get_expected_sensible_heat_flux', &
            'not implemented yet.', &
             source, revision, revdate)

end subroutine get_expected_sensible_heat_flux


subroutine get_expected_net_C_production(state, state_time, ens_index, location, obs_key, obs_val, istatus)
! the routine must return values for:
! obs_val -- the computed forward operator value
! istatus -- return code: 0=ok, > 0 is error, < 0 reserved for system use
real(r8),            intent(in)  :: state(:)
type(time_type),     intent(in)  :: state_time
integer,             intent(in)  :: ens_index
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_key
real(r8),            intent(out) :: obs_val
integer,             intent(out) :: istatus

obs_val = MISSING_R8
istatus = 1

call error_handler(E_ERR, 'get_expected_net_C_production', &
            'not implemented yet.', &
             source, revision, revdate)

end subroutine get_expected_net_C_production



end module obs_def_tower_mod

! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------
