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
!TOWER_GLOBAL_RADIATION,         KIND_RADIATION,                COMMON_CODE
!TOWER_NET_CARBON_FLUX,          KIND_NET_CARBON_FLUX,          COMMON_CODE
!TOWER_LATENT_HEAT_FLUX,         KIND_LATENT_HEAT_FLUX
!TOWER_SENSIBLE_HEAT_FLUX,       KIND_SENSIBLE_HEAT_FLUX
!TOWER_NETC_ECO_EXCHANGE,        KIND_NET_CARBON_PRODUCTION
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

use        types_mod, only : r4, r8, MISSING_R8, PI, deg2rad
use     location_mod, only : location_type, get_location
use time_manager_mod, only : time_type
use    utilities_mod, only : register_module, E_ERR, E_MSG, error_handler, &
                             check_namelist_read, find_namelist_in_file,   &
                             nmlfileunit, do_output, do_nml_file, do_nml_term, &
                             nc_check
use typesizes
use netcdf

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

integer,            allocatable, dimension(:) :: ncid
character(len=129), allocatable, dimension(:) :: fname
real(r8),           allocatable, dimension(:) :: lon, lat, time
integer :: nlon, nlat, ntime, ens_size

! namelist items
character(len=256) :: casename = 'clm_tim'
logical            :: debug = .true.

namelist /obs_def_tower_nml/ casename, debug

contains

!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Start of executable routines
!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine initialize_module

! Called once to set values and allocate space

integer :: iunit, io, rc, i
integer :: dimid, varid

! Prevent multiple calls from executing this code more than once.
if (module_initialized) return

module_initialized = .true.

! Log the version of this source file.
call register_module(source, revision, revdate)

! Read the namelist entry.
call find_namelist_in_file("input.nml", "obs_def_tower_nml", iunit)
read(iunit, nml = obs_def_tower_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_def_tower_nml")

! Record the namelist values used for the run ... 
if (do_nml_file()) write(nmlfileunit, nml=obs_def_tower_nml)
if (do_nml_term()) write(     *     , nml=obs_def_tower_nml)

! FIXME: Figure out how many files (i.e. ensemble size)
ens_size = 4

allocate(fname(ens_size),ncid(ens_size))

ENSEMBLE : do i = 1,ens_size
   write(fname(i),'(A,''.clm2_'',I4.4,''.h.nc'')')trim(casename),i
enddo ENSEMBLE

ncid = 0

i = 1
if (debug) write(*,*)'Reading observations from ',trim(fname(i))

! Harvest information from the first observation file.
! FIXME All other files will be opened to make sure they have the same dimensions.

if (ncid(i) == 0) then ! we need to open it
   call nc_check(nf90_open(trim(fname(i)), nf90_nowrite, ncid(i)), &
       'obs_def_tower_mod:initialize_routine','open '//trim(fname(i)))
endif

call nc_check(nf90_inq_dimid(ncid(i), 'lon', dimid), &
       'obs_def_tower_mod:initialize_routine','inq_dimid lon '//trim(fname(i)))
call nc_check(nf90_inquire_dimension(ncid(i), dimid, len=nlon), &
       'obs_def_tower_mod:initialize_routine','inquire_dimension lon '//trim(fname(i)))

call nc_check(nf90_inq_dimid(ncid(i), 'lat', dimid), &
       'obs_def_tower_mod:initialize_routine','inq_dimid lat '//trim(fname(i)))
call nc_check(nf90_inquire_dimension(ncid(i), dimid, len=nlat), &
       'obs_def_tower_mod:initialize_routine','inquire_dimension lat '//trim(fname(i)))

call nc_check(nf90_inq_dimid(ncid(i), 'time', dimid), &
       'obs_def_tower_mod:initialize_routine','inq_dimid time '//trim(fname(i)))
call nc_check(nf90_inquire_dimension(ncid(i), dimid, len=ntime), &
       'obs_def_tower_mod:initialize_routine','inquire_dimension lat '//trim(fname(i)))

allocate(lon(nlon),lat(nlat),time(ntime))

call nc_check(nf90_inq_varid(ncid(i), 'lon', varid), &
       'obs_def_tower_mod:initialize_routine','inq_varid lon '//trim(fname(i)))
call nc_check(nf90_get_var(ncid(i), varid, lon),     &
       'obs_def_tower_mod:initialize_routine', 'get_var')

call nc_check(nf90_inq_varid(ncid(i), 'lat', varid), &
       'obs_def_tower_mod:initialize_routine','inq_varid lat '//trim(fname(i)))
call nc_check(nf90_get_var(ncid(i), varid, lat),     &
       'obs_def_tower_mod:initialize_routine', 'get_var')

call nc_check(nf90_inq_varid(ncid(i), 'time', varid), &
       'obs_def_tower_mod:initialize_routine','inq_varid time '//trim(fname(i)))
call nc_check(nf90_get_var(ncid(i), varid, time),     &
       'obs_def_tower_mod:initialize_routine', 'get_var')

if (debug) write(*,*)'obs_def_tower  lon',lon
if (debug) write(*,*)'obs_def_tower  lat',lat
if (debug) write(*,*)'obs_def_tower time',time

! FIXME
! check all other ensemble member history files to make sure metadata is the same.

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

if ( .not. module_initialized ) call initialize_module

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

if ( .not. module_initialized ) call initialize_module

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
!
! float NEE(time, lat, lon) ;
!          NEE:long_name = "net ecosystem exchange of carbon, blah, blah, blah" ;
!          NEE:units = "gC/m^2/s" ;
!          NEE:cell_methods = "time: mean" ;
!          NEE:_FillValue = 1.e+36f ;
!          NEE:missing_value = 1.e+36f ;

real(r8),            intent(in)  :: state(:)
type(time_type),     intent(in)  :: state_time
integer,             intent(in)  :: ens_index
type(location_type), intent(in)  :: location
integer,             intent(in)  :: obs_key
real(r8),            intent(out) :: obs_val
integer,             intent(out) :: istatus

integer,  dimension(NF90_MAX_VAR_DIMS) :: dimids
real(r8), dimension(3) :: loc
integer,  dimension(3) :: ncstart, nccount
integer,  dimension(1) :: loninds, latinds
integer                :: gridloni, gridlatj
integer                :: varid, xtype, ndims, natts, dimlen
integer                :: io1, io2
real(r8)               :: loc_lon, loc_lat
real(r4), dimension(1) :: hyperslab
real(r4)               :: spvalR4
real(r8)               :: scale_factor, add_offset

if ( .not. module_initialized ) call initialize_module

obs_val = MISSING_R8
istatus = 1

if (ens_index > ens_size) then
   write(string1,*)'believed to have ',ens_size,'ensemble members for observation operator.'
   write(string2,*)'asking to use operator for ensemble member ',ens_index
   call error_handler(E_ERR, 'get_expected_net_C_production', &
              string1, source, revision, revdate, text2=string2)
endif

! bombproofing ... make sure the netcdf file is open.

! bombproofing ... make sure the variable is the shape and size we expect

call nc_check(nf90_inq_varid(ncid(ens_index), 'NEE', varid), &
              'get_expected_net_C_production', 'inq_varid NEE ')
call nc_check(nf90_inquire_variable(ncid(ens_index), varid, xtype=xtype, ndims=ndims, &
            dimids=dimids, natts=natts),'get_expected_net_C_production','inquire variable NEE ')

if (ndims /= 3) then
   write(string1,*)'NEE is supposed to have 3 dimensions, it has',ndims
   call error_handler(E_ERR, 'get_expected_net_C_production', &
              string1, source, revision, revdate)
endif

! If the variable is not a NF90_FLOAT, then the assumptions for processing
! the missing_value, _FillValue, etc., may not be correct.
if (xtype /= NF90_FLOAT) then
   write(string1,*)'NEE supposed to be a 32 bit real. xtype = ',NF90_FLOAT,' it is ',xtype 
   call error_handler(E_ERR, 'get_expected_net_C_production', &
              string1, source, revision, revdate)
endif

! Dimension 1 is longitude
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(1), len=dimlen), &
              'get_expected_net_C_production', 'inquire_dimension NEE 1')
if (dimlen /= nlon) then
   write(string1,*)'LON has length',nlon,'NEE has ',dimlen,'longitudes.'
   call error_handler(E_ERR, 'get_expected_net_C_production', &
              string1, source, revision, revdate)
endif

! Dimension 2 is latitude
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(2), len=dimlen), &
              'get_expected_net_C_production', 'inquire_dimension NEE 2')
if (dimlen /= nlat) then
   write(string1,*)'LAT has length',nlat,'NEE has ',dimlen,'latitudes.'
   call error_handler(E_ERR, 'get_expected_net_C_production', &
              string1, source, revision, revdate)
endif

! Dimension 3 is time
call nc_check(nf90_inquire_dimension(ncid(ens_index), dimids(3), len=dimlen), &
              'get_expected_net_C_production', 'inquire_dimension NEE 3')
if (dimlen /= ntime) then
   write(string1,*)'TIME has length',ntime,'NEE has ',dimlen,'times.'
   call error_handler(E_ERR, 'get_expected_net_C_production', &
              string1, source, revision, revdate)
endif

! Find the grid cell of interest 
! Get the individual locations values

loc      = get_location(location)       ! loc is in DEGREES
loc_lon  = loc(1)
loc_lat  = loc(2)
latinds  = minloc(abs(lat - loc_lat))   ! these return 'arrays' ...
loninds  = minloc(abs(lon - loc_lon))   ! these return 'arrays' ...
gridlatj = latinds(1)
gridloni = loninds(1)

if (debug .and. do_output()) then
   write(*,*)'get_expected_net_C_production:targetlon, lon, lon index is ', &
                                           loc_lon,lon(gridloni),gridloni
   write(*,*)'get_expected_net_C_production:targetlat, lat, lat index is ', &
                                           loc_lat,lat(gridlatj),gridlatj
endif

! For now, just grab the last timestep ... 

ncstart = (/ gridloni, gridlatj, ntime /)
nccount = (/        1,        1,     1 /)

call nc_check(nf90_get_var(ncid(ens_index), varid, hyperslab, start=ncstart, count=nccount), &
     'get_expected_net_C_production', 'get_var')

obs_val = hyperslab(1)

! Apply any netCDF attributes ...

io1 = nf90_get_att(ncid(ens_index), varid, '_FillValue' , spvalR4)
if ((io1 == NF90_NOERR) .and. (hyperslab(1) == spvalR4)) obs_val = MISSING_R8

io2 = nf90_get_att(ncid(ens_index), varid, 'missing_value' , spvalR4)
if ((io2 == NF90_NOERR) .and. (hyperslab(1) == spvalR4)) obs_val = MISSING_R8

io1 = nf90_get_att(ncid(ens_index), varid, 'scale_factor', scale_factor)
io2 = nf90_get_att(ncid(ens_index), varid, 'add_offset'  , add_offset)

if ( (io1 == NF90_NOERR) .and. (io2 == NF90_NOERR) ) then
   if (obs_val /= MISSING_R8) obs_val = obs_val * scale_factor + add_offset
elseif (io1 == NF90_NOERR) then
   if (obs_val /= MISSING_R8) obs_val = obs_val * scale_factor
elseif (io2 == NF90_NOERR) then
   if (obs_val /= MISSING_R8) obs_val = obs_val + add_offset
endif

if (obs_val /= MISSING_R8) istatus = 0

end subroutine get_expected_net_C_production



end module obs_def_tower_mod

! END DART PREPROCESS MODULE CODE
!-----------------------------------------------------------------------------
