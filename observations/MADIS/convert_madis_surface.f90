!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_madis_surface - program that reads a MADIS netCDF land 
!                           surface observation file and writes a text 
!                           file of observations within the analysis 
!                           time.  The text file can be used in other 
!                           programs that write obs_seq.out files.
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program convert_madis_surface

use        types_mod, only : r8
use time_manager_mod, only : time_type, set_calendar_type, set_date, &
                             increment_time, get_time, GREGORIAN, operator(-)
use     location_mod, only : VERTISSURFACE
use obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                             static_init_obs_sequence, init_obs, write_obs_seq, & 
                             append_obs_to_seq, init_obs_sequence, get_num_obs, & 
                             set_copy_meta_data, set_qc_meta_data
use       meteor_mod, only : sat_vapor_pressure, specific_humidity, & 
                             wind_dirspd_to_uv, invert_altimeter
use ncep_obs_err_mod, only : ncep_land_temp_error, ncep_land_wind_error, &
                             ncep_land_altim_error, ncep_land_moist_error
use     obs_kind_mod, only : LAND_SFC_U_WIND_COMPONENT, LAND_SFC_V_WIND_COMPONENT, &
                             LAND_SFC_TEMPERATURE, LAND_SFC_SPECIFIC_HUMIDITY, & 
                             LAND_SFC_ALTIMETER                
use           netcdf

implicit none

character(len=16),  parameter :: surface_netcdf_file = 'surface_input.nc'
character(len=129), parameter :: surface_out_file    = 'obs_seq.land_sfc'

integer, parameter   :: dsecobs    = 420, &   ! observation window
                        num_copies = 1,   &   ! number of copies in sequence
                        num_qc     = 1        ! number of QC entries

character (len=129) :: meta_data
character (len=80)  :: name
character (len=19)  :: datestr
integer :: rcode, ncid, varid, nobs, n, i, oday, osec, dday, &
           dsec, nused, iyear, imonth, iday, ihour, imin, isec
logical :: file_exist
real(r8) :: alti_miss, tair_miss, tdew_miss, wdir_miss, wspd_miss, uwnd, &
            vwnd, qobs, qsat, oerr, stn_pres, qc

integer,  allocatable :: tobs(:)
real(r8), allocatable :: lat(:), lon(:), elev(:), alti(:), tair(:), & 
                         tdew(:), wdir(:), wspd(:), latu(:), lonu(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs
type(time_type)         :: comp_day0, time_obs, time_anal

print*,'Enter the analysis time (yyyy-mm-dd_hh:mm:ss)'
read*,datestr

! put the analysis date into DART format
call set_calendar_type(GREGORIAN)
read(datestr(1:4),   fmt='(i4)') iyear
read(datestr(6:7),   fmt='(i2)') imonth
read(datestr(9:10),  fmt='(i2)') iday
read(datestr(12:13), fmt='(i2)') ihour
read(datestr(15:16), fmt='(i2)') imin
read(datestr(18:19), fmt='(i2)') isec
time_anal = set_date(iyear, imonth, iday, ihour, imin, isec)
call get_time(time_anal, osec, oday)
comp_day0 = set_date(1970, 1, 1, 0, 0, 0)

rcode = nf90_open(surface_netcdf_file, nf90_nowrite, ncid)
call check( nf90_inq_dimid(ncid, "recNum", varid) )
call check( nf90_inquire_dimension(ncid, varid, name, nobs) )

allocate( lat(nobs))  ;  allocate( lon(nobs))
allocate(latu(nobs))  ;  allocate(lonu(nobs))
allocate(elev(nobs))  ;  allocate(alti(nobs))
allocate(tair(nobs))  ;  allocate(tdew(nobs))
allocate(wdir(nobs))  ;  allocate(wspd(nobs))
allocate(tobs(nobs))

! read the latitude array
call check( nf90_inq_varid(ncid, "latitude", varid) )
call check( nf90_get_var(ncid, varid, lat) )

! read the latitude array
call check( nf90_inq_varid(ncid, "longitude", varid) )
call check( nf90_get_var(ncid, varid, lon) )

! read the elevation array
call check( nf90_inq_varid(ncid, "elevation", varid) )
call check( nf90_get_var(ncid, varid, elev) )

! read the altimeter setting array
call check( nf90_inq_varid(ncid, "altimeter", varid) )
call check( nf90_get_var(ncid, varid, alti) )
call check( nf90_get_att(ncid, varid, '_FillValue', alti_miss) )

! read the air temperature array
call check( nf90_inq_varid(ncid, "temperature", varid) )
call check( nf90_get_var(ncid, varid, tair) )
call check( nf90_get_att(ncid, varid, '_FillValue', tair_miss) )

! read the dew-point temperature array
call check( nf90_inq_varid(ncid, "dewpoint", varid) )
call check( nf90_get_var(ncid, varid, tdew) )
call check( nf90_get_att(ncid, varid, '_FillValue', tdew_miss) )

! read the wind direction array
call check( nf90_inq_varid(ncid, "windDir", varid) )
call check( nf90_get_var(ncid, varid, wdir) )
call check( nf90_get_att(ncid, varid, '_FillValue', wdir_miss) )

! read the wind speed array
call check( nf90_inq_varid(ncid, "windSpeed", varid) )
call check( nf90_get_var(ncid, varid, wspd) )
call check( nf90_get_att(ncid, varid, '_FillValue', wspd_miss) )

! read the observation time array
call check( nf90_inq_varid(ncid, "timeObs", varid) )
call check( nf90_get_var(ncid, varid, tobs) )

call check( nf90_close(ncid) )

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs, num_copies, num_qc)
inquire(file=surface_out_file, exist=file_exist)
if ( file_exist ) then

  call read_obs_seq(surface_out_file, 0, 0, 5*nobs, obs_seq)

else

  call init_obs_sequence(obs_seq, num_copies, num_qc, 5*nobs)
  do i = 1, num_copies
    meta_data = 'NCEP BUFR observation'
    call set_copy_meta_data(obs_seq, i, meta_data)
  end do
  do i = 1, num_qc
    meta_data = 'NCEP QC index'
    call set_qc_meta_data(obs_seq, i, meta_data)
  end do

end if

nused = 0
obsloop: do n = 1, nobs

  ! determine whether observation is close to analysis time
  time_obs = increment_time(comp_day0, mod(tobs(n),86400), tobs(n) / 86400)
  call get_time((time_anal - time_obs), dsec, dday)
  if ( dsec > dsecobs .or. dday > 0 ) cycle obsloop
  if ( lon(n) < 0.0_r8 )  lon(n) = lon(n) + 360.0_r8

  ! check to make sure this observation has not been used
  do i = 1, nused
    if ( lon(n) == lon(i) .and. lat(n) == lat(i) ) cycle obsloop
  end do
  qc = 1

  ! add altimeter data to text file
  if ( alti(n) /= alti_miss ) then

    call create_obs_type(lat(n), lon(n), elev(n), VERTISSURFACE, alti(n) * 0.01_r8, & 
                         LAND_SFC_ALTIMETER, ncep_land_altim_error, oday, osec, qc, obs)
    call append_obs_to_seq(obs_seq, obs)

  end if

  ! add wind component data to text file
  if ( wdir(n) .ne. wdir_miss .and. wspd(n) .ne. wspd_miss ) then

    call wind_dirspd_to_uv(wdir(n), wspd(n), uwnd, vwnd)
    if ( abs(uwnd) < 150.0_r8 .and. abs(vwnd) < 150.0_r8 ) then

      call create_obs_type(lat(n), lon(n), elev(n), VERTISSURFACE, uwnd, &
                           LAND_SFC_U_WIND_COMPONENT, ncep_land_wind_error, & 
                           oday, osec, qc, obs)
      call append_obs_to_seq(obs_seq, obs)
      call create_obs_type(lat(n), lon(n), elev(n), VERTISSURFACE, vwnd, &
                           LAND_SFC_V_WIND_COMPONENT, ncep_land_wind_error, &
                           oday, osec, qc, obs)
      call append_obs_to_seq(obs_seq, obs)

    end if

  end if

  ! add air temperature data to text file
  if ( tair(n) .ne. tair_miss ) then 

    call create_obs_type(lat(n), lon(n), elev(n), VERTISSURFACE, tair(n), &
                         LAND_SFC_TEMPERATURE, ncep_land_temp_error, oday, osec, qc, obs)
    call append_obs_to_seq(obs_seq, obs)

  end if

  ! add dew-point temperature data to text file, but as specific humidity
  if ( tair(n) .ne. tair_miss .and. tdew(n) .ne. tdew_miss .and. alti(n) .ne. alti_miss ) then

    stn_pres = invert_altimeter(alti(n) * 0.01_r8, elev(n)) * 100.0_r8
    qobs = specific_humidity(sat_vapor_pressure(tdew(n)), stn_pres)
    qsat = specific_humidity(sat_vapor_pressure(tair(n)), stn_pres)
    oerr = max(ncep_land_moist_error * qsat, 0.0001_r8)

    if ( abs(qobs) < 100.0_r8 ) then
      call create_obs_type(lat(n), lon(n), elev(n), VERTISSURFACE, qobs, &
                           LAND_SFC_SPECIFIC_HUMIDITY, oerr, oday, osec, qc, obs)
      call append_obs_to_seq(obs_seq, obs)
    end if

  end if

  nused = nused + 1
  latu(nused) = lat(n)
  lonu(nused) = lon(n)

end do obsloop

if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, surface_out_file)

stop
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   check - subroutine that checks the flag from a netCDF function.  If 
!           there is an error, the message is displayed.
!
!    istatus - netCDF output flag
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check( istatus ) 

use netcdf

implicit none

integer, intent (in) :: istatus

if(istatus /= nf90_noerr) print*,'Netcdf error: ',trim(nf90_strerror(istatus))

end subroutine check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   create_obs_type - subroutine that is used to create an observation
!                     type from observation data.
!
!    lat   - latitude of observation
!    lon   - longitude of observation
!    pres  - pressure of observation
!    vcord - vertical coordinate
!    obsv  - observation value
!    okind - observation kind
!    oerr  - observation error
!    day   - gregorian day
!    sec   - gregorian second
!    qc    - quality control value
!    obs   - observation type
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_obs_type(lat, lon, pres, vcord, obsv, okind, oerr, day, sec, qc, obs)

use types_mod,        only : r8
use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location
use     location_mod, only : location_type, set_location
use time_manager_mod, only : time_type, set_time

implicit none

integer, intent(in)         :: okind, vcord, day, sec
real(r8), intent(in)        :: lat, lon, pres, obsv, oerr, qc
type(obs_type), intent(inout) :: obs

real(r8)              :: obs_val(1), qc_val(1)
type(obs_def_type)    :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, pres, vcord))
call set_obs_def_kind(obs_def, okind)
call set_obs_def_time(obs_def, set_time(sec, day))
call set_obs_def_error_variance(obs_def, oerr * oerr)
call set_obs_def(obs, obs_def)

obs_val(1) = obsv
call set_obs_values(obs, obs_val)
qc_val(1)  = qc
call set_qc(obs, qc_val)

return
end subroutine create_obs_type
