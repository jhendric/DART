!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_madis_marine - program that reads a MADIS netCDF marine
!                          surface observation file and writes a text
!                          file of observations within the analysis
!                          time.  The text file can be used in other
!                          programs that write obs_seq.out files.
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program convert_madis_marine

use             types_mod, only : r8
use      time_manager_mod, only : time_type, set_calendar_type, set_date, &
                                  increment_time, get_time, operator(-), GREGORIAN
use          location_mod, only : VERTISSURFACE
use      obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                                  static_init_obs_sequence, init_obs, write_obs_seq, &
                                  append_obs_to_seq, init_obs_sequence, get_num_obs, &
                                  set_copy_meta_data, set_qc_meta_data
use            meteor_mod, only : sat_vapor_pressure, specific_humidity, & 
                                  wind_dirspd_to_uv
use      ncep_obs_err_mod, only : ncep_marine_temp_error, ncep_marine_moist_error, &
                                  ncep_marine_wind_error, ncep_marine_altim_error
use          obs_kind_mod, only : MARINE_SFC_U_WIND_COMPONENT, MARINE_SFC_V_WIND_COMPONENT, &
                                  MARINE_SFC_TEMPERATURE, MARINE_SFC_SPECIFIC_HUMIDITY, &
                                  MARINE_SFC_ALTIMETER
use obs_def_altimeter_mod, only : compute_altimeter
use                netcdf

implicit none

character(len=15),  parameter :: marine_netcdf_file = 'marine_input.nc'
character(len=129), parameter :: marine_out_file    = 'obs_seq.marine_sfc'

integer, parameter   :: dsecobs = 2700,   &   ! observation window
                        num_copies = 1,   &   ! number of copies in sequence
                        num_qc     = 1        ! number of QC entries

real(r8), parameter :: def_elev = 0.0_r8

character (len=129) :: meta_data
character (len=80)  :: name
character (len=19)  :: datestr
integer :: rcode, ncid, varid, nobs, n, i, dday, dsec, oday, &
           osec, nused, iyear, imonth, iday, ihour, imin, isec
logical :: file_exist
real(r8) :: sfcp_miss, tair_miss, tdew_miss, wdir_miss, wspd_miss, uwnd, &
            vwnd, altim, qobs, qsat, qobserr, slp_miss, elev_miss, qc

integer,  allocatable :: tobs(:)
real(r8), allocatable :: lat(:), lon(:), elev(:), sfcp(:), tair(:), slp(:), & 
                         tdew(:), wdir(:), wspd(:), latu(:), lonu(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs
type(time_type)         :: comp_day0, time_anal, time_obs

print*,'Enter the observation date (yyyy-mm-dd_hh:mm:ss)'
read*,datestr

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

rcode = nf90_open(marine_netcdf_file, nf90_nowrite, ncid)

call check( nf90_inq_dimid(ncid, "recNum", varid) )
call check( nf90_inquire_dimension(ncid, varid, name, nobs) )

allocate( lat(nobs))  ;  allocate( lon(nobs))
allocate(latu(nobs))  ;  allocate(lonu(nobs))
allocate(elev(nobs))  ;  allocate(sfcp(nobs))
allocate(tair(nobs))  ;  allocate(tdew(nobs))
allocate(wdir(nobs))  ;  allocate(wspd(nobs))
allocate(tobs(nobs))  ;  allocate(slp(nobs))

! read the latitude array
call check( nf90_inq_varid(ncid, "latitude", varid) )
call check( nf90_get_var(ncid, varid, lat) )

! read the longitude array
call check( nf90_inq_varid(ncid, "longitude", varid) )
call check( nf90_get_var(ncid, varid, lon) )

! read the elevation array
call check( nf90_inq_varid(ncid, "elevation", varid) )
call check( nf90_get_var(ncid, varid, elev) )
call check( nf90_get_att(ncid, varid, '_FillValue', elev_miss) )

! read the altimeter setting array
call check( nf90_inq_varid(ncid, "stationPress", varid) )
call check( nf90_get_var(ncid, varid, sfcp) )
call check( nf90_get_att(ncid, varid, '_FillValue', sfcp_miss) )

! read the sea-level pressure array
call check( nf90_inq_varid(ncid, "seaLevelPress", varid) )
call check( nf90_get_var(ncid, varid, slp) )
call check( nf90_get_att(ncid, varid, '_FillValue', slp_miss) )

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
inquire(file=marine_out_file, exist=file_exist)
if ( file_exist ) then

  call read_obs_seq(marine_out_file, 0, 0, 5*nobs, obs_seq)

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

  ! check the lat/lon values to see if they are ok
  if ( lat(n) >  90.0_r8 .or. lat(n) <  -90.0_r8 ) cycle obsloop
  if ( lon(n) > 180.0_r8 .or. lon(n) < -180.0_r8 ) cycle obsloop

  ! determine whether observation is close to the analysis time
  time_obs = increment_time(comp_day0, mod(tobs(n),86400), tobs(n) / 86400)
  call get_time((time_anal - time_obs), dsec, dday)
  if ( dsec > dsecobs .or. dday > 0 ) cycle obsloop
  if ( lon(n) < 0.0_r8 )  lon(n) = lon(n) + 360.0_r8

  do i = 1, nused
    if ( lon(n) == lon(i) .and. lat(n) == lat(i) ) cycle obsloop
  end do
  qc = 1.0_r8

  ! add altimeter data to obs_seq
  if ( sfcp(n) /= sfcp_miss .and. elev(n) /= elev_miss ) then

    altim = compute_altimeter(sfcp(n) * 0.01_r8, elev(n))
    call create_obs_type(lat(n), lon(n), elev(n), VERTISSURFACE, altim, &
                         MARINE_SFC_ALTIMETER, ncep_marine_altim_error, &
                         oday, osec, qc, obs)
    call append_obs_to_seq(obs_seq, obs)

  !  if surface pressure and elevation do not exist, use SLP.
  else if ( slp(n) /= slp_miss ) then

    altim = compute_altimeter(slp(n) * 0.01_r8, 0.0_r8) 
    call create_obs_type(lat(n), lon(n), def_elev, VERTISSURFACE, altim, &
                         MARINE_SFC_ALTIMETER, ncep_marine_altim_error, &
                         oday, osec, qc, obs)
    call append_obs_to_seq(obs_seq, obs)

  end if
  if ( elev(n) == elev_miss )  elev(n) = def_elev

  ! add wind component data to obs. sequence
  if ( wdir(n) .ne. wdir_miss .and. wspd(n) .ne. wspd_miss ) then

    call wind_dirspd_to_uv(wdir(n), wspd(n), uwnd, vwnd)
    if ( abs(uwnd) < 150.0_r8 .and. abs(vwnd) < 150.0_r8 ) then

      call create_obs_type(lat(n), lon(n), elev(n), VERTISSURFACE, uwnd, &
                           MARINE_SFC_U_WIND_COMPONENT, ncep_marine_wind_error, &
                           oday, osec, qc, obs)
      call append_obs_to_seq(obs_seq, obs)
      call create_obs_type(lat(n), lon(n), elev(n), VERTISSURFACE, vwnd, &
                           MARINE_SFC_V_WIND_COMPONENT, ncep_marine_wind_error, &
                           oday, osec, qc, obs)
      call append_obs_to_seq(obs_seq, obs)

    end if

  end if

  ! add air temperature data to obs. sequence
  if ( tair(n) .ne. tair_miss ) then 
    
    call create_obs_type(lat(n), lon(n), def_elev, VERTISSURFACE, tair(n), &
                         MARINE_SFC_TEMPERATURE, ncep_marine_temp_error, &
                         oday, osec, qc, obs)
    call append_obs_to_seq(obs_seq, obs)
 
  end if

  ! add dew-point temperature data to obs. sequence, but as specific humidity
  if ( tair(n) .ne. tair_miss .and. tdew(n) .ne. tdew_miss .and. sfcp(n) .ne. sfcp_miss ) then

    qobs    = specific_humidity(sat_vapor_pressure(tdew(n)), sfcp(n))
    qsat    = specific_humidity(sat_vapor_pressure(tair(n)), sfcp(n))
    qobserr = max(ncep_marine_moist_error * qsat, 0.0001_r8)

    if ( abs(qobs) < 100.0_r8 ) then

      call create_obs_type(lat(n), lon(n), def_elev, VERTISSURFACE, qobs, &
                           MARINE_SFC_SPECIFIC_HUMIDITY, qobserr, oday, osec, qc, obs)
      call append_obs_to_seq(obs_seq, obs)

    end if

  end if

  nused = nused + 1
  latu(nused) = lat(n)
  lonu(nused) = lon(n)

end do obsloop

if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, marine_out_file)

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

if(istatus /= nf90_noerr) then 
  print*,'Netcdf error: ',trim(nf90_strerror(istatus))
  stop
end if

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
