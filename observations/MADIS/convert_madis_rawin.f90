! DART software - Copyright � 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program convert_madis_rawin

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   convert_madis_rawin - program that reads a netCDF file from the 
!                         MADIS database that contains rawinsonde data 
!                         and writes a DART obs_seq file using the DART 
!                         library routines.
!
!     created Dec. 2007 Ryan Torn, NCAR/MMM
!     modified Dec. 2008 Soyoung Ha and David Dowell, NCAR/MMM
!     - added dewpoint as an output variable
!     - added relative humidity as an output variable
!
!     modified to use a common set of utilities, better netcdf error checks,
!     able to insert obs with any time correctly (not only monotonically
!     increasing times)    nancy collins,  ncar/image   11 march 2010
!     
!     keep original obs times, make source for all converters as similar
!     as possbile.   nancy collins,  ncar/image   26 march 2010
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, missing_r8
use     utilities_mod, only : nc_check, initialize_utilities, finalize_utilities
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                                  increment_time, get_time, operator(-), GREGORIAN
use      location_mod, only : VERTISSURFACE, VERTISPRESSURE, VERTISHEIGHT
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data
use        meteor_mod, only : sat_vapor_pressure, specific_humidity, & 
                              wind_dirspd_to_uv, pres_alt_to_pres, &
                              temp_and_dewpoint_to_rh
use       obs_err_mod, only : rawin_temp_error, rawin_wind_error, &
                              rawin_pres_error, rawin_rel_hum_error
use dewpoint_obs_err_mod, only : dewpt_error_from_rh_and_temp, &
                                 rh_error_from_dewpt_and_temp
use obs_def_altimeter_mod, only : compute_altimeter
use          obs_kind_mod, only : RADIOSONDE_U_WIND_COMPONENT,      & 
                                  RADIOSONDE_V_WIND_COMPONENT,      & 
                                  RADIOSONDE_TEMPERATURE,           & 
                                  RADIOSONDE_SPECIFIC_HUMIDITY,     &
                                  RADIOSONDE_RELATIVE_HUMIDITY,     &
                                  RADIOSONDE_DEWPOINT,              &
                                  RADIOSONDE_SURFACE_ALTIMETER 
use     obs_utilities_mod, only : add_obs_to_seq, create_3d_obs, &
                                  getdimlen, getvar_int, set_missing_name

use           netcdf

implicit none

character(len=19),  parameter :: rawin_in_file  = 'rawin_input.nc'
character(len=129), parameter :: rawin_out_file = 'obs_seq.rawin'

! the following logical parameters control which water-vapor variables appear in the output file,
! whether to use the NCEP error or Lin and Hubbard (2004) moisture error model, and if the
! input file has data quality control fields, whether to use or ignore them.
logical, parameter :: LH_err                    = .false.
logical, parameter :: include_specific_humidity = .true.
logical, parameter :: include_relative_humidity = .false.
logical, parameter :: include_dewpoint          = .false.
logical, parameter :: use_input_qc              = .true. 

integer, parameter :: num_copies = 1,   &   ! number of copies in sequence
                      num_qc     = 1        ! number of QC entries

character(len=129) :: meta_data

integer :: oday, osec, nman, nsig, nsound, nused, & 
           nmaxml, nmaxsw, nmaxst, maxobs, nvars_man, nvars_sigt, k, n, i, ncid

integer, allocatable :: obscnt(:)

logical :: fexist, sigwnd, sigtmp, first_obs

real(r8) :: otime, lat, lon, elev, uwnd, vwnd, qobs, qsat, dptk, oerr, &
            pres_miss, wdir_miss, wspd_miss, tair_miss, tdew_miss, prespa, & 
            qc, altim, rh, qerr  ! , time_miss
real(r8), allocatable :: latu(:), lonu(:)

real(r8), allocatable :: pres(:), wdir(:), wspd(:), tair(:), tdew(:)
integer,  allocatable :: qc_pres(:), qc_wdir(:), qc_wspd(:), qc_tair(:), qc_tdew(:)

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: comp_day0, time_obs, prev_time


call initialize_utilities('convert_madis_rawin')

! prompt for optional significant level values
print*,'Include significant level winds, temperature?: '
read*, sigwnd, sigtmp

! shouldn't need to do this now - the nf90 open will catch it below
inquire(file=trim(rawin_in_file), exist=fexist) ! check for drop file 
if ( .NOT. fexist ) then
  print*,'Rawinsonde file ',rawin_in_file,' does not exist, exiting'
  stop
endif

! put the reference date into DART format
call set_calendar_type(GREGORIAN)
comp_day0 = set_date(1970, 1, 1, 0, 0, 0)

first_obs = .true.


call nc_check( nf90_open(rawin_in_file, nf90_nowrite, ncid), &
               'convert_madis_rawin', 'opening file '//trim(rawin_in_file))

call getdimlen(ncid, "recNum", nsound)
call set_missing_name("missing_value")

allocate(obscnt(nsound))  
allocate(latu(nsound))  ; allocate(lonu(nsound))
nmaxml = 0  ;  nmaxsw = 0  ;  nmaxst = 0

call getvar_int(ncid, "numMand", obscnt)
do n = 1, nsound 
  if ( obscnt(n) > nmaxml .and. obscnt(n) < 25 )  nmaxml = obscnt(n)
end do

if ( sigwnd ) then
  call getvar_int(ncid, "numSigW", obscnt)
  do n = 1, nsound
    if ( obscnt(n) > nmaxsw .and. obscnt(n) < 150 )  nmaxsw = obscnt(n)
  end do
endif

if ( sigtmp ) then
  call getvar_int(ncid, "numSigT", obscnt)
  do n = 1, nsound
    if ( obscnt(n) > nmaxst .and. obscnt(n) < 150 )  nmaxst = obscnt(n)
  end do
endif

nvars_man = 4
nvars_sigt = 1
if (include_specific_humidity) then
  nvars_man  = nvars_man  + 1
  nvars_sigt = nvars_sigt + 1
endif
if (include_relative_humidity) then
  nvars_man  = nvars_man  + 1
  nvars_sigt = nvars_sigt + 1
endif
if (include_dewpoint) then
  nvars_man  = nvars_man  + 1
  nvars_sigt = nvars_sigt + 1
endif

maxobs = nsound * (nvars_man * nmaxml + 2 * nmaxsw + nvars_sigt * nmaxst + 1)
deallocate(obscnt)

!  either read existing obs_seq or create a new one
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

inquire(file=rawin_out_file, exist=fexist)
if ( fexist ) then

  ! existing file found, append to it
  call read_obs_seq(rawin_out_file, 0, 0, maxobs, obs_seq)

else

  ! create a new one
  call init_obs_sequence(obs_seq, num_copies, num_qc, maxobs)
  do n = 1, num_copies
    meta_data = 'MADIS observation'
    call set_copy_meta_data(obs_seq, n, meta_data)
  end do
  do n = 1, num_qc
    meta_data = 'Data QC'
    call set_qc_meta_data(obs_seq, n, meta_data)
  end do

endif

! Set the DART data quality control.  Be consistent with NCEP codes;
! 0 is 'must use', 1 is good, no reason not to use it.
qc = 1.0_r8

nused = 0
sondeloop : do n = 1, nsound !  loop over all soundings in the file 

  call getvar_real_1d_1val(ncid, "staLat",  n, lat  )
  call getvar_real_1d_1val(ncid, "staLon",  n, lon  )
  call getvar_real_1d_1val(ncid, "staElev", n, elev )
  call getvar_int_1d_1val (ncid, "numMand", n, nman )
  call getvar_real_1d_1val(ncid, "synTime", n, otime)  !! , time_miss)
  ! the original code had a line to get the fill value here but it
  ! was commented out.  is there one?  do we need it?

  if (nman <= 0 .or. nman > nmaxml) cycle sondeloop

  if ( otime < 0.0_r8 ) cycle sondeloop

  ! compute time of observation
  time_obs = increment_time(comp_day0, nint(otime))

  ! check the lat/lon values to see if they are ok
  if ( lat >  90.0_r8 .or. lat <  -90.0_r8 ) cycle sondeloop
  if ( lon > 180.0_r8 .or. lon < -180.0_r8 ) cycle sondeloop

  ! change lon from -180 to 180 into 0-360
  if (lon < 0.0_r8) lon = lon + 360.0_r8

  ! Check for duplicate observations
  do i = 1, nused
    if ( lon == lonu(i) .and. &
         lat == latu(i) ) cycle sondeloop
  end do


  ! extract actual time of observation in file into oday, osec.
  call get_time(time_obs, osec, oday)

  allocate(pres(nman))  ;  allocate(tair(nman))  ;  allocate(tdew(nman))
  allocate(wdir(nman))  ;  allocate(wspd(nman))

  allocate(qc_pres(nman))  ;  allocate(qc_tair(nman))  ;  allocate(qc_tdew(nman))
  allocate(qc_wdir(nman))  ;  allocate(qc_wspd(nman))

  call getvar_real_2d(ncid, "prMan", n, nman, pres, pres_miss)
  call getvar_real_2d(ncid, "tpMan", n, nman, tair, tair_miss)
  call getvar_real_2d(ncid, "tdMan", n, nman, tdew, tdew_miss)
  call getvar_real_2d(ncid, "wdMan", n, nman, wdir, wdir_miss)
  call getvar_real_2d(ncid, "wsMan", n, nman, wspd, wspd_miss)

  ! if user says to use QC, read them in or fill if not there
  if (use_input_qc) then
     call get_or_fill_QC_2d(ncid, "prManQCR", n, nman, qc_pres)
     call get_or_fill_QC_2d(ncid, "tpManQCR", n, nman, qc_tair)
     call get_or_fill_QC_2d(ncid, "tdManQCR", n, nman, qc_tdew)
     call get_or_fill_QC_2d(ncid, "wdManQCR", n, nman, qc_wdir)
     call get_or_fill_QC_2d(ncid, "wsManQCR", n, nman, qc_wspd)
  else
     qc_pres = 0
     qc_tair = 0 ;  qc_tdew = 0
     qc_wdir = 0 ;  qc_wspd = 0
  endif

  if ( pres(1) /= pres_miss .and. qc_pres(1) == 0 ) then

    altim = compute_altimeter(pres(1), elev)
    oerr  = rawin_pres_error(pres_alt_to_pres(elev) * 0.01_r8)
    if ( altim >= 890.0_r8 .and. altim <= 1100.0_r8 .and. oerr /= missing_r8 ) then

      call create_3d_obs(lat, lon, elev, VERTISSURFACE, altim, &
                         RADIOSONDE_SURFACE_ALTIMETER, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

    endif

  endif

  do k = 2, nman   ! obtain the mandatory level data

    prespa = pres(k) * 100.0_r8

    if ( wdir(k) /= wdir_miss .and. qc_wdir(k) == 0 .and. &
         wspd(k) /= wspd_miss .and. qc_wspd(k) == 0  ) then

      call wind_dirspd_to_uv(wdir(k), wspd(k), uwnd, vwnd)
      oerr = rawin_wind_error(pres(k))
      if ( abs(uwnd) <= 150.0_r8 .and. & 
           abs(vwnd) <= 150.0_r8 .and. oerr /= missing_r8 ) then

        call create_3d_obs(lat, lon, prespa, VERTISPRESSURE, uwnd, &
                           RADIOSONDE_U_WIND_COMPONENT, oerr, oday, osec, qc, obs)
        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
 
        call create_3d_obs(lat, lon, prespa, VERTISPRESSURE, vwnd, &
                           RADIOSONDE_V_WIND_COMPONENT, oerr, oday, osec, qc, obs)
        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

      endif

    endif

    if ( tair(k) /= tair_miss .and. qc_tair(k) == 0 ) then

      oerr = rawin_temp_error(pres(k))
      if ( tair(k) >= 180.0_r8 .and. &
           tair(k) <= 330.0_r8 .and. oerr /= missing_r8 ) then

        call create_3d_obs(lat, lon, prespa, VERTISPRESSURE, tair(k), &
                           RADIOSONDE_TEMPERATURE, oerr, oday, osec, qc, obs)
        call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

      endif

    endif
  
    ! if the air and dewpoint obs are both ok, then see which of the possible
    ! three types of moisture obs to generate.
    if ( tair(k) /= tair_miss .and. qc_tair(k) == 0 .and. &
         tdew(k) /= tdew_miss .and. qc_tdew(k) == 0  ) then

      ! tdew is the dewpoint depression
      dptk = tair(k) - tdew(k)

      if ( include_specific_humidity ) then

        qobs = specific_humidity(sat_vapor_pressure(dptk),    prespa)
        qsat = specific_humidity(sat_vapor_pressure(tair(k)), prespa)
        if ( LH_err ) then
          qerr = rh_error_from_dewpt_and_temp(tair(k), dptk)
        else
          qerr = rawin_rel_hum_error(pres(k), tair(k), qobs / qsat)
        endif
        oerr = max(qerr * qsat, 0.0001_r8)
  
        if ( qobs >  0.0_r8  .and. &
             qobs <= 0.07_r8 .and. qerr /= missing_r8 ) then
  
          call create_3d_obs(lat, lon, prespa, VERTISPRESSURE, qobs, &
                             RADIOSONDE_SPECIFIC_HUMIDITY, oerr, oday, osec, qc, obs)
          call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
        endif
  
      endif
  
      if ( include_relative_humidity ) then
  
        rh = temp_and_dewpoint_to_rh(tair(k), dptk)
        if ( LH_err ) then
          oerr = rh_error_from_dewpt_and_temp(tair(k), dptk)
        else
          oerr = rawin_rel_hum_error(pres(k), tair(k), rh)
        endif
  
        if ( rh >  0.0_r8 .and. &
             rh <= 1.5_r8 .and. oerr /= missing_r8 ) then

          call create_3d_obs(lat, lon, prespa, VERTISPRESSURE, rh, &
                             RADIOSONDE_RELATIVE_HUMIDITY, oerr, oday, osec, qc, obs)
          call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
        endif

      endif
  
      if ( include_dewpoint ) then

        rh = temp_and_dewpoint_to_rh(tair(k), dptk)
        oerr = dewpt_error_from_rh_and_temp(tair(k), rh)
  
        if ( rh >  0.0_r8 .and. &
             rh <= 1.5_r8 .and. oerr /= missing_r8 ) then

          call create_3d_obs(lat, lon, prespa, VERTISPRESSURE, dptk, &
                             RADIOSONDE_DEWPOINT, oerr, oday, osec, qc, obs)
          call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

        endif

      endif

    endif  ! quality control/missing check on tair, tdew

  end do
  deallocate(pres, wdir, wspd, tair, tdew, qc_pres, qc_wdir, qc_wspd, qc_tair, qc_tdew)

  !  If desired, read the significant-level temperature data, write to obs_seq
  call getvar_int_1d_1val(ncid, "numSigT", n, nsig )

  if ( sigtmp .and. nsig <= nmaxst ) then

    allocate(pres(nsig))     ;  allocate(tair(nsig))     ;  allocate(tdew(nsig))
    allocate(qc_pres(nsig))  ;  allocate(qc_tair(nsig))  ;  allocate(qc_tdew(nsig))

    !  read significant level data
    call getvar_real_2d(ncid, "prSigT", n, nsig, pres, pres_miss)
    call getvar_real_2d(ncid, "tpSigT", n, nsig, tair, tair_miss)
    call getvar_real_2d(ncid, "tdSigT", n, nsig, tdew, tdew_miss)

    if (use_input_qc) then
       call get_or_fill_QC_2d(ncid, "prSigTQCR", n, nsig, qc_pres)
       call get_or_fill_QC_2d(ncid, "tpSigTQCR", n, nsig, qc_tair)
       call get_or_fill_QC_2d(ncid, "tdSigTQCR", n, nsig, qc_tdew)
    else
       qc_pres = 0
       qc_tair = 0
       qc_tdew = 0
    endif
  
    do k = 1, nsig

      prespa = pres(k) * 100.0_r8

      if ( tair(k) /= tair_miss .and. qc_tair(k) == 0 ) then
  
        oerr = rawin_temp_error(pres(k))
        if ( tair(k) >= 180.0_r8 .and. &
             tair(k) <= 330.0_r8 .and. oerr /= missing_r8 ) then

          call create_3d_obs(lat, lon, prespa, VERTISPRESSURE, tair(k), &
                             RADIOSONDE_TEMPERATURE, oerr, oday, osec, qc, obs)
          call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
        endif
    
      endif
  
      ! if the air and dewpoint obs are both ok, then see which of the possible
      ! three types of moisture obs to generate.
      if ( tair(k) /= tair_miss .and. qc_tair(k) == 0 .and. &
           tdew(k) /= tdew_miss .and. qc_tdew(k) == 0  ) then

        ! tdew is the dewpoint depression
        dptk = tair(k) - tdew(k)

        if ( include_specific_humidity ) then

          qobs = specific_humidity(sat_vapor_pressure(dptk),    prespa)
          qsat = specific_humidity(sat_vapor_pressure(tair(k)), prespa)
          if ( LH_err ) then
            qerr = rh_error_from_dewpt_and_temp(tair(k), dptk)
          else
            qerr = rawin_rel_hum_error(pres(k), tair(k), qobs / qsat)
          endif
          oerr = max(qerr * qsat, 0.0001_r8)
          if ( qobs >  0.0_r8  .and. &
               qobs <= 0.07_r8 .and. qerr /= missing_r8 ) then
  
            call create_3d_obs(lat, lon, prespa, VERTISPRESSURE, qobs, &
                               RADIOSONDE_SPECIFIC_HUMIDITY, oerr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
          endif
  
        endif
  
        if ( include_relative_humidity ) then

          rh = temp_and_dewpoint_to_rh(tair(k), dptk)
          if ( LH_err ) then
            oerr = rh_error_from_dewpt_and_temp(tair(k), dptk)
          else
            oerr = rawin_rel_hum_error(pres(k), tair(k), rh)
          endif
          
          if ( rh >  0.0_r8 .and. &
               rh <= 1.5_r8 .and. oerr /= missing_r8 ) then

            call create_3d_obs(lat, lon, prespa, VERTISPRESSURE, rh, &
                               RADIOSONDE_RELATIVE_HUMIDITY, oerr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
          endif

        endif
  
        if ( include_dewpoint ) then
  
          rh = temp_and_dewpoint_to_rh(tair(k), dptk)
          oerr = dewpt_error_from_rh_and_temp(tair(k), rh)
  
          if ( rh >  0.0_r8 .and. &
               rh <= 1.5_r8 .and. oerr /= missing_r8 ) then

            call create_3d_obs(lat, lon, prespa, VERTISPRESSURE, dptk, &
                               RADIOSONDE_DEWPOINT, oerr, oday, osec, qc, obs)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
  
          endif

        endif

      endif  ! quality control/missing check on tair and tdew

    end do
    deallocate(pres, tair, tdew, qc_pres, qc_tair, qc_tdew)

  endif

  !  If desired, read the significant-level wind data, write to obs_seq
  call getvar_int_1d_1val(ncid, "numSigW", n, nsig )
  
  if ( sigwnd .and. nsig <= nmaxsw ) then

    allocate(pres(nsig))     ;  allocate(wdir(nsig))     ;  allocate(wspd(nsig))
    allocate(qc_pres(nsig))  ;  allocate(qc_wdir(nsig))  ;  allocate(qc_wspd(nsig))

    !  read significant level data
    call getvar_real_2d(ncid, "htSigW", n, nsig, pres, pres_miss)
    call getvar_real_2d(ncid, "wdSigW", n, nsig, wdir, wdir_miss)
    call getvar_real_2d(ncid, "wsSigW", n, nsig, wspd, wspd_miss)

    if (use_input_qc) then
       call get_or_fill_QC_2d(ncid, "wdSigTQCR", n, nsig, qc_wdir)
       call get_or_fill_QC_2d(ncid, "wsSigTQCR", n, nsig, qc_wspd)
    else
       qc_wdir = 0
       qc_wspd = 0
    endif

    do k = 1, nsig

      !  add data to the observation sequence here.
      if ( wdir(k) /= wdir_miss .and. qc_wdir(k) == 0 .and. &
           wspd(k) /= wspd_miss .and. qc_wspd(k) == 0  ) then

        call wind_dirspd_to_uv(wdir(k), wspd(k), uwnd, vwnd)
        oerr = rawin_wind_error(pres_alt_to_pres(pres(k)) * 0.01_r8)
        if ( abs(uwnd) <= 150.0_r8 .and. & 
             abs(vwnd) <= 150.0_r8 .and. oerr /= missing_r8 ) then

          call create_3d_obs(lat, lon, pres(k), VERTISHEIGHT, uwnd, &
                             RADIOSONDE_U_WIND_COMPONENT, oerr, oday, osec, qc, obs)
          call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

          call create_3d_obs(lat, lon, pres(k), VERTISHEIGHT, vwnd, &
                             RADIOSONDE_V_WIND_COMPONENT, oerr, oday, osec, qc, obs)
          call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)

        endif

      endif

    end do
    deallocate(pres, wdir, wspd, qc_pres, qc_wdir, qc_wspd)

  endif

  nused = nused + 1
  latu(nused) = lat
  lonu(nused) = lon

enddo sondeloop

! have to close at end of loop, unlike other versions of the madis converters
call nc_check( nf90_close(ncid), &
               'convert_madis_rawin', 'closing file '//trim(rawin_in_file))

! if we added any obs to the sequence, write it now.
if ( get_num_obs(obs_seq) > 0 )  call write_obs_seq(obs_seq, rawin_out_file)

!end of main program
call finalize_utilities()

contains

! specialized versions of the netcdf get routines that seem to be
! pretty specific to this version of the code, so i didn't put them
! in the general observations utilities file.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   getvar_real_1d_1val - subroutine that inquires, gets the variable, and fills 
!            in the missing value attribute if that arg is present.
!            takes a single start, uses count=1, returns a scalar
!
!      ncid - open netcdf file handle
!      varname - string name of netcdf variable
!      start - starting index in the 1d array
!      dout - output value.  real(r8)
!      dmiss - value that signals a missing value   real(r8), optional
!
!     created 11 Mar 2010,  nancy collins,  ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getvar_real_1d_1val(ncid, varname, start, dout, dmiss)
 integer,            intent(in)   :: ncid
 character(len = *), intent(in)   :: varname
 integer,            intent(in)   :: start
 real(r8),           intent(out)  :: dout
 real(r8), optional, intent(out)  :: dmiss

integer :: varid

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'getvar_real', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, dout, start = (/ start /) ), &
               'getvar_real', 'getting var '// trim(varname))

if (present(dmiss)) &
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', dmiss), &
               'getvar_real', 'getting attr "_FillValue" for '//trim(varname))

end subroutine getvar_real_1d_1val

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   getvar_int_1d_1val - subroutine that inquires, gets the variable, and fills 
!            in the missing value attribute if that arg is present.
!            takes a single start, uses count=1, returns a scalar
!
!      ncid - open netcdf file handle
!      varname - string name of netcdf variable
!      start - starting index in the 1d array
!      dout - output value.  int
!      dmiss - value that signals a missing value   int, optional
!
!     created 11 Mar 2010,  nancy collins,  ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getvar_int_1d_1val(ncid, varname, start, dout, dmiss)
 integer,            intent(in)   :: ncid
 character(len = *), intent(in)   :: varname
 integer,            intent(in)   :: start
 integer,            intent(out)  :: dout
 integer,  optional, intent(out)  :: dmiss

integer :: varid

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'getvar_int_1d_1val', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, dout, start = (/ start /) ), &
               'getvar_int_1d_1val', 'getting var '// trim(varname))

if (present(dmiss)) &
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', dmiss), &
               'getvar_int_1d_1val', 'getting attr "_FillValue" for '//trim(varname))

end subroutine getvar_int_1d_1val

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   getvar_real_2d - subroutine that inquires, gets the variable, and fills 
!            in the missing value attribute if that arg is present.
!     SPECIALIZED for this use - assumes start = (/ 1, n /) and count = (/ m, 1 /)
!           so takes a scalar start, count, returns a 1d_array
!
!      ncid - open netcdf file handle
!      varname - string name of netcdf variable
!      start - starting index in the 2d array.  integer
!      count - nitems to get. integer
!      darray - output array.  real(r8)
!      dmiss - value that signals a missing value   real(r8), optional
!
!     created 11 Mar 2010,  nancy collins,  ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine getvar_real_2d(ncid, varname, start, count, darray, dmiss)
 integer,            intent(in)   :: ncid
 character(len = *), intent(in)   :: varname
 integer,            intent(in)   :: start
 integer,            intent(in)   :: count
 real(r8),           intent(out)  :: darray(:)
 real(r8), optional, intent(out)  :: dmiss

integer :: varid

! read the data for the requested array, and get the fill value
call nc_check( nf90_inq_varid(ncid, varname, varid), &
               'getvar_real_2d', 'inquire var '// trim(varname))
call nc_check( nf90_get_var(ncid, varid, darray, &
                start=(/ 1, start /), count=(/ count, 1 /) ), &
               'getvar_real_2d', 'getting var '// trim(varname))

if (present(dmiss)) &
   call nc_check( nf90_get_att(ncid, varid, '_FillValue', dmiss), &
               'getvar_real_2d', 'getting attr "_FillValue" for '//trim(varname))

end subroutine getvar_real_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   get_or_fill_QC_2d - subroutine which gets the requested netcdf variable
!           but if it isn't there, it fills the array with 0s.  not an
!           error if it's not present.  assumes integer data array
!     SPECIALIZED for this use - assumes start = (/ 1, n /) and count = (/ m, 1 /)
!           so takes a scalar start, count, returns a 1d_array
!           also prints out a message if fill used.
!
!      ncid - open netcdf file handle
!      varname - string name of netcdf variable
!      start - starting index in the 2d array.  integer
!      count - nitems to get. integer
!      darray - output array.  integer
!
!     created Mar 8, 2010    nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_or_fill_QC_2d(ncid, varname, start, count, darray)
 integer,            intent(in)    :: ncid
 character(len = *), intent(in)    :: varname
 integer,            intent(in)    :: start
 integer,            intent(in)    :: count
 integer,            intent(inout) :: darray(:)

integer :: varid, nfrc

! test to see if variable is present.  if yes, read it in.
! otherwise, set to fill value, or 0 if none given.

nfrc = nf90_inq_varid(ncid, varname, varid)
if (nfrc == NF90_NOERR) then
   call nc_check( nf90_get_var(ncid, varid, darray, &
                  start=(/ 1, start /), count=(/ count, 1 /) ), &
                  'get_or_fill_int_2d', 'reading '//trim(varname) )
else
   darray = 0
   if (start == 1) & 
     print *, 'QC field named ' // trim(varname) // ' was not found in input, 0 used instead'
endif

end subroutine get_or_fill_QC_2d

end program convert_madis_rawin
