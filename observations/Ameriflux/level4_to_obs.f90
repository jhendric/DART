! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program level4_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   level4_to_obs - a program that only needs minor customization to read
!      in a text-based dataset - either white-space separated values or
!      fixed-width column data.
!
!     created 3 May 2012   Tim Hoar NCAR/IMAGe
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use         types_mod, only : r8, MISSING_R8

use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              register_module, error_handler, E_MSG, E_ERR, &
                              open_file, close_file, do_nml_file, do_nml_term, &
                              check_namelist_read, find_namelist_in_file, &
                              nmlfileunit

use  time_manager_mod, only : time_type, set_calendar_type, GREGORIAN, &
                              set_date, set_time, get_time, print_time, &
                              print_date, operator(-), operator(+), operator(>), &
                              operator(<), operator(==), operator(<=), operator(>=)

use      location_mod, only : VERTISHEIGHT

use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, &
                              init_obs_sequence, get_num_obs, &
                              set_copy_meta_data, set_qc_meta_data

use      obs_kind_mod, only : TOWER_SENSIBLE_HEAT_FLUX, &
                              TOWER_NETC_ECO_EXCHANGE,  &
                              TOWER_LATENT_HEAT_FLUX

implicit none

!-----------------------------------------------------------------------
! version controlled file description for error handling, do not edit
!-----------------------------------------------------------------------

character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!-----------------------------------------------------------------------
! Namelist with default values
!-----------------------------------------------------------------------

character(len=128) :: text_input_file = 'textdata.input'
character(len=128) :: obs_out_file    = 'obs_seq.out'
integer            :: year
real(r8)           :: timezoneoffset
real(r8)           :: latitude
real(r8)           :: longitude
real(r8)           :: elevation
real(r8)           :: flux_height
real(r8)           :: maxgoodqc       = 3
logical            :: verbose         = .false.

namelist /level4_to_obs_nml/ text_input_file, obs_out_file, year, &
             timezoneoffset, latitude, longitude, elevation, &
             flux_height, maxgoodqc, verbose

!-----------------------------------------------------------------------
! globally-scoped variables
!-----------------------------------------------------------------------

character(len=256)      :: input_line, string1, string2, string3
integer                 :: iline, nlines
logical                 :: file_exist, first_obs
integer                 :: n, i, oday, osec, rcio, iunit
integer                 :: num_copies, num_qc, max_obs
real(r8)                :: oerr, qc
type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: time_obs, prev_time, offset
real(r8), parameter     :: umol_to_gC = (1.0_r8/1000000.0_r8) * 12.0_r8

type towerdata
  type(time_type)   :: time_obs
  character(len=20) :: monthstring = 'month'
  character(len=20) :: daystring   = 'day'
  character(len=20) :: hourstring  = 'hour'
  character(len=20) :: doystring   = 'doy'
  character(len=20) :: neestring   = 'nee_or_fmds'
  character(len=20) :: neeQCstring = 'nee_or_fmdsqc'
  character(len=20) :: lestring    = 'le_f'
  character(len=20) :: leQCstring  = 'le_fqc'
  character(len=20) :: hstring     = 'h_f'
  character(len=20) :: hQCstring   = 'h_fqc'
  integer  :: monthindex
  integer  :: dayindex
  integer  :: hourindex
  integer  :: doyindex
  integer  :: neeindex
  integer  :: neeQCindex
  integer  :: leindex
  integer  :: leQCindex
  integer  :: hindex
  integer  :: hQCindex
  integer  :: month
  integer  :: day
  real(r8) :: hour
  real(r8) :: doy
  real(r8) :: nee
  integer  :: neeQC
  real(r8) :: le
  integer  :: leQC
  real(r8) :: h
  integer  :: hQC
end type towerdata

type(towerdata) :: tower

!-----------------------------------------------------------------------
! start of executable code
!-----------------------------------------------------------------------

call initialize_utilities('level4_to_obs')

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "level4_to_obs_nml", iunit)
read(iunit, nml = level4_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, "level4_to_obs_nml")

! Record the namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=level4_to_obs_nml)
if (do_nml_term()) write(     *     , nml=level4_to_obs_nml)

! time setup
call set_calendar_type(GREGORIAN)
offset    = set_time(nint(abs(timezoneoffset)*3600.0_r8),0)
prev_time = set_time(0, 0)

if (verbose) print *, 'tower located at lat, lon, elev  =', latitude, longitude, elevation
if (verbose) print *, 'flux observations taken at       =', flux_height,'m'

! check the lat/lon values to see if they are ok
if (longitude < 0.0_r8) longitude = longitude + 360.0_r8

if (( latitude > 90.0_r8 .or. latitude  <  -90.0_r8 ) .or. &
    (longitude <  0.0_r8 .or. longitude >  360.0_r8 )) then

   write (string2,*)'latitude  should be [-90, 90] but is ',latitude
   write (string3,*)'longitude should be [  0,360] but is ',longitude

   string1 ='tower location error in input.nml&level4_to_obs_nml'
   call error_handler(E_ERR,'level4_to_obs', string1, source, revision, &
                      revdate, text2=string2,text3=string3)

endif

! We need to know the maximum number of observations in the input file.
! Each line has info for the 3 observations we want.
! The max possible number of obs needs to be specified but it
! will only write out the actual number created.
! Each observation in this series will have a single
! observation value and a quality control flag.  
! Initialize two empty observations - one to track location
! in observation sequence - the other is for the new observation.

iunit = open_file(text_input_file, 'formatted', 'read')
if (verbose) print *, 'opened input file ' // trim(text_input_file)

nlines     = count_file_lines(iunit)
max_obs    = 3*nlines
num_copies = 1
num_qc     = 1
first_obs  = .true.

call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

! the first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call set_qc_meta_data(  obs_seq, 1, 'Ameriflux QC')

! The first line describes all the fields ... column headers, if you will

rewind(iunit)
call decode_header(iunit)

obsloop: do iline = 2,nlines

   ! read in entire text line into a buffer
   read(iunit,'(A)',iostat=rcio) input_line
   if (rcio < 0) exit obsloop
   if (rcio > 0) then
      write (string1,'(''Cannot read (error '',i3,'') line '',i8,'' in '',A)') &
                    rcio, iline, trim(text_input_file)
      call error_handler(E_ERR,'main', string1, source, revision, revdate)
   endif

   ! parse the line into the tower structure (including the observation time)
   call stringparse(input_line, iline)

   if (iline <= 2) then
      write(*,*)''
      write(*,*)'Check of the first observation: (column,string,value)'
      write(*,*)tower%monthindex, tower%monthstring , tower%month
      write(*,*)tower%dayindex  , tower%daystring   , tower%day
      write(*,*)tower%hourindex , tower%hourstring  , tower%hour
      write(*,*)tower%doyindex  , tower%doystring   , tower%doy
      write(*,*)tower%hindex    , tower%hstring     , tower%h
      write(*,*)tower%hQCindex  , tower%hQCstring   , tower%hQC
      write(*,*)tower%leindex   , tower%lestring    , tower%le
      write(*,*)tower%leQCindex , tower%leQCstring  , tower%leQC
      write(*,*)tower%neeindex  , tower%neestring   , tower%nee
      write(*,*)tower%neeQCindex, tower%neeQCstring , tower%neeQC
      call print_date(tower%time_obs, 'observation date is')
      call print_time(tower%time_obs, 'observation time is')
   end if

   if (verbose) call print_date(tower%time_obs, 'obs time is')

   call get_time(tower%time_obs, osec, oday)

   ! make an obs derived type, and then add it to the sequence
   ! If the QC value is good, use the observation.
   ! Increasingly larger QC values are more questionable quality data.
   ! The observation errors are from page 183, Table 7.1(A) in 
   ! Chapter 7 of a book by A.D. Richardson et al. via Andy Fox.  

   if (tower%hQC <= maxgoodqc) then   ! Sensible Heat Flux [W m-2]
      oerr = 10.0_r8 + abs(tower%h)*0.22_r8
      qc   = real(tower%hQC,r8)
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%h, &
                         TOWER_LATENT_HEAT_FLUX, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)
   endif

   if (tower%leQC <= maxgoodqc) then   ! Latent Heat Flux [W m-2]
      oerr = 10.0_r8 + abs(tower%le)*0.32_r8
      qc   = real(tower%leQC,r8)
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%le, &
                         TOWER_SENSIBLE_HEAT_FLUX, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)
   endif

   if (tower%neeQC <= maxgoodqc) then   ! Net Ecosystem Exchange  [umol m-2 s-1]
      if (tower%nee <= 0) then
         oerr = (2.0_r8 + abs(tower%nee)*0.1_r8) * umol_to_gC
      else
         oerr = (2.0_r8 + abs(tower%nee)*0.4_r8) * umol_to_gC
      endif
      tower%nee = -tower%nee * umol_to_gC ! to match convention in CLM [gC m-2 s-1]
      qc        = real(tower%neeQC,r8)
      call create_3d_obs(latitude, longitude, flux_height, VERTISHEIGHT, tower%nee, &
                         TOWER_NETC_ECO_EXCHANGE, oerr, oday, osec, qc, obs)
      call add_obs_to_seq(obs_seq, obs, tower%time_obs, prev_obs, prev_time, first_obs)
   endif

end do obsloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (verbose) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   create_3d_obs - subroutine that is used to create an observation
!                   type from observation data.
!
!       NOTE: assumes the code is using the threed_sphere locations module,
!             that the observation has a single data value and a single
!             qc value, and that this obs type has no additional required
!             data (e.g. gps and radar obs need additional data per obs)
!
!    lat   - latitude of observation
!    lon   - longitude of observation
!    vval  - vertical coordinate
!    vkind - kind of vertical coordinate (pressure, level, etc)
!    obsv  - observation value
!    okind - observation kind
!    oerr  - observation error (in units of standard deviation)
!    day   - gregorian day
!    sec   - gregorian second
!    qc    - quality control value
!    obs   - observation type
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!     adapted for more generic use 11 Mar 2010, nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_3d_obs(lat, lon, vval, vkind, obsv, okind, oerr, day, sec, qc, obs)
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location
use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
use time_manager_mod, only : time_type, set_time
use     location_mod, only : set_location

 integer,        intent(in)    :: okind, vkind, day, sec
 real(r8),       intent(in)    :: lat, lon, vval, obsv, oerr, qc
 type(obs_type), intent(inout) :: obs

real(r8)           :: obs_val(1), qc_val(1)
type(obs_def_type) :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
call set_obs_def_kind(obs_def, okind)
call set_obs_def_time(obs_def, set_time(sec, day))
call set_obs_def_error_variance(obs_def, oerr * oerr)
call set_obs_def(obs, obs_def)

obs_val(1) = obsv
call set_obs_values(obs, obs_val)
qc_val(1)  = qc
call set_qc(obs, qc_val)

end subroutine create_3d_obs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   add_obs_to_seq -- adds an observation to a sequence.  inserts if first
!           obs, inserts with a prev obs to save searching if that's possible.
!
!     seq - observation sequence to add obs to
!     obs - observation, already filled in, ready to add
!     obs_time - time of this observation, in dart time_type format
!     prev_obs - the previous observation that was added to this sequence
!                (will be updated by this routine)
!     prev_time - the time of the previously added observation (will also
!                be updated by this routine)
!     first_obs - should be initialized to be .true., and then will be
!                updated by this routine to be .false. after the first obs
!                has been added to this sequence.
!
!     created Mar 8, 2010   nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine add_obs_to_seq(seq, obs, obs_time, prev_obs, prev_time, first_obs)

use obs_sequence_mod, only : obs_sequence_type, obs_type, insert_obs_in_seq
use time_manager_mod, only : time_type, operator(>=)

type(obs_sequence_type), intent(inout) :: seq
type(obs_type),          intent(inout) :: obs, prev_obs
type(time_type),         intent(in)    :: obs_time
type(time_type),         intent(inout) :: prev_time
logical,                 intent(inout) :: first_obs

! insert(seq,obs) always works (i.e. it inserts the obs in
! proper time format) but it can be slow with a long file.
! supplying a previous observation that is older (or the same
! time) as the new one speeds up the searching a lot.

if(first_obs) then    ! for the first observation, no prev_obs
   call insert_obs_in_seq(seq, obs)
   first_obs = .false.
else
   if(obs_time >= prev_time) then  ! same time or later than previous obs
      call insert_obs_in_seq(seq, obs, prev_obs)
   else                            ! earlier, search from start of seq
      call insert_obs_in_seq(seq, obs)
   endif
endif

! update for next time
prev_obs  = obs
prev_time = obs_time

end subroutine add_obs_to_seq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   count_file_lines --
!           count the lines in a text file.
!           rewinds the unit after counting.
!
!     iunit - handle to the already-open text file
!
!     created May 2, 2012   Tim Hoar, NCAR/IMAGe
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function count_file_lines(iunit)

integer, intent(in) :: iunit
integer :: count_file_lines

integer :: i
character(len=128) :: oneline

integer, parameter :: tenmillion = 10000000
rewind(iunit)

count_file_lines = 0
countloop : do i = 1,tenmillion

   read(iunit,'(A)',iostat=rcio) oneline

   if (rcio < 0) exit countloop ! end of file
   if (rcio > 0) then
      write (string1,'('' read around line '',i8)')i
      call error_handler(E_ERR,'count_file_lines', string1, &
                         source, revision, revdate)
   endif
   count_file_lines = count_file_lines + 1

enddo countloop
rewind(iunit)

if (count_file_lines >= tenmillion) then
   write (string1,'('' suspiciously large number of lines '',i8)')count_file_lines
   call error_handler(E_MSG,'count_file_lines', string1, &
                         source, revision, revdate)
endif

end function count_file_lines




subroutine decode_header(iunit)
! Reads the first line of the header and parses the information.
! FIXME ... decode the header ... do not assume ...
integer, intent(in) :: iunit

read(iunit,'(A)',iostat=rcio) input_line
if (rcio /= 0) then
  write(string1,*)'Cannot parse header. Begins <',trim(input_line(1:40)),'>'
  call error_handler(E_ERR,'decode_header',string1, source, revision, revdate)
endif

call error_handler(E_MSG,'decode_header','hardcoding values for now ... dangerous', &
                     source, revision, revdate)

tower%monthindex = 1
tower%dayindex   = 2
tower%hourindex  = 3
tower%doyindex   = 4
tower%hindex     = 15
tower%hQCindex   = 16
tower%leindex    = 17
tower%leQCindex  = 18
tower%neeindex   = 26
tower%neeQCindex = 27

end subroutine decode_header



subroutine stringparse(str1,linenum)
! just declare everything as reals and chunk it

character(len=*), intent(in) :: str1
integer         , intent(in) :: linenum

real(r8), dimension(34) :: values
integer :: iday, ihour, imin, isec, seconds
type(time_type) :: time0, time1, time2

values = MISSING_R8

read(str1,*,iostat=rcio) values
if (rcio /= 0) then
  write(string1,*)'Cannot parse line',linenum,'. Begins <',trim(str1(1:40)),'>'
  call error_handler(E_ERR,'stringparse',string1, source, revision, revdate)
endif

! Stuff what we want into the tower structure
!
! Convert to 'CLM-friendly' units AFTER we determine observation error variance.
! That happens in the main routine.
!
! NEE_or_fMDS    has units     [umolCO2 m-2 s-1] 
! H_f            has units     [W m-2]
! LE_f           has units     [W m-2]
!
! (CLM) NEE      has units     [gC m-2 s-1]

tower%month = nint(values(tower%monthindex))
tower%day   = nint(values(tower%dayindex  ))
tower%hour  =      values(tower%hourindex )
tower%doy   =      values(tower%doyindex  )
tower%nee   =      values(tower%neeindex  )
tower%neeQC = nint(values(tower%neeQCindex))
tower%le    =      values(tower%leindex   )
tower%leQC  = nint(values(tower%leQCindex ))
tower%h     =      values(tower%hindex    )
tower%hQC   = nint(values(tower%hQCindex  ))

! decode the time pieces ... two times ...
! The LAST line of these files is knackered ... and we have to check that
! if the doy is greater than the ymd ...
! 12,31,23.5,366.979    N-1
!  1, 1, 0.0,367.000    N

iday    = int(tower%doy)
ihour   = int(tower%hour)
seconds = nint((tower%hour - real(ihour,r8))*3600)
imin    = seconds / 60
isec    = seconds - imin * 60
time0   = set_date(year, tower%month, tower%day, ihour, imin, isec)

isec    = ihour*3600 + imin*60 + isec
time1   = set_date(year,1,1,0,0,0) + set_time(isec, (iday-1))
time2   = time0 - time1

call get_time(time2, isec, iday)

if ( iday > 0 ) then
   ! we need to change the day ...

   tower%time_obs = time1

   if (verbose) then
      write(string1,*)'converting ',tower%month,tower%day,tower%hour,tower%doy
      write(string2,*)'the day-of-year indicates we should amend the month/day values.' 
      call error_handler(E_MSG,'stringparse', string1, source, revision, &
                      revdate, text2=string2)

      call print_date(time0, 'stringparse: using ymd date is')
      call print_date(time1, 'stringparse: using doy date is')
      call print_time(time0, 'stringparse: using ymd time is')
      call print_time(time1, 'stringparse: using doy time is')
      call print_time(time2, 'stringparse: difference     is')
   endif
else

   tower%time_obs = time0

endif

! 8AM East Coast is 1PM Greenwich 
if (timezoneoffset < 0.0_r8) then
   tower%time_obs = tower%time_obs + offset
else
   tower%time_obs = tower%time_obs - offset
endif

! The QC values can be 'missing' ... in which case the values are too

if (tower%neeQC < 0) tower%neeQC = maxgoodqc + 1000 
if (tower%leQC  < 0) tower%leQC  = maxgoodqc + 1000
if (tower%hQC   < 0) tower%hQC   = maxgoodqc + 1000

! if (tower%neeQC < maxgoodqc) then
!    write(*,*)'nee umol m-2 s-1 ',tower%nee
!    write(*,*)'nee   gC m-2 s-1 ',tower%nee*umol_to_gC
! endif
 
end subroutine stringparse



end program level4_to_obs


! LEVEL 4 VARIABLE DESCRIPTION
!
! Variables description:
! Level 4 data are obtained from the level 3 products, data are ustar filtered,
! gap-filled using different methods and partitioned.
! Datasets are also aggregated from daily to monthly.
! Flags with information regarding quality of the original and gapfilled data are added.
!
! Half hourly dataset variables description:
!
! - Month          : from 1 to 12
! - Day            : day of the month
! - Hour           : from 0 to 23.5, indicates the end of the half hour of measurement
! - DoY            : decimal day of the year
! - Rg_f           : global radiation filled [W m-2]
! - Rg_fqc         : global radiation quality flags:
!                    0 = original, 1 = A (most reliable), 2 = B (medium), 3 = C (least reliable).
!                    (Refer to Reichstein et al. 2005 Global Change Biology )
! - Ta_f           : air temperature filled [�C]
! - Ta_fqc         : air temperature quality flags:
!                    0 = original, 1 = A (most reliable), 2 = B (medium), 3 = C (least reliable).
!                    (Refer to Reichstein et al. 2005 Global Change Biology )
! - VPD_f          : vapour pressure deficit [hPa]
! - VPD_fqc        : vapour pressure deficit quality flags:
!                    0 = original, 1 = A (most reliable), 2 = B (medium), 3 = C (least reliable).
!                    (Refer to Reichstein et al. 2005 Global Change Biology )
! - Ts_f           : soil temperature filled [�C]
! - Ts_fqc         : soil temperature quality flags:
!                    0 = original, 1 = A (most reliable), 2 = B (medium), 3 = C (least reliable).
!                    (Refer to Reichstein et al. 2005 Global Change Biology )
! - Precip         : precipitation [mm]
! - SWC            : soil water content [%vol]
! - H_f            : sensible heat flux filled [W m-2]
! - H_fqc          : sensible heat flux quality flags:
!                    0 = original, 1 = A (most reliable), 2 = B (medium), 3 = C (least reliable).
!                    (Refer to Reichstein et al. 2005 Global Change Biology )
! - LE_f           : latent heat flux filled [W m-2]
! - LE_fqc         : latent heat flux quality flags:
!                    0 = original, 1 = A (most reliable), 2 = B (medium), 3 = C (least reliable).
!                    (Refer to Reichstein et al. 2005 Global Change Biology )
! - qf_NEE_st      : fluxes quality flags as defined in the Level3 product
! - qf_NEE_or      : fluxes quality flags as defined in the Level3 product
! - Reco_st        : Estimated ecosystem respiration according to the short-term temperature
!                    response of night-time fluxes based on NEE_st
!                    (Refer to Reichstein et al. 2005 Global Change Biology )
!                    [umolCO2 m-2 s-1]
! - Reco_or        : Estimated ecosystem respiration according to the short-term temperature
!                    response of night-time fluxes based on NEE_or
!                    (Refer to Reichstein et al. 2005 Global Change Biology )
!                    [umolCO2 m-2 s-1]
! - NEE_st_fMDS    : NEE_st filled using the Marginal Distribution Sampling method
!                    (Refer to Reichstein et al. 2005 Global Change Biology )
!                    [umolCO2 m-2 s-1]
! - NEE_st_fMDSqc  : NEE_st_fMDS quality flags:
!                    0 = original, 1 = A (most reliable), 2 = B (medium), 3 = C (least reliable).
!                    (Refer to Reichstein et al. 2005 Global Change Biology )
! - GPP_st_MDS     : Gross Primary Production calculated as GPP_st_MDS = Reco_st - NEE_st_MDS
!                    [umolCO2 m-2 s-1]
! - NEE_or_fMDS    : NEE_or filled using the Marginal Distribution Sampling method
!                    (Refer to Reichstein et al. 2005 Global Change Biology )
!                    [umolCO2 m-2 s-1]
! - NEE_or_fMDSqc  : NEE_or_fMDS quality flags:
!                    0 = original, 1 = A (most reliable), 2 = B (medium), 3 = C (least reliable).
!                    (Refer to Reichstein et al. 2005 Global Change Biology )
! - GPP_or_MDS     : Gross Primary Production calculated as GPP_or_MDS = Reco_or - NEE_or_MDS
!                    [umolCO2 m-2 s-1]
! - NEE_st_fANN    : NEE_st filled using the Artificial Neural Network method
!                    (Refer to Papale et al. 2003 Global Change Biology and to the Other Information section in this document)
!                    [umolCO2 m-2 s-1]
! - NEE_st_fANNqc  : NEE_st_fANN quality flags:
!                    0 = original, 1 = filled using original meteorological inputs or filled with qc=1,
!                    2 = filled using filled meteorological inputs with qc=2 or 3,
!                    3 = not filled using ANN due to one or more input missed but filled with the MDS method
! - GPP_st_ANN     : Gross Primary Production calculated as GPP_st_ ANN = Reco_st - NEE_st_ ANN
!                    [umolCO2 m-2 s-1]
! - NEE_or_f ANN   : NEE_or filled using the Artificial Neural Network method
!                    (Refer to Papale et al. 2003 Global Change Biology and to the Other Information section in this document)
!                    [umolCO2 m-2 s-1]
! - NEE_or_f ANNqc : NEE_or_fANN quality flags:
!                    0 = original, 1 = filled using original meteorological inputs or filled with qc=1,
!                    2 = filled using filled meteorological inputs with qc=2 or 3,
!                    3 = not filled using ANN due to one or more input missed but filled with the MDS method
! - GPP_or_ ANN    : Gross Primary Production calculated as GPP_or_ ANN = Reco_or - NEE_or_ ANN
!                    [umolCO2 m-2 s-1]
