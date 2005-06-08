! Data Assimilation Research Testbed -- DART
! Copyright 2004, 2005, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 
module model_mod

! <next five lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
! $Name$ 
!
! This is a non-divergent barotropic model on the sphere. Currently makes
! use of NAG based transforms which are not available on NCAR systems.

use types_mod, only : r8
use transforms_mod
use ncd_file_mod
use nag_wrap_mod,     only : g05ddf_wrap
use loc_and_dist_mod, only : loc_type, get_dist, set_loc

implicit none
private

public :: init_model, get_model_size, lat_max, num_lon, init_conditions, & 
   adv_1step, advance, &
   output, barot_to_dp, dp_to_barot, delta_t, adv_true_state, &
   dp_to_grid, lon, lat, model_state_location, diag_output_index, &
   num_fourier, num_spherical, model_output, trans_spherical_to_grid, &
   get_close_pts, grid_to_dp, state_loc, model_get_close_states

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
   source   = "$Source$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

! Flag for using real data or perfect model
logical, parameter :: use_real_data = .false.

! Truncation for T42 follows:
! Reduced grid
!integer, parameter :: num_windows = 1, lat_max = 44, num_lon = 84
 integer, parameter :: num_windows = 1, lat_max = 64, num_lon = 128
 integer, parameter :: num_fourier = 42, fourier_inc = 1
 integer, parameter :: num_spherical = 43

! Truncation for T21 follows:
!integer, parameter :: num_windows = 1, lat_max = 32, num_lon = 64
!integer, parameter :: num_windows = 1, lat_max = 22, num_lon = 42
!integer, parameter :: num_fourier = 21, fourier_inc = 1
!integer, parameter :: num_spherical = 22 

! Truncation for T10 follows
!integer, parameter :: num_windows = 1, lat_max = 16, num_lon = 32
!integer, parameter :: num_fourier = 10, fourier_inc = 1
!integer, parameter :: num_spherical = 11 

real, parameter :: radius = 6.4e6, omega = 7.292e-5

real :: dif_days, delta_t, real_time
complex :: force(0:num_fourier, 0:num_spherical)
integer :: fourier_lim, spherical_lim


! Following is for standard dynamical systems interface; physical space
integer, parameter :: model_size = lat_max * num_lon

! Definitions of lats and lons
real(r8) :: lat(lat_max), lon(num_lon)

! Define the location of the state variables in module storage
type(loc_type) :: state_loc(model_size)

! Define output indices for diagnostics
integer :: diag_output_index(9) 

contains



  subroutine output(x, time)
!---------------------------------------------------------------------
! subroutine output(x, time)

implicit none

real, intent(in) :: x(model_size)
real, intent(in) :: time

end subroutine output



  function model_state_location()
!---------------------------------------------------------------------
! function model_state_location()

implicit none

type (loc_type) :: model_state_location(model_size)

integer :: i, j, index
real    :: rlat(lat_max), rlon(num_lon)

! Compute the lat and lons for the Gaussian grid and put them in storage

call get_deg_lat(rlat)
call get_deg_lon(rlon)

lat = rlat
lon = rlon

! Load these into structure

do i = 1, num_lon
   do j = 1, lat_max
      index = j + (i - 1)*lat_max
      call set_loc(model_state_location(index), lon(i), lat(j))
   end do
end do

end function model_state_location



  subroutine barot_init(dif_days_in, delt, force_in, fourier_lim_in, spherical_lim_in)
!---------------------------------------------------------------------
! subroutine barot_init(dif_days_in, delt, force_in, fourier_lim_in, spherical_lim_in)
!
! Calls the initialization routines for the spherical harmonic transforms.
! Sets del8 diffusion time on smallest wave, forcing coefficient and
! forcing patter.

implicit none

real,    intent(in) :: dif_days_in, delt
complex, intent(in) :: force_in(0:num_fourier, 0:num_spherical)
integer, intent(in) :: fourier_lim_in, spherical_lim_in

real    :: rlat(lat_max), rlon(num_lon)
integer :: i, j

call initialize_transforms(radius, num_windows, lat_max, num_lon, &
   num_fourier, fourier_inc, num_spherical, .false., .true., .true., 0.0)

dif_days      = dif_days_in
delta_t       = delt
real_time     = 0.0
force         = force_in
fourier_lim   = fourier_lim_in
spherical_lim = spherical_lim_in

! Compute the lat and lons for the Gaussian grid and put them in storage

call get_deg_lat(rlat)
call get_deg_lon(rlon)

lat = rlat
lon = rlon

! Quick temporary output of model grid
!do i = 1, num_lon
!   do j = 1, lat_max
!      if(j /= 1 .and. j /= lat_max) then
!         write(*, *) rlon(i), rlat(j), 1e11
!      else
!         write(*, *) rlon(i), 0.99 * rlat(j), 1e11
!      endif
!   end do
!end do
!if(1 == 1) stop

end subroutine barot_init



  function forwrd(psisp, nsteps)
!---------------------------------------------------------------------
! function forwrd(psisp, nsteps)
!
! Integrates barotropic model forward for nsteps from initial spectral
! streamfunction field psisp.  Timestep is delta_t, force is an
! additional forcing field, and dif_days is a timescal for del8
! diffusion on smallest retained wave vorticity. The first timestep is
!  a standard forward, the second a leapfrog, and all others are third
!  order adams-bashforth steps.

complex, intent(in), dimension(0:num_fourier, 0:num_spherical) :: psisp
integer, intent(in) :: nsteps
complex             :: forwrd(0:num_fourier, 0:num_spherical)

complex, dimension(0:num_fourier, 0:num_spherical)    :: vsp, vold, psi_temp
complex, dimension(0:num_fourier, 0:num_spherical, 3) :: dt

integer :: i

! Forced components stay put

psi_temp = psisp
do n = 0, min(num_fourier, spherical_lim)
   do m = 0, min(n, fourier_lim)
      if(m <= fourier_lim .and. n <= spherical_lim) then
         psi_temp(m, n-m) = force(m, n-m)
      end if
   end do
end do

vsp = compute_laplacian(psisp, 1)                !  compute spectral vort
dt(:, :, 1) = fordt(vsp)                         !  compute dT
call ftstpwv(vsp, vold, dt(:, :, 1), delta_t)    !  do forward time step

!do n = 0, num_fourier
!   do m = 0, n
do n = 0, min(num_fourier, spherical_lim)
   do m = 0, min(n, fourier_lim)
      if(m <= fourier_lim .and. n <= spherical_lim) then
         forwrd(m, n-m) = force(m, n-m)
      end if
   end do
end do

if(nsteps == 1) return                           !  quit here if only one step

dt(:, :, 2) = fordt(vsp)                         !  do second leapfrog step; compute dt
call lfrogwv(vsp, vold, dt(:, :, 2), delta_t)
do n = 0, min(num_fourier, spherical_lim)
   do m = 0, min(n, fourier_lim)
      if(m <= fourier_lim .and. n <= spherical_lim) then
         forwrd(m, n-m) = force(m, n-m)
      end if
   end do
end do


do i = 3, nsteps                           !  do loop of adams-bashforth 3rd order steps

   index = (i - 1) - (i - 1) / 3 * 3 + 1   !  compute dt in rotating storage
   dt(:, :, index) = fordt(vsp)

   call adams3(vsp, vold, dt, index, delta_t) ! compute timestep using rotating dt storage

   do n = 0, min(num_fourier, spherical_lim)
      do m = 0, min(n, fourier_lim)
!        if(m <= 4 .and. n <= 10) then
         if(m <= 1 .and. n <= 20) then
            forwrd(m, n-m) = force(m, n-m)
         end if
      end do
   end do
end do

forwrd = compute_laplacian(vsp, -1)          ! Copy into result

do n = 0, min(num_fourier, spherical_lim)
   do m = 0, min(n, fourier_lim)
      if(m <= fourier_lim .and. n <= spherical_lim) then

!         forwrd(m, n-m) = force(m, n-m)  
! for no forcing model at any resolution
!         forwrd(m, n-m) = forwrd(m, n-m)
! Used for t21
!         forwrd(m, n-m) = forwrd(m, n-m) + 0.2*(force(m, n-m) - forwrd(m, n-m))
! Used for t42
         forwrd(m, n-m) = forwrd(m, n-m) + 0.06*(force(m, n-m) - forwrd(m, n-m))
      end if
   end do
end do

end function forwrd



  function fordt(vsp)
!---------------------------------------------------------------------
! function fordt(vsp)

implicit none

!  computes the time derivative of vorticity in spectral model

complex, dimension(0:num_fourier, 0:num_spherical), intent(in) :: vsp

complex, dimension(0:num_fourier, 0:num_spherical) :: qsp, psisp, current_state
complex, dimension(0:num_fourier, 0:num_spherical) :: fordt, diffusion

real    :: dtph(num_lon, lat_max), dif2_coef
integer :: i, n, m

psisp = compute_laplacian(vsp, -1)                 !  get spectral psi (inverse laplacian)
qsp   = vsp                                        !  get absolute vorticity
qsp(0, 1) = qsp(0, 1) + 2. * omega * sqrt(2.0/3.0)

dtph = jacob(psisp, qsp)                           !  compute physical space dt

!  convert to spectral: ??? Do trunctation at this point???

call trans_grid_to_spherical(dtph, fordt, .true.)

! Use del8 diffusion with dif_days timescale on smallest retained wave

!dif8_coef = ((num_spherical - 1) * num_spherical / radius**2)**4 * (86400 * dif_days) 

! Compute del8 diffusion by repeated calls to avoid overflow
! Use dif_days timescale diffusion on smallest retained wave

dif2_coef = ((num_spherical - 1) * num_spherical / radius**2)
diffusion = vsp
do i = 1, 4
   diffusion = compute_laplacian(diffusion, 1) / dif2_coef
end do

diffusion = diffusion / (86400.0 * dif_days)      ! Add in the timescale normalization


! Get sign, add del8 diffusion
!fordt = -1.0 * fordt - compute_laplacian(vsp, 4) / dif8_coef

fordt = -1.0 * fordt - diffusion

end function fordt



  function jacob(asp, bsp)
!---------------------------------------------------------------------
! function jacob(asp, bsp)
!
! Computes jacobian of two spectral fields and returns it in physical
! space.
!

implicit none


complex, intent(in) :: asp(:, :), bsp(:, :)
real                :: jacob(num_lon, lat_max)

real,    dimension(num_lon, lat_max) :: almph, amuph, blmph, bmuph
complex, dimension(0:num_fourier, 0:num_spherical) :: da_dlon, da_dlat, &
                                                      db_dlon, db_dlat

integer :: i

! transform lambda/mu derivatives of b

call compute_gradient_cos(asp, da_dlon, da_dlat)
call compute_gradient_cos(bsp, db_dlon, db_dlat)

call trans_spherical_to_grid(da_dlon, almph)
call trans_spherical_to_grid(da_dlat, amuph)
call trans_spherical_to_grid(db_dlon, blmph)
call trans_spherical_to_grid(db_dlat, bmuph)

jacob = (almph * bmuph - amuph * blmph)

call divide_by_cos2(jacob)        ! Get rid of the factor of cos2

end function jacob



  subroutine ftstpwv(lap, lapold, dt, deltat)
!-----------------------------------------------------------------------
! subroutine ftstpwv(lap, lapold, dt, deltat)

implicit none

! Does a single forward timestep

real,    intent(in)    :: deltat
complex, intent(in)    ::  dt(:, :)
complex, intent(inout) :: lap(:, :)
complex, intent(out)   :: lapold(:, :)

lapold = lap
lap    = lap + dt * deltat

end subroutine ftstpwv



  subroutine lfrogwv(lap, lapold, dt, deltat)
!----------------------------------------------------------------------
! subroutine lfrogwv(lap, lapold, dt, deltat)
!
!  does leapfrog timestepping
!

implicit none

complex, intent(inout) :: lap(:, :), lapold(:, :)
complex, intent(in)    :: dt(:, :)
real,    intent(in)    :: deltat

complex :: temp(size(lap, 1), size(lap, 2))

temp   = lapold + dt * 2 * deltat
lapold = lap
lap    = temp

end subroutine lfrogwv



  subroutine adams3(lap, lapold, dt, index, deltat)
!-----------------------------------------------------------------------
! subroutine adams3(lap, lapold, dt, index, deltat)
!
!  compute third order adams-bashforth timestep.  dt contains a 
!  rotating set of the last three time derivatives and index 
!  points to the element with the most recently updated.
!

implicit none

complex, intent(in)    :: dt(:, :, :)
integer, intent(in)    :: index
real,    intent(in)    :: deltat
complex, intent(inout) :: lap(:, :), lapold(:, :)

integer :: inm1, inm2, inm3

! Rotating storage for dt, last dimension must be 3

if(size(dt, 3) /= 3) then
   write(*, *) 'Error in subroutine adams3 in timstp; bad dt index'
   stop
endif

!  compute indices to the most recent, second and third most recent dts

inm1 = index
inm2 = (index + 1) - ((index + 1) / 3 * 3) + 1
inm3 = index - (index / 3 * 3) + 1

!  do the adams-bashforth step; see durran, 1991 for details

lapold = lap
lap    = lap + (deltat / 12.0) * (23.0 * dt(:, :, inm1) - 16.0 * &
                  dt(:, :, inm2) + 5.0 * dt(:, :, inm3))

end subroutine adams3



  subroutine model_output(x, time, field_name)
!-------------------------------------------------------------------------
! subroutine model_output(x, time, field_name)

!use ncd_file_mod

implicit none

integer,     intent(in) :: time
real,        intent(in) :: x(num_lon, lat_max)
character*4, intent(in) :: field_name

real, save        :: lons(num_lon), lats(lat_max)
character(len=30) :: my_axes(2)
real              :: valid_range(2)
logical, save     :: first = .true.
integer           :: start(2)

if(first) then
   write(*, *) 'ENTERING NCD_INIT LOOP'
   call get_deg_lon(lons)
   call get_deg_lat(lats)
   call def_axis('lon', lons, 'degrees_E', 'x', 'longitude', 1)
   call def_axis('lat', lats, 'degrees_N', 'y', 'latitude', 1)
   call def_time_axis('time', 'days',0.0, 'time extent', 1)
   call init_ncd_file('out.nc', 'time', 'stream function')
   my_axes(1) = 'lon'
   my_axes(2) = 'lat'
   valid_range(1) = -1e10
   valid_range(2) = 1e10
   call init_diag_field(field_name, 'out.nc', (/'lon','lat'/), &
   'percent', .FALSE., field_name, valid_range, -1.0, .FALSE., 1)
   first = .false.
end if

start = 1
call sum_diag(field_name, 'out.nc', x, start)
call output_diag('out.nc', 1.0*time)

end subroutine model_output



subroutine init_model()
!-------------------------------------------------------------------------
!
! For historical reasons, init_conditions does the model initialization
! for the barotropic model. This should be rewritten at some point.
! WARNING: If barotropic model is given a run-time resolution
! setting capability this will have to be changed.

end subroutine init_model




!-------------------------------------------------------------------------
! Following subroutines are for standard dynamical systems interface.
! WARNING: The dynamical systems routines all use real(r8),
! the barot model uses real.
!-------------------------------------------------------------------------



  subroutine init_conditions(x)
!-------------------------------------------------------------------------
! subroutine init_conditions(x)

implicit none

real(r8), intent(out) :: x(model_size)

complex, dimension(0:num_fourier, 0:num_spherical) :: psisp, force_in
complex :: temp
integer :: m, n
real    :: delt, dif_days_in

! Define the interesting indexes for variables to do diag output; span lats

do m = 1, 9
   diag_output_index(m) = (m - 1) * (lat_max / 9.0) + 1
end do

! Let's read in one of the old format files to act as an initial condition
!open(unit = 10, file = '/home/jla/psi/t21psijan80sp')
!open(unit = 11, file = '/home/jla/psi/t21psijan14sp')
!open(unit = 11, file = '/net/jla/t21psijan14sp')

! Unit 81 is used for real data runs;

open(unit = 81, file = '/t90/jla/assim/work/nov_to_mar')

! WARNING: Remember that old model uses TOTAL wavenumber

force_in = 0.0; psisp = 0.0
do n = 0, 21
   do m = 0, n
!      read(10, 31) temp
 31   format(1x, 2(e10.4, 1x))
!      if(n <= num_fourier) force_in(m, n - m) = temp
!      read(11, 31) temp
      read(81, 31) temp
      if(n <= num_fourier) psisp(m, n - m) = temp
   end do
end do

!close(unit = 11)

force_in = psisp

!  add in some del8 diffusion and forcing
!dif_days_in = 1.0
!dif_days_in = 0.95
!dif_days_in = 0.80
! Used for t21 and t42 reduced grid with forcing
dif_days_in = 2.0
! Used for t21 or t42 unforced for real data
!dif_days_in = 100.0


! Set timestep
! For t21 reduced grid
!delt = 3600.0
! For t42 reduced grid, t42 unforced or t21 unforced with real data
delt = 1800.0

!  initialize the triangular model; limits on fourier and spherical
call barot_init(dif_days_in, delt, force_in, 4, 10)

! Convert the psisp field for ics to format for dynamical systems
x = barot_to_dp(psisp)

end subroutine init_conditions



function get_model_size()
!-----------------------------------------------------------------------

implicit none

integer :: get_model_size

get_model_size = model_size

end function get_model_size



function dp_to_barot(x)
!-------------------------------------------------------------------------

implicit none

! Converts from real(r8) physical space array to barotropic 
! spherical harmonic stream function

real(r8), intent(in) :: x(model_size)
complex              :: dp_to_barot(0:num_fourier, 0:num_spherical)

real :: phys(num_lon, lat_max)

phys = dp_to_grid(x)      ! First copy into real

call trans_grid_to_spherical(phys, dp_to_barot)

end function dp_to_barot



  function dp_to_grid(x)
!-------------------------------------------------------------------------
! function dp_to_grid(x)

implicit none

! Converts from real(r8) physical space flat array to physical
! space array

real(r8), intent(in) :: x(model_size)
real                 :: dp_to_grid(num_lon, lat_max)

integer :: i, j, index

! First copy into real
do i = 1, num_lon
   do j = 1, lat_max
      index = j + (i - 1) * lat_max
      dp_to_grid(i, j) = real(x(index))
   end do
end do

end function dp_to_grid



  function barot_to_dp(psi)
!-------------------------------------------------------------------------
! function barot_to_dp(psi)
!
! Converts from spherical harmonic stream function to real(r8)
! physical space array.

implicit none


complex, intent(in) :: psi(0:num_fourier, 0:num_spherical)
real(r8)            :: barot_to_dp(model_size)

real    :: phys(num_lon, lat_max)
integer :: i, j, index

call trans_spherical_to_grid(psi, phys)    ! Transform to physical space

barot_to_dp = grid_to_dp(phys)             ! Copy into dp array

end function barot_to_dp



  function grid_to_dp(x)
!-------------------------------------------------------------------------
! function grid_to_dp(x)
!
! Converts from physical space array to flat double precison array
!

implicit none

real, intent(in) :: x(num_lon, lat_max)
real(r8)         :: grid_to_dp(model_size)

integer :: i, j, index

do i = 1, num_lon
   do j = 1, lat_max
      index = j + (i - 1) * lat_max
      grid_to_dp(index) = x(i, j)      ! Copy into dp array
   end do
end do

end function grid_to_dp



  subroutine adv_true_state(x)
!-------------------------------------------------------------------------
! subroutine adv_true_state(x)

implicit none

real(r8), intent(inout) :: x(model_size)

integer  :: m, n
real(r8) :: x_grid(num_lon, lat_max)
complex, dimension(0:num_fourier, 0:num_spherical) :: psisp

if(.not. use_real_data) then
   call adv_1step(x)
else

   ! Have real data every day

   real_time = real_time + delta_t * 48.0

   write(*, *) 'real time is ', real_time

   if(int(real_time / (24. * 3600.)) * (24. * 3600.)  == real_time) then
      write(*, *) 'updating observed state'

      ! WARNING: Remember that old model uses TOTAL wavenumber

      psisp = 0.0
      do n = 0, 21
         do m = 0, n
 31         format(1x, 2(e10.4, 1x))
            read(81, 31) psisp(m, n - m)
         end do
      end do
   
      write(*, *) 'Observed 4, 4 is ', psisp(4, 4)

      ! Next need to convert this to current model resolution single dimension state

      x_grid = dble(dp_to_grid(barot_to_dp(psisp)))

      ! RETURN X as FULL STATE SPACE STATE

      x = barot_to_dp(psisp)
   endif
endif

end subroutine adv_true_state



  subroutine adv_1step(x)
!-------------------------------------------------------------------------
! subroutine adv_1step(x)
!
! For now advances for a single 12-hour step

implicit none

real(r8), intent(inout) :: x(model_size)

complex :: psisp(0:num_fourier, 0:num_spherical)
integer :: i


psisp = dp_to_barot(x)      ! Convert to spectral, advance for 24 steps

! TEMPORARY KLUGE TO GET 24 hours for real data at T42: additional 24 steps

if(use_real_data) then
   psisp = forwrd(psisp, 48)
else
   psisp = forwrd(psisp, 24)
endif

x = barot_to_dp(psisp)      ! Convert back to grid dp

end subroutine adv_1step



  subroutine advance(x, nsteps, x_new)
!-------------------------------------------------------------------------
! subroutine advance(x, nsteps, x_new)

implicit none
integer,  intent(in)  :: nsteps
real(r8), intent(in)  :: x(model_size)
real(r8), intent(out) :: x_new(model_size)

integer :: i

x_new = x
do i = 1, nsteps
   call adv_1step(x_new)
end do

end subroutine advance



  subroutine filter(x)
!-------------------------------------------------------------------------
! subroutine filter(x)
!
! Does a spherical harmonic filtering truncation on a field

implicit none

real(r8), intent(inout) :: x(model_size)

integer, parameter :: filter_lim = 21
complex            :: psisp(0:num_fourier, 0:num_spherical)

psisp = dp_to_barot(x)

psisp(filter_lim-1:num_fourier, :)   = 0.0
psisp(:, filter_lim-1:num_spherical) = 0.0

x = barot_to_dp(psisp)        ! Convert back to grid dp

end subroutine filter


!subroutine get_close_pts(list, num)
!-------------------------------------------------------------------------
!subroutine get_close_pts(list, num)

! In the long run, this is too big for big models, will need another way?

!implicit none

!integer, intent(in) :: num
!integer, intent(out) :: list(model_size, num)
!real(r8):: dist(num)
!integer :: lon_id(num), lat_id(num)
!integer :: i, j, k, ii, jj, ind, ind2, j_lo, j_hi, base_index, index
!real(r8) :: tdist
!real :: rlat(lat_max), rlon(num_lon)
!type (loc_type) :: a, b


! Get lats and longs
!call get_deg_lat(rlat)
!call get_deg_lon(rlon)
!lat = dble(rlat)
!lon = dble(rlon)

! Get distances for each latitude row, lons are just uniform offset
!do j = 1, lat_max

! Initialize the temporary list
!   dist(:) = huge(dist)
!   lon_id = 0
!   lat_id = 0

!   a%lon = 0.0
!   a%lat = lat(j)
!   do ii = 0, num_lon - 1
! For efficiency, limit number of latitudes searched
!      j_lo = j - (sqrt(1.0 * num) / 2.0 + 1.0)
!      if(j_lo < 1) j_lo = 1
!      j_hi = j + (sqrt(1.0 * num) / 2.0 + 1.0)
!      if(j_hi > lat_max) j_hi = lat_max
      
!      do jj = j_lo, j_hi
!         b%lon = lon(ii + 1)
!         b%lat = lat(jj)
!         tdist = get_dist(a, b)
! Insert this distance into the list that holds num
!         do ind = 1, num
!            if(tdist < dist(ind)) then
!               do ind2 = num, ind + 1, -1
!                  dist(ind2) = dist(ind2 - 1)
!                  lon_id(ind2) = lon_id(ind2 - 1)
!                  lat_id(ind2) = lat_id(ind2 - 1)
!               end do
!               dist(ind) = tdist
!               lon_id(ind) = ii
!               lat_id(ind) = jj
!               goto 10
!            endif
!         end do
! 10   end do
!   end do
! Now load up close points for each lon point in this lat row
!   do i = 1, num_lon
!      base_index = j + (i - 1) * lat_max
!      write(*, *) 'base i, j, index ', i, j, base_index
!      do k = 1, num
!         ii = i + lon_id(k)
!         if(ii > num_lon) ii = ii - num_lon
!         jj = lat_id(k)
!         index = jj + (ii - 1) * lat_max
!!         write(*, *) 'neighbor ', k, ' lon lat ind ', ii, jj, index
!         list(base_index, k) = index
!      end do
!   end do
!end do

!end subroutine get_close_pts





subroutine model_get_close_states(o_loc, radius, number, indices, dist, x)
!--------------------------------------------------------------------
! 
! Stub for computation of get close states

implicit none

type(location_type), intent(in) :: o_loc
real(r8), intent(in) :: radius
integer, intent(out) :: number, indices(:)
real(r8), intent(out) :: dist(:)
real(r8), intent(in) :: x(:)

! Because of F90 limits this stub must be here telling assim_model
! to do exhaustive search (number = -1 return)
number = -1

end subroutine model_get_close_states



!-------------------------------------------------------------------------
! End of model_mod.f90
!-------------------------------------------------------------------------

end module model_mod
