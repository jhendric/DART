! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use        types_mod, only : r8
use time_manager_mod, only : time_type, set_time, get_time
use     location_mod, only : location_type, set_location, get_location, &
                             LocationDims, LocationName, LocationLName, &
                             get_close_maxdist_init, get_close_obs_init, get_close_obs

use    utilities_mod, only : register_module, error_handler, E_ERR, E_MSG, logfileunit, &
                             find_namelist_in_file, check_namelist_read

implicit none
private

public :: get_model_size, &
          adv_1step, &
          get_state_meta_data, &
          model_interpolate, &
          get_model_time_step, &
          end_model, &
          static_init_model, &
          init_time, &
          init_conditions, &
          model_get_close_states, &
          nc_write_model_atts, &
          nc_write_model_vars, &
          nc_read_model_vars, &
          pert_model_state, &
          get_close_maxdist_init, get_close_obs_init, get_close_obs, ens_mean_for_model

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!---------------------------------------------------------------
! Basic model parameters controlled by nameslist; have defaults
!
integer  :: model_size       = 960
real(r8) :: forcing          = 15.00_r8
real(r8) :: delta_t          = 0.001_r8
real(r8) :: space_time_scale = 10.00_r8
real(r8) :: coupling         = 3.00_r8
integer  :: K                = 32
integer  :: smooth_steps     = 16
integer  :: time_step_days    = 0
integer  :: time_step_seconds = 3600
integer  :: model_number      = 3 ! (2 for single scale, 3 for 2-scale, Lorenz 05)

namelist /model_nml/ model_size, forcing, delta_t, space_time_scale, coupling, K, &
       smooth_steps, time_step_days, time_step_seconds, model_number

!---------------------------------------------------------------- 
! Define the location of the state variables in module storage
type(location_type), allocatable :: state_loc(:)
type(time_type) :: time_step

! Define the averaging function for the production of x from z (L2k4)
real(r8), allocatable :: a(:)

! Define some parameters for computational efficiency
integer  :: H
integer  :: K2
integer  :: K4
integer  :: ss2
real(r8) :: sts2

contains

!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------
! Initializes class data for this model. For now, simply outputs the
! identity info, sets the location of the state variables, and initializes
! the time type for the time stepping (is this general enough for time???)

real(r8) :: x_loc
real(r8) :: ri
real(r8) :: alpha, beta
integer  :: i, iunit, io
integer  :: j

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the namelist entry
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the namelist values used for the run ...
call error_handler(E_MSG,'static_init_model','model_nml values are',' ',' ',' ')
write(logfileunit, nml=model_nml)
write(     *     , nml=model_nml)

! Create storage for locations
allocate(state_loc(model_size))

! Define the locations of the model state variables
do i = 1, model_size
   x_loc = (i - 1.0_r8) / model_size
   state_loc(i) =  set_location(x_loc)
end do

! The time_step in terms of a time type must also be initialized. Need
! to determine appropriate non-dimensionalization conversion for L96 from
! Shree Khare.
time_step = set_time(time_step_seconds, time_step_days)

! Generate the alpha and beta parameters for the calculation of "a"
alpha = (3.0_r8*(smooth_steps**2) + 3.0_r8) &
      / (2.0_r8*(smooth_steps**3) + 4.0_r8*smooth_steps)
beta  = (2.0_r8*(smooth_steps**2) + 1.0_r8) &
      / (1.0_r8*(smooth_steps**4) + 2.0_r8*(smooth_steps**2))

! The "a" vector is a smoothing filter for the production of x and y from z
! in L2k4. Apologies for the "ri" and "j" construct
allocate(a(2*smooth_steps + 1))
ri = - smooth_steps - 1.0_r8
j = 0
do i = - smooth_steps, smooth_steps
   j = j + 1
   ri = ri + 1.0_r8
   a(j) = alpha - beta*abs(ri)
end do
a(1)=a(1)/2.00_r8
a(2*smooth_steps+1)=a(2*smooth_steps+1)/2.00_r8

! defining parameters to help reduce the number of operations in the calculation
! of dz/dt
H    = K/2
K2   = 2*K
K4   = 4*K
ss2  = 2*smooth_steps
sts2 = space_time_scale**2   

end subroutine static_init_model


subroutine comp_dt(z, dt) 
!------------------------------------------------------------------
! subroutine comp_dt(z, dt)
! 
! Computes the time tendency of the lorenz 2004 model given current state.
!
! The model equations are given by
! 
! Model 2 (II)
!      dX_i
!      ---- = [X,X]_{K,i} -  X_i + F 
!       dt                
!                         
!
! Model 3 (III)
!      dZ_i
!      ---- = [X,X]_{K,i} + b^2 (-Y_{i-2}Y_{i-1} + Y_{i-1}Y_{i+1})
!       dt                +  c  (-Y_{i-2}X_{i-1} + Y_{i-1}X_{i+1})
!                         -  X_i - b Y_i + F,
!
! where
!
!     [X,X]_{K,i} = -W_{i-2K}W_{i-K} 
!                 +  sumprime_{j=-(K/2)}^{K/2} W_{i-K+j}X_{i+K+j}/K,
!
!      W_i =  sumprime_{j=-(K/2)}^{K/2} X_{i-j}/K,
!
! and sumprime denotes a special kind of summation where the first
! and last terms are divided by 2.
!
! NOTE: The equations above are only valid for K even.  If K is odd,
! then sumprime is replaced by the traditional sum, and the K/2 limits
! of summation are replaced by (K-1)/2. THIS CODE ONLY IMPLEMENTS THE
! K EVEN SOLUTION!!!
!
! The variable that is integrated is X (model II) or Z (model III), 
! but the integration of Z requires
! the variables X and Y.  For model III they are obtained by
!
!      X_i = sumprime_{j= -J}^{J} a_j Z_{i+j}
!      Y_i = Z_i - X_i.
!
! The "a" coefficients are given by
!
!      a_j = alpha - beta |j|,
! 
! where
!
!      alpha = (3J^2 + 3)/(2J^3 + 4J)
!      beta  = (2J^2 + 1)/(1J^4 + 2J^2).
!
! This choice of alpha and beta ensures that X_i will equal Z_i
! when Z_i varies quadratically over the interval 2J.   This choice
! of alpha and beta means that sumprime a_j = 1 and 
! sumprime (j^2) a_j = 0.
!
! Note that the impact of this filtering is to put large-scale
! variations into the X variable, and small-scale variations into
! the Y variable.
! 
! The parameter names above are based on those that appear in
! Lorenz 04.  To map to the code below, set:
!
!       F = forcing 
!       b = space_time_scale
!       c = coupling
!       K = K
!       J = smooth_steps

real(r8), intent( in)        ::  z(:)
real(r8), intent(out)        :: dt(:)
real(r8), dimension(size(z)) :: x, y
real(r8)                     :: xwrap(- K4:model_size + K4)
real(r8)                     :: ywrap(- K4:model_size + K4)
real(r8)                     ::    wx(- K4:model_size + K4)
real(r8)                     :: xx
integer                      :: i, j

! could branch this differently for more effecient model II

if ( model_number == 3 ) then
   ! Decompose z into x and y
   call z2xy(z,x,y)
elseif ( model_number == 2 ) then
   x = z
   y = 0.0_r8   ! just a dummy
else
   call error_handler(E_ERR,'comp_dt',&
         'Do not know that model number', source, revision, revdate)
endif

! Deal with cyclic boundary conditions using buffers
do i = 1, model_size
   xwrap(i) = x(i)
   ywrap(i) = y(i)
end do

! Fill the xwrap and ywrap buffers
do i = 1, K4
   xwrap(- K4 + i)       = xwrap(model_size - K4 +i)
   xwrap(model_size + i) = xwrap(i)
   ywrap(- K4 + i)       = ywrap(model_size - K4 +i)
   ywrap(model_size + i) = ywrap(i)
end do

! Calculate the W's
do i = 1, model_size
   wx(i) = xwrap(i - (-H))/2.00_r8
   do j = - H + 1, H - 1
      wx(i) = wx(i) + xwrap(i - j)
   end do
   wx(i) = wx(i) + xwrap(i - H)/2.00_r8
   wx(i) = wx(i)/K
end do

! Fill the W buffers
do i = 1, K4
   wx(- K4 + i)       = wx(model_size - K4 + i)
   wx(model_size + i) = wx(i)
end do

! Generate dz/dt
do i = 1, model_size
   xx = wx(i - K + (-H))*xwrap(i + K + (-H))/2.00_r8
   do j = - H + 1, H - 1
      xx = xx + wx(i - K + j)*xwrap(i + K + j)
   end do
   xx = xx + wx(i - K + H)*xwrap(i + K + H)/2.00_r8
   xx = - wx(i - K2)*wx(i - K) + xx/K
      
   if ( model_number == 3 ) then
     dt(i) = xx + (sts2)*( - ywrap(i - 2)*ywrap(i - 1) &
         + ywrap(i - 1)*ywrap(i + 1)) + coupling*( - ywrap(i - 2)*xwrap(i - 1) &
         + ywrap(i - 1)*xwrap(i + 1)) - xwrap(i) - space_time_scale*ywrap(i) &
         + forcing
   else ! must be model II
     dt(i) = xx - xwrap(i) + forcing
   endif

end do

end subroutine comp_dt


subroutine z2xy(z,x,y)
!------------------------------------------------------------------
! subroutine z2xy(z,x,y)
!
! Decomposes z into x and y for L2k4

integer :: i, j, ia
real(r8), intent( in) :: z(:)
real(r8), intent(out) :: x(:)
real(r8), intent(out) :: y(:)
real(r8)              :: zwrap(- ss2:model_size + ss2)

! Fill zwrap
do i = 1, model_size
   zwrap(i) = z(i)
end do
zwrap( - ss2) = zwrap(model_size - ss2)
do i = 1, ss2
   zwrap( - ss2 + i) = zwrap(model_size - ss2 + i)
   zwrap(model_size + i) = zwrap(i)
end do

! Generate the x variables
do i = 1, model_size
   ia = 1
   x(i) = a(ia)*zwrap(i - ( - smooth_steps))/2.00_r8
   do j = - smooth_steps + 1, smooth_steps - 1
      ia = ia + 1
      x(i) = x(i) + a(ia)*zwrap(i - j)
   end do
   ia = ia + 1
   x(i) = x(i) + a(ia)*zwrap(i - smooth_steps)/2.00_r8
end do

! Generate the y variables
do i = 1, model_size
   y(i) = z(i) - x(i)
end do

end subroutine z2xy


subroutine init_conditions(x)
!------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Initial conditions for lorenz 04
! It is assumed that this is called before any other routines in this
! module. Should probably make that more formal and perhaps enforce for
! more comprehensive models.


real(r8), intent(out) :: x(:)

x    = forcing
x(1) = 1.001_r8 * forcing

end subroutine init_conditions



subroutine adv_1step(x, time)
!------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
! Does single time step advance for lorenz 04 model
! using four-step rk time step
! The Time argument is needed for compatibility with more complex models
! that need to know the time to compute their time tendency and is not
! used in L04. Is there a better way to do this in F90 than to just hang
! this argument out everywhere?


real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time

real(r8), dimension(size(x)) :: x1, x2, x3, x4, dx, inter
real(r8), dimension(size(x)) :: dxt

call comp_dt(x, dx)    !  Compute the first intermediate step
x1    = delta_t * dx
inter = x + x1 / 2.0_r8

call comp_dt(inter, dx)!  Compute the second intermediate step
x2    = delta_t * dx
inter = x + x2 / 2.0_r8

call comp_dt(inter, dx)!  Compute the third intermediate step
x3    = delta_t * dx
inter = x + x3

call comp_dt(inter, dx)!  Compute fourth intermediate step
x4 = delta_t * dx

!  Compute new value for x

dxt = x1/6.0_r8 + x2/3.0_r8 + x3/3.0_r8 + x4/6.0_r8
x = x + dxt

end subroutine adv_1step



function get_model_size()
!------------------------------------------------------------------
! function get_model_size()
!
! Returns size of model

integer :: get_model_size

get_model_size = model_size

end function get_model_size



subroutine init_time(time)
!------------------------------------------------------------------
!
! Gets the initial time for a state from the model. Where should this info
! come from in the most general case?

type(time_type), intent(out) :: time

! For now, just set to 0
time = set_time(0, 0)

end subroutine init_time



subroutine model_interpolate(x, location, itype, obs_val, istatus)
!------------------------------------------------------------------
!
! Interpolates from state vector x to the location. It's not particularly
! happy dumping all of this straight into the model. Eventually some
! concept of a grid underlying models but above locations is going to
! be more general. May want to wait on external infrastructure projects
! for this?

! Argument itype is not used here because there is only one type of variable.
! Type is needed to allow swap consistency with more complex models.
   

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus
 
integer :: lower_index, upper_index, i
real(r8) :: lctn, lctnfrac

! All forward operators supported   
istatus = 0

! Convert location to real
lctn = get_location(location)
! Multiply by model size assuming domain is [0, 1] cyclic
lctn = model_size * lctn

lower_index = int(lctn) + 1
upper_index = lower_index + 1
if(lower_index > model_size) lower_index = lower_index - model_size
if(upper_index > model_size) upper_index = upper_index - model_size

lctnfrac = lctn - int(lctn)
obs_val = (1.0_r8 - lctnfrac) * x(lower_index) + lctnfrac * x(upper_index)

if(1 == 1) return

!!!obs_val = obs_val ** 2
!!!if(1 == 1) return

! Temporarily add on an observation from the other side of the domain, too
lower_index = lower_index + model_size / 2
if(lower_index > model_size) lower_index = lower_index - model_size
upper_index = upper_index + model_size / 2
if(upper_index > model_size) upper_index = upper_index - model_size
obs_val = obs_val + &
   lctnfrac * x(lower_index) + (1.0_r8 - lctnfrac) * x(upper_index)
if(1 == 1) return


! Next one does an average over a range of points
obs_val = 0.0_r8
lower_index = lower_index - 7
upper_index = upper_index - 7
if(lower_index < 1) lower_index = lower_index + model_size
if(upper_index < 1) upper_index = upper_index + model_size

do i = 1, 15
   if(lower_index > model_size) lower_index = lower_index - model_size
   if(upper_index > model_size) upper_index = upper_index - model_size
   obs_val = obs_val + &
      (1.0_r8 - lctnfrac) * x(lower_index) + lctnfrac * x(upper_index)
   lower_index = lower_index + 1
   upper_index = upper_index + 1
end do

end subroutine model_interpolate




function get_model_time_step()
!------------------------------------------------------------------
! function get_model_time_step()
!
! Returns the the time step of the model. In the long run should be repalced
! by a more general routine that returns details of a general time-stepping
! capability.

type(time_type) :: get_model_time_step

get_model_time_step = time_step

end function get_model_time_step



subroutine get_state_meta_data(index_in, location, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. This is not a function because the more general
! form of the call has a second intent(out) optional argument kind.
! Maybe a functional form should be added?


integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type                                      

location = state_loc(index_in)
if (present(var_type)) var_type = 1    ! default variable type

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Nothing for L04 for now.


end subroutine end_model


subroutine model_get_close_states(o_loc, radius, inum, indices, dist, x)
!------------------------------------------------------------------
!
! Stub for computation of get close states
   
type(location_type), intent( in) :: o_loc
real(r8),            intent( in) :: radius  
integer,             intent(out) :: inum
integer,             intent( in) :: indices(:)
real(r8),            intent( in) :: dist(:)
real(r8),            intent( in) :: x(:)

! Because of F90 limits this stub must be here telling assim_model
! to do exhaustive search (inum = -1 return)
inum = -1

end subroutine model_get_close_states


function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH Jan 24 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_04 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

use typeSizes
use netcdf

integer, intent(in)  :: ncFileID      ! netCDF file identifier
integer              :: ierr          ! return value of function

!--------------------------------------------------------------------
! General netCDF variables
!--------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID

!--------------------------------------------------------------------
! netCDF variables for Location
!--------------------------------------------------------------------

integer :: LocationVarID
integer :: StateVarDimID, StateVarVarID
integer :: StateVarID, MemberDimID, TimeDimID

!--------------------------------------------------------------------
! local variables
!--------------------------------------------------------------------

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer             :: i
type(location_type) :: lctn
ierr = 0                      ! assume normal termination

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!--------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(nf90_sync(ncFileID)) ! Ensure netCDF file is current
call check(nf90_Redef(ncFileID))

!--------------------------------------------------------------------
! Determine ID's from stuff already in the netCDF file
!--------------------------------------------------------------------

! make sure time is unlimited dimid

call check(nf90_inq_dimid(ncFileID,"copy",dimid=MemberDimID))
call check(nf90_inq_dimid(ncFileID,"time",dimid=TimeDimID))

!--------------------------------------------------------------------
! Write Global Attributes
!--------------------------------------------------------------------
call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source", source ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision", revision ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate", revdate ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model", "Lorenz_04"))
if ( model_number == 2 ) then
   call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_scale", "single"))
else if ( model_number == 3 ) then
   call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_scale", "2-scale"))
endif
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_forcing", forcing ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_delta_t", delta_t ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "space_time_scale", space_time_scale ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "coupling", coupling ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "K", K ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "smooth_steps", smooth_steps ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "time_step_days", time_step_days ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "time_step_seconds", time_step_seconds ))

!--------------------------------------------------------------------
! Define the model size, state variable dimension ... whatever ...
!--------------------------------------------------------------------

call check(nf90_def_dim(ncid=ncFileID, name="StateVariable", &
                        len=model_size, dimid = StateVarDimID))

!--------------------------------------------------------------------
! Define the Location Variable and add Attributes
! Some of the atts come from location_mod (via the USE: stmnt)
! CF standards for Locations:
! http://www.cgd.ucar.edu/cms/eaton/netcdf/CF-working.html#ctype
!--------------------------------------------------------------------

call check(NF90_def_var(ncFileID, name=trim(adjustl(LocationName)), xtype=nf90_double, &
              dimids = StateVarDimID, varid=LocationVarID) )
call check(nf90_put_att(ncFileID, LocationVarID, "long_name", trim(adjustl(LocationLName))))
call check(nf90_put_att(ncFileID, LocationVarID, "dimension", LocationDims ))
call check(nf90_put_att(ncFileID, LocationVarID, "units", "nondimensional"))
call check(nf90_put_att(ncFileID, LocationVarID, "valid_range", (/ 0.0_r8, 1.0_r8 /)))

!--------------------------------------------------------------------
! Define either the "state vector" variables -OR- the "prognostic" variables.
!--------------------------------------------------------------------

! Define the state vector coordinate variable
call check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
           dimids=StateVarDimID, varid=StateVarVarID))
call check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"))
call check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical") )
call check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, model_size /)))

! Define the actual state vector
call check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_double, &
           dimids = (/ StateVarDimID, MemberDimID, TimeDimID /), varid=StateVarID))
call check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"))

! Leave define mode so we can fill
call check(nf90_enddef(ncfileID))

! Fill the state variable coordinate variable
call check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ))

!--------------------------------------------------------------------
! Fill the location variable
!--------------------------------------------------------------------

do i = 1,model_size
   call get_state_meta_data(i,lctn)
   call check(nf90_put_var(ncFileID, LocationVarID, get_location(lctn), (/ i /) ))
enddo

!--------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!--------------------------------------------------------------------
call check(nf90_sync(ncFileID))

write (*,*)'Model attributes written, netCDF file synched ...'

contains
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)
      integer, intent ( in) :: istatus
      if(istatus /= nf90_noerr) call error_handler(E_ERR,'nc_write_model_atts',&
         trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check
end function nc_write_model_atts




function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH 23 May 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_04 model, each state variable is at a separate location.
! that's all the model-specific attributes I can think of ...
!
! assim_model_mod:init_diag_output uses information from the location_mod
!     to define the location dimension and variable ID. All we need to do
!     is query, verify, and fill ...
!
! Typical sequence for adding new dimensions,variables,attributes:
! NF90_OPEN             ! open existing netCDF dataset
!    NF90_redef         ! put into define mode
!    NF90_def_dim       ! define additional dimensions (if any)
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset

use typeSizes
use netcdf

integer,                intent(in) :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(in) :: statevec
integer,                intent(in) :: copyindex
integer,                intent(in) :: timeindex
integer                            :: ierr          ! return value of function

!-------------------------------------------------------------------------------
! General netCDF variables
!-------------------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID

!-------------------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------------------

ierr = 0                      ! assume normal termination

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!-------------------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

! no matter the value of "output_state_vector", we only do one thing.

call check(NF90_inq_varid(ncFileID, "state", StateVarID) )
call check(NF90_put_var(ncFileID, StateVarID, statevec,  &
             start=(/ 1, copyindex, timeindex /)))

! write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID))
! write (*,*)'netCDF file is synched ...'

contains
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)
   integer, intent ( in) :: istatus
      if(istatus /= nf90_noerr) call error_handler(E_ERR,'nc_write_model_vars',&
         trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check
end function nc_write_model_vars




function nc_read_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
! Reads the model-specific variables from a netCDF file

use typeSizes
use netcdf

integer,                intent(in)  :: ncFileID      ! netCDF file identifier
real(r8), dimension(:), intent(out) :: statevec
integer,                intent(in)  :: copyindex
integer,                intent(in)  :: timeindex
integer                             :: ierr          ! return value of function

!-------------------------------------------------------------------------------
! General netCDF variables
!-------------------------------------------------------------------------------

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: StateVarID

!-------------------------------------------------------------------------------
! local variables
!-------------------------------------------------------------------------------

ierr = 0                      ! assume normal termination

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!-------------------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))

! no matter the value of "output_state_vector", we only do one thing.

call check(NF90_inq_varid(ncFileID, "state", StateVarID) )
call check(NF90_get_var(ncFileID, StateVarID, statevec,  &
             start=(/ 1, copyindex, timeindex /)))

! write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID))
! write (*,*)'netCDF file is synched ...'

contains
   ! Internal subroutine - checks error status after each netcdf, prints
   !                       text message each time an error code is returned.
   subroutine check(istatus)
   integer, intent ( in) :: istatus
      if(istatus /= nf90_noerr) call error_handler(E_ERR,'nc_write_model_vars',&
         trim(nf90_strerror(istatus)), source, revision, revdate)
   end subroutine check

end function nc_read_model_vars



subroutine pert_model_state(state, pert_state, interf_provided)
!------------------------------------------------------------------
! subroutine pert_model_state(state, pert_state, interf_provided)
!
! Perturbs a model state for generating initial ensembles
! Returning interf_provided means go ahead and do this with uniform
! small independent perturbations.

real(r8), intent(in)  :: state(:)
real(r8), intent(in)  :: pert_state(:)
logical,  intent(out) :: interf_provided

interf_provided = .false.

end subroutine pert_model_state



subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! Not used in low-order models

real(r8), intent(in) :: ens_mean(:)

end subroutine ens_mean_for_model


!===================================================================
! End of model_mod
!===================================================================
end module model_mod