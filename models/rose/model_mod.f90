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

!-----------------------------------------------------------------------
!
!     Interface for ROSE model
!
!-----------------------------------------------------------------------

! DART Modules 
use        types_mod, only : r8, pi
use time_manager_mod, only : time_type,set_time,print_time
use     location_mod, only : location_type, set_location, get_location,&
                             query_location, get_dist  
use    utilities_mod, only : file_exist, open_file, close_file, &       
                             error_handler, E_ERR, E_MSG, E_WARN, logfileunit, &
                             find_namelist_in_file, check_namelist_read
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
! ROSE Modules
use params, only : nx, ny, nz, nbcon
 
!-----------------------------------------------------------------------

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
          pert_model_state, &
          !ROSE specific routines!
          model_type, &
          prog_var_to_vector, &
          vector_to_prog_var, &
          read_ROSE_restart, &
          update_ROSE_restart, &
          init_model_instance, &
          end_model_instance

!-----------------------------------------------------------------------

type model_type
  real(r8), pointer :: vars_3d(:,:,:,:)
end type model_type

!----------------------------------------------------------------------

! Global storage for describing ROSE model class
integer :: model_size 
type(time_type) :: Time_step_ROSE

! Arrays to store lats, lons
real(r8) :: lons(nx), lats(ny), levs(nz)

! Random sequence and init for pert_model_state
logical :: first_pert_call = .true.
type(random_seq_type)   :: random_seq

!----------------------------------------------------------------------
! Nameslist variables with default values follow
! Namelist variables for defining state vector, and default values
integer :: state_num_3d = 9             ! # of 3d fields to read from file
namelist /model_nml/ state_num_3d

integer :: ntime = 8       
logical :: output_prog_diag = .false.                          !NOT USED
character (len=50) :: input_dir = '../DAinput/'                !NOT USED
character (len=50) :: out_dir   = '../DAoutput/'               !NOT USED
character (len=30) :: ncep_file = 'nmc_lbc.02.nc'              !NOT USED
character (len=30) :: restart_file = 'NMC_SOC.day151_2002.dat' !NOT USED
real(kind=r8) :: h_tune = pi                                   !NOT USED  
real(kind=r8) :: z_tune = 1.0                                  !NOT USED
real(kind=r8) :: target_time = 0.125 ! e.g., 168.0 = 7 days * 24 [hr]        !NOT USED

namelist /rose_nml/ target_time, &
                    input_dir, out_dir, ncep_file, restart_file,&
                    output_prog_diag, &
                    h_tune, z_tune, &
                    ntime


integer, parameter :: TYPE_U0 = 0, TYPE_V0 = 1, TYPE_T0 = 2, & 
                      TYPE_U = 3, TYPE_V = 4, TYPE_T = 5, & 
                      TYPE_Q_H = 6, TYPE_Q_OH = 7, TYPE_Q_O = 8
!----------------------------------------------------------------------

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

contains

!==================================================================


subroutine static_init_model()
!------------------------------------------------------------------
!
! Called to do one time initialization of the model. As examples,
! might define information about the model size or model timestep.
! In models that require pre-computed static data, for instance
! spherical harmonic weights, these would also be computed here.
! Can be a NULL INTERFACE for the simplest models.
!

integer :: iunit, ierr, io, i, j, k
integer :: Time_step_seconds, Time_step_days
integer :: seconds_of_day = 86400
real(r8) :: d_lat, d_lon
real(r8) :: z_m
real(r8) :: dz = 2500.0_r8, zbot = 17500.0_r8
character(len=129) :: err_string, nml_string

! Read the namelist rose_nml from the file rose.nml
call find_namelist_in_file("rose.nml", "rose_nml", iunit)
read(iunit, nml = rose_nml, iostat = io)
call check_namelist_read(iunit, io, "rose_nml")

call error_handler(E_MSG,'static_init_model','rose_nml values are',source,revision,revdate)
write(logfileunit,nml=rose_nml)
write(     *     ,nml=rose_nml)


! Read the namelist entry for model_mod from file input.nml
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

call error_handler(E_MSG,'static_init_model','model_nml values are',source,revision,revdate)
write(logfileunit,nml=model_nml)
write(     *     ,nml=model_nml)

! Compute overall model size and put in global storage
model_size = nx * ny * nz * state_num_3d

! Set the model minimum time step from the namelist seconds and days input
! ROSE control variable "ntime" specify the number of steps per hour

Time_step_seconds = nint(3600.0_r8 / real(ntime))
Time_step_days = 0
write(*, *) 'time step secs days ' , Time_step_seconds, Time_step_days
if (Time_step_seconds > seconds_of_day) then
   Time_step_days = floor(real(Time_step_seconds) / real(seconds_of_day))
   Time_step_seconds = mod(Time_step_seconds, Time_step_days * seconds_of_day) 
endif
Time_step_ROSE = set_time(Time_step_seconds, Time_step_days)

call print_time(Time_step_ROSE)

! lon: long_name = "geographic longitude", units = "degrees" ;
d_lon = 360.0_r8/real(nx)
do i = 1, nx
   lons(i) = (i-1)*d_lon
enddo

! lat: long_name = "geographic latitude",  units = "degrees" ;
d_lat = 180.0_r8/real(ny)
do j = 1, ny
   lats(j) = -90 + d_lat/2. + (j-1)*d_lat
enddo

! Pressure levels ! from msetfix.f
! 1013 [mb] exp(-z/7km) with a fixed scale height (7km) 
! (2.5 km increment from 17.5 km)
do k = 1, nz
z_m = (k-1)*dz + zbot
levs(k) = 1013._r8*exp(-z_m/7.e3)   ![hPa]
enddo

end subroutine static_init_model


subroutine init_conditions(x)
!------------------------------------------------------------------
! subroutine init_conditions(x)
!
! Returns a model state vector, x, that is some sort of appropriate
! initial condition for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

real(r8), intent(out) :: x(:)


end subroutine init_conditions



subroutine adv_1step(x, time)
!------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to 
! compute a timestep, for instance for radiation compuations.
! This interface is only called if the namelist parameter
! async is set to 0 in perfect_model_obs of filter or if the 
! program integrate_model is to be used to advance the model
! state as a separate executable. If one of these options
! is not going to be used (the model will only be advanced as
! a separate model-specific executable), this can be a 
! NULL INTERFACE.

real(r8), intent(inout) :: x(:)
type(time_type), intent(in) :: time


end subroutine adv_1step



function get_model_size()
!------------------------------------------------------------------
!
! Returns the size of the model as an integer. Required for all
! applications.

integer :: get_model_size

get_model_size = model_size

end function get_model_size



subroutine init_time(time)
!------------------------------------------------------------------
!
! Companion interface to init_conditions. Returns a time that is somehow 
! appropriate for starting up a long integration of the model.
! At present, this is only used if the namelist parameter 
! start_from_restart is set to .false. in the program perfect_model_obs.
! If this option is not to be used in perfect_model_obs, or if no 
! synthetic data experiments using perfect_model_obs are planned, 
! this can be a NULL INTERFACE.

type(time_type), intent(out) :: time


end subroutine init_time



subroutine model_interpolate(x, location, itype, obs_val, istatus)
!------------------------------------------------------------------
!
! Given a state vector, a location, and a model state variable type,
! interpolates the state variable field`to that location and returns
! the value in obs_val. The istatus variable should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is a model specific integer that specifies the type of field (for
! instance temperature, zonal wind component, etc.). In low order
! models that have no notion of types of variables, this argument can
! be ingored. For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observerd), this can be a NULL INTERFACE.

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

! Default for successful return
istatus = 0

end subroutine model_interpolate



function get_model_time_step()
!------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.
!
! Limited by ROSE's fixed time step (also by 2-step restart file)

type(time_type) :: get_model_time_step

get_model_time_step =  Time_step_ROSE

end function get_model_time_step



subroutine get_state_meta_data(index_in, location, var_type)
!------------------------------------------------------------------
!
! Given an integer index into the state vector structure, returns the
! associated location. A second intent(out) optional argument kind
! can be returned if the model has more than one type of field (for
! instance temperature and zonal wind component). This interface is
! required for all filter applications as it is required for computing
! the distance between observations and state variables.

integer,             intent(in)  :: index_in
type(location_type), intent(out) :: location
integer,             intent(out), optional :: var_type


integer  :: indx, num_per_col, col_num, col_elem, lon_index,&
    & lat_index, lev_index
real(r8) :: lon, lat, lev
integer  :: local_var_type, var_type_temp

! Easier to compute with a 0 to size - 1 index
indx = index_in - 1

! Compute number of items per column
num_per_col = nz * state_num_3d 

! What column is this index in
col_num  = indx / num_per_col 
col_elem = indx - col_num * num_per_col

! What lon and lat index for this column
lon_index = col_num / ny      ! ny is number of ROSE latitude grid
lat_index = col_num - lon_index * ny

! Get actual lon lat values from static_init arrays ???
lon = lons(lon_index + 1)
lat = lats(lat_index + 1)

! Now figure out which beast in column this is
lev_index = col_elem / state_num_3d 
lev = levs(lev_index + 1)

! Find which var_type this element is
var_type_temp = mod(col_elem, state_num_3d )
if(var_type_temp == 0) then
  local_var_type = TYPE_U0
else if(var_type_temp == 1) then
  local_var_type = TYPE_V0
else if(var_type_temp == 2) then
  local_var_type = TYPE_T0
else if(var_type_temp == 3) then
  local_var_type = TYPE_U
else if(var_type_temp == 4) then
  local_var_type = TYPE_V
else if(var_type_temp == 5) then
  local_var_type = TYPE_T
else if(var_type_temp == 6) then
  local_var_type = TYPE_Q_H
else if(var_type_temp == 7) then 
  local_var_type = TYPE_Q_OH
else
  local_var_type = TYPE_Q_O
endif

write(*, '(1x,3(f6.2,1x),i3)') lon, lat, lev, local_var_type

location = set_location(lon, lat, lev, 2)  ! 2 == pressure (hPa)

! If the type is wanted, return it
if(present(var_type)) var_type = local_var_type

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storgae, etc.

end subroutine end_model



!subroutine model_get_close_states(o_loc, radius, nfound, indices, dist, x)
subroutine model_get_close_states(o_loc, radius, inum, indices, dist, x)
!------------------------------------------------------------------
! 
! Computes a list of model state variable indices that are within 
! distance radius of a given location, o_loc. The units of the radius
! and the metric for computing distances is defined by the location module
! that is in  use. The number of state variables that are within radius
! of o_loc is returned in inum. The indices of each of these is 
! returned in indices and the corresponding distance in dist. The model
! state is available in x because it is sometimes required to determine
! the distance (for instance, the current model surface pressure field
! is required to compute the location of state variables in a sigma
! vertical coordinate model). A model can choose to do no computation
! here and return a value of -1 in inum. If this happens, the filter
! will do a naive search through ALL state variables for close states.
! This can work fine in low-order models, but can be far too expensive
! in large models.

type(location_type), intent(in) :: o_loc
real(r8), intent(in) :: radius
integer, intent(out) :: inum, indices(:)
real(r8), intent(out) :: dist(:)
real(r8), intent(in) :: x(:)

type(location_type) :: s_loc
real(r8) :: loc_array(3), sloc_array(3), t_dist
integer,  allocatable :: lon_ind(:), lat_ind(:)
real(r8), allocatable :: close_dist(:)
integer  :: num_per_col, col_base_index
integer  :: i, j,  num, max_size


! Simplest interface just sets inum to -1 and returns
! inum = -1


loc_array = get_location(o_loc)

! Number found starts at 0
inum = 0
! Num of close horizontal grid points starts at 0, too
num = 0
! For now, just allocate enough space for all grid points, may want
! to make this smaller at some point for big models.
max_size = nx * ny  ! number_of_longitude * number_of_latitude
allocate(lon_ind(max_size), lat_ind(max_size), close_dist(max_size))
! Look for close grid points in horizontal only
call grid_close_states2(o_loc, lons, lats, nx, ny, radius, &
   num, lon_ind, lat_ind, close_dist)

! Compute size of grid storage for full levels
num_per_col = nz * state_num_3d 
! For vertical localization need the vertical pressure structure for this column
do i = 1, num
   col_base_index = ((lon_ind(i) - 1) * nx + lat_ind(i) - 1) * num_per_col
   do j = 1, num_per_col
      ! Added for vertical localization, 17 May, 2004
      call get_state_meta_data(col_base_index + j, s_loc) 
      t_dist = get_dist(s_loc, o_loc)

       if(t_dist < radius) then
         inum = inum + 1
         if(inum <= size(indices)) indices(inum) = col_base_index + j
         if(inum <= size(dist)) dist(inum) = t_dist
      endif

   enddo
enddo

deallocate(lon_ind, lat_ind, close_dist)

end subroutine model_get_close_states



function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! Writes the model-specific attributes to a netCDF file
! TJH Jan 24 2003
!
! TJH 29 July 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
!
! For the lorenz_96 model, each state variable is at a separate location.
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

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID


integer :: StateVarDimID
integer :: MemberDimID, TimeDimID
integer :: lonDimID, latDimID, levDimID
integer :: lonVarID, latVarID, levVarID
integer :: u1VarID, v1VarID, t1VarID, uVarID, vVarID, tVarID
integer :: qnHVarID, qnOHVarID, qnOVarID

character(len=129) :: errstring
character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

!------------------------------------------------------------------

ierr = 0  ! assume normal termination 

!--------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file
!--------------------------------------------------------------------

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(nf90_Redef(ncFileID))

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies
!-------------------------------------------------------------------------------
                                                                                                           
call check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID))
call check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid=  TimeDimID))
                                                                                                           
if ( TimeDimID /= unlimitedDimId ) then
   write(errstring,*)'Time Dimension ID ',TimeDimID,' should equal Unlimited Dimension ID',unlimitedDimID
   call error_handler(E_ERR,'nc_write_model_atts', errstring, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Write Global Attributes
!-------------------------------------------------------------------------------
                                                                                                           
call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date",str1))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source",source))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate",revdate))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","ROSE"))

!-------------------------------------------------------------------------------
! Define the dimensions IDs
!-------------------------------------------------------------------------------

call check(nf90_def_dim(ncid=ncFileID, name="lon",   len = nx,   dimid =   lonDimID))
call check(nf90_def_dim(ncid=ncFileID, name="lat",   len = ny,   dimid =   latDimID))
call check(nf90_def_dim(ncid=ncFileID, name="lev",   len = nz,   dimid =   levDimID))

!-------------------------------------------------------------------------------
! Create the (empty) Variables and the Attributes
!-------------------------------------------------------------------------------

call check(nf90_def_var(ncFileID, name="lon", xtype=nf90_double, &
                                                dimids=lonDimID, varid=lonVarID) ) 
call check(nf90_put_att(ncFileID, lonVarID, "long_name", "longitude"))
call check(nf90_put_att(ncFileID, lonVarID, "cartesian_axis", "X"))
call check(nf90_put_att(ncFileID, lonVarID, "units", "degrees_east")) ! check
call check(nf90_put_att(ncFileID, lonVarID, "valid_range", (/ 0.0_r8, 360.0_r8 /)))

call check(nf90_def_var(ncFileID, name="lat", xtype=nf90_double, &
                                                dimids=latDimID, varid=latVarID) ) 
call check(nf90_put_att(ncFileID, latVarID, "long_name", "latitude"))
call check(nf90_put_att(ncFileID, latVarID, "cartesian_axis", "Y"))
call check(nf90_put_att(ncFileID, latVarID, "units", "degrees_north")) !check
call check(nf90_put_att(ncFileID, latVarID, "valid_range", (/ -90.0_r8, 90.0_r8 /)))

call check(nf90_def_var(ncFileID, name="lev", xtype=nf90_double, &
                                                dimids=levDimID, varid=levVarID) ) 
call check(nf90_put_att(ncFileID, levVarID, "long_name", "level"))
call check(nf90_put_att(ncFileID, levVarID, "cartesian_axis", "Z"))
call check(nf90_put_att(ncFileID, levVarID, "units", "hPa"))
call check(nf90_put_att(ncFileID, levVarID, "positive", "down")) !check

!----------------------------------------------------------------------------
! Create attributes for the state vector
!----------------------------------------------------------------------------

call check(nf90_def_var(ncid=ncFileID, name="u1", xtype=nf90_real, &
       dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
       varid  = u1VarID))
call check(nf90_put_att(ncFileID, u1VarID, "long_name", "zonal wind component")) !check
call check(nf90_put_att(ncFileID, u1VarID, "units", "m/s"))
                                                                                                           
call check(nf90_def_var(ncid=ncFileID, name="v1", xtype=nf90_real, &
       dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
       varid  = v1VarID))
call check(nf90_put_att(ncFileID, v1VarID, "long_name", "meridional wind component")) !check
call check(nf90_put_att(ncFileID, v1VarID, "units", "m/s"))

call check(nf90_def_var(ncid=ncFileID, name="t1", xtype=nf90_real, &
       dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
       varid  = t1VarID))
call check(nf90_put_att(ncFileID, t1VarID, "long_name", "temperature")) !check
call check(nf90_put_att(ncFileID, t1VarID, "units", "degrees Kelvin"))


call check(nf90_def_var(ncid=ncFileID, name="u", xtype=nf90_real, &
       dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
       varid  = uVarID))
call check(nf90_put_att(ncFileID, uVarID, "long_name", "zonal wind component"))
call check(nf90_put_att(ncFileID, uVarID, "units", "m/s"))
                                                                                                           
call check(nf90_def_var(ncid=ncFileID, name="v", xtype=nf90_real, &
       dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
       varid  = vVarID))
call check(nf90_put_att(ncFileID, vVarID, "long_name", "meridional wind component"))
call check(nf90_put_att(ncFileID, vVarID, "units", "m/s"))

call check(nf90_def_var(ncid=ncFileID, name="t", xtype=nf90_real, &
       dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
       varid  = tVarID))
call check(nf90_put_att(ncFileID, tVarID, "long_name", "temperature"))
call check(nf90_put_att(ncFileID, tVarID, "units", "degrees Kelvin"))

call check(nf90_def_var(ncid=ncFileID, name="qnH", xtype=nf90_real, &
       dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
       varid  = qnHVarID))
call check(nf90_put_att(ncFileID, qnHVarID, "long_name", "mixing ratio H"))
call check(nf90_put_att(ncFileID, qnHVarID, "units", "?"))

call check(nf90_def_var(ncid=ncFileID, name="qnOH", xtype=nf90_real, &
       dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
       varid  = qnOHVarID))
call check(nf90_put_att(ncFileID, qnOHVarID, "long_name", "mixing ratio OH"))
call check(nf90_put_att(ncFileID, qnOHVarID, "units", "?"))

call check(nf90_def_var(ncid=ncFileID, name="qnO", xtype=nf90_real, &
       dimids = (/ levDimID, lonDimID, latDimID, MemberDimID, unlimitedDimID /), &
       varid  = qnOVarID))
call check(nf90_put_att(ncFileID, qnOVarID, "long_name", "mixing ratio O"))
call check(nf90_put_att(ncFileID, qnOVarID, "units", "?"))

call check(nf90_enddef(ncfileID))

!-------------------------------------------------------------------------------
! Fill the variables
!-------------------------------------------------------------------------------

call check(nf90_put_var(ncFileID, lonVarID, lons))
call check(nf90_put_var(ncFileID, latVarID, lats))
call check(nf90_put_var(ncFileID, levVarID, levs))

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call check(nf90_sync(ncFileID))
write (*,*)'nc_write_model_atts: netCDF file ',ncFileID,' is synched ...'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains
                                                                                                           
! Internal subroutine - checks error status after each netcdf, prints
!                       text message each time an error code is returned.

subroutine check(istatus)
    integer, intent (in) :: istatus
    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'nc_write_model_atts', &
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
! For the lorenz_96 model, each state variable is at a separate location.
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

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID
integer :: u1VarID, v1VarID, t1VarID, uVarID, vVarID, tVarID
integer :: qnHVarID, qnOHVarID, qnOVarID
type(model_type)                   :: Var

ierr = 0  ! assume normal termination

call check(nf90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))


call init_model_instance(Var)
call vector_to_prog_var(statevec, Var)

call check(NF90_inq_varid(ncFileID,  "u1",  u1VarID))
call check(nf90_put_var( ncFileID,  u1VarId, var%vars_3d(:,:,:, 1), &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
                                                                                                                                     
call check(NF90_inq_varid(ncFileID,  "v1",  v1VarID))
call check(nf90_put_var( ncFileID,  v1VarId, var%vars_3d(:,:,:, 2), &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))
                                                                                                                                     
call check(NF90_inq_varid(ncFileID,  "t1",  t1VarID))
call check(nf90_put_var( ncFileID,  t1VarId, var%vars_3d(:,:,:, 3), &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))

call check(NF90_inq_varid(ncFileID,  "u",  uVarID))
call check(nf90_put_var( ncFileID,  uVarId, var%vars_3d(:,:,:, 4), &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))

call check(NF90_inq_varid(ncFileID,  "v",  vVarID))
call check(nf90_put_var( ncFileID,  vVarId, var%vars_3d(:,:,:, 5), &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))

call check(NF90_inq_varid(ncFileID,  "t",  tVarID))
call check(nf90_put_var( ncFileID,  tVarId, var%vars_3d(:,:,:, 6), &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))

call check(NF90_inq_varid(ncFileID,  "qnH",  qnHVarID))
call check(nf90_put_var( ncFileID,  qnHVarId, var%vars_3d(:,:,:, 7), &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))

call check(NF90_inq_varid(ncFileID,  "qnOH",  qnOHVarID))
call check(nf90_put_var( ncFileID,  qnOHVarId, var%vars_3d(:,:,:, 8), &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))

call check(NF90_inq_varid(ncFileID,  "qnH",  qnOVarID))
call check(nf90_put_var( ncFileID,  qnOVarId, var%vars_3d(:,:,:, 9), &
           start=(/ 1, 1, 1, copyindex, timeindex /) ))

write (*,*)'Finished filling variables ...'
call check(nf90_sync(ncFileID))
write (*,*)'netCDF file is synched ...'

call end_model_instance(Var)   ! should avoid any memory leaking

contains
                                                                                                         
 ! Internal subroutine - checks error status after each netcdf, prints
 !                       text message each time an error code is returned.
 subroutine check(istatus)
  integer, intent ( in) :: istatus
    if(istatus /= nf90_noerr) call error_handler(E_ERR, 'nc_write_model_vars', &
          trim(nf90_strerror(istatus)), source, revision, revdate)
                                                                                                         
 end subroutine check

end function nc_write_model_vars


subroutine pert_model_state(state, pert_state, interf_provided)
!------------------------------------------------------------------
!
! Perturbs a model state for generating initial ensembles.
! The perturbed state is returned in pert_state.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the inter_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding an O(0.1) magnitude perturbation to each
! model state variable independently. The interf_provided argument
! should be returned as .true. if the model want to do its own
! perturbing of states.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

integer                 :: i, variable_type
type(location_type)     :: temp_loc

! An interface is provided
interf_provided = .true.

! If first call initialize random sequence
if(first_pert_call) then
   call init_random_seq(random_seq)
   first_pert_call = .false.
endif

do i = 1, get_model_size()
   call get_state_meta_data(i, temp_loc, variable_type)
   if(variable_type == TYPE_U .or. variable_type == TYPE_V .or. variable_type == TYPE_T .or. &
   variable_type == TYPE_U0 .or. variable_type == TYPE_V0 .or. variable_type == TYPE_T0) then
      pert_state(i) = random_gaussian(random_seq, state(i), 0.5_r8)
   else
      pert_state(i) = state(i)
   endif
end do

end subroutine pert_model_state


subroutine grid_close_states2(o_loc, lons, lats, nlon, nlat, radius, &
                   num, close_lon_ind, close_lat_ind, close_dist)
!=======================================================================
! subroutine grid_close_states2(o_loc, lons, lats, nlon, nlat, radius, &
!                  num, close_lon_ind, close_lat_ind, close_dist)
!
!
! Finds close state points from a particular grid; Just uses horizontal
! distance by setting pressure of state location to same as observation.

type(location_type), intent(in)    :: o_loc
integer,             intent(in)    :: nlon, nlat
real(r8),            intent(in)    :: lons(nlon), lats(nlat), radius
integer,             intent(inout) :: num
integer,             intent(inout) :: close_lon_ind(:), close_lat_ind(:)
real(r8),            intent(out)   :: close_dist(:)

real(r8) :: glat, glon, loc_array(3), o_lon, o_lat, o_lev
real(r8) :: gdist, diff, row_dist(nlon)
integer  :: blat_ind, blon_ind, i, j, lat_ind, lon_ind
integer  :: row_lon_ind(nlon), row_num
real(r8), parameter :: glev = 1.0_r8
type(location_type) :: loc

! Get the lat and lon from the loc
loc_array = get_location(o_loc)
o_lon = loc_array(1)
o_lat = loc_array(2)

! Get index to closest lat and lon for this observation
blat_ind = get_closest_lat_index(o_lat, lats, nlat)
!write(*, *) 'closest latitude in grid is ', blat_ind, lats(blat_ind)
blon_ind = get_closest_lon_index(o_lon, lons, nlon)
!write(*, *) 'closest longitude in grid is ', blon_ind, lons(blon_ind)

! Begin a search along the latitude axis in the positive direction
do lat_ind = blat_ind, nlat
   glat = lats(lat_ind)
   ! Take care of storage round-off
   if(glat < -90.0_r8) glat =  0.0_r8
   if(glat >  90.0_r8) glat = 90.0_r8

   ! Search all the contiguous close longitudes around the base longitude
   call lon_search(glat, glev, blon_ind, o_loc, radius, lons, &
      row_lon_ind, row_dist, row_num)

   ! If none are found, it's time to search in the negative latitude direction
   if(row_num == 0) goto 11

   ! Copy the points found in the row into summary storage
   close_lon_ind(num+1 : num+row_num) = row_lon_ind(1:row_num)
   close_lat_ind(num+1 : num+row_num) = lat_ind
   close_dist(num+1 : num+row_num) = row_dist(1:row_num)
   num = num + row_num
end do

! Search in the negative lat direction
11 continue
do lat_ind = blat_ind - 1, 1, -1
   glat = lats(lat_ind)
   ! Take care of storage round-off
   if(glat < -90.0_r8) glat =  0.0_r8
   if(glat >  90.0_r8) glat = 90.0_r8

   ! Search all the contiguous close longitudes around the base longitude
   call lon_search(glat, glev, blon_ind, o_loc, radius, lons, &
      row_lon_ind, row_dist, row_num)

   ! If none are found, it's time to give up
   if(row_num == 0) return

   ! Copy the points found in the row into summary storage
   close_lon_ind(num+1 : num+row_num) = row_lon_ind(1:row_num)
   close_lat_ind(num+1 : num+row_num) = lat_ind
   close_dist(num+1 : num+row_num) = row_dist(1:row_num)
   num = num + row_num
end do

end subroutine grid_close_states2

subroutine lon_search(glat, glev, blon_ind, o_loc, radius, lons, &
                      close_lon_ind, close_dist, num)
!=======================================================================
!  subroutine lon_search(glat, glev, blon_ind, o_loc, radius, lons, &
!                     close_lon_ind, close_dist, num)
!
! Given an observation location and radius and a latitude row from a grid,
! searches to find all longitude points in this row that are within radius
! of the observation location and returns their latitude index, longitude
! index, and the distance between them and the observation.

real(r8),            intent(in)  :: glat, glev, radius, lons(:)
integer,             intent(in)  :: blon_ind
type(location_type), intent(in)  :: o_loc
integer,             intent(out) :: close_lon_ind(:), num
real(r8),            intent(out) :: close_dist(:)

type(location_type) :: loc
integer  :: nlon, j, max_pos, lon_ind, which_vert
real(r8) :: glon, gdist, olev

! Total number found is 0 at start
num = 0
nlon = size(lons)

! Search as far as possible in the positive direction
do j = 0, nlon - 1
   max_pos = j
   lon_ind = blon_ind + j
   if(lon_ind > nlon) lon_ind = lon_ind - nlon
   glon = lons(lon_ind)

   ! Correct for longitude storage round-off
   if(glon > 360.0_r8) glon = 360.0_r8
   if(glon <   0.0_r8) glon =   0.0_r8

   ! Use same vertical "philosophy" as the existing location object.
   ! As of April, 2004 -- the vertical is (still) ignored in get_dist.
   which_vert = nint(query_location(o_loc))
   ! Only looking for horizontal distance here, use same level as obs.
   olev = query_location(o_loc, 'VLOC')  
   loc        = set_location(glon, glat, olev, which_vert)
   gdist      = get_dist(loc, o_loc)
   if(gdist <= radius) then
      num = num + 1
      close_lon_ind(num) = lon_ind
      close_dist(num) = gdist
      ! If radius is too far for closest longitude, 
      ! no need to search further or to search other side
   else if (j == 0) then
      return
   else
      ! Look in negative longitude offset direction next
      goto 21
   endif
end do
! Falling off end means the whole longitude circle has been searched; move along
return

! Search around the other way
21 continue
do j = 1, nlon - 1 - max_pos

   lon_ind = blon_ind - j
   if(lon_ind < 1) lon_ind = nlon + lon_ind
   glon = lons(lon_ind)

   ! Correct for longitude storage round-off
   if(glon > 360.0_r8) glon = 360.0_r8
   if(glon <   0.0_r8) glon =   0.0_r8

   ! Use same vertical "philosophy" as the existing location object.
   ! As of April, 2004 -- the vertical is (still) ignored in get_dist.
   which_vert = nint(query_location(o_loc))
   ! Only looking for horizontal distance here, use same level as obs.
   olev = query_location(o_loc, 'VLOC')  
   loc        = set_location(glon, glat, olev, which_vert)
   gdist      = get_dist(loc, o_loc)

   if(gdist <= radius) then
      num = num + 1
      close_lon_ind(num) = lon_ind
      close_dist(num) = gdist
   else
      ! No more longitudes in negative direction
      return
   endif
end do

end subroutine lon_search


function get_closest_lat_index(o_lat, lats, nlat)
!=======================================================================
! function get_closest_lat_index(o_lat, lats, nlat)
!

integer, intent(in) :: nlat
real(r8), intent(in) :: o_lat, lats(nlat)
integer :: get_closest_lat_index

real(r8) :: lat_bot, lat_top, lat_int, diff
integer :: lower_ind

! Find closest lat
lat_bot = lats(1)
lat_top = lats(nlat)
lat_int = lats(2) - lats(1)
if(o_lat <= lat_bot) then
   get_closest_lat_index = 1
else if(o_lat >= lat_top) then
   get_closest_lat_index = nlat
else
   diff = (o_lat - lat_bot) / lat_int
   lower_ind = int(diff) + 1
   if(diff - int(diff) < 0.5) then
      get_closest_lat_index = lower_ind
   else
      get_closest_lat_index = lower_ind + 1
   endif
endif

end function get_closest_lat_index



function get_closest_lon_index(o_lon, lons, nlon)
!=======================================================================
! function get_closest_lon_index(o_lon, lons, nlon)

integer, intent(in) :: nlon
real(r8), intent(in) :: o_lon, lons(nlon)
integer :: get_closest_lon_index

real(r8) :: diff, lon_bot, lon_top, lon_int
integer :: lower_ind, blon_ind

! Find closest longitude on grid to given longitude
lon_bot = lons(1)
lon_top = lons(nlon)
lon_int = lons(2) - lons(1)
if(o_lon <= lon_bot) then
   diff = (lon_bot - o_lon) / lon_int
   if(diff > 0.5) then
      get_closest_lon_index = nlon
   else
      get_closest_lon_index = 1
   end if
else if(o_lon >= lon_top) then
   diff = (o_lon - lon_top) / lon_int
   if(diff > 0.5) then
      get_closest_lon_index = 1
   else
      get_closest_lon_index = nlon
   end if
else
   diff = (o_lon - lon_bot) / lon_int
   lower_ind = int(diff) + 1
   if(diff - int(diff) < 0.5) then
      get_closest_lon_index = lower_ind
   else
      get_closest_lon_index = lower_ind + 1
   end if
end if

end function get_closest_lon_index


subroutine prog_var_to_vector(var, x)
!=======================================================================
! based on prog_var_to_vector in model_mod for CAM.
!

type(model_type), intent(in) :: var
real(r8), intent(out) :: x(:)

integer :: i, j, k, nf, indx
character(len=129) :: errstring

! Do order as ps, t, u, v, q, tracers to be consistent with b-grid

! Start copying fields to straight vector
indx = 0
do i = 1, nx !longitude
   do j = 1, ny !latitude
! u,v,t,q, and tracers at successively lower levels
      do k = 1, nz ! height
         do nf= 1, state_num_3d
            indx = indx + 1
            x(indx) = var%vars_3d(k, i, j, nf)
         end do
      end do
   end do
end do

! Temporary check
if(indx /= model_size) then
   write(errstring, *) 'indx ',indx,' model_size ',model_size,' must be equal '
   call error_handler(E_ERR, 'prog_var_to_vector', errstring, source, revision, revdate)
endif

end subroutine prog_var_to_vector


subroutine vector_to_prog_var(x, var) 
!=======================================================================
! subroutine vector_to_prog_var(x, var) 
!

real(r8), intent(in) :: x(:)
type(model_type), intent(out) :: var

integer :: i, j, k, nf, indx
character(len=129) :: errstring

! Start copying fields from straight vector
indx = 0
do i = 1, nx ! longitude
   do j = 1, ny ! latitude
! u,v,t,q  and tracers at successive levels
      do k = 1, nz ! height
         do nf = 1, state_num_3d
            indx = indx + 1
            var%vars_3d(k, i, j, nf) = x(indx)
         end do 
      end do
   end do
end do

! Temporary check
if(indx /= model_size) then
   write(errstring, *) 'indx ',indx,' model_size ',model_size,' must be equal '
   call error_handler(E_ERR, 'vector_to_prog_var', errstring, source, revision, revdate)
endif

end subroutine vector_to_prog_var

subroutine init_model_instance(var)
!=======================================================================
! subroutine init_model_instance(var)
!
! Initializes an instance of a ROSE model state variable

type(model_type), intent(out) :: var

! Initialize the storage space and return

allocate(var%vars_3d(nz, nx, ny, state_num_3d))

end subroutine init_model_instance

subroutine end_model_instance(var)
!=======================================================================
! subroutine end_model_instance(var)
!
! Ends an instance of a rose model state variable
                                                                                                         
type(model_type), intent(inout) :: var
                                                                                                         
deallocate(var%vars_3d)
                                                                                                         
end subroutine end_model_instance

subroutine update_ROSE_restart(file_name, var)
!=======================================================================
! update ROSE restart file fields
!
  character (len = *), intent(in) :: file_name
  type(model_type), intent(in) :: var
  integer :: iunit, ierr, iftype

! ROSE restart file prognostic variables
! from a3d.chem.f
  integer :: doy                   ! day of year
  integer :: utsec                 ! universal time of model in seconds
  integer :: iyear                 ! number of years model has run
! from params.mod.f
  integer :: year0                 ! start year
  integer :: day0                  ! start day of year (day 1 = Jan. 1)
  integer :: ut0                   ! start universal time (sec) 
! from dynam.mod.f
  real, dimension (nz,nx,ny) :: un1, vn1, tn1, un0, vn0, tn0
  real, dimension (nz)       :: tref                                 
  real                       :: treflb, trefub
! from chem.mod.f ... constituent mixing ratios
  real, dimension (nz,nx,ny,nbcon) :: qn1
  real, dimension (nz,ny) :: q_o2, q_n2

  character(len=129) :: errstring
!
  integer :: dummy_rose(3)
  integer :: daynum
  real    :: gmt_frac        ! fraction of day since 0Z
!====================================================================

if(file_exist(file_name)) then

   ! Determine if we have a style 1 (old school) or not restart file.
   iunit = open_file(file_name, form = 'unformatted', action = 'read')
   read(iunit, iostat=ierr) iftype
   if ( ierr /= 0 ) then
      write(errstring,*)'restart file is ',trim(adjustl(file_name))
      call error_handler(E_MSG,'update_ROSE_restart',errstring,source,revision,revdate)
      call error_handler(E_ERR,'update_ROSE_restart','unknown restart file type',&
           source,revision,revdate)
   endif

   if (iftype == 1) then ! style 1 (old school) restart file.

      read(iunit, iostat=ierr) dummy_rose, gmt_frac, daynum, tref, treflb, &
                  trefub, un1, vn1, tn1, un0, vn0, tn0, &
                  qn1(:,:,:,1:nbcon-1), q_o2, q_n2

      if ( ierr /= 0 ) then
         write(errstring,*)'type 1 restart file is ',trim(adjustl(file_name))
         call error_handler(E_MSG,'update_ROSE_restart',errstring,source,revision,revdate)
         call error_handler(E_ERR,'update_ROSE_restart','read error',&
              source,revision,revdate)
      endif

      doy   = mod(dummy_rose(3), 365) 
      utsec = int(gmt_frac * 24.0 * 3600.0)

   else ! anything else is 'newish' 

      read(iunit, iostat=ierr) iyear, doy, utsec, year0, day0, ut0, tref, & 
                  treflb, trefub, un1, vn1, tn1, un0, vn0, tn0, &
                  qn1, q_o2, q_n2
      if ( ierr /= 0 ) then
         write(errstring,*)'type ',iunit,' restart file is ',trim(adjustl(file_name))
         call error_handler(E_MSG,'update_ROSE_restart',errstring,source,revision,revdate)
         call error_handler(E_ERR,'update_ROSE_restart','read error',&
              source,revision,revdate)
      endif

   endif

   call close_file(iunit)
else
   write(errstring,*) trim(adjustl(file_name)),' not available.'
   call error_handler(E_ERR,'update_ROSE_restart',errstring,source,revision,revdate)
endif

un1 = var%vars_3d(:,:,:,1)
vn1 = var%vars_3d(:,:,:,2)
tn1 = var%vars_3d(:,:,:,3)
un0 = var%vars_3d(:,:,:, 4)
vn0 = var%vars_3d(:,:,:, 5)
tn0 = var%vars_3d(:,:,:, 6)
qn1(:,:,:,7)  = var%vars_3d(:,:,:,7) ! H
qn1(:,:,:,8)  = var%vars_3d(:,:,:,8) ! OH
qn1(:,:,:,18) = var%vars_3d(:,:,:,9) ! O


iunit = open_file(file_name, form = 'unformatted', action = 'write')

if (iftype == 1) then ! style 1 (old school) restart file.
  write(iunit) iftype
  write(iunit) dummy_rose, gmt_frac, daynum, tref, treflb, &
                  trefub, un1, vn1, tn1, un0, vn0, tn0, &
                  qn1(:,:,:,1:nbcon-1), q_o2, q_n2
else
  write(iunit) iftype
  write(iunit) iyear, doy, utsec, year0, day0, ut0, &
             tref, treflb, trefub, un1, vn1, tn1, &
             un0, vn0, tn0, qn1, q_o2, q_n2
endif

call close_file(iunit)

end subroutine update_ROSE_restart


subroutine read_ROSE_restart(file_name, var, model_time)
!=======================================================================
! read ROSE restart file fields that have been updated
!
  character (len = *), intent(in) :: file_name
  type(model_type),   intent(out) :: var
  type(time_type),   intent(out)  :: model_time
  integer :: iunit, ierr, iftype
  integer :: model_seconds, model_days

! ROSE restart file prognostic variables
! from a3d.chem.f
  integer :: doy                   ! day of year
  integer :: utsec                 ! universal time of model in seconds
  integer :: iyear                 ! number of years model has run
! from params.mod.f
  integer :: year0                 ! start year
  integer :: day0                  ! start day of year (day 1 = Jan. 1)
  integer :: ut0                   ! start universal time (sec) 
! from dynam.mod.f
  real, dimension (nz,nx,ny) :: un1, vn1, tn1, un0, vn0, tn0
  real, dimension (nz)       :: tref                                 
  real                       :: treflb, trefub
! from chem.mod.f ... constituent mixing ratios
  real, dimension (nz,nx,ny,nbcon) :: qn1
  real, dimension (nz,ny)          :: q_o2, q_n2
!
  integer :: dummy_rose(3)
  integer :: daynum
  real    :: gmt_frac        ! fraction of day since 0Z

  character(len=129) :: errstring
!
!====================================================================

if(file_exist(file_name)) then

   iunit = open_file(file_name, form = 'unformatted', action = 'read')

   ! Determine if we have a style 1 (old school) or not restart file.
   read(iunit, iostat=ierr) iftype
   if ( ierr /= 0 ) then
      write(errstring,*)'restart file is ',trim(adjustl(file_name))
      call error_handler(E_MSG,'read_ROSE_restart',errstring,source,revision,revdate)
      call error_handler(E_ERR,'read_ROSE_restart','unknown restart file type',&
           source,revision,revdate)
   endif

   if (iftype == 1) then ! style 1 (old school) restart file.

      read(iunit, iostat=ierr) dummy_rose, gmt_frac, daynum, tref, treflb, &
                  trefub, un1, vn1, tn1, un0, vn0, tn0, &
                  qn1(:,:,:,1:nbcon-1), q_o2, q_n2

      if ( ierr /= 0 ) then
         write(errstring,*)'type 1 restart file is ',trim(adjustl(file_name))
         call error_handler(E_MSG,'read_ROSE_restart',errstring,source,revision,revdate)
         call error_handler(E_ERR,'read_ROSE_restart','read error',&
              source,revision,revdate)
      endif
      
      doy   = mod(dummy_rose(3), 365) 
      utsec = int(gmt_frac * 24.0 * 3600.)
                                                                                                              
   else ! anything else is 'newish' 

      read(iunit, iostat=ierr) iyear, doy, utsec, year0, day0, ut0, tref, & 
                  treflb, trefub, un1, vn1, tn1, un0, vn0, tn0, &
                  qn1, q_o2, q_n2
      if ( ierr /= 0 ) then
         write(errstring,*)'type ', iftype ,' restart file is ',trim(adjustl(file_name))
         call error_handler(E_MSG,'read_ROSE_restart',errstring,source,revision,revdate)
         call error_handler(E_ERR,'read_ROSE_restart','read error',&
              source,revision,revdate)
      endif

   endif

   call close_file(iunit)
else
   call error_handler(E_ERR,'read_ROSE_restart','rose_restart.dat not available',source,revision,revdate)
endif

var%vars_3d(:,:,:,1)  = un1 
var%vars_3d(:,:,:,2)  = vn1 
var%vars_3d(:,:,:,3)  = tn1 
var%vars_3d(:,:,:,4) = un0           ! at model_time
var%vars_3d(:,:,:,5) = vn0           ! at model_time
var%vars_3d(:,:,:,6) = tn0           ! at model_time
var%vars_3d(:,:,:,7) = qn1(:,:,:,7)  ! H
var%vars_3d(:,:,:,8) = qn1(:,:,:,8)  ! OH
var%vars_3d(:,:,:,9) = qn1(:,:,:,18) ! O

print*, "read_ROSE_restart: BEFORE model_time: utsec = ", utsec
print*, "read_ROSE_restart: BEFORE model_time: doy = ", doy
model_time = set_time(utsec, doy)
print*, "read_ROSE_restart: AFTER model_time :"
call print_time(model_time)


end subroutine read_ROSE_restart


!===================================================================
! End of model_mod
!===================================================================
end module model_mod
