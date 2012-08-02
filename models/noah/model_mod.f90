! DART software - Copyright 2004 - 2011 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

module model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

! This is a noah_1d showing the interfaces required for a model to be compliant
! with the DART data assimilation infrastructure. The public interfaces listed
! must all be supported with the argument lists as indicated. Many of the interfaces
! are not required for minimal implementation (see the discussion of each
! interface and look for NULL INTERFACE). 

! Modules that are absolutely required for use are listed
use        types_mod, only : r8, MISSING_R8, obstypelength
use time_manager_mod, only : time_type, set_time, set_date, get_time,          &
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=)

use     location_mod, only : location_type, get_dist, query_location,          &
                             get_close_maxdist_init, get_close_type,           &
                             set_location, get_location, horiz_dist_only,      &
                             vert_is_undef,    VERTISUNDEF,                    &
                             vert_is_surface,  VERTISSURFACE,                  &
                             vert_is_level,    VERTISLEVEL,                    &
                             vert_is_pressure, VERTISPRESSURE,                 &
                             vert_is_height,   VERTISHEIGHT,                   &
                             get_close_obs_init, get_close_obs,                &
                             set_location_missing, write_location

use    utilities_mod, only : register_module, error_handler, nc_check,         &
                             E_ERR, E_MSG, logfileunit, get_unit,              &
                             nmlfileunit, do_output, do_nml_file, do_nml_term, &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, file_exist, find_textfile_dims,        &
                             file_to_text

use     obs_kind_mod, only : KIND_SOIL_TEMPERATURE,   &
                             KIND_LIQUID_WATER,       &
                             KIND_ICE,                &
                             KIND_SNOWCOVER_FRAC,     &
                             KIND_SNOW_THICKNESS,     &
                             KIND_LEAF_CARBON,        &
                             KIND_WATER_TABLE_DEPTH,  &
                             paramname_length,        &
                             get_raw_obs_kind_index

use mpi_utilities_mod, only: my_task_id
use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use typesizes
use netcdf

implicit none
private

! required by DART code - will be called from filter and other
! DART executables.  interfaces to these routines are fixed and
! cannot be changed in any way.
public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          end_model,              &
          static_init_model,      &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_state,       &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          ens_mean_for_model

! not required by DART but for larger models can be useful for
! utility programs that are tightly tied to the other parts of
! the model_mod code.
public :: noah1d_to_dart_vector, &
          dart_vector_to_model_file, &
          get_noah1D_restart_filename

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

!------------------------------------------------------------------
! The variables in the noah restart file that are used to create the
! DART state vector are specified in the input.nml:model_nml namelist.
!
!    noah_variables  = 'STC',    'KIND_SOIL_TEMPERATURE',
!                      'SMC',    'KIND_SOIL_MOISTURE',
!                      'SH2O',   'KIND_LIQUID_SOIL_MOISTURE',
!                      'T1',     'KIND_SKIN_TEMPERATURE',
!                      'SNOWH',  'KIND_SNOW_DEPTH',
!                      'SNEQV',  'KIND_LIQUID_EQUIVALENT',
!                      'CMC',    'KIND_CANOPY_WATER',
!------------------------------------------------------------------

integer :: nfields
integer, parameter :: max_state_variables = 40
integer, parameter :: num_state_table_columns = 2
character(len=obstypelength) :: variable_table(max_state_variables, num_state_table_columns)

!------------------------------------------------------------------
! things which can/should be in the DART model_nml
!------------------------------------------------------------------

character(len=128)    :: noah_netcdf_filename   = 'OUTPUT.NC'
character(len=128)    :: noah_namelist_filename = 'somelocation.dat'
integer               :: assimilation_period_days     = 0
integer               :: assimilation_period_seconds  = 60
real(r8)              :: model_perturbation_amplitude = 0.2
logical               :: output_state_vector          = .true.
character(len=32)     :: calendar = 'Gregorian'
integer               :: debug    = 0  ! turn up for more and more debug messages
character(len=obstypelength) :: noah_variables(max_state_variables*num_state_table_columns) = ' '

namelist /model_nml/ noah_netcdf_filename, noah_namelist_filename, &
          assimilation_period_days, assimilation_period_seconds,   &
          model_perturbation_amplitude, output_state_vector,       &
          calendar, debug, noah_variables

!------------------------------------------------------------------
! Everything needed to recreate the NOAH METADTA_NAMELIST
!
! To restart the file, we write a new namelist.
! DART needs to write a NOAH-compatible namelist. 
!------------------------------------------------------------------

integer, parameter :: nSoilLayers = 4

character(len=12) :: startdate
character(len=12) :: enddate
logical  :: loop_for_a_while
real(r8) :: Latitude
real(r8) :: Longitude
integer  :: Forcing_Timestep
integer  :: Noahlsm_Timestep
logical  :: Sea_ice_point
real(r8), dimension(nSoilLayers) :: Soil_layer_thickness
real(r8), dimension(nSoilLayers) :: Soil_Temperature
real(r8), dimension(nSoilLayers) :: Soil_Moisture
real(r8), dimension(nSoilLayers) :: Soil_Liquid
real(r8) :: Skin_Temperature
real(r8) :: Canopy_water
real(r8) :: Snow_depth
real(r8) :: Snow_equivalent
real(r8) :: Deep_Soil_Temperature
character(len=256) :: Landuse_dataset
integer  :: Soil_type_index
integer  :: Vegetation_type_index
integer  :: Urban_veg_category
integer  :: glacial_veg_category
integer  :: Slope_type_index
real(r8) :: Max_snow_albedo
real(r8) :: Air_temperature_level
real(r8) :: Wind_level
real(r8) :: Green_Vegetation_Min
real(r8) :: Green_Vegetation_Max
logical  :: Usemonalb
logical  :: Rdlai2d
integer  :: sfcdif_option
integer  :: iz0tlnd
real(r8), dimension(12) :: Albedo_monthly
real(r8), dimension(12) :: Shdfac_monthly
real(r8), dimension(12) ::    lai_monthly
real(r8), dimension(12) ::  Z0brd_monthly

namelist /METADATA_NAMELIST/ startdate, enddate, loop_for_a_while,   &
         Latitude, Longitude, Forcing_Timestep, Noahlsm_Timestep,    &
         Sea_ice_point, Soil_layer_thickness, Soil_Temperature,      &
         Soil_Moisture, Soil_Liquid, Skin_Temperature, Canopy_water, &
         Snow_depth, Snow_equivalent, Deep_Soil_Temperature, Landuse_dataset, &
         Soil_type_index, Vegetation_type_index, Urban_veg_category, &
         glacial_veg_category, Slope_type_index, Max_snow_albedo,    &
         Air_temperature_level, Wind_level, Green_Vegetation_Min,    &
         Green_Vegetation_Max, Usemonalb, Rdlai2d, sfcdif_option,    &
         iz0tlnd, Albedo_monthly, Shdfac_monthly, lai_monthly, Z0brd_monthly

! We are going to create a DART state vector out of the following 17 items

type noahtype
   private
   real(r8), dimension(nSoilLayers) :: Soil_Temperature
   real(r8), dimension(nSoilLayers) :: Soil_Moisture
   real(r8), dimension(nSoilLayers) :: Soil_Liquid
   real(r8) :: Skin_Temperature
   real(r8) :: Canopy_water
   real(r8) :: Snow_depth
   real(r8) :: Snow_equivalent
   real(r8) :: Deep_Soil_Temperature
end type noahtype

! define model parameters here
type(time_type)     :: time_step
type(location_type) :: state_loc(0:nSoilLayers)
type(noahtype)      :: noah1d

! Everything needed to describe a variable

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=obstypelength), dimension(NF90_MAX_VAR_DIMS) :: dimnames
   integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
   integer :: numdims
   integer :: maxlevels
   integer :: xtype
   integer :: varsize     ! prod(dimlens(1:numdims))
   integer :: index1      ! location in dart state vector of first occurrence
   integer :: indexN      ! location in dart state vector of last  occurrence
   integer :: dart_kind
   character(len=paramname_length) :: kind_string
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

!------------------------------------------------------------------------------
! These are the metadata arrays that are the same size as the state vector.

real(r8), allocatable, dimension(:) :: ens_mean ! may be needed for forward ops
real(r8), allocatable, dimension(:) :: levels   ! depth

!------------------------------------------------------------------
! module storage
!------------------------------------------------------------------

integer            :: model_size       ! the state vector length
type(time_type)    :: model_time       ! valid time of the model state
type(time_type)    :: model_time_step  ! smallest time to adv model
character(len=256) :: string1, string2, string3
logical, save      :: module_initialized = .false.
real(r8), dimension(nSoilLayers) :: soil_depths

!==================================================================
contains
!==================================================================



subroutine static_init_model()
!------------------------------------------------------------------
! one time initialization of the model

! Local variables - all the important ones have module scope

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
character(len=NF90_MAX_NAME)          :: varname
character(len=obstypelength)          :: dimname
character(len=paramname_length)       :: kind_string

integer  :: VarID, dimlen, varsize 
integer  :: iunit, io, ivar
integer  :: i, index1, nLayers

if ( module_initialized ) return ! only need to do this once.

! Since this routine calls other routines that could call this routine
! we'll say we've been initialized pretty dang early.
module_initialized = .true.

! Print module information to log file and stdout.
call register_module(source, revision, revdate)

! Read the DART namelist
call find_namelist_in_file("input.nml", "model_nml", iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, "model_nml")

! Record the DART namelist values used for the run ...
if (do_nml_file()) write(nmlfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

! Read the NOAH namelist
call find_namelist_in_file(trim(noah_namelist_filename), "METADATA_NAMELIST", iunit)
read(iunit, nml = METADATA_NAMELIST, iostat = io)
call check_namelist_read(iunit, io, "METADATA_NAMELIST")

! Record the NOAH namelist
if (do_nml_file()) write(nmlfileunit, nml=METADATA_NAMELIST)
if (do_nml_term()) write(     *     , nml=METADATA_NAMELIST)

! Check to make sure the required input files exist
if ( .not. file_exist(noah_netcdf_filename) ) then
   write(string1,*) 'cannot open file ', trim(noah_netcdf_filename),' for reading.'
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif
if ( .not. file_exist(noah_namelist_filename) ) then
   write(string1,*) 'cannot open file ', trim(noah_namelist_filename),' for reading.'
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif

! convert the [-180,180] longitudes to [0,360)

if (Longitude < 0.0_r8) Longitude = Longitude + 360.0_r8

! The time_step in terms of a time type must also be initialized.

call set_calendar_type( calendar )

call nc_check(nf90_open(adjustl(noah_netcdf_filename), NF90_NOWRITE, iunit), &
                   'static_init_model', 'open '//trim(noah_netcdf_filename))

model_time      = get_state_time(iunit, trim(noah_netcdf_filename))
model_time_step = set_time(assimilation_period_seconds, assimilation_period_days)

if (debug > 0) then
   call print_date(model_time     ,'static_init_model:model date')
   call print_time(model_time     ,'static_init_model:model time')
   call print_time(model_time_step,'static_init_model:model timestep')
endif

! Make sure the number of soil layers is as we expect

call nc_check(nf90_inq_dimid(iunit, 'num_soil_layers', dimIDs(1)), &
                  'static_init_model','inq_dimid num_soil_layers '//trim(noah_netcdf_filename))
call nc_check(nf90_inquire_dimension(iunit, dimIDs(1), len=nLayers), &
                  'static_init_model','inquire_dimension Time '//trim(noah_netcdf_filename))

if (nSoilLayers /= nLayers) then
   write(string1,*) 'Expected ',nSoilLayers,' soil layers ', &
                       trim(noah_netcdf_filename),' has ',nLayers
   call error_handler(E_ERR,'static_init_model',string1,source,revision,revdate)
endif

! convert soil thicknesses to depths
! closer to the center of the earth is an increasingly large negative number
soil_depths(1) = Soil_layer_thickness(1) 
do i = 2,nSoilLayers
   soil_depths(i) = soil_depths(i-1) + Soil_layer_thickness(i)
enddo
soil_depths = -1.0_r8 * soil_depths

! there are only nSoilLayers + 1 different locations

state_loc(0) = set_location(Longitude, Latitude, 0.0_r8, VERTISHEIGHT)
do i = 1,nSoilLayers
   state_loc(i) = set_location(Longitude, Latitude, soil_depths(i), VERTISHEIGHT)
enddo

!---------------------------------------------------------------
! Compile the list of NOAH variables to use in the creation
! of the DART state vector. Required to determine model_size.
!
! Verify all variables are in the NOAH netcdf file.
! Compute the offsets into the state vector for each variable type.
! Record the extent of the variable type in the state vector.

call verify_state_variables( noah_variables, iunit, noah_netcdf_filename, &
                             nfields, variable_table )

index1  = 1
FILL_PROGVAR : do ivar = 1, nfields

   varname                   = trim(variable_table(ivar,1))
   kind_string               = trim(variable_table(ivar,2))
   progvar(ivar)%varname     = varname
   progvar(ivar)%kind_string = kind_string
   progvar(ivar)%dart_kind   = get_raw_obs_kind_index( progvar(ivar)%kind_string )
   progvar(ivar)%dimlens     = 0
   progvar(ivar)%dimnames    = ' '
   progvar(ivar)%maxlevels   = 0

   string2 = trim(noah_netcdf_filename)//' '//trim(varname)

   call nc_check(nf90_inq_varid(iunit, trim(varname), VarID), &
            'static_init_model', 'inq_varid '//trim(string2))

   call nc_check(nf90_inquire_variable(iunit, VarID, dimids=dimIDs, &
                 ndims=progvar(ivar)%numdims, xtype=progvar(ivar)%xtype), &
            'static_init_model', 'inquire '//trim(string2))

   ! If the long_name and/or units attributes are set, get them.
   ! They are not REQUIRED to exist but are nice to use if they are present.

   if( nf90_inquire_attribute(    iunit, VarID, 'long_name') == NF90_NOERR ) then
      call nc_check( nf90_get_att(iunit, VarID, 'long_name' , progvar(ivar)%long_name), &
                  'static_init_model', 'get_att long_name '//trim(string2))
   else
      progvar(ivar)%long_name = varname
   endif

   if( nf90_inquire_attribute(    iunit, VarID, 'units') == NF90_NOERR )  then
      call nc_check( nf90_get_att(iunit, VarID, 'units' , progvar(ivar)%units), &
                  'static_init_model', 'get_att units '//trim(string2))
   else
      progvar(ivar)%units = 'unknown'
   endif

   ! These variables have a Time dimension. We only want the most recent time.

   varsize = 1
   dimlen  = 1
   DimensionLoop : do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string2)
      call nc_check(nf90_inquire_dimension(iunit, dimIDs(i), name=dimname, len=dimlen), &
                                          'static_init_model', string1)

      if ((trim(dimname) == 'Time') .or. (trim(dimname) == 'time')) dimlen = 1
      progvar(ivar)%dimlens( i) = dimlen
      progvar(ivar)%dimnames(i) = dimname
      varsize = varsize * dimlen

   enddo DimensionLoop

   progvar(ivar)%varsize     = varsize
   progvar(ivar)%index1      = index1
   progvar(ivar)%indexN      = index1 + varsize - 1
   index1                    = index1 + varsize      ! sets up for next variable

   if ((debug > 8) .and. do_output()) then
      write(logfileunit,*)
      write(logfileunit,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(logfileunit,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(logfileunit,*) '  units       ',trim(progvar(ivar)%units)
      write(logfileunit,*) '  xtype       ',progvar(ivar)%xtype
      write(logfileunit,*) '  dimnames    ',progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
      write(logfileunit,*) '  dimlens     ',progvar(ivar)%dimlens( 1:progvar(ivar)%numdims)
      write(logfileunit,*) '  numdims     ',progvar(ivar)%numdims
      write(logfileunit,*) '  varsize     ',progvar(ivar)%varsize
      write(logfileunit,*) '  index1      ',progvar(ivar)%index1
      write(logfileunit,*) '  indexN      ',progvar(ivar)%indexN
      write(logfileunit,*) '  dart_kind   ',progvar(ivar)%dart_kind
      write(logfileunit,*) '  kind_string ',progvar(ivar)%kind_string

      write(     *     ,*)
      write(     *     ,*) trim(progvar(ivar)%varname),' variable number ',ivar
      write(     *     ,*) '  long_name   ',trim(progvar(ivar)%long_name)
      write(     *     ,*) '  units       ',trim(progvar(ivar)%units)
      write(     *     ,*) '  xtype       ',progvar(ivar)%xtype
      write(     *     ,*) '  dimnames    ',progvar(ivar)%dimnames(1:progvar(ivar)%numdims)
      write(     *     ,*) '  dimlens     ',progvar(ivar)%dimlens( 1:progvar(ivar)%numdims)
      write(     *     ,*) '  numdims     ',progvar(ivar)%numdims
      write(     *     ,*) '  varsize     ',progvar(ivar)%varsize
      write(     *     ,*) '  index1      ',progvar(ivar)%index1
      write(     *     ,*) '  indexN      ',progvar(ivar)%indexN
      write(     *     ,*) '  dart_kind   ',progvar(ivar)%dart_kind
      write(     *     ,*) '  kind_string ',progvar(ivar)%kind_string
   endif

enddo FILL_PROGVAR

call nc_check(nf90_close(iunit), 'static_init_model', 'close '//trim(noah_netcdf_filename))

model_size = progvar(nfields)%indexN

if ((debug > 8) .and. do_output()) then

   write(*,*)
   do i=1,nSoilLayers
      write(*,*)'soil layer',i,soil_depths(i)
   enddo

   write(*,*)
   do i=0,nSoilLayers
      call write_location(iunit,state_loc(i),charstring=string1)
      write(*,*)'location ',i,' is ',trim(string1)
   enddo

endif

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

if ( .not. module_initialized ) call static_init_model

x = MISSING_R8

end subroutine init_conditions



subroutine adv_1step(x, time)
!------------------------------------------------------------------
! subroutine adv_1step(x, time)
!
! Does a single timestep advance of the model. The input value of
! the vector x is the starting condition and x is updated to reflect
! the changed state after a timestep. The time argument is intent
! in and is used for models that need to know the date/time to 
! compute a timestep, for instance for radiation computations.
! This interface is only called if the namelist parameter
! async is set to 0 in perfect_model_obs of filter or if the 
! program integrate_model is to be used to advance the model
! state as a separate executable. If one of these options
! is not going to be used (the model will only be advanced as
! a separate model-specific executable), this can be a 
! NULL INTERFACE.

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

if ( .not. module_initialized ) call static_init_model

end subroutine adv_1step



function get_model_size()
!------------------------------------------------------------------
!
! Returns the size of the model as an integer. Required for all
! applications.

integer :: get_model_size

if ( .not. module_initialized ) call static_init_model

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

if ( .not. module_initialized ) call static_init_model

time = model_time

end subroutine init_time



subroutine model_interpolate(x, location, itype, obs_val, istatus)
!------------------------------------------------------------------
!
! Given a state vector, a location, and a model state variable type,
! interpolates the state variable field to that location and returns
! the value in obs_val. The istatus variable should be returned as
! 0 unless there is some problem in computing the interpolation in
! which case an alternate value should be returned. The itype variable
! is a model specific integer that specifies the type of field (for
! instance temperature, zonal wind component, etc.). In low order
! models that have no notion of types of variables, this argument can
! be ignored. For applications in which only perfect model experiments
! with identity observations (i.e. only the value of a particular
! state variable is observed), this can be a NULL INTERFACE.

real(r8),            intent(in) :: x(:)
type(location_type), intent(in) :: location
integer,             intent(in) :: itype
real(r8),           intent(out) :: obs_val
integer,            intent(out) :: istatus

if ( .not. module_initialized ) call static_init_model

! This should be the result of the interpolation of a
! given kind (itype) of variable at the given location.
obs_val = MISSING_R8

! The return code for successful return should be 0. 
! Any positive number is an error.
! Negative values are reserved for use by the DART framework.
! Using distinct positive values for different types of errors can be
! useful in diagnosing problems.
istatus = 1

end subroutine model_interpolate



function get_model_time_step()
!------------------------------------------------------------------
!
! Returns the the time step of the model; the smallest increment
! in time that the model is capable of advancing the state in a given
! implementation. This interface is required for all applications.

type(time_type) :: get_model_time_step

if ( .not. module_initialized ) call static_init_model

get_model_time_step = model_time_step

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

integer,             intent(in)            :: index_in
type(location_type), intent(out)           :: location
integer,             intent(out), optional :: var_type

integer :: n, layer, ivar

if ( .not. module_initialized ) call static_init_model

layer = -1

FindIndex : do n = 1,nfields
   if( (progvar(n)%index1 <= index_in) .and. (index_in <= progvar(n)%indexN) ) then
      layer    = index_in - progvar(n)%index1 + 1
      var_type = progvar(n)%dart_kind
      ivar     = n
      exit FindIndex
   endif
enddo FindIndex

if ((debug > 0) .and. do_output()) then
   write(*,*)'get_state_meta_data: index_in is ',index_in
   write(*,*)'get_state_meta_data: ivar     is ',ivar
   write(*,*)'get_state_meta_data: layer    is ',layer
endif

if( layer == -1 ) then
     write(string1,*) 'Problem, cannot find base_offset, index_in is: ', index_in
     call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
endif

if (progvar(ivar)%varsize == 1) layer = 0

location = state_loc(layer)

end subroutine get_state_meta_data



subroutine end_model()
!------------------------------------------------------------------
!
! Does any shutdown and clean-up needed for model. Can be a NULL
! INTERFACE if the model has no need to clean up storage, etc.

! good style ... perhaps you could deallocate stuff (from static_init_model?).
! deallocate(state_loc)
if ( .not. module_initialized ) call static_init_model

end subroutine end_model



function nc_write_model_atts( ncFileID ) result (ierr)
!------------------------------------------------------------------
! TJH 24 Oct 2006 -- Writes the model-specific attributes to a netCDF file.
!     This includes coordinate variables and some metadata, but NOT
!     the model state vector. We do have to allocate SPACE for the model
!     state vector, but that variable gets filled as the model advances.
!
! As it stands, this routine will work for ANY model, with no modification.
!
! The simplest possible netCDF file would contain a 3D field
! containing the state of 'all' the ensemble members. This requires
! three coordinate variables -- one for each of the dimensions 
! [model_size, ensemble_member, time]. A little metadata is useful, 
! so we can also create some 'global' attributes. 
! This is what is implemented here.
!
! Once the simplest case is working, this routine (and nc_write_model_vars)
! can be extended to create a more logical partitioning of the state vector,
! fundamentally creating a netCDF file with variables that are easily 
! plotted. The bgrid model_mod is perhaps a good one to view, keeping
! in mind it is complicated by the fact it has two coordinate systems. 
! There are stubs in this template, but they are only stubs.
!
! TJH 29 Jul 2003 -- for the moment, all errors are fatal, so the
! return code is always '0 == normal', since the fatal errors stop execution.
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

integer :: StateVarDimID    ! netCDF pointer to state variable dimension (model size)
integer :: MemberDimID      ! netCDF pointer to dimension of ensemble    (ens_size)
integer :: TimeDimID        ! netCDF pointer to time dimension           (unlimited)
integer :: nSoilLayersDimID !   ..     ..                                (# soil layers)

integer :: StateVarVarID   ! netCDF pointer to state variable coordinate array
integer :: StateVarID      ! netCDF pointer to 3D [state,copy,time] array

!----------------------------------------------------------------------
! variables for the namelist output
!----------------------------------------------------------------------

character(len=129), allocatable, dimension(:) :: textblock
integer :: LineLenDimID, nlinesDimID, nmlVarID
integer :: nlines, linelen
logical :: has_ncommas_namelist

!----------------------------------------------------------------------
! we are going to need these to record the creation date in the netCDF file.
! This is entirely optional, but nice.
!----------------------------------------------------------------------

character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic
character(len=NF90_MAX_NAME) :: str1

integer :: i

if ( .not. module_initialized ) call static_init_model

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
! and then put into define mode.
!-------------------------------------------------------------------------------

ierr = -1 ! assume things go poorly

call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                     "nc_write_model_atts", "inquire")
call nc_check(nf90_redef(ncFileID), "nc_write_model_atts", "redef")

!-------------------------------------------------------------------------------
! We need the dimension ID for the number of copies/ensemble members, and
! we might as well check to make sure that Time is the Unlimited dimension. 
! Our job is create the 'model size' dimension.
!-------------------------------------------------------------------------------

call nc_check(nf90_inq_dimid(ncid=ncFileID, name="copy", dimid=MemberDimID), &
                            "nc_write_model_atts", "inq_dimid copy")
call nc_check(nf90_inq_dimid(ncid=ncFileID, name="time", dimid= TimeDimID), &
                            "nc_write_model_atts", "inq_dimid time")

if ( TimeDimID /= unlimitedDimId ) then
   write(string1,*)"Time Dimension ID ",TimeDimID, &
                     " should equal Unlimited Dimension ID",unlimitedDimID
   call error_handler(E_ERR,"nc_write_model_atts", string1, source, revision, revdate)
endif

!-------------------------------------------------------------------------------
! Define the model size / state variable dimension / whatever ...
!-------------------------------------------------------------------------------

call nc_check(nf90_def_dim(ncid=ncFileID, name="StateVariable",  &
                           len=model_size, dimid=StateVarDimID), &
                           "nc_write_model_atts", "def_dim state")

call nc_check(nf90_def_dim(ncid=ncFileID, name="nSoilLayers",  &
                           len=nSoilLayers, dimid=nSoilLayersDimID), &
                           "nc_write_model_atts", "def_dim nSoilLayers")

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call DATE_AND_TIME(crdate,crtime,crzone,values)
write(str1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
                  values(1), values(2), values(3), values(5), values(6), values(7)

call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "creation_date" ,str1), &
                          "nc_write_model_atts", "put_att creation_date")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_source"  ,source), &
                          "nc_write_model_atts", "put_att model_source")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revision",revision), &
                          "nc_write_model_atts", "put_att model_revision")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model_revdate" ,revdate), &
                          "nc_write_model_atts", "put_att model_revdate")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "model","noah_1d"), &
                          "nc_write_model_atts", "put_att model")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "calendar",trim(calendar)), &
                          "nc_write_model_atts", "put_att calendar")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "Longitude",Longitude), &
                          "nc_write_model_atts", "put_att Longitude")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "Latitude",Latitude), &
                          "nc_write_model_atts", "put_att Latitude")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "Forcing_Timestep",Forcing_Timestep), &
                          "nc_write_model_atts", "put_att Forcing_Timestep")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "Noahlsm_Timestep",Noahlsm_Timestep), &
                          "nc_write_model_atts", "put_att Noahlsm_Timestep")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "Soil_type_index",Soil_type_index), &
                          "nc_write_model_atts", "put_att Soil_type_index")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "Vegetation_type_index",Vegetation_type_index), &
                          "nc_write_model_atts", "put_att Vegetation_type_index")
call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, "Urban_veg_category",Urban_veg_category), &
                          "nc_write_model_atts", "put_att Urban_veg_category")

!-------------------------------------------------------------------------------
! Determine shape of most important namelist
!-------------------------------------------------------------------------------

call find_textfile_dims('ncommas_vars.nml', nlines, linelen)
if (nlines > 0) then
  has_ncommas_namelist = .true.
else
  has_ncommas_namelist = .false.
endif

if (debug > 0)    print *, 'ncommas namelist: nlines, linelen = ', nlines, linelen

if (has_ncommas_namelist) then
   allocate(textblock(nlines))
   textblock = ''

   call nc_check(nf90_def_dim(ncid=ncFileID, name='nlines', &
                 len = nlines, dimid = nlinesDimID), &
                 'nc_write_model_atts', 'def_dim nlines ')

   call nc_check(nf90_def_var(ncFileID,name='ncommas_in', xtype=nf90_char,    &
                 dimids = (/ linelenDimID, nlinesDimID /),  varid=nmlVarID), &
                 'nc_write_model_atts', 'def_var ncommas_in')
   call nc_check(nf90_put_att(ncFileID, nmlVarID, 'long_name',       &
                 'contents of ncommas_in namelist'), 'nc_write_model_atts', 'put_att ncommas_in')

endif

!-------------------------------------------------------------------------------
! Here is the extensible part. The simplest scenario is to output the state vector,
! parsing the state vector into model-specific parts is complicated, and you need
! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
! complicated part.
!-------------------------------------------------------------------------------

if ( output_state_vector ) then

   !----------------------------------------------------------------------------
   ! Create a variable for the state vector
   !----------------------------------------------------------------------------

  ! Define the state vector coordinate variable and some attributes.
   call nc_check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=NF90_INT, &
                              dimids=StateVarDimID, varid=StateVarVarID), &
                             "nc_write_model_atts", "def_var StateVariable")
   call nc_check(nf90_put_att(ncFileID, StateVarVarID,"long_name","State Variable ID"), &
                             "nc_write_model_atts", "put_att StateVariable long_name")
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "units",     "indexical"), &
                             "nc_write_model_atts", "put_att StateVariable units")
   call nc_check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, model_size /)), &
                             "nc_write_model_atts", "put_att StateVariable valid_range")

   ! Define the actual (3D) state vector, which gets filled as time goes on ... 
   call nc_check(nf90_def_var(ncid=ncFileID, name="state", xtype=NF90_REAL, &
                 dimids = (/ StateVarDimID, MemberDimID, unlimitedDimID /), &
                 varid=StateVarID), "nc_write_model_atts", "def_var state")
   call nc_check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"), &
                             "nc_write_model_atts", "put_att state long_name")

   ! Leave define mode so we can fill the coordinate variable.
   call nc_check(nf90_enddef(ncfileID),"nc_write_model_atts", "state_vector enddef")

   ! Fill the state variable coordinate variable
   call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /)), &
                                    "nc_write_model_atts", "put_var state")

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------

   ! This block is a stub for something more complicated.
   ! Usually, the control for the execution of this block is a namelist variable.
   ! Take a peek at the bgrid model_mod.f90 for a (rather complicated) example.

   call nc_check(nf90_enddef(ncfileID), "nc_write_model_atts", "prognostic enddef")

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------
call nc_check(nf90_sync(ncFileID),"nc_write_model_atts", "sync")

ierr = 0 ! If we got here, things went well.

end function nc_write_model_atts



function nc_write_model_vars( ncFileID, statevec, copyindex, timeindex ) result (ierr)         
!------------------------------------------------------------------
! TJH 24 Oct 2006 -- Writes the model variables to a netCDF file.
!
! TJH 29 Jul 2003 -- for the moment, all errors are fatal, so the
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

integer :: StateVarID

if ( .not. module_initialized ) call static_init_model

!-------------------------------------------------------------------------------
! make sure ncFileID refers to an open netCDF file, 
!-------------------------------------------------------------------------------

ierr = -1 ! assume things go poorly

call nc_check(nf90_inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID), &
                          "nc_write_model_vars", "inquire")

if ( output_state_vector ) then

   call nc_check(nf90_inq_varid(ncFileID, "state", StateVarID), &
                               "nc_write_model_vars", "inq_varid state" )
   call nc_check(nf90_put_var(ncFileID, StateVarID, statevec,  &
                              start=(/ 1, copyindex, timeindex /)), &
                             "nc_write_model_vars", "put_var state")                   

else

   !----------------------------------------------------------------------------
   ! We need to process the prognostic variables.
   !----------------------------------------------------------------------------

   ! This block is a stub for something more complicated.
   ! Usually, the control for the execution of this block is a namelist variable.
   ! Take a peek at the bgrid model_mod.f90 for a (rather complicated) example.
   !
   ! Generally, it is necessary to take the statevec and decompose it into 
   ! the separate prognostic variables. In this (commented out) example,
   ! global_Var is a user-defined type that has components like:
   ! global_Var%ps, global_Var%t, ... etc. Each of those can then be passed
   ! directly to the netcdf put_var routine. This may cause a huge storage
   ! hit, so large models may want to avoid the duplication if possible.

   ! call vector_to_prog_var(statevec, get_model_size(), global_Var)

   ! the 'start' array is crucial. In the following example, 'ps' is a 2D
   ! array, and the netCDF variable "ps" is a 4D array [lat,lon,copy,time]

   ! call nc_check(nf90_inq_varid(ncFileID, "ps", psVarID), &
   !                             "nc_write_model_vars",  "inq_varid ps")
   ! call nc_check(nf90_put_var( ncFileID, psVarID, global_Var%ps, &
   !                             start=(/ 1, 1, copyindex, timeindex /)), &
   !                            "nc_write_model_vars", "put_var ps")

endif

!-------------------------------------------------------------------------------
! Flush the buffer and leave netCDF file open
!-------------------------------------------------------------------------------

call nc_check(nf90_sync(ncFileID), "nc_write_model_vars", "sync")

ierr = 0 ! If we got here, things went well.

end function nc_write_model_vars



subroutine pert_model_state(state, pert_state, interf_provided)
!------------------------------------------------------------------
!
! Perturbs a model state for generating initial ensembles.
! The perturbed state is returned in pert_state.
! A model may choose to provide a NULL INTERFACE by returning
! .false. for the interf_provided argument. This indicates to
! the filter that if it needs to generate perturbed states, it
! may do so by adding an O(0.1) magnitude perturbation to each
! model state variable independently. The interf_provided argument
! should be returned as .true. if the model wants to do its own
! perturbing of states.  The returned pert_state should in any
! case be valid, since it will be read by filter even if 
! interf_provided is .false.

real(r8), intent(in)  :: state(:)
real(r8), intent(out) :: pert_state(:)
logical,  intent(out) :: interf_provided

if ( .not. module_initialized ) call static_init_model

pert_state      = state
interf_provided = .false.

end subroutine pert_model_state




subroutine ens_mean_for_model(ens_mean)
!------------------------------------------------------------------
! Not used in low-order models

real(r8), intent(in) :: ens_mean(:)

if ( .not. module_initialized ) call static_init_model

end subroutine ens_mean_for_model


!==================================================================
! PUBLIC interfaces that aren't required by the DART code but are
! generally useful for other related utility programs.
! (less necessary for small models; generally used for larger models
! with predefined file formats and control structures.)
!==================================================================


subroutine verify_state_variables( state_variables, ncid, filename, ngood, table )
!------------------------------------------------------------------

character(len=*), dimension(:),   intent(in)  :: state_variables
integer,                          intent(in)  :: ncid
character(len=*),                 intent(in)  :: filename
integer,                          intent(out) :: ngood
character(len=*), dimension(:,:), intent(out) :: table

integer :: nrows, ncols, i, VarID
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr

if ( .not. module_initialized ) call static_init_model

nrows = size(table,1)
ncols = size(table,2)

ngood = 0
MyLoop : do i = 1, nrows

   varname    = trim(state_variables(2*i -1))
   dartstr    = trim(state_variables(2*i   ))
   table(i,1) = trim(varname)
   table(i,2) = trim(dartstr)

   if ( table(i,1) == ' ' .and. table(i,2) == ' ' ) exit MyLoop ! Found end of list.

   if ( table(i,1) == ' ' .or. table(i,2) == ' ' ) then
      string1 = 'model_nml:clm_state_variables not fully specified'
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Make sure variable exists in netCDF file

   write(string1,'(''there is no variable '',a,'' in '',a)') trim(varname), trim(filename)
   call nc_check(NF90_inq_varid(ncid, trim(varname), VarID), &
                 'verify_state_variables', trim(string1))

   ! Make sure DART kind is valid

   if( get_raw_obs_kind_index(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'verify_state_variables',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector

   if ((debug > 3) .and. do_output()) then
      write(logfileunit,*)'variable ',i,' is ',trim(table(i,1)), '   ', trim(table(i,2))
      write(     *     ,*)'variable ',i,' is ',trim(table(i,1)), '   ', trim(table(i,2))
   endif

   ngood = ngood + 1
enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'verify_state_variables',string1,source,revision,revdate,text2=string2)
endif

end subroutine verify_state_variables




subroutine get_noah1D_restart_filename( noah1D_restart_filename )
!------------------------------------------------------------------
character(len=*), intent(out) :: noah1d_restart_filename

if ( .not. module_initialized ) call static_init_model

noah1d_restart_filename = noah_netcdf_filename

end subroutine get_noah1D_restart_filename




subroutine noah1d_to_dart_vector(filename, state_vector, restart_time)
!------------------------------------------------------------------
! Reads the current time and state variables from a model data
! file and packs them into a dart state vector.

character(len=*), intent(in)    :: filename
real(r8),         intent(inout) :: state_vector(:)
type(time_type),  intent(out)   :: restart_time

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, ncstart, nccount
character(len=NF90_MAX_NAME) :: varname

integer  :: year,month,day,hours,minutes
integer  :: ncid, ncNdims, dimlen, VarID
integer  :: i, j, indx, ivar, ntimes

real(r8), allocatable, dimension(:)   :: data_1d_array
real(r8), allocatable, dimension(:,:) :: data_2d_array

if ( .not. module_initialized ) call static_init_model

state_vector(:) = MISSING_R8

! Check that the input file exists ...

if ( .not. file_exist(filename) ) then
   write(string1,*) 'file <', trim(filename),'> does not exist.'
   call error_handler(E_ERR,'noah1d_to_dart_vector',string1,source,revision,revdate)
endif

call nc_check(nf90_open(adjustl(filename), NF90_NOWRITE, ncid), &
                   'noah1d_to_dart_vector', 'open '//trim(filename))

restart_time = get_state_time(ncid,filename)

if (do_output()) call print_time(restart_time,'time in restart file '//trim(filename))
if (do_output()) call print_date(restart_time,'date in restart file '//trim(filename))

! Start counting and filling the state vector one item at a time,
! repacking the Nd arrays into a single 1d list of numbers.

do ivar=1, nfields

   ntimes     = -1
   ncstart(:) = -1
   nccount(:) = -1
   varname    = trim(progvar(ivar)%varname)
   string3    = trim(filename)//' '//trim(varname)

   call nc_check(nf90_inq_varid(ncid, varname, VarID), &
            'noah1d_to_dart_vector', 'inq_varid '//trim(string3))

   call nc_check(nf90_inquire_variable(ncid,VarID,dimids=dimIDs,ndims=ncNdims), &
            'noah1d_to_dart_vector', 'inquire '//trim(string3))

   ! Check the rank of the variable

   if ( ncNdims /= progvar(ivar)%numdims ) then
      write(string1, *) 'netCDF rank of '//trim(varname)//' does not match derived type knowledge'
      write(string2, *) 'netCDF rank is ',ncNdims,' expected ',progvar(ivar)%numdims
      call error_handler(E_ERR,'noah1d_to_dart_vector', string1, &
                        source,revision,revdate,text2=string2)
   endif

   ! Check the shape of the variable
   ! making sure we only compare the last timestep ...

   do i = 1,progvar(ivar)%numdims

      write(string1,'(''inquire dimension'',i2,A)') i,trim(string3)
      call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlen), &
            'noah1d_to_dart_vector', string1)

      ncstart(i) = 1
      nccount(i) = dimlen

      if ( progvar(ivar)%dimnames(i) == 'Time' ) then
         ntimes     = dimlen
         ncstart(i) = dimlen
         nccount(i) = 1
      elseif ( dimlen /= progvar(ivar)%dimlens(i) ) then
         write(string1,*) trim(string3),' dim/dimlen ',i,dimlen,' not ',progvar(ivar)%dimlens(i)
         call error_handler(E_ERR,'noah1d_to_dart_vector',string1,source,revision,revdate)
      endif

   enddo

   ! Pack the variable into the DART state vector

   indx = progvar(ivar)%index1

   if (ncNdims == 1) then

      dimlen = progvar(ivar)%dimlens(1)
      allocate(data_1d_array(dimlen))

      call nc_check(nf90_get_var(ncid, VarID, data_1d_array,  &
                     start=ncstart(1:1), count=nccount(1:1)), &
                  'noah1d_to_dart_vector', 'get_var '//trim(string3))

      do i = 1, dimlen
         state_vector(indx) = data_1d_array(i)
         indx = indx + 1
      enddo
      deallocate(data_1d_array)

   elseif (ncNdims == 2) then

      dimlen = progvar(ivar)%dimlens(1)
      allocate(data_2d_array(progvar(ivar)%dimlens(1), progvar(ivar)%dimlens(2)))

      call nc_check(nf90_get_var(ncid, VarID, data_2d_array,  &
                     start=ncstart(1:2), count=nccount(1:2)), &
                  'noah1d_to_dart_vector', 'get_var '//trim(string3))

      do j = 1, progvar(ivar)%dimlens(2)
      do i = 1, progvar(ivar)%dimlens(1)
         state_vector(indx) = data_2d_array(i,j)
         indx = indx + 1
      enddo
      enddo
      deallocate(data_2d_array)

   else

      write(string1, *)'Variable '//trim(varname)//' has ',ncNdims,' dimensions.'
      write(string2, *)'cannot handle that.'
      call error_handler(E_ERR,'noah1d_to_dart_vector', string1, &
                        source,revision,revdate,text2=string2)

   endif

   indx = indx - 1
   if ( indx /= progvar(ivar)%indexN ) then
      write(string1, *)'Variable '//trim(varname)//' filled wrong.'
      write(string2, *)'Should have ended at ',progvar(ivar)%indexN,' actually ended at ',indx
      call error_handler(E_ERR,'noah1d_to_dart_vector', string1, &
                        source,revision,revdate,text2=string2)
   endif

enddo

if (do_output() .and. (debug > 8)) then
   write(*,*)progvar(ivar)%varname, 'newest time is ',ntimes
   do i = 1,size(state_vector)
      write(*,*)'state vector(',i,') is',state_vector(i)
   enddo
endif

end subroutine noah1d_to_dart_vector


 
subroutine dart_vector_to_model_file(state_vector, filename, statedate, adv_to_time)
!------------------------------------------------------------------
! Writes the current time and state variables from a dart state
! vector (1d array) into a ncommas netcdf restart file.
!
real(r8),         intent(in) :: state_vector(:)
character(len=*), intent(in) :: filename
type(time_type),  intent(in) :: statedate
type(time_type),  intent(in), optional :: adv_to_time

integer :: iunit

if ( .not. module_initialized ) call static_init_model

iunit = open_file(trim(filename), form='formatted', action='write')

! FIXME TJH
! unpack state_vector into namelist variables
! convert statedate to namelist variable
! convert adv_to_time to namelist variable if needed

write(iunit,nml=METADATA_NAMELIST)
close(iunit)

end subroutine dart_vector_to_model_file



function get_state_time(ncid, filename)
!------------------------------------------------------------------
! The restart netcdf files have the time of the state.
! We are always using the 'most recent' which is, by defn, the last one.
!
!        Time = UNLIMITED ; // (35039 currently)
!        num_soil_layers = 4 ;
!        DatStrLen = 12 ;
!variables:
!        char Times(Time, DatStrLen) ;
!                Times:description = "UTC time of data output" ;
!                Times:units = "YYYYMMDD HH:mm" ;

type(time_type) :: get_state_time
integer,          intent(in) :: ncid
character(len=*), intent(in) :: filename

character(len=12), allocatable, dimension(:) :: datestring
integer               :: year, month, day, hour, minute
integer               :: DimID, VarID, strlen, ntimes
integer, dimension(2) :: ncstart, nccount

if ( .not. module_initialized ) call static_init_model

! Get the dimensions for the strings of times

call nc_check(nf90_inq_dimid(ncid, 'Time', DimID), &
                  'get_state_time','inq_dimid Time '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=ntimes), &
                  'get_state_time','inquire_dimension Time '//trim(filename))

call nc_check(nf90_inq_dimid(ncid, 'DatStrLen', DimID), &
                  'get_state_time','inq_dimid DatStrLen '//trim(filename))
call nc_check(nf90_inquire_dimension(ncid, DimID, len=strlen), &
                  'get_state_time','inquire_dimension DatStrLen '//trim(filename))

if (strlen /= 12) then
   write(string1,*)"DatStrLen string length ",strlen," /= 12 "
   call error_handler(E_ERR,"get_state_time", string1, source, revision, revdate)
endif

! Get the last Time string

call nc_check(nf90_inq_varid(ncid, 'Times', VarID), &
                   'get_state_time', 'inq_varid Times '//trim(filename))

! for some reason, this did not work ... its a bit more compact than reading the whole thing.
!ncstart = (/ 1, ntimes /)
!nccount = (/ 1,      1 /)
!call nc_check(nf90_get_var(ncid, VarID, datestring, start=ncstart, count=nccount), &
!                  'get_state_time', 'get_var Times '//trim(filename))

allocate(datestring(ntimes))
call nc_check(nf90_get_var(ncid, VarID, datestring), &
                   'get_state_time', 'get_var Times '//trim(filename))

if (debug > 0) write(*,*)'Last time is '//trim(datestring(ntimes))

read(datestring(ntimes),'(i4,i2,i2,i2,i2)')year, month, day, hour, minute

get_state_time = set_date(year, month, day, hours=hour, minutes=minute, seconds=0)

deallocate(datestring)

end function get_state_time



!===================================================================
! End of model_mod
!===================================================================
end module model_mod
