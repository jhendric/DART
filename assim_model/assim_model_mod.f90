module assim_model_mod
!
! <next four lines automatically updated by CVS, do not edit>
! $Source$ 
! $Revision$ 
! $Date$ 
! $Author$ 
!

! This module is used to wrap around the basic portions of existing dynamical models to
! add capabilities needed by the standard assimilation methods.

! NEED TO ADD ON ONLY CLAUSES
use location_mod, only : location_type, get_dist, write_location, read_location, &
                         LocationDims, LocationName, LocationLName
! I've had a problem with putting in the only for time_manager on the pgf90 compiler (JLA).
use time_manager_mod
use utilities_mod, only : get_unit
use types_mod
use model_mod, only : get_model_size, static_init_model, get_state_meta_data, &
   get_model_time_step, model_interpolate, init_conditions, init_time, adv_1step, &
   end_model, model_get_close_states, nc_write_model_atts

private

public static_init_assim_model, init_diag_output, get_model_size, get_closest_state_time_to, &
   get_initial_condition, get_state_meta_data, get_close_states, get_num_close_states, &
   get_model_time, get_model_state_vector, copy_assim_model, advance_state, interpolate, &
   set_model_time, set_model_state_vector, write_state_restart, read_state_restart, &
   output_diagnostics, end_assim_model, assim_model_type, init_diag_input, input_diagnostics, &
   get_diag_input_copy_meta_data, init_assim_model, get_state_vector_ptr


! Eventually need to be very careful to implement this to avoid state vector copies which
! will be excruciatingly costly (storage at least) in big models.
type assim_model_type
!   private
   real(r8), pointer :: state_vector(:)
   type(time_type) :: time
   integer :: model_size       ! TJH request
   integer :: copyID           ! TJH request
! would like to include character string to indicate which netCDF variable --
! replace "state" in output_diagnostics ...
end type assim_model_type

! Permanent class storage for model_size
integer :: model_size


type(time_type) :: time_step

! Everybody needs to know about this ... TJH Jan 28, 2003

character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"


contains

!======================================================================


subroutine init_assim_model(state)
!----------------------------------------------------------------------
!
! Allocates storage for an instance of an assim_model_type. With this
! implementation, need to be VERY careful about assigment and maintaining
! permanent storage locations. Need to revisit the best way to do 
! assim_model_copy below.

implicit none

type(assim_model_type), intent(inout) :: state

! Get the model_size from the model
model_size = get_model_size()

allocate(state%state_vector(model_size))
state%model_size = model_size

end subroutine init_assim_model




subroutine static_init_assim_model()
!----------------------------------------------------------------------
! subroutine static_init_assim_model()
!
! Initializes class data for the assim_model. Also calls the static
! initialization for the underlying model. So far, this simply 
! is initializing the position of the state variables as location types.

implicit none

! Change output to diagnostic output block ... 

write(*,*)'assim_model attributes:'
write(*,*)'   ',trim(adjustl(source))
write(*,*)'   ',trim(adjustl(revision))
write(*,*)'   ',trim(adjustl(revdate))

! Call the underlying model's static initialization
call static_init_model()

end subroutine static_init_assim_model



function init_diag_output(FileName, global_meta_data, &
                  copies_of_field_per_time, meta_data_per_copy) result(ncFileID)
!--------------------------------------------------------------------------------
!
! Typical sequence:
! NF90_OPEN             ! create netCDF dataset: enter define mode
!    NF90_def_dim       ! define dimenstions: from name and length
!    NF90_def_var       ! define variables: from name, type, and dims
!    NF90_put_att       ! assign attribute values
! NF90_ENDDEF           ! end definitions: leave define mode
!    NF90_put_var       ! provide values for variable
! NF90_CLOSE            ! close: save updated netCDF dataset
!
! Time is a funny beast ... 
! Many packages decode the time:units attribute to convert the offset to a calendar
! date/time format. Using an offset simplifies many operations, but is not the
! way we like to see stuff plotted. The "approved" calendars are:
! gregorian or standard 
!      Mixed Gregorian/Julian calendar as defined by Udunits. This is the default. 
!  noleap   Modern calendar without leap years, i.e., all years are 365 days long. 
!  360_day  All years are 360 days divided into 30 day months. 
!  julian   Julian calendar. 
!  none     No calendar. 
!
! location is another one ...
!

use typeSizes
use netcdf
implicit none

character(len=*), intent(in) :: FileName, global_meta_data
integer,          intent(in) :: copies_of_field_per_time
character(len=*), intent(in) :: meta_data_per_copy(copies_of_field_per_time)
integer                      :: ncFileID

integer             :: i, model_size, metadata_length
type(location_type) :: state_loc

integer :: StateVarDimID, StateVarVarID     ! for each model parameter/State Variable
integer ::   MemberDimID,   MemberVarID     ! for each "copy" or ensemble member
integer ::     TimeDimID,     TimeVarID
integer :: LocationDimID, LocationVarID
integer :: MetadataDimID, MetadataVarID
integer ::                   StateVarID     ! for ENTIRE Model State

if(.not. byteSizesOK()) then
   print *, "Compiler does not appear to support required kinds of variables."
   stop
end if

model_size      = get_model_size()
metadata_length = LEN(meta_data_per_copy(1))

! Create the file
call check(nf90_create(path = trim(FileName)//".nc", cmode = nf90_clobber, ncid = ncFileID))

! Define the dimensions
call check(nf90_def_dim(ncid=ncFileID, &
             name="metadatalength", len = metadata_length,        dimid = metadataDimID))

call check(nf90_def_dim(ncid=ncFileID, &
             name="StateVariable",  len=model_size,               dimid = StateVarDimID))

call check(nf90_def_dim(ncid=ncFileID, &
             name="locationrank",   len = LocationDims,           dimid = LocationDimID))

call check(nf90_def_dim(ncid=ncFileID, &
             name="copy",           len=copies_of_field_per_time, dimid = MemberDimID))

call check(nf90_def_dim(ncid=ncFileID, &
             name="time",           len = nf90_unlimited,         dimid = TimeDimID))

!-------------------------------------------------------------------------------
! Write Global Attributes 
!-------------------------------------------------------------------------------

call check(nf90_put_att(ncFileID, NF90_GLOBAL, "title", global_meta_data))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "assim_model_source", source ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "assim_model_revision", revision ))
call check(nf90_put_att(ncFileID, NF90_GLOBAL, "assim_model_revdate", revdate ))

!-------------------------------------------------------------------------------
! Create variables and attributes
!-------------------------------------------------------------------------------

!    State ID
call check(nf90_def_var(ncid=ncFileID,name="StateVariable", xtype=nf90_int, &
                                     dimids=StateVarDimID, varid=StateVarVarID))
call check(nf90_put_att(ncFileID, StateVarVarID, "long_name", "State Variable ID"))
call check(nf90_put_att(ncFileID, StateVarVarID, "units",     "nondimensional") )
call check(nf90_put_att(ncFileID, StateVarVarID, "valid_range", (/ 1, model_size /)))


!    Copy ID
call check(nf90_def_var(ncid=ncFileID, name="copy", xtype=nf90_int, dimids=MemberDimID, &
                                                                    varid=MemberVarID))
call check(nf90_put_att(ncFileID, MemberVarID, "long_name", "ensemble member or copy"))
call check(nf90_put_att(ncFileID, MemberVarID, "units",     "nondimensional") )
call check(nf90_put_att(ncFileID, MemberVarID, "valid_range", (/ 1, copies_of_field_per_time /)))


!    Metadata for each Copy
call check(nf90_def_var(ncid=ncFileID,name="CopyMetaData", xtype=nf90_char,    &
                        dimids = (/ metadataDimID, MemberDimID /),  varid=metadataVarID))
call check(nf90_put_att(ncFileID, metadataVarID, "long_name",       &
                        "Metadata for each copy/member"))


!    Time -- the unlimited dimension
call check(nf90_def_var(ncFileID, name="time", xtype=nf90_double, dimids=TimeDimID, &
                                                                  varid =TimeVarID) )
call check(nf90_put_att(ncFileID, TimeVarID, "long_name", "time"))
call check(nf90_put_att(ncFileID, TimeVarID, "calendar", "gregorian"))
call check(nf90_put_att(ncFileID, TimeVarID, "cartesian_axis", "T"))
call check(nf90_put_att(ncFileID, TimeVarID, "axis", "T"))
call check(nf90_put_att(ncFileID, TimeVarID, "units", "days since 0000-00-00 00:00:00"))


!    The State 
call check(nf90_def_var(ncid=ncFileID, name="state", xtype=nf90_double, &
                        dimids = (/ StateVarDimID, MemberDimID, TimeDimID /), varid=StateVarID))
call check(nf90_put_att(ncFileID, StateVarID, "long_name", "model state or fcopy"))

!    The Locations -- are part of the model (some models have multiple grids).
!    They are written by model_mod:nc_write_model_atts

! Leave define mode
call check(nf90_enddef(ncfileID))

!-------------------------------------------------------------------------------
! Fill the dimension variables
! The time variable is filled as time progresses.
! The state variable is filled similarly ...
!-------------------------------------------------------------------------------

call check(nf90_put_var(ncFileID, MemberVarID,   (/ (i,i=1,copies_of_field_per_time) /) ))
call check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ))
call check(nf90_put_var(ncFileID, metadataVarID, meta_data_per_copy ))

i =  nc_write_model_atts(ncFileID)
if ( i < 0 ) then
   print *,'assim_model_mod:nc_write_model_atts  bombed ', i
else if ( i > 0 ) then
   print *,'assim_model_mod:nc_write_model_atts  bombed ', i
endif

!-------------------------------------------------------------------------------
call check(nf90_sync(ncFileID))               ! sync to disk, but leave open
!-------------------------------------------------------------------------------

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      print *,'assim_model_mod:init_diag_output'
      stop
    end if
  end subroutine check  

end function init_diag_output



function init_diag_input(file_name, global_meta_data, model_size, copies_of_field_per_time)
!--------------------------------------------------------------------------
!
! Initializes a model state diagnostic file for input. A file id is
! returned which for now is just an integer unit number.

implicit none

integer :: init_diag_input
character(len = *), intent(in) :: file_name
character(len = *), intent(out) ::  global_meta_data
integer, intent(out) :: model_size, copies_of_field_per_time

integer :: i

init_diag_input = get_unit()
open(unit = init_diag_input, file = file_name)
read(init_diag_input, *) global_meta_data

! Read the model size
read(init_diag_input, *) model_size

! Read the number of copies of field per time
read(init_diag_input, *) copies_of_field_per_time

end function init_diag_input



subroutine get_diag_input_copy_meta_data(file_id, model_size_out, num_copies, &
   location, meta_data_per_copy)
!-------------------------------------------------------------------------
!
! Returns the meta data associated with each copy of data in
! a diagnostic input file. Should be called immediately after 
! function init_diag_input.

implicit none

integer, intent(in) :: file_id, model_size_out, num_copies
type(location_type), intent(out) :: location(model_size_out)
character(len = *) :: meta_data_per_copy(num_copies)

character(len=129) :: header
integer :: i, j

! Should have space checks, etc here
! Read the meta data associated with each copy
do i = 1, num_copies
   read(file_id, *) j, meta_data_per_copy(i)
end do

! Will need other metadata, too; Could be as simple as writing locations
read(file_id, *) header
if(header /= 'locat') then
   write(*, *) 'Error: get_diag_input_copy_meta_data expected to read "locat"'
   stop
endif

! Read in the locations
do i = 1, model_size_out
   location(i) =  read_location(file_id)
end do

end subroutine get_diag_input_copy_meta_data




  function get_closest_state_time_to(assim_model, time)
!----------------------------------------------------------------------
!
! Returns the time closest to the given time that the model can reach
! with its state. Initial implementation just assumes fixed timestep.
! Need to describe potentially more general time-stepping capabilities
! from the underlying model in the long run.

implicit none

type(assim_model_type), intent(in) :: assim_model
type(time_type), intent(in) :: time
type(time_type) :: get_closest_state_time_to

type(time_type) :: model_time, delta_time, time_step

! CAREFUL WITH FLOATING POINT DIVISION AND COMPARISONS

! Get the model time step capabilities
time_step = get_model_time_step()

model_time = assim_model%time
if(model_time > time) then
   write(*, *) 'Error in get_closest_state_to_time: model time > time'
   stop
endif

delta_time = time - model_time

get_closest_state_time_to = (delta_time / time_step) * time_step + model_time

end function get_closest_state_time_to



subroutine get_initial_condition(x)
!----------------------------------------------------------------------
! function get_initial_condition()
!
! Initial conditions . This returns an initial assim_model_type
! which includes both a state vector and a time. Design of exactly where this 
! stuff should come from is still evolving (12 July, 2002) but for now can 
! start at time offset 0 with the initial state .
! Need to carefully coordinate this with the times for observations.

implicit none

type(assim_model_type), intent(inout) :: x

call init_conditions(x%state_vector)

call init_time(x%time)

end subroutine get_initial_condition




subroutine get_close_states(location, radius, number, indices, dist)
!---------------------------------------------------------------------
! subroutine get_close_states(location, radius, number, indices)
!
! Returns a list of indices for model state vector points that are
! within distance radius of the location. Might want to add an option
! to return the distances, too. This is written in a model independent
! form at present, hence it is in assim_model_mod. HOWEVER, for
! efficiency in large models, this will have to be model specific at
! some point. At that time, need a way to test to see if this 
! generic form should be over loaded (how to do this in F90 ) by 
! some model specific method.

implicit none

type(location_type), intent(in) :: location
real(r8), intent(in) :: radius
integer, intent(out) :: number, indices(:)
real(r8), intent(out) :: dist(:)

type(location_type) :: state_loc
integer :: index, i
real(r8) :: this_dist

! If model provides a working get_close_states, use it; otherwise search
! Direct use of model dependent stuff, needs to be automated (F90 can't do this
call model_get_close_states(location, radius, number, indices, dist)

! If number returns as -1, not implemented
if(number == -1) then
   index = 0
   model_size = get_model_size()
   do i = 1, model_size
      call get_state_meta_data(i, state_loc)
      this_dist = get_dist(location, state_loc)
      if(this_dist < radius) then
         index = index + 1
         if(index <= size(indices)) indices(index) = i
         if(index <= size(dist)) dist(index) = this_dist
      end if
   end do
   number = index
endif

! If size has overflowed, indicate this with negative size return
if(number > size(indices) .or. number > size(dist)) then
   number = -1 * number
end if

end subroutine get_close_states



function get_num_close_states(location, radius)
!-----------------------------------------------------------------------
!
! Returns number of state vector points located within distance radius
! of the location.

implicit none

integer :: get_num_close_states
type(location_type), intent(in) :: location
real(r8), intent(in) :: radius

type(location_type) :: state_loc
integer :: i, indices(1)
real(r8) :: dist(1)


! call direct model get close with storage that is too 
! small and get size from this
! model_get_close_states returns -1 if it is not implemented
call model_get_close_states(location, radius, get_num_close_states, indices, dist)

if(get_num_close_states == -1) then
! Do exhaustive search
   get_num_close_states = 0
   do i = 1, model_size
      call get_state_meta_data(i, state_loc)
! INTERESTING NOTE: Because of floating point round-off in comps
! this can give a 'variable' number of num close for certain obs, should fix
      if(get_dist(location, state_loc) < radius) get_num_close_states= get_num_close_states + 1
   end do

endif
   
end function get_num_close_states



function get_model_time(assim_model)
!-----------------------------------------------------------------------
!
! Returns the time component of a assim_model extended state.

implicit none

type(time_type) :: get_model_time
type(assim_model_type), intent(in) :: assim_model

get_model_time = assim_model%time

end function get_model_time



function get_state_vector_ptr(assim_model)
!------------------------------------------------------------------------
!
! Returns a pointer directly into the assim_model state vector storage.

implicit none

real(r8), pointer :: get_state_vector_ptr(:)
type(assim_model_type), intent(in) :: assim_model

get_state_vector_ptr => assim_model%state_vector

end function get_state_vector_ptr





subroutine copy_assim_model(model_out, model_in)
!-------------------------------------------------------------------------
!
! Does a copy of assim_model, should be overloaded to =? Still need to be
! very careful about trying to limit copies of the potentially huge state
! vectors for big models.  Interaction with pointer storage?

implicit none

type(assim_model_type), intent(out) :: model_out
type(assim_model_type), intent(in)  :: model_in

integer :: i

! Need to make sure to copy the actual storage and not just the pointer (verify)
model_out%time       = model_in%time
model_out%model_size = model_in%model_size

do i = 1, model_in%model_size
   model_out%state_vector(i) = model_in%state_vector(i)
end do

end subroutine copy_assim_model





subroutine advance_state(assim_model, num, target_time, asynch)
!-----------------------------------------------------------------------
!
! Advances the model extended state until time is equal (within roundoff?)
! of the target_time. For L96 this is relatively straightforward with 
! fixed time steps, etc.

implicit none

integer, intent(in) :: num
type(assim_model_type), intent(inout) :: assim_model(num)
type(time_type), intent(in) :: target_time
logical, intent(in) :: asynch

type(time_type) :: model_time, time_step

integer :: seconds, days, i, len, control_unit, ic_file_unit, ud_file_unit

character(len = 26), dimension(num) :: ic_file_name, ud_file_name 
character(len = 128) :: input_string

! NEED TO BE CAREFUL ABOUT FLOATING POINT TESTS: Being sloppy here

! If none of these needs advancing just return
do i = 1, num
   model_time = get_model_time(assim_model(i))
   if(model_time /= target_time) goto 10
end do
return

! Loop through each model state and advance
10 do i = 1, num
   write(*, *) 'advancing model state ', i
   model_time = get_model_time(assim_model(i))

! Check for time error; use error handler when available
! TJH -- cannot write time_type to stdout -- they have private
!        components and the compiler balks. time_manager_mod
!        provides a routine to return the seconds and days
!        from a time type. 



   if(model_time > target_time) then
      write(*, *) 'Error in advance_state, target_time before model_time'
      call get_time(model_time,seconds,days)
      write(*, *) 'Model_time  (days, seconds) ', days, seconds
      call get_time(target_time,seconds,days)
      write(*, *) 'Target time (days, seconds) ', days, seconds
      stop
   endif

! Two blocks here for now, one for single executable, one for asynch multiple execs

!------------- Block for single executable ----------------------------
! At some point probably need to push the determination of the time back
! into the model itself and out of assim_model
   if(.not. asynch) then
      time_step = get_model_time_step()
      call get_time(time_step, seconds, days)
      do while(model_time < target_time)
         call adv_1step(assim_model(i)%state_vector, model_time)
         model_time = model_time + time_step
         call get_time(model_time, seconds, days)
      end do

! Set the time to updated value
      assim_model(i)%time = model_time
   else
!-------------- End single executable block ---------------------------


!-------------- Block for multiple asynch executables -----------------

! Loop to write out state for each member to separate file, preface with target time
      if(i < 10) then
         write(ic_file_name(i), 11) 'assim_model_state_ic', i
         write(ud_file_name(i), 11) 'assim_model_state_ud', i
      else if(i < 100) then
         write(ic_file_name(i), 21) 'assim_model_state_ic', i
         write(ud_file_name(i), 21) 'assim_model_state_ud', i
      else if(i < 1000) then
         write(ic_file_name(i), 31) 'assim_model_state_ic', i
         write(ud_file_name(i), 31) 'assim_model_state_ud', i
      else if(i < 10000) then
         write(ic_file_name(i), 41) 'assim_model_state_ic', i
         write(ud_file_name(i), 41) 'assim_model_state_ud', i
      else 
         write(*, *) 'advance_state in assim_model_mod needs ensemble size < 10000'
         stop
      endif

 11   format(a21, i1)
 21   format(a21, i2)
 31   format(a21, i3)
 41   format(a21, i4)
      write(*, *) 'ic and ud files ', i, ic_file_name(i), ud_file_name(i)

! Output the destination time followed by assim_model_state 
      ic_file_unit = get_unit()
      open(unit = ic_file_unit, file = ic_file_name(i))
! Write the time to which to advance
      call write_time(ic_file_unit, target_time)
! Write the assim model extended state   
      call write_state_restart(assim_model(i), ic_file_unit)
      close(ic_file_unit)
   endif

!-------------- End of multiple async executables block ---------------

end do


! Also need synchronization block at the end for the asynch
if(asynch) then
! Write out the file names to a control file
   control_unit = get_unit()
   open(unit = control_unit, file = 'filter_control')
   write(control_unit, *) num
   do i = 1, num
      write(control_unit, 51) ic_file_name(i)
      write(control_unit, 51) ud_file_name(i)
 51   format(a26)
   end do
   close(control_unit)

! Create the file async_may_go to allow the async model integrations
   control_unit = get_unit()
   open(unit = control_unit, file = 'async_may_go')
   write(control_unit, *) 1
   close(control_unit)

! Suspend on a read from standard in for integer value
   do 
      read(*, *) input_string
      if(trim(input_string) == 'All_done:Please_proceed') exit
! Following line can allow diagnostic pass through of output
      write(*, *) 'ECHO:', input_string
   end do


write(*, *) 'got clearance to proceed in advance_state'

! All should be done, read in the states and proceed
   do i = 1, num
      ud_file_unit = get_unit()
      open(unit = ud_file_unit, file = ud_file_name(i))
      call read_state_restart(assim_model(i), ud_file_unit)
      close(ud_file_unit)
   end do

end if


end subroutine advance_state



function interpolate(x, location, type)
!---------------------------------------------------------------------
!
! Interpolates from the state vector in an assim_model_type to the
! location. Will need to be generalized for more complex state vector
! types. It might be better to be passing an assim_model_type with
! the associated time through here, but that requires changing the
! entire observation side of the class tree. Reconsider this at a 
! later date (JLA, 15 July, 2002). Type for now is an integer that
! specifies what sort of variable from the model should be interpolated.

implicit none

real(r8) :: interpolate
real(r8), intent(in) :: x(:)
type(location_type), intent(in) :: location
integer, intent(in) :: type

interpolate = model_interpolate(x, location, type)

end function interpolate



subroutine set_model_time(assim_model, time)
!-----------------------------------------------------------------------
!
! Sets the time in an assim_model type

implicit none

type(assim_model_type), intent(inout) :: assim_model
type(time_type), intent(in) :: time

assim_model%time = time

end subroutine set_model_time



subroutine set_model_state_vector(assim_model, state)
!-----------------------------------------------------------------------
!
! Sets the state vector part of an assim_model_type

implicit none

type(assim_model_type), intent(inout) :: assim_model
real(r8), intent(in) :: state(:)

! Check the size for now
if(size(state) /= get_model_size()) then
   write(*, *) 'Error: Input state vector is wrong size in set_model_state_vector'
   stop
endif

assim_model%state_vector = state

end subroutine set_model_state_vector



subroutine write_state_restart(assim_model, file)
!----------------------------------------------------------------------
!
! Write a restart file given a model extended state and a unit number 
! opened to the restart file. (Need to reconsider what is passed to 
! identify file or if file can even be opened within this routine).

implicit none

type (assim_model_type), intent(in) :: assim_model
integer, intent(in) :: file

integer :: days, seconds

! This needs to be done more carefully, consider this an extended stub
! Write time first
call write_time(file, assim_model%time)

! Write the state vector
write(file, *) assim_model%state_vector

end subroutine write_state_restart



subroutine read_state_restart(assim_model, file)
!----------------------------------------------------------------------
!
! Read a restart file given a unit number (see write_state_restart)

implicit none

type(assim_model_type), intent(out) :: assim_model
integer, intent(in) :: file

integer :: seconds, days

! Read the time
assim_model%time = read_time(file)

! Read the state vector
read(file, *) assim_model%state_vector

end subroutine read_state_restart



subroutine output_diagnostics(ncFileID, state, copy_index)
!-------------------------------------------------------------------
! Outputs the "state" to the supplied netCDF file. 
!
! the time, and an optional index saying which
! copy of the metadata this state is associated with.
!
! ncFileID       the netCDF file identifier
! state          the copy of the state vector
! copy_index     which copy of the state vector (ensemble member ID)
!
! TJH 28 Aug 2002 original netCDF implementation 
! TJH  7 Feb 2003 [created time_manager_mod:nc_get_tindex] 
!     substantially modified to handle time in a much better manner
!      

use typeSizes
use netcdf
implicit none

integer,                intent(in) :: ncFileID
type(assim_model_type), intent(in) :: state
integer, optional,      intent(in) :: copy_index

integer :: nDimensions, nVariables, nAttributes, unlimitedDimID, StateVarID
integer :: i,ierr, timeindex, copyindex

if (.not. present(copy_index) ) then     ! we are dependent on the fact
   copyindex = 1                         ! there is a copyindex == 1
else                                     ! if the optional argument is
   copyindex = copy_index                ! not specified, we'd better
endif                                    ! have a backup plan

timeindex = nc_get_tindex(ncFileID, state%time)
if ( timeindex < 0 ) then
   write(*,*)'ERROR: ',trim(adjustl(source))
   write(*,*)'ERROR: output_diagnostics: model%time not in netcdf file'
   write(*,*)'ERROR: model%time : '
   call write_time(6,state%time) ! hardwired unit will change with error handler
   write(*,*)'ERROR: netcdf file ID ',ncFileID
   stop
endif

call check(NF90_Inquire(ncFileID, nDimensions, nVariables, nAttributes, unlimitedDimID))
call check(NF90_inq_varid(ncFileID, "state", StateVarID)) ! Get state Variable ID
call check(NF90_put_var(ncFileID, StateVarID, state%state_vector, start=(/ 1, copyindex, timeindex /)))

call check(NF90_sync(ncFileID))

contains

  ! Internal subroutine - checks error status after each netcdf, prints 
  !                       text message each time an error code is returned. 
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
    end if
  end subroutine check  

end subroutine output_diagnostics



subroutine input_diagnostics(file_id, state, copy_index)
!------------------------------------------------------------------
!
! Reads in diagnostic state output from file_id for copy_index
! copy. Need to make this all more rigorously enforced.

implicit none

integer, intent(in) :: file_id
! MAYBE SHOULDN'T use assim model type here, but just state and time ?
type(assim_model_type), intent(inout) :: state
integer, intent(out) :: copy_index

character(len=5) :: header

! Read in the time
state%time = read_time(file_id)

! Read in the copy index
read(file_id, *) header
if(header /= 'fcopy')  then
   write(*, *) 'Error: expected "copy" in input_diagnostics'
   write(*, *) 'header read was: ', header
   stop
endif

read(file_id, *) copy_index

! Read in the state vector
read(file_id, *) state%state_vector

end subroutine input_diagnostics




subroutine end_assim_model()
!--------------------------------------------------------------------
!
! Closes down assim_model; nothing to do for L96

implicit none

call end_model()

end subroutine end_assim_model



function get_model_state_vector(assim_model)
!--------------------------------------------------------------------
!
! Returns the state vector component of an assim_model extended state.

real(r8) :: get_model_state_vector(model_size)
type(assim_model_type), intent(in) :: assim_model

get_model_state_vector = assim_model%state_vector

end function get_model_state_vector


!
!===================================================================
! End of assim_model_mod
!===================================================================
!
end module assim_model_mod
