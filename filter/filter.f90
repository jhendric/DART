! Data Assimilation Research Testbed -- DART
! Copyright 2004, Data Assimilation Initiative, University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

program filter

! <next four lines automatically updated by CVS, do not edit>
! $Source$
! $Revision$
! $Date$
! $Author$
!

use        types_mod, only : r8, missing_r
use obs_sequence_mod, only : obs_sequence_type, write_obs_sequence, &
   read_obs_sequence, get_num_obs_sets, get_obs_sequence_time, &
   get_num_obs_in_set, get_expected_obs, get_diag_obs_err_cov, &
   get_obs_values, obs_sequence_def_copy, inc_num_obs_copies, &
   set_single_obs_value, get_num_close_states, get_close_states
use time_manager_mod, only : time_type, set_time, print_time, operator(/=), &
   operator(>)
use    utilities_mod, only :  get_unit, open_file, close_file, register_module, &
                              check_nml_error, file_exist, error_handler, E_ERR, &
                              logfileunit, initialize_utilities, finalize_utilities
use  assim_model_mod, only : assim_model_type, static_init_assim_model, &
   get_model_size, get_closest_state_time_to, &
   advance_state, set_model_time, get_model_time, init_diag_output, &
   output_diagnostics, finalize_diag_output, init_assim_model, get_state_vector_ptr, &
   write_state_restart, read_state_restart, get_state_meta_data, &
   binary_restart_files, aoutput_diagnostics, aread_state_restart, &
   aget_closest_state_time_to, awrite_state_restart, Aadvance_state, pert_model_state
use   random_seq_mod, only : random_seq_type, init_random_seq, random_gaussian
use  assim_tools_mod, only : obs_increment, update_from_obs_inc, assim_tools_init
use   cov_cutoff_mod, only : comp_cov_factor
use     location_mod, only : location_type
use   reg_factor_mod, only : comp_reg_factor
use         sort_mod, only : sort

implicit none

! CVS Generated file description for error handling, do not edit
character(len=128) :: &
source   = "$Source$", &
revision = "$Revision$", &
revdate  = "$Date$"

type(obs_sequence_type) :: seq, prior_seq, posterior_seq
type(time_type)         :: time1, time2
type(random_seq_type)   :: random_seq


integer :: i, j, k, ind, iunit, prior_obs_unit, posterior_obs_unit, io
integer :: prior_state_unit, posterior_state_unit, num_obs_in_set, ierr
integer :: PriorStateUnit, PosteriorStateUnit
integer :: lji, meta_data_size
integer :: output_ens_mean_index, output_ens_spread_index
integer :: model_size, num_obs_sets
integer :: grp_size, grp_bot, grp_top, group, num_greater_1
real(r8) :: reg_factor, median
real(r8), allocatable :: sum_reg_factor(:, :), reg_factor_series(:, :, :)
real(r8), allocatable :: regress(:), a_returned(:)

! Storage for direct access to ensemble state vectors
real(r8),        allocatable :: ens(:, :), ens_mean(:), ens_spread(:), x(:)
type(time_type), allocatable :: ens_time(:)
type(time_type)              :: ens_mean_time, ens_spread_time, x_time

real(r8), allocatable  :: obs_inc(:), ens_inc(:), ens_obs(:), swath(:)
real(r8), allocatable  :: obs_err_cov(:), obs(:)
real(r8)               :: cov_factor, mean_inc, sd_ratio
character(len = 129), allocatable   :: ens_copy_meta_data(:)

logical :: interf_provided

! Storage with fixed size for observation space diagnostics
real(r8) :: ges(1), anl(1)

! Set a reasonable upper bound on number of close states, will be increased if needed
integer, parameter    :: first_num_close = 100000
integer               :: num_close_ptr(1) 
integer,  allocatable :: close_ptr(:, :)         ! First element size should be 1
real(r8), allocatable ::  dist_ptr(:, :)         ! First element size should be 1

! Test storage for variance ratio
real(r8) :: var_ratio_sum, var_ratio

! Temporary storage to allow decent CAM initial ensemble perturbations
integer :: var_type
type(location_type) :: location

! Temporary storage to test adaptive error capability; should be moved to assim_tools
integer :: slope_index = 0

!----------------------------------------------------------------
! Namelist input with default values
!
integer  :: async = 0, ens_size = 20
real(r8) :: cutoff      = 0.2_r8
real(r8) :: cov_inflate = 1.0_r8
logical  :: start_from_restart = .false., output_restart = .false.
! if init_time_days and seconds are negative initial time is 0, 0
! for no restart or comes from restart if restart exists
integer  :: init_time_days    = 0
integer  :: init_time_seconds = 0
! Control diagnostic output for state variables
logical  :: output_state_ens_mean = .true., output_state_ens_spread = .true.
integer  :: num_output_ens_members = 0
integer  :: output_interval = 1
integer  :: num_groups = 1
real(r8) :: confidence_slope = 0.0_r8
logical  :: get_mean_reg = .false., get_median_reg = .false.

character(len = 129) :: obs_sequence_file_name = "obs_seq.in", &
                        restart_in_file_name = 'filter_ics', &
                        restart_out_file_name = 'filter_restart'

logical :: output_obs_diagnostics = .false.

namelist /filter_nml/async, ens_size, cutoff, cov_inflate, &
   start_from_restart, output_restart, &
   obs_sequence_file_name, restart_in_file_name, restart_out_file_name, &
   init_time_days, init_time_seconds, output_state_ens_mean, &
   output_state_ens_spread, num_output_ens_members, output_interval, &
   num_groups, confidence_slope, output_obs_diagnostics, get_mean_reg, get_median_reg

!----------------------------------------------------------------
! Start of the routine
!----------------------------------------------------------------

call initialize_utilities
call register_module(source,revision,revdate)
write(logfileunit,*)'STARTING filter ...'
call assim_tools_init()

! Begin by reading the namelist input
if(file_exist('input.nml')) then
   iunit = open_file('input.nml', action = 'read')
   ierr = 1
   do while(ierr /= 0)
      read(iunit, nml = filter_nml, iostat = io, end = 11)
      ierr = check_nml_error(io, 'filter_nml')
   enddo
 11 continue
   call close_file(iunit)
endif
write(logfileunit, nml=filter_nml)

! Now know the ensemble size; allocate all the storage
write(*, *) 'the ensemble size is ', ens_size
allocate(obs_inc(ens_size), ens_inc(ens_size), ens_obs(ens_size), swath(ens_size), &
   ens_copy_meta_data(ens_size + 2), regress(num_groups), a_returned(num_groups))

! Set an initial size for the close state pointers
allocate(close_ptr(1, first_num_close), dist_ptr(1, first_num_close))

! Input the obs_sequence
iunit = get_unit()
open(unit = iunit, file = obs_sequence_file_name)
seq = read_obs_sequence(iunit)
close(iunit)

! Count of number of sets in the sequence
num_obs_sets = get_num_obs_sets(seq)

! Copy just the definitions part of the sequence to the two output obs sequences
if(output_obs_diagnostics) then
   call obs_sequence_def_copy(prior_seq, seq)
   call obs_sequence_def_copy(posterior_seq, seq)
endif

! Set up the metadata for the output ensemble observations
do i = 1, num_output_ens_members
   if(i < 10000) then
      write(ens_copy_meta_data(i), '(a15,1x,i6)') 'ensemble member', i
   else
      write(*, *) 'output metadata in filter needs ensemble size < 10000'
      call error_handler(E_ERR,'filter', 'output metadata in filter needs ensemble size < 10000', &
           source, revision, revdate)
   endif
end do

meta_data_size = num_output_ens_members
if(output_state_ens_mean) then
   meta_data_size = meta_data_size + 1
   ens_copy_meta_data(meta_data_size) = 'ensemble mean'
   output_ens_mean_index = meta_data_size
endif
if(output_state_ens_spread) then
   meta_data_size = meta_data_size + 1
   ens_copy_meta_data(meta_data_size) = 'ensemble spread'
   output_ens_spread_index = meta_data_size
endif

!  to add the ens_mean to the output in addition to the ensemble members
if(output_obs_diagnostics) then
   call inc_num_obs_copies(prior_seq, ens_size + 1, ens_copy_meta_data)
   call inc_num_obs_copies(posterior_seq, ens_size + 1, ens_copy_meta_data)
endif

! Initialize the model class data now that obs_sequence is all set up
call static_init_assim_model()
model_size = get_model_size()
allocate(ens(ens_size, model_size), ens_time(ens_size), ens_mean(model_size))

! Allocate storage for ensemble spread if needed
if(output_state_ens_spread) allocate(ens_spread(model_size))

! Set up diagnostic output for model state, if output is desired
if(  output_state_ens_spread .or. output_state_ens_mean .or. &
    ( num_output_ens_members > 0 ) ) then
   PriorStateUnit     = init_diag_output('Prior_Diag', &
                           'prior ensemble state', meta_data_size, ens_copy_meta_data)
   PosteriorStateUnit = init_diag_output('Posterior_Diag', &
                           'posterior ensemble state', meta_data_size, ens_copy_meta_data)
endif

! Set a time type for initial time if namelist inputs are not negative
if(init_time_days >= 0) then
   time1 = set_time(init_time_seconds, init_time_days)
else
   time1 = set_time(0, 0)
endif

!------------------- Read restart if requested ----------------------

if(start_from_restart) then
   iunit = get_unit()
   if (binary_restart_files ) then
      open(unit = iunit, file = restart_in_file_name, form = "unformatted")
   else
      open(unit = iunit, file = restart_in_file_name)
   endif

   do i = 1, ens_size
      write(*, *) 'trying to read restart ', i
      if (binary_restart_files ) then
         call aread_state_restart(ens_time(i), ens(i, :), iunit, "unformatted")
      else
         call aread_state_restart(ens_time(i), ens(i, :), iunit)
      endif

      ! If init_time_days and init_time_seconds are not < 0, set time to them
      if(init_time_days >= 0) ens_time(i) = time1
   end do
   close(iunit)

   !-----------------  Restart read in --------------------------------
else
   !-----  Block to do cold start initialization of ensembles ----------
   ! Initialize the control and ensemble states and set up direct pointers

   ! WARNING: THIS IS COUNTERINTUITIVE: IF START FROM RESTART IS FALSE,
   ! STILL USE A RESTART FILE TO GET SINGLE CONTROL RUN TO PERTURB AROUND.
   allocate(x(model_size))
   iunit = get_unit()
   if (binary_restart_files ) then
      open(unit = iunit, file = restart_in_file_name, form = "unformatted")
   else
      open(unit = iunit, file = restart_in_file_name)
   endif

   ! Get the initial condition
   if (binary_restart_files ) then
      call aread_state_restart(x_time, x, iunit, "unformatted")
   else
      call aread_state_restart(x_time, x, iunit)
   endif
   close(iunit)

   ! Initialize a repeatable random sequence for perturbations
   call init_random_seq(random_seq)

   ! Perturb for ensembles; 
   do i = 1, ens_size
      call pert_model_state(x, ens(i, :), interf_provided)
      ! If model does not provide a perturbing interface, do it here with uniform 0.002
      if(.not. interf_provided) then
         do j = 1, model_size
            ens(i, j) = random_gaussian(random_seq, x(j), 0.002_r8) 
         end do
      endif
      ! Set time to 0, 0 if none specified, otherwise to specified
      ens_time(i) = time1
   end do
   deallocate(x)
   !-------------------- End of cold start ensemble initialization block ------
endif

! Temporary print of initial model time
write(logfileunit, *) 'initial model time is '
write(     *     , *) 'initial model time is '
call print_time(ens_time(1))

! Advance the model and ensemble to the closest time to the next
! available observations (need to think hard about these model time interfaces).

AdvanceTime : do i = 1, num_obs_sets

   call get_obs_sequence_time(seq, i, time1)
   write(*, *) ' '
   write(*, *) 'time of obs set ', i
   call print_time(time1)

   ! If the model time is past the obs set time, just need to skip???
   if(ens_time(1) > time1) cycle AdvanceTime

   time2 = aget_closest_state_time_to(ens_time(1), time1)

   ! Advance all the ensembles (to the time of the first ensemble)
   if(time2 /= ens_time(1)) call Aadvance_state(ens_time, ens, ens_size, time2, async)

   ! Tag the ensemble mean and spread with the current time
   ens_mean_time   = ens_time(1)
   ens_spread_time = ens_time(1)

   ! Inflate each group separately;  Divide ensemble into num_groups groups
   grp_size = ens_size / num_groups
   Group_inflate: do group = 1, num_groups
      grp_bot = (group - 1) * grp_size + 1
      grp_top = grp_bot + grp_size - 1
      do k = 1, model_size
         ens_mean(k) = sum(ens(grp_bot:grp_top, k)) / grp_size
         do j = grp_bot, grp_top
            ens(j, k) = ens_mean(k) + (ens(j, k) - ens_mean(k)) * sqrt(cov_inflate)
         end do
      end do
   end do Group_inflate

   ! Need global ensemble mean for diagnostics after this
   do k = 1, model_size
      ens_mean(k) = sum(ens(:, k)) / ens_size
   end do

   ! Output state diagnostics as required: NOTE: Prior has been inflated
   if(i / output_interval * output_interval == i) then
      do j = 1, num_output_ens_members
         call aoutput_diagnostics( PriorStateUnit, ens_time(j), ens(j, :), j)
      end do
   end if

   ! Output ensemble mean if requested
   if(output_state_ens_mean .and. i / output_interval * output_interval == i) &
      call aoutput_diagnostics(PriorStateUnit, ens_mean_time, ens_mean, output_ens_mean_index)

   ! Compute and output ensemble spread if requested
   if(output_state_ens_spread  .and. i / output_interval * output_interval == i) then
      do k = 1, model_size
         ens_spread(k) = get_ens_spread(ens, ens_mean(k), ens_size, k)
      end do
      call aoutput_diagnostics(PriorStateUnit, ens_spread_time, ens_spread, output_ens_spread_index)
   endif

   ! How many observations in this set
   num_obs_in_set = get_num_obs_in_set(seq, i)

   ! Allocate storage for the ensemble priors for this number of observations
   allocate(obs_err_cov(num_obs_in_set), obs(num_obs_in_set)) 

   ! Storage for diagnosing fixed obs set reg_factor
   if(i == 1) then
      if(get_mean_reg) then
         allocate(sum_reg_factor(num_obs_in_set, model_size))
         sum_reg_factor = 0.0_r8
      endif
      if(get_median_reg) then
         allocate(reg_factor_series(num_obs_in_set, model_size, num_obs_sets))
         reg_factor_series = 0.0_r8
      endif
   endif

   ! Get the observational error covariance (diagonal at present)
   call get_diag_obs_err_cov(seq, i, obs_err_cov)

   ! Get the observations; from copy 1 for now
   call get_obs_values(seq, i, obs, 1)

   ! Output the ensemble mean and all ensemble members at observation space 
   if(output_obs_diagnostics) then
      do j = 1, num_obs_in_set

         ! then output each ensemble member
         do k = 1, ens_size
            call get_expected_obs(seq, i, ens(k, :), ges(1:1), j)
            call set_single_obs_value(prior_seq, i, j, ges(1), k)
         enddo

         ! output ensemble mean last
         k = ens_size + 1
         call get_expected_obs(seq, i, ens_mean, ges(1:1), j)
         call set_single_obs_value(prior_seq, i, j, ges(1), k)

      end do
   endif

   ! Loop through each observation in the set
   Observations : do j = 1, num_obs_in_set
      ! Compute the ensemble prior for this ob
      do k = 1, ens_size
         call get_expected_obs(seq, i, ens(k, :), ens_obs(k:k), j)
      end do

      ! Divide ensemble into num_groups groups
      grp_size = ens_size / num_groups
      Group1: do group = 1, num_groups
         grp_bot = (group - 1) * grp_size + 1
         grp_top = grp_bot + grp_size - 1

         ! Call obs_increment to do observation space
         call obs_increment(ens_obs(grp_bot:grp_top), ens_size/num_groups, obs(j), &
            obs_err_cov(j), obs_inc(grp_bot:grp_top), confidence_slope, a_returned(group))
      end do Group1

      ! Getting close states for each scalar observation for now
222   call get_close_states(seq, i, 2.0*cutoff, num_close_ptr, close_ptr, dist_ptr, j)
      if(num_close_ptr(1) < 0) then
         deallocate(close_ptr, dist_ptr)
         allocate(close_ptr(1, -1 * num_close_ptr(1)), dist_ptr(1, -1 * num_close_ptr(1)))
         goto 222
      endif

      do k = 1, ens_size
         if (ens_obs(k) == missing_r) num_close_ptr(1) = 0
      end do
      if (obs(j) == missing_r) num_close_ptr(1) = 0


      ! Now loop through each close state variable for this observation
      do k = 1, num_close_ptr(1)
         ind = close_ptr(1, k)
         ! Compute distance dependent envelope
         cov_factor = comp_cov_factor(dist_ptr(1, k), cutoff)

         ! Get the ensemble elements for this state variable and do regression
          swath = ens(:, ind)

         ! Loop through the groups
         Group2: do group = 1, num_groups
            grp_bot = (group - 1) * grp_size + 1
            grp_top = grp_bot + grp_size - 1
            call update_from_obs_inc(ens_obs(grp_bot:grp_top), &
               obs_inc(grp_bot:grp_top), swath(grp_bot:grp_top), ens_size/num_groups, &
               a_returned(group), ens_inc(grp_bot:grp_top), regress(group))
         end do Group2

         ! Compute an information factor for impact of this observation on this state
         if(num_groups == 0) then
             reg_factor = 1.0_r8
         else
            reg_factor = comp_reg_factor(num_groups, regress, i, j, ind)
            if(get_mean_reg) sum_reg_factor(j, k) = sum_reg_factor(j, k) + reg_factor
            if(get_median_reg) reg_factor_series(j, k, i) = reg_factor
         endif

         ! COV_FACTOR NOW COMES IN HERE FOR USE WITH GROUPS; JLA 11/19/03
         reg_factor = min(reg_factor, cov_factor)

         ! Do the final update for this state variable
         ens(:, ind) = ens(:, ind) + reg_factor * ens_inc(:)
      end do

   end do Observations

   ! Free up the storage for this obs set
   deallocate(obs_err_cov, obs)

   ! Output posterior diagnostics

   ! Output state diagnostics as requested
   if(i / output_interval * output_interval == i) then
      do j = 1, num_output_ens_members
          call aoutput_diagnostics(     PosteriorStateUnit, ens_time(j), ens(j, :), j)
      end do
   end if

   ! Compute ensemble mean if either mean or spread to be output
   if(output_state_ens_mean .or. output_state_ens_spread  &
      .and. i / output_interval * output_interval == i) then
      do k = 1, model_size
         ens_mean(k) = sum(ens(:, k)) / ens_size
      end do
   endif

   ! Output the ensemble mean of analysis at the observation space
   if(output_obs_diagnostics) then
       do j = 1, num_obs_in_set

          ! write out the obs, guess, and analysis at OBS space
          ! then output each ensemble member
          do k = 1, ens_size
             call get_expected_obs(seq, i, ens(k, :), anl(1:1), j)
             call set_single_obs_value(posterior_seq, i, j, anl(1), k)
          enddo

          ! output ensemble mean last
          k = ens_size + 1
          call get_expected_obs(seq, i, ens_mean, anl(1:1), j)
          call set_single_obs_value(posterior_seq, i, j, anl(1), k)

       end do 
   endif

   ! Output an ensemble mean if requested
   if(output_state_ens_mean  .and. i / output_interval * output_interval == i) &
      call aoutput_diagnostics(PosteriorStateUnit, ens_mean_time, ens_mean, output_ens_mean_index)

   ! Compute and output state_ens_spread if requested
   if(output_state_ens_spread  .and. i / output_interval * output_interval == i) then
      do k = 1, model_size
         ens_spread(k) = get_ens_spread(ens, ens_mean(k), ens_size, k)
      end do
      call aoutput_diagnostics(PosteriorStateUnit, ens_spread_time, ens_spread, output_ens_spread_index)
   endif

end do AdvanceTime

! properly dispose of the diagnostics files

ierr = finalize_diag_output(PriorStateUnit)
ierr = finalize_diag_output(PosteriorStateUnit)

! Initialize the model state space diagnostic output files

! Output the observation space diagnostic files

if(output_obs_diagnostics) then

   prior_obs_unit = get_unit()
   open(unit = prior_obs_unit, file = 'prior_obs_diagnostics')
   call write_obs_sequence(prior_obs_unit, prior_seq)
   close(prior_obs_unit)

   posterior_obs_unit = get_unit()
   open(unit = posterior_obs_unit, file = 'posterior_obs_diagnostics')
   call write_obs_sequence(posterior_obs_unit, posterior_seq)
   close(posterior_obs_unit)
 
endif

! Output a restart file if requested
if(output_restart) then
   iunit = get_unit()
   if (binary_restart_files ) then
      open(unit = iunit, file = restart_out_file_name, form = "unformatted")
      do i = 1, ens_size
         call awrite_state_restart(ens_time(i), ens(i, :), iunit, "unformatted")
      end do
   else
      open(unit = iunit, file = restart_out_file_name)
      do i = 1, ens_size
         call awrite_state_restart(ens_time(i), ens(i, :), iunit)
      end do
   endif
   close(iunit)
endif

! Output the regression factor means if requested
if(get_mean_reg) then
   iunit = get_unit()
   open(unit = iunit, file = 'time_mean_reg')
   write(iunit, *) num_obs_in_set, model_size
   do j = 1, num_obs_in_set
      do i = 1, model_size
         write(iunit, *) j, i, sum_reg_factor(j, i) / num_obs_sets
      end do
   end do
   close(iunit)
endif

! Output the regression factor medians if requested
if(get_median_reg) then
   iunit = get_unit()
   open(unit = iunit, file = 'time_median_reg')
   write(iunit, *) num_obs_in_set, model_size
   do j = 1, num_obs_in_set
      do i = 1, model_size
         reg_factor_series(j, i, :) = sort(reg_factor_series(j, i, :))
         write(iunit, *) j, i, reg_factor_series(j, i, num_obs_sets / 2)
      end do
   end do
   close(iunit)
endif

!===========================================================

write(logfileunit,*)'FINISHED filter.'
write(logfileunit,*)

call finalize_utilities ! closes the log file.

contains



function get_ens_spread(ens, ens_mean, ens_size, index)
!----------------------------------------------------------
!
! Computes the ensemble mean of the index element of the
! state vector.

implicit none

integer,  intent(in) :: ens_size, index
real(r8), intent(in) :: ens_mean
real(r8), intent(in) :: ens(:, :)
real(r8)             :: get_ens_spread

get_ens_spread = sum((ens(:, index) - ens_mean)**2)
get_ens_spread = sqrt(get_ens_spread / (ens_size - 1))

end function get_ens_spread


end program filter
